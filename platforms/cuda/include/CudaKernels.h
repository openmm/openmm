#ifndef OPENMM_CUDAKERNELS_H_
#define OPENMM_CUDAKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "CudaPlatform.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaFFT3D.h"
#include "CudaParameterSet.h"
#include "CudaSort.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/CustomIntegratorUtilities.h"
#include "lepton/CompiledExpression.h"
#include "lepton/ExpressionProgram.h"
#include <cufft.h>

namespace OpenMM {

/**
 * This abstract class defines an interface for code that can compile CUDA kernels.  This allows a plugin to take advantage of runtime compilation
 * when running on recent versions of CUDA.
 */
class CudaCompilerKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CudaCompilerKernel";
    }
    CudaCompilerKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Compile a kernel to PTX.
     *
     * @param source     the source code for the kernel
     * @param options    the flags to be passed to the compiler
     * @param cu         the CudaContext for which the kernel is being compiled
     */
    virtual std::string createModule(const std::string& source, const std::string& flags, CudaContext& cu) = 0;
};

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class CudaCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    CudaCalcForcesAndEnergyKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcForcesAndEnergyKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * This is called at the beginning of each force/energy computation, before calcForcesAndEnergy() has been called on
     * any ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     */
    void beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups);
    /**
     * This is called at the end of each force/energy computation, after calcForcesAndEnergy() has been called on
     * every ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     * @param valid         the method may set this to false to indicate the results are invalid and the force/energy
     *                      calculation should be repeated
     * @return the potential energy of the system.  This value is added to all values returned by ForceImpls'
     * calcForcesAndEnergy() methods.  That is, each force kernel may <i>either</i> return its contribution to the
     * energy directly, <i>or</i> add it to an internal buffer so that it will be included here.
     */
    double finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid);
private:
   CudaContext& cu;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class CudaUpdateStateDataKernel : public UpdateStateDataKernel {
public:
    CudaUpdateStateDataKernel(std::string name, const Platform& platform, CudaContext& cu) : UpdateStateDataKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Get the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    double getTime(const ContextImpl& context) const;
    /**
     * Set the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    void setTime(ContextImpl& context, double time);
    /**
     * Get the positions of all particles.
     *
     * @param positions  on exit, this contains the particle positions
     */
    void getPositions(ContextImpl& context, std::vector<Vec3>& positions);
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containg the particle positions
     */
    void setPositions(ContextImpl& context, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles.
     *
     * @param velocities  on exit, this contains the particle velocities
     */
    void getVelocities(ContextImpl& context, std::vector<Vec3>& velocities);
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containg the particle velocities
     */
    void setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities);
    /**
     * Get the current forces on all particles.
     *
     * @param forces  on exit, this contains the forces
     */
    void getForces(ContextImpl& context, std::vector<Vec3>& forces);
    /**
     * Get the current derivatives of the energy with respect to context parameters.
     *
     * @param derivs  on exit, this contains the derivatives
     */
    void getEnergyParameterDerivatives(ContextImpl& context, std::map<std::string, double>& derivs);
    /**
     * Get the current periodic box vectors.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the current periodic box vectors.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * Create a checkpoint recording the current state of the Context.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(ContextImpl& context, std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(ContextImpl& context, std::istream& stream);
private:
    CudaContext& cu;
};

/**
 * This kernel modifies the positions of particles to enforce distance constraints.
 */
class CudaApplyConstraintsKernel : public ApplyConstraintsKernel {
public:
    CudaApplyConstraintsKernel(std::string name, const Platform& platform, CudaContext& cu) : ApplyConstraintsKernel(name, platform),
            cu(cu), hasInitializedKernel(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Update particle positions to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the distance tolerance within which constraints must be satisfied.
     */
    void apply(ContextImpl& context, double tol);
    /**
     * Update particle velocities to enforce constraints.
     *
     * @param context    the context in which to execute this kernel
     * @param tol        the velocity tolerance within which constraints must be satisfied.
     */
    void applyToVelocities(ContextImpl& context, double tol);
private:
    CudaContext& cu;
    bool hasInitializedKernel;
    CUfunction applyDeltasKernel;
};

/**
 * This kernel recomputes the positions of virtual sites.
 */
class CudaVirtualSitesKernel : public VirtualSitesKernel {
public:
    CudaVirtualSitesKernel(std::string name, const Platform& platform, CudaContext& cu) : VirtualSitesKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Compute the virtual site locations.
     *
     * @param context    the context in which to execute this kernel
     */
    void computePositions(ContextImpl& context);
private:
    CudaContext& cu;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    CudaCalcHarmonicBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcHarmonicBondForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicBondForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force);
private:
    class ForceInfo;
    int numBonds;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    CudaCalcCustomBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomBondForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcCustomBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomBondForce& force);
private:
    class ForceInfo;
    int numBonds;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaParameterSet* params;
    CudaArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    CudaCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcHarmonicAngleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the HarmonicAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
public:
    CudaCalcCustomAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomAngleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcCustomAngleForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomAngleForce this kernel will be used for
     */
    void initialize(const System& system, const CustomAngleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomAngleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomAngleForce& force);
private:
    class ForceInfo;
    int numAngles;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaParameterSet* params;
    CudaArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    CudaCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcPeriodicTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the PeriodicTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const PeriodicTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the PeriodicTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaArray params;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    CudaCalcRBTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcRBTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RBTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const RBTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RBTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const RBTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaArray params1;
    CudaArray params2;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
public:
    CudaCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCMAPTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CMAPTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const CMAPTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CMAPTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    std::vector<int2> mapPositionsVec;
    CudaArray coefficients;
    CudaArray mapPositions;
    CudaArray torsionMaps;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
public:
    CudaCalcCustomTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcCustomTorsionForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const CustomTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomTorsionForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force);
private:
    class ForceInfo;
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaParameterSet* params;
    CudaArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CudaCalcNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcNonbondedForceKernel(name, platform),
            cu(cu), hasInitializedFFT(false), sort(NULL), dispersionFft(NULL), fft(NULL), pmeio(NULL), usePmeStream(false) {
    }
    ~CudaCalcNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param includeDirect  true if direct space interactions should be included
     * @param includeReciprocal  true if reciprocal space interactions should be included
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the NonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const NonbondedForce& force);
    /**
     * Get the parameters being used for PME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the dispersion parameters being used for the dispersion term in LJPME.
     * 
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
    class SortTrait : public CudaSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "(-2147483647-1)";}
        const char* getMaxKey() const {return "2147483647";}
        const char* getMaxValue() const {return "make_int2(2147483647, 2147483647)";}
        const char* getSortKey() const {return "value.y";}
    };
    class ForceInfo;
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncStreamPreComputation;
    class SyncStreamPostComputation;
    CudaContext& cu;
    ForceInfo* info;
    bool hasInitializedFFT;
    CudaArray charges;
    CudaArray sigmaEpsilon;
    CudaArray exceptionParams;
    CudaArray exclusionAtoms;
    CudaArray exclusionParams;
    CudaArray baseParticleParams;
    CudaArray baseExceptionParams;
    CudaArray particleParamOffsets;
    CudaArray exceptionParamOffsets;
    CudaArray particleOffsetIndices;
    CudaArray exceptionOffsetIndices;
    CudaArray globalParams;
    CudaArray cosSinSums;
    CudaArray pmeGrid1;
    CudaArray pmeGrid2;
    CudaArray pmeBsplineModuliX;
    CudaArray pmeBsplineModuliY;
    CudaArray pmeBsplineModuliZ;
    CudaArray pmeDispersionBsplineModuliX;
    CudaArray pmeDispersionBsplineModuliY;
    CudaArray pmeDispersionBsplineModuliZ;
    CudaArray pmeAtomGridIndex;
    CudaArray pmeEnergyBuffer;
    CudaSort* sort;
    Kernel cpuPme;
    PmeIO* pmeio;
    CUstream pmeStream;
    CUevent pmeSyncEvent, paramsSyncEvent;
    CudaFFT3D* fft;
    cufftHandle fftForward;
    cufftHandle fftBackward;
    CudaFFT3D* dispersionFft;
    cufftHandle dispersionFftForward;
    cufftHandle dispersionFftBackward;
    CUfunction computeParamsKernel, computeExclusionParamsKernel;
    CUfunction ewaldSumsKernel;
    CUfunction ewaldForcesKernel;
    CUfunction pmeGridIndexKernel;
    CUfunction pmeDispersionGridIndexKernel;
    CUfunction pmeSpreadChargeKernel;
    CUfunction pmeDispersionSpreadChargeKernel;
    CUfunction pmeFinishSpreadChargeKernel;
    CUfunction pmeDispersionFinishSpreadChargeKernel;
    CUfunction pmeEvalEnergyKernel;
    CUfunction pmeEvalDispersionEnergyKernel;
    CUfunction pmeConvolutionKernel;
    CUfunction pmeDispersionConvolutionKernel;
    CUfunction pmeInterpolateForceKernel;
    CUfunction pmeInterpolateDispersionForceKernel;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;
    double ewaldSelfEnergy, dispersionCoefficient, alpha, dispersionAlpha;
    int interpolateForceThreads;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    bool hasCoulomb, hasLJ, usePmeStream, useCudaFFT, doLJPME, usePosqCharges, recomputeParams, hasOffsets;
    NonbondedMethod nonbondedMethod;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    CudaCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomNonbondedForceKernel(name, platform),
            cu(cu), params(NULL), forceCopy(NULL), system(system), hasInitializedKernel(false) {
    }
    ~CudaCalcCustomNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const CustomNonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomNonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force);
private:
    class ForceInfo;
    void initInteractionGroups(const CustomNonbondedForce& force, const std::string& interactionSource, const std::vector<std::string>& tableTypes);
    CudaContext& cu;
    ForceInfo* info;
    CudaParameterSet* params;
    CudaArray globals;
    CudaArray interactionGroupData, filteredGroupData, numGroupTiles;
    CUfunction interactionGroupKernel, prepareNeighborListKernel, buildNeighborListKernel;
    std::vector<void*> interactionGroupArgs, prepareNeighborListArgs, buildNeighborListArgs;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    double longRangeCoefficient;
    std::vector<double> longRangeCoefficientDerivs;
    bool hasInitializedLongRangeCorrection, hasInitializedKernel, hasParamDerivs, useNeighborList;
    int numGroupThreadBlocks;
    CustomNonbondedForce* forceCopy;
    const System& system;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CudaCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    CudaCalcGBSAOBCForceKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcGBSAOBCForceKernel(name, platform), cu(cu),
            hasCreatedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    void initialize(const System& system, const GBSAOBCForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the GBSAOBCForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force);
private:
    class ForceInfo;
    double prefactor, surfaceAreaFactor, cutoff;
    bool hasCreatedKernels;
    int maxTiles;
    CudaContext& cu;
    ForceInfo* info;
    CudaArray charges;
    CudaArray params;
    CudaArray bornSum;
    CudaArray bornRadii;
    CudaArray bornForce;
    CudaArray obcChain;
    CUfunction computeBornSumKernel;
    CUfunction reduceBornSumKernel;
    CUfunction force1Kernel;
    CUfunction reduceBornForceKernel;
    std::vector<void*> computeSumArgs, force1Args;
};

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system.
 */
class CudaCalcCustomGBForceKernel : public CalcCustomGBForceKernel {
public:
    CudaCalcCustomGBForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomGBForceKernel(name, platform),
            hasInitializedKernels(false), cu(cu), params(NULL), computedValues(NULL), energyDerivs(NULL), energyDerivChain(NULL), system(system) {
    }
    ~CudaCalcCustomGBForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomGBForce this kernel will be used for
     */
    void initialize(const System& system, const CustomGBForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomGBForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomGBForce& force);
private:
    class ForceInfo;
    double cutoff;
    bool hasInitializedKernels, needParameterGradient, needEnergyParamDerivs;
    int maxTiles, numComputedValues;
    CudaContext& cu;
    ForceInfo* info;
    CudaParameterSet* params;
    CudaParameterSet* computedValues;
    CudaParameterSet* energyDerivs;
    CudaParameterSet* energyDerivChain;
    std::vector<CudaParameterSet*> dValuedParam;
    std::vector<CudaArray> dValue0dParam;
    CudaArray longEnergyDerivs;
    CudaArray globals;
    CudaArray valueBuffers;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    std::vector<bool> pairValueUsesParam, pairEnergyUsesParam, pairEnergyUsesValue;
    const System& system;
    CUfunction pairValueKernel, perParticleValueKernel, pairEnergyKernel, perParticleEnergyKernel, gradientChainRuleKernel;
    std::vector<void*> pairValueArgs, perParticleValueArgs, pairEnergyArgs, perParticleEnergyArgs, gradientChainRuleArgs;
    std::string pairValueSrc, pairEnergySrc;
    std::map<std::string, std::string> pairValueDefines, pairEnergyDefines;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
public:
    CudaCalcCustomExternalForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomExternalForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcCustomExternalForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomExternalForce this kernel will be used for
     */
    void initialize(const System& system, const CustomExternalForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomExternalForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomExternalForce& force);
private:
    class ForceInfo;
    int numParticles;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    const System& system;
    CudaParameterSet* params;
    CudaArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system.
 */
class CudaCalcCustomHbondForceKernel : public CalcCustomHbondForceKernel {
public:
    CudaCalcCustomHbondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomHbondForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), donorParams(NULL), acceptorParams(NULL), system(system) {
    }
    ~CudaCalcCustomHbondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomHbondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomHbondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomHbondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomHbondForce& force);
private:
    class ForceInfo;
    int numDonors, numAcceptors;
    bool hasInitializedKernel;
    CudaContext& cu;
    ForceInfo* info;
    CudaParameterSet* donorParams;
    CudaParameterSet* acceptorParams;
    CudaArray globals;
    CudaArray donors;
    CudaArray acceptors;
    CudaArray donorExclusions;
    CudaArray acceptorExclusions;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    std::vector<void*> donorArgs, acceptorArgs;
    const System& system;
    CUfunction donorKernel, acceptorKernel;
};

/**
 * This kernel is invoked by CustomCentroidBondForce to calculate the forces acting on the system.
 */
class CudaCalcCustomCentroidBondForceKernel : public CalcCustomCentroidBondForceKernel {
public:
    CudaCalcCustomCentroidBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomCentroidBondForceKernel(name, platform),
            cu(cu), params(NULL), system(system) {
    }
    ~CudaCalcCustomCentroidBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCentroidBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomCentroidBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCentroidBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force);

private:
    class ForceInfo;
    int numGroups, numBonds;
    bool needEnergyParamDerivs;
    CudaContext& cu;
    ForceInfo* info;
    CudaParameterSet* params;
    CudaArray globals;
    CudaArray groupParticles;
    CudaArray groupWeights;
    CudaArray groupOffsets;
    CudaArray groupForces;
    CudaArray bondGroups;
    CudaArray centerPositions;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    std::vector<void*> groupForcesArgs;
    CUfunction computeCentersKernel, groupForcesKernel, applyForcesKernel;
    const System& system;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
 */
class CudaCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
public:
    CudaCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomCompoundBondForceKernel(name, platform),
            cu(cu), params(NULL), system(system) {
    }
    ~CudaCalcCustomCompoundBondForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCompoundBondForce this kernel will be used for
     */
    void initialize(const System& system, const CustomCompoundBondForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCompoundBondForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force);

private:
    class ForceInfo;
    int numBonds;
    CudaContext& cu;
    ForceInfo* info;
    CudaParameterSet* params;
    CudaArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    const System& system;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system.
 */
class CudaCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    CudaCalcCustomManyParticleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomManyParticleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), params(NULL), system(system) {
    }
    ~CudaCalcCustomManyParticleForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomManyParticleForce this kernel will be used for
     */
    void initialize(const System& system, const CustomManyParticleForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomManyParticleForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force);

private:
    class ForceInfo;
    CudaContext& cu;
    ForceInfo* info;
    bool hasInitializedKernel;
    NonbondedMethod nonbondedMethod;
    int maxNeighborPairs, forceWorkgroupSize, findNeighborsWorkgroupSize;
    CudaParameterSet* params;
    CudaArray particleTypes;
    CudaArray orderIndex;
    CudaArray particleOrder;
    CudaArray exclusions;
    CudaArray exclusionStartIndex;
    CudaArray blockCenter;
    CudaArray blockBoundingBox;
    CudaArray neighborPairs;
    CudaArray numNeighborPairs;
    CudaArray neighborStartIndex;
    CudaArray numNeighborsForAtom;
    CudaArray neighbors;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray> tabulatedFunctions;
    std::vector<void*> forceArgs, blockBoundsArgs, neighborsArgs, startIndicesArgs, copyPairsArgs;
    const System& system;
    CUfunction forceKernel, blockBoundsKernel, neighborsKernel, startIndicesKernel, copyPairsKernel;
    CUdeviceptr globalsPtr;
    CUevent event;
};

/**
 * This kernel is invoked by GayBerneForce to calculate the forces acting on the system.
 */
class CudaCalcGayBerneForceKernel : public CalcGayBerneForceKernel {
public:
    CudaCalcGayBerneForceKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcGayBerneForceKernel(name, platform), cu(cu),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GayBerneForce this kernel will be used for
     */
    void initialize(const System& system, const GayBerneForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the GayBerneForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const GayBerneForce& force);
private:
    class ForceInfo;
    class ReorderListener;
    void sortAtoms();
    CudaContext& cu;
    ForceInfo* info;
    bool hasInitializedKernels;
    int numRealParticles, numExceptions, maxNeighborBlocks;
    GayBerneForce::NonbondedMethod nonbondedMethod;
    CudaArray sortedParticles;
    CudaArray axisParticleIndices;
    CudaArray sigParams;
    CudaArray epsParams;
    CudaArray scale;
    CudaArray exceptionParticles;
    CudaArray exceptionParams;
    CudaArray aMatrix;
    CudaArray bMatrix;
    CudaArray gMatrix;
    CudaArray exclusions;
    CudaArray exclusionStartIndex;
    CudaArray blockCenter;
    CudaArray blockBoundingBox;
    CudaArray neighbors;
    CudaArray neighborIndex;
    CudaArray neighborBlockCount;
    CudaArray sortedPos;
    CudaArray torque;
    std::vector<bool> isRealParticle;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::pair<int, int> > excludedPairs;
    std::vector<void*> framesArgs, blockBoundsArgs, neighborsArgs, forceArgs, torqueArgs;
    CUfunction framesKernel, blockBoundsKernel, neighborsKernel, forceKernel, torqueKernel;
    CUevent event;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomCVForceKernel : public CalcCustomCVForceKernel {
public:
    CudaCalcCustomCVForceKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcCustomCVForceKernel(name, platform),
            cu(cu), hasInitializedListeners(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomCVForce this kernel will be used for
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    void initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy);
    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the CustomCVForce for computing collective variables
     */
    void copyState(ContextImpl& context, ContextImpl& innerContext);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomCVForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomCVForce& force);
private:
    class ForceInfo;
    class ReorderListener;
    CudaContext& cu;
    bool hasInitializedListeners;
    Lepton::ExpressionProgram energyExpression;
    std::vector<std::string> variableNames, paramDerivNames, globalParameterNames;
    std::vector<Lepton::ExpressionProgram> variableDerivExpressions;
    std::vector<Lepton::ExpressionProgram> paramDerivExpressions;
    std::vector<CudaArray> cvForces;
    CudaArray invAtomOrder;
    CudaArray innerInvAtomOrder;
    CUfunction copyStateKernel, copyForcesKernel, addForcesKernel;
};

/**
 * This kernel is invoked by RMSDForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcRMSDForceKernel : public CalcRMSDForceKernel {
public:
    CudaCalcRMSDForceKernel(std::string name, const Platform& platform, CudaContext& cu) : CalcRMSDForceKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the RMSDForce this kernel will be used for
     */
    void initialize(const System& system, const RMSDForce& force);
    /**
     * Record the reference positions and particle indices.
     */
    void recordParameters(const RMSDForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * This is the internal implementation of execute(), templatized on whether we're
     * using single or double precision.
     */
    template <class REAL>
    double executeImpl(ContextImpl& context);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RMSDForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const RMSDForce& force);
private:
    class ForceInfo;
    CudaContext& cu;
    ForceInfo* info;
    double sumNormRef;
    CudaArray referencePos;
    CudaArray particles;
    CudaArray buffer;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class CudaIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    CudaIntegrateVerletStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateVerletStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const VerletIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator);
private:
    CudaContext& cu;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CudaIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    CudaIntegrateLangevinStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateLangevinStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    void initialize(const System& system, const LangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const LangevinIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator);
private:
    CudaContext& cu;
    double prevTemp, prevFriction, prevStepSize;
    CudaArray params;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class CudaIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    CudaIntegrateBrownianStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateBrownianStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the BrownianIntegrator this kernel will be used for
     */
    void initialize(const System& system, const BrownianIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const BrownianIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator);
private:
    CudaContext& cu;
    double prevTemp, prevFriction, prevStepSize;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by VariableVerletIntegrator to take one time step.
 */
class CudaIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
public:
    CudaIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateVariableVerletStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableVerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableVerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    double execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableVerletIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator);
private:
    CudaContext& cu;
    int blockSize;
    CUfunction kernel1, kernel2, selectSizeKernel;
};

/**
 * This kernel is invoked by VariableLangevinIntegrator to take one time step.
 */
class CudaIntegrateVariableLangevinStepKernel : public IntegrateVariableLangevinStepKernel {
public:
    CudaIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateVariableLangevinStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VariableLangevinIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableLangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     * @return the size of the step that was taken
     */
    double execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VariableLangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator);
private:
    CudaContext& cu;
    int blockSize;
    CudaArray params;
    CUfunction kernel1, kernel2, selectSizeKernel;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by CustomIntegrator to take one time step.
 */
class CudaIntegrateCustomStepKernel : public IntegrateCustomStepKernel {
public:
    enum GlobalTargetType {DT, VARIABLE, PARAMETER};
    CudaIntegrateCustomStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateCustomStepKernel(name, platform), cu(cu),
            hasInitializedKernels(false), needsEnergyParamDerivs(false) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the CustomIntegrator this kernel will be used for
     */
    void initialize(const System& system, const CustomIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the CustomIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated.
     *                       On exit, this should specify whether the cached forces are valid at the
     *                       end of the step.
     */
    void execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the CustomIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated.
     *                       On exit, this should specify whether the cached forces are valid at the
     *                       end of the step.
     */
    double computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
    /**
     * Get the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    on exit, this contains the values
     */
    void getGlobalVariables(ContextImpl& context, std::vector<double>& values) const;
    /**
     * Set the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    a vector containing the values
     */
    void setGlobalVariables(ContextImpl& context, const std::vector<double>& values);
    /**
     * Get the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    on exit, this contains the values
     */
    void getPerDofVariable(ContextImpl& context, int variable, std::vector<Vec3>& values) const;
    /**
     * Set the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    a vector containing the values
     */
    void setPerDofVariable(ContextImpl& context, int variable, const std::vector<Vec3>& values);
private:
    class ReorderListener;
    class GlobalTarget;
    class DerivFunction;
    std::string createPerDofComputation(const std::string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator,
        const std::string& forceName, const std::string& energyName, std::vector<const TabulatedFunction*>& functions,
        std::vector<std::pair<std::string, std::string> >& functionNames);
    void prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
    Lepton::ExpressionTreeNode replaceDerivFunctions(const Lepton::ExpressionTreeNode& node, OpenMM::ContextImpl& context);
    void findExpressionsForDerivs(const Lepton::ExpressionTreeNode& node, std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& variableNodes);
    void recordGlobalValue(double value, GlobalTarget target, CustomIntegrator& integrator);
    void recordChangedParameters(ContextImpl& context);
    bool evaluateCondition(int step);
    CudaContext& cu;
    double energy;
    float energyFloat;
    int numGlobalVariables, sumWorkGroupSize;
    bool hasInitializedKernels, deviceGlobalsAreCurrent, modifiesParameters, keNeedsForce, hasAnyConstraints, needsEnergyParamDerivs;
    std::vector<bool> deviceValuesAreCurrent;
    mutable std::vector<bool> localValuesAreCurrent; 
    CudaArray globalValues;
    CudaArray sumBuffer;
    CudaArray summedValue;
    CudaArray uniformRandoms;
    CudaArray randomSeed;
    CudaArray perDofEnergyParamDerivs;
    std::vector<CudaArray> tabulatedFunctions, perDofValues;
    std::map<int, double> savedEnergy;
    std::map<int, CudaArray> savedForces;
    std::set<int> validSavedForces;
    mutable std::vector<std::vector<float4> > localPerDofValuesFloat;
    mutable std::vector<std::vector<double4> > localPerDofValuesDouble;
    std::map<std::string, double> energyParamDerivs;
    std::vector<std::string> perDofEnergyParamDerivNames;
    std::vector<double> localPerDofEnergyParamDerivs;
    std::vector<double> localGlobalValues;
    std::vector<double> initialGlobalVariables;
    std::vector<std::vector<CUfunction> > kernels;
    std::vector<std::vector<std::vector<void*> > > kernelArgs;
    std::vector<void*> kineticEnergyArgs;
    CUfunction randomKernel, kineticEnergyKernel, sumKineticEnergyKernel;
    std::vector<CustomIntegrator::ComputationType> stepType;
    std::vector<CustomIntegratorUtilities::Comparison> comparisons;
    std::vector<std::vector<Lepton::CompiledExpression> > globalExpressions;
    CompiledExpressionSet expressionSet;
    std::vector<bool> needsGlobals;
    std::vector<bool> needsForces;
    std::vector<bool> needsEnergy;
    std::vector<bool> computeBothForceAndEnergy;
    std::vector<bool> invalidatesForces;
    std::vector<bool> merged;
    std::vector<int> forceGroupFlags;
    std::vector<int> blockEnd;
    std::vector<int> requiredGaussian;
    std::vector<int> requiredUniform;
    std::vector<int> stepEnergyVariableIndex;
    std::vector<int> globalVariableIndex;
    std::vector<int> parameterVariableIndex;
    int gaussianVariableIndex, uniformVariableIndex, dtVariableIndex;
    std::vector<std::string> parameterNames;
    std::vector<GlobalTarget> stepTarget;
};

class CudaIntegrateCustomStepKernel::GlobalTarget {
public:
    CudaIntegrateCustomStepKernel::GlobalTargetType type;
    int variableIndex;
    GlobalTarget() {
    }
    GlobalTarget(CudaIntegrateCustomStepKernel::GlobalTargetType type, int variableIndex) : type(type), variableIndex(variableIndex) {
    }
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class CudaApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    CudaApplyAndersenThermostatKernel(std::string name, const Platform& platform, CudaContext& cu) : ApplyAndersenThermostatKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param thermostat the AndersenThermostat this kernel will be used for
     */
    void initialize(const System& system, const AndersenThermostat& thermostat);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     */
    void execute(ContextImpl& context);
private:
    CudaContext& cu;
    int randomSeed;
    CudaArray atomGroups;
    CUfunction kernel;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class CudaApplyMonteCarloBarostatKernel : public ApplyMonteCarloBarostatKernel {
public:
    CudaApplyMonteCarloBarostatKernel(std::string name, const Platform& platform, CudaContext& cu) : ApplyMonteCarloBarostatKernel(name, platform), cu(cu),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param barostat   the MonteCarloBarostat this kernel will be used for
     */
    void initialize(const System& system, const Force& barostat);
    /**
     * Attempt a Monte Carlo step, scaling particle positions (or cluster centers) by a specified value.
     * This version scales the x, y, and z positions independently.
     * This is called BEFORE the periodic box size is modified.  It should begin by translating each particle
     * or cluster into the first periodic box, so that coordinates will still be correct after the box size
     * is changed.
     *
     * @param context    the context in which to execute this kernel
     * @param scaleX     the scale factor by which to multiply particle x-coordinate
     * @param scaleY     the scale factor by which to multiply particle y-coordinate
     * @param scaleZ     the scale factor by which to multiply particle z-coordinate
     */
    void scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ);
    /**
     * Reject the most recent Monte Carlo step, restoring the particle positions to where they were before
     * scaleCoordinates() was last called.
     *
     * @param context    the context in which to execute this kernel
     */
    void restoreCoordinates(ContextImpl& context);
private:
    CudaContext& cu;
    bool hasInitializedKernels;
    int numMolecules;
    CudaArray savedPositions;
    CudaArray savedForces;
    CudaArray moleculeAtoms;
    CudaArray moleculeStartIndex;
    CUfunction kernel;
    std::vector<int> lastAtomOrder;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class CudaRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    CudaRemoveCMMotionKernel(std::string name, const Platform& platform, CudaContext& cu) : RemoveCMMotionKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel, setting up the particle masses.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CMMotionRemover this kernel will be used for
     */
    void initialize(const System& system, const CMMotionRemover& force);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     */
    void execute(ContextImpl& context);
private:
    CudaContext& cu;
    int frequency;
    CudaArray cmMomentum;
    CUfunction kernel1, kernel2;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAKERNELS_H_*/

