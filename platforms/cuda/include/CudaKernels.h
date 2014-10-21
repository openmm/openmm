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
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
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
#include "CudaParameterSet.h"
#include "CudaSort.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include <cufft.h>

namespace OpenMM {

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
     * @return the potential energy of the system.  This value is added to all values returned by ForceImpls'
     * calcForcesAndEnergy() methods.  That is, each force kernel may <i>either</i> return its contribution to the
     * energy directly, <i>or</i> add it to an internal buffer so that it will be included here.
     */
    double finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups);
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
    void setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) const;
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
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcHarmonicBondForceKernel();
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
    int numBonds;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    CudaCalcCustomBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomBondForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL), globals(NULL) {
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
    int numBonds;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaParameterSet* params;
    CudaArray* globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    CudaCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcHarmonicAngleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcHarmonicAngleForceKernel();
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
    int numAngles;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
public:
    CudaCalcCustomAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomAngleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL), globals(NULL) {
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
    int numAngles;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaParameterSet* params;
    CudaArray* globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    CudaCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcPeriodicTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }
    ~CudaCalcPeriodicTorsionForceKernel();
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
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaArray* params;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    CudaCalcRBTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcRBTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params1(NULL), params2(NULL) {
    }
    ~CudaCalcRBTorsionForceKernel();
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
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaArray* params1;
    CudaArray* params2;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
public:
    CudaCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCMAPTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), coefficients(NULL), mapPositions(NULL), torsionMaps(NULL) {
    }
    ~CudaCalcCMAPTorsionForceKernel();
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
private:
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaArray* coefficients;
    CudaArray* mapPositions;
    CudaArray* torsionMaps;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
public:
    CudaCalcCustomTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), system(system), params(NULL), globals(NULL) {
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
    int numTorsions;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaParameterSet* params;
    CudaArray* globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CudaCalcNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcNonbondedForceKernel(name, platform),
            cu(cu), hasInitializedFFT(false), sigmaEpsilon(NULL), exceptionParams(NULL), cosSinSums(NULL), directPmeGrid(NULL), reciprocalPmeGrid(NULL),
            pmeBsplineModuliX(NULL), pmeBsplineModuliY(NULL), pmeBsplineModuliZ(NULL),  pmeAtomRange(NULL), pmeAtomGridIndex(NULL), sort(NULL), pmeio(NULL) {
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
private:
    class SortTrait : public CudaSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "INT_MIN";}
        const char* getMaxKey() const {return "INT_MAX";}
        const char* getMaxValue() const {return "make_int2(INT_MAX, INT_MAX)";}
        const char* getSortKey() const {return "value.y";}
    };
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncStreamPreComputation;
    class SyncStreamPostComputation;
    CudaContext& cu;
    bool hasInitializedFFT;
    CudaArray* sigmaEpsilon;
    CudaArray* exceptionParams;
    CudaArray* cosSinSums;
    CudaArray* directPmeGrid;
    CudaArray* reciprocalPmeGrid;
    CudaArray* pmeBsplineModuliX;
    CudaArray* pmeBsplineModuliY;
    CudaArray* pmeBsplineModuliZ;
    CudaArray* pmeAtomRange;
    CudaArray* pmeAtomGridIndex;
    CudaSort* sort;
    Kernel cpuPme;
    PmeIO* pmeio;
    CUstream pmeStream;
    CUevent pmeSyncEvent;
    cufftHandle fftForward;
    cufftHandle fftBackward;
    CUfunction ewaldSumsKernel;
    CUfunction ewaldForcesKernel;
    CUfunction pmeGridIndexKernel;
    CUfunction pmeSpreadChargeKernel;
    CUfunction pmeFinishSpreadChargeKernel;
    CUfunction pmeEvalEnergyKernel;
    CUfunction pmeConvolutionKernel;
    CUfunction pmeInterpolateForceKernel;
    std::map<std::string, std::string> pmeDefines;
    std::vector<std::pair<int, int> > exceptionAtoms;
    double ewaldSelfEnergy, dispersionCoefficient, alpha;
    int interpolateForceThreads;
    bool hasCoulomb, hasLJ, usePmeStream;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    CudaCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomNonbondedForceKernel(name, platform),
            cu(cu), params(NULL), globals(NULL), interactionGroupData(NULL), forceCopy(NULL), system(system), hasInitializedKernel(false) {
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
    void initInteractionGroups(const CustomNonbondedForce& force, const std::string& interactionSource);
    CudaContext& cu;
    CudaParameterSet* params;
    CudaArray* globals;
    CudaArray* interactionGroupData;
    CUfunction interactionGroupKernel;
    std::vector<void*> interactionGroupArgs;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray*> tabulatedFunctions;
    double longRangeCoefficient;
    bool hasInitializedLongRangeCorrection, hasInitializedKernel;
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
            hasCreatedKernels(false), params(NULL), bornSum(NULL), bornRadii(NULL), bornForce(NULL), obcChain(NULL) {
    }
    ~CudaCalcGBSAOBCForceKernel();
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
    double prefactor, surfaceAreaFactor;
    bool hasCreatedKernels;
    int maxTiles;
    CudaContext& cu;
    CudaArray* params;
    CudaArray* bornSum;
    CudaArray* bornRadii;
    CudaArray* bornForce;
    CudaArray* obcChain;
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
            hasInitializedKernels(false), cu(cu), params(NULL), computedValues(NULL), energyDerivs(NULL), energyDerivChain(NULL), longEnergyDerivs(NULL), globals(NULL),
            valueBuffers(NULL), system(system) {
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
    bool hasInitializedKernels, needParameterGradient;
    int maxTiles, numComputedValues;
    CudaContext& cu;
    CudaParameterSet* params;
    CudaParameterSet* computedValues;
    CudaParameterSet* energyDerivs;
    CudaParameterSet* energyDerivChain;
    CudaArray* longEnergyDerivs;
    CudaArray* globals;
    CudaArray* valueBuffers;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray*> tabulatedFunctions;
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
            hasInitializedKernel(false), cu(cu), system(system), params(NULL), globals(NULL) {
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
    int numParticles;
    bool hasInitializedKernel;
    CudaContext& cu;
    const System& system;
    CudaParameterSet* params;
    CudaArray* globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
};

/**
 * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system.
 */
class CudaCalcCustomHbondForceKernel : public CalcCustomHbondForceKernel {
public:
    CudaCalcCustomHbondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomHbondForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), donorParams(NULL), acceptorParams(NULL), donors(NULL), acceptors(NULL),
            globals(NULL), donorExclusions(NULL), acceptorExclusions(NULL), system(system) {
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
    int numDonors, numAcceptors;
    bool hasInitializedKernel;
    CudaContext& cu;
    CudaParameterSet* donorParams;
    CudaParameterSet* acceptorParams;
    CudaArray* globals;
    CudaArray* donors;
    CudaArray* acceptors;
    CudaArray* donorExclusions;
    CudaArray* acceptorExclusions;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray*> tabulatedFunctions;
    std::vector<void*> donorArgs, acceptorArgs;
    const System& system;
    CUfunction donorKernel, acceptorKernel;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
 */
class CudaCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
public:
    CudaCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomCompoundBondForceKernel(name, platform),
            cu(cu), params(NULL), globals(NULL), system(system) {
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
    int numBonds;
    CudaContext& cu;
    CudaParameterSet* params;
    CudaArray* globals;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray*> tabulatedFunctions;
    const System& system;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system.
 */
class CudaCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    CudaCalcCustomManyParticleForceKernel(std::string name, const Platform& platform, CudaContext& cu, const System& system) : CalcCustomManyParticleForceKernel(name, platform),
            hasInitializedKernel(false), cu(cu), params(NULL), particleTypes(NULL), orderIndex(NULL), particleOrder(NULL), exclusions(NULL),
            exclusionStartIndex(NULL), blockCenter(NULL), blockBoundingBox(NULL), neighborPairs(NULL), numNeighborPairs(NULL), neighborStartIndex(NULL),
            numNeighborsForAtom(NULL), neighbors(NULL), system(system) {
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
    CudaContext& cu;
    bool hasInitializedKernel;
    NonbondedMethod nonbondedMethod;
    int maxNeighborPairs, forceWorkgroupSize, findNeighborsWorkgroupSize;
    CudaParameterSet* params;
    CudaArray* particleTypes;
    CudaArray* orderIndex;
    CudaArray* particleOrder;
    CudaArray* exclusions;
    CudaArray* exclusionStartIndex;
    CudaArray* blockCenter;
    CudaArray* blockBoundingBox;
    CudaArray* neighborPairs;
    CudaArray* numNeighborPairs;
    CudaArray* neighborStartIndex;
    CudaArray* numNeighborsForAtom;
    CudaArray* neighbors;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<CudaArray*> tabulatedFunctions;
    std::vector<void*> forceArgs, blockBoundsArgs, neighborsArgs, startIndicesArgs, copyPairsArgs;
    const System& system;
    CUfunction forceKernel, blockBoundsKernel, neighborsKernel, startIndicesKernel, copyPairsKernel;
    CUdeviceptr globalsPtr;
    CUevent event;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class CudaIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    CudaIntegrateVerletStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateVerletStepKernel(name, platform), cu(cu) {
    }
    ~CudaIntegrateVerletStepKernel();
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
    double prevStepSize;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CudaIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    CudaIntegrateLangevinStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateLangevinStepKernel(name, platform), cu(cu), params(NULL) {
    }
    ~CudaIntegrateLangevinStepKernel();
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
    CudaArray* params;
    CUfunction kernel1, kernel2;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class CudaIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    CudaIntegrateBrownianStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateBrownianStepKernel(name, platform), cu(cu) {
    }
    ~CudaIntegrateBrownianStepKernel();
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
    ~CudaIntegrateVariableVerletStepKernel();
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
    CudaIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateVariableLangevinStepKernel(name, platform),
            cu(cu), params(NULL) {
    }
    ~CudaIntegrateVariableLangevinStepKernel();
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
    CudaArray* params;
    CUfunction kernel1, kernel2, selectSizeKernel;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by CustomIntegrator to take one time step.
 */
class CudaIntegrateCustomStepKernel : public IntegrateCustomStepKernel {
public:
    CudaIntegrateCustomStepKernel(std::string name, const Platform& platform, CudaContext& cu) : IntegrateCustomStepKernel(name, platform), cu(cu),
            hasInitializedKernels(false), localValuesAreCurrent(false), globalValues(NULL), contextParameterValues(NULL), sumBuffer(NULL), potentialEnergy(NULL),
            kineticEnergy(NULL), uniformRandoms(NULL), randomSeed(NULL), perDofValues(NULL) {
    }
    ~CudaIntegrateCustomStepKernel();
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
    std::string createGlobalComputation(const std::string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator, const std::string& energyName);
    std::string createPerDofComputation(const std::string& variable, const Lepton::ParsedExpression& expr, int component, CustomIntegrator& integrator, const std::string& forceName, const std::string& energyName);
    void prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid);
    void recordChangedParameters(ContextImpl& context);
    CudaContext& cu;
    double prevStepSize;
    int numGlobalVariables;
    bool hasInitializedKernels, deviceValuesAreCurrent, modifiesParameters, keNeedsForce;
    mutable bool localValuesAreCurrent;
    CudaArray* globalValues;
    CudaArray* contextParameterValues;
    CudaArray* sumBuffer;
    CudaArray* potentialEnergy;
    CudaArray* kineticEnergy;
    CudaArray* uniformRandoms;
    CudaArray* randomSeed;
    std::map<int, CudaArray*> savedForces;
    std::set<int> validSavedForces;
    CudaParameterSet* perDofValues;
    mutable std::vector<std::vector<float> > localPerDofValuesFloat;
    mutable std::vector<std::vector<double> > localPerDofValuesDouble;
    std::vector<float> contextValuesFloat;
    std::vector<double> contextValuesDouble;
    std::vector<std::vector<CUfunction> > kernels;
    std::vector<std::vector<std::vector<void*> > > kernelArgs;
    std::vector<void*> kineticEnergyArgs;
    CUfunction sumPotentialEnergyKernel, randomKernel, kineticEnergyKernel, sumKineticEnergyKernel;
    std::vector<CustomIntegrator::ComputationType> stepType;
    std::vector<bool> needsForces;
    std::vector<bool> needsEnergy;
    std::vector<bool> invalidatesForces;
    std::vector<bool> merged;
    std::vector<int> forceGroup;
    std::vector<int> requiredGaussian;
    std::vector<int> requiredUniform;
    std::vector<std::string> parameterNames;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class CudaApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    CudaApplyAndersenThermostatKernel(std::string name, const Platform& platform, CudaContext& cu) : ApplyAndersenThermostatKernel(name, platform), cu(cu),
            atomGroups(NULL) {
    }
    ~CudaApplyAndersenThermostatKernel();
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
    CudaArray* atomGroups;
    CUfunction kernel;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class CudaApplyMonteCarloBarostatKernel : public ApplyMonteCarloBarostatKernel {
public:
    CudaApplyMonteCarloBarostatKernel(std::string name, const Platform& platform, CudaContext& cu) : ApplyMonteCarloBarostatKernel(name, platform), cu(cu),
            hasInitializedKernels(false), savedPositions(NULL), moleculeAtoms(NULL), moleculeStartIndex(NULL) {
    }
    ~CudaApplyMonteCarloBarostatKernel();
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
    CudaArray* savedPositions;
    CudaArray* moleculeAtoms;
    CudaArray* moleculeStartIndex;
    CUfunction kernel;
    std::vector<int> lastAtomOrder;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class CudaRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    CudaRemoveCMMotionKernel(std::string name, const Platform& platform, CudaContext& cu) : RemoveCMMotionKernel(name, platform), cu(cu), cmMomentum(NULL) {
    }
    ~CudaRemoveCMMotionKernel();
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
    CudaArray* cmMomentum;
    CUfunction kernel1, kernel2;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAKERNELS_H_*/

