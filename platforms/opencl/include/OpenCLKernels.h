#ifndef OPENMM_OPENCLKERNELS_H_
#define OPENMM_OPENCLKERNELS_H_

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

#include "OpenCLPlatform.h"
#include "OpenCLArray.h"
#include "OpenCLContext.h"
#include "OpenCLFFT3D.h"
#include "OpenCLParameterSet.h"
#include "OpenCLSort.h"
#include "openmm/kernels.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/CustomIntegratorUtilities.h"
#include "lepton/CompiledExpression.h"
#include "lepton/ExpressionProgram.h"
#include "openmm/System.h"

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class OpenCLCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    OpenCLCalcForcesAndEnergyKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcForcesAndEnergyKernel(name, platform), cl(cl) {
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
   OpenCLContext& cl;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class OpenCLUpdateStateDataKernel : public UpdateStateDataKernel {
public:
    OpenCLUpdateStateDataKernel(std::string name, const Platform& platform, OpenCLContext& cl) : UpdateStateDataKernel(name, platform), cl(cl) {
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
    OpenCLContext& cl;
};

/**
 * This kernel modifies the positions of particles to enforce distance constraints.
 */
class OpenCLApplyConstraintsKernel : public ApplyConstraintsKernel {
public:
    OpenCLApplyConstraintsKernel(std::string name, const Platform& platform, OpenCLContext& cl) : ApplyConstraintsKernel(name, platform),
            cl(cl), hasInitializedKernel(false) {
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
    OpenCLContext& cl;
    bool hasInitializedKernel;
    cl::Kernel applyDeltasKernel;
};

/**
 * This kernel recomputes the positions of virtual sites.
 */
class OpenCLVirtualSitesKernel : public VirtualSitesKernel {
public:
    OpenCLVirtualSitesKernel(std::string name, const Platform& platform, OpenCLContext& cl) : VirtualSitesKernel(name, platform), cl(cl) {
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
    OpenCLContext& cl;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    OpenCLCalcHarmonicBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcHarmonicBondForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLArray params;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    OpenCLCalcCustomBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomBondForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL) {
    }
    ~OpenCLCalcCustomBondForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    OpenCLCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcHarmonicAngleForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLArray params;
};

/**
 * This kernel is invoked by CustomAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomAngleForceKernel : public CalcCustomAngleForceKernel {
public:
    OpenCLCalcCustomAngleForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomAngleForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL) {
    }
    ~OpenCLCalcCustomAngleForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    OpenCLCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcPeriodicTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLArray params;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    OpenCLCalcRBTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcRBTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLArray params;
};

/**
 * This kernel is invoked by CMAPTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCMAPTorsionForceKernel : public CalcCMAPTorsionForceKernel {
public:
    OpenCLCalcCMAPTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCMAPTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    std::vector<mm_int2> mapPositionsVec;
    OpenCLArray coefficients;
    OpenCLArray mapPositions;
    OpenCLArray torsionMaps;
};

/**
 * This kernel is invoked by CustomTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomTorsionForceKernel : public CalcCustomTorsionForceKernel {
public:
    OpenCLCalcCustomTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL) {
    }
    ~OpenCLCalcCustomTorsionForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class OpenCLCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    OpenCLCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcNonbondedForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), sort(NULL), fft(NULL), dispersionFft(NULL), pmeio(NULL), usePmeQueue(false) {
    }
    ~OpenCLCalcNonbondedForceKernel();
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
     * Get the parameters being used for the dispersion term in LJPME.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
private:
    class SortTrait : public OpenCLSort::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "INT_MIN";}
        const char* getMaxKey() const {return "INT_MAX";}
        const char* getMaxValue() const {return "(int2) (INT_MAX, INT_MAX)";}
        const char* getSortKey() const {return "value.y";}
    };
    class ForceInfo;
    class PmeIO;
    class PmePreComputation;
    class PmePostComputation;
    class SyncQueuePreComputation;
    class SyncQueuePostComputation;
    OpenCLContext& cl;
    ForceInfo* info;
    bool hasInitializedKernel;
    OpenCLArray charges;
    OpenCLArray sigmaEpsilon;
    OpenCLArray exceptionParams;
    OpenCLArray exclusionAtoms;
    OpenCLArray exclusionParams;
    OpenCLArray baseParticleParams;
    OpenCLArray baseExceptionParams;
    OpenCLArray particleParamOffsets;
    OpenCLArray exceptionParamOffsets;
    OpenCLArray particleOffsetIndices;
    OpenCLArray exceptionOffsetIndices;
    OpenCLArray globalParams;
    OpenCLArray cosSinSums;
    OpenCLArray pmeGrid1;
    OpenCLArray pmeGrid2;
    OpenCLArray pmeBsplineModuliX;
    OpenCLArray pmeBsplineModuliY;
    OpenCLArray pmeBsplineModuliZ;
    OpenCLArray pmeDispersionBsplineModuliX;
    OpenCLArray pmeDispersionBsplineModuliY;
    OpenCLArray pmeDispersionBsplineModuliZ;
    OpenCLArray pmeBsplineTheta;
    OpenCLArray pmeAtomRange;
    OpenCLArray pmeAtomGridIndex;
    OpenCLArray pmeEnergyBuffer;
    OpenCLSort* sort;
    cl::CommandQueue pmeQueue;
    cl::Event pmeSyncEvent;
    OpenCLFFT3D* fft;
    OpenCLFFT3D* dispersionFft;
    Kernel cpuPme;
    PmeIO* pmeio;
    SyncQueuePostComputation* syncQueue;
    cl::Kernel computeParamsKernel, computeExclusionParamsKernel;
    cl::Kernel ewaldSumsKernel;
    cl::Kernel ewaldForcesKernel;
    cl::Kernel pmeAtomRangeKernel;
    cl::Kernel pmeDispersionAtomRangeKernel;
    cl::Kernel pmeZIndexKernel;
    cl::Kernel pmeDispersionZIndexKernel;
    cl::Kernel pmeUpdateBsplinesKernel;
    cl::Kernel pmeDispersionUpdateBsplinesKernel;
    cl::Kernel pmeSpreadChargeKernel;
    cl::Kernel pmeDispersionSpreadChargeKernel;
    cl::Kernel pmeFinishSpreadChargeKernel;
    cl::Kernel pmeDispersionFinishSpreadChargeKernel;
    cl::Kernel pmeConvolutionKernel;
    cl::Kernel pmeDispersionConvolutionKernel;
    cl::Kernel pmeEvalEnergyKernel;
    cl::Kernel pmeDispersionEvalEnergyKernel;
    cl::Kernel pmeInterpolateForceKernel;
    cl::Kernel pmeDispersionInterpolateForceKernel;
    std::map<std::string, std::string> pmeDefines;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;
    double ewaldSelfEnergy, dispersionCoefficient, alpha, dispersionAlpha;
    int gridSizeX, gridSizeY, gridSizeZ;
    int dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ;
    bool hasCoulomb, hasLJ, usePmeQueue, doLJPME, usePosqCharges, recomputeParams, hasOffsets;
    NonbondedMethod nonbondedMethod;
    static const int PmeOrder = 5;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    OpenCLCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomNonbondedForceKernel(name, platform),
            cl(cl), params(NULL), forceCopy(NULL), system(system), hasInitializedKernel(false) {
    }
    ~OpenCLCalcCustomNonbondedForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    OpenCLArray interactionGroupData, filteredGroupData, numGroupTiles;
    cl::Kernel interactionGroupKernel, prepareNeighborListKernel, buildNeighborListKernel;
    std::vector<void*> interactionGroupArgs;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
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
class OpenCLCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    OpenCLCalcGBSAOBCForceKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcGBSAOBCForceKernel(name, platform), cl(cl),
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLArray params;
    OpenCLArray charges;
    OpenCLArray bornSum;
    OpenCLArray longBornSum;
    OpenCLArray bornRadii;
    OpenCLArray bornForce;
    OpenCLArray longBornForce;
    OpenCLArray obcChain;
    cl::Kernel computeBornSumKernel;
    cl::Kernel reduceBornSumKernel;
    cl::Kernel force1Kernel;
    cl::Kernel reduceBornForceKernel;
};

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomGBForceKernel : public CalcCustomGBForceKernel {
public:
    OpenCLCalcCustomGBForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomGBForceKernel(name, platform),
            hasInitializedKernels(false), cl(cl), params(NULL), computedValues(NULL), energyDerivs(NULL), energyDerivChain(NULL), system(system) {
    }
    ~OpenCLCalcCustomGBForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLParameterSet* params;
    OpenCLParameterSet* computedValues;
    OpenCLParameterSet* energyDerivs;
    OpenCLParameterSet* energyDerivChain;
    std::vector<OpenCLParameterSet*> dValuedParam;
    std::vector<OpenCLArray> dValue0dParam;
    OpenCLArray longEnergyDerivs;
    OpenCLArray globals;
    OpenCLArray valueBuffers;
    OpenCLArray longValueBuffers;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
    std::vector<bool> pairValueUsesParam, pairEnergyUsesParam, pairEnergyUsesValue;
    const System& system;
    cl::Kernel pairValueKernel, perParticleValueKernel, pairEnergyKernel, perParticleEnergyKernel, gradientChainRuleKernel;
    std::string pairValueSrc, pairEnergySrc;
    std::map<std::string, std::string> pairValueDefines, pairEnergyDefines;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
public:
    OpenCLCalcCustomExternalForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomExternalForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL) {
    }
    ~OpenCLCalcCustomExternalForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    const System& system;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
};

/**
 * This kernel is invoked by CustomHbondForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomHbondForceKernel : public CalcCustomHbondForceKernel {
public:
    OpenCLCalcCustomHbondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomHbondForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), donorParams(NULL), acceptorParams(NULL), system(system) {
    }
    ~OpenCLCalcCustomHbondForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLParameterSet* donorParams;
    OpenCLParameterSet* acceptorParams;
    OpenCLArray globals;
    OpenCLArray donors;
    OpenCLArray acceptors;
    OpenCLArray donorBufferIndices;
    OpenCLArray acceptorBufferIndices;
    OpenCLArray donorExclusions;
    OpenCLArray acceptorExclusions;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
    const System& system;
    cl::Kernel donorKernel, acceptorKernel;
};

/**
 * This kernel is invoked by CustomCentroidBondForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomCentroidBondForceKernel : public CalcCustomCentroidBondForceKernel {
public:
    OpenCLCalcCustomCentroidBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomCentroidBondForceKernel(name, platform),
            cl(cl), params(NULL), system(system) {
    }
    ~OpenCLCalcCustomCentroidBondForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    OpenCLArray groupParticles;
    OpenCLArray groupWeights;
    OpenCLArray groupOffsets;
    OpenCLArray groupForces;
    OpenCLArray bondGroups;
    OpenCLArray centerPositions;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
    cl::Kernel computeCentersKernel, groupForcesKernel, applyForcesKernel;
    const System& system;
};

/**
 * This kernel is invoked by CustomCompoundBondForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomCompoundBondForceKernel : public CalcCustomCompoundBondForceKernel {
public:
    OpenCLCalcCustomCompoundBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomCompoundBondForceKernel(name, platform),
            cl(cl), params(NULL), system(system) {
    }
    ~OpenCLCalcCustomCompoundBondForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
    const System& system;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    OpenCLCalcCustomManyParticleForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, const System& system) : CalcCustomManyParticleForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), params(NULL), system(system) {
    }
    ~OpenCLCalcCustomManyParticleForceKernel();
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
    OpenCLContext& cl;
    ForceInfo* info;
    bool hasInitializedKernel;
    NonbondedMethod nonbondedMethod;
    int maxNeighborPairs, forceWorkgroupSize, findNeighborsWorkgroupSize;
    OpenCLParameterSet* params;
    OpenCLArray globals;
    OpenCLArray particleTypes;
    OpenCLArray orderIndex;
    OpenCLArray particleOrder;
    OpenCLArray exclusions;
    OpenCLArray exclusionStartIndex;
    OpenCLArray blockCenter;
    OpenCLArray blockBoundingBox;
    OpenCLArray neighborPairs;
    OpenCLArray numNeighborPairs;
    OpenCLArray neighborStartIndex;
    OpenCLArray numNeighborsForAtom;
    OpenCLArray neighbors;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<OpenCLArray> tabulatedFunctions;
    const System& system;
    cl::Kernel forceKernel, blockBoundsKernel, neighborsKernel, startIndicesKernel, copyPairsKernel;
};

/**
 * This kernel is invoked by GayBerneForce to calculate the forces acting on the system.
 */
class OpenCLCalcGayBerneForceKernel : public CalcGayBerneForceKernel {
public:
    OpenCLCalcGayBerneForceKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcGayBerneForceKernel(name, platform), cl(cl),
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
    OpenCLContext& cl;
    ForceInfo* info;
    bool hasInitializedKernels;
    int numRealParticles, maxNeighborBlocks;
    GayBerneForce::NonbondedMethod nonbondedMethod;
    OpenCLArray sortedParticles;
    OpenCLArray axisParticleIndices;
    OpenCLArray sigParams;
    OpenCLArray epsParams;
    OpenCLArray scale;
    OpenCLArray exceptionParticles;
    OpenCLArray exceptionParams;
    OpenCLArray aMatrix;
    OpenCLArray bMatrix;
    OpenCLArray gMatrix;
    OpenCLArray exclusions;
    OpenCLArray exclusionStartIndex;
    OpenCLArray blockCenter;
    OpenCLArray blockBoundingBox;
    OpenCLArray neighbors;
    OpenCLArray neighborIndex;
    OpenCLArray neighborBlockCount;
    OpenCLArray sortedPos;
    OpenCLArray torque;
    std::vector<bool> isRealParticle;
    std::vector<std::pair<int, int> > exceptionAtoms;
    std::vector<std::pair<int, int> > excludedPairs;
    cl::Kernel framesKernel, blockBoundsKernel, neighborsKernel, forceKernel, torqueKernel;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomCVForceKernel : public CalcCustomCVForceKernel {
public:
    OpenCLCalcCustomCVForceKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcCustomCVForceKernel(name, platform),
            cl(cl), hasInitializedKernels(false) {
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    Lepton::ExpressionProgram energyExpression;
    std::vector<std::string> variableNames, paramDerivNames, globalParameterNames;
    std::vector<Lepton::ExpressionProgram> variableDerivExpressions;
    std::vector<Lepton::ExpressionProgram> paramDerivExpressions;
    std::vector<OpenCLArray> cvForces;
    OpenCLArray invAtomOrder;
    OpenCLArray innerInvAtomOrder;
    cl::Kernel copyStateKernel, copyForcesKernel, addForcesKernel;
};

/**
 * This kernel is invoked by RMSDForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcRMSDForceKernel : public CalcRMSDForceKernel {
public:
    OpenCLCalcRMSDForceKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcRMSDForceKernel(name, platform), cl(cl) {
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
    OpenCLContext& cl;
    ForceInfo* info;
    double sumNormRef;
    OpenCLArray referencePos;
    OpenCLArray particles;
    OpenCLArray buffer;
    cl::Kernel kernel1, kernel2;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class OpenCLIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    OpenCLIntegrateVerletStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateVerletStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
    }
    ~OpenCLIntegrateVerletStepKernel();
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    cl::Kernel kernel1, kernel2;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class OpenCLIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    OpenCLIntegrateLangevinStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateLangevinStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
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
    OpenCLContext& cl;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    OpenCLArray params;
    cl::Kernel kernel1, kernel2;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class OpenCLIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    OpenCLIntegrateBrownianStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateBrownianStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false), prevTemp(-1), prevFriction(-1), prevStepSize(-1) {
    }
    ~OpenCLIntegrateBrownianStepKernel();
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
    OpenCLContext& cl;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    cl::Kernel kernel1, kernel2;
};

/**
 * This kernel is invoked by VariableVerletIntegrator to take one time step.
 */
class OpenCLIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
public:
    OpenCLIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateVariableVerletStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
    }
    ~OpenCLIntegrateVariableVerletStepKernel();
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    int blockSize;
    cl::Kernel kernel1, kernel2, selectSizeKernel;
};

/**
 * This kernel is invoked by VariableLangevinIntegrator to take one time step.
 */
class OpenCLIntegrateVariableLangevinStepKernel : public IntegrateVariableLangevinStepKernel {
public:
    OpenCLIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateVariableLangevinStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    int blockSize;
    OpenCLArray params;
    cl::Kernel kernel1, kernel2, selectSizeKernel;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by CustomIntegrator to take one time step.
 */
class OpenCLIntegrateCustomStepKernel : public IntegrateCustomStepKernel {
public:
    enum GlobalTargetType {DT, VARIABLE, PARAMETER};
    OpenCLIntegrateCustomStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateCustomStepKernel(name, platform), cl(cl),
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
    OpenCLContext& cl;
    double energy;
    float energyFloat;
    int numGlobalVariables, sumWorkGroupSize;
    bool hasInitializedKernels, deviceGlobalsAreCurrent, modifiesParameters, keNeedsForce, hasAnyConstraints, needsEnergyParamDerivs;
    std::vector<bool> deviceValuesAreCurrent;
    mutable std::vector<bool> localValuesAreCurrent;
    OpenCLArray globalValues;
    OpenCLArray sumBuffer;
    OpenCLArray summedValue;
    OpenCLArray uniformRandoms;
    OpenCLArray randomSeed;
    OpenCLArray perDofEnergyParamDerivs;
    std::vector<OpenCLArray> tabulatedFunctions, perDofValues;
    std::map<int, double> savedEnergy;
    std::map<int, OpenCLArray> savedForces;
    std::set<int> validSavedForces;
    mutable std::vector<std::vector<mm_float4> > localPerDofValuesFloat;
    mutable std::vector<std::vector<mm_double4> > localPerDofValuesDouble;
    std::map<std::string, double> energyParamDerivs;
    std::vector<std::string> perDofEnergyParamDerivNames;
    std::vector<cl_double> localPerDofEnergyParamDerivs;
    std::vector<double> localGlobalValues;
    std::vector<double> initialGlobalVariables;
    std::vector<std::vector<cl::Kernel> > kernels;
    cl::Kernel randomKernel, kineticEnergyKernel, sumKineticEnergyKernel;
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

class OpenCLIntegrateCustomStepKernel::GlobalTarget {
public:
    OpenCLIntegrateCustomStepKernel::GlobalTargetType type;
    int variableIndex;
    GlobalTarget() {
    }
    GlobalTarget(OpenCLIntegrateCustomStepKernel::GlobalTargetType type, int variableIndex) : type(type), variableIndex(variableIndex) {
    }
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class OpenCLApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    OpenCLApplyAndersenThermostatKernel(std::string name, const Platform& platform, OpenCLContext& cl) : ApplyAndersenThermostatKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    int randomSeed;
    OpenCLArray atomGroups;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by MonteCarloBarostat to adjust the periodic box volume
 */
class OpenCLApplyMonteCarloBarostatKernel : public ApplyMonteCarloBarostatKernel {
public:
    OpenCLApplyMonteCarloBarostatKernel(std::string name, const Platform& platform, OpenCLContext& cl) : ApplyMonteCarloBarostatKernel(name, platform), cl(cl),
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
    OpenCLContext& cl;
    bool hasInitializedKernels;
    int numMolecules;
    OpenCLArray savedPositions;
    OpenCLArray savedForces;
    OpenCLArray moleculeAtoms;
    OpenCLArray moleculeStartIndex;
    cl::Kernel kernel;
    std::vector<int> lastAtomOrder;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class OpenCLRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    OpenCLRemoveCMMotionKernel(std::string name, const Platform& platform, OpenCLContext& cl) : RemoveCMMotionKernel(name, platform), cl(cl) {
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
    OpenCLContext& cl;
    int frequency;
    OpenCLArray cmMomentum;
    cl::Kernel kernel1, kernel2;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLKERNELS_H_*/
