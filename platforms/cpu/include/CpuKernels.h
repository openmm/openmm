#ifndef OPENMM_CPUKERNELS_H_
#define OPENMM_CPUKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CpuBondForce.h"
#include "CpuCustomGBForce.h"
#include "CpuCustomManyParticleForce.h"
#include "CpuCustomNonbondedForce.h"
#include "CpuGBSAOBCForce.h"
#include "CpuLangevinDynamics.h"
#include "CpuNeighborList.h"
#include "CpuNonbondedForce.h"
#include "CpuPlatform.h"
#include "openmm/kernels.h"
#include "openmm/System.h"

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class CpuCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    class InitForceTask;
    class SumForceTask;
    CpuCalcForcesAndEnergyKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data, ContextImpl& context);
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
    CpuPlatform::PlatformData& data;
    Kernel referenceKernel;
    std::vector<RealVec> lastPositions;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CpuCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    CpuCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcHarmonicAngleForceKernel(name, platform), data(data), angleIndexArray(NULL), angleParamArray(NULL), usePeriodic(false) {
    }
    ~CpuCalcHarmonicAngleForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numAngles;
    int **angleIndexArray;
    RealOpenMM **angleParamArray;
    CpuBondForce bondForce;
    bool usePeriodic;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CpuCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    CpuCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcPeriodicTorsionForceKernel(name, platform), data(data), torsionIndexArray(NULL), torsionParamArray(NULL), usePeriodic(false) {
    }
    ~CpuCalcPeriodicTorsionForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numTorsions;
    int **torsionIndexArray;
    RealOpenMM **torsionParamArray;
    CpuBondForce bondForce;
    bool usePeriodic;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CpuCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    CpuCalcRBTorsionForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcRBTorsionForceKernel(name, platform), data(data), torsionIndexArray(NULL), torsionParamArray(NULL), usePeriodic(false) {
    }
    ~CpuCalcRBTorsionForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numTorsions;
    int **torsionIndexArray;
    RealOpenMM **torsionParamArray;
    CpuBondForce bondForce;
    bool usePeriodic;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CpuCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CpuCalcNonbondedForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data);
    ~CpuCalcNonbondedForceKernel();
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
private:
    class PmeIO;
    CpuPlatform::PlatformData& data;
    int numParticles, num14;
    int **bonded14IndexArray;
    double **bonded14ParamArray;
    double nonbondedCutoff, switchingDistance, rfDielectric, ewaldAlpha, ewaldSelfEnergy, dispersionCoefficient;
    int kmax[3], gridSize[3];
    bool useSwitchingFunction, useOptimizedPme, hasInitializedPme;
    std::vector<std::set<int> > exclusions;
    std::vector<std::pair<float, float> > particleParams;
    NonbondedMethod nonbondedMethod;
    CpuNonbondedForce* nonbonded;
    Kernel optimizedPme;
    CpuBondForce bondForce;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class CpuCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    CpuCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data);
    ~CpuCalcCustomNonbondedForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numParticles;
    double **particleParamArray;
    double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
    bool useSwitchingFunction, hasInitializedLongRangeCorrection;
    CustomNonbondedForce* forceCopy;
    std::map<std::string, double> globalParamValues;
    std::vector<std::set<int> > exclusions;
    std::vector<std::string> parameterNames, globalParameterNames;
    std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;
    NonbondedMethod nonbondedMethod;
    CpuCustomNonbondedForce* nonbonded;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CpuCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    CpuCalcGBSAOBCForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) : CalcGBSAOBCForceKernel(name, platform),
            data(data) {
    }
    ~CpuCalcGBSAOBCForceKernel();
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
    CpuPlatform::PlatformData& data;
    std::vector<std::pair<float, float> > particleParams;
    CpuGBSAOBCForce obc;
};

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system.
 */
class CpuCalcCustomGBForceKernel : public CalcCustomGBForceKernel {
public:
    CpuCalcCustomGBForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcCustomGBForceKernel(name, platform), data(data) {
    }
    ~CpuCalcCustomGBForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numParticles;
    bool isPeriodic;
    RealOpenMM **particleParamArray;
    RealOpenMM nonbondedCutoff;
    CpuCustomGBForce* ixn;
    std::vector<std::set<int> > exclusions;
    std::vector<std::string> particleParameterNames, globalParameterNames, valueNames;
    std::vector<OpenMM::CustomGBForce::ComputationType> valueTypes;
    std::vector<OpenMM::CustomGBForce::ComputationType> energyTypes;
    NonbondedMethod nonbondedMethod;
};

/**
 * This kernel is invoked by CustomManyParticleForce to calculate the forces acting on the system and the energy of the system.
 */
class CpuCalcCustomManyParticleForceKernel : public CalcCustomManyParticleForceKernel {
public:
    CpuCalcCustomManyParticleForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) : CalcCustomManyParticleForceKernel(name, platform),
            data(data), ixn(NULL) {
    }
    ~CpuCalcCustomManyParticleForceKernel();
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
    CpuPlatform::PlatformData& data;
    int numParticles;
    RealOpenMM cutoffDistance;
    RealOpenMM **particleParamArray;
    CpuCustomManyParticleForce* ixn;
    std::vector<std::string> globalParameterNames;
    NonbondedMethod nonbondedMethod;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CpuIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    CpuIntegrateLangevinStepKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) : IntegrateLangevinStepKernel(name, platform),
            data(data), dynamics(NULL) {
    }
    ~CpuIntegrateLangevinStepKernel();
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
    CpuPlatform::PlatformData& data;
    CpuLangevinDynamics* dynamics;
    std::vector<RealOpenMM> masses;
    double prevTemp, prevFriction, prevStepSize;
};

} // namespace OpenMM

#endif /*OPENMM_CPUKERNELS_H_*/

