#ifndef OPENMM_REFERENCEKERNELS_H_
#define OPENMM_REFERENCEKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "kernels.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include "SimTKReference/ReferenceNeighborList.h"

class CpuObc;
class ReferenceAndersenThermostat;
class ReferenceBrownianDynamics;
class ReferenceStochasticDynamics;
class ReferenceShakeAlgorithm;
class ReferenceVerletDynamics;

namespace OpenMM {

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    ReferenceCalcHarmonicBondForceKernel(std::string name, const Platform& platform) : CalcHarmonicBondForceKernel(name, platform) {
    }
    ~ReferenceCalcHarmonicBondForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicBondForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicBondForce& force);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicBondForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    int numBonds;
    int **bondIndexArray;
    RealOpenMM **bondParamArray;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    ReferenceCalcHarmonicAngleForceKernel(std::string name, const Platform& platform) : CalcHarmonicAngleForceKernel(name, platform) {
    }
    ~ReferenceCalcHarmonicAngleForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    void initialize(const System& system, const HarmonicAngleForce& force);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicAngleForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    int numAngles;
    int **angleIndexArray;
    RealOpenMM **angleParamArray;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    ReferenceCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform) : CalcPeriodicTorsionForceKernel(name, platform) {
    }
    ~ReferenceCalcPeriodicTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the PeriodicTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const PeriodicTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the PeriodicTorsionForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    int numTorsions;
    int **torsionIndexArray;
    RealOpenMM **torsionParamArray;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    ReferenceCalcRBTorsionForceKernel(std::string name, const Platform& platform) : CalcRBTorsionForceKernel(name, platform) {
    }
    ~ReferenceCalcRBTorsionForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the RBTorsionForce this kernel will be used for
     */
    void initialize(const System& system, const RBTorsionForce& force);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the RBTorsionForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    int numTorsions;
    int **torsionIndexArray;
    RealOpenMM **torsionParamArray;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class ReferenceCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    ReferenceCalcNonbondedForceKernel(std::string name, const Platform& platform) : CalcNonbondedForceKernel(name, platform) {
    }
    ~ReferenceCalcNonbondedForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     * @param exclusions the i'th element lists the indices of all particles with which the i'th particle should not interact through
     *                   nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
     *                   the standard nonbonded calculation.
     */
    void initialize(const System& system, const NonbondedForce& force, const std::vector<std::set<int> >& exclusions);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the NonbondedForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    int numParticles, num14;
    int **exclusionArray, **bonded14IndexArray;
    RealOpenMM **particleParamArray, **bonded14ParamArray;
    RealOpenMM nonbondedCutoff, periodicBoxSize[3];
    std::vector<std::set<int> > exclusions;
    NonbondedMethod nonbondedMethod;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked by GBSAOBCForceField to calculate the forces acting on the system.
 */
class ReferenceCalcGBSAOBCForceFieldKernel : public CalcGBSAOBCForceFieldKernel {
public:
    ReferenceCalcGBSAOBCForceFieldKernel(std::string name, const Platform& platform) : CalcGBSAOBCForceFieldKernel(name, platform) {
    }
    ~ReferenceCalcGBSAOBCForceFieldKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForceField this kernel will be used for
     */
    void initialize(const System& system, const GBSAOBCForceField& force);
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCForceField
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    CpuObc* obc;
    std::vector<RealOpenMM> charges;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class ReferenceIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    ReferenceIntegrateVerletStepKernel(std::string name, const Platform& platform) : IntegrateVerletStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateVerletStepKernel();
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
    void execute(OpenMMContextImpl& context, const VerletIntegrator& integrator);
private:
    ReferenceVerletDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevStepSize;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class ReferenceIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    ReferenceIntegrateLangevinStepKernel(std::string name, const Platform& platform) : IntegrateLangevinStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateLangevinStepKernel();
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
    void execute(OpenMMContextImpl& context, const LangevinIntegrator& integrator);
private:
    ReferenceStochasticDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class ReferenceIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    ReferenceIntegrateBrownianStepKernel(std::string name, const Platform& platform) : IntegrateBrownianStepKernel(name, platform),
        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
    }
    ~ReferenceIntegrateBrownianStepKernel();
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
    void execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator);
private:
    ReferenceBrownianDynamics* dynamics;
    ReferenceShakeAlgorithm* shake;
    RealOpenMM* masses;
    RealOpenMM** shakeParameters;
    int** constraintIndices;
    int numConstraints;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class ReferenceApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    ReferenceApplyAndersenThermostatKernel(std::string name, const Platform& platform) : ApplyAndersenThermostatKernel(name, platform), thermostat(0) {
    }
    ~ReferenceApplyAndersenThermostatKernel();
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
    void execute(OpenMMContextImpl& context);
private:
    ReferenceAndersenThermostat* thermostat;
    RealOpenMM* masses;
};

/**
 * This kernel is invoked to calculate the kinetic energy of the system.
 */
class ReferenceCalcKineticEnergyKernel : public CalcKineticEnergyKernel {
public:
    ReferenceCalcKineticEnergyKernel(std::string name, const Platform& platform) : CalcKineticEnergyKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    double execute(OpenMMContextImpl& context);
private:
    std::vector<double> masses;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class ReferenceRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    ReferenceRemoveCMMotionKernel(std::string name, const Platform& platform) : RemoveCMMotionKernel(name, platform) {
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
    void execute(OpenMMContextImpl& context);
private:
    std::vector<double> masses;
    int frequency;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEKERNELS_H_*/
