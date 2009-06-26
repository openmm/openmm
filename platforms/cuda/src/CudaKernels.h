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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
#include "openmm/kernels.h"
#include "kernels/gputypes.h"
#include "openmm/System.h"

class CudaAndersenThermostat;
class CudaBrownianDynamics;
class CudaStochasticDynamics;
class CudaShakeAlgorithm;
class CudaVerletDynamics;

namespace OpenMM {

/**
 * This kernel is invoked at the start of each force evaluation to clear the forces.
 */
class CudaInitializeForcesKernel : public InitializeForcesKernel {
public:
    CudaInitializeForcesKernel(std::string name, const Platform& platform) : InitializeForcesKernel(name, platform) {
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
    void execute(OpenMMContextImpl& context);
};

/**
 * This kernel is invoked to get or set the current time.
 */
class CudaUpdateTimeKernel : public UpdateTimeKernel {
public:
    CudaUpdateTimeKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : UpdateTimeKernel(name, platform), data(data) {
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
    double getTime(const OpenMMContextImpl& context) const;
    /**
     * Set the current time (in picoseconds).
     *
     * @param context    the context in which to execute this kernel
     */
    void setTime(OpenMMContextImpl& context, double time);
private:
    CudaPlatform::PlatformData& data;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    CudaCalcHarmonicBondForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcHarmonicBondForceKernel(name, platform), data(data), system(system) {
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
    CudaPlatform::PlatformData& data;
    System& system;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    CudaCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcHarmonicAngleForceKernel(name, platform), data(data), system(system) {
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
    CudaPlatform::PlatformData& data;
    System& system;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    CudaCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcPeriodicTorsionForceKernel(name, platform), data(data), system(system) {
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
    CudaPlatform::PlatformData& data;
    System& system;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    CudaCalcRBTorsionForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcRBTorsionForceKernel(name, platform), data(data), system(system) {
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
    CudaPlatform::PlatformData& data;
    System& system;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CudaCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CudaCalcNonbondedForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcNonbondedForceKernel(name, platform), data(data), system(system) {
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
    CudaPlatform::PlatformData& data;
    int numParticles, num14;
    System& system;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class CudaCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    CudaCalcGBSAOBCForceKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : CalcGBSAOBCForceKernel(name, platform), data(data) {
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
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    void executeForces(OpenMMContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCForce
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    CudaPlatform::PlatformData& data;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class CudaIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
public:
    CudaIntegrateVerletStepKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : IntegrateVerletStepKernel(name, platform), data(data) {
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
    void execute(OpenMMContextImpl& context, const VerletIntegrator& integrator);
private:
    CudaPlatform::PlatformData& data;
    double prevStepSize;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CudaIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    CudaIntegrateLangevinStepKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : IntegrateLangevinStepKernel(name, platform), data(data) {
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
    void execute(OpenMMContextImpl& context, const LangevinIntegrator& integrator);
private:
    CudaPlatform::PlatformData& data;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class CudaIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    CudaIntegrateBrownianStepKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : IntegrateBrownianStepKernel(name, platform), data(data) {
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
    void execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator);
private:
    CudaPlatform::PlatformData& data;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked by VariableVerletIntegrator to take one time step.
 */
class CudaIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
public:
    CudaIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : IntegrateVariableVerletStepKernel(name, platform), data(data) {
    }
    ~CudaIntegrateVariableVerletStepKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableVerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    void execute(OpenMMContextImpl& context, const VariableVerletIntegrator& integrator);
private:
    CudaPlatform::PlatformData& data;
    double prevErrorTol;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class CudaApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    CudaApplyAndersenThermostatKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : ApplyAndersenThermostatKernel(name, platform), data(data) {
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
    void execute(OpenMMContextImpl& context);
private:
    CudaPlatform::PlatformData& data;
    double prevTemp, prevFrequency, prevStepSize;
};

/**
 * This kernel is invoked to calculate the kinetic energy of the system.
 */
class CudaCalcKineticEnergyKernel : public CalcKineticEnergyKernel {
public:
    CudaCalcKineticEnergyKernel(std::string name, const Platform& platform) : CalcKineticEnergyKernel(name, platform) {
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
class CudaRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    CudaRemoveCMMotionKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : RemoveCMMotionKernel(name, platform), data(data) {
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
    CudaPlatform::PlatformData& data;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAKERNELS_H_*/
