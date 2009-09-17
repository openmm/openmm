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

#include "OpenCLPlatform.h"
#include "openmm/kernels.h"
#include "openmm/System.h"

namespace OpenMM {

///**
// * This kernel is invoked at the start of each force evaluation to clear the forces.
// */
//class OpenCLInitializeForcesKernel : public InitializeForcesKernel {
//public:
//    OpenCLInitializeForcesKernel(std::string name, const Platform& platform) : InitializeForcesKernel(name, platform) {
//    }
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     */
//    void initialize(const System& system);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void execute(ContextImpl& context);
//};
//
///**
// * This kernel is invoked to get or set the current time.
// */
//class OpenCLUpdateTimeKernel : public UpdateTimeKernel {
//public:
//    OpenCLUpdateTimeKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : UpdateTimeKernel(name, platform), data(data) {
//    }
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     */
//    void initialize(const System& system);
//    /**
//     * Get the current time (in picoseconds).
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    double getTime(const ContextImpl& context) const;
//    /**
//     * Set the current time (in picoseconds).
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void setTime(ContextImpl& context, double time);
//private:
//    OpenCLPlatform::PlatformData& data;
//};
//
///**
// * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
// */
//class OpenCLCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
//public:
//    OpenCLCalcHarmonicBondForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcHarmonicBondForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcHarmonicBondForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the HarmonicBondForce this kernel will be used for
//     */
//    void initialize(const System& system, const HarmonicBondForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the HarmonicBondForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    int numBonds;
//    OpenCLPlatform::PlatformData& data;
//    System& system;
//};
//
///**
// * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
// */
//class OpenCLCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
//public:
//    OpenCLCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcHarmonicAngleForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcHarmonicAngleForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the HarmonicAngleForce this kernel will be used for
//     */
//    void initialize(const System& system, const HarmonicAngleForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the HarmonicAngleForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    int numAngles;
//    OpenCLPlatform::PlatformData& data;
//    System& system;
//};
//
///**
// * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
// */
//class OpenCLCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
//public:
//    OpenCLCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcPeriodicTorsionForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcPeriodicTorsionForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the PeriodicTorsionForce this kernel will be used for
//     */
//    void initialize(const System& system, const PeriodicTorsionForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the PeriodicTorsionForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    int numTorsions;
//    OpenCLPlatform::PlatformData& data;
//    System& system;
//};
//
///**
// * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
// */
//class OpenCLCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
//public:
//    OpenCLCalcRBTorsionForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcRBTorsionForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcRBTorsionForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the RBTorsionForce this kernel will be used for
//     */
//    void initialize(const System& system, const RBTorsionForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the RBTorsionForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    int numTorsions;
//    OpenCLPlatform::PlatformData& data;
//    System& system;
//};
//
///**
// * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
// */
//class OpenCLCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
//public:
//    OpenCLCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcNonbondedForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcNonbondedForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the NonbondedForce this kernel will be used for
//     */
//    void initialize(const System& system, const NonbondedForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the NonbondedForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    OpenCLPlatform::PlatformData& data;
//    int numParticles;
//    System& system;
//};
//
///**
// * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
// */
//class OpenCLCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
//public:
//    OpenCLCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data, System& system) : CalcCustomNonbondedForceKernel(name, platform), data(data), system(system) {
//    }
//    ~OpenCLCalcCustomNonbondedForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the CustomNonbondedForce this kernel will be used for
//     */
//    void initialize(const System& system, const CustomNonbondedForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the CustomNonbondedForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    void updateGlobalParams(ContextImpl& context);
//    OpenCLPlatform::PlatformData& data;
//    int numParticles;
//    std::vector<std::string> globalParamNames;
//    std::vector<float> globalParamValues;
//    System& system;
//};
//
///**
// * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
// */
//class OpenCLCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
//public:
//    OpenCLCalcGBSAOBCForceKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : CalcGBSAOBCForceKernel(name, platform), data(data) {
//    }
//    ~OpenCLCalcGBSAOBCForceKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the GBSAOBCForce this kernel will be used for
//     */
//    void initialize(const System& system, const GBSAOBCForce& force);
//    /**
//     * Execute the kernel to calculate the forces.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void executeForces(ContextImpl& context);
//    /**
//     * Execute the kernel to calculate the energy.
//     *
//     * @param context    the context in which to execute this kernel
//     * @return the potential energy due to the GBSAOBCForce
//     */
//    double executeEnergy(ContextImpl& context);
//private:
//    OpenCLPlatform::PlatformData& data;
//};
//
///**
// * This kernel is invoked by VerletIntegrator to take one time step.
// */
//class OpenCLIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
//public:
//    OpenCLIntegrateVerletStepKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : IntegrateVerletStepKernel(name, platform), data(data) {
//    }
//    ~OpenCLIntegrateVerletStepKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param integrator the VerletIntegrator this kernel will be used for
//     */
//    void initialize(const System& system, const VerletIntegrator& integrator);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     * @param integrator the VerletIntegrator this kernel is being used for
//     */
//    void execute(ContextImpl& context, const VerletIntegrator& integrator);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevStepSize;
//};
//
///**
// * This kernel is invoked by LangevinIntegrator to take one time step.
// */
//class OpenCLIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
//public:
//    OpenCLIntegrateLangevinStepKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : IntegrateLangevinStepKernel(name, platform), data(data) {
//    }
//    ~OpenCLIntegrateLangevinStepKernel();
//    /**
//     * Initialize the kernel, setting up the particle masses.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param integrator the LangevinIntegrator this kernel will be used for
//     */
//    void initialize(const System& system, const LangevinIntegrator& integrator);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     * @param integrator the LangevinIntegrator this kernel is being used for
//     */
//    void execute(ContextImpl& context, const LangevinIntegrator& integrator);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevTemp, prevFriction, prevStepSize;
//};
//
///**
// * This kernel is invoked by BrownianIntegrator to take one time step.
// */
//class OpenCLIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
//public:
//    OpenCLIntegrateBrownianStepKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : IntegrateBrownianStepKernel(name, platform), data(data) {
//    }
//    ~OpenCLIntegrateBrownianStepKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param integrator the BrownianIntegrator this kernel will be used for
//     */
//    void initialize(const System& system, const BrownianIntegrator& integrator);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     * @param integrator the BrownianIntegrator this kernel is being used for
//     */
//    void execute(ContextImpl& context, const BrownianIntegrator& integrator);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevTemp, prevFriction, prevStepSize;
//};
//
///**
// * This kernel is invoked by VariableVerletIntegrator to take one time step.
// */
//class OpenCLIntegrateVariableVerletStepKernel : public IntegrateVariableVerletStepKernel {
//public:
//    OpenCLIntegrateVariableVerletStepKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : IntegrateVariableVerletStepKernel(name, platform), data(data) {
//    }
//    ~OpenCLIntegrateVariableVerletStepKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param integrator the VerletIntegrator this kernel will be used for
//     */
//    void initialize(const System& system, const VariableVerletIntegrator& integrator);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     * @param integrator the VerletIntegrator this kernel is being used for
//     * @param maxTime    the maximum time beyond which the simulation should not be advanced
//     */
//    void execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevErrorTol;
//};
//
///**
// * This kernel is invoked by VariableLangevinIntegrator to take one time step.
// */
//class OpenCLIntegrateVariableLangevinStepKernel : public IntegrateVariableLangevinStepKernel {
//public:
//    OpenCLIntegrateVariableLangevinStepKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : IntegrateVariableLangevinStepKernel(name, platform), data(data) {
//    }
//    ~OpenCLIntegrateVariableLangevinStepKernel();
//    /**
//     * Initialize the kernel, setting up the particle masses.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param integrator the VariableLangevinIntegrator this kernel will be used for
//     */
//    void initialize(const System& system, const VariableLangevinIntegrator& integrator);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     * @param integrator the VariableLangevinIntegrator this kernel is being used for
//     * @param maxTime    the maximum time beyond which the simulation should not be advanced
//     */
//    void execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevTemp, prevFriction, prevErrorTol;
//};
//
///**
// * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
// */
//class OpenCLApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
//public:
//    OpenCLApplyAndersenThermostatKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : ApplyAndersenThermostatKernel(name, platform), data(data) {
//    }
//    ~OpenCLApplyAndersenThermostatKernel();
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param thermostat the AndersenThermostat this kernel will be used for
//     */
//    void initialize(const System& system, const AndersenThermostat& thermostat);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void execute(ContextImpl& context);
//private:
//    OpenCLPlatform::PlatformData& data;
//    double prevTemp, prevFrequency, prevStepSize;
//};
//
///**
// * This kernel is invoked to calculate the kinetic energy of the system.
// */
//class OpenCLCalcKineticEnergyKernel : public CalcKineticEnergyKernel {
//public:
//    OpenCLCalcKineticEnergyKernel(std::string name, const Platform& platform) : CalcKineticEnergyKernel(name, platform) {
//    }
//    /**
//     * Initialize the kernel.
//     *
//     * @param system     the System this kernel will be applied to
//     */
//    void initialize(const System& system);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    double execute(ContextImpl& context);
//private:
//    std::vector<double> masses;
//};
//
///**
// * This kernel is invoked to remove center of mass motion from the system.
// */
//class OpenCLRemoveCMMotionKernel : public RemoveCMMotionKernel {
//public:
//    OpenCLRemoveCMMotionKernel(std::string name, const Platform& platform, OpenCLPlatform::PlatformData& data) : RemoveCMMotionKernel(name, platform), data(data) {
//    }
//    /**
//     * Initialize the kernel, setting up the particle masses.
//     *
//     * @param system     the System this kernel will be applied to
//     * @param force      the CMMotionRemover this kernel will be used for
//     */
//    void initialize(const System& system, const CMMotionRemover& force);
//    /**
//     * Execute the kernel.
//     *
//     * @param context    the context in which to execute this kernel
//     */
//    void execute(ContextImpl& context);
//private:
//    OpenCLPlatform::PlatformData& data;
//};

} // namespace OpenMM

#endif /*OPENMM_OPENCLKERNELS_H_*/
