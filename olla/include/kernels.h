#ifndef OPENMM_KERNELS_H_
#define OPENMM_KERNELS_H_

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

#include "AndersenThermostat.h"
#include "BrownianIntegrator.h"
#include "CMMotionRemover.h"
#include "GBSAOBCForce.h"
#include "GBVIForce.h"
#include "HarmonicAngleForce.h"
#include "HarmonicBondForce.h"
#include "KernelImpl.h"
#include "LangevinIntegrator.h"
#include "PeriodicTorsionForce.h"
#include "RBTorsionForce.h"
#include "NonbondedForce.h"
#include "Stream.h"
#include "System.h"
#include "VerletIntegrator.h"
#include <set>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked at the start of each force evaluation to clear the forces.
 */
class InitializeForcesKernel : public KernelImpl {
public:
    static std::string Name() {
        return "InitializeForces";
    }
    InitializeForcesKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void execute(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcHarmonicBondForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcHarmonicBondForce";
    }
    CalcHarmonicBondForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicBondForce this kernel will be used for
     */
    virtual void initialize(const System& system, const HarmonicBondForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicBondForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcHarmonicAngleForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcHarmonicAngleForce";
    }
    CalcHarmonicAngleForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HarmonicAngleForce this kernel will be used for
     */
    virtual void initialize(const System& system, const HarmonicAngleForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicAngleForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcPeriodicTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcPeriodicTorsionForce";
    }
    CalcPeriodicTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the PeriodicTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const PeriodicTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the PeriodicTorsionForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcRBTorsionForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcRBTorsionForce";
    }
    CalcRBTorsionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the RBTorsionForce this kernel will be used for
     */
    virtual void initialize(const System& system, const RBTorsionForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the RBTorsionForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcNonbondedForceKernel : public KernelImpl {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2,
        Ewald = 3
    };
    static std::string Name() {
        return "CalcNonbondedForce";
    }
    CalcNonbondedForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    virtual void initialize(const System& system, const NonbondedForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the NonbondedForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGBSAOBCForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGBSAOBCForces";
    }
    CalcGBSAOBCForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const System& system, const GBSAOBCForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by GBVIForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGBVIForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGBVIForces";
    }
    CalcGBVIForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system      the System this kernel will be applied to
     * @param force       the GBVIForce this kernel will be used for
     * @param scaledRadii scaled radii
     */
    virtual void initialize(const System& system, const GBVIForce& force, const std::vector<double>& scaledRadii) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(OpenMMContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBVIForce
     */
    virtual double executeEnergy(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked by VerletIntegrator to take one time step.
 */
class IntegrateVerletStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateVerletStep";
    }
    IntegrateVerletStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const VerletIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     */
    virtual void execute(OpenMMContextImpl& context, const VerletIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class IntegrateLangevinStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateLangevinStep";
    }
    IntegrateLangevinStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the LangevinIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const LangevinIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the LangevinIntegrator this kernel is being used for
     */
    virtual void execute(OpenMMContextImpl& context, const LangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class IntegrateBrownianStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateBrownianStep";
    }
    IntegrateBrownianStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the BrownianIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const BrownianIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the BrownianIntegrator this kernel is being used for
     */
    virtual void execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class ApplyAndersenThermostatKernel : public KernelImpl {
public:
    static std::string Name() {
        return "ApplyAndersenThermostat";
    }
    ApplyAndersenThermostatKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param thermostat the AndersenThermostat this kernel will be used for
     */
    virtual void initialize(const System& system, const AndersenThermostat& thermostat) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void execute(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked to calculate the kinetic energy of the system.
 */
class CalcKineticEnergyKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcKineticEnergy";
    }
    CalcKineticEnergyKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual double execute(OpenMMContextImpl& context) = 0;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class RemoveCMMotionKernel : public KernelImpl {
public:
    static std::string Name() {
        return "RemoveCMMotion";
    }
    RemoveCMMotionKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the CMMotionRemover this kernel will be used for
     */
    virtual void initialize(const System& system, const CMMotionRemover& force) = 0;
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void execute(OpenMMContextImpl& context) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_KERNELS_H_*/
