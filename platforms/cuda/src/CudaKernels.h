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

#include "CudaPlatform.h"
#include "kernels.h"
#include "kernels/gpuTypes.h"
#include "System.h"

class CudaAndersenThermostat;
class CudaBrownianDynamics;
class CudaStochasticDynamics;
class CudaShakeAlgorithm;
class CudaVerletDynamics;

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class CudaCalcStandardMMForceFieldKernel : public CalcStandardMMForceFieldKernel {
public:
    CudaCalcStandardMMForceFieldKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data, System& system) : CalcStandardMMForceFieldKernel(name, platform), data(data), system(system) {
    }
    ~CudaCalcStandardMMForceFieldKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the StandardMMForceField this kernel will be used for
     * @param exclusions the i'th element lists the indices of all atoms with which the i'th atom should not interact through
     *                   nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
     *                   the standard nonbonded calculation.
     */
    void initialize(const System& system, const StandardMMForceField& force, const std::vector<std::set<int> >& exclusions);
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
     * @return the potential energy due to the StandardMMForceField
     */
    double executeEnergy(OpenMMContextImpl& context);
private:
    CudaPlatform::PlatformData& data;
    int numAtoms, numBonds, numAngles, numPeriodicTorsions, numRBTorsions, num14;
    System& system;
};

/**
 * This kernel is invoked by GBSAOBCForceField to calculate the forces acting on the system.
 */
class CudaCalcGBSAOBCForceFieldKernel : public CalcGBSAOBCForceFieldKernel {
public:
    CudaCalcGBSAOBCForceFieldKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : CalcGBSAOBCForceFieldKernel(name, platform), data(data) {
    }
    ~CudaCalcGBSAOBCForceFieldKernel();
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
    CudaPlatform::PlatformData& data;
};
//
///**
// * This kernel is invoked by VerletIntegrator to take one time step.
// */
//class CudaIntegrateVerletStepKernel : public IntegrateVerletStepKernel {
//public:
//    CudaIntegrateVerletStepKernel(std::string name, const Platform& platform) : IntegrateVerletStepKernel(name, platform),
//        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
//    }
//    ~CudaIntegrateVerletStepKernel();
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
//    void execute(OpenMMContextImpl& context, const VerletIntegrator& integrator);
//private:
//    CudaVerletDynamics* dynamics;
//    CudaShakeAlgorithm* shake;
//    RealOpenMM* masses;
//    RealOpenMM** shakeParameters;
//    int** constraintIndices;
//    int numConstraints;
//    double prevStepSize;
//};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class CudaIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    CudaIntegrateLangevinStepKernel(std::string name, const Platform& platform, CudaPlatform::PlatformData& data) : IntegrateLangevinStepKernel(name, platform), data(data) {
    }
    ~CudaIntegrateLangevinStepKernel();
    /**
     * Initialize the kernel, setting up the atomic masses.
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
//
///**
// * This kernel is invoked by BrownianIntegrator to take one time step.
// */
//class CudaIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
//public:
//    CudaIntegrateBrownianStepKernel(std::string name, const Platform& platform) : IntegrateBrownianStepKernel(name, platform),
//        dynamics(0), shake(0), masses(0), shakeParameters(0), constraintIndices(0) {
//    }
//    ~CudaIntegrateBrownianStepKernel();
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
//    void execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator);
//private:
//    CudaBrownianDynamics* dynamics;
//    CudaShakeAlgorithm* shake;
//    RealOpenMM* masses;
//    RealOpenMM** shakeParameters;
//    int** constraintIndices;
//    int numConstraints;
//    double prevTemp, prevFriction, prevStepSize;
//};
//
///**
// * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the atom velocities.
// */
//class CudaApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
//public:
//    CudaApplyAndersenThermostatKernel(std::string name, const Platform& platform) : ApplyAndersenThermostatKernel(name, platform), thermostat(0) {
//    }
//    ~CudaApplyAndersenThermostatKernel();
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
//    void execute(OpenMMContextImpl& context);
//private:
//    CudaAndersenThermostat* thermostat;
//    RealOpenMM* masses;
//};

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
     * Initialize the kernel, setting up the atomic masses.
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
