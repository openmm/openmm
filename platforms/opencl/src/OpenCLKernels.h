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
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "OpenCLParameterSet.h"
#include "openmm/kernels.h"
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
     * This is called at the beginning of each force computation, before calcForces() has been called on
     * any ForceImpl.
     *
     * @param context    the context in which to execute this kernel
     */
    void beginForceComputation(ContextImpl& context);
    /**
     * This is called at the end of each force computation, after calcForces() has been called on
     * every ForceImpl.
     *
     * @param context    the context in which to execute this kernel
     */
    void finishForceComputation(ContextImpl& context);
    /**
     * This is called at the beginning of each energy computation, before calcEnergy() has been called on
     * any ForceImpl.
     *
     * @param context    the context in which to execute this kernel
     */
    void beginEnergyComputation(ContextImpl& context);
    /**
     * This is called at the end of each energy computation, after calcEnergy() has been called on
     * every ForceImpl.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy of the system.  This value is added to all values returned by ForceImpls'
     * calcEnergy() methods.  That is, each force kernel may <i>either</i> return its contribution to the
     * energy directly, <i>or</i> add it to an internal buffer so that it will be included here.
     */
    double finishEnergyComputation(ContextImpl& context);
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
private:
    OpenCLContext& cl;
};

/**
 * This kernel is invoked by HarmonicBondForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcHarmonicBondForceKernel : public CalcHarmonicBondForceKernel {
public:
    OpenCLCalcHarmonicBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcHarmonicBondForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL), indices(NULL) {
    }
    ~OpenCLCalcHarmonicBondForceKernel();
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
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicBondForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numBonds;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLArray<mm_float2>* params;
    OpenCLArray<mm_int4>* indices;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by CustomBondForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomBondForceKernel : public CalcCustomBondForceKernel {
public:
    OpenCLCalcCustomBondForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcCustomBondForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL), indices(NULL), globals(NULL) {
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
     * Execute the kernel to calculate the forces.
     *
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the CustomBondForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numBonds;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLParameterSet* params;
    OpenCLArray<mm_int4>* indices;
    OpenCLArray<cl_float>* globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by HarmonicAngleForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcHarmonicAngleForceKernel : public CalcHarmonicAngleForceKernel {
public:
    OpenCLCalcHarmonicAngleForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcHarmonicAngleForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
    }
    ~OpenCLCalcHarmonicAngleForceKernel();
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
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the HarmonicAngleForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numAngles;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLArray<mm_float2>* params;
    OpenCLArray<mm_int8>* indices;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by PeriodicTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcPeriodicTorsionForceKernel : public CalcPeriodicTorsionForceKernel {
public:
    OpenCLCalcPeriodicTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcPeriodicTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
    }
    ~OpenCLCalcPeriodicTorsionForceKernel();
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
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the PeriodicTorsionForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numTorsions;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLArray<mm_float4>* params;
    OpenCLArray<mm_int8>* indices;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by RBTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcRBTorsionForceKernel : public CalcRBTorsionForceKernel {
public:
    OpenCLCalcRBTorsionForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcRBTorsionForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system) {
    }
    ~OpenCLCalcRBTorsionForceKernel();
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
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the RBTorsionForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numTorsions;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLArray<mm_float8>* params;
    OpenCLArray<mm_int8>* indices;
    cl::Kernel kernel;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class OpenCLCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    OpenCLCalcNonbondedForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcNonbondedForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), sigmaEpsilon(NULL), exceptionParams(NULL), exceptionIndices(NULL), cosSinSums(NULL) {
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
     * Execute the kernel to calculate the forces.
     *
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the NonbondedForce
     */
    double executeEnergy(ContextImpl& context);
private:
    OpenCLContext& cl;
    bool hasInitializedKernel;
    OpenCLArray<mm_float2>* sigmaEpsilon;
    OpenCLArray<mm_float4>* exceptionParams;
    OpenCLArray<mm_int4>* exceptionIndices;
    OpenCLArray<mm_float2>* cosSinSums;
    cl::Kernel exceptionsKernel;
    cl::Kernel ewaldSumsKernel;
    cl::Kernel ewaldForcesKernel;
    double cutoffSquared, ewaldSelfEnergy;
};

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class OpenCLCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    OpenCLCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcCustomNonbondedForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), params(NULL), globals(NULL), tabulatedFunctionParams(NULL), system(system) {
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
     * Execute the kernel to calculate the forces.
     *
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the CustomNonbondedForce
     */
    double executeEnergy(ContextImpl& context);
private:
    bool hasInitializedKernel;
    OpenCLContext& cl;
    OpenCLParameterSet* params;
    OpenCLArray<cl_float>* globals;
    OpenCLArray<mm_float4>* tabulatedFunctionParams;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray<mm_float4>*> tabulatedFunctions;
    System& system;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class OpenCLCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {
public:
    OpenCLCalcGBSAOBCForceKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcGBSAOBCForceKernel(name, platform), cl(cl), hasCreatedKernels(false) {
    }
    ~OpenCLCalcGBSAOBCForceKernel();
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
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCForce
     */
    double executeEnergy(ContextImpl& context);
private:
    double prefactor;
    bool hasCreatedKernels;
    OpenCLContext& cl;
    OpenCLArray<mm_float2>* params;
    OpenCLArray<cl_float>* bornSum;
    OpenCLArray<cl_float>* bornRadii;
    OpenCLArray<cl_float>* bornForce;
    OpenCLArray<cl_float>* obcChain;
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
    OpenCLCalcCustomGBForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcCustomGBForceKernel(name, platform),
            hasInitializedKernels(false), cl(cl), params(NULL), globals(NULL), valueBuffers(NULL), tabulatedFunctionParams(NULL), system(system) {
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
     * Execute the kernel to calculate the forces.
     *
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the CustomGBForce
     */
    double executeEnergy(ContextImpl& context);
private:
    bool hasInitializedKernels;
    OpenCLContext& cl;
    OpenCLParameterSet* params;
    OpenCLParameterSet* computedValues;
    OpenCLArray<cl_float>* globals;
    OpenCLArray<cl_float>* valueBuffers;
    OpenCLArray<mm_float4>* tabulatedFunctionParams;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    std::vector<OpenCLArray<mm_float4>*> tabulatedFunctions;
    System& system;
    cl::Kernel pairValueKernel, perParticleValueKernel, pairEnergyKernel, perParticleEnergyKernel;
};

/**
 * This kernel is invoked by CustomExternalForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcCustomExternalForceKernel : public CalcCustomExternalForceKernel {
public:
    OpenCLCalcCustomExternalForceKernel(std::string name, const Platform& platform, OpenCLContext& cl, System& system) : CalcCustomExternalForceKernel(name, platform),
            hasInitializedKernel(false), cl(cl), system(system), params(NULL), indices(NULL), globals(NULL) {
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
     * Execute the kernel to calculate the forces.
     *
     * @param context    the context in which to execute this kernel
     */
    void executeForces(ContextImpl& context);
    /**
     * Execute the kernel to calculate the energy.
     *
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the CustomExternalForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numParticles;
    bool hasInitializedKernel;
    OpenCLContext& cl;
    System& system;
    OpenCLParameterSet* params;
    OpenCLArray<cl_int>* indices;
    OpenCLArray<cl_float>* globals;
    std::vector<std::string> globalParamNames;
    std::vector<cl_float> globalParamValues;
    cl::Kernel kernel;
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
private:
    OpenCLContext& cl;
    double prevStepSize;
    bool hasInitializedKernels;
    cl::Kernel kernel1, kernel2;
};

/**
 * This kernel is invoked by LangevinIntegrator to take one time step.
 */
class OpenCLIntegrateLangevinStepKernel : public IntegrateLangevinStepKernel {
public:
    OpenCLIntegrateLangevinStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateLangevinStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false), params(NULL), xVector(NULL), vVector(NULL) {
    }
    ~OpenCLIntegrateLangevinStepKernel();
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
private:
    OpenCLContext& cl;
    double prevTemp, prevFriction, prevStepSize;
    bool hasInitializedKernels;
    OpenCLArray<cl_float>* params;
    OpenCLArray<mm_float4>* xVector;
    OpenCLArray<mm_float4>* vVector;
    cl::Kernel kernel1, kernel2, kernel3;
};

/**
 * This kernel is invoked by BrownianIntegrator to take one time step.
 */
class OpenCLIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {
public:
    OpenCLIntegrateBrownianStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) : IntegrateBrownianStepKernel(name, platform), cl(cl) {
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
     * @param integrator the VerletIntegrator this kernel will be used for
     */
    void initialize(const System& system, const VariableVerletIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the VerletIntegrator this kernel is being used for
     * @param maxTime    the maximum time beyond which the simulation should not be advanced
     */
    void execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime);
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
    ~OpenCLIntegrateVariableLangevinStepKernel();
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
     */
    void execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime);
private:
    OpenCLContext& cl;
    bool hasInitializedKernels;
    int blockSize;
    OpenCLArray<cl_float>* params;
    OpenCLArray<mm_float4>* xVector;
    OpenCLArray<mm_float4>* vVector;
    cl::Kernel kernel1, kernel2, kernel3, selectSizeKernel;
    double prevTemp, prevFriction, prevErrorTol;
};

/**
 * This kernel is invoked by AndersenThermostat at the start of each time step to adjust the particle velocities.
 */
class OpenCLApplyAndersenThermostatKernel : public ApplyAndersenThermostatKernel {
public:
    OpenCLApplyAndersenThermostatKernel(std::string name, const Platform& platform, OpenCLContext& cl) : ApplyAndersenThermostatKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
    }
    ~OpenCLApplyAndersenThermostatKernel();
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
    cl::Kernel kernel;
    double prevTemp, prevFriction, prevStepSize;
};

/**
 * This kernel is invoked to calculate the kinetic energy of the system.
 */
class OpenCLCalcKineticEnergyKernel : public CalcKineticEnergyKernel {
public:
    OpenCLCalcKineticEnergyKernel(std::string name, const Platform& platform, OpenCLContext& cl) : CalcKineticEnergyKernel(name, platform), cl(cl) {
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
    double execute(ContextImpl& context);
private:
    OpenCLContext& cl;
    std::vector<double> masses;
};

/**
 * This kernel is invoked to remove center of mass motion from the system.
 */
class OpenCLRemoveCMMotionKernel : public RemoveCMMotionKernel {
public:
    OpenCLRemoveCMMotionKernel(std::string name, const Platform& platform, OpenCLContext& cl) : RemoveCMMotionKernel(name, platform), cl(cl), cmMomentum(NULL) {
    }
    ~OpenCLRemoveCMMotionKernel();
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
    OpenCLArray<mm_float4>* cmMomentum;
    cl::Kernel kernel1, kernel2;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLKERNELS_H_*/
