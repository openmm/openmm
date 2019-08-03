#ifndef CUDA_DRUDE_KERNELS_H_
#define CUDA_DRUDE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
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

#include "openmm/DrudeKernels.h"
#include "CudaContext.h"
#include "CudaArray.h"
#include "lbfgs.h"

namespace OpenMM {

/**
 * This kernel is invoked by DrudeForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcDrudeForceKernel : public CalcDrudeForceKernel {
public:
    CudaCalcDrudeForceKernel(std::string name, const Platform& platform, CudaContext& cu) :
            CalcDrudeForceKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the DrudeForce this kernel will be used for
     */
    void initialize(const System& system, const DrudeForce& force);
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
     * @param force      the DrudeForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const DrudeForce& force);
private:
    CudaContext& cu;
    CudaArray particleParams;
    CudaArray pairParams;
};

/**
 * This kernel is invoked by DrudeLangevinIntegrator to take one time step
 */
class CudaIntegrateDrudeLangevinStepKernel : public IntegrateDrudeLangevinStepKernel {
public:
    CudaIntegrateDrudeLangevinStepKernel(std::string name, const Platform& platform, CudaContext& cu) :
            IntegrateDrudeLangevinStepKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the DrudeLangevinIntegrator this kernel will be used for
     * @param force      the DrudeForce to get particle parameters from
     */
    void initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the DrudeLangevinIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context     the context in which to execute this kernel
     * @param integrator  the DrudeLangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator);
private:
    CudaContext& cu;
    double prevStepSize;
    CudaArray normalParticles;
    CudaArray pairParticles;
    CUfunction kernel1, kernel2, hardwallKernel;
};

/**
 * This kernel is invoked by DrudeSCFIntegrator to take one time step
 */
class CudaIntegrateDrudeSCFStepKernel : public IntegrateDrudeSCFStepKernel {
public:
    CudaIntegrateDrudeSCFStepKernel(std::string name, const Platform& platform, CudaContext& cu) :
            IntegrateDrudeSCFStepKernel(name, platform), cu(cu), minimizerPos(NULL) {
    }
    ~CudaIntegrateDrudeSCFStepKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the DrudeSCFIntegrator this kernel will be used for
     * @param force      the DrudeForce to get particle parameters from
     */
    void initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the DrudeSCFIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const DrudeSCFIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context     the context in which to execute this kernel
     * @param integrator  the DrudeSCFIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator);
private:
    void minimize(ContextImpl& context, double tolerance);
    CudaContext& cu;
    double prevStepSize;
    std::vector<int> drudeParticles;
    lbfgsfloatval_t *minimizerPos;
    lbfgs_parameter_t minimizerParams;
    CUfunction kernel1, kernel2;
};

} // namespace OpenMM

#endif /*CUDA_DRUDE_KERNELS_H_*/
