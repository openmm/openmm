#ifndef DRUDE_KERNELS_H_
#define DRUDE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "openmm/DrudeForce.h"
#include "openmm/DrudeLangevinIntegrator.h"
#include "openmm/DrudeSCFIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by DrudeForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcDrudeForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcDrudeForce";
    }
    CalcDrudeForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the DrudeForce this kernel will be used for
     */
    virtual void initialize(const System& system, const DrudeForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the DrudeForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const DrudeForce& force) = 0;
};

/**
 * This kernel is invoked by DrudeLangevinIntegrator to take one time step.
 */
class IntegrateDrudeLangevinStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateDrudeLangevinStep";
    }
    IntegrateDrudeLangevinStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the DrudeLangevinIntegrator this kernel will be used for
     * @param force      the DrudeForce to get particle parameters from
     */
    virtual void initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the DrudeLangevinIntegrator this kernel is being used for
     */
    virtual void execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by DrudeSCFIntegrator to take one time step.
 */
class IntegrateDrudeSCFStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateDrudeSCFStep";
    }
    IntegrateDrudeSCFStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the DrudeSCFIntegrator this kernel will be used for
     * @param force      the DrudeForce to get particle parameters from
     */
    virtual void initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the DrudeSCFIntegrator this kernel is being used for
     */
    virtual void execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) = 0;
};

} // namespace OpenMM

#endif /*DRUDE_KERNELS_H_*/
