#ifndef RPMD_KERNELS_H_
#define RPMD_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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

#include "openmm/RPMDIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by RPMDIntegrator to take one time step, and to get and
 * set the state of system copies.
 */
class IntegrateRPMDStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateRPMDStep";
    }
    IntegrateRPMDStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the RPMDIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const RPMDIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    virtual void execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) = 0;
    /**
     * Get the positions of all particles in one copy of the system.
     */
    virtual void setPositions(int copy, const std::vector<Vec3>& positions) = 0;
    /**
     * Get the velocities of all particles in one copy of the system.
     */
    virtual void setVelocities(int copy, const std::vector<Vec3>& velocities) = 0;
    /**
     * Copy positions and velocities for one copy into the context.
     */
    virtual void copyToContext(int copy, ContextImpl& context) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) = 0;
};

} // namespace OpenMM

#endif /*RPMD_KERNELS_H_*/
