#ifndef OPENCL_RPMD_KERNELS_H_
#define OPENCL_RPMD_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2018 Stanford University and the Authors.      *
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

#include "openmm/RpmdKernels.h"
#include "OpenCLContext.h"
#include "OpenCLArray.h"
#include <map>
namespace OpenMM {

/**
 * This kernel is invoked by RPMDIntegrator to take one time step, and to get and
 * set the state of system copies.
 */
class OpenCLIntegrateRPMDStepKernel : public IntegrateRPMDStepKernel {
public:
    OpenCLIntegrateRPMDStepKernel(std::string name, const Platform& platform, OpenCLContext& cl) :
            IntegrateRPMDStepKernel(name, platform), cl(cl), hasInitializedKernel(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the RPMDIntegrator this kernel will be used for
     */
    void initialize(const System& system, const RPMDIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    void execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid);
    /**
     * Compute the kinetic energy.
     * 
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator);
    /**
     * Get the positions of all particles in one copy of the system.
     */
    void setPositions(int copy, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles in one copy of the system.
     */
    void setVelocities(int copy, const std::vector<Vec3>& velocities);
    /**
     * Copy positions and velocities for one copy into the context.
     */
    void copyToContext(int copy, ContextImpl& context);
private:
    void initializeKernels(ContextImpl& context);
    void computeForces(ContextImpl& context);
    std::string createFFT(int size, const std::string& variable, bool forward);
    OpenCLContext& cl;
    bool hasInitializedKernel;
    int numCopies, numParticles, workgroupSize;
    std::map<int, int> groupsByCopies;
    int groupsNotContracted;
    OpenCLArray forces;
    OpenCLArray positions;
    OpenCLArray velocities;
    OpenCLArray contractedForces;
    OpenCLArray contractedPositions;
    cl::Kernel pileKernel, stepKernel, velocitiesKernel, copyToContextKernel, copyFromContextKernel, translateKernel;
    std::map<int, cl::Kernel> positionContractionKernels;
    std::map<int, cl::Kernel> forceContractionKernels;
};

} // namespace OpenMM

#endif /*OPENCL_RPMD_KERNELS_H_*/
