#ifndef REFERENCE_DRUDE_KERNELS_H_
#define REFERENCE_DRUDE_KERNELS_H_

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

#include "ReferencePlatform.h"
#include "openmm/DrudeKernels.h"
#include "openmm/Vec3.h"
#include "lbfgs.h"
#include <utility>
#include <vector>

class ReferenceConstraintAlgorithm;

namespace OpenMM {

/**
 * This kernel is invoked by DrudeForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcDrudeForceKernel : public CalcDrudeForceKernel {
public:
    ReferenceCalcDrudeForceKernel(std::string name, const Platform& platform) : CalcDrudeForceKernel(name, platform) {
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
    std::vector<int> particle, particle1, particle2, particle3, particle4;
    std::vector<double> charge, polarizability, aniso12, aniso34;
    std::vector<int> pair1, pair2;
    std::vector<double> pairThole;
};

/**
 * This kernel is invoked by DrudeLangevinIntegrator to take one time step
 */
class ReferenceIntegrateDrudeLangevinStepKernel : public IntegrateDrudeLangevinStepKernel {
public:
    ReferenceIntegrateDrudeLangevinStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) :
        IntegrateDrudeLangevinStepKernel(name, platform), data(data) {
    }
    ~ReferenceIntegrateDrudeLangevinStepKernel();
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
    ReferencePlatform::PlatformData& data;
    std::vector<int> normalParticles;
    std::vector<std::pair<int, int> > pairParticles;
    std::vector<double> particleMass;
    std::vector<double> particleInvMass;
    std::vector<double> pairInvTotalMass;
    std::vector<double> pairInvReducedMass;
};

/**
 * This kernel is invoked by DrudeSCFIntegrator to take one time step
 */
class ReferenceIntegrateDrudeSCFStepKernel : public IntegrateDrudeSCFStepKernel {
public:
    ReferenceIntegrateDrudeSCFStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) :
        IntegrateDrudeSCFStepKernel(name, platform), data(data), minimizerPos(NULL) {
    }
    ~ReferenceIntegrateDrudeSCFStepKernel();
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
    ReferencePlatform::PlatformData& data;
    std::vector<int> drudeParticles;
    std::vector<double> particleInvMass;
    lbfgsfloatval_t *minimizerPos;
    lbfgs_parameter_t minimizerParams;
    double maxDrudeDistance;
};

} // namespace OpenMM

#endif /*REFERENCE_DRUDE_KERNELS_H_*/
