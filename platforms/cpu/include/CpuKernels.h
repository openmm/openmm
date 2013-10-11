#ifndef OPENMM_CPUKERNELS_H_
#define OPENMM_CPUKERNELS_H_

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

#include "CpuPlatform.h"
#include "CpuNeighborList.h"
#include "openmm/kernels.h"
#include "openmm/System.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class CpuCalcNonbondedForceKernel : public CalcNonbondedForceKernel {
public:
    CpuCalcNonbondedForceKernel(std::string name, const Platform& platform) : CalcNonbondedForceKernel(name, platform), bonded14IndexArray(NULL), bonded14ParamArray(NULL) {
    }
    ~CpuCalcNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @param includeDirect  true if direct space interactions should be included
     * @param includeReciprocal  true if reciprocal space interactions should be included
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the NonbondedForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const NonbondedForce& force);
private:
    int numParticles, num14;
    int **bonded14IndexArray;
    double **bonded14ParamArray;
    double nonbondedCutoff, switchingDistance, rfDielectric, ewaldAlpha, ewaldSelfEnergy, dispersionCoefficient;
    int kmax[3], gridSize[3];
    bool useSwitchingFunction;
    std::vector<std::set<int> > exclusions;
    std::vector<std::pair<float, float> > particleParams;
    std::vector<float> posq;
    std::vector<float> forces;
    NonbondedMethod nonbondedMethod;
    CpuNeighborList neighborList;
};

} // namespace OpenMM

#endif /*OPENMM_CPUKERNELS_H_*/

