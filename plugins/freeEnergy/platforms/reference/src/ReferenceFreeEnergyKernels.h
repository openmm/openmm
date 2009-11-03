#ifndef OPENMM_REFERENCE_FREE_ENERGY_KERNELS_H_
#define OPENMM_REFERENCE_FREE_ENERGY_KERNELS_H_

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

#include "ReferencePlatform.h"
#include "openmm/freeEnergyKernels.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include "SimTKReference/ReferenceNeighborList.h"
#include "gbsa/CpuGBVISoftcore.h"
#include "gbsa/CpuObcSoftcore.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedSoftcoreForce to calculate the forces acting on the system.
 */
class ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel : public CalcNonbondedSoftcoreForceKernel {
public:
    ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel(std::string name, const Platform& platform) : CalcNonbondedSoftcoreForceKernel(name, platform) {
    }
    ~ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedSoftcoreForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedSoftcoreForce& force);
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
     * @return the potential energy due to the NonbondedSoftcoreForce
     */
    double executeEnergy(ContextImpl& context);
private:
    int numParticles, num14;
    int **exclusionArray, **bonded14IndexArray;
    RealOpenMM **particleParamArray, **bonded14ParamArray;
    RealOpenMM nonbondedCutoff, periodicBoxSize[3], rfDielectric, ewaldAlpha;
    int kmax[3];
    std::vector<std::set<int> > exclusions;
    NonbondedSoftcoreMethod nonbondedMethod;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked by GBSAOBCForce to calculate the forces acting on the system.
 */
class ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel : public CalcGBSAOBCSoftcoreForceKernel {
public:
    ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel(std::string name, const Platform& platform) : CalcGBSAOBCSoftcoreForceKernel(name, platform) {
    }
    ~ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCSoftcoreForce this kernel will be used for
     */
    void initialize(const System& system, const GBSAOBCSoftcoreForce& force);
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
     * @return the potential energy due to the GBSAOBCSoftcoreForce
     */
    double executeEnergy(ContextImpl& context);
private:
    CpuObcSoftcore* obc;
    std::vector<RealOpenMM> charges;
};

/**
 * This kernel is invoked by GBVISoftcoreForce to calculate the forces acting on the system.
 */
class ReferenceFreeEnergyCalcGBVISoftcoreForceKernel : public CalcGBVISoftcoreForceKernel {
public:
    ReferenceFreeEnergyCalcGBVISoftcoreForceKernel(std::string name, const Platform& platform) : CalcGBVISoftcoreForceKernel(name, platform) {
    }
    ~ReferenceFreeEnergyCalcGBVISoftcoreForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system       the System this kernel will be applied to
     * @param force        the GBVISoftcoreForce this kernel will be used for
     * @param scaled radii the scaled radii (Eq. 5 of Labute paper)
     */
    void initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double> & scaledRadii);
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
     * @return the potential energy due to the GBVISoftcoreForce
     */
    double executeEnergy(ContextImpl& context);
private:
    CpuGBVISoftcore* gbviSoftcore;
    std::vector<RealOpenMM> charges;
};


} // namespace OpenMM

#endif /*OPENMM_REFERENCE_FREE_ENERGY_KERNELS_H_*/
