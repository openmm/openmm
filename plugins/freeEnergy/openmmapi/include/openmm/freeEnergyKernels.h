#ifndef OPENMM_FREE_ENERGY_KERNELS_H_
#define OPENMM_FREE_ENERGY_KERNELS_H_

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

#include "openmm/GBSAOBCSoftcoreForce.h"
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/NonbondedSoftcoreForce.h"
#include "openmm/System.h"
#include "openmm/KernelImpl.h"
#include <set>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by GBSAOBCSoftcoreForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGBSAOBCSoftcoreForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGBSAOBCSoftcoreForce";
    }
    CalcGBSAOBCSoftcoreForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCSoftcoreForce this kernel will be used for
     */
    virtual void initialize(const System& system, const GBSAOBCSoftcoreForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(ContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBSAOBCSoftcoreForce
     */
    virtual double executeEnergy(ContextImpl& context) = 0;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcNonbondedSoftcoreForceKernel : public KernelImpl {
public:
    enum NonbondedSoftcoreMethod {
        NoCutoff          = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic    = 2,
        Ewald             = 3,
        PME               = 4
    };
    static std::string Name() {
        return "CalcNonbondedSoftcoreForce";
    }
    CalcNonbondedSoftcoreForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedSoftcoreForce this kernel will be used for
     */
    virtual void initialize(const System& system, const NonbondedSoftcoreForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(ContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the NonbondedSoftcoreForce
     */
    virtual double executeEnergy(ContextImpl& context) = 0;
};

/**
 * This kernel is invoked by GBVISoftcoreForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcGBVISoftcoreForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcGBVISoftcoreForce";
    }
    CalcGBVISoftcoreForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system      the System this kernel will be applied to
     * @param force       the GBVISoftcoreForce this kernel will be used for
     * @param scaledRadii scaled radii
     */
    virtual void initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double>& scaledRadii) = 0;
    /**
     * Execute the kernel to calculate the forces.
     * 
     * @param context    the context in which to execute this kernel
     */
    virtual void executeForces(ContextImpl& context) = 0;
    /**
     * Execute the kernel to calculate the energy.
     * 
     * @param context    the context in which to execute this kernel
     * @return the potential energy due to the GBVISoftcoreForce
     */
    virtual double executeEnergy(ContextImpl& context) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_FREE_ENERGY_KERNELS_H_*/
