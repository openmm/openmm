#ifndef OPENMM_HIPKERNELS_H_
#define OPENMM_HIPKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Portions copyright (c) 2020-2022 Advanced Micro Devices, Inc.              *
 * Authors: Peter Eastman, Nicholas Curtis                                    *
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

#include "HipPlatform.h"
#include "HipArray.h"
#include "HipContext.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/common/CommonCalcNonbondedForce.h"
#include "openmm/common/ComputeSort.h"
#include "openmm/common/FFT3D.h"

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class HipCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    HipCalcForcesAndEnergyKernel(std::string name, const Platform& platform, HipContext& cu) : CalcForcesAndEnergyKernel(name, platform), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     */
    void initialize(const System& system);
    /**
     * This is called at the beginning of each force/energy computation, before calcForcesAndEnergy() has been called on
     * any ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     */
    void beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups);
    /**
     * This is called at the end of each force/energy computation, after calcForcesAndEnergy() has been called on
     * every ForceImpl.
     *
     * @param context       the context in which to execute this kernel
     * @param includeForce  true if forces should be computed
     * @param includeEnergy true if potential energy should be computed
     * @param groups        a set of bit flags for which force groups to include
     * @param valid         the method may set this to false to indicate the results are invalid and the force/energy
     *                      calculation should be repeated
     * @return the potential energy of the system.  This value is added to all values returned by ForceImpls'
     * calcForcesAndEnergy() methods.  That is, each force kernel may <i>either</i> return its contribution to the
     * energy directly, <i>or</i> add it to an internal buffer so that it will be included here.
     */
    double finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid);
private:
   HipContext& cu;
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class HipCalcNonbondedForceKernel : public CommonCalcNonbondedForceKernel {
public:
    HipCalcNonbondedForceKernel(std::string name, const Platform& platform, HipContext& cu, const System& system) :
            CommonCalcNonbondedForceKernel(name, platform, cu, system), cu(cu) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedForce& force);
private:
    HipContext& cu;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class HipCalcCustomCVForceKernel : public CommonCalcCustomCVForceKernel {
public:
    HipCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcCustomCVForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<HipPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

/**
 * This kernel is invoked by ATMForce to calculate the forces acting on the system and the energy of the system.
 */
class HipCalcATMForceKernel : public CommonCalcATMForceKernel {
public:
    HipCalcATMForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<HipPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

} // namespace OpenMM

#endif /*OPENMM_HIPKERNELS_H_*/
