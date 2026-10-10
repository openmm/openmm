#ifndef OPENMM_METALKERNELS_H_
#define OPENMM_METALKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
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

#include "MetalPlatform.h"
#include "MetalArray.h"
#include "MetalContext.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/common/CommonCalcNonbondedForce.h"
#include "openmm/common/CommonCalcConstantPotentialForce.h"
#include "openmm/common/ComputeSort.h"
#include "openmm/common/FFT3D.h"

namespace OpenMM {

/**
 * This kernel is invoked at the beginning and end of force and energy computations.  It gives the
 * Platform a chance to clear buffers and do other initialization at the beginning, and to do any
 * necessary work at the end to determine the final results.
 */
class MetalCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    MetalCalcForcesAndEnergyKernel(std::string name, const Platform& platform, MetalContext& mc) : CalcForcesAndEnergyKernel(name, platform), mc(mc) {
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
   MetalContext& mc;
};

/**
 * This kernel provides methods for setting and retrieving various state data: time, positions,
 * velocities, and forces.
 */
class MetalUpdateStateDataKernel : public CommonUpdateStateDataKernel {
public:
    MetalUpdateStateDataKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonUpdateStateDataKernel(name, platform, cc) {
    }
    /**
     * Set the positions of all particles.
     *
     * @param positions  a vector containg the particle positions
     */
    void setPositions(ContextImpl& context, const std::vector<Vec3>& positions);
    /**
     * Set the velocities of all particles.
     *
     * @param velocities  a vector containg the particle velocities
     */
    void setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities);
};

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class MetalCalcNonbondedForceKernel : public CommonCalcNonbondedForceKernel {
public:
    MetalCalcNonbondedForceKernel(std::string name, const Platform& platform, MetalContext& mc, const System& system) :
            CommonCalcNonbondedForceKernel(name, platform, mc, system), mc(mc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the NonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const NonbondedForce& force);
private:
    MetalContext& mc;
};

/**
 * This kernel is invoked by ConstantPotentialForce to calculate the forces acting on the system.
 */
class MetalCalcConstantPotentialForceKernel : public CommonCalcConstantPotentialForceKernel {
public:
    MetalCalcConstantPotentialForceKernel(std::string name, const Platform& platform, MetalContext& mc, const System& system) :
            CommonCalcConstantPotentialForceKernel(name, platform, mc, system), mc(mc) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the ConstantPotentialForce this kernel will be used for
     */
    void initialize(const System& system, const ConstantPotentialForce& force);
private:
    MetalContext& mc;
};

/**
 * This kernel is invoked by CustomCVForce to calculate the forces acting on the system and the energy of the system.
 */
class MetalCalcCustomCVForceKernel : public CommonCalcCustomCVForceKernel {
public:
    MetalCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcCustomCVForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<MetalPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

class MetalCalcATMForceKernel : public CommonCalcATMForceKernel {
public:
    MetalCalcATMForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<MetalPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

} // namespace OpenMM

#endif /*OPENMM_METALKERNELS_H_*/
