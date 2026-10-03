/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/include/OpenCLKernels.h                          *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALKERNELS_H_
#define OPENMM_METALKERNELS_H_

#include "MetalContext.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/common/CommonCalcNonbondedForce.h"
#include "openmm/common/CommonCalcConstantPotentialForce.h"
#include "openmm/common/CommonMinimizeKernel.h"

namespace OpenMM {

/** @brief Read forces using the currently selected accumulator format. */
class MetalUpdateStateDataKernel : public CommonUpdateStateDataKernel {
public:
    /** @brief Retain the context for floating-accumulator readback. */
    MetalUpdateStateDataKernel(std::string name, const Platform& platform, MetalContext& context) :
            CommonUpdateStateDataKernel(name, platform, context), context(context) {
    }
    /** @brief Use ordinary Common readback except while floating accumulators are selected. */
    void getForces(ContextImpl& context, std::vector<Vec3>& forces) override;
private:
    MetalContext& context;
};

/**
 * @brief Run Common minimization without CPU fallback.
 *
 * The original Context, ForceImpl objects, applied parameters, and callbacks
 * remain active. When OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS is enabled,
 * accumulator pipelines change temporarily to avoid the Q32.32 range limit.
 * Disabling it retains the ordinary fixed-point path and Common reductions.
 * That path diagnoses out-of-range individual contributions, not every possible
 * overflow caused by summing many individually representable contributions.
 */
class MetalMinimizeKernel : public CommonMinimizeKernel {
public:
    /** @brief Select optional floating accumulators and always prohibit CPU fallback. */
    MetalMinimizeKernel(std::string name, const Platform& platform, MetalContext& context);
    /** @brief Minimize on the original Context and restore its accumulator mode on every exit. */
    void execute(ContextImpl& context, double tolerance, int maxIterations, MinimizationReporter* reporter) override;
private:
    MetalContext& context;
};

/** @brief Clear buffers, schedule shared interactions, and finish a force evaluation. */
class MetalCalcForcesAndEnergyKernel : public CalcForcesAndEnergyKernel {
public:
    /** @brief Borrow the compute context used by this kernel. */
    MetalCalcForcesAndEnergyKernel(std::string name, const Platform& platform, MetalContext& context) :
            CalcForcesAndEnergyKernel(name, platform), context(context) {
    }
    /** @brief Perform no additional initialization; the context owns shared utilities. */
    void initialize(const System& system) override;
    /** @brief Clear accumulators, update parameters, and prepare nonbonded interactions. */
    void beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) override;
    /**
     * @brief Evaluate shared interactions and finish the force and energy reductions.
     * @param valid Set false when a neighbor-list overflow requires another evaluation.
     * @return Potential energy accumulated in the shared buffers and post-computations.
     */
    double finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy,
            int groups, bool& valid) override;
private:
    MetalContext& context;
};

/** @brief Configure the shared NonbondedForce implementation for GPU-only execution. */
class MetalCalcNonbondedForceKernel : public CommonCalcNonbondedForceKernel {
public:
    /** @brief Forward the simulation resources to the Common Compute implementation. */
    MetalCalcNonbondedForceKernel(std::string name, const Platform& platform, MetalContext& context, const System& system) :
            CommonCalcNonbondedForceKernel(name, platform, context, system) {
    }
    /** @brief Use the main GPU queue and independently selectable PME/LJPME grid formats; never CPU PME. */
    void initialize(const System& system, const NonbondedForce& force) override;
};

/** @brief Configure the shared ConstantPotentialForce implementation for the GPU. */
class MetalCalcConstantPotentialForceKernel : public CommonCalcConstantPotentialForceKernel {
public:
    /** @brief Forward the simulation resources to the Common Compute implementation. */
    MetalCalcConstantPotentialForceKernel(std::string name, const Platform& platform, MetalContext& context, const System& system) :
            CommonCalcConstantPotentialForceKernel(name, platform, context, system) {
    }
    /** @brief Use GPU execution with the independently selectable ConstantPotential grid format. */
    void initialize(const System& system, const ConstantPotentialForce& force) override;
};

/** @brief Connect CustomCVForce's shared implementation to its inner Metal Context. */
class MetalCalcCustomCVForceKernel : public CommonCalcCustomCVForceKernel {
public:
    /** @brief Forward the compute context to the shared implementation. */
    MetalCalcCustomCVForceKernel(std::string name, const Platform& platform, ComputeContext& context) :
            CommonCalcCustomCVForceKernel(name, platform, context) {
    }
    /** @return The compute context owned by the linked inner simulation Context. */
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) override;
};

/** @brief Connect ATMForce's shared implementation to its inner Metal Context. */
class MetalCalcATMForceKernel : public CommonCalcATMForceKernel {
public:
    /** @brief Forward the compute context to the shared implementation. */
    MetalCalcATMForceKernel(std::string name, const Platform& platform, ComputeContext& context) :
            CommonCalcATMForceKernel(name, platform, context) {
    }
    /** @return The compute context owned by the linked inner simulation Context. */
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) override;
};

} // namespace OpenMM

#endif // OPENMM_METALKERNELS_H_
