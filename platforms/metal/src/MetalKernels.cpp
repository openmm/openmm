/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/src/OpenCLKernels.cpp                            *
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

#include "MetalKernels.h"
#include "MetalPlatform.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"

using namespace OpenMM;
using namespace std;

#ifndef OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
#define OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS 1
#endif

namespace {

#if OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
/** @brief Restore the original accumulator mode even if minimization or a reporter throws. */
class FloatingAccumulatorScope {
public:
    explicit FloatingAccumulatorScope(MetalContext& context) : context(context),
            previous(context.getUseFloatingPointAccumulators()), active(true) {
        try {
            context.setUseFloatingPointAccumulators(true);
        }
        catch (...) {
            context.restoreFloatingPointAccumulators(previous);
            throw;
        }
    }
    ~FloatingAccumulatorScope() noexcept {
        if (active)
            context.restoreFloatingPointAccumulators(previous);
    }
    /** @brief Restore normally while allowing GPU synchronization errors to propagate. */
    void finish() {
        context.setUseFloatingPointAccumulators(previous);
        active = false;
    }
private:
    MetalContext& context;
    bool previous, active;
};
#else
/** @brief Check Q32.32 contributions without leaving linked Context diagnostics enabled. */
class FixedPointRangeScope {
public:
    explicit FixedPointRangeScope(MetalContext& context) : context(context), active(true) {
        try {
            context.setCheckFixedPointRange(true);
        }
        catch (...) {
            context.restoreCheckFixedPointRange();
            throw;
        }
    }
    ~FixedPointRangeScope() noexcept {
        if (active)
            context.restoreCheckFixedPointRange();
    }
    /** @brief Report a recorded range error, then disable checking on normal exit. */
    void finish() {
        context.checkFixedPointRange();
        context.setCheckFixedPointRange(false);
        active = false;
    }
private:
    MetalContext& context;
    bool active;
};
#endif

} // namespace

void MetalUpdateStateDataKernel::getForces(ContextImpl& simulation, vector<Vec3>& forces) {
    if (!context.getUseFloatingPointAccumulators()) {
        CommonUpdateStateDataKernel::getForces(simulation, forces);
        return;
    }
    ContextSelector selector(context);
    vector<double> values;
    context.downloadFixedPointBuffer(context.getLongForceBuffer(), values);
    const vector<int>& order = context.getAtomIndex();
    const int count = simulation.getSystem().getNumParticles();
    const int padded = context.getPaddedNumAtoms();
    forces.resize(count);
    for (int i = 0; i < count; i++)
        forces[order[i]] = Vec3(values[i], values[i+padded], values[i+2*padded]);
}

MetalMinimizeKernel::MetalMinimizeKernel(string name, const Platform& platform, MetalContext& context) :
        CommonMinimizeKernel(name, platform, context, false, OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS != 0),
        context(context) {
}

void MetalMinimizeKernel::execute(ContextImpl& simulation, double tolerance, int maxIterations,
        MinimizationReporter* reporter) {
#if OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
    FloatingAccumulatorScope scope(context);
    CommonMinimizeKernel::execute(simulation, tolerance, maxIterations, reporter);
#else
    FixedPointRangeScope scope(context);
    try {
        CommonMinimizeKernel::execute(simulation, tolerance, maxIterations, reporter);
    }
    catch (...) {
        // Prefer a recorded GPU range error over secondary failures produced
        // by an invalid gradient. Scope destruction disables checking either way.
        context.checkFixedPointRange();
        throw;
    }
#endif
    scope.finish();
}

void MetalCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void MetalCalcForcesAndEnergyKernel::beginComputation(ContextImpl& simulation, bool includeForces, bool includeEnergy, int groups) {
    context.setForcesValid(true);
    context.clearAutoclearBuffers();
    context.updateGlobalParamValues();
    for (auto computation : context.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    context.setComputeForceCount(context.getComputeForceCount()+1);
    context.getNonbondedUtilities().prepareInteractions(groups);
    map<string, double>& derivatives = context.getEnergyParamDerivWorkspace();
    for (auto& parameter : simulation.getParameters())
        derivatives[parameter.first] = 0;
}

double MetalCalcForcesAndEnergyKernel::finishComputation(ContextImpl& simulation, bool includeForces, bool includeEnergy,
        int groups, bool& valid) {
    context.getBondedUtilities().computeInteractions(groups);
    context.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double energy = 0;
    for (auto computation : context.getPostComputations())
        energy += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    context.reduceForces();
    context.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        energy += context.reduceEnergy();
    if (!context.getForcesValid())
        valid = false;
    return energy;
}

void MetalCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    commonInitialize(system, force, false, false, true, false);
}

void MetalCalcConstantPotentialForceKernel::initialize(const System& system, const ConstantPotentialForce& force) {
    commonInitialize(system, force, false, true);
}

ComputeContext& MetalCalcCustomCVForceKernel::getInnerComputeContext(ContextImpl& innerContext) {
    return *static_cast<MetalPlatform::PlatformData*>(innerContext.getPlatformData())->computeContext;
}

ComputeContext& MetalCalcATMForceKernel::getInnerComputeContext(ContextImpl& innerContext) {
    return *static_cast<MetalPlatform::PlatformData*>(innerContext.getPlatformData())->computeContext;
}
