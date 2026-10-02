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

#include "MetalKernels.h"
#include "openmm/internal/ContextImpl.h"

using namespace OpenMM;
using namespace std;

static void getMetalPmeParameters(MetalContext& mc, bool& usePmeQueue, bool& useFixedPointChargeSpreading) {
    usePmeQueue = (!mc.getPlatformData().disablePmeStream && !mc.getPlatformData().useCpuPme);
    useFixedPointChargeSpreading = mc.getPlatformData().deterministicForces;
}

void MetalCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void MetalCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    mc.setForcesValid(true);
    mc.clearAutoclearBuffers();
    mc.updateGlobalParamValues();
    for (auto computation : mc.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    MetalNonbondedUtilities& nb = mc.getNonbondedUtilities();
    mc.setComputeForceCount(mc.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = mc.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
    mc.flushQueue();
}

double MetalCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    mc.getBondedUtilities().computeInteractions(groups);
    mc.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double sum = 0.0;
    for (auto computation : mc.getPostComputations())
        sum += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    mc.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        sum += mc.reduceEnergy();
    if (!mc.getForcesValid())
        valid = false;
    mc.flushQueue();
    return sum;
}

void MetalCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    bool usePmeQueue, useFixedPointChargeSpreading;
    getMetalPmeParameters(mc, usePmeQueue, useFixedPointChargeSpreading);
    commonInitialize(system, force, usePmeQueue, false, useFixedPointChargeSpreading, mc.getPlatformData().useCpuPme);
}

void MetalCalcConstantPotentialForceKernel::initialize(const System& system, const ConstantPotentialForce& force) {
    bool usePmeQueue, useFixedPointChargeSpreading;
    getMetalPmeParameters(mc, usePmeQueue, useFixedPointChargeSpreading);
    commonInitialize(system, force, false, useFixedPointChargeSpreading);
}
