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

#include "HipKernels.h"
#include "HipForceInfo.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/common/ContextSelector.h"
#include "CommonKernelSources.h"
#include "HipBondedUtilities.h"
#include "HipExpressionUtilities.h"
#include "HipIntegrationUtilities.h"
#include "HipNonbondedUtilities.h"
#include "HipKernelSources.h"
#include "HipQueue.h"
#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMUtilities.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <assert.h>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result, prefix) \
    if (result != hipSuccess) { \
        std::stringstream m; \
        m<<prefix<<": "<<HipContext::getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

void HipCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void HipCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.setForcesValid(true);
    ContextSelector selector(cu);
    cu.clearAutoclearBuffers();
    for (auto computation : cu.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    HipNonbondedUtilities& nb = cu.getNonbondedUtilities();
    cu.setComputeForceCount(cu.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = cu.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
}

double HipCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    ContextSelector selector(cu);
    cu.getBondedUtilities().computeInteractions(groups);
    cu.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double sum = 0.0;
    for (auto computation : cu.getPostComputations())
        sum += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    cu.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        sum += cu.reduceEnergy();
    if (!cu.getForcesValid())
        valid = false;
    return sum;
}

void HipCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    bool usePmeQueue = (!cu.getPlatformData().disablePmeStream && !cu.getPlatformData().useCpuPme);
    bool useFixedPointChargeSpreading = cu.getUseDoublePrecision() || !cu.getSupportsHardwareFloatGlobalAtomicAdd() || cu.getPlatformData().deterministicForces;
    commonInitialize(system, force, usePmeQueue, false, useFixedPointChargeSpreading, cu.getPlatformData().useCpuPme);
}
