/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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

#include "CudaKernels.h"
#include "CudaForceInfo.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/common/ContextSelector.h"
#include "CommonKernelSources.h"
#include "CudaBondedUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaIntegrationUtilities.h"
#include "CudaNonbondedUtilities.h"
#include "CudaKernelSources.h"
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
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

void CudaCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void CudaCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.setForcesValid(true);
    ContextSelector selector(cu);
    cu.clearAutoclearBuffers();
    for (auto computation : cu.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    cu.setComputeForceCount(cu.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = cu.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
}

double CudaCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
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

void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    bool usePmeQueue = (!cu.getPlatformData().disablePmeStream && !cu.getPlatformData().useCpuPme);
    bool useFixedPointChargeSpreading = cu.getUseDoublePrecision() || cu.getPlatformData().deterministicForces;
    commonInitialize(system, force, usePmeQueue, false, useFixedPointChargeSpreading, cu.getPlatformData().useCpuPme);
}
