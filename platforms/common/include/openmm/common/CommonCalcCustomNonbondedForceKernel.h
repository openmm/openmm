#ifndef OPENMM_COMMONCALCCUSTOMNONBONDEDFORCEKERNEL_H_
#define OPENMM_COMMONCALCCUSTOMNONBONDEDFORCEKERNEL_H_

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

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeParameterSet.h"
#include "openmm/Platform.h"
#include "openmm/kernels.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by CustomNonbondedForce to calculate the forces acting on the system.
 */
class CommonCalcCustomNonbondedForceKernel : public CalcCustomNonbondedForceKernel {
public:
    CommonCalcCustomNonbondedForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomNonbondedForceKernel(name, platform),
            cc(cc), params(NULL), computedValues(NULL), forceCopy(NULL), system(system), hasInitializedKernel(false) {
    }
    ~CommonCalcCustomNonbondedForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const CustomNonbondedForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the CustomNonbondedForce to copy the parameters from
     * @param firstParticle  the index of the first particle whose parameters might have changed
     * @param lastParticle   the index of the last particle whose parameters might have changed
     */
    void copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle);
private:
    class ForceInfo;
    class LongRangePostComputation;
    class LongRangeTask;
    void initInteractionGroups(const CustomNonbondedForce& force, const std::string& interactionSource, const std::vector<std::string>& tableTypes);
    ComputeContext& cc;
    ForceInfo* info;
    ComputeParameterSet* params;
    ComputeParameterSet* computedValues;
    ComputeArray globals, interactionGroupData, filteredGroupData, numGroupTiles;
    ComputeKernel interactionGroupKernel, prepareNeighborListKernel, buildNeighborListKernel, computedValuesKernel;
    std::vector<void*> interactionGroupArgs;
    std::vector<std::string> globalParamNames;
    std::vector<float> globalParamValues;
    std::vector<ComputeArray> tabulatedFunctionArrays;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    std::vector<std::string> paramNames, computedValueNames;
    std::vector<ComputeParameterInfo> paramBuffers, computedValueBuffers;
    double longRangeCoefficient;
    std::vector<double> longRangeCoefficientDerivs;
    bool hasInitializedLongRangeCorrection, hasInitializedKernel, hasParamDerivs, useNeighborList;
    int numGroupThreadBlocks;
    CustomNonbondedForce* forceCopy;
    CustomNonbondedForceImpl::LongRangeCorrectionData longRangeCorrectionData;
    const System& system;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONCALCCUSTOMNONBONDEDFORCEKERNEL_H_*/
