#ifndef OPENMM_COMMONCALCCUSTOMGBFORCEKERNEL_H_
#define OPENMM_COMMONCALCCUSTOMGBFORCEKERNEL_H_

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
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by CustomGBForce to calculate the forces acting on the system.
 */
class CommonCalcCustomGBForceKernel : public CalcCustomGBForceKernel {
public:
    CommonCalcCustomGBForceKernel(std::string name, const Platform& platform, ComputeContext& cc, const System& system) : CalcCustomGBForceKernel(name, platform),
            hasInitializedKernels(false), cc(cc), params(NULL), computedValues(NULL), energyDerivs(NULL), energyDerivChain(NULL), system(system) {
    }
    ~CommonCalcCustomGBForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the CustomGBForce this kernel will be used for
     */
    void initialize(const System& system, const CustomGBForce& force);
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
     * @param force      the CustomGBForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const CustomGBForce& force);
private:
    class ForceInfo;
    double cutoff;
    bool hasInitializedKernels, needGlobalParams, needParameterGradient, needEnergyParamDerivs;
    int maxTiles, numComputedValues;
    ComputeContext& cc;
    ForceInfo* info;
    ComputeParameterSet* params;
    ComputeParameterSet* computedValues;
    ComputeParameterSet* energyDerivs;
    ComputeParameterSet* energyDerivChain;
    std::vector<ComputeParameterSet*> dValuedParam;
    std::vector<ComputeArray> dValue0dParam;
    ComputeArray longEnergyDerivs, valueBuffers;
    std::vector<ComputeArray> tabulatedFunctionArrays;
    std::map<std::string, int> tabulatedFunctionUpdateCount;
    std::vector<bool> pairValueUsesParam, pairEnergyUsesParam, pairEnergyUsesValue;
    const System& system;
    ComputeKernel pairValueKernel, perParticleValueKernel, pairEnergyKernel, perParticleEnergyKernel, gradientChainRuleKernel;
    std::string pairValueSrc, pairEnergySrc;
    std::map<std::string, std::string> pairValueDefines, pairEnergyDefines;
};

} // namespace OpenMM

#endif /*OPENMM_COMMONCALCCUSTOMGBFORCEKERNEL_H_*/
