#ifndef OPENMM_CUSTOMNONBONDEDFORCEIMPL_H_
#define OPENMM_CUSTOMNONBONDEDFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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

#include "ForceImpl.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/Kernel.h"
#include "openmm/internal/ThreadPool.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CompiledVectorExpression.h"
#include <utility>
#include <map>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of CustomNonbondedForce.
 */

class OPENMM_EXPORT CustomNonbondedForceImpl : public ForceImpl {
public:
    class LongRangeCorrectionData;
    CustomNonbondedForceImpl(const CustomNonbondedForce& owner);
    ~CustomNonbondedForceImpl();
    void initialize(ContextImpl& context);
    const CustomNonbondedForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle);
    /**
     * Prepare for computing the long range correction.  This function pre-computes anything
     * that depends only on the Force (such as particle parameters) but not on information in
     * the Context (such as global parameters).  This allows the coefficient to be updated
     * more quickly when global parameters change.
     */
    static LongRangeCorrectionData prepareLongRangeCorrection(const CustomNonbondedForce& force, int numThreads);
    /**
     * Compute the coefficient which, when divided by the periodic box volume, gives the
     * long range correction to the energy.  If the Force computes parameter derivatives,
     * also compute the corresponding derivatives of the correction.
     */
    static void calcLongRangeCorrection(const CustomNonbondedForce& force, LongRangeCorrectionData& data, const Context& context, double& coefficient, std::vector<double>& derivatives, ThreadPool& threads);
private:
    static double integrateInteraction(Lepton::CompiledVectorExpression& expression, const std::vector<double>& params1, const std::vector<double>& params2,
            const std::vector<double>& computedValues1, const std::vector<double>& computedValues2, const CustomNonbondedForce& force, const Context& context,
            const std::vector<std::string>& paramNames, const std::vector<std::string>& computedValueNames);
    const CustomNonbondedForce& owner;
    Kernel kernel;
};

class CustomNonbondedForceImpl::LongRangeCorrectionData {
public:
    CustomNonbondedForce::NonbondedMethod method;
    std::vector<std::vector<double> > classes;
    std::vector<std::string> paramNames, computedValueNames;
    std::map<std::pair<int, int>, long long int> interactionCount;
    std::vector<Lepton::CompiledVectorExpression> energyExpression;
    std::vector<std::vector<Lepton::CompiledVectorExpression> > derivExpressions;
    std::vector<Lepton::CompiledExpression> computedValueExpressions;
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMNONBONDEDFORCEIMPL_H_*/
