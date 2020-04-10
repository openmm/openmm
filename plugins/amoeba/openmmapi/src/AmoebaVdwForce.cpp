/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/AmoebaVdwForce.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"

using namespace OpenMM;
using std::string;
using std::vector;

AmoebaVdwForce::AmoebaVdwForce() : nonbondedMethod(NoCutoff), sigmaCombiningRule("CUBIC-MEAN"), epsilonCombiningRule("HHG"), cutoff(1.0e+10), useDispersionCorrection(true), alchemicalMethod(None), n(5), alpha(0.7), usesVdwpr(false), usesLJ(false) {
}

int AmoebaVdwForce::addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor, bool isAlchemical) {
    parameters.push_back(VdwInfo(parentIndex, sigma, epsilon, reductionFactor, isAlchemical));
    return parameters.size()-1;
}

void AmoebaVdwForce::getParticleParameters(int particleIndex, int& parentIndex,
                                           double& sigma, double& epsilon, double& reductionFactor, bool& isAlchemical) const {
    parentIndex     = parameters[particleIndex].parentIndex;
    sigma           = parameters[particleIndex].sigma;
    epsilon         = parameters[particleIndex].epsilon;
    reductionFactor = parameters[particleIndex].reductionFactor;
    isAlchemical    = parameters[particleIndex].isAlchemical;
}

void AmoebaVdwForce::setParticleParameters(int particleIndex, int parentIndex,
                                           double sigma, double epsilon, double reductionFactor, bool isAlchemical) {
    parameters[particleIndex].parentIndex     = parentIndex;
    parameters[particleIndex].sigma           = sigma;
    parameters[particleIndex].epsilon         = epsilon;
    parameters[particleIndex].reductionFactor = reductionFactor;
    parameters[particleIndex].isAlchemical    = isAlchemical;
}

void AmoebaVdwForce::setSigmaCombiningRule(const std::string& inputSigmaCombiningRule) {
    sigmaCombiningRule = inputSigmaCombiningRule;
}

const std::string& AmoebaVdwForce::getSigmaCombiningRule() const {
    return sigmaCombiningRule;
}

void AmoebaVdwForce::setEpsilonCombiningRule(const std::string& inputEpsilonCombiningRule) {
    epsilonCombiningRule = inputEpsilonCombiningRule;
}

const std::string& AmoebaVdwForce::getEpsilonCombiningRule() const {
    return epsilonCombiningRule;
}

int AmoebaVdwForce::getNumCondensedTypes() const {
    return (int) sigEpsTable.size();
}

int AmoebaVdwForce::addCondensedType(int type) {
    condensedTypes.push_back(type);
    return condensedTypes.size() - 1;
}

void AmoebaVdwForce::getCondensedType(int particleIndex, int& type) const {
    type = condensedTypes[particleIndex];
}

void AmoebaVdwForce::setPairSigmaEpsilon(int type1, int type2, double sig, double eps) {
    if (!usesVdwpr) {
        usesVdwpr = true;
        std::set<int> uniqueTypes(condensedTypes.begin(), condensedTypes.end());
        size_t nunique = uniqueTypes.size();
        sigEpsTable.resize(nunique);
        for (size_t i = 0; i < nunique; ++i) {
            sigEpsTable[i].resize(nunique);
        }
    }

    PairSigEps& pr1 = sigEpsTable[type1][type2];
    pr1.sig = sig;
    pr1.eps = eps;
    if (type1 != type2) {
        PairSigEps& pr2 = sigEpsTable[type2][type1];
        pr2.sig = sig;
        pr2.eps = eps;
    }
}

void AmoebaVdwForce::getPairSigmaEpsilon(int type1, int type2, double& sig, double& eps) const {
    const PairSigEps& pr = sigEpsTable[type1][type2];
    sig = pr.sig;
    eps = pr.eps;
}

bool AmoebaVdwForce::getUsePairwiseVdw() const {
    return usesVdwpr;
}

    
void AmoebaVdwForce::setUseLennardJones() {
    usesLJ = true;
}


bool AmoebaVdwForce::getUseLennardJones() const {
    return usesLJ;
}


void AmoebaVdwForce::setParticleExclusions(int particleIndex, const std::vector< int >& inputExclusions) {

   if (exclusions.size() < parameters.size()) {
       exclusions.resize(parameters.size());
   }
   if (static_cast<int>(exclusions.size()) < particleIndex) {
       exclusions.resize(particleIndex + 10);
   }
   for (unsigned int ii = 0; ii < inputExclusions.size(); ii++) {
       exclusions[particleIndex].push_back(inputExclusions[ii]);
   }
}

void AmoebaVdwForce::getParticleExclusions(int particleIndex, std::vector< int >& outputExclusions) const {

   if (particleIndex < static_cast<int>(exclusions.size())) {
       outputExclusions.resize(exclusions[particleIndex].size());
       for (unsigned int ii = 0; ii < exclusions[particleIndex].size(); ii++) {
           outputExclusions[ii] = exclusions[particleIndex][ii];
       }
   }

}

double AmoebaVdwForce::getCutoffDistance() const {
    return cutoff;
}

void AmoebaVdwForce::setCutoffDistance(double inputCutoff) {
    cutoff = inputCutoff;
}

void AmoebaVdwForce::setCutoff(double inputCutoff) {
    setCutoffDistance(inputCutoff);
}

double AmoebaVdwForce::getCutoff() const {
    return getCutoffDistance();
}

AmoebaVdwForce::NonbondedMethod AmoebaVdwForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void AmoebaVdwForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 1)
        throw OpenMMException("AmoebaVdwForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

AmoebaVdwForce::AlchemicalMethod AmoebaVdwForce::getAlchemicalMethod() const {
    return alchemicalMethod;
}

void AmoebaVdwForce::setAlchemicalMethod(AlchemicalMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("AmoebaVdwForce: Illegal value for alchemical method");
    alchemicalMethod = method;
}

void AmoebaVdwForce::setSoftcorePower(int power) {
    n = power;
}

int AmoebaVdwForce::getSoftcorePower() const {
    return n;
}

void AmoebaVdwForce::setSoftcoreAlpha(double a) {
    alpha = a;
}

double AmoebaVdwForce::getSoftcoreAlpha() const {
    return alpha;
}


ForceImpl* AmoebaVdwForce::createImpl() const {
    return new AmoebaVdwForceImpl(*this);
}

void AmoebaVdwForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaVdwForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
