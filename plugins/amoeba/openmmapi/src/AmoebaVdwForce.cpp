/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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
#include "openmm/internal/AssertionUtilities.h"

using namespace OpenMM;
using std::string;
using std::vector;

AmoebaVdwForce::AmoebaVdwForce() : nonbondedMethod(NoCutoff), potentialFunction(Buffered147),
        sigmaCombiningRule("CUBIC-MEAN"), epsilonCombiningRule("HHG"), cutoff(1.0e+10), useDispersionCorrection(true),
        useTypes(false), alchemicalMethod(None), n(5), alpha(0.7) {
}

int AmoebaVdwForce::addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor, bool isAlchemical) {
    if (useTypes)
        throw OpenMMException("AmoebaVdwForce: must use the same version of addParticle() for all particles");
    parameters.push_back(VdwInfo(parentIndex, sigma, epsilon, -1, reductionFactor, isAlchemical));
    return parameters.size()-1;
}

int AmoebaVdwForce::addParticle(int parentIndex, int typeIndex, double reductionFactor, bool isAlchemical) {
    if (parameters.size() > 0 && !useTypes)
        throw OpenMMException("AmoebaVdwForce: must use the same version of addParticle() for all particles");
    useTypes = true;
    parameters.push_back(VdwInfo(parentIndex, 1.0, 0.0, typeIndex, reductionFactor, isAlchemical));
    return parameters.size()-1;
}

void AmoebaVdwForce::getParticleParameters(int particleIndex, int& parentIndex,
                                           double& sigma, double& epsilon, double& reductionFactor, bool& isAlchemical, int& typeIndex) const {
    ASSERT_VALID_INDEX(particleIndex, parameters);
    parentIndex     = parameters[particleIndex].parentIndex;
    sigma           = parameters[particleIndex].sigma;
    epsilon         = parameters[particleIndex].epsilon;
    reductionFactor = parameters[particleIndex].reductionFactor;
    isAlchemical    = parameters[particleIndex].isAlchemical;
    typeIndex       = parameters[particleIndex].typeIndex;
}

void AmoebaVdwForce::setParticleParameters(int particleIndex, int parentIndex,
                                           double sigma, double epsilon, double reductionFactor, bool isAlchemical, int typeIndex) {
    ASSERT_VALID_INDEX(particleIndex, parameters);
    parameters[particleIndex].parentIndex     = parentIndex;
    parameters[particleIndex].sigma           = sigma;
    parameters[particleIndex].epsilon         = epsilon;
    parameters[particleIndex].reductionFactor = reductionFactor;
    parameters[particleIndex].isAlchemical    = isAlchemical;
    parameters[particleIndex].typeIndex       = typeIndex;
}

int AmoebaVdwForce::addParticleType(double sigma, double epsilon) {
    types.push_back(ParticleTypeInfo(sigma, epsilon));
    return types.size()-1;
}

void AmoebaVdwForce::getParticleTypeParameters(int typeIndex, double& sigma, double& epsilon) const {
    ASSERT_VALID_INDEX(typeIndex, types);
    sigma = types[typeIndex].sigma;
    epsilon = types[typeIndex].epsilon;
}

void AmoebaVdwForce::setParticleTypeParameters(int typeIndex, double sigma, double epsilon) {
    ASSERT_VALID_INDEX(typeIndex, types);
    types[typeIndex].sigma = sigma;
    types[typeIndex].epsilon = epsilon;
}

int AmoebaVdwForce::addTypePair(int type1, int type2, double sigma, double epsilon) {
    pairs.push_back(TypePairInfo(type1, type2, sigma, epsilon));
    return pairs.size()-1;
}

void AmoebaVdwForce::getTypePairParameters(int pairIndex, int& type1, int& type2, double& sigma, double& epsilon) const {
    ASSERT_VALID_INDEX(pairIndex, pairs);
    type1 = pairs[pairIndex].type1;
    type2 = pairs[pairIndex].type2;
    sigma = pairs[pairIndex].sigma;
    epsilon = pairs[pairIndex].epsilon;
}

void AmoebaVdwForce::setTypePairParameters(int pairIndex, int type1, int type2, double sigma, double epsilon) {
    ASSERT_VALID_INDEX(pairIndex, pairs);
    pairs[pairIndex].type1 = type1;
    pairs[pairIndex].type2 = type2;
    pairs[pairIndex].sigma = sigma;
    pairs[pairIndex].epsilon = epsilon;
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
