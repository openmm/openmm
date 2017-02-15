/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include <cmath>

using namespace OpenMM;
using std::string;
using std::vector;

AmoebaVdwForce::AmoebaVdwForce()
    : nonbondedMethod(NoCutoff)
    , sigmaCombiningRule("CUBIC-MEAN")
    , epsilonCombiningRule("HHG")
    , functionalForm("BUFFERED-14-7")
    , cutoff(1.0e+10)
    , useDispersionCorrection(false)
    , numVdwprTypes(0) {}

void AmoebaVdwForce::setParticleParameters(int particleIndex,
                                           int parentIndex,
                                           int vdwprType,
                                           double sigma,
                                           double epsilon,
                                           double reductionFactor,
                                           double lambda) {
    parameters[particleIndex].parentIndex = parentIndex;
    parameters[particleIndex].vdwprType = vdwprType;
    parameters[particleIndex].sigma = sigma;
    parameters[particleIndex].epsilon = epsilon;
    parameters[particleIndex].reductionFactor = reductionFactor;
    parameters[particleIndex].lambda = lambda;
}

void AmoebaVdwForce::getParticleParameters(int particleIndex,
                                           int& parentIndex,
                                           int& vdwprType,
                                           double& sigma,
                                           double& epsilon,
                                           double& reductionFactor,
                                           double& lambda) const {
    parentIndex = parameters[particleIndex].parentIndex;
    vdwprType = parameters[particleIndex].vdwprType;
    sigma = parameters[particleIndex].sigma;
    epsilon = parameters[particleIndex].epsilon;
    reductionFactor = parameters[particleIndex].reductionFactor;
    lambda = parameters[particleIndex].lambda;
}

int AmoebaVdwForce::addParticle(int parentIndex,
                                int vdwprType,
                                double sigma,
                                double epsilon,
                                double reductionFactor,
                                double lambda) {
    parameters.push_back(VdwInfo(parentIndex,
                                 vdwprType,
                                 sigma,
                                 epsilon,
                                 reductionFactor,
                                 lambda));
    return static_cast<int>(parameters.size() - 1);
}

void AmoebaVdwForce::computeCombinedSigmaEpsilon() {
    std::map<int, int> m_typeMap; // m_typeMap.at(old) = new; typeMap.at(new) = old;
    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        m_typeMap[parameters[i].vdwprType] = -1;
    }

    numVdwprTypes = static_cast<int>(m_typeMap.size());
    typeMap.resize(numVdwprTypes);
    arguments.resize(numVdwprTypes * numVdwprTypes);

    std::map<int, int>::iterator it = m_typeMap.begin();
    for (int i = 0; i < numVdwprTypes; ++i) {
        typeMap[i] = it->first;
        it->second = i;
        ++it;
    }

    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        int tp = parameters[i].vdwprType; // tp = old
        parameters[i].vdwprType = m_typeMap.at(tp);
    }

    std::vector<double> sigmaByType(numVdwprTypes);
    std::vector<double> epsilonByType(numVdwprTypes);
    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        int tp = parameters[i].vdwprType; // tp = new
        sigmaByType[tp] = parameters[i].sigma;
        epsilonByType[tp] = parameters[i].epsilon;
    }

    // compute combined sigma and epsilon

    for (int i = 0; i < numVdwprTypes; ++i) {
        int p = i * numVdwprTypes + i;
        VdwprInfo para(sigmaByType[i] * 2, epsilonByType[i]);
        arguments[p] = para;
    }

    for (int i = 0; i < numVdwprTypes; ++i) {
        for (int j = i + 1; j < numVdwprTypes; ++j) {
            int p = i * numVdwprTypes + j;

            double sigma1 = sigmaByType[i];
            double sigma2 = sigmaByType[j];
            double epsilon1 = epsilonByType[i];
            double epsilon2 = epsilonByType[j];

            VdwprInfo para;

            if (sigmaCombiningRule == "ARITHMETIC") {
                para.combinedSigma = sigma1 + sigma2;
            }

            if (sigmaCombiningRule == "GEOMETRIC") {
                para.combinedSigma = 2 * std::sqrt(sigma1 * sigma2);
            }

            if (sigmaCombiningRule == "CUBIC-MEAN") {
                double sigma12 = sigma1 * sigma1;
                double sigma22 = sigma2 * sigma2;
                double sigmasum = sigma12 + sigma22;
                double sigmacub = sigma12 * sigma1 + sigma22 * sigma2;
                para.combinedSigma = (sigmasum == 0.0 ? 0.0 : 2 * sigmacub / sigmasum);
            }

            if (epsilonCombiningRule == "ARITHMETIC") {
                para.combinedEpsilon = 0.5 * (epsilon1 + epsilon2);
            }

            if (epsilonCombiningRule == "GEOMETRIC") {
                para.combinedEpsilon = std::sqrt(epsilon1 * epsilon2);
            }

            if (epsilonCombiningRule == "HARMONIC") {
                double epssum = epsilon1 + epsilon2;
                para.combinedEpsilon = (epssum == 0.0 ? 0.0 : 2 * epsilon1 * epsilon2 / epssum);
            }

            if (epsilonCombiningRule == "HHG") {
                double eps_s = std::sqrt(epsilon1) + std::sqrt(epsilon2);
                para.combinedEpsilon = (eps_s == 0.0 ? 0.0 : 4 * epsilon1 * epsilon2 / (eps_s * eps_s));
            }

            arguments[p] = para;
            p = j * numVdwprTypes + i;
            arguments[p] = para;
        }
    }
}

void AmoebaVdwForce::setNumVdwprTypes(int newNum) {
    numVdwprTypes = newNum;
}

int AmoebaVdwForce::getNewVdwprType(int oldType) const {
    for (int i = 0; i < numVdwprTypes; ++i) {
        if (typeMap[i] == oldType) {
            return i;
        }
    }
    return -1;
}

int AmoebaVdwForce::getOldVdwprType(int newType) const {
    return typeMap[newType];
}

void AmoebaVdwForce::setOldVdwprType(int newType, int oldType) {
    typeMap[newType] = oldType;
}

void AmoebaVdwForce::resize(int newSize) {
    arguments.resize(newSize * newSize);
    typeMap.resize(newSize);
}

void AmoebaVdwForce::setVdwprParametersByOldTypes(int oldtype1, int oldtype2, double combinedSigma, double combinedEpsilon) {
    VdwprInfo para = VdwprInfo(combinedSigma, combinedEpsilon);
    int ntype1 = getNewVdwprType(oldtype1);
    int ntype2 = getNewVdwprType(oldtype2);
    if (ntype1 == -1 || ntype2 == -1) {
        fprintf(stdout, "\n VDWPR (%d, %d) are ignored.\n", oldtype1, oldtype2);
        return;
    }
    int v1 = ntype1 * numVdwprTypes + ntype2;
    int v2 = ntype2 * numVdwprTypes + ntype1;
    arguments[v1] = para;
    arguments[v2] = para;
}

int AmoebaVdwForce::addVdwprByOldTypes(int oldtype1, int oldtype2, double combinedSigma, double combinedEpsilon) {
    this->setVdwprParametersByOldTypes(oldtype1, oldtype2, combinedSigma, combinedEpsilon);
    return static_cast<int>(arguments.size() - 1);
}

void AmoebaVdwForce::setVdwprParameters(int ntype1, int ntype2, double combinedSigma, double combinedEpsilon) {
    VdwprInfo para = VdwprInfo(combinedSigma, combinedEpsilon);
    int v1 = ntype1 * numVdwprTypes + ntype2;
    int v2 = ntype2 * numVdwprTypes + ntype1;
    arguments[v1] = para;
    arguments[v2] = para;
}

void AmoebaVdwForce::getVdwprParameters(int ntype1, int ntype2, double& combinedSigma, double& combinedEpsilon) const {
    int v1 = ntype1 * numVdwprTypes + ntype2;
    VdwprInfo para = arguments.at(v1);
    combinedSigma = para.combinedSigma;
    combinedEpsilon = para.combinedEpsilon;
}

int AmoebaVdwForce::addVdwpr(int ntype1, int ntype2, double combinedSigma, double combinedEpsilon) {
    this->setVdwprParameters(ntype1, ntype2, combinedSigma, combinedEpsilon);
    return static_cast<int>(arguments.size() - 1);
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

void AmoebaVdwForce::setFunctionalForm(const std::string& inputFuncForm) {
    functionalForm = inputFuncForm;
}

const std::string& AmoebaVdwForce::getFunctionalForm() const {
    return functionalForm;
}

void AmoebaVdwForce::setParticleExclusions(int particleIndex, const std::vector<int>& inputExclusions) {
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

void AmoebaVdwForce::getParticleExclusions(int particleIndex, std::vector<int>& outputExclusions) const {
    if (particleIndex < static_cast<int>(exclusions.size())) {
        outputExclusions.resize(exclusions[particleIndex].size());
        for (unsigned int ii = 0; ii < exclusions[particleIndex].size(); ii++) {
            outputExclusions[ii] = exclusions[particleIndex][ii];
        }
    }
}

void AmoebaVdwForce::setCutoff(double inputCutoff) {
    cutoff = inputCutoff;
}

double AmoebaVdwForce::getCutoff() const {
    return cutoff;
}

AmoebaVdwForce::NonbondedMethod AmoebaVdwForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void AmoebaVdwForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

ForceImpl* AmoebaVdwForce::createImpl() const {
    return new AmoebaVdwForceImpl(*this);
}

void AmoebaVdwForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaVdwForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
