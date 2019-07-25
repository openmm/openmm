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
#include "openmm/AmoebaCTForce.h"
#include "openmm/internal/AmoebaCTForceImpl.h"
#include <cmath>

using namespace OpenMM;
using std::string;
using std::vector;

AmoebaCTForce::AmoebaCTForce()
    : nonbondedMethod(NoCutoff)
    , apreCombiningRule("GEOMETRIC")
    , bexpCombiningRule("ARITHMETIC")
    , cutoff(1.0e+10)
    , numCTprTypes(0) {}

void AmoebaCTForce::setParticleParameters(int particleIndex,
                                          int CTprType,
                                          double apre,
                                          double bexp) {
    parameters[particleIndex].CTprType = CTprType;
    parameters[particleIndex].apre = apre;
    parameters[particleIndex].bexp = bexp;
}

void AmoebaCTForce::getParticleParameters(int particleIndex,
                                          int& CTprType,
                                          double& apre,
                                          double& bexp) const {
    CTprType = parameters[particleIndex].CTprType;
    apre = parameters[particleIndex].apre;
    bexp = parameters[particleIndex].bexp;
}

int AmoebaCTForce::addParticle(int CTprType,
                               double apre,
                               double bexp) {
    parameters.push_back(CTInfo(CTprType,
                                apre,
                                bexp));
    return static_cast<int>(parameters.size() - 1);
}

void AmoebaCTForce::computeCombinedApreBexp() {
    std::map<int, int> m_typeMap; // m_typeMap.at(old) = new; typeMap.at(new) = old;
    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        m_typeMap[parameters[i].CTprType] = -1;
    }

    numCTprTypes = static_cast<int>(m_typeMap.size());
    typeMap.resize(numCTprTypes);
    arguments.resize(numCTprTypes * numCTprTypes);

    std::map<int, int>::iterator it = m_typeMap.begin();
    for (int i = 0; i < numCTprTypes; ++i) {
        typeMap[i] = it->first;
        it->second = i;
        ++it;
    }

    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        int tp = parameters[i].CTprType; // tp = old
        parameters[i].CTprType = m_typeMap.at(tp);
    }

    std::vector<double> apreByType(numCTprTypes);
    std::vector<double> bexpByType(numCTprTypes);
    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        int tp = parameters[i].CTprType; // tp = new
        apreByType[tp] = parameters[i].apre;
        bexpByType[tp] = parameters[i].bexp;

        //fprintf(stdout, "\n apre checkpoint1 %f\n", parameters[i].apre); //PASSED
        //fprintf(stdout, "\n bexp %f\n", parameters[i].bexp);//PASSED
    }

    // compute combined apre and bexp

    for (int i = 0; i < numCTprTypes; ++i) {
        int p = i * numCTprTypes + i;
        CTprInfo para(apreByType[i], bexpByType[i]);
        arguments[p] = para;
    }

    for (int i = 0; i < numCTprTypes; ++i) {
        for (int j = i + 1; j < numCTprTypes; ++j) {
            int p = i * numCTprTypes + j;

            double apre1 = apreByType[i];
            double apre2 = apreByType[j];
            double bexp1 = bexpByType[i];
            double bexp2 = bexpByType[j];
            CTprInfo para;

            if (apreCombiningRule == "ARITHMETIC") {
                para.combinedApre = 0.5*(apre1 + apre2);
            }
            
            if (apreCombiningRule == "GEOMETRIC") {
                para.combinedApre = std::sqrt(apre1 * apre2);
            }


            if (bexpCombiningRule == "ARITHMETIC") {
                para.combinedBexp = 0.5*(bexp1 + bexp2);
            }

            if (bexpCombiningRule == "GEOMETRIC") {
                para.combinedBexp = std::sqrt(bexp1 * bexp2);
            }
            
            //fprintf(stdout, "\n combinedApre %f %f %f\n", para.combinedApre, apre1, apre2); //PASSED
            //fprintf(stdout, "\n combinedApre %f\n", para.combinedApre); //PASSED
            //fprintf(stdout, "\n combinedBexp %f\n", para.combinedBexp); //PASSED

            arguments[p] = para;
            p = j * numCTprTypes + i;
            arguments[p] = para;
        }
    }
}

void AmoebaCTForce::setNumCTprTypes(int newNum) {
    numCTprTypes = newNum;
}

int AmoebaCTForce::getNewCTprType(int oldType) const {
    for (int i = 0; i < numCTprTypes; ++i) {
        if (typeMap[i] == oldType) {
            return i;
        }
    }
    return -1;
}

int AmoebaCTForce::getOldCTprType(int newType) const {
    return typeMap[newType];
}

void AmoebaCTForce::setOldCTprType(int newType, int oldType) {
    typeMap[newType] = oldType;
}

void AmoebaCTForce::resize(int newSize) {
    arguments.resize(newSize * newSize);
    typeMap.resize(newSize);
}

void AmoebaCTForce::setCTprParametersByOldTypes(int oldtype1, int oldtype2, double combinedApre, double combinedBexp) {
    CTprInfo para = CTprInfo(combinedApre, combinedBexp);
    int ntype1 = getNewCTprType(oldtype1);
    int ntype2 = getNewCTprType(oldtype2);
    if (ntype1 == -1 || ntype2 == -1) {
        fprintf(stdout, "\n CTPR (%d, %d) are ignored.\n", oldtype1, oldtype2);
        return;
    }
    int ct1 = ntype1 * numCTprTypes + ntype2;
    int ct2 = ntype2 * numCTprTypes + ntype1;
    arguments[ct1] = para;
    arguments[ct2] = para;
}

int AmoebaCTForce::addCTprByOldTypes(int oldtype1, int oldtype2, double combinedApre, double combinedBexp) {
    this->setCTprParametersByOldTypes(oldtype1, oldtype2, combinedApre, combinedBexp);
    return static_cast<int>(arguments.size() - 1);
}

void AmoebaCTForce::setCTprParameters(int ntype1, int ntype2, double combinedApre, double combinedBexp) {
    CTprInfo para = CTprInfo(combinedApre, combinedBexp);
    int ct1 = ntype1 * numCTprTypes + ntype2;
    int ct2 = ntype2 * numCTprTypes + ntype1;
    arguments[ct1] = para;
    arguments[ct2] = para;
}

void AmoebaCTForce::getCTprParameters(int ntype1, int ntype2, double& combinedApre, double& combinedBexp) const {
    int ct1 = ntype1 * numCTprTypes + ntype2;
    CTprInfo para = arguments.at(ct1);
    combinedApre = para.combinedApre;
    combinedBexp = para.combinedBexp;
}

int AmoebaCTForce::addCTpr(int ntype1, int ntype2, double combinedApre, double combinedBexp) {
    this->setCTprParameters(ntype1, ntype2, combinedApre, combinedBexp);
    return static_cast<int>(arguments.size() - 1);
}

void AmoebaCTForce::setApreCombiningRule(const std::string& inputApreCombiningRule) {
    apreCombiningRule = inputApreCombiningRule;
}

const std::string& AmoebaCTForce::getApreCombiningRule() const {
    return apreCombiningRule;
}

void AmoebaCTForce::setBexpCombiningRule(const std::string& inputBexpCombiningRule) {
    bexpCombiningRule = inputBexpCombiningRule;
}

const std::string& AmoebaCTForce::getBexpCombiningRule() const {
    return bexpCombiningRule;
}

void AmoebaCTForce::setParticleExclusions(int particleIndex, const std::vector<int>& inputExclusions) {
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

void AmoebaCTForce::getParticleExclusions(int particleIndex, std::vector<int>& outputExclusions) const {
    if (particleIndex < static_cast<int>(exclusions.size())) {
        outputExclusions.resize(exclusions[particleIndex].size());
        for (unsigned int ii = 0; ii < exclusions[particleIndex].size(); ii++) {
            outputExclusions[ii] = exclusions[particleIndex][ii];
        }
    }
}

void AmoebaCTForce::setCutoffDistance(double inputCutoff) {
    cutoff = inputCutoff;
}

double AmoebaCTForce::getCutoffDistance() const {
    return cutoff;
}

void AmoebaCTForce::setCutoff(double inputCutoff) {
    setCutoffDistance(inputCutoff);
}

double AmoebaCTForce::getCutoff() const {
    return getCutoffDistance();
}

AmoebaCTForce::NonbondedMethod AmoebaCTForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void AmoebaCTForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

ForceImpl* AmoebaCTForce::createImpl() const {
    return new AmoebaCTForceImpl(*this);
}

void AmoebaCTForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaCTForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
