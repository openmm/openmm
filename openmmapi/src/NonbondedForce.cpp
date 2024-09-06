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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

NonbondedForce::NonbondedForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), switchingDistance(-1.0), rfDielectric(78.3),
        ewaldErrorTol(5e-4), alpha(0.0), dalpha(0.0), useSwitchingFunction(false), useDispersionCorrection(true), exceptionsUsePeriodic(false), recipForceGroup(-1),
        includeDirectSpace(true), nx(0), ny(0), nz(0), dnx(0), dny(0), dnz(0), numContexts(0) {
}

NonbondedForce::NonbondedMethod NonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void NonbondedForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 5)
        throw OpenMMException("NonbondedForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double NonbondedForce::getCutoffDistance() const {
    return cutoffDistance;
}

void NonbondedForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

bool NonbondedForce::getUseSwitchingFunction() const {
    return useSwitchingFunction;
}

void NonbondedForce::setUseSwitchingFunction(bool use) {
    useSwitchingFunction = use;
}

double NonbondedForce::getSwitchingDistance() const {
    return switchingDistance;
}

void NonbondedForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}

double NonbondedForce::getReactionFieldDielectric() const {
    return rfDielectric;
}

void NonbondedForce::setReactionFieldDielectric(double dielectric) {
    rfDielectric = dielectric;
}

double NonbondedForce::getEwaldErrorTolerance() const {
    return ewaldErrorTol;
}

void NonbondedForce::setEwaldErrorTolerance(double tol) {
    ewaldErrorTol = tol;
}

void NonbondedForce::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->alpha;
    nx = this->nx;
    ny = this->ny;
    nz = this->nz;
}

void NonbondedForce::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->dalpha;
    nx = this->dnx;
    ny = this->dny;
    nz = this->dnz;
}

void NonbondedForce::setPMEParameters(double alpha, int nx, int ny, int nz) {
    this->alpha = alpha;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
}

void NonbondedForce::setLJPMEParameters(double alpha, int nx, int ny, int nz) {
    this->dalpha = alpha;
    this->dnx = nx;
    this->dny = ny;
    this->dnz = nz;
}

void NonbondedForce::getPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const NonbondedForceImpl&>(getImplInContext(context)).getPMEParameters(alpha, nx, ny, nz);
}

void NonbondedForce::getLJPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const NonbondedForceImpl&>(getImplInContext(context)).getLJPMEParameters(alpha, nx, ny, nz);
}

int NonbondedForce::addParticle(double charge, double sigma, double epsilon) {
    particles.push_back(ParticleInfo(charge, sigma, epsilon));
    return particles.size()-1;
}

void NonbondedForce::getParticleParameters(int index, double& charge, double& sigma, double& epsilon) const {
    ASSERT_VALID_INDEX(index, particles);
    charge = particles[index].charge;
    sigma = particles[index].sigma;
    epsilon = particles[index].epsilon;
}

void NonbondedForce::setParticleParameters(int index, double charge, double sigma, double epsilon) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].charge = charge;
    particles[index].sigma = sigma;
    particles[index].epsilon = epsilon;
    if (numContexts > 0) {
        firstChangedParticle = min(index, firstChangedParticle);
        lastChangedParticle = max(index, lastChangedParticle);
    }}

int NonbondedForce::addException(int particle1, int particle2, double chargeProd, double sigma, double epsilon, bool replace) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "NonbondedForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}
void NonbondedForce::getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd, double& sigma, double& epsilon) const {
    ASSERT_VALID_INDEX(index, exceptions);
    particle1 = exceptions[index].particle1;
    particle2 = exceptions[index].particle2;
    chargeProd = exceptions[index].chargeProd;
    sigma = exceptions[index].sigma;
    epsilon = exceptions[index].epsilon;
}

void NonbondedForce::setExceptionParameters(int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon) {
    ASSERT_VALID_INDEX(index, exceptions);
    exceptions[index].particle1 = particle1;
    exceptions[index].particle2 = particle2;
    exceptions[index].chargeProd = chargeProd;
    exceptions[index].sigma = sigma;
    exceptions[index].epsilon = epsilon;
    if (numContexts > 0) {
        firstChangedException = min(index, firstChangedException);
        lastChangedException = max(index, lastChangedException);
    }}

ForceImpl* NonbondedForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to particles and exceptions.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
        firstChangedException = exceptions.size();
        lastChangedException = -1;
    }
    numContexts++;
    return new NonbondedForceImpl(*this);
}

void NonbondedForce::createExceptionsFromBonds(const vector<pair<int, int> >& bonds, double coulomb14Scale, double lj14Scale) {
    for (auto& bond : bonds)
        if (bond.first < 0 || bond.second < 0 || bond.first >= particles.size() || bond.second >= particles.size())
            throw OpenMMException("createExceptionsFromBonds: Illegal particle index in list of bonds");

    // Find particles separated by 1, 2, or 3 bonds.

    vector<set<int> > exclusions(particles.size());
    vector<set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        bonded12[bond.first].insert(bond.second);
        bonded12[bond.second].insert(bond.first);
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        addExclusionsToSet(bonded12, exclusions[i], i, i, 2);

    // Find particles separated by 1 or 2 bonds and create the exceptions.

    for (int i = 0; i < (int) exclusions.size(); ++i) {
        set<int> bonded13;
        addExclusionsToSet(bonded12, bonded13, i, i, 1);
        for (int j : exclusions[i]) {
            if (j < i) {
                if (bonded13.find(j) == bonded13.end()) {
                    // This is a 1-4 interaction.

                    const ParticleInfo& particle1 = particles[j];
                    const ParticleInfo& particle2 = particles[i];
                    const double chargeProd = coulomb14Scale*particle1.charge*particle2.charge;
                    const double sigma = 0.5*(particle1.sigma+particle2.sigma);
                    const double epsilon = lj14Scale*std::sqrt(particle1.epsilon*particle2.epsilon);
                    addException(j, i, chargeProd, sigma, epsilon);
                }
                else {
                    // This interaction should be completely excluded.

                    addException(j, i, 0.0, 1.0, 0.0);
                }
            }
        }
    }
}

void NonbondedForce::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const {
    for (int i : bonded12[fromParticle]) {
        if (i != baseParticle)
            exclusions.insert(i);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, i, currentLevel-1);
    }
}

int NonbondedForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& NonbondedForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void NonbondedForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double NonbondedForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void NonbondedForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

int NonbondedForce::getGlobalParameterIndex(const std::string& parameter) const {
    for (int i = 0; i < globalParameters.size(); i++)
        if (globalParameters[i].name == parameter)
            return i;
    throw OpenMMException("There is no global parameter called '"+parameter+"'");
}

int NonbondedForce::addParticleParameterOffset(const std::string& parameter, int particleIndex, double chargeScale, double sigmaScale, double epsilonScale) {
    particleOffsets.push_back(ParticleOffsetInfo(getGlobalParameterIndex(parameter), particleIndex, chargeScale, sigmaScale, epsilonScale));
    return particleOffsets.size()-1;
}

void NonbondedForce::getParticleParameterOffset(int index, std::string& parameter, int& particleIndex, double& chargeScale, double& sigmaScale, double& epsilonScale) const {
    ASSERT_VALID_INDEX(index, particleOffsets);
    parameter = globalParameters[particleOffsets[index].parameter].name;
    particleIndex = particleOffsets[index].particle;
    chargeScale = particleOffsets[index].chargeScale;
    sigmaScale = particleOffsets[index].sigmaScale;
    epsilonScale = particleOffsets[index].epsilonScale;
}

void NonbondedForce::setParticleParameterOffset(int index, const std::string& parameter, int particleIndex, double chargeScale, double sigmaScale, double epsilonScale) {
    ASSERT_VALID_INDEX(index, particleOffsets);
    particleOffsets[index].parameter = getGlobalParameterIndex(parameter);
    particleOffsets[index].particle = particleIndex;
    particleOffsets[index].chargeScale = chargeScale;
    particleOffsets[index].sigmaScale = sigmaScale;
    particleOffsets[index].epsilonScale = epsilonScale;
}

int NonbondedForce::addExceptionParameterOffset(const std::string& parameter, int exceptionIndex, double chargeProdScale, double sigmaScale, double epsilonScale) {
    exceptionOffsets.push_back(ExceptionOffsetInfo(getGlobalParameterIndex(parameter), exceptionIndex, chargeProdScale, sigmaScale, epsilonScale));
    return exceptionOffsets.size()-1;
}

void NonbondedForce::getExceptionParameterOffset(int index, std::string& parameter, int& exceptionIndex, double& chargeProdScale, double& sigmaScale, double& epsilonScale) const {
    ASSERT_VALID_INDEX(index, exceptionOffsets);
    parameter = globalParameters[exceptionOffsets[index].parameter].name;
    exceptionIndex = exceptionOffsets[index].exception;
    chargeProdScale = exceptionOffsets[index].chargeProdScale;
    sigmaScale = exceptionOffsets[index].sigmaScale;
    epsilonScale = exceptionOffsets[index].epsilonScale;
}

void NonbondedForce::setExceptionParameterOffset(int index, const std::string& parameter, int exceptionIndex, double chargeProdScale, double sigmaScale, double epsilonScale) {
    ASSERT_VALID_INDEX(index, exceptionOffsets);
    exceptionOffsets[index].parameter = getGlobalParameterIndex(parameter);
    exceptionOffsets[index].exception = exceptionIndex;
    exceptionOffsets[index].chargeProdScale = chargeProdScale;
    exceptionOffsets[index].sigmaScale = sigmaScale;
    exceptionOffsets[index].epsilonScale = epsilonScale;
}

int NonbondedForce::getReciprocalSpaceForceGroup() const {
    return recipForceGroup;
}

void NonbondedForce::setReciprocalSpaceForceGroup(int group) {
    if (group < -1 || group > 31)
        throw OpenMMException("Force group must be between -1 and 31");
    recipForceGroup = group;
}

bool NonbondedForce::getIncludeDirectSpace() const {
    return includeDirectSpace;
}

void NonbondedForce::setIncludeDirectSpace(bool include) {
    includeDirectSpace = include;
}

void NonbondedForce::updateParametersInContext(Context& context) {
    dynamic_cast<NonbondedForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context),
            firstChangedParticle, lastChangedParticle, firstChangedException, lastChangedException);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed particles and exceptions.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
        firstChangedException = exceptions.size();
        lastChangedException = -1;
    }
}

bool NonbondedForce::getExceptionsUsePeriodicBoundaryConditions() const {
    return exceptionsUsePeriodic;
}

void NonbondedForce::setExceptionsUsePeriodicBoundaryConditions(bool periodic) {
    exceptionsUsePeriodic = periodic;
}
