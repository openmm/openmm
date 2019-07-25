/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

NonbondedForce::NonbondedForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), switchingDistance(-1.0), rfDielectric(78.3),
        ewaldErrorTol(5e-4), alpha(0.0), useSwitchingFunction(false), useDispersionCorrection(true), recipForceGroup(-1), nx(0), ny(0), nz(0) {
}

NonbondedForce::NonbondedMethod NonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void NonbondedForce::setNonbondedMethod(NonbondedMethod method) {
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

void NonbondedForce::setPMEParameters(double alpha, int nx, int ny, int nz) {
    this->alpha = alpha;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
}

void NonbondedForce::getPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const NonbondedForceImpl&>(getImplInContext(context)).getPMEParameters(alpha, nx, ny, nz);
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
}

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
}

ForceImpl* NonbondedForce::createImpl() const {
    return new NonbondedForceImpl(*this);
}

void NonbondedForce::createExceptionsFromBonds(const vector<pair<int, int> >& bonds, double coulomb14Scale, double lj14Scale) {
    for (int i = 0; i < (int) bonds.size(); ++i)
        if (bonds[i].first < 0 || bonds[i].second < 0 || bonds[i].first >= particles.size() || bonds[i].second >= particles.size())
            throw OpenMMException("createExceptionsFromBonds: Illegal particle index in list of bonds");

    // Find particles separated by 1, 2, or 3 bonds.

    vector<set<int> > exclusions(particles.size());
    vector<set<int> > bonded12(exclusions.size());
    for (int i = 0; i < (int) bonds.size(); ++i) {
        bonded12[bonds[i].first].insert(bonds[i].second);
        bonded12[bonds[i].second].insert(bonds[i].first);
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        addExclusionsToSet(bonded12, exclusions[i], i, i, 2);

    // Find particles separated by 1 or 2 bonds and create the exceptions.

    for (int i = 0; i < (int) exclusions.size(); ++i) {
        set<int> bonded13;
        addExclusionsToSet(bonded12, bonded13, i, i, 1);
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            if (*iter < i) {
                if (bonded13.find(*iter) == bonded13.end()) {
                    // This is a 1-4 interaction.

                    const ParticleInfo& particle1 = particles[*iter];
                    const ParticleInfo& particle2 = particles[i];
                    const double chargeProd = coulomb14Scale*particle1.charge*particle2.charge;
                    const double sigma = 0.5*(particle1.sigma+particle2.sigma);
                    const double epsilon = lj14Scale*std::sqrt(particle1.epsilon*particle2.epsilon);
                    addException(*iter, i, chargeProd, sigma, epsilon);
                }
                else {
                    // This interaction should be completely excluded.

                    addException(*iter, i, 0.0, 1.0, 0.0);
                }
            }
    }
}

void NonbondedForce::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const {
    for (set<int>::const_iterator iter = bonded12[fromParticle].begin(); iter != bonded12[fromParticle].end(); ++iter) {
        if (*iter != baseParticle)
            exclusions.insert(*iter);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, *iter, currentLevel-1);
    }
}

int NonbondedForce::getReciprocalSpaceForceGroup() const {
    return recipForceGroup;
}

void NonbondedForce::setReciprocalSpaceForceGroup(int group) {
    if (group < -1 || group > 31)
        throw OpenMMException("Force group must be between -1 and 31");
    recipForceGroup = group;
}

void NonbondedForce::updateParametersInContext(Context& context) {
    dynamic_cast<NonbondedForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
