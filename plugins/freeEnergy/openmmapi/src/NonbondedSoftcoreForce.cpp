/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "openmm/NonbondedSoftcoreForce.h"
#include "openmm/internal/NonbondedSoftcoreForceImpl.h"

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

NonbondedSoftcoreForce::NonbondedSoftcoreForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0),  rfDielectric(78.3), ewaldErrorTol(1e-4) {
}

NonbondedSoftcoreForce::NonbondedSoftcoreMethod NonbondedSoftcoreForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void NonbondedSoftcoreForce::setNonbondedMethod(NonbondedSoftcoreMethod method) {
    nonbondedMethod = method;
}

double NonbondedSoftcoreForce::getCutoffDistance() const {
    return cutoffDistance;
}

void NonbondedSoftcoreForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

double NonbondedSoftcoreForce::getReactionFieldDielectric( void ) const {
    return rfDielectric;
}

void NonbondedSoftcoreForce::setReactionFieldDielectric(double dielectric) {
    rfDielectric = dielectric;
}

int NonbondedSoftcoreForce::addParticle(double charge, double sigma, double epsilon, double softcoreLJLambda) {
    particles.push_back(ParticleInfo(charge, sigma, epsilon, softcoreLJLambda));
    return particles.size()-1;
}

void NonbondedSoftcoreForce::getParticleParameters(int index, double& charge, double& sigma, double& epsilon) const {
    charge   = particles[index].charge;
    sigma    = particles[index].sigma;
    epsilon  = particles[index].epsilon;
}

void NonbondedSoftcoreForce::getParticleParameters(int index, double& charge, double& sigma, double& epsilon, double& softcoreLJLambda) const {
    charge           = particles[index].charge;
    sigma            = particles[index].sigma;
    epsilon          = particles[index].epsilon;
    softcoreLJLambda = particles[index].softcoreLJLambda;
}

void NonbondedSoftcoreForce::setParticleParameters(int index, double charge, double sigma, double epsilon, double softcoreLJLambda) {
    particles[index].charge            = charge;
    particles[index].sigma             = sigma;
    particles[index].epsilon           = epsilon;
    particles[index].softcoreLJLambda  = softcoreLJLambda;
}

int NonbondedSoftcoreForce::addException(int particle1, int particle2, double chargeProd, double sigma, double epsilon, bool replace, double softcoreLJLambda ) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "NonbondedSoftcoreForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}
void NonbondedSoftcoreForce::getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd, double& sigma, double& epsilon ) const {
    particle1        = exceptions[index].particle1;
    particle2        = exceptions[index].particle2;
    chargeProd       = exceptions[index].chargeProd;
    sigma            = exceptions[index].sigma;
    epsilon          = exceptions[index].epsilon;
}

void NonbondedSoftcoreForce::getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd, double& sigma, double& epsilon, double& softcoreLJLambda ) const {
    particle1        = exceptions[index].particle1;
    particle2        = exceptions[index].particle2;
    chargeProd       = exceptions[index].chargeProd;
    sigma            = exceptions[index].sigma;
    epsilon          = exceptions[index].epsilon;
    softcoreLJLambda = exceptions[index].softcoreLJLambda;
}

void NonbondedSoftcoreForce::setExceptionParameters(int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon, double softcoreLJLambda) {
    exceptions[index].particle1           = particle1;
    exceptions[index].particle2           = particle2;
    exceptions[index].chargeProd          = chargeProd;
    exceptions[index].sigma               = sigma;
    exceptions[index].epsilon             = epsilon;
    exceptions[index].softcoreLJLambda    = softcoreLJLambda;
}

ForceImpl* NonbondedSoftcoreForce::createImpl() {
    return new NonbondedSoftcoreForceImpl(*this);
}

void NonbondedSoftcoreForce::createExceptionsFromBonds(const vector<pair<int, int> >& bonds, double coulomb14Scale, double lj14Scale) {

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

void NonbondedSoftcoreForce::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const {
    for (set<int>::const_iterator iter = bonded12[fromParticle].begin(); iter != bonded12[fromParticle].end(); ++iter) {
        if (*iter != baseParticle)
            exclusions.insert(*iter);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, *iter, currentLevel-1);
    }
}
