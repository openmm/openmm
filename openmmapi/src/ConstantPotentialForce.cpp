/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "openmm/ConstantPotentialForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ConstantPotentialForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

ConstantPotentialForce::ConstantPotentialForce() : constantPotentialMethod(CG),
        cutoffDistance(1.0), ewaldErrorTol(5e-4), alpha(0.0), cgErrorTol(1e-4),
        chargeTarget(0.0), exceptionsUsePeriodic(false),
        useChargeConstraint(false), usePreconditioner(true),
        nx(0), ny(0), nz(0), numContexts(0) {
}

double ConstantPotentialForce::getCutoffDistance() const {
    return cutoffDistance;
}

void ConstantPotentialForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

double ConstantPotentialForce::getEwaldErrorTolerance() const {
    return ewaldErrorTol;
}

void ConstantPotentialForce::setEwaldErrorTolerance(double tol) {
    ewaldErrorTol = tol;
}

void ConstantPotentialForce::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = this->alpha;
    nx = this->nx;
    ny = this->ny;
    nz = this->nz;
}

void ConstantPotentialForce::setPMEParameters(double alpha, int nx, int ny, int nz) {
    this->alpha = alpha;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
}

void ConstantPotentialForce::getPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const {
    dynamic_cast<const ConstantPotentialForceImpl&>(getImplInContext(context)).getPMEParameters(alpha, nx, ny, nz);
}

int ConstantPotentialForce::addParticle(double charge) {
    particles.push_back(ParticleInfo(charge));
    return particles.size()-1;
}

void ConstantPotentialForce::getParticleParameters(int index, double& charge) const {
    ASSERT_VALID_INDEX(index, particles);
    charge = particles[index].charge;
}

void ConstantPotentialForce::setParticleParameters(int index, double charge) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].charge = charge;
    if (numContexts > 0) {
        firstChangedParticle = min(index, firstChangedParticle);
        lastChangedParticle = max(index, lastChangedParticle);
    }
}

int ConstantPotentialForce::addException(int particle1, int particle2, double chargeProd, bool replace) {
    map<pair<int, int>, int>::iterator iter = exceptionMap.find(pair<int, int>(particle1, particle2));
    int newIndex;
    if (iter == exceptionMap.end())
        iter = exceptionMap.find(pair<int, int>(particle2, particle1));
    if (iter != exceptionMap.end()) {
        if (!replace) {
            stringstream msg;
            msg << "ConstantPotentialForce: There is already an exception for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[iter->second] = ExceptionInfo(particle1, particle2, chargeProd);
        newIndex = iter->second;
        exceptionMap.erase(iter->first);
    }
    else {
        exceptions.push_back(ExceptionInfo(particle1, particle2, chargeProd));
        newIndex = exceptions.size()-1;
    }
    exceptionMap[pair<int, int>(particle1, particle2)] = newIndex;
    return newIndex;
}
void ConstantPotentialForce::getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd) const {
    ASSERT_VALID_INDEX(index, exceptions);
    particle1 = exceptions[index].particle1;
    particle2 = exceptions[index].particle2;
    chargeProd = exceptions[index].chargeProd;
}

void ConstantPotentialForce::setExceptionParameters(int index, int particle1, int particle2, double chargeProd) {
    ASSERT_VALID_INDEX(index, exceptions);
    exceptions[index].particle1 = particle1;
    exceptions[index].particle2 = particle2;
    exceptions[index].chargeProd = chargeProd;
    if (numContexts > 0) {
        firstChangedException = min(index, firstChangedException);
        lastChangedException = max(index, lastChangedException);
    }}

ForceImpl* ConstantPotentialForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to particles and exceptions.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
        firstChangedException = exceptions.size();
        lastChangedException = -1;
        firstChangedElectrode = electrodes.size();
        lastChangedElectrode = -1;
    }
    numContexts++;
    return new ConstantPotentialForceImpl(*this);
}

void ConstantPotentialForce::createExceptionsFromBonds(const vector<pair<int, int> >& bonds, double coulomb14Scale) {
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
                    addException(j, i, chargeProd);
                }
                else {
                    // This interaction should be completely excluded.

                    addException(j, i, 0.0);
                }
            }
        }
    }
}

void ConstantPotentialForce::addExclusionsToSet(const vector<set<int> >& bonded12, set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const {
    for (int i : bonded12[fromParticle]) {
        if (i != baseParticle)
            exclusions.insert(i);
        if (currentLevel > 0)
            addExclusionsToSet(bonded12, exclusions, baseParticle, i, currentLevel-1);
    }
}

void ConstantPotentialForce::updateParametersInContext(Context& context) {
    dynamic_cast<ConstantPotentialForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context),
            firstChangedParticle, lastChangedParticle, firstChangedException, lastChangedException, firstChangedElectrode, lastChangedElectrode);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can
        // reset the tracking of changed particles, exceptions, and electrodes.
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
        firstChangedException = exceptions.size();
        lastChangedException = -1;
        firstChangedElectrode = electrodes.size();
        lastChangedElectrode = -1;
    }
}

bool ConstantPotentialForce::getExceptionsUsePeriodicBoundaryConditions() const {
    return exceptionsUsePeriodic;
}

void ConstantPotentialForce::setExceptionsUsePeriodicBoundaryConditions(bool periodic) {
    exceptionsUsePeriodic = periodic;
}

ConstantPotentialForce::ConstantPotentialMethod ConstantPotentialForce::getConstantPotentialMethod() const {
    return constantPotentialMethod;
}

void ConstantPotentialForce::setConstantPotentialMethod(ConstantPotentialMethod method) {
    if (method < 0 || method > 1) {
        throw OpenMMException("ConstantPotentialForce: Illegal value for constant potential method");
    }
    constantPotentialMethod = method;
}

bool ConstantPotentialForce::getUsePreconditioner() const {
    return usePreconditioner;
}

void ConstantPotentialForce::setUsePreconditioner(bool use) {
    usePreconditioner = use;
}

double ConstantPotentialForce::getCGErrorTolerance() const {
    return cgErrorTol;
}

void ConstantPotentialForce::setCGErrorTolerance(double tol) {
    cgErrorTol = tol;
}

int ConstantPotentialForce::addElectrode(const std::set<int>& electrodeParticles, double potential, double gaussianWidth, double thomasFermiScale) {
    electrodes.push_back(ElectrodeInfo(electrodeParticles, potential, gaussianWidth, thomasFermiScale));
    return electrodes.size() - 1;
}

void ConstantPotentialForce::getElectrodeParameters(int index, std::set<int>& electrodeParticles, double& potential, double& gaussianWidth, double& thomasFermiScale) const {
    ASSERT_VALID_INDEX(index, electrodes);
    electrodeParticles = electrodes[index].particles;
    potential = electrodes[index].potential;
    gaussianWidth = electrodes[index].gaussianWidth;
    thomasFermiScale = electrodes[index].thomasFermiScale;
}

void ConstantPotentialForce::setElectrodeParameters(int index, const std::set<int>& electrodeParticles, double potential, double gaussianWidth, double thomasFermiScale) {
    ASSERT_VALID_INDEX(index, electrodes);
    electrodes[index].particles = electrodeParticles;
    electrodes[index].potential = potential;
    electrodes[index].gaussianWidth = gaussianWidth;
    electrodes[index].thomasFermiScale = thomasFermiScale;
    if (numContexts > 0) {
        firstChangedElectrode = min(index, firstChangedElectrode);
        lastChangedElectrode = max(index, lastChangedElectrode);
    }
}

bool ConstantPotentialForce::getUseChargeConstraint() const {
    return useChargeConstraint;
}

void ConstantPotentialForce::setUseChargeConstraint(bool use) {
    useChargeConstraint = use;
}

double ConstantPotentialForce::getChargeConstraintTarget() const {
    return chargeTarget;
}

void ConstantPotentialForce::setChargeConstraintTarget(double charge) {
    chargeTarget = charge;
}

void ConstantPotentialForce::getExternalField(Vec3& field) const {
    field = externalField;
}

void ConstantPotentialForce::setExternalField(const Vec3& field) {
    externalField = field;
}

void ConstantPotentialForce::getCharges(Context& context, std::vector<double>& charges) const {
    dynamic_cast<ConstantPotentialForceImpl&>(getImplInContext(context)).getCharges(getContextImpl(context), charges);
}
