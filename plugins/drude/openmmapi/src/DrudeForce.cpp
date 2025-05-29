/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2024 Stanford University and the Authors.      *
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
#include "openmm/DrudeForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/DrudeForceImpl.h"

using namespace OpenMM;
using namespace std;

DrudeForce::DrudeForce() : usePeriodic(false) {
}

int DrudeForce::addParticle(int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34) {
    if (polarizability <= 0)
        throw OpenMMException("Polarizability must be positive");
    if ((aniso12 <= 0 && particle2 != -1) || (aniso34 <= 0 && particle3 != -1 && particle4 != -1))
        throw OpenMMException("Anisotropy factors must be positive");
    particles.push_back(ParticleInfo(particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34));
    return particles.size()-1;
}

void DrudeForce::getParticleParameters(int index, int& particle, int& particle1, int& particle2, int& particle3, int& particle4, double& charge, double& polarizability, double& aniso12, double& aniso34) const {
    ASSERT_VALID_INDEX(index, particles);
    particle = particles[index].particle;
    particle1 = particles[index].particle1;
    particle2 = particles[index].particle2;
    particle3 = particles[index].particle3;
    particle4 = particles[index].particle4;
    charge = particles[index].charge;
    polarizability = particles[index].polarizability;
    aniso12 = particles[index].aniso12;
    aniso34 = particles[index].aniso34;
}

void DrudeForce::setParticleParameters(int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34) {
    ASSERT_VALID_INDEX(index, particles);
    if (polarizability <= 0)
        throw OpenMMException("Polarizability must be positive");
    if ((aniso12 <= 0 && particle2 != -1) || (aniso34 <= 0 && particle3 != -1 && particle4 != -1))
        throw OpenMMException("Anisotropy factors must be positive");
    particles[index].particle = particle;
    particles[index].particle1 = particle1;
    particles[index].particle2 = particle2;
    particles[index].particle3 = particle3;
    particles[index].particle4 = particle4;
    particles[index].charge = charge;
    particles[index].polarizability = polarizability;
    particles[index].aniso12 = aniso12;
    particles[index].aniso34 = aniso34;
}

int DrudeForce::addScreenedPair(int particle1, int particle2, double thole) {
    screenedPairs.push_back(ScreenedPairInfo(particle1, particle2, thole));
    return screenedPairs.size()-1;
}
void DrudeForce::getScreenedPairParameters(int index, int& particle1, int& particle2, double& thole) const {
    ASSERT_VALID_INDEX(index, screenedPairs);
    particle1 = screenedPairs[index].particle1;
    particle2 = screenedPairs[index].particle2;
    thole = screenedPairs[index].thole;
}

void DrudeForce::setScreenedPairParameters(int index, int particle1, int particle2, double thole) {
    ASSERT_VALID_INDEX(index, screenedPairs);
    screenedPairs[index].particle1 = particle1;
    screenedPairs[index].particle2 = particle2;
    screenedPairs[index].thole = thole;
}

ForceImpl* DrudeForce::createImpl() const {
    return new DrudeForceImpl(*this);
}

void DrudeForce::updateParametersInContext(Context& context) {
    dynamic_cast<DrudeForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

void DrudeForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usePeriodic = periodic;
}

bool DrudeForce::usesPeriodicBoundaryConditions() const {
    return usePeriodic;
}
