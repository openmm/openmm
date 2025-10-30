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
#include "openmm/LCPOForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/LCPOForceImpl.h"

using namespace OpenMM;
using namespace std;

LCPOForce::LCPOForce() : usesPeriodic(false), numContexts(0) {
}

void LCPOForce::setUsesPeriodicBoundaryConditions(bool periodic) {
    usesPeriodic = periodic;
}

int LCPOForce::getNumParticles() const {
    return particles.size();
}

int LCPOForce::addParticle(double radius, double p1, double p2, double p3, double p4) {
    particles.push_back(ParticleInfo(radius, p1, p2, p3, p4));
    return particles.size() - 1;
}

void LCPOForce::getParticleParameters(int index, double& radius, double& p1, double& p2, double& p3, double& p4) const {
    ASSERT_VALID_INDEX(index, particles);
    const ParticleInfo& info = particles[index];
    radius = info.radius;
    p1 = info.p1;
    p2 = info.p2;
    p3 = info.p3;
    p4 = info.p4;
}

void LCPOForce::setParticleParameters(int index, double radius, double p1, double p2, double p3, double p4) {
    ASSERT_VALID_INDEX(index, particles);
    ParticleInfo& info = particles[index];
    info.radius = radius;
    info.p1 = p1;
    info.p2 = p2;
    info.p3 = p3;
    info.p4 = p4;
    if (numContexts > 0) {
        firstChangedParticle = min(index, firstChangedParticle);
        lastChangedParticle = max(index, lastChangedParticle);
    }
}

void LCPOForce::updateParametersInContext(Context& context) {
    dynamic_cast<LCPOForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedParticle, lastChangedParticle);
    if (numContexts == 1) {
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
    }
}

bool LCPOForce::usesPeriodicBoundaryConditions() const {
    return usesPeriodic;
}

ForceImpl* LCPOForce::createImpl() const {
    if (numContexts == 0) {
        firstChangedParticle = particles.size();
        lastChangedParticle = -1;
    }
    numContexts++;
    return new LCPOForceImpl(*this);
}
