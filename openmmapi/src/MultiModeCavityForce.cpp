/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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

#include "openmm/MultiModeCavityForce.h"
#include "openmm/internal/MultiModeCavityForceImpl.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;

MultiModeCavityForce::MultiModeCavityForce(int numModes, double omega1, double lambda1,
                                           double cavityLength, double moleculeZ,
                                           double photonMass) :
        numModes(numModes), omega1(omega1), lambda1(lambda1),
        cavityLength(cavityLength), moleculeZ(moleculeZ),
        photonMass(photonMass), dsePrefactor(0.0) {
    if (numModes < 1)
        throw OpenMMException("MultiModeCavityForce: numModes must be >= 1");
    if (omega1 <= 0)
        throw OpenMMException("MultiModeCavityForce: omega1 must be positive");
    if (photonMass <= 0)
        throw OpenMMException("MultiModeCavityForce: photonMass must be positive");
    if (cavityLength <= 0)
        throw OpenMMException("MultiModeCavityForce: cavityLength must be positive");
    
    // Precompute spatial profiles and DSE prefactor
    recomputeDerivedQuantities();
}

void MultiModeCavityForce::setOmega1(double omega1) {
    if (omega1 <= 0)
        throw OpenMMException("MultiModeCavityForce: omega1 must be positive");
    this->omega1 = omega1;
    recomputeDerivedQuantities();
}

void MultiModeCavityForce::setLambda1(double lambda1) {
    this->lambda1 = lambda1;
    recomputeDerivedQuantities();
}

void MultiModeCavityForce::setPhotonMass(double mass) {
    if (mass <= 0)
        throw OpenMMException("MultiModeCavityForce: photonMass must be positive");
    this->photonMass = mass;
    recomputeDerivedQuantities();
}

int MultiModeCavityForce::addCavityParticle(int particleIndex) {
    if ((int)cavityParticleIndices.size() >= numModes) {
        std::stringstream msg;
        msg << "MultiModeCavityForce: already added " << numModes << " cavity particles (one per mode)";
        throw OpenMMException(msg.str());
    }
    cavityParticleIndices.push_back(particleIndex);
    return cavityParticleIndices.size() - 1;
}

int MultiModeCavityForce::getCavityParticleIndex(int modeIndex) const {
    if (modeIndex < 0 || modeIndex >= (int)cavityParticleIndices.size()) {
        std::stringstream msg;
        msg << "MultiModeCavityForce: mode index " << modeIndex << " out of range [0, " 
            << cavityParticleIndices.size() << ")";
        throw OpenMMException(msg.str());
    }
    return cavityParticleIndices[modeIndex];
}

void MultiModeCavityForce::recomputeDerivedQuantities() {
    // Precompute spatial profiles: f_n(z0) = sin(n * pi * z0 / L)
    // n is 1-indexed, so for mode index i (0-based), n = i+1
    spatialProfiles.resize(numModes);
    for (int i = 0; i < numModes; i++) {
        int n = i + 1;
        spatialProfiles[i] = std::sin(n * M_PI * moleculeZ / cavityLength);
    }
    
    // Precompute DSE prefactor: (1/2) * sum_n( eps_n^2 / K_n * f_n^2 )
    // eps_n = lambda_n * omega_n = sqrt(n) * lambda1 * n * omega1 = n^(3/2) * lambda1 * omega1
    // K_n = photonMass * omega_n^2 = photonMass * n^2 * omega1^2
    // eps_n^2 / K_n = (n^3 * lambda1^2 * omega1^2) / (photonMass * n^2 * omega1^2) = n * lambda1^2 / photonMass
    // So DSE prefactor = (1/2) * sum_n( n * lambda1^2 / photonMass * f_n^2 )
    // = (lambda1^2 / (2 * photonMass)) * sum_n( n * f_n^2 )
    dsePrefactor = 0.0;
    double prefactorBase = lambda1 * lambda1 / (2.0 * photonMass);
    for (int i = 0; i < numModes; i++) {
        int n = i + 1;
        double fn = spatialProfiles[i];
        dsePrefactor += n * fn * fn;
    }
    dsePrefactor *= prefactorBase;
}

double MultiModeCavityForce::getHarmonicEnergy(const Context& context) const {
    return dynamic_cast<const MultiModeCavityForceImpl&>(getImplInContext(context)).getHarmonicEnergy();
}

double MultiModeCavityForce::getCouplingEnergy(const Context& context) const {
    return dynamic_cast<const MultiModeCavityForceImpl&>(getImplInContext(context)).getCouplingEnergy();
}

double MultiModeCavityForce::getDipoleSelfEnergy(const Context& context) const {
    return dynamic_cast<const MultiModeCavityForceImpl&>(getImplInContext(context)).getDipoleSelfEnergy();
}

double MultiModeCavityForce::getTotalCavityEnergy(const Context& context) const {
    return getHarmonicEnergy(context) + getCouplingEnergy(context) + getDipoleSelfEnergy(context);
}

void MultiModeCavityForce::updateParametersInContext(Context& context) {
    dynamic_cast<MultiModeCavityForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* MultiModeCavityForce::createImpl() const {
    return new MultiModeCavityForceImpl(*this);
}
