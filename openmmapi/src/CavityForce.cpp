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

#include "openmm/CavityForce.h"
#include "openmm/internal/CavityForceImpl.h"
#include "openmm/OpenMMException.h"
#include <algorithm>

using namespace OpenMM;

CavityForce::CavityForce(int cavityParticleIndex, double omegac, double lambdaCoupling, double photonMass) :
        cavityParticleIndex(cavityParticleIndex), omegac(omegac), 
        lambdaCoupling(lambdaCoupling), photonMass(photonMass),
        couplingOnStep(-1), couplingOnValue(0.0),
        cavityDriveEnabled(false), cavityDriveAmplitude(0.0), cavityDriveFrequency(0.0), 
        cavityDrivePhase(0.0), cavityDriveEnvelopeType(0), cavityDriveEnvParam1(0.0), cavityDriveEnvParam2(0.0),
        directLaserEnabled(false), directLaserAmplitude(0.0), directLaserFrequency(0.0),
        directLaserPhase(0.0), directLaserEnvelopeType(0), directLaserEnvParam1(0.0), directLaserEnvParam2(0.0) {
    if (omegac <= 0)
        throw OpenMMException("CavityForce: omegac must be positive");
    if (photonMass <= 0)
        throw OpenMMException("CavityForce: photonMass must be positive");
}

void CavityForce::setCouplingOnStep(int step, double value) {
    couplingOnStep = step;
    couplingOnValue = value;
}

void CavityForce::setLambdaCouplingSchedule(const std::vector<std::pair<int, double>>& schedule) {
    couplingSchedule = schedule;
    // Sort by timestep
    std::sort(couplingSchedule.begin(), couplingSchedule.end(),
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.first < b.first;
              });
}

double CavityForce::getLambdaCouplingAtStep(int step) const {
    // First check the simple on-step mechanism
    if (couplingOnStep >= 0) {
        if (step >= couplingOnStep)
            return couplingOnValue;
        else
            return lambdaCoupling;
    }
    
    // Fall back to the full schedule if set
    if (couplingSchedule.empty())
        return lambdaCoupling;
    
    // Find the schedule entry that applies at this step
    double value = lambdaCoupling;
    for (const auto& entry : couplingSchedule) {
        if (entry.first <= step)
            value = entry.second;
        else
            break;
    }
    return value;
}

double CavityForce::getHarmonicEnergy(const Context& context) const {
    return dynamic_cast<const CavityForceImpl&>(getImplInContext(context)).getHarmonicEnergy();
}

double CavityForce::getCouplingEnergy(const Context& context) const {
    return dynamic_cast<const CavityForceImpl&>(getImplInContext(context)).getCouplingEnergy();
}

double CavityForce::getDipoleSelfEnergy(const Context& context) const {
    return dynamic_cast<const CavityForceImpl&>(getImplInContext(context)).getDipoleSelfEnergy();
}

double CavityForce::getTotalCavityEnergy(const Context& context) const {
    return getHarmonicEnergy(context) + getCouplingEnergy(context) + getDipoleSelfEnergy(context);
}

// Case (2): Cavity-mode driving methods
void CavityForce::setCavityDriveAmplitude(double f0) {
    cavityDriveAmplitude = f0;
}

double CavityForce::getCavityDriveAmplitude() const {
    return cavityDriveAmplitude;
}

void CavityForce::setCavityDriveFrequency(double omega_d) {
    cavityDriveFrequency = omega_d;
}

double CavityForce::getCavityDriveFrequency() const {
    return cavityDriveFrequency;
}

void CavityForce::setCavityDrivePhase(double phi) {
    cavityDrivePhase = phi;
}

double CavityForce::getCavityDrivePhase() const {
    return cavityDrivePhase;
}

void CavityForce::setCavityDriveEnvelope(const std::string& type, double param1, double param2) {
    if (type == "constant") {
        cavityDriveEnvelopeType = 0;
    } else if (type == "gaussian") {
        cavityDriveEnvelopeType = 1;
    } else if (type == "square") {
        cavityDriveEnvelopeType = 2;
    } else if (type == "exponential") {
        cavityDriveEnvelopeType = 3;
    } else {
        throw OpenMMException("CavityForce: Unknown envelope type: " + type);
    }
    cavityDriveEnvParam1 = param1;
    cavityDriveEnvParam2 = param2;
}

int CavityForce::getCavityDriveEnvelopeType() const {
    return cavityDriveEnvelopeType;
}

double CavityForce::getCavityDriveEnvelopeParam1() const {
    return cavityDriveEnvParam1;
}

double CavityForce::getCavityDriveEnvelopeParam2() const {
    return cavityDriveEnvParam2;
}

void CavityForce::setCavityDriveEnabled(bool enable) {
    cavityDriveEnabled = enable;
}

bool CavityForce::getCavityDriveEnabled() const {
    return cavityDriveEnabled;
}

// Case (1): Direct molecule-laser coupling methods
void CavityForce::setDirectLaserAmplitude(double E0) {
    directLaserAmplitude = E0;
}

double CavityForce::getDirectLaserAmplitude() const {
    return directLaserAmplitude;
}

void CavityForce::setDirectLaserFrequency(double omega_L) {
    directLaserFrequency = omega_L;
}

double CavityForce::getDirectLaserFrequency() const {
    return directLaserFrequency;
}

void CavityForce::setDirectLaserPhase(double phi) {
    directLaserPhase = phi;
}

double CavityForce::getDirectLaserPhase() const {
    return directLaserPhase;
}

void CavityForce::setDirectLaserEnvelope(const std::string& type, double param1, double param2) {
    if (type == "constant") {
        directLaserEnvelopeType = 0;
    } else if (type == "gaussian") {
        directLaserEnvelopeType = 1;
    } else if (type == "square") {
        directLaserEnvelopeType = 2;
    } else if (type == "exponential") {
        directLaserEnvelopeType = 3;
    } else {
        throw OpenMMException("CavityForce: Unknown envelope type: " + type);
    }
    directLaserEnvParam1 = param1;
    directLaserEnvParam2 = param2;
}

int CavityForce::getDirectLaserEnvelopeType() const {
    return directLaserEnvelopeType;
}

double CavityForce::getDirectLaserEnvelopeParam1() const {
    return directLaserEnvParam1;
}

double CavityForce::getDirectLaserEnvelopeParam2() const {
    return directLaserEnvParam2;
}

void CavityForce::setDirectLaserCouplingEnabled(bool enable) {
    directLaserEnabled = enable;
}

bool CavityForce::getDirectLaserCouplingEnabled() const {
    return directLaserEnabled;
}

double CavityForce::getCavityDriveEnergy(const Context& context) const {
    return dynamic_cast<const CavityForceImpl&>(getImplInContext(context)).getCavityDriveEnergy();
}

double CavityForce::getDirectLaserEnergy(const Context& context) const {
    return dynamic_cast<const CavityForceImpl&>(getImplInContext(context)).getDirectLaserEnergy();
}

void CavityForce::updateParametersInContext(Context& context) {
    dynamic_cast<CavityForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* CavityForce::createImpl() const {
    return new CavityForceImpl(*this);
}
