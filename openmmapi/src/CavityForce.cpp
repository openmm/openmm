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
        couplingOnStep(-1), couplingOnValue(0.0) {
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

void CavityForce::updateParametersInContext(Context& context) {
    dynamic_cast<CavityForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* CavityForce::createImpl() const {
    return new CavityForceImpl(*this);
}
