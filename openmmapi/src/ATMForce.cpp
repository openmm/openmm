/* -------------------------------------------------------------------------- *
 *                    OpenMM's Alchemical Transfer Force                      *
 * -------------------------------------------------------------------------- *
 * This is a Force of the OpenMM molecular simulation toolkit                 *
 * that implements the Alchemical Transfer Potential                          *
 * for absolute and relative binding free energy estimation                   *
 * (https://doi.org/10.1021/acs.jcim.1c01129). The code is derived from the   *
 * ATMMetaForce plugin                                                        *
 * https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin               *
 * with support from the National Science Foundation CAREER 1750511           *
 *                                                                            *
 * Portions copyright (c) 2021-2023 by the Authors                            *
 * Authors: Emilio Gallicchio                                                 *
 * Contributors: Peter Eastman                                                *
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

#include "openmm/ATMForce.h"
#include "openmm/Force.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/internal/ATMForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>
#include <cmath>
#include <string>

using namespace OpenMM;
using namespace std;

ATMForce::ATMForce(const string& energy) : energyExpression(energy) {
}

ATMForce::ATMForce(double lambda1, double lambda2, double alpha, double uh, double w0, double umax, double ubcore, double acore, double direction) {
    if (alpha < 0)
        throw OpenMMException("ATMForce: alpha cannot be negative");
    if (lambda1 != lambda2 && alpha == 0)
        throw OpenMMException("ATMForce: alpha must be positive when lambda1 and lambda2 are different");
    if (umax < ubcore)
        throw OpenMMException("ATMForce: umax cannot be less than ubcore");
    if (acore < 0)
        throw OpenMMException("ATMForce: acore cannot be negative");
    if (direction != 1.0 && direction != -1.0)
        throw OpenMMException("ATMForce: direction must be either 1 or -1");
    string referencePotExpression = "select(step(Direction), u0, u1) + ";
    string alchemicalPotExpression = "select(Lambda2-Lambda1 , ((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0, Lambda2*usc + W0);";
    string softCoreExpression = "usc = select(Acore, select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);"
                                "fsc = (z^Acore-1)/(z^Acore+1);"
                                "z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;"
                                "y = (u-Ubcore)/(Umax-Ubcore);"
	                        "u = select(step(Direction), 1, -1)*(u1-u0)";
    setEnergyFunction(referencePotExpression + alchemicalPotExpression + softCoreExpression);
    addGlobalParameter(Lambda1(), lambda1);
    addGlobalParameter(Lambda2(), lambda2);
    addGlobalParameter(Alpha(), alpha);
    addGlobalParameter(Uh(), uh);
    addGlobalParameter(W0(), w0);
    addGlobalParameter(Umax(), umax);
    addGlobalParameter(Ubcore(), ubcore);
    addGlobalParameter(Acore(), acore);
    addGlobalParameter(Direction(), direction);
}

ATMForce::~ATMForce() {
    for (Force* force : forces)
        delete force;
}

const string& ATMForce::getEnergyFunction() const {
    return energyExpression;
}

void ATMForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

int ATMForce::addParticle(const Vec3& displacement1, const Vec3& displacement0) {
    particles.push_back(ParticleInfo(particles.size(), displacement1, displacement0));
    return particles.size()-1;
}

void ATMForce::getParticleParameters(int index, Vec3& displacement1, Vec3& displacement0) const {
    ASSERT_VALID_INDEX(index, particles);
    displacement1 = particles[index].displacement1;
    displacement0 = particles[index].displacement0;
}

void ATMForce::setParticleParameters(int index, const Vec3& displacement1, const Vec3& displacement0) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].displacement1 = displacement1;
    particles[index].displacement0 = displacement0;
}

int ATMForce::addForce(Force* force) {
    forces.push_back(force);
    return forces.size()-1;
}

Force& ATMForce::getForce(int index) const {
    ASSERT_VALID_INDEX(index, forces);
    return *forces[index];
}

int ATMForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& ATMForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void ATMForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double ATMForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void ATMForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void ATMForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& ATMForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

ForceImpl* ATMForce::createImpl() const {
    return new ATMForceImpl(*this);
}

void ATMForce::updateParametersInContext(OpenMM::Context& context) {
    dynamic_cast<ATMForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

bool ATMForce::usesPeriodicBoundaryConditions() const {
    for (auto& force : forces)
        if (force->usesPeriodicBoundaryConditions())
            return true;
    return false;
}

void ATMForce::getPerturbationEnergy(OpenMM::Context& context, double& u0, double& u1, double& energy) {
    dynamic_cast<ATMForceImpl&>(getImplInContext(context)).getPerturbationEnergy(getContextImpl(context), u0, u1, energy);
}

