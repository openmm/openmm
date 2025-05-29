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
 * Portions copyright (c) 2021-2024 by the Authors                            *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/ATMForceImpl.h"

#include "openmm/NonbondedForce.h"
#include "openmm/kernels.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/Vec3.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>

#include <iostream>

using namespace OpenMM;
using namespace std;

ATMForceImpl::ATMForceImpl(const ATMForce& owner) : owner(owner), innerIntegrator0(1.0), innerIntegrator1(1.0),
        innerContext0(NULL), innerContext1(NULL) {
    Lepton::ParsedExpression expr = Lepton::Parser::parse(owner.getEnergyFunction()).optimize();
    energyExpression = expr.createCompiledExpression();
    u0DerivExpression = expr.differentiate("u0").createCompiledExpression();
    u1DerivExpression = expr.differentiate("u1").createCompiledExpression();
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(owner.getGlobalParameterName(i));
    globalValues.resize(globalParameterNames.size());
    map<string, double*> variableLocations;
    variableLocations["u0"] = &state0Energy;
    variableLocations["u1"] = &state1Energy;
    for (int i = 0; i < globalParameterNames.size(); i++)
        variableLocations[globalParameterNames[i]] = &globalValues[i];
    energyExpression.setVariableLocations(variableLocations);
    u0DerivExpression.setVariableLocations(variableLocations);
    u1DerivExpression.setVariableLocations(variableLocations);
    for (int i = 0; i < owner.getNumEnergyParameterDerivatives(); i++) {
        string name = owner.getEnergyParameterDerivativeName(i);
        paramDerivNames.push_back(name);
        paramDerivExpressions.push_back(expr.differentiate(name).createCompiledExpression());
        paramDerivExpressions[i].setVariableLocations(variableLocations);
    }
}

ATMForceImpl::~ATMForceImpl() {
    if (innerContext0 != NULL)
        delete innerContext0;
    if (innerContext1 != NULL)
        delete innerContext1;
}

void ATMForceImpl::copySystem(ContextImpl& context, const OpenMM::System& system, OpenMM::System& innerSystem) {
    //copy particles
    for (int i = 0; i < system.getNumParticles(); i++)
        innerSystem.addParticle(system.getParticleMass(i));

    //copy periodic box dimensions
    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    innerSystem.setDefaultPeriodicBoxVectors(a, b, c);

    // Add forces to the inner contexts
    for (int i = 0; i < owner.getNumForces(); i++) {
        const Force &force = owner.getForce(i);
        innerSystem.addForce(XmlSerializer::clone<Force>(force));
    }
}

void ATMForceImpl::initialize(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();

    copySystem(context, system, innerSystem0);
    copySystem(context, system, innerSystem1);

    // Create the inner context.

    innerContext0 = context.createLinkedContext(innerSystem0, innerIntegrator0);
    innerContext1 = context.createLinkedContext(innerSystem1, innerIntegrator1);

    // Create the kernel.

    kernel = context.getPlatform().createKernel(CalcATMForceKernel::Name(), context);
    kernel.getAs<CalcATMForceKernel>().initialize(context.getSystem(), owner);
}

double ATMForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) == 0)
        return 0.0;

    ContextImpl& innerContextImpl0 = getContextImpl(*innerContext0);
    ContextImpl& innerContextImpl1 = getContextImpl(*innerContext1);

    // Copy the coordinates etc. from the context to the inner contexts

    kernel.getAs<CalcATMForceKernel>().copyState(context, innerContextImpl0, innerContextImpl1);

    // Evaluate energy and forces for the two systems

    state0Energy = innerContextImpl0.calcForcesAndEnergy(includeForces, true);
    state1Energy = innerContextImpl1.calcForcesAndEnergy(includeForces, true);

    // set global parameters for energy expression

    for (int i = 0; i < globalParameterNames.size(); i++)
        globalValues[i] = context.getParameter(globalParameterNames[i]);

    // Protect against overflow when the hybrid potential function does
    // not depend on u0 or u1 and their values are unbounded; typically at the endstates

    double dEdu0 = u0DerivExpression.evaluate();
    double dEdu1 = u1DerivExpression.evaluate();
    double epsi = std::numeric_limits<float>::min();
    double maxEnergy = std::numeric_limits<float>::max();
    if(fabs(dEdu0) < epsi && (isnan(state0Energy) || isinf(state0Energy)))
	state0Energy = maxEnergy;
    if(fabs(dEdu1) < epsi && (isnan(state1Energy) || isinf(state1Energy)))
	state1Energy = maxEnergy;

    // Compute the alchemical energy and forces.

    combinedEnergy = energyExpression.evaluate();
    if (includeForces) {
        map<string, double> energyParamDerivs;
        for (int i = 0; i < paramDerivExpressions.size(); i++)
            energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate();
        kernel.getAs<CalcATMForceKernel>().applyForces(context, innerContextImpl0, innerContextImpl1, dEdu0, dEdu1, energyParamDerivs);
    }

    return (includeEnergy ? combinedEnergy : 0.0);
}

std::map<std::string, double> ATMForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters.insert(innerContext0->getParameters().begin(), innerContext0->getParameters().end());
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

std::vector<std::string> ATMForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcATMForceKernel::Name());
    return names;
}

vector<pair<int, int> > ATMForceImpl::getBondedParticles() const {
    vector<pair<int, int> > bonds;
    const ContextImpl& innerContextImpl = getContextImpl(*innerContext0);
    for (auto& impl : innerContextImpl.getForceImpls()) {
        for (auto& bond : impl->getBondedParticles())
            bonds.push_back(bond);
    }
    return bonds;
}

void ATMForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcATMForceKernel>().copyParametersToContext(context, owner);
}

void ATMForceImpl::getPerturbationEnergy(ContextImpl& context, double& u1, double& u0, double& energy) {
    calcForcesAndEnergy(context, false, true, -1);
    u0 = state0Energy;
    u1 = state1Energy;
    energy = combinedEnergy;
}
