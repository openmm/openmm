/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2017 Stanford University and the Authors.      *
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

#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCVForceImpl.h"
#include "openmm/kernels.h"
#include "openmm/serialization/XmlSerializer.h"
#include <map>

using namespace OpenMM;
using namespace std;

CustomCVForceImpl::CustomCVForceImpl(const CustomCVForce& owner) : owner(owner), innerIntegrator(1.0),
        innerContext(NULL) {
}

CustomCVForceImpl::~CustomCVForceImpl() {
    if (innerContext != NULL)
        delete innerContext;
}

void CustomCVForceImpl::initialize(ContextImpl& context) {
    // Construct the inner system used to evaluate collective variables.
    
    const System& system = context.getSystem();
    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    innerSystem.setDefaultPeriodicBoxVectors(a, b, c);
    for (int i = 0; i < system.getNumParticles(); i++)
        innerSystem.addParticle(system.getParticleMass(i));
    for (int i = 0; i < owner.getNumCollectiveVariables(); i++) {
        Force* variable = XmlSerializer::clone<Force>(owner.getCollectiveVariable(i));
        variable->setForceGroup(i);
        NonbondedForce* nonbonded = dynamic_cast<NonbondedForce*>(variable);
        if (nonbonded != NULL)
            nonbonded->setReciprocalSpaceForceGroup(-1);
        innerSystem.addForce(variable);
    }
    
    // Create the inner context.
    
    innerContext = context.createLinkedContext(innerSystem, innerIntegrator);
    vector<Vec3> positions(system.getNumParticles(), Vec3());
    innerContext->setPositions(positions);
    
    // Create the kernel.
    
    kernel = context.getPlatform().createKernel(CalcCustomCVForceKernel::Name(), context);
    kernel.getAs<CalcCustomCVForceKernel>().initialize(context.getSystem(), owner, getContextImpl(*innerContext));
}

double CustomCVForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomCVForceKernel>().execute(context, getContextImpl(*innerContext), includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomCVForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomCVForceKernel::Name());
    return names;
}

map<string, double> CustomCVForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters.insert(innerContext->getParameters().begin(), innerContext->getParameters().end());
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

void CustomCVForceImpl::getCollectiveVariableValues(ContextImpl& context, vector<double>& values) {
    kernel.getAs<CalcCustomCVForceKernel>().copyState(context, getContextImpl(*innerContext));
    values.clear();
    for (int i = 0; i < innerSystem.getNumForces(); i++) {
        double value = innerContext->getState(State::Energy, false, 1<<i).getPotentialEnergy();
        values.push_back(value);
    }
}

Context& CustomCVForceImpl::getInnerContext() {
    return *innerContext;
}

void CustomCVForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomCVForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
