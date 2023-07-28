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
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include <iostream>

using namespace OpenMM;
using namespace std;

ATMForceImpl::ATMForceImpl(const ATMForce& owner) : owner(owner), innerIntegrator1(1.0), innerIntegrator2(1.0) {
}

ATMForceImpl::~ATMForceImpl() {
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

    copySystem(context, system, innerSystem1);
    copySystem(context, system, innerSystem2);

    // Create the inner context.

    innerContext1 = context.createLinkedContext(innerSystem1, innerIntegrator1);
    innerContext2 = context.createLinkedContext(innerSystem2, innerIntegrator2);
    vector<Vec3> positions(system.getNumParticles(), Vec3());
    innerContext1->setPositions(positions);
    innerContext2->setPositions(positions);

    // Create the kernel.

    kernel = context.getPlatform().createKernel(CalcATMForceKernel::Name(), context);
    kernel.getAs<CalcATMForceKernel>().initialize(context.getSystem(), owner);
}

double ATMForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << owner.getForceGroup())) == 0)
        return 0.0;

    ContextImpl& innerContextImpl1 = getContextImpl(*innerContext1);
    ContextImpl& innerContextImpl2 = getContextImpl(*innerContext2);

    //copies the coordinates etc. from the context to the inner contexts
    kernel.getAs<CalcATMForceKernel>().copyState(context, innerContextImpl1, innerContextImpl2);

    //evaluate variable energy and forces for original system
    double state1Energy = innerContextImpl1.calcForcesAndEnergy(true, true);

    //evaluate variable energy and force for the displaced system
    double state2Energy = innerContextImpl2.calcForcesAndEnergy(true, true);

    //evaluate the alchemical energy
    double energy = kernel.getAs<CalcATMForceKernel>().execute(context,
            innerContextImpl1, innerContextImpl2,
            state1Energy, state2Energy,
            includeForces, includeEnergy);

    //retrieve the perturbation energy
    perturbationEnergy = kernel.getAs<CalcATMForceKernel>().getPerturbationEnergy();

    return (includeEnergy ? energy : 0.0);
}

std::map<std::string, double> ATMForceImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    parameters[ATMForce::Lambda1()] = getOwner().getDefaultLambda1();
    parameters[ATMForce::Lambda2()] = getOwner().getDefaultLambda2();
    parameters[ATMForce::Alpha()] = getOwner().getDefaultAlpha();
    parameters[ATMForce::U0()] = getOwner().getDefaultU0();
    parameters[ATMForce::W0()] = getOwner().getDefaultW0();
    parameters[ATMForce::Umax()] = getOwner().getDefaultUmax();
    parameters[ATMForce::Ubcore()] = getOwner().getDefaultUbcore();
    parameters[ATMForce::Acore()] = getOwner().getDefaultAcore();
    parameters[ATMForce::Direction()] = getOwner().getDefaultDirection();
    return parameters;
}

std::vector<std::string> ATMForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcATMForceKernel::Name());
    return names;
}

vector<pair<int, int> > ATMForceImpl::getBondedParticles() const {
    vector<pair<int, int> > bonds;
    const ContextImpl& innerContextImpl = getContextImpl(*innerContext1);
    for (auto& impl : innerContextImpl.getForceImpls()) {
        for (auto& bond : impl->getBondedParticles())
            bonds.push_back(bond);
    }
    return bonds;
}

void ATMForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcATMForceKernel>().copyParametersToContext(context, owner);
}
