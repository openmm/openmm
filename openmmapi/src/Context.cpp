/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ForceImpl.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <cmath>
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

Context::Context(const System& system, Integrator& integrator, ContextImpl& linked) : properties(linked.getOwner().properties) {
    // This is used by ContextImpl::createLinkedContext().
    impl = new ContextImpl(*this, system, integrator, &linked.getPlatform(), properties, &linked);
    impl->initialize();
}

Context::Context(const System& system, Integrator& integrator) : properties(map<string, string>()) {
    impl = new ContextImpl(*this, system, integrator, 0, properties);
    impl->initialize();
}

Context::Context(const System& system, Integrator& integrator, Platform& platform) : properties(map<string, string>()) {
    impl = new ContextImpl(*this, system, integrator, &platform, properties);
    impl->initialize();
}

Context::Context(const System& system, Integrator& integrator, Platform& platform, const map<string, string>& properties) : properties(properties) {
    impl = new ContextImpl(*this, system, integrator, &platform, properties);
    impl->initialize();
}

Context::~Context() {
    delete impl;
}

const System& Context::getSystem() const {
    return impl->getSystem();

}

const Integrator& Context::getIntegrator() const {
    return impl->getIntegrator();
}

Integrator& Context::getIntegrator() {
    return impl->getIntegrator();
}

const Platform& Context::getPlatform() const {
    return impl->getPlatform();
}

Platform& Context::getPlatform() {
    return impl->getPlatform();
}

State Context::getState(int types, bool enforcePeriodicBox, int groups) const {
    State::StateBuilder builder(impl->getTime());
    Vec3 periodicBoxSize[3];
    impl->getPeriodicBoxVectors(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2]);
    builder.setPeriodicBoxVectors(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2]);
    bool includeForces = types&State::Forces;
    bool includeEnergy = types&State::Energy;
    bool includeParameterDerivs = types&State::ParameterDerivatives;
    if (includeForces || includeEnergy || includeParameterDerivs) {
        double energy = impl->calcForcesAndEnergy(includeForces || includeEnergy || includeParameterDerivs, includeEnergy, groups);
        if (includeEnergy)
            builder.setEnergy(impl->calcKineticEnergy(), energy);
        if (includeForces) {
            vector<Vec3> forces;
            impl->getForces(forces);
            builder.setForces(forces);
        }
    }
    if (types&State::Parameters) {
        map<string, double> params;
        for (auto& param : impl->parameters)
            params[param.first] = param.second;
        builder.setParameters(params);
    }
    if (types&State::ParameterDerivatives) {
        map<string, double> derivs;
        impl->getEnergyParameterDerivatives(derivs);
        builder.setEnergyParameterDerivatives(derivs);
    }
    if (types&State::Positions) {
        vector<Vec3> positions;
        impl->getPositions(positions);
        if (enforcePeriodicBox) {
            const vector<vector<int> >& molecules = impl->getMolecules();
            for (auto& mol : molecules) {
                // Find the molecule center.

                Vec3 center;
                for (int j : mol)
                    center += positions[j];
                center *= 1.0/mol.size();

                // Find the displacement to move it into the first periodic box.
                Vec3 diff;
                diff += periodicBoxSize[2]*floor(center[2]/periodicBoxSize[2][2]);
                diff += periodicBoxSize[1]*floor((center[1]-diff[1])/periodicBoxSize[1][1]);
                diff += periodicBoxSize[0]*floor((center[0]-diff[0])/periodicBoxSize[0][0]);

                // Translate all the particles in the molecule.
                for (int j : mol)
                    positions[j] -= diff;
            }
        }
        builder.setPositions(positions);
    }
    if (types&State::Velocities) {
        vector<Vec3> velocities;
        impl->getVelocities(velocities);
        builder.setVelocities(velocities);
    }
    return builder.getState();
}

void Context::setState(const State& state) {
    setTime(state.getTime());
    Vec3 a, b, c;
    state.getPeriodicBoxVectors(a, b, c);
    setPeriodicBoxVectors(a, b, c);
    if ((state.getDataTypes()&State::Positions) != 0)
        setPositions(state.getPositions());
    if ((state.getDataTypes()&State::Velocities) != 0)
        setVelocities(state.getVelocities());
    if ((state.getDataTypes()&State::Parameters) != 0)
        for (auto& param : state.getParameters())
            setParameter(param.first, param.second);
}

void Context::setTime(double time) {
    impl->setTime(time);
}

void Context::setPositions(const vector<Vec3>& positions) {
    if ((int) positions.size() != impl->getSystem().getNumParticles())
        throw OpenMMException("Called setPositions() on a Context with the wrong number of positions");
    impl->setPositions(positions);
}

void Context::setVelocities(const vector<Vec3>& velocities) {
    if ((int) velocities.size() != impl->getSystem().getNumParticles())
        throw OpenMMException("Called setVelocities() on a Context with the wrong number of velocities");
    impl->setVelocities(velocities);
}

void Context::setVelocitiesToTemperature(double temperature, int randomSeed) {
    const System& system = impl->getSystem();
    
    // Generate the list of Gaussian random numbers.
    
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(randomSeed, sfmt);
    vector<double> randoms;
    while (randoms.size() < system.getNumParticles()*3) {
        double x, y, r2;
        do {
            x = 2.0*genrand_real2(sfmt)-1.0;
            y = 2.0*genrand_real2(sfmt)-1.0;
            r2 = x*x + y*y;
        } while (r2 >= 1.0 || r2 == 0.0);
        double multiplier = sqrt((-2.0*log(r2))/r2);
        randoms.push_back(x*multiplier);
        randoms.push_back(y*multiplier);
    }
    
    // Assign the velocities.
    
    vector<Vec3> velocities(system.getNumParticles(), Vec3());
    int nextRandom = 0;
    for (int i = 0; i < system.getNumParticles(); i++) {
        double mass = system.getParticleMass(i);
        if (mass != 0) {
            double velocityScale = sqrt(BOLTZ*temperature/mass);
            velocities[i] = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*velocityScale;
        }
    }
    setVelocities(velocities);
    impl->applyVelocityConstraints(1e-5);
}

const map<string, double>& Context::getParameters() const {
    return impl->getParameters();
}

double Context::getParameter(const string& name) const {
    return impl->getParameter(name);
}

void Context::setParameter(const string& name, double value) {
    impl->setParameter(name, value);
}

void Context::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    impl->setPeriodicBoxVectors(a, b, c);
}

void Context::applyConstraints(double tol) {
    impl->applyConstraints(tol);
}

void Context::applyVelocityConstraints(double tol) {
    impl->applyVelocityConstraints(tol);
}

void Context::computeVirtualSites() {
    impl->computeVirtualSites();
}

void Context::reinitialize(bool preserveState) {
    const System& system = impl->getSystem();
    Integrator& integrator = impl->getIntegrator();
    Platform& platform = impl->getPlatform();
    stringstream checkpoint(ios_base::out | ios_base::in | ios_base::binary);
    if (preserveState)
        createCheckpoint(checkpoint);
    integrator.cleanup();
    delete impl;
    impl = new ContextImpl(*this, system, integrator, &platform, properties);
    impl->initialize();
    if (preserveState)
        loadCheckpoint(checkpoint);
}

void Context::createCheckpoint(ostream& stream) {
    impl->createCheckpoint(stream);
}

void Context::loadCheckpoint(istream& stream) {
    impl->loadCheckpoint(stream);
}

ContextImpl& Context::getImpl() {
    return *impl;
}

const ContextImpl& Context::getImpl() const {
    return *impl;
}

const vector<vector<int> >& Context::getMolecules() const {
    return impl->getMolecules();
}
