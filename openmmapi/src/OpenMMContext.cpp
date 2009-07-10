/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "openmm/OpenMMContext.h"
#include "openmm/internal/OpenMMContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace std;

OpenMMContext::OpenMMContext(System& system, Integrator& integrator) {
    impl = new OpenMMContextImpl(*this, system, integrator, 0);
}

OpenMMContext::OpenMMContext(System& system, Integrator& integrator, Platform& platform) {
    impl = new OpenMMContextImpl(*this, system, integrator, &platform);
}

OpenMMContext::~OpenMMContext() {
	delete impl;
}

const System& OpenMMContext::getSystem() const {
    return impl->getSystem();

}

System& OpenMMContext::getSystem() {
    return impl->getSystem();
}

const Integrator& OpenMMContext::getIntegrator() const {
    return impl->getIntegrator();
}

Integrator& OpenMMContext::getIntegrator() {
    return impl->getIntegrator();
}

const Platform& OpenMMContext::getPlatform() const {
    return impl->getPlatform();
}

Platform& OpenMMContext::getPlatform() {
    return impl->getPlatform();
}

State OpenMMContext::getState(int types) const {
    State state(impl->getTime(), impl->getSystem().getNumParticles(), State::DataType(types));
    if (types&State::Energy)
        state.setEnergy(impl->calcKineticEnergy(), impl->calcPotentialEnergy());
    if (types&State::Forces) {
        impl->calcForces();
        impl->getForces().saveToArray(&state.updForces()[0]);
    }
    if (types&State::Parameters) {
        for (map<string, double>::const_iterator iter = impl->parameters.begin(); iter != impl->parameters.end(); iter++)
            state.updParameters()[iter->first] = iter->second;
    }
    if (types&State::Positions)
        impl->getPositions().saveToArray(&state.updPositions()[0]);
    if (types&State::Velocities)
        impl->getVelocities().saveToArray(&state.updVelocities()[0]);
    return state;
}

void OpenMMContext::setTime(double time) {
    impl->setTime(time);
}

void OpenMMContext::setPositions(const vector<Vec3>& positions) {
    if ((int) positions.size() != impl->getSystem().getNumParticles())
        throw OpenMMException("Called setPositions() on an OpenMMContext with the wrong number of positions");
    impl->getPositions().loadFromArray(&positions[0]);
}

void OpenMMContext::setVelocities(const vector<Vec3>& velocities) {
    if ((int) velocities.size() != impl->getSystem().getNumParticles())
        throw OpenMMException("Called setVelocities() on an OpenMMContext with the wrong number of velocities");
    impl->getVelocities().loadFromArray(&velocities[0]);
}

double OpenMMContext::getParameter(const string& name) {
    return impl->getParameter(name);
}

void OpenMMContext::setParameter(const string& name, double value) {
    impl->setParameter(name, value);
}

void OpenMMContext::reinitialize() {
    System& system = impl->getSystem();
    Integrator& integrator = impl->getIntegrator();
    Platform& platform = impl->getPlatform();
    delete impl;
    impl = new OpenMMContextImpl(*this, system, integrator, &platform);
}
