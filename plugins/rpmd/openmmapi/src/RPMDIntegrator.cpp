/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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

#include "openmm/RPMDIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/RpmdKernels.h"
#include "openmm/RPMDUpdater.h"
#include "SimTKOpenMMRealType.h"
#include <cmath>
#include <string>

using namespace OpenMM;
using namespace std;

RPMDIntegrator::RPMDIntegrator(int numCopies, double temperature, double frictionCoeff, double stepSize, const map<int, int>& contractions) :
        numCopies(numCopies), applyThermostat(true), contractions(contractions), forcesAreValid(false), hasSetPosition(false), hasSetVelocity(false), isFirstStep(true) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
}

RPMDIntegrator::RPMDIntegrator(int numCopies, double temperature, double frictionCoeff, double stepSize) :
        numCopies(numCopies), applyThermostat(true), forcesAreValid(false), hasSetPosition(false), hasSetVelocity(false), isFirstStep(true) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
}

void RPMDIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    if (contextRef.getSystem().getNumConstraints() > 0)
        throw OpenMMException("RPMDIntegrator cannot be used with Systems that include constraints");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateRPMDStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateRPMDStepKernel>().initialize(contextRef.getSystem(), *this);
}

void RPMDIntegrator::cleanup() {
    kernel = Kernel();
}

void RPMDIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> RPMDIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateRPMDStepKernel::Name());
    return names;
}

void RPMDIntegrator::setPositions(int copy, const vector<Vec3>& positions) {
    kernel.getAs<IntegrateRPMDStepKernel>().setPositions(copy, positions);
    forcesAreValid = false;
    hasSetPosition = true;
}

void RPMDIntegrator::setVelocities(int copy, const vector<Vec3>& velocities) {
    kernel.getAs<IntegrateRPMDStepKernel>().setVelocities(copy, velocities);
    hasSetVelocity = true;
}

State RPMDIntegrator::getState(int copy, int types, bool enforcePeriodicBox, int groups) {
    if (isFirstStep) {
        // Call setPositions() on the Context so it doesn't think the user is trying to
        // run a simulation without setting positions first.  These positions will
        // immediately get overwritten by the ones stored in this integrator.

        vector<Vec3> p(context->getSystem().getNumParticles(), Vec3());
        context->getOwner().setPositions(p);
        isFirstStep = false;
    }
    kernel.getAs<IntegrateRPMDStepKernel>().copyToContext(copy, *context);
    State state = context->getOwner().getState(types, enforcePeriodicBox && copy == 0, groups);
    if (enforcePeriodicBox && copy > 0 && (types&State::Positions) != 0) {
        // Apply periodic boundary conditions based on copy 0.  Otherwise, molecules might end
        // up in different places for different copies.

        kernel.getAs<IntegrateRPMDStepKernel>().copyToContext(0, *context);
        State state2 = context->getOwner().getState(State::Positions, false, groups);
        vector<Vec3> positions = state.getPositions();
        const vector<Vec3>& refPos = state2.getPositions();
        const vector<vector<int> >& molecules = context->getMolecules();
        Vec3 periodicBoxSize[3];
        state2.getPeriodicBoxVectors(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2]);
        for (auto& mol : molecules) {
            // Find the molecule center.

            Vec3 center;
            for (int j : mol)
                center += refPos[j];
            center *= 1.0/mol.size();

            // Find the displacement to move it into the first periodic box.
            Vec3 diff;
            diff += periodicBoxSize[2]*floor(center[2]/periodicBoxSize[2][2]);
            diff += periodicBoxSize[1]*floor((center[1]-diff[1])/periodicBoxSize[1][1]);
            diff += periodicBoxSize[0]*floor((center[0]-diff[0])/periodicBoxSize[0][0]);

            // Translate all the particles in the molecule.
            for (int j : mol) {
                Vec3& pos = positions[j];
                pos -= diff;
            }
        }

        // Construct the new State.

        State::StateBuilder builder(state.getTime());
        builder.setPositions(positions);
        builder.setPeriodicBoxVectors(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2]);
        if (types&State::Velocities)
            builder.setVelocities(state.getVelocities());
        if (types&State::Forces)
            builder.setForces(state.getForces());
        if (types&State::Parameters)
            builder.setParameters(state.getParameters());
        if (types&State::Energy)
            builder.setEnergy(state.getKineticEnergy(), state.getPotentialEnergy());
        state = builder.getState();
    }
    return state;
}

double RPMDIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateRPMDStepKernel>().computeKineticEnergy(*context, *this);
}

void RPMDIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    if (!hasSetPosition) {
        // Initialize the positions from the context.

        State s = context->getOwner().getState(State::Positions);
        for (int i = 0; i < numCopies; i++)
            setPositions(i, s.getPositions());
    }
    if (!hasSetVelocity) {
        // Initialize the velocities from the context.

        State s = context->getOwner().getState(State::Velocities);
        for (int i = 0; i < numCopies; i++)
            setVelocities(i, s.getVelocities());
    }
    if (isFirstStep) {
        // Call setPositions() on the Context so it doesn't think the user is trying to
        // run a simulation without setting positions first.  These positions will
        // immediately get overwritten by the ones stored in this integrator.

        vector<Vec3> p(context->getSystem().getNumParticles(), Vec3());
        context->getOwner().setPositions(p);
        isFirstStep = false;
    }
    for (auto impl : context->getForceImpls()) {
        RPMDUpdater* updater = dynamic_cast<RPMDUpdater*>(impl);
        if (updater != NULL)
            updater->updateRPMDState(*context);
    }
    for (int i = 0; i < steps; ++i) {
        kernel.getAs<IntegrateRPMDStepKernel>().execute(*context, *this, forcesAreValid);
        forcesAreValid = true;
    }
}

double RPMDIntegrator::getTotalEnergy() {
    const System& system = owner->getSystem();
    int numParticles = system.getNumParticles();
    double energy = 0.0;
    const double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    const double wn = numCopies*BOLTZ*temperature/hbar;
    State prevState = getState(numCopies-1, State::Positions);
    for (int i = 0; i < numCopies; i++) {
        // Add the energy of this copy.

        State state = getState(i, State::Positions | State::Energy);
        energy += state.getKineticEnergy()+state.getPotentialEnergy();

        // Add the energy from the springs connecting it to the previous copy.

        for (int j = 0; j < numParticles; j++) {
            Vec3 delta = state.getPositions()[j]-prevState.getPositions()[j];
            energy += 0.5*wn*wn*system.getParticleMass(j)*delta.dot(delta);
        }
        prevState = state;
    }
    return energy;
}
