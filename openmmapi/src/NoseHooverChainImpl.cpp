/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
 * Authors: Andreas Krämer and Andrew C. Simmonett                            *
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

#include "openmm/internal/NoseHooverChainImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <string>
#include <vector>

using namespace OpenMM;
using std::vector;

NoseHooverChainImpl::NoseHooverChainImpl(const NoseHooverChain& owner) : owner(owner) {
}

void NoseHooverChainImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(PropagateNoseHooverChainKernel::Name(), context);
    kernel.getAs<PropagateNoseHooverChainKernel>().initialize(context.getSystem(), owner);
}

void NoseHooverChainImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    // This kernel updates the Nose-Hoover particles when invoked explicitly
    // via its propagate(), which is called from the integrator.
}

double NoseHooverChainImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    // This kernel doesn't compute an energy or force.  Instead it is invoked
    // directly by the integrator via it's propagate() method.
    return 0;
}

std::map<std::string, double> NoseHooverChainImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    const auto &owner = getOwner();
    double frequency = owner.getDefaultCollisionFrequency();
    int chainLength = owner.getDefaultChainLength();
    double T = owner.getDefaultTemperature();
    double kT = BOLTZ * T;
    int DOFs = owner.getDefaultNumDegreesOfFreedom();
    parameters[owner.Temperature()] = T;
    parameters[owner.CollisionFrequency()] = frequency;
    parameters[owner.ChainLength()] = chainLength;
    parameters[owner.NumYoshidaSuzukiTimeSteps()] = owner.getDefaultNumYoshidaSuzukiTimeSteps();
    parameters[owner.NumMultiTimeSteps()] = owner.getDefaultNumMultiTimeSteps();
    parameters[owner.NumDegreesOfFreedom()] = DOFs;
    for(int i = 0; i < chainLength; ++i) {
        parameters[owner.Force(i)] = 0;
        parameters[owner.Position(i)] = 0;
        parameters[owner.Mass(i)] = kT / (frequency * frequency);
    }
    parameters[owner.Mass(0)] *= DOFs;

    // Set the velocities to the appropriate Boltzmann distribution;
    // this is copied from Context::setVelocitiesToTemperature()
    OpenMM_SFMT::SFMT sfmt;
    int randomSeed = 0; //TODO figure out where / how this should be handled
    init_gen_rand(randomSeed, sfmt);
    vector<double> randoms;
    while (randoms.size() < chainLength) {
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
    int nextRandom = 0;
    for (int i = 0; i < chainLength; i++) {
        double velocityScale = 1 / frequency; // = sqrt(kT / Q) = sqrt(kT / [ kT frequency^-2])
        parameters[owner.Velocity(i)] = velocityScale * randoms[nextRandom++];
    }
    // N.B. the zeroth entry is computed as a function of the instantaneous KE at the start of propagate
    for(int i = 1; i < chainLength-1; ++i) {
        const double & v = parameters[owner.Velocity(i)];
        parameters[owner.Force(i+1)] = (parameters[owner.Mass(i)] * v * v - kT) / parameters[owner.Mass(i+1)];
    }


    return parameters;
}

std::vector<std::string> NoseHooverChainImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(PropagateNoseHooverChainKernel::Name());
    return names;
}
