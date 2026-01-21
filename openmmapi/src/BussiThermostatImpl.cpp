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

#include "openmm/internal/BussiThermostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include <vector>
#include <set>

using namespace OpenMM;
using std::vector;
using std::set;

BussiThermostatImpl::BussiThermostatImpl(const BussiThermostat& owner) : owner(owner), stepCount(0) {
}

void BussiThermostatImpl::initialize(ContextImpl& context) {
    if (owner.getDefaultTemperature() < 0)
        throw OpenMMException("BussiThermostat: temperature cannot be negative");
    if (owner.getDefaultTau() <= 0)
        throw OpenMMException("BussiThermostat: tau must be positive");
    
    // Determine which particles to apply the thermostat to
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    
    if (owner.getApplyToAllParticles()) {
        // Apply to all particles with non-zero mass
        for (int i = 0; i < numParticles; i++) {
            if (system.getParticleMass(i) > 0)
                particleIndices.push_back(i);
        }
    } else {
        // Use the specified particle set
        const set<int>& particles = owner.getParticles();
        for (int particle : particles) {
            if (particle < 0 || particle >= numParticles)
                throw OpenMMException("BussiThermostat: particle index out of range");
            if (system.getParticleMass(particle) > 0)
                particleIndices.push_back(particle);
        }
    }
    
    if (particleIndices.empty())
        throw OpenMMException("BussiThermostat: no particles to apply thermostat to");
    
    // Create the kernel
    kernel = context.getPlatform().createKernel(ApplyBussiThermostatKernel::Name(), context);
    kernel.getAs<ApplyBussiThermostatKernel>().initialize(context.getSystem(), owner, particleIndices);
}

void BussiThermostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    // Only apply the thermostat at the specified frequency
    if (stepCount % owner.getFrequency() == 0) {
        kernel.getAs<ApplyBussiThermostatKernel>().execute(context);
    }
    stepCount++;
}

std::map<std::string, double> BussiThermostatImpl::getDefaultParameters() {
    return {
        {BussiThermostat::Temperature(), owner.getDefaultTemperature()},
        {BussiThermostat::Tau(), owner.getDefaultTau()},
        {BussiThermostat::ReservoirEnergyTranslational(), 0.0},
        {BussiThermostat::ReservoirEnergyRotational(), 0.0}
    };
}

std::vector<std::string> BussiThermostatImpl::getKernelNames() {
    return {ApplyBussiThermostatKernel::Name()};
}
