/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/Force.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include <map>
#include <vector>

using namespace OpenMM;
using std::map;
using std::vector;
using std::string;

ContextImpl::ContextImpl(Context& owner, System& system, Integrator& integrator, Platform* platform) :
				owner(owner), system(system), 
				integrator(integrator), platform(platform), 
				platformData(NULL)
{
    vector<string> kernelNames;
    kernelNames.push_back(CalcKineticEnergyKernel::Name());
    kernelNames.push_back(InitializeForcesKernel::Name());
    kernelNames.push_back(UpdateTimeKernel::Name());
    for (int i = 0; i < system.getNumForces(); ++i) {
        forceImpls.push_back(system.getForce(i).createImpl());
        map<string, double> forceParameters = forceImpls[forceImpls.size()-1]->getDefaultParameters();
        parameters.insert(forceParameters.begin(), forceParameters.end());
        vector<string> forceKernels = forceImpls[forceImpls.size()-1]->getKernelNames();
        kernelNames.insert(kernelNames.begin(), forceKernels.begin(), forceKernels.end());
    }
    vector<string> integratorKernels = integrator.getKernelNames();
    kernelNames.insert(kernelNames.begin(), integratorKernels.begin(), integratorKernels.end());
    if (platform == 0)
        this->platform = platform = &Platform::findPlatform(kernelNames);
    else if (!platform->supportsKernels(kernelNames))
        throw OpenMMException("Specified a Platform for a Context which does not support all required kernels");
    platform->contextCreated(*this);
    initializeForcesKernel = platform->createKernel(InitializeForcesKernel::Name(), *this);
    dynamic_cast<InitializeForcesKernel&>(initializeForcesKernel.getImpl()).initialize(system);
    kineticEnergyKernel = platform->createKernel(CalcKineticEnergyKernel::Name(), *this);
    dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).initialize(system);
    updateTimeKernel = platform->createKernel(UpdateTimeKernel::Name(), *this);
    dynamic_cast<UpdateTimeKernel&>(updateTimeKernel.getImpl()).initialize(system);
    for (size_t i = 0; i < forceImpls.size(); ++i)
        forceImpls[i]->initialize(*this);
    integrator.initialize(*this);
    positions = platform->createStream("particlePositions", system.getNumParticles(), Stream::Double3, *this);
    velocities = platform->createStream("particleVelocities", system.getNumParticles(), Stream::Double3, *this);
    forces = platform->createStream("particleForces", system.getNumParticles(), Stream::Double3, *this);
    double zero[] = {0.0, 0.0, 0.0};
    velocities.fillWithValue(&zero);
}

ContextImpl::~ContextImpl() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        delete forceImpls[i];
    platform->contextDestroyed(*this);
}

double ContextImpl::getTime() const {
    return dynamic_cast<const UpdateTimeKernel&>(updateTimeKernel.getImpl()).getTime(*this);
}

void ContextImpl::setTime(double t) {
    dynamic_cast<UpdateTimeKernel&>(updateTimeKernel.getImpl()).setTime(*this, t);
}

double ContextImpl::getParameter(std::string name) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called getParameter() with invalid parameter name");
    return parameters[name];
}

void ContextImpl::setParameter(std::string name, double value) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called setParameter() with invalid parameter name");
    parameters[name] = value;
}

void ContextImpl::calcForces() {
    dynamic_cast<InitializeForcesKernel&>(initializeForcesKernel.getImpl()).execute(*this);
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->calcForces(*this, forces);
}

double ContextImpl::calcKineticEnergy() {
    return dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).execute(*this);
}

double ContextImpl::calcPotentialEnergy() {
    double energy = 0.0;
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        energy += forceImpls[i]->calcEnergy(*this);
    return energy;
}

void ContextImpl::updateContextState() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->updateContextState(*this);
}

void* ContextImpl::getPlatformData() {
    return platformData;
}

const void* ContextImpl::getPlatformData() const {
    return platformData;
}

void ContextImpl::setPlatformData(void* data) {
    platformData = data;
}
