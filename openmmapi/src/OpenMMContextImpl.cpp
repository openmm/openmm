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

#include "Force.h"
#include "Integrator.h"
#include "OpenMMException.h"
#include "System.h"
#include "kernels.h"
#include "internal/ForceImpl.h"
#include "internal/OpenMMContextImpl.h"
#include <map>
#include <vector>

using namespace OpenMM;
using std::map;
using std::vector;
using std::string;

OpenMMContextImpl::OpenMMContextImpl(OpenMMContext& owner, System& system, Integrator& integrator, Platform* platform) :
				owner(owner), system(system), 
				integrator(integrator), platform(platform), 
				platformData(NULL),
				time(0)
{
    vector<string> kernelNames;
    kernelNames.push_back(CalcKineticEnergyKernel::Name());
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
        throw OpenMMException("Specified a Platform for an OpenMMContext which does not support all required kernels");
    platform->contextCreated(*this);
    kineticEnergyKernel = platform->createKernel(CalcKineticEnergyKernel::Name(), *this);
    vector<double> masses(system.getNumAtoms());
    for (size_t i = 0; i < masses.size(); ++i)
        masses[i] = system.getAtomMass(i);
    dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).initialize(system);
    for (size_t i = 0; i < forceImpls.size(); ++i)
        forceImpls[i]->initialize(*this);
    integrator.initialize(*this);
    positions = platform->createStream("atomPositions", system.getNumAtoms(), Stream::Double3, *this);
    velocities = platform->createStream("atomVelocities", system.getNumAtoms(), Stream::Double3, *this);
    forces = platform->createStream("atomForces", system.getNumAtoms(), Stream::Double3, *this);
    double zero[] = {0.0, 0.0, 0.0};
    velocities.fillWithValue(&zero);
}

OpenMMContextImpl::~OpenMMContextImpl() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        delete forceImpls[i];
    platform->contextDestroyed(*this);
}

double OpenMMContextImpl::getParameter(std::string name) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called getParameter() with invalid parameter name");
    return parameters[name];
}

void OpenMMContextImpl::setParameter(std::string name, double value) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called setParameter() with invalid parameter name");
    parameters[name] = value;
}

void OpenMMContextImpl::calcForces() {
    double zero[] = {0.0, 0.0, 0.0};
    forces.fillWithValue(zero);
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->calcForces(*this, forces);
}

double OpenMMContextImpl::calcKineticEnergy() {
    return dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).execute(*this);
}

double OpenMMContextImpl::calcPotentialEnergy() {
    double energy = 0.0;
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        energy += forceImpls[i]->calcEnergy(*this);
    return energy;
}

void OpenMMContextImpl::updateContextState() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->updateContextState(*this);
}

void OpenMMContextImpl::reinitialize() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        delete forceImpls[i];
    forceImpls.resize(0);
    for (int i = 0; i < system.getNumForces(); ++i) {
        forceImpls.push_back(system.getForce(i).createImpl());
        forceImpls[i]->initialize(*this);
    }
    integrator.initialize(*this);
}

void* OpenMMContextImpl::getPlatformData() {
    return platformData;
}

void OpenMMContextImpl::setPlatformData(void* data) {
    platformData = data;
}
