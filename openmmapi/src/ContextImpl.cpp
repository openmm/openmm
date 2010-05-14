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
#include <utility>
#include <vector>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::string;

ContextImpl::ContextImpl(Context& owner, System& system, Integrator& integrator, Platform* platform, const map<string, string>& properties) :
         owner(owner), system(system), integrator(integrator), hasInitializedForces(false), platform(platform), platformData(NULL) {
    vector<string> kernelNames;
    kernelNames.push_back(CalcKineticEnergyKernel::Name());
    kernelNames.push_back(CalcForcesAndEnergyKernel::Name());
    kernelNames.push_back(UpdateStateDataKernel::Name());
    for (int i = 0; i < system.getNumForces(); ++i) {
        forceImpls.push_back(system.getForce(i).createImpl());
        map<string, double> forceParameters = forceImpls[forceImpls.size()-1]->getDefaultParameters();
        parameters.insert(forceParameters.begin(), forceParameters.end());
        vector<string> forceKernels = forceImpls[forceImpls.size()-1]->getKernelNames();
        kernelNames.insert(kernelNames.begin(), forceKernels.begin(), forceKernels.end());
    }
    hasInitializedForces = true;
    vector<string> integratorKernels = integrator.getKernelNames();
    kernelNames.insert(kernelNames.begin(), integratorKernels.begin(), integratorKernels.end());
    if (platform == 0)
        this->platform = platform = &Platform::findPlatform(kernelNames);
    else if (!platform->supportsKernels(kernelNames))
        throw OpenMMException("Specified a Platform for a Context which does not support all required kernels");
    platform->contextCreated(*this, properties);
    initializeForcesKernel = platform->createKernel(CalcForcesAndEnergyKernel::Name(), *this);
    dynamic_cast<CalcForcesAndEnergyKernel&>(initializeForcesKernel.getImpl()).initialize(system);
    kineticEnergyKernel = platform->createKernel(CalcKineticEnergyKernel::Name(), *this);
    dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).initialize(system);
    updateStateDataKernel = platform->createKernel(UpdateStateDataKernel::Name(), *this);
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).initialize(system);
    for (size_t i = 0; i < forceImpls.size(); ++i)
        forceImpls[i]->initialize(*this);
    integrator.initialize(*this);
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setVelocities(*this, vector<Vec3>(system.getNumParticles()));
}

ContextImpl::~ContextImpl() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        delete forceImpls[i];
    platform->contextDestroyed(*this);
}

double ContextImpl::getTime() const {
    return dynamic_cast<const UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getTime(*this);
}

void ContextImpl::setTime(double t) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setTime(*this, t);
}

void ContextImpl::getPositions(std::vector<Vec3>& positions) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getPositions(*this, positions);
}

void ContextImpl::setPositions(const std::vector<Vec3>& positions) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setPositions(*this, positions);
}

void ContextImpl::getVelocities(std::vector<Vec3>& velocities) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getVelocities(*this, velocities);
}

void ContextImpl::setVelocities(const std::vector<Vec3>& velocities) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setVelocities(*this, velocities);
}

void ContextImpl::getForces(std::vector<Vec3>& forces) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getForces(*this, forces);
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
    CalcForcesAndEnergyKernel& kernel = dynamic_cast<CalcForcesAndEnergyKernel&>(initializeForcesKernel.getImpl());
    kernel.beginForceComputation(*this);
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->calcForces(*this);
    kernel.finishForceComputation(*this);
}

double ContextImpl::calcKineticEnergy() {
    return dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).execute(*this);
}

double ContextImpl::calcPotentialEnergy() {
    CalcForcesAndEnergyKernel& kernel = dynamic_cast<CalcForcesAndEnergyKernel&>(initializeForcesKernel.getImpl());
    kernel.beginEnergyComputation(*this);
    double energy = 0.0;
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        energy += forceImpls[i]->calcEnergy(*this);
    energy += kernel.finishEnergyComputation(*this);
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

const vector<vector<int> >& ContextImpl::getMolecules() const {
    if (!hasInitializedForces)
        throw OpenMMException("ContextImpl: getMolecules() cannot be called until all ForceImpls have been initialized");
    if (molecules.size() > 0 || system.getNumParticles() == 0)
        return molecules;

    // First make a list of bonds and constraints.

    vector<pair<int, int> > bonds;
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        bonds.push_back(std::make_pair(particle1, particle2));
    }
    for (int i = 0; i < (int) forceImpls.size(); i++) {
        vector<pair<int, int> > forceBonds = forceImpls[i]->getBondedParticles();
        bonds.insert(bonds.end(), forceBonds.begin(), forceBonds.end());
    }

    // Make a list of every other particle to which each particle is connected

    int numParticles = system.getNumParticles();
    vector<vector<int> > particleBonds(numParticles);
    for (int i = 0; i < (int) bonds.size(); i++) {
        particleBonds[bonds[i].first].push_back(bonds[i].second);
        particleBonds[bonds[i].second].push_back(bonds[i].first);
    }

    // Now tag particles by which molecule they belong to.

    vector<int> particleMolecule(numParticles, -1);
    int numMolecules = 0;
    for (int i = 0; i < numParticles; i++)
        if (particleMolecule[i] == -1)
            tagParticlesInMolecule(i, numMolecules++, particleMolecule, particleBonds);
    molecules.resize(numMolecules);
    for (int i = 0; i < numParticles; i++)
        molecules[particleMolecule[i]].push_back(i);
    return molecules;
}

void ContextImpl::tagParticlesInMolecule(int particle, int molecule, vector<int>& particleMolecule, vector<vector<int> >& particleBonds) {
    // Recursively tag particles as belonging to a particular molecule.

    particleMolecule[particle] = molecule;
    for (int i = 0; i < (int) particleBonds[particle].size(); i++)
        if (particleMolecule[particleBonds[particle][i]] == -1)
            tagParticlesInMolecule(particleBonds[particle][i], molecule, particleMolecule, particleBonds);
}
