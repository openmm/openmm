/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/State.h"
#include "openmm/VirtualSite.h"
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
    if (system.getNumParticles() == 0)
        throw OpenMMException("Cannot create a Context for a System with no particles");
    
    // Check for errors in virtual sites and massless particles.
    
    for (int i = 0; i < system.getNumParticles(); i++) {
        if (system.isVirtualSite(i)) {
            if (system.getParticleMass(i) != 0.0)
                throw OpenMMException("Virtual site has nonzero mass");
            const VirtualSite& site = system.getVirtualSite(i);
            for (int j = 0; j < site.getNumParticles(); j++)
                if (system.isVirtualSite(site.getParticle(j)))
                    throw OpenMMException("A virtual site cannot depend on another virtual site");
        }
    }
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        if (system.getParticleMass(particle1) == 0.0 || system.getParticleMass(particle2) == 0.0)
            throw OpenMMException("A constraint cannot involve a massless particle");
    }
    
    // Find the list of kernels required.
    
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
    
    // Create and initialize kernels and other objects.
    
    platform->contextCreated(*this, properties);
    initializeForcesKernel = platform->createKernel(CalcForcesAndEnergyKernel::Name(), *this);
    dynamic_cast<CalcForcesAndEnergyKernel&>(initializeForcesKernel.getImpl()).initialize(system);
    kineticEnergyKernel = platform->createKernel(CalcKineticEnergyKernel::Name(), *this);
    dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).initialize(system);
    updateStateDataKernel = platform->createKernel(UpdateStateDataKernel::Name(), *this);
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).initialize(system);
    applyConstraintsKernel = platform->createKernel(ApplyConstraintsKernel::Name(), *this);
    dynamic_cast<ApplyConstraintsKernel&>(applyConstraintsKernel.getImpl()).initialize(system);
    virtualSitesKernel = platform->createKernel(VirtualSitesKernel::Name(), *this);
    dynamic_cast<VirtualSitesKernel&>(virtualSitesKernel.getImpl()).initialize(system);
    Vec3 periodicBoxVectors[3];
    system.getDefaultPeriodicBoxVectors(periodicBoxVectors[0], periodicBoxVectors[1], periodicBoxVectors[2]);
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setPeriodicBoxVectors(*this, periodicBoxVectors[0], periodicBoxVectors[1], periodicBoxVectors[2]);
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
    integrator.stateChanged(State::Positions);
}

void ContextImpl::getVelocities(std::vector<Vec3>& velocities) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getVelocities(*this, velocities);
}

void ContextImpl::setVelocities(const std::vector<Vec3>& velocities) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setVelocities(*this, velocities);
    integrator.stateChanged(State::Velocities);
}

void ContextImpl::getForces(std::vector<Vec3>& forces) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getForces(*this, forces);
}

const std::map<std::string, double>& ContextImpl::getParameters() const {
    return parameters;
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
    integrator.stateChanged(State::Parameters);
}

void ContextImpl::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) {
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).getPeriodicBoxVectors(*this, a, b, c);
}

void ContextImpl::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw OpenMMException("First periodic box vector must be parallel to x.");
    if (b[0] != 0.0 || b[2] != 0.0)
        throw OpenMMException("Second periodic box vector must be parallel to y.");
    if (c[0] != 0.0 || c[1] != 0.0)
        throw OpenMMException("Third periodic box vector must be parallel to z.");
    dynamic_cast<UpdateStateDataKernel&>(updateStateDataKernel.getImpl()).setPeriodicBoxVectors(*this, a, b, c);
}

void ContextImpl::applyConstraints(double tol) {
    dynamic_cast<ApplyConstraintsKernel&>(applyConstraintsKernel.getImpl()).apply(*this, tol);
}

void ContextImpl::computeVirtualSites() {
    dynamic_cast<VirtualSitesKernel&>(virtualSitesKernel.getImpl()).computePositions(*this);
}

double ContextImpl::calcForcesAndEnergy(bool includeForces, bool includeEnergy) {
    CalcForcesAndEnergyKernel& kernel = dynamic_cast<CalcForcesAndEnergyKernel&>(initializeForcesKernel.getImpl());
    double energy = 0.0;
    kernel.beginComputation(*this, includeForces, includeEnergy);
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        energy += forceImpls[i]->calcForcesAndEnergy(*this, includeForces, includeEnergy);
    energy += kernel.finishComputation(*this, includeForces, includeEnergy);
    return energy;
}

double ContextImpl::calcKineticEnergy() {
    return dynamic_cast<CalcKineticEnergyKernel&>(kineticEnergyKernel.getImpl()).execute(*this);
}

void ContextImpl::updateContextState() {
    for (int i = 0; i < (int) forceImpls.size(); ++i)
        forceImpls[i]->updateContextState(*this);
}

const vector<ForceImpl*>& ContextImpl::getForceImpls() const {
    return forceImpls;
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
    for (int i = 0; i < system.getNumParticles(); i++) {
        if (system.isVirtualSite(i)) {
            const VirtualSite& site = system.getVirtualSite(i);
            for (int j = 0; j < site.getNumParticles(); j++)
                bonds.push_back(std::make_pair(i, site.getParticle(j)));
        }
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
