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

#include "openmm/Force.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/kernels.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/State.h"
#include "openmm/VirtualSite.h"
#include "openmm/Context.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <string.h>

using namespace OpenMM;
using namespace std;
const static char CHECKPOINT_MAGIC_BYTES[] = "OpenMM Binary Checkpoint\n";


ContextImpl::ContextImpl(Context& owner, const System& system, Integrator& integrator, Platform* platform, const map<string, string>& properties, ContextImpl* originalContext) :
        owner(owner), system(system), integrator(integrator), hasInitializedForces(false), hasSetPositions(false), integratorIsDeleted(false),
        lastForceGroups(-1), platform(platform), platformData(NULL) {
    int numParticles = system.getNumParticles();
    if (numParticles == 0)
        throw OpenMMException("Cannot create a Context for a System with no particles");
    
    // Check for errors in virtual sites and massless particles.
    
    for (int i = 0; i < numParticles; i++) {
        if (system.isVirtualSite(i)) {
            if (system.getParticleMass(i) != 0.0)
                throw OpenMMException("Virtual site has nonzero mass");
            const VirtualSite& site = system.getVirtualSite(i);
            for (int j = 0; j < site.getNumParticles(); j++)
                if (system.isVirtualSite(site.getParticle(j)))
                    throw OpenMMException("A virtual site cannot depend on another virtual site");
        }
    }
    set<pair<int, int> > constraintAtoms;
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        if (particle1 == particle2)
            throw OpenMMException("A constraint cannot connect a particle to itself");
        if (particle1 < 0 || particle2 < 0 || particle1 >= numParticles || particle2 >= numParticles)
            throw OpenMMException("Illegal particle index in constraint");
        double mass1 = system.getParticleMass(particle1);
        double mass2 = system.getParticleMass(particle2);
        if ((mass1 == 0.0 && mass2 != 0.0) || (mass2 == 0.0 && mass1 != 0.0))
            throw OpenMMException("A constraint cannot involve a massless particle");
        pair<int, int> atoms = make_pair(min(particle1, particle2), max(particle1, particle2));
        if (constraintAtoms.find(atoms) != constraintAtoms.end())
            throw OpenMMException("The System has two constraints between the same atoms.  This will produce a singular constraint matrix.");
        constraintAtoms.insert(atoms);
    }
    
    // Validate the list of properties.

    const vector<string>& platformProperties = platform->getPropertyNames();
    map<string, string> validatedProperties;
    for (auto& prop : properties) {
        string property = prop.first;
        if (platform->deprecatedPropertyReplacements.find(property) != platform->deprecatedPropertyReplacements.end())
            property = platform->deprecatedPropertyReplacements[property];
        bool valid = false;
        for (auto& p : platformProperties)
            if (p == property) {
                valid = true;
                break;
            }
        if (!valid)
            throw OpenMMException("Illegal property name: "+prop.first);
        validatedProperties[property] = prop.second;
    }
    
    // Find the list of kernels required.
    
    vector<string> kernelNames;
    kernelNames.push_back(CalcForcesAndEnergyKernel::Name());
    kernelNames.push_back(UpdateStateDataKernel::Name());
    kernelNames.push_back(ApplyConstraintsKernel::Name());
    kernelNames.push_back(VirtualSitesKernel::Name());
    for (int i = 0; i < system.getNumForces(); ++i) {
        forceImpls.push_back(system.getForce(i).createImpl());
        vector<string> forceKernels = forceImpls[forceImpls.size()-1]->getKernelNames();
        kernelNames.insert(kernelNames.begin(), forceKernels.begin(), forceKernels.end());
    }
    hasInitializedForces = true;
    vector<string> integratorKernels = integrator.getKernelNames();
    kernelNames.insert(kernelNames.begin(), integratorKernels.begin(), integratorKernels.end());
    
    // Select a platform to use.
    
    vector<pair<double, Platform*> > candidatePlatforms;
    if (platform == NULL) {
        char* defaultPlatform = getenv("OPENMM_DEFAULT_PLATFORM");
        if (defaultPlatform != NULL)
            platform = &Platform::getPlatformByName(string(defaultPlatform));
    }
    if (platform == NULL) {
        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
            Platform& p = Platform::getPlatform(i);
            if (p.supportsKernels(kernelNames))
                candidatePlatforms.push_back(make_pair(p.getSpeed(), &p));
        }
        if (candidatePlatforms.size() == 0)
            throw OpenMMException("No Platform supports all the requested kernels");
        sort(candidatePlatforms.begin(), candidatePlatforms.end());
    }
    else {
        if (!platform->supportsKernels(kernelNames))
            throw OpenMMException("Specified a Platform for a Context which does not support all required kernels");
        candidatePlatforms.push_back(make_pair(platform->getSpeed(), platform));
    }
    for (int i = candidatePlatforms.size()-1; i >= 0; i--) {
        try {
            this->platform = platform = candidatePlatforms[i].second;
            if (originalContext == NULL)
                platform->contextCreated(*this, validatedProperties);
            else
                platform->linkedContextCreated(*this, *originalContext);
            break;
        }
        catch (...) {
            if (i > 0)
                continue;
            throw;
        }
    }
}

void ContextImpl::initialize() {
    // Create and initialize kernels and other objects.
    
    initializeForcesKernel = platform->createKernel(CalcForcesAndEnergyKernel::Name(), *this);
    initializeForcesKernel.getAs<CalcForcesAndEnergyKernel>().initialize(system);
    updateStateDataKernel = platform->createKernel(UpdateStateDataKernel::Name(), *this);
    updateStateDataKernel.getAs<UpdateStateDataKernel>().initialize(system);
    applyConstraintsKernel = platform->createKernel(ApplyConstraintsKernel::Name(), *this);
    applyConstraintsKernel.getAs<ApplyConstraintsKernel>().initialize(system);
    virtualSitesKernel = platform->createKernel(VirtualSitesKernel::Name(), *this);
    virtualSitesKernel.getAs<VirtualSitesKernel>().initialize(system);
    Vec3 periodicBoxVectors[3];
    system.getDefaultPeriodicBoxVectors(periodicBoxVectors[0], periodicBoxVectors[1], periodicBoxVectors[2]);
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setPeriodicBoxVectors(*this, periodicBoxVectors[0], periodicBoxVectors[1], periodicBoxVectors[2]);
    for (size_t i = 0; i < forceImpls.size(); ++i) {
        forceImpls[i]->initialize(*this);
        map<string, double> forceParameters = forceImpls[i]->getDefaultParameters();
        parameters.insert(forceParameters.begin(), forceParameters.end());
    }
    integrator.initialize(*this);
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setVelocities(*this, vector<Vec3>(system.getNumParticles()));
}

ContextImpl::~ContextImpl() {
    for (auto force : forceImpls)
        delete force;
    
    // Make sure all kernels get properly deleted before contextDestroyed() is called.
    
    initializeForcesKernel = Kernel();
    updateStateDataKernel = Kernel();
    applyConstraintsKernel = Kernel();
    virtualSitesKernel = Kernel();
    if (!integratorIsDeleted) {
        // The Context is being deleted before the Integrator, so call cleanup() on it now.
        
        integrator.cleanup();
        integrator.context = NULL;
    }
    platform->contextDestroyed(*this);
}

double ContextImpl::getTime() const {
    return updateStateDataKernel.getAs<const UpdateStateDataKernel>().getTime(*this);
}

void ContextImpl::setTime(double t) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setTime(*this, t);
}

void ContextImpl::getPositions(std::vector<Vec3>& positions) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().getPositions(*this, positions);
}

void ContextImpl::setPositions(const std::vector<Vec3>& positions) {
    hasSetPositions = true;
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setPositions(*this, positions);
    integrator.stateChanged(State::Positions);
}

void ContextImpl::getVelocities(std::vector<Vec3>& velocities) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().getVelocities(*this, velocities);
}

void ContextImpl::setVelocities(const std::vector<Vec3>& velocities) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setVelocities(*this, velocities);
    integrator.stateChanged(State::Velocities);
}

void ContextImpl::getForces(std::vector<Vec3>& forces) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().getForces(*this, forces);
}

const std::map<std::string, double>& ContextImpl::getParameters() const {
    return parameters;
}

double ContextImpl::getParameter(std::string name) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called getParameter() with invalid parameter name: "+name);
    return parameters[name];
}

void ContextImpl::setParameter(std::string name, double value) {
    if (parameters.find(name) == parameters.end())
        throw OpenMMException("Called setParameter() with invalid parameter name: "+name);
    parameters[name] = value;
    integrator.stateChanged(State::Parameters);
}

void ContextImpl::getEnergyParameterDerivatives(std::map<std::string, double>& derivs) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().getEnergyParameterDerivatives(*this, derivs);
}

void ContextImpl::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) {
    updateStateDataKernel.getAs<UpdateStateDataKernel>().getPeriodicBoxVectors(*this, a, b, c);
}

void ContextImpl::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw OpenMMException("First periodic box vector must be parallel to x.");
    if (b[2] != 0.0)
        throw OpenMMException("Second periodic box vector must be in the x-y plane.");
    if (a[0] <= 0.0 || b[1] <= 0.0 || c[2] <= 0.0 || a[0] < 2*fabs(b[0]) || a[0] < 2*fabs(c[0]) || b[1] < 2*fabs(c[1]))
        throw OpenMMException("Periodic box vectors must be in reduced form.");
    updateStateDataKernel.getAs<UpdateStateDataKernel>().setPeriodicBoxVectors(*this, a, b, c);
}

void ContextImpl::applyConstraints(double tol) {
    if (!hasSetPositions)
        throw OpenMMException("Particle positions have not been set");
    applyConstraintsKernel.getAs<ApplyConstraintsKernel>().apply(*this, tol);
}

void ContextImpl::applyVelocityConstraints(double tol) {
    if (!hasSetPositions)
        throw OpenMMException("Particle positions have not been set");
    applyConstraintsKernel.getAs<ApplyConstraintsKernel>().applyToVelocities(*this, tol);
}

void ContextImpl::computeVirtualSites() {
    virtualSitesKernel.getAs<VirtualSitesKernel>().computePositions(*this);
}

double ContextImpl::calcForcesAndEnergy(bool includeForces, bool includeEnergy, int groups) {
    if (!hasSetPositions)
        throw OpenMMException("Particle positions have not been set");
    lastForceGroups = groups;
    CalcForcesAndEnergyKernel& kernel = initializeForcesKernel.getAs<CalcForcesAndEnergyKernel>();
    while (true) {
        double energy = 0.0;
        kernel.beginComputation(*this, includeForces, includeEnergy, groups);
        for (auto force : forceImpls)
            energy += force->calcForcesAndEnergy(*this, includeForces, includeEnergy, groups);
        bool valid = true;
        energy += kernel.finishComputation(*this, includeForces, includeEnergy, groups, valid);
        if (valid)
            return energy;
    }
}

int& ContextImpl::getLastForceGroups() {
    return lastForceGroups;
}

double ContextImpl::calcKineticEnergy() {
    return integrator.computeKineticEnergy();
}

bool ContextImpl::updateContextState() {
    bool forcesInvalid = false;
    for (auto force : forceImpls)
        force->updateContextState(*this, forcesInvalid);
    return forcesInvalid;
}

const vector<ForceImpl*>& ContextImpl::getForceImpls() const {
    return forceImpls;
}

vector<ForceImpl*>& ContextImpl::getForceImpls() {
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
    for (auto force : forceImpls) {
        vector<pair<int, int> > forceBonds = force->getBondedParticles();
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
    for (auto& bond : bonds) {
        particleBonds[bond.first].push_back(bond.second);
        particleBonds[bond.second].push_back(bond.first);
    }

    // Now identify particles by which molecule they belong to.

    molecules = findMolecules(numParticles, particleBonds);
    return molecules;
}

vector<vector<int> > ContextImpl::findMolecules(int numParticles, vector<vector<int> >& particleBonds) {
    // This is essentially a recursive algorithm, but it is reformulated as a loop to avoid
    // stack overflows.  It selects a particle, marks it as a new molecule, then recursively
    // marks every particle bonded to it as also being in that molecule.
    
    vector<int> particleMolecule(numParticles, -1);
    int numMolecules = 0;
    for (int i = 0; i < numParticles; i++)
        if (particleMolecule[i] == -1) {
            // Start a new molecule.
            
            vector<int> particleStack;
            vector<int> neighborStack;
            particleStack.push_back(i);
            neighborStack.push_back(0);
            int molecule = numMolecules++;
            
            // Recursively tag all the bonded particles.
            
            while (particleStack.size() > 0) {
                int particle = particleStack.back();
                particleMolecule[particle] = molecule;
                int& neighbor = neighborStack.back();
                while (neighbor < particleBonds[particle].size() && particleMolecule[particleBonds[particle][neighbor]] != -1)
                    neighbor++;
                if (neighbor < particleBonds[particle].size()) {
                    particleStack.push_back(particleBonds[particle][neighbor]);
                    neighborStack.push_back(0);
                }
                else {
                    particleStack.pop_back();
                    neighborStack.pop_back();
                }
            }
        }
    
    // Build the final output vector.
    
    vector<vector<int> > molecules(numMolecules);
    for (int i = 0; i < numParticles; i++)
        molecules[particleMolecule[i]].push_back(i);
    return molecules;
}

static void writeString(ostream& stream, string str) {
    int length = str.size();
    stream.write((char*) &length, sizeof(int));
    stream.write((char*) &str[0], length);
}

static string readString(istream& stream) {
    int length;
    stream.read((char*) &length, sizeof(int));
    string str(length, ' ');
    stream.read((char*) &str[0], length);
    return str;
}

void ContextImpl::createCheckpoint(ostream& stream) {
    stream.write(CHECKPOINT_MAGIC_BYTES, sizeof(CHECKPOINT_MAGIC_BYTES)/sizeof(CHECKPOINT_MAGIC_BYTES[0]));
    writeString(stream, getPlatform().getName());
    int numParticles = getSystem().getNumParticles();
    stream.write((char*) &numParticles, sizeof(int));
    int numParameters = parameters.size();
    stream.write((char*) &numParameters, sizeof(int));
    for (auto& param : parameters) {
        writeString(stream, param.first);
        stream.write((char*) &param.second, sizeof(double));
    }
    updateStateDataKernel.getAs<UpdateStateDataKernel>().createCheckpoint(*this, stream);
    stream.flush();
}

void ContextImpl::loadCheckpoint(istream& stream) {
    static const int magiclength = sizeof(CHECKPOINT_MAGIC_BYTES)/sizeof(CHECKPOINT_MAGIC_BYTES[0]);
    char magicbytes[magiclength];
    stream.read(magicbytes, magiclength);
    if (memcmp(magicbytes, CHECKPOINT_MAGIC_BYTES, magiclength) != 0)
        throw OpenMMException("loadCheckpoint: Checkpoint header was not correct");

    string platformName = readString(stream);
    if (platformName != getPlatform().getName())
        throw OpenMMException("loadCheckpoint: Checkpoint was created with a different Platform: "+platformName);
    int numParticles;
    stream.read((char*) &numParticles, sizeof(int));
    if (numParticles != getSystem().getNumParticles())
        throw OpenMMException("loadCheckpoint: Checkpoint contains the wrong number of particles");
    int numParameters;
    stream.read((char*) &numParameters, sizeof(int));
    for (int i = 0; i < numParameters; i++) {
        string name = readString(stream);
        double value;
        stream.read((char*) &value, sizeof(double));
        parameters[name] = value;
    }
    updateStateDataKernel.getAs<UpdateStateDataKernel>().loadCheckpoint(*this, stream);
    hasSetPositions = true;
}

void ContextImpl::systemChanged() {
    integrator.stateChanged(State::Energy);
}

Context* ContextImpl::createLinkedContext(const System& system, Integrator& integrator) {
    return new Context(system, integrator, *this);
}
