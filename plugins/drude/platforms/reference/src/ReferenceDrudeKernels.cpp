/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2013 Stanford University and the Authors.      *
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

#include "ReferenceDrudeKernels.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "SimTKReference/ReferenceCCMAAlgorithm.h"
#include "SimTKReference/ReferenceVirtualSites.h"
#include <set>

using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static void findAnglesForCCMA(const System& system, vector<ReferenceCCMAAlgorithm::AngleInfo>& angles) {
    for (int i = 0; i < system.getNumForces(); i++) {
        const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
        if (force != NULL) {
            for (int j = 0; j < force->getNumAngles(); j++) {
                int atom1, atom2, atom3;
                double angle, k;
                force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, (RealOpenMM)angle));
            }
        }
    }
}

static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& inverseMasses, double timeShift, ReferenceConstraintAlgorithm* constraints) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    
    // Compute the shifted velocities.
    
    vector<RealVec> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (inverseMasses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift*inverseMasses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.
    
    if (constraints != NULL) {
        constraints->setTolerance(1e-4);
        constraints->applyToVelocities(numParticles, posData, shiftedVel, inverseMasses);
    }
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (inverseMasses[i] > 0)
            energy += (shiftedVel[i].dot(shiftedVel[i]))/inverseMasses[i];
    return 0.5*energy;
}


void ReferenceCalcDrudeForceKernel::initialize(const System& system, const DrudeForce& force) {
    // Initialize particle parameters.
    
    int numParticles = force.getNumParticles();
    particle.resize(numParticles);
    particle1.resize(numParticles);
    particle2.resize(numParticles);
    particle3.resize(numParticles);
    particle4.resize(numParticles);
    charge.resize(numParticles);
    polarizability.resize(numParticles);
    aniso12.resize(numParticles);
    aniso34.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
        force.getParticleParameters(i, particle[i], particle1[i], particle2[i], particle3[i], particle4[i], charge[i], polarizability[i], aniso12[i], aniso34[i]);
    
    // Initialize screened pair parameters.
    
    int numPairs = force.getNumScreenedPairs();
    pair1.resize(numPairs);
    pair2.resize(numPairs);
    pairThole.resize(numPairs);
    for (int i = 0; i < numPairs; i++)
        force.getScreenedPairParameters(i, pair1[i], pair2[i], pairThole[i]);
}

double ReferenceCalcDrudeForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    int numParticles = particle.size();
    double energy = 0;
    
    // Compute the interactions from the harmonic springs.
    
    for (int i = 0; i < numParticles; i++) {
        int p = particle[i];
        int p1 = particle1[i];
        int p2 = particle2[i];
        int p3 = particle3[i];
        int p4 = particle4[i];
        
        RealOpenMM a1 = (p2 == -1 ? 1 : aniso12[i]);
        RealOpenMM a2 = (p3 == -1 || p4 == -1 ? 1 : aniso34[i]);
        RealOpenMM a3 = 3-a1-a2;
        RealOpenMM k3 = charge[i]*charge[i]/(polarizability[i]*a3);
        RealOpenMM k1 = charge[i]*charge[i]/(polarizability[i]*a1) - k3;
        RealOpenMM k2 = charge[i]*charge[i]/(polarizability[i]*a2) - k3;
        
        // Compute the isotropic force.
        
        RealVec delta = pos[p]-pos[p1];
        RealOpenMM r2 = delta.dot(delta);
        energy += 0.5*k3*r2;
        force[p] -= delta*k3;
        force[p1] += delta*k3;
        
        // Compute the first anisotropic force.
        
        if (p2 != -1) {
            RealVec dir = pos[p1]-pos[p2];
            RealOpenMM invDist = 1.0/sqrt(dir.dot(dir));
            dir *= invDist;
            RealOpenMM rprime = dir.dot(delta);
            energy += 0.5*k1*rprime*rprime;
            RealVec f1 = dir*(k1*rprime); 
            RealVec f2 = (delta-dir*rprime)*(k1*rprime*invDist);
            force[p] -= f1;
            force[p1] += f1-f2;
            force[p2] += f2;
        }
        
        // Compute the second anisotropic force.
        
        if (p3 != -1 && p4 != -1) {
            RealVec dir = pos[p3]-pos[p4];
            RealOpenMM invDist = 1.0/sqrt(dir.dot(dir));
            dir *= invDist;
            RealOpenMM rprime = dir.dot(delta);
            energy += 0.5*k2*rprime*rprime;
            RealVec f1 = dir*(k2*rprime);
            RealVec f2 = (delta-dir*rprime)*(k2*rprime*invDist);
            force[p] -= f1;
            force[p1] += f1;
            force[p3] -= f2;
            force[p4] += f2;
        }
    }
    
    // Compute the screened interaction between bonded dipoles.
    
    int numPairs = pair1.size();
    for (int i = 0; i < numPairs; i++) {
        int dipole1 = pair1[i];
        int dipole2 = pair2[i];
        int dipole1Particles[] = {particle[dipole1], particle1[dipole1]};
        int dipole2Particles[] = {particle[dipole2], particle1[dipole2]};
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                int p1 = dipole1Particles[j];
                int p2 = dipole2Particles[k];
                RealOpenMM chargeProduct = charge[dipole1]*charge[dipole2]*(j == k ? 1 : -1);
                RealVec delta = pos[p1]-pos[p2];
                RealOpenMM r = sqrt(delta.dot(delta));
                RealOpenMM u = r*pairThole[i]/pow(polarizability[dipole1]*polarizability[dipole2], 1.0/6.0);
                RealOpenMM screening = 1.0 - (1.0+0.5*u)*exp(-u);
                energy += ONE_4PI_EPS0*chargeProduct*screening/r;
                RealVec f = delta*(ONE_4PI_EPS0*chargeProduct*screening/(r*r*r));
                force[p1] += f;
                force[p2] -= f;
            }
    }
    return energy;
}

void ReferenceCalcDrudeForceKernel::copyParametersToContext(ContextImpl& context, const DrudeForce& force) {
    if (force.getNumParticles() != particle.size())
        throw OpenMMException("updateParametersInContext: The number of Drude particles has changed");
    if (force.getNumScreenedPairs() != pair1.size())
        throw OpenMMException("updateParametersInContext: The number of screened pairs has changed");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge[i], polarizability[i], aniso12[i], aniso34[i]);
        if (p != particle[i] || p1 != particle1[i] || p2 != particle2[i] || p3 != particle3[i] || p4 != particle4[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
    for (int i = 0; i < force.getNumScreenedPairs(); i++) {
        int p1, p2;
        force.getScreenedPairParameters(i, p1, p2, pairThole[i]);
        if (p1 != pair1[i] || p2 != pair2[i])
            throw OpenMMException("updateParametersInContext: A particle index for a screened pair has changed");
    }
}

ReferenceIntegrateDrudeLangevinStepKernel::~ReferenceIntegrateDrudeLangevinStepKernel() {
    if (constraints != NULL)
        delete constraints;
}

void ReferenceIntegrateDrudeLangevinStepKernel::initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) {
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    // Identify particle pairs and ordinary particles.
    
    set<int> particles;
    vector<RealOpenMM> particleMass;
    for (int i = 0; i < system.getNumParticles(); i++) {
        particles.insert(i);
        double mass = system.getParticleMass(i);
        particleMass.push_back(mass);
        particleInvMass.push_back(mass == 0.0 ? 0.0 : 1.0/mass);
    }
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        particles.erase(p);
        particles.erase(p1);
        pairParticles.push_back(make_pair(p, p1));
        double m1 = system.getParticleMass(p);
        double m2 = system.getParticleMass(p1);
        pairInvTotalMass.push_back(1.0/(m1+m2));
        pairInvReducedMass.push_back((m1+m2)/(m1*m2));
    }
    normalParticles.insert(normalParticles.begin(), particles.begin(), particles.end());
    
    // Prepare constraints.
    
    int numConstraints = system.getNumConstraints();
    if (numConstraints > 0) {
        vector<pair<int, int> > constraintIndices(numConstraints);
        vector<RealOpenMM> constraintDistances(numConstraints);
        for (int i = 0; i < numConstraints; ++i) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            constraintIndices[i].first = particle1;
            constraintIndices[i].second = particle2;
            constraintDistances[i] = static_cast<RealOpenMM>(distance);
        }
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(system, angles);
        constraints = new ReferenceCCMAAlgorithm(system.getNumParticles(), numConstraints, constraintIndices, constraintDistances, particleMass, angles, (RealOpenMM)integrator.getConstraintTolerance());
    }
}

void ReferenceIntegrateDrudeLangevinStepKernel::execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& vel = extractVelocities(context);
    vector<RealVec>& force = extractForces(context);
    
    // Update velocities of ordinary particles.
    
    const RealOpenMM vscale = exp(-integrator.getStepSize()*integrator.getFriction());
    const RealOpenMM fscale = (1-vscale)/integrator.getFriction();
    const RealOpenMM kT = BOLTZ*integrator.getTemperature();
    const RealOpenMM noisescale = sqrt(2*kT*integrator.getFriction())*sqrt(0.5*(1-vscale*vscale)/integrator.getFriction());
    for (int i = 0; i < (int) normalParticles.size(); i++) {
        int index = normalParticles[i];
        RealOpenMM invMass = particleInvMass[index];
        if (invMass != 0.0) {
            RealOpenMM sqrtInvMass = sqrt(invMass);
            for (int j = 0; j < 3; j++)
                vel[index][j] = vscale*vel[index][j] + fscale*invMass*force[index][j] + noisescale*sqrtInvMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        }
    }
    
    // Update velocities of Drude particle pairs.
    
    const RealOpenMM vscaleDrude = exp(-integrator.getStepSize()*integrator.getDrudeFriction());
    const RealOpenMM fscaleDrude = (1-vscaleDrude)/integrator.getDrudeFriction();
    const RealOpenMM kTDrude = BOLTZ*integrator.getDrudeTemperature();
    const RealOpenMM noisescaleDrude = sqrt(2*kTDrude*integrator.getDrudeFriction())*sqrt(0.5*(1-vscaleDrude*vscaleDrude)/integrator.getDrudeFriction());
    for (int i = 0; i < (int) pairParticles.size(); i++) {
        int p1 = pairParticles[i].first;
        int p2 = pairParticles[i].second;
        RealOpenMM mass1fract = pairInvTotalMass[i]/particleInvMass[p1];
        RealOpenMM mass2fract = pairInvTotalMass[i]/particleInvMass[p2];
        RealOpenMM sqrtInvTotalMass = sqrt(pairInvTotalMass[i]);
        RealOpenMM sqrtInvReducedMass = sqrt(pairInvReducedMass[i]);
        RealVec cmVel = vel[p1]*mass1fract+vel[p2]*mass2fract;
        RealVec relVel = vel[p2]-vel[p1];
        RealVec cmForce = force[p1]+force[p2];
        RealVec relForce = force[p2]*mass1fract - force[p1]*mass2fract;
        for (int j = 0; j < 3; j++) {
            cmVel[j] = vscale*cmVel[j] + fscale*pairInvTotalMass[i]*cmForce[j] + noisescale*sqrtInvTotalMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
            relVel[j] = vscaleDrude*relVel[j] + fscaleDrude*pairInvReducedMass[i]*relForce[j] + noisescaleDrude*sqrtInvReducedMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        }
        vel[p1] = cmVel-relVel*mass2fract;
        vel[p2] = cmVel+relVel*mass1fract;
    }

    // Update the particle positions.
    
    int numParticles = particleInvMass.size();
    vector<RealVec> xPrime(numParticles);
    RealOpenMM dt = integrator.getStepSize();
    for (int i = 0; i < numParticles; i++)
        if (particleInvMass[i] != 0.0)
            xPrime[i] = pos[i]+vel[i]*dt;
    
    // Apply constraints.
    
    if (constraints != NULL)
        constraints->apply(numParticles, pos, xPrime, particleInvMass);
    
    // Record the constrained positions and velocities.
    
    RealOpenMM dtInv = 1.0/dt;
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] = (xPrime[i]-pos[i])*dtInv;
            pos[i] = xPrime[i];
        }
    }
    ReferenceVirtualSites::computePositions(context.getSystem(), pos);
    data.time += integrator.getStepSize();
    data.stepCount++;
}

double ReferenceIntegrateDrudeLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, particleInvMass, 0.5*integrator.getStepSize(), constraints);
}

ReferenceIntegrateDrudeSCFStepKernel::~ReferenceIntegrateDrudeSCFStepKernel() {
    if (constraints != NULL)
        delete constraints;
    if (minimizerPos != NULL)
        lbfgs_free(minimizerPos);
}

void ReferenceIntegrateDrudeSCFStepKernel::initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) {
    // Identify Drude particles.
    
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        drudeParticles.push_back(p);
    }

    // Record particle masses.

    vector<RealOpenMM> particleMass;
    for (int i = 0; i < system.getNumParticles(); i++) {
        double mass = system.getParticleMass(i);
        particleMass.push_back(mass);
        particleInvMass.push_back(mass == 0.0 ? 0.0 : 1.0/mass);
    }
    
    // Prepare constraints.
    
    int numConstraints = system.getNumConstraints();
    if (numConstraints > 0) {
        vector<pair<int, int> > constraintIndices(numConstraints);
        vector<RealOpenMM> constraintDistances(numConstraints);
        for (int i = 0; i < numConstraints; ++i) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            constraintIndices[i].first = particle1;
            constraintIndices[i].second = particle2;
            constraintDistances[i] = static_cast<RealOpenMM>(distance);
        }
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(system, angles);
        constraints = new ReferenceCCMAAlgorithm(system.getNumParticles(), numConstraints, constraintIndices, constraintDistances, particleMass, angles, (RealOpenMM)integrator.getConstraintTolerance());
    }
    
    // Initialize the energy minimizer.
    
    minimizerPos = lbfgs_malloc(drudeParticles.size()*3);
    if (minimizerPos == NULL)
        throw OpenMMException("DrudeSCFIntegrator: Failed to allocate memory");
    lbfgs_parameter_init(&minimizerParams);
    minimizerParams.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    if (sizeof(RealOpenMM) < 8)
        minimizerParams.xtol = 1e-7;
}

void ReferenceIntegrateDrudeSCFStepKernel::execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& vel = extractVelocities(context);
    vector<RealVec>& force = extractForces(context);
    
    // Update the positions and velocities.
    
    int numParticles = particleInvMass.size();
    vector<RealVec> xPrime(numParticles);
    RealOpenMM dt = integrator.getStepSize();
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] += force[i]*particleInvMass[i]*dt;
            xPrime[i] = pos[i]+vel[i]*dt;
        }
    }
        
    // Apply constraints.
    
    if (constraints != NULL)
        constraints->apply(numParticles, pos, xPrime, particleInvMass);
    
    // Record the constrained positions and velocities.
    
    RealOpenMM dtInv = 1.0/dt;
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] = (xPrime[i]-pos[i])*dtInv;
            pos[i] = xPrime[i];
        }
    }
    
    // Update the positions of virtual sites and Drude particles.
    
    ReferenceVirtualSites::computePositions(context.getSystem(), pos);
    minimize(context, integrator.getMinimizationErrorTolerance());
    data.time += integrator.getStepSize();
    data.stepCount++;
}

double ReferenceIntegrateDrudeSCFStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, particleInvMass, 0.5*integrator.getStepSize(), constraints);
}

struct MinimizerData {
    ContextImpl& context;
    vector<int>& drudeParticles;
    MinimizerData(ContextImpl& context, vector<int>& drudeParticles) : context(context), drudeParticles(drudeParticles) {}
};

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    MinimizerData* data = reinterpret_cast<MinimizerData*>(instance);
    ContextImpl& context = data->context;
    vector<int>& drudeParticles = data->drudeParticles;
    int numDrudeParticles = drudeParticles.size();

    // Compute the force and energy for this configuration.

    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    for (int i = 0; i < numDrudeParticles; i++)
        pos[drudeParticles[i]] = RealVec(x[3*i], x[3*i+1], x[3*i+2]);
    double energy = context.calcForcesAndEnergy(true, true);
    for (int i = 0; i < numDrudeParticles; i++) {
        RealVec f = force[drudeParticles[i]];
        g[3*i] = -f[0];
        g[3*i+1] = -f[1];
        g[3*i+2] = -f[2];
    }
    return energy;
}

void ReferenceIntegrateDrudeSCFStepKernel::minimize(ContextImpl& context, double tolerance) {
    // Record the initial positions and determine a normalization constant for scaling the tolerance.

    vector<RealVec>& pos = extractPositions(context);
    int numDrudeParticles = drudeParticles.size();
    double norm = 0.0;
    for (int i = 0; i < numDrudeParticles; i++) {
        RealVec p = pos[drudeParticles[i]];
        minimizerPos[3*i] = p[0];
        minimizerPos[3*i+1] = p[1];
        minimizerPos[3*i+2] = p[2];
        norm += p.dot(p);
    }
    norm /= numDrudeParticles;
    norm = (norm < 1 ? 1 : sqrt(norm));
    minimizerParams.epsilon = tolerance/norm;
    
    // Perform the minimization.

    lbfgsfloatval_t fx;
    MinimizerData data(context, drudeParticles);
    lbfgs(numDrudeParticles*3, minimizerPos, &fx, evaluate, NULL, &data, &minimizerParams);
}