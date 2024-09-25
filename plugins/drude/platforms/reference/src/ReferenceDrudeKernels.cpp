/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2024 Stanford University and the Authors.      *
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
#include "SimTKOpenMMUtilities.h"
#include "ReferenceConstraints.h"
#include "ReferenceVirtualSites.h"
#include "ReferenceForce.h"
#include <set>

using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->positions;
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->velocities;
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->forces;
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return data->periodicBoxVectors;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->constraints;
}

static const ReferenceVirtualSites& extractVirtualSites(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->virtualSites;
}

static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& inverseMasses, double timeShift) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    vector<Vec3> shiftedVel(numParticles);
    context.computeShiftedVelocities(timeShift, shiftedVel);
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
    periodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcDrudeForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);
    int numParticles = particle.size();
    double energy = 0;
    
    // Compute the interactions from the harmonic springs.
    
    for (int i = 0; i < numParticles; i++) {
        int p = particle[i];
        int p1 = particle1[i];
        int p2 = particle2[i];
        int p3 = particle3[i];
        int p4 = particle4[i];
        
        double a1 = (p2 == -1 ? 1 : aniso12[i]);
        double a2 = (p3 == -1 || p4 == -1 ? 1 : aniso34[i]);
        double a3 = 3-a1-a2;
        double k3 = ONE_4PI_EPS0*charge[i]*charge[i]/(polarizability[i]*a3);
        double k1 = ONE_4PI_EPS0*charge[i]*charge[i]/(polarizability[i]*a1) - k3;
        double k2 = ONE_4PI_EPS0*charge[i]*charge[i]/(polarizability[i]*a2) - k3;
        
        // Compute the isotropic force.
        
        Vec3 delta = pos[p]-pos[p1];
        double r2 = delta.dot(delta);
        energy += 0.5*k3*r2;
        force[p] -= delta*k3;
        force[p1] += delta*k3;
        
        // Compute the first anisotropic force.
        
        if (p2 != -1) {
            Vec3 dir = pos[p1]-pos[p2];
            double invDist = 1.0/sqrt(dir.dot(dir));
            dir *= invDist;
            double rprime = dir.dot(delta);
            energy += 0.5*k1*rprime*rprime;
            Vec3 f1 = dir*(k1*rprime); 
            Vec3 f2 = (delta-dir*rprime)*(k1*rprime*invDist);
            force[p] -= f1;
            force[p1] += f1-f2;
            force[p2] += f2;
        }
        
        // Compute the second anisotropic force.
        
        if (p3 != -1 && p4 != -1) {
            Vec3 dir = pos[p3]-pos[p4];
            double invDist = 1.0/sqrt(dir.dot(dir));
            dir *= invDist;
            double rprime = dir.dot(delta);
            energy += 0.5*k2*rprime*rprime;
            Vec3 f1 = dir*(k2*rprime);
            Vec3 f2 = (delta-dir*rprime)*(k2*rprime*invDist);
            force[p] -= f1;
            force[p1] += f1;
            force[p3] -= f2;
            force[p4] += f2;
        }
    }
    
    // Compute the screened interaction between bonded dipoles.
    
    int numPairs = pair1.size();
    Vec3* boxVectors = extractBoxVectors(context);
    for (int i = 0; i < numPairs; i++) {
        int dipole1 = pair1[i];
        int dipole2 = pair2[i];
        int dipole1Particles[] = {particle[dipole1], particle1[dipole1]};
        int dipole2Particles[] = {particle[dipole2], particle1[dipole2]};
        double uscale = pairThole[i]/pow(polarizability[dipole1]*polarizability[dipole2], 1.0/6.0);
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                int p1 = dipole1Particles[j];
                int p2 = dipole2Particles[k];
                double chargeProduct = charge[dipole1]*charge[dipole2]*(j == k ? 1 : -1);
                double deltaR[ReferenceForce::LastDeltaRIndex];
                if (periodic)
                    ReferenceForce::getDeltaRPeriodic(pos[p2], pos[p1], boxVectors, deltaR);
                else
                    ReferenceForce::getDeltaR(pos[p2], pos[p1], deltaR);
                Vec3 delta(deltaR[ReferenceForce::XIndex], deltaR[ReferenceForce::YIndex], deltaR[ReferenceForce::ZIndex]);
                double r = deltaR[ReferenceForce::RIndex];
                double u = r*uscale;
                double screening = 1.0 - (1.0+0.5*u)*exp(-u);
                energy += ONE_4PI_EPS0*chargeProduct*screening/r;
                Vec3 f = delta*(ONE_4PI_EPS0*chargeProduct/(r*r))*(screening/r-0.5*(1+u)*exp(-u)*uscale);
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
}

void ReferenceIntegrateDrudeLangevinStepKernel::initialize(const System& system, const DrudeLangevinIntegrator& integrator, const DrudeForce& force) {
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    // Identify particle pairs and ordinary particles.
    
    set<int> particles;
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
}

void ReferenceIntegrateDrudeLangevinStepKernel::execute(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& vel = extractVelocities(context);
    vector<Vec3>& force = extractForces(context);
    
    // Update velocities of ordinary particles.
    
    const double vscale = exp(-integrator.getStepSize()*integrator.getFriction());
    const double fscale = (1-vscale)/integrator.getFriction();
    const double kT = BOLTZ*integrator.getTemperature();
    const double noisescale = sqrt(2*kT*integrator.getFriction())*sqrt(0.5*(1-vscale*vscale)/integrator.getFriction());
    for (int index : normalParticles) {
        double invMass = particleInvMass[index];
        if (invMass != 0.0) {
            double sqrtInvMass = sqrt(invMass);
            for (int j = 0; j < 3; j++)
                vel[index][j] = vscale*vel[index][j] + fscale*invMass*force[index][j] + noisescale*sqrtInvMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        }
    }
    
    // Update velocities of Drude particle pairs.
    
    const double vscaleDrude = exp(-integrator.getStepSize()*integrator.getDrudeFriction());
    const double fscaleDrude = (1-vscaleDrude)/integrator.getDrudeFriction();
    const double kTDrude = BOLTZ*integrator.getDrudeTemperature();
    const double noisescaleDrude = sqrt(2*kTDrude*integrator.getDrudeFriction())*sqrt(0.5*(1-vscaleDrude*vscaleDrude)/integrator.getDrudeFriction());
    for (int i = 0; i < (int) pairParticles.size(); i++) {
        int p1 = pairParticles[i].first;
        int p2 = pairParticles[i].second;
        double mass1fract = pairInvTotalMass[i]/particleInvMass[p1];
        double mass2fract = pairInvTotalMass[i]/particleInvMass[p2];
        double sqrtInvTotalMass = sqrt(pairInvTotalMass[i]);
        double sqrtInvReducedMass = sqrt(pairInvReducedMass[i]);
        Vec3 cmVel = vel[p1]*mass1fract+vel[p2]*mass2fract;
        Vec3 relVel = vel[p2]-vel[p1];
        Vec3 cmForce = force[p1]+force[p2];
        Vec3 relForce = force[p2]*mass1fract - force[p1]*mass2fract;
        for (int j = 0; j < 3; j++) {
            cmVel[j] = vscale*cmVel[j] + fscale*pairInvTotalMass[i]*cmForce[j] + noisescale*sqrtInvTotalMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
            relVel[j] = vscaleDrude*relVel[j] + fscaleDrude*pairInvReducedMass[i]*relForce[j] + noisescaleDrude*sqrtInvReducedMass*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        }
        vel[p1] = cmVel-relVel*mass2fract;
        vel[p2] = cmVel+relVel*mass1fract;
    }

    // Update the particle positions.
    
    int numParticles = particleInvMass.size();
    vector<Vec3> xPrime(numParticles);
    double dt = integrator.getStepSize();
    for (int i = 0; i < numParticles; i++)
        if (particleInvMass[i] != 0.0)
            xPrime[i] = pos[i]+vel[i]*dt;
    
    // Apply constraints.
    
    extractConstraints(context).apply(pos, xPrime, particleInvMass, integrator.getConstraintTolerance());
    
    // Record the constrained positions and velocities.
    
    double dtInv = 1.0/dt;
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] = (xPrime[i]-pos[i])*dtInv;
            pos[i] = xPrime[i];
        }
    }

    // Apply hard wall constraints.

    const double maxDrudeDistance = integrator.getMaxDrudeDistance();
    if (maxDrudeDistance > 0) {
        const double hardwallscaleDrude = sqrt(kTDrude);
        for (int i = 0; i < (int) pairParticles.size(); i++) {
            int p1 = pairParticles[i].first;
            int p2 = pairParticles[i].second;
            Vec3 delta = pos[p1]-pos[p2];
            double r = sqrt(delta.dot(delta));
            double rInv = 1/r;
            if (rInv*maxDrudeDistance < 1.0) {
                // The constraint has been violated, so make the inter-particle distance "bounce"
                // off the hard wall.
                
                if (rInv*maxDrudeDistance < 0.5)
                    throw OpenMMException("Drude particle moved too far beyond hard wall constraint");
                Vec3 bondDir = delta*rInv;
                Vec3 vel1 = vel[p1];
                Vec3 vel2 = vel[p2];
                double mass1 = particleMass[p1];
                double mass2 = particleMass[p2];
                double deltaR = r-maxDrudeDistance;
                double deltaT = dt;
                double dotvr1 = vel1.dot(bondDir);
                Vec3 vb1 = bondDir*dotvr1;
                Vec3 vp1 = vel1-vb1;
                if (mass2 == 0) {
                    // The parent particle is massless, so move only the Drude particle.

                    if (dotvr1 != 0.0)
                        deltaT = deltaR/abs(dotvr1);
                    if (deltaT > dt)
                        deltaT = dt;
                    dotvr1 = -dotvr1*hardwallscaleDrude/(abs(dotvr1)*sqrt(mass1));
                    double dr = -deltaR + deltaT*dotvr1;
                    pos[p1] += bondDir*dr;
                    vel[p1] = vp1 + bondDir*dotvr1;
                }
                else {
                    // Move both particles.

                    double invTotalMass = pairInvTotalMass[i];
                    double dotvr2 = vel2.dot(bondDir);
                    Vec3 vb2 = bondDir*dotvr2;
                    Vec3 vp2 = vel2-vb2;
                    double vbCMass = (mass1*dotvr1 + mass2*dotvr2)*invTotalMass;
                    dotvr1 -= vbCMass;
                    dotvr2 -= vbCMass;
                    if (dotvr1 != dotvr2)
                        deltaT = deltaR/abs(dotvr1-dotvr2);
                    if (deltaT > dt)
                        deltaT = dt;
                    double vBond = hardwallscaleDrude/sqrt(mass1);
                    dotvr1 = -dotvr1*vBond*mass2*invTotalMass/abs(dotvr1);
                    dotvr2 = -dotvr2*vBond*mass1*invTotalMass/abs(dotvr2);
                    double dr1 = -deltaR*mass2*invTotalMass + deltaT*dotvr1;
                    double dr2 = deltaR*mass1*invTotalMass + deltaT*dotvr2;
                    dotvr1 += vbCMass;
                    dotvr2 += vbCMass;
                    pos[p1] += bondDir*dr1;
                    pos[p2] += bondDir*dr2;
                    vel[p1] = vp1 + bondDir*dotvr1;
                    vel[p2] = vp2 + bondDir*dotvr2;
                }
            }
        }
    }
    extractVirtualSites(context).computePositions(context.getSystem(), pos);
    data.time += integrator.getStepSize();
    data.stepCount++;
}

double ReferenceIntegrateDrudeLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeLangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, particleInvMass, 0.5*integrator.getStepSize());
}

ReferenceIntegrateDrudeSCFStepKernel::~ReferenceIntegrateDrudeSCFStepKernel() {
}

void ReferenceIntegrateDrudeSCFStepKernel::initialize(const System& system, const DrudeSCFIntegrator& integrator, const DrudeForce& force) {
    // Record information about the drude particles.

    int numDrude = force.getNumParticles();
    particle.resize(numDrude);
    particle1.resize(numDrude);
    particle2.resize(numDrude);
    particle3.resize(numDrude);
    particle4.resize(numDrude);
    k1.resize(numDrude);
    k2.resize(numDrude);
    k3.resize(numDrude);
    for (int i = 0; i < numDrude; i++) {
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, particle[i], particle1[i], particle2[i], particle3[i], particle4[i], charge, polarizability, aniso12, aniso34);
        double a1 = (particle2[i] == -1 ? 1 : aniso12);
        double a2 = (particle3[i] == -1 || particle4[i] == -1 ? 1 : aniso34);
        double a3 = 3-a1-a2;
        k3[i] = ONE_4PI_EPS0*charge*charge/(polarizability*a3);
        k1[i] = ONE_4PI_EPS0*charge*charge/(polarizability*a1) - k3[i];
        k2[i] = ONE_4PI_EPS0*charge*charge/(polarizability*a2) - k3[i];
    }

    // Record particle masses.

    vector<double> particleMass;
    for (int i = 0; i < system.getNumParticles(); i++) {
        double mass = system.getParticleMass(i);
        particleMass.push_back(mass);
        particleInvMass.push_back(mass == 0.0 ? 0.0 : 1.0/mass);
    }
}

void ReferenceIntegrateDrudeSCFStepKernel::execute(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& vel = extractVelocities(context);
    vector<Vec3>& force = extractForces(context);
    
    // Update the positions and velocities.
    
    int numParticles = particleInvMass.size();
    vector<Vec3> xPrime(numParticles);
    double dt = integrator.getStepSize();
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] += force[i]*particleInvMass[i]*dt;
            xPrime[i] = pos[i]+vel[i]*dt;
        }
    }
        
    // Apply constraints.
    
    extractConstraints(context).apply(pos, xPrime, particleInvMass, integrator.getConstraintTolerance());
    
    // Record the constrained positions and velocities.
    
    double dtInv = 1.0/dt;
    for (int i = 0; i < numParticles; i++) {
        if (particleInvMass[i] != 0.0) {
            vel[i] = (xPrime[i]-pos[i])*dtInv;
            pos[i] = xPrime[i];
        }
    }
    
    // Update the positions of virtual sites and Drude particles.
    
    extractVirtualSites(context).computePositions(context.getSystem(), pos);
    minimize(context, integrator.getMinimizationErrorTolerance());
    data.time += integrator.getStepSize();
    data.stepCount++;
}

double ReferenceIntegrateDrudeSCFStepKernel::computeKineticEnergy(ContextImpl& context, const DrudeSCFIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, particleInvMass, 0.5*integrator.getStepSize());
}

void ReferenceIntegrateDrudeSCFStepKernel::minimize(ContextImpl& context, double tolerance) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);
    int numDrude = particle.size();
    double lastForce = 0;
    for (int iteration = 0; iteration < 50; iteration++) {
        context.calcForcesAndEnergy(true, false, context.getIntegrator().getIntegrationForceGroups());
        double totalForce = 0;
        for (int i = 0; i < numDrude; i++) {
            int p = particle[i];
            int p1 = particle1[i];
            int p2 = particle2[i];
            int p3 = particle3[i];
            int p4 = particle4[i];
            Vec3 fscale(k3[i], k3[i], k3[i]);
            if (p2 != -1) {
                Vec3 dir = pos[p1]-pos[p2];
                dir /= sqrt(dir.dot(dir));;
                fscale += k1[i]*dir;
            }
            if (p3 != -1 && p4 != -1) {
                Vec3 dir = pos[p3]-pos[p4];
                dir /= sqrt(dir.dot(dir));;
                fscale += k2[i]*dir;
            }
            Vec3 f = force[p];
            double f2 = f.dot(f);
            totalForce += f2;
            double damping = (sqrt(f2) > 10*tolerance ? 0.5 : 1.0);
            for (int i = 0; i < 3; i++)
                pos[p][i] += damping*f[i]/fscale[i];
        }
        if (sqrt(totalForce/(3*numDrude)) < tolerance || (iteration > 0 && totalForce > 0.9*lastForce))
            break;
        lastForce = totalForce;
    }
}
