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

#include "ReferenceRpmdKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"

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

ReferenceIntegrateRPMDStepKernel::~ReferenceIntegrateRPMDStepKernel() {
    if (fft != NULL)
        fftpack_destroy(fft);
}
void ReferenceIntegrateRPMDStepKernel::initialize(const System& system, const RPMDIntegrator& integrator) {
    int numCopies = integrator.getNumCopies();
    int numParticles = system.getNumParticles();
    positions.resize(numCopies);
    velocities.resize(numCopies);
    forces.resize(numCopies);
    for (int i = 0; i < numCopies; i++) {
        positions[i].resize(numParticles);
        velocities[i].resize(numParticles);
        forces[i].resize(numParticles);
    }
    fftpack_init_1d(&fft, numCopies);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    const int numCopies = positions.size();
    const int numParticles = positions[0].size();
    const RealOpenMM dt = integrator.getStepSize();
    const RealOpenMM halfdt = 0.5*dt;
    const System& system = context.getSystem();
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& vel = extractVelocities(context);
    vector<RealVec>& f = extractForces(context);
    
    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid) {
        for (int i = 0; i < numCopies; i++) {
            pos = positions[i];
            vel = velocities[i];
            context.computeVirtualSites();
            context.calcForcesAndEnergy(true, false);
            forces[i] = f;
        }
    }

    // Apply the PILE-L thermostat.
    
    vector<t_complex> v(numCopies);
    vector<t_complex> q(numCopies);
    const RealOpenMM hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    const RealOpenMM scale = 1.0/sqrt((RealOpenMM) numCopies);
    const RealOpenMM nkT = numCopies*BOLTZ*integrator.getTemperature();
    const RealOpenMM twown = 2.0*nkT/hbar;
    const RealOpenMM c1_0 = exp(-halfdt*integrator.getFriction());
    const RealOpenMM c2_0 = sqrt(1.0-c1_0*c1_0);
    for (int particle = 0; particle < numParticles; particle++) {
        if (system.getParticleMass(particle) == 0.0)
            continue;
        const RealOpenMM c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
        for (int component = 0; component < 3; component++) {
            for (int k = 0; k < numCopies; k++)
                v[k] = t_complex(scale*velocities[k][particle][component], 0.0);
            fftpack_exec_1d(fft, FFTPACK_FORWARD, &v[0], &v[0]);
            
            // Apply a local Langevin thermostat to the centroid mode.

            v[0].re = v[0].re*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

            // Use critical damping white noise for the remaining modes.
            
            for (int k = 1; k <= numCopies/2; k++) {
                const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                const RealOpenMM wk = twown*sin(k*M_PI/numCopies);
                const RealOpenMM c1 = exp(-2.0*wk*halfdt);
                const RealOpenMM c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                const RealOpenMM c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                RealOpenMM rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                RealOpenMM rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                v[k] = v[k]*c1 + t_complex(rand1, rand2);
                if (k < numCopies-k)
                    v[numCopies-k] = v[numCopies-k]*c1 + t_complex(rand1, -rand2);
            }
            fftpack_exec_1d(fft, FFTPACK_BACKWARD, &v[0], &v[0]);
            for (int k = 0; k < numCopies; k++)
                velocities[k][particle][component] = scale*v[k].re;
        }
    }

    // Update velocities.
    
    for (int i = 0; i < numCopies; i++)
        for (int j = 0; j < numParticles; j++)
            if (system.getParticleMass(j) != 0.0)
                velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    for (int particle = 0; particle < numParticles; particle++) {
        if (system.getParticleMass(particle) == 0.0)
            continue;
        for (int component = 0; component < 3; component++) {
            for (int k = 0; k < numCopies; k++) {
                q[k] = t_complex(scale*positions[k][particle][component], 0.0);
                v[k] = t_complex(scale*velocities[k][particle][component], 0.0);
            }
            fftpack_exec_1d(fft, FFTPACK_FORWARD, &q[0], &q[0]);
            fftpack_exec_1d(fft, FFTPACK_FORWARD, &v[0], &v[0]);
            q[0] += v[0]*dt;
            for (int k = 1; k < numCopies; k++) {
                const RealOpenMM wk = twown*sin(k*M_PI/numCopies);
                const RealOpenMM wt = wk*dt;
                const RealOpenMM coswt = cos(wt);
                const RealOpenMM sinwt = sin(wt);
                const RealOpenMM wm = wk*system.getParticleMass(particle);
                const t_complex vprime = v[k]*coswt - q[k]*(wk*sinwt); // Advance velocity from t to t+dt
                q[k] = v[k]*(sinwt/wk) + q[k]*coswt; // Advance position from t to t+dt
                v[k] = vprime;
            }
            fftpack_exec_1d(fft, FFTPACK_BACKWARD, &q[0], &q[0]);
            fftpack_exec_1d(fft, FFTPACK_BACKWARD, &v[0], &v[0]);
            for (int k = 0; k < numCopies; k++) {
                positions[k][particle][component] = scale*q[k].re;
                velocities[k][particle][component] = scale*v[k].re;
            }
        }
    }
    
    // Calculate forces based on the updated positions.
    
    for (int i = 0; i < numCopies; i++) {
        pos = positions[i];
        vel = velocities[i];
        context.computeVirtualSites();
        context.updateContextState();
        positions[i] = pos;
        velocities[i] = vel;
        context.calcForcesAndEnergy(true, false);
        forces[i] = f;
    }

    // Update velocities.
    
    for (int i = 0; i < numCopies; i++)
        for (int j = 0; j < numParticles; j++)
            if (system.getParticleMass(j) != 0.0)
                velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));

    // Apply the PILE-L thermostat again.
    
    for (int particle = 0; particle < numParticles; particle++) {
        if (system.getParticleMass(particle) == 0.0)
            continue;
        const RealOpenMM c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
        for (int component = 0; component < 3; component++) {
            for (int k = 0; k < numCopies; k++)
                v[k] = t_complex(scale*velocities[k][particle][component], 0.0);
            fftpack_exec_1d(fft, FFTPACK_FORWARD, &v[0], &v[0]);
            
            // Apply a local Langevin thermostat to the centroid mode.

            v[0].re = v[0].re*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

            // Use critical damping white noise for the remaining modes.
            
            for (int k = 1; k <= numCopies/2; k++) {
                const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                const RealOpenMM wk = twown*sin(k*M_PI/numCopies);
                const RealOpenMM c1 = exp(-2.0*wk*halfdt);
                const RealOpenMM c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                const RealOpenMM c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                RealOpenMM rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                RealOpenMM rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                v[k] = v[k]*c1 + t_complex(rand1, rand2);
                if (k < numCopies-k)
                    v[numCopies-k] = v[numCopies-k]*c1 + t_complex(rand1, -rand2);
            }
            fftpack_exec_1d(fft, FFTPACK_BACKWARD, &v[0], &v[0]);
            for (int k = 0; k < numCopies; k++)
                velocities[k][particle][component] = scale*v[k].re;
        }
    }
    
    // Update the time.
    
    context.setTime(context.getTime()+dt);
}

double ReferenceIntegrateRPMDStepKernel::computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    vector<RealVec>& velData = extractVelocities(context);
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double mass = system.getParticleMass(i);
        if (mass > 0) {
            RealVec v = velData[i];
            energy += mass*(v.dot(v));
        }
    }
    return 0.5*energy;
}

void ReferenceIntegrateRPMDStepKernel::setPositions(int copy, const vector<Vec3>& pos) {
    int numParticles = positions[copy].size();
    for (int i = 0; i < numParticles; i++)
        positions[copy][i] = pos[i];
}

void ReferenceIntegrateRPMDStepKernel::setVelocities(int copy, const vector<Vec3>& vel) {
    int numParticles = velocities[copy].size();
    for (int i = 0; i < numParticles; i++)
        velocities[copy][i] = vel[i];
}

void ReferenceIntegrateRPMDStepKernel::copyToContext(int copy, ContextImpl& context) {
    extractPositions(context) = positions[copy];
    extractVelocities(context) = velocities[copy];
}
