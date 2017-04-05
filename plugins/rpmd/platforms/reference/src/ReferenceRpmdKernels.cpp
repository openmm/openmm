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
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKOpenMMUtilities.h"

using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->velocities);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

ReferenceIntegrateRPMDStepKernel::~ReferenceIntegrateRPMDStepKernel() {
    if (fft != NULL)
        fftpack_destroy(fft);
    for (auto& c : contractionFFT)
        if (c.second != NULL)
            fftpack_destroy(c.second);
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
    
    // Build a list of contractions.
    
    groupsNotContracted = -1;
    const map<int, int>& contractions = integrator.getContractions();
    int maxContractedCopies = 0;
    for (auto& c : contractions) {
        int group = c.first;
        int copies = c.second;
        if (group < 0 || group > 31)
            throw OpenMMException("RPMDIntegrator: Force group must be between 0 and 31");
        if (copies < 0 || copies > numCopies)
            throw OpenMMException("RPMDIntegrator: Number of copies for contraction cannot be greater than the total number of copies being simulated");
        if (copies != numCopies) {
            if (groupsByCopies.find(copies) == groupsByCopies.end()) {
                groupsByCopies[copies] = 1<<group;
                contractionFFT[copies] = NULL;
                fftpack_init_1d(&contractionFFT[copies], copies);
                if (copies > maxContractedCopies)
                    maxContractedCopies = copies;
            }
            else
                groupsByCopies[copies] |= 1<<group;
            groupsNotContracted -= 1<<group;
        }
    }
    
    // Create workspace for doing contractions.
    
    contractedPositions.resize(maxContractedCopies);
    contractedForces.resize(maxContractedCopies);
    for (int i = 0; i < maxContractedCopies; i++) {
        contractedPositions[i].resize(numParticles);
        contractedForces[i].resize(numParticles);
    }
}

void ReferenceIntegrateRPMDStepKernel::execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid) {
    const int numCopies = positions.size();
    const int numParticles = positions[0].size();
    const double dt = integrator.getStepSize();
    const double halfdt = 0.5*dt;
    const System& system = context.getSystem();
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& vel = extractVelocities(context);
    vector<Vec3>& f = extractForces(context);
    
    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid)
        computeForces(context, integrator);

    // Apply the PILE-L thermostat.
    
    vector<t_complex> v(numCopies);
    vector<t_complex> q(numCopies);
    const double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    const double scale = 1.0/sqrt((double) numCopies);
    const double nkT = numCopies*BOLTZ*integrator.getTemperature();
    const double twown = 2.0*nkT/hbar;
    const double c1_0 = exp(-halfdt*integrator.getFriction());
    const double c2_0 = sqrt(1.0-c1_0*c1_0);
    if (integrator.getApplyThermostat()) {
        for (int particle = 0; particle < numParticles; particle++) {
            if (system.getParticleMass(particle) == 0.0)
                continue;
            const double c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
            for (int component = 0; component < 3; component++) {
                for (int k = 0; k < numCopies; k++)
                    v[k] = t_complex(scale*velocities[k][particle][component], 0.0);
                fftpack_exec_1d(fft, FFTPACK_FORWARD, &v[0], &v[0]);

                // Apply a local Langevin thermostat to the centroid mode.

                v[0].re = v[0].re*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

                // Use critical damping white noise for the remaining modes.

                for (int k = 1; k <= numCopies/2; k++) {
                    const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                    const double wk = twown*sin(k*M_PI/numCopies);
                    const double c1 = exp(-2.0*wk*halfdt);
                    const double c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                    const double c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                    double rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                    double rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                    v[k] = v[k]*c1 + t_complex(rand1, rand2);
                    if (k < numCopies-k)
                        v[numCopies-k] = v[numCopies-k]*c1 + t_complex(rand1, -rand2);
                }
                fftpack_exec_1d(fft, FFTPACK_BACKWARD, &v[0], &v[0]);
                for (int k = 0; k < numCopies; k++)
                    velocities[k][particle][component] = scale*v[k].re;
            }
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
                const double wk = twown*sin(k*M_PI/numCopies);
                const double wt = wk*dt;
                const double coswt = cos(wt);
                const double sinwt = sin(wt);
                const double wm = wk*system.getParticleMass(particle);
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
    
    computeForces(context, integrator);

    // Update velocities.
    
    for (int i = 0; i < numCopies; i++)
        for (int j = 0; j < numParticles; j++)
            if (system.getParticleMass(j) != 0.0)
                velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));

    // Apply the PILE-L thermostat again.
    
    if (integrator.getApplyThermostat()) {
        for (int particle = 0; particle < numParticles; particle++) {
            if (system.getParticleMass(particle) == 0.0)
                continue;
            const double c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
            for (int component = 0; component < 3; component++) {
                for (int k = 0; k < numCopies; k++)
                    v[k] = t_complex(scale*velocities[k][particle][component], 0.0);
                fftpack_exec_1d(fft, FFTPACK_FORWARD, &v[0], &v[0]);

                // Apply a local Langevin thermostat to the centroid mode.

                v[0].re = v[0].re*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

                // Use critical damping white noise for the remaining modes.

                for (int k = 1; k <= numCopies/2; k++) {
                    const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                    const double wk = twown*sin(k*M_PI/numCopies);
                    const double c1 = exp(-2.0*wk*halfdt);
                    const double c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                    const double c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                    double rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                    double rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                    v[k] = v[k]*c1 + t_complex(rand1, rand2);
                    if (k < numCopies-k)
                        v[numCopies-k] = v[numCopies-k]*c1 + t_complex(rand1, -rand2);
                }
                fftpack_exec_1d(fft, FFTPACK_BACKWARD, &v[0], &v[0]);
                for (int k = 0; k < numCopies; k++)
                    velocities[k][particle][component] = scale*v[k].re;
            }
        }
    }
    
    // Update the time.
    
    context.setTime(context.getTime()+dt);
}

void ReferenceIntegrateRPMDStepKernel::computeForces(ContextImpl& context, const RPMDIntegrator& integrator) {
    const int totalCopies = positions.size();
    const int numParticles = positions[0].size();
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& vel = extractVelocities(context);
    vector<Vec3>& f = extractForces(context);
    
    // Compute forces from all groups that didn't have a specified contraction.
    
    for (int i = 0; i < totalCopies; i++) {
        pos = positions[i];
        vel = velocities[i];
        context.computeVirtualSites();
        Vec3 initialBox[3];
        context.getPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
        context.updateContextState();
        Vec3 finalBox[3];
        context.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
        if (initialBox[0] != finalBox[0] || initialBox[1] != finalBox[1] || initialBox[2] != finalBox[2])
            throw OpenMMException("Standard barostats cannot be used with RPMDIntegrator.  Use RPMDMonteCarloBarostat instead.");
        positions[i] = pos;
        velocities[i] = vel;
        context.calcForcesAndEnergy(true, false, groupsNotContracted);
        forces[i] = f;
    }
    
    // Now loop over contractions and compute forces from them.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        int groupFlags = g.second;
        fftpack* shortFFT = contractionFFT[copies];
        
        // Find the contracted positions.
        
        vector<t_complex> q(totalCopies);
        const double scale1 = 1.0/totalCopies;
        for (int particle = 0; particle < numParticles; particle++) {
            for (int component = 0; component < 3; component++) {
                // Transform to the frequency domain, set high frequency components to zero, and transform back.
                
                for (int k = 0; k < totalCopies; k++)
                    q[k] = t_complex(positions[k][particle][component], 0.0);
                fftpack_exec_1d(fft, FFTPACK_FORWARD, &q[0], &q[0]);
                if (copies > 1) {
                    int start = (copies+1)/2;
                    int end = totalCopies-copies+start;
                    for (int k = end; k < totalCopies; k++)
                        q[k-(totalCopies-copies)] = q[k];
                    fftpack_exec_1d(shortFFT, FFTPACK_BACKWARD, &q[0], &q[0]);
                }
                for (int k = 0; k < copies; k++)
                    contractedPositions[k][particle][component] = scale1*q[k].re;
            }
        }
        
        // Compute forces.

        for (int i = 0; i < copies; i++) {
            pos = contractedPositions[i];
            context.computeVirtualSites();
            context.calcForcesAndEnergy(true, false, groupFlags);
            contractedForces[i] = f;
        }
        
        // Apply the forces to the original copies.
        
        const double scale2 = 1.0/copies;
        for (int particle = 0; particle < numParticles; particle++) {
            for (int component = 0; component < 3; component++) {
                // Transform to the frequency domain, pad with zeros, and transform back.
                
                for (int k = 0; k < copies; k++)
                    q[k] = t_complex(contractedForces[k][particle][component], 0.0);
                if (copies > 1)
                    fftpack_exec_1d(shortFFT, FFTPACK_FORWARD, &q[0], &q[0]);
                int start = (copies+1)/2;
                int end = totalCopies-copies+start;
                for (int k = end; k < totalCopies; k++)
                    q[k] = q[k-(totalCopies-copies)];
                for (int k = start; k < end; k++)
                    q[k] = t_complex(0, 0);
                fftpack_exec_1d(fft, FFTPACK_BACKWARD, &q[0], &q[0]);
                for (int k = 0; k < totalCopies; k++)
                    forces[k][particle][component] += scale2*q[k].re;
            }
        }
    }
}

double ReferenceIntegrateRPMDStepKernel::computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    vector<Vec3>& velData = extractVelocities(context);
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double mass = system.getParticleMass(i);
        if (mass > 0) {
            Vec3 v = velData[i];
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
