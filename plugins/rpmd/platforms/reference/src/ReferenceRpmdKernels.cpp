/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org.                                        *
 *                                                                            *
 * Portions copyright (c) 2011-2022 Stanford University and the Authors.      *
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
#include <set>
#ifdef _MSC_VER
  #define POCKETFFT_NO_VECTORS
#endif
#include "pocketfft_hdronly.h"
#include <complex>
#include <fstream>
#include <sstream>

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
                if (copies > maxContractedCopies)
                    maxContractedCopies = copies;
            }
            else
                groupsByCopies[copies] |= 1<<group;
            groupsNotContracted -= 1<<group;
        }
    }
    groupsNotContracted &= integrator.getIntegrationForceGroups();

    const map<int, int>& particleTypes = integrator.getParticleTypes();
    const set<int>& quantumTypes = integrator.getQuantumParticleTypes();
    bool defaultQuantum = integrator.getDefaultQuantum();
    isQuantumParticle.assign(numParticles, true);
    hybridMode = false;
    if (!particleTypes.empty() || !quantumTypes.empty() || !defaultQuantum) {
        for (int i = 0; i < numParticles; i++) {
            int type = 0;
            auto it = particleTypes.find(i);
            if (it != particleTypes.end())
                type = it->second;
            isQuantumParticle[i] = (type == 0) ? defaultQuantum : (quantumTypes.count(type) > 0);
            if (!isQuantumParticle[i])
                hybridMode = true;
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
    const RPMDIntegrator::ThermostatType thermostatType = integrator.getThermostatType();
    int nInner = integrator.getNumInnerSteps();
    if (nInner < 1)
        throw OpenMMException("RPMDIntegrator: numInnerSteps must be at least 1.");
    int innerMask = integrator.getInnerForceGroups() & integrator.getIntegrationForceGroups();
    int outerMask = integrator.getIntegrationForceGroups() & ~innerMask;
    bool useMTS = (nInner > 1);
    if (useMTS && (innerMask == 0 || outerMask == 0))
        throw OpenMMException("RPMDIntegrator: MTS requires non-empty inner and outer force groups (within integrationForceGroups).");

    // Loop over copies and compute the force on each one.
    
    if (!forcesAreValid) {
        if (!useMTS)
            computeForces(context, integrator);
        else
            computeForcesWithMask(context, integrator, outerMask);
    }

    // Apply the thermostat (first half).
    
    vector<complex<double>> v(numCopies);
    vector<complex<double>> q(numCopies);
    const double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    const double scale = 1.0/sqrt((double) numCopies);
    const double nkT = numCopies*BOLTZ*integrator.getTemperature();
    const double twown = 2.0*nkT/hbar;
    const double c1_0 = exp(-halfdt*integrator.getFriction());
    const double c2_0 = sqrt(1.0-c1_0*c1_0);
    
    // Bussi thermostat parameters for centroid (PILE_G mode)
    const double c1_bussi = exp(-halfdt*integrator.getCentroidFriction());
    const bool useFFL = (thermostatType == RPMDIntegrator::FastForwardLangevin);
    vector<vector<Vec3> > velSavedFFL;
    
    if (integrator.getApplyThermostat() && thermostatType != RPMDIntegrator::NoneThermo) {
        
        // For PILE_G mode, apply Bussi thermostat to centroid first
        if (thermostatType == RPMDIntegrator::PileG) {
            applyBussiCentroidThermostat(system, integrator, numCopies, numParticles, scale, nkT, c1_bussi);
        }
        if (useFFL)
            velSavedFFL = velocities;
        
        for (int particle = 0; particle < numParticles; particle++) {
            if (system.getParticleMass(particle) == 0.0)
                continue;
            if (hybridMode && !isQuantumParticle[particle])
                continue;
            const double c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
            for (int component = 0; component < 3; component++) {
                for (int k = 0; k < numCopies; k++)
                    v[k] = complex<double>(scale*velocities[k][particle][component], 0.0);
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, v.data(), v.data(), 1.0, 0);

                if (thermostatType == RPMDIntegrator::Pile || thermostatType == RPMDIntegrator::FastForwardLangevin) {
                    v[0].real(v[0].real()*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                }

                // Use critical damping white noise for the remaining (internal) modes.

                for (int k = 1; k <= numCopies/2; k++) {
                    const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                    const double wk = twown*sin(k*M_PI/numCopies);
                    const double c1 = exp(-2.0*wk*halfdt);
                    const double c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                    const double c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                    double rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                    double rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                    v[k] = v[k]*c1 + complex<double>(rand1, rand2);
                    if (k < numCopies-k)
                        v[numCopies-k] = v[numCopies-k]*c1 + complex<double>(rand1, -rand2);
                }
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, v.data(), v.data(), 1.0, 0);
                for (int k = 0; k < numCopies; k++)
                    velocities[k][particle][component] = scale*v[k].real();
            }
        }
        if (useFFL) {
            for (int ii = 0; ii < numCopies; ii++)
                for (int jj = 0; jj < numParticles; jj++) {
                    if (hybridMode && !isQuantumParticle[jj])
                        continue;
                    for (int c = 0; c < 3; c++)
                        if (velocities[ii][jj][c] * velSavedFFL[ii][jj][c] < 0.0)
                            velocities[ii][jj][c] = -velocities[ii][jj][c];
                }
        }
    }

    auto applyThermostatSecondHalf = [&]() {
        if (integrator.getApplyThermostat() && thermostatType != RPMDIntegrator::NoneThermo) {
            if (thermostatType == RPMDIntegrator::PileG) {
                applyBussiCentroidThermostat(system, integrator, numCopies, numParticles, scale, nkT, c1_bussi);
            }
            if (useFFL)
                velSavedFFL = velocities;
            for (int particle = 0; particle < numParticles; particle++) {
                if (system.getParticleMass(particle) == 0.0)
                    continue;
                if (hybridMode && !isQuantumParticle[particle])
                    continue;
                const double c3_0 = c2_0*sqrt(nkT/system.getParticleMass(particle));
                for (int component = 0; component < 3; component++) {
                    for (int k = 0; k < numCopies; k++)
                        v[k] = complex<double>(scale*velocities[k][particle][component], 0.0);
                    pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, v.data(), v.data(), 1.0, 0);
                    if (thermostatType == RPMDIntegrator::Pile || thermostatType == RPMDIntegrator::FastForwardLangevin) {
                        v[0].real(v[0].real()*c1_0 + c3_0*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                    }
                    for (int k = 1; k <= numCopies/2; k++) {
                        const bool isCenter = (numCopies%2 == 0 && k == numCopies/2);
                        const double wk = twown*sin(k*M_PI/numCopies);
                        const double c1 = exp(-2.0*wk*halfdt);
                        const double c2 = sqrt((1.0-c1*c1)/2) * (isCenter ? sqrt(2.0) : 1.0);
                        const double c3 = c2*sqrt(nkT/system.getParticleMass(particle));
                        double rand1 = c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
                        double rand2 = (isCenter ? 0.0 : c3*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
                        v[k] = v[k]*c1 + complex<double>(rand1, rand2);
                        if (k < numCopies-k)
                            v[numCopies-k] = v[numCopies-k]*c1 + complex<double>(rand1, -rand2);
                    }
                    pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, v.data(), v.data(), 1.0, 0);
                    for (int k = 0; k < numCopies; k++)
                        velocities[k][particle][component] = scale*v[k].real();
                }
            }
            if (useFFL) {
                for (int ii = 0; ii < numCopies; ii++)
                    for (int jj = 0; jj < numParticles; jj++) {
                        if (hybridMode && !isQuantumParticle[jj])
                            continue;
                        for (int c = 0; c < 3; c++)
                            if (velocities[ii][jj][c] * velSavedFFL[ii][jj][c] < 0.0)
                                velocities[ii][jj][c] = -velocities[ii][jj][c];
                    }
            }
        }
    };

    auto evolveRingPolymer = [&](double stepDt) {
        for (int particle = 0; particle < numParticles; particle++) {
            if (system.getParticleMass(particle) == 0.0)
                continue;
            if (hybridMode && !isQuantumParticle[particle]) {
                Vec3 vcm(0, 0, 0);
                for (int k = 0; k < numCopies; k++)
                    vcm += velocities[k][particle];
                vcm *= 1.0/numCopies;
                for (int k = 0; k < numCopies; k++) {
                    positions[k][particle] += vcm * stepDt;
                    velocities[k][particle] = vcm;
                }
                continue;
            }
            for (int component = 0; component < 3; component++) {
                for (int k = 0; k < numCopies; k++) {
                    q[k] = complex<double>(scale*positions[k][particle][component], 0.0);
                    v[k] = complex<double>(scale*velocities[k][particle][component], 0.0);
                }
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, q.data(), q.data(), 1.0, 0);
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, v.data(), v.data(), 1.0, 0);
                q[0] += v[0]*stepDt;
                for (int k = 1; k < numCopies; k++) {
                    const double wk = twown*sin(k*M_PI/numCopies);
                    const double wt = wk*stepDt;
                    const double coswt = cos(wt);
                    const double sinwt = sin(wt);
                    const complex<double> vprime = v[k]*coswt - q[k]*(wk*sinwt);
                    q[k] = v[k]*(sinwt/wk) + q[k]*coswt;
                    v[k] = vprime;
                }
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, q.data(), q.data(), 1.0, 0);
                pocketfft::c2c({(size_t) numCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, v.data(), v.data(), 1.0, 0);
                for (int k = 0; k < numCopies; k++) {
                    positions[k][particle][component] = scale*q[k].real();
                    velocities[k][particle][component] = scale*v[k].real();
                }
            }
        }
    };

    if (!useMTS) {
    // Update velocities.
    
    for (int i = 0; i < numCopies; i++)
        for (int j = 0; j < numParticles; j++)
            if (system.getParticleMass(j) != 0.0)
                velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));
    
    // Evolve the free ring polymer by transforming to the frequency domain.

    evolveRingPolymer(dt);
    
    // Calculate forces based on the updated positions.
    
    computeForces(context, integrator);

    // Update velocities.
    
    for (int i = 0; i < numCopies; i++)
        for (int j = 0; j < numParticles; j++)
            if (system.getParticleMass(j) != 0.0)
                velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));

    applyThermostatSecondHalf();
    } else {
        const double dtInner = dt / nInner;
        const double halfdtInner = 0.5 * dtInner;
        for (int i = 0; i < numCopies; i++)
            for (int j = 0; j < numParticles; j++)
                if (system.getParticleMass(j) != 0.0)
                    velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));
        for (int inner = 0; inner < nInner; inner++) {
            computeForcesWithMask(context, integrator, innerMask);
            for (int i = 0; i < numCopies; i++)
                for (int j = 0; j < numParticles; j++)
                    if (system.getParticleMass(j) != 0.0)
                        velocities[i][j] += forces[i][j]*(halfdtInner/system.getParticleMass(j));
            evolveRingPolymer(dtInner);
            if (hybridMode) {
                for (int particle = 0; particle < numParticles; particle++) {
                    if (isQuantumParticle[particle] || system.getParticleMass(particle) == 0.0)
                        continue;
                    Vec3 velAvg(0, 0, 0);
                    for (int k = 0; k < numCopies; k++)
                        velAvg += velocities[k][particle];
                    velAvg *= 1.0/numCopies;
                    for (int k = 0; k < numCopies; k++) {
                        positions[k][particle] = positions[0][particle];
                        velocities[k][particle] = velAvg;
                    }
                }
            }
            computeForcesWithMask(context, integrator, innerMask);
            for (int i = 0; i < numCopies; i++)
                for (int j = 0; j < numParticles; j++)
                    if (system.getParticleMass(j) != 0.0)
                        velocities[i][j] += forces[i][j]*(halfdtInner/system.getParticleMass(j));
        }
        computeForcesWithMask(context, integrator, outerMask);
        for (int i = 0; i < numCopies; i++)
            for (int j = 0; j < numParticles; j++)
                if (system.getParticleMass(j) != 0.0)
                    velocities[i][j] += forces[i][j]*(halfdt/system.getParticleMass(j));
        applyThermostatSecondHalf();
    }

    if (hybridMode) {
        for (int particle = 0; particle < numParticles; particle++) {
            if (isQuantumParticle[particle] || system.getParticleMass(particle) == 0.0)
                continue;
            Vec3 velAvg(0, 0, 0);
            for (int k = 0; k < numCopies; k++)
                velAvg += velocities[k][particle];
            velAvg *= 1.0/numCopies;
            for (int k = 0; k < numCopies; k++) {
                positions[k][particle] = positions[0][particle];
                velocities[k][particle] = velAvg;
            }
        }
    }
    
    // Update the time.
    
    context.setTime(context.getTime()+dt);
}

void ReferenceIntegrateRPMDStepKernel::applyBussiCentroidThermostat(const System& system, const RPMDIntegrator& integrator,
                                                                    int numCopies, int numParticles, double scale,
                                                                    double nkT, double c1) {
    // Mirror i-PI ThermoSVR.step() for the centroid mode used by ThermoPILE_G.
    // This leaves the critical-damped internal PILE modes untouched and applies
    // the same stochastic velocity rescaling update to the centroid momentum.
    const double kT = BOLTZ * integrator.getTemperature();
    const double kPerDof = 0.5 * kT;
    
    // Step 1: Compute current centroid kinetic energy and store centroid velocities
    // Centroid velocity = (1/numCopies) * sum of bead velocities
    double centroidKE = 0.0;
    int ndof = 0;  // Number of degrees of freedom for centroid
    vector<Vec3> centroidVel(numParticles, Vec3(0.0, 0.0, 0.0));
    
    for (int particle = 0; particle < numParticles; particle++) {
        double mass = system.getParticleMass(particle);
        if (mass == 0.0)
            continue;
        if (hybridMode && !isQuantumParticle[particle])
            continue;
        
        // Compute centroid velocity for this particle
        for (int copy = 0; copy < numCopies; copy++) {
            centroidVel[particle] += velocities[copy][particle];
        }
        centroidVel[particle] *= (1.0 / numCopies);
        
        centroidKE += 0.5 * mass * centroidVel[particle].dot(centroidVel[particle]);
        ndof += 3;
    }
    
    if (ndof == 0)
        return;
    
    if (centroidKE <= 0.0)
        return;

    double r1 = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
    double rg = 0.0;
    for (int i = 0; i < ndof - 1; i++) {
        double rnd = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        rg += rnd * rnd;
    }
    double alpha2 = c1 + (kPerDof / centroidKE) * (1.0 - c1) * (r1 * r1 + rg)
                  + 2.0 * r1 * sqrt((kPerDof / centroidKE) * c1 * (1.0 - c1));
    if (alpha2 < 0.0)
        alpha2 = 0.0;
    double alpha = sqrt(alpha2);
    if ((r1 + sqrt(2.0 * centroidKE / kPerDof * c1 / (1.0 - c1))) < 0.0)
        alpha *= -1.0;

    double deltaAlpha = alpha - 1.0;
    for (int particle = 0; particle < numParticles; particle++) {
        if (system.getParticleMass(particle) == 0.0)
            continue;
        if (hybridMode && !isQuantumParticle[particle])
            continue;
        Vec3 deltaCentroid = centroidVel[particle] * deltaAlpha;
        for (int copy = 0; copy < numCopies; copy++)
            velocities[copy][particle] += deltaCentroid;
    }
}

void ReferenceIntegrateRPMDStepKernel::computeForces(ContextImpl& context, const RPMDIntegrator& integrator) {
    computeForcesWithMask(context, integrator, integrator.getIntegrationForceGroups());
}

void ReferenceIntegrateRPMDStepKernel::computeForcesWithMask(ContextImpl& context, const RPMDIntegrator& integrator, int groupMask) {
    const int totalCopies = positions.size();
    const int numParticles = positions[0].size();
    const System& system = context.getSystem();
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& vel = extractVelocities(context);
    vector<Vec3>& f = extractForces(context);
    int integrationMask = integrator.getIntegrationForceGroups();
    int nonContracted = groupsNotContracted & groupMask & integrationMask;
    
    // Compute forces from all groups that didn't have a specified contraction.
    
    for (int i = 0; i < totalCopies; i++) {
        pos = positions[i];
        // Hybrid classical particles: use bead-0 coordinates for force evaluation (i-PI frozen-ring / mean-field convention).
        if (hybridMode) {
            for (int j = 0; j < numParticles; j++) {
                if (isQuantumParticle[j] || system.getParticleMass(j) == 0.0)
                    continue;
                pos[j] = positions[0][j];
            }
        }
        vel = velocities[i];
        context.computeVirtualSites();
        Vec3 initialBox[3];
        context.getPeriodicBoxVectors(initialBox[0], initialBox[1], initialBox[2]);
        context.updateContextState();
        Vec3 finalBox[3];
        context.getPeriodicBoxVectors(finalBox[0], finalBox[1], finalBox[2]);
        if (initialBox[0] != finalBox[0] || initialBox[1] != finalBox[1] || initialBox[2] != finalBox[2])
            throw OpenMMException("Standard barostats cannot be used with RPMDIntegrator.  Use RPMDMonteCarloBarostat or RPMDStochasticCellRescalingBarostat instead.");
        positions[i] = pos;
        velocities[i] = vel;
        context.calcForcesAndEnergy(true, false, nonContracted);
        forces[i] = f;
    }
    
    // Now loop over contractions and compute forces from them.
    
    for (auto& g : groupsByCopies) {
        int copies = g.first;
        int groupFlags = g.second & groupMask & integrationMask;
        if (groupFlags == 0)
            continue;
        
        // Find the contracted positions.
        
        vector<complex<double>> q(totalCopies);
        const double scale1 = 1.0/totalCopies;
        for (int particle = 0; particle < numParticles; particle++) {
            for (int component = 0; component < 3; component++) {
                // Transform to the frequency domain, set high frequency components to zero, and transform back.
                
                for (int k = 0; k < totalCopies; k++)
                    q[k] = complex<double>(positions[k][particle][component], 0.0);
                pocketfft::c2c({(size_t) totalCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, q.data(), q.data(), 1.0, 0);
                if (copies > 1) {
                    int start = (copies+1)/2;
                    int end = totalCopies-copies+start;
                    for (int k = end; k < totalCopies; k++)
                        q[k-(totalCopies-copies)] = q[k];
                    pocketfft::c2c({(size_t) copies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, q.data(), q.data(), 1.0, 0);
                }
                for (int k = 0; k < copies; k++)
                    contractedPositions[k][particle][component] = scale1*q[k].real();
            }
        }
        
        // Compute forces.

        for (int i = 0; i < copies; i++) {
            pos = contractedPositions[i];
            if (hybridMode) {
                for (int j = 0; j < numParticles; j++) {
                    if (isQuantumParticle[j] || system.getParticleMass(j) == 0.0)
                        continue;
                    pos[j] = contractedPositions[0][j];
                }
            }
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
                    q[k] = complex<double>(contractedForces[k][particle][component], 0.0);
                if (copies > 1)
                    pocketfft::c2c({(size_t) copies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, true, q.data(), q.data(), 1.0, 0);
                int start = (copies+1)/2;
                int end = totalCopies-copies+start;
                for (int k = end; k < totalCopies; k++)
                    q[k] = q[k-(totalCopies-copies)];
                for (int k = start; k < end; k++)
                    q[k] = complex<double>(0, 0);
                pocketfft::c2c({(size_t) totalCopies}, {sizeof(complex<double>)}, {sizeof(complex<double>)}, {0}, false, q.data(), q.data(), 1.0, 0);
                for (int k = 0; k < totalCopies; k++)
                    forces[k][particle][component] += scale2*q[k].real();
            }
        }
    }

    if (integrator.getSuzukiChinEnabled() && groupMask == integrator.getIntegrationForceGroups()) {
        const System& system = context.getSystem();
        double beta = 1.0 / (BOLTZ * integrator.getTemperature());
        double betaP = beta / totalCopies;
        double hbar = 1.054571628e-34 * AVOGADRO / (1000 * 1e-12);
        double prefactor = (1.0/3.0) * betaP * betaP * hbar * hbar / 12.0;
        double epsilon = integrator.getSuzukiChinEpsilon();
        for (int copy = 0; copy < totalCopies; copy++) {
            vector<Vec3> origPos = positions[copy];
            vector<Vec3> origForces = forces[copy];
            for (int j = 0; j < numParticles; j++) {
                double mass = system.getParticleMass(j);
                if (mass == 0.0)
                    continue;
                Vec3 F = origForces[j];
                double fn = sqrt(F.dot(F));
                if (fn < 1e-20)
                    continue;
                Vec3 d = F * (1.0 / fn);
                positions[copy][j] = origPos[j] + d * epsilon;
            }
            pos = positions[copy];
            vel = velocities[copy];
            context.calcForcesAndEnergy(true, false, integrator.getIntegrationForceGroups());
            vector<Vec3> fPlus = f;
            for (int j = 0; j < numParticles; j++) {
                double mass = system.getParticleMass(j);
                if (mass == 0.0)
                    continue;
                Vec3 F = origForces[j];
                double fn = sqrt(F.dot(F));
                if (fn < 1e-20)
                    continue;
                Vec3 d = F * (1.0 / fn);
                positions[copy][j] = origPos[j] - d * epsilon;
            }
            pos = positions[copy];
            context.calcForcesAndEnergy(true, false, integrator.getIntegrationForceGroups());
            vector<Vec3> fMinus = f;
            positions[copy] = origPos;
            for (int j = 0; j < numParticles; j++) {
                double mass = system.getParticleMass(j);
                if (mass == 0.0)
                    continue;
                double coeff = -prefactor / mass / (2.0 * epsilon);
                forces[copy][j] = origForces[j] + (fPlus[j] - fMinus[j]) * coeff;
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
