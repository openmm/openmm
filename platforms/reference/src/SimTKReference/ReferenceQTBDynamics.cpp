/* Portions copyright (c) 2025 Stanford University and the Authors.
 * Authors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKOpenMMUtilities.h"
#include "ReferenceQTBDynamics.h"
#include "ReferencePlatform.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <cmath>
#include <map>

#ifdef _MSC_VER
  #define POCKETFFT_NO_VECTORS
#endif
#include "pocketfft_hdronly.h"

using namespace OpenMM;
using namespace std;

ReferenceQTBDynamics::ReferenceQTBDynamics(const System& system, const QTBIntegrator& integrator) :
            ReferenceDynamics(system.getNumParticles(), integrator.getStepSize(), integrator.getTemperature()), friction(integrator.getFriction()), stepIndex(0) {
    segmentLength = (int) ceil(integrator.getSegmentLength()/integrator.getStepSize());
    int numParticles = system.getNumParticles();
    xPrime.resize(numParticles);
    oldx.resize(numParticles);
    inverseMasses.resize(numParticles);
    noise.resize(3*3*segmentLength*numParticles);
    randomForce.resize(segmentLength*numParticles);
    segmentVelocity.resize(segmentLength*numParticles);
    particleType.resize(numParticles);
    for (int i = 0; i < noise.size(); i++)
        noise[i] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

    // Record information about groups defined by particle types.

    map<int, int> typeIndex;
    map<int, double> massTable;
    const auto& types = integrator.getParticleTypes();
    double defaultAdaptationRate = integrator.getDefaultAdaptationRate();
    for (auto particle : types) {
        int type = particle.second;
        double mass = system.getParticleMass(particle.first);
        if (typeIndex.find(type) == typeIndex.end()) {
            typeIndex[type] = typeIndex.size();
            double rate = defaultAdaptationRate;
            const auto& typeRates = integrator.getTypeAdaptationRates();
            if (typeRates.find(type) != typeRates.end())
                rate = typeRates.at(type);
            typeAdaptationRate.push_back(rate);
            typeParticles.push_back(vector<int>());
            typeMass.push_back(mass);
            massTable[type] = mass;
        }
        if (mass != massTable[type])
            throw OpenMMException("QTBIntegrator: All particles of the same type must have the same mass");
        particleType[particle.first] = typeIndex[type];
        typeParticles[type].push_back(particle.first);
    }
    for (int i = 0; i < system.getNumParticles(); i++)
        if (types.find(i) == types.end()) {
            // This particle's type isn't set, so define a new type for it.
            particleType[i] = typeParticles.size();
            typeAdaptationRate.push_back(defaultAdaptationRate);
            typeParticles.push_back({i});
            typeMass.push_back(system.getParticleMass(i));
        }

    // Calculate the target energy distribution.

    numFreq = (3*segmentLength+1)/2;
    theta.resize(numFreq);
    cutoffFunction.resize(numFreq);
    double cutoff = integrator.getCutoffFrequency();
    double cutoffWidth = cutoff/100;
    for (int i = 1; i < numFreq; i++) {
        double w = M_PI*i/(numFreq*getDeltaT());
        cutoffFunction[i] = 1.0/(1.0+exp((w-cutoff)/cutoffWidth));
    }

    // Allocate space for adaptation.

    int numTypes = typeParticles.size();
    adaptedFriction.resize(numTypes, vector<double>(numFreq, friction));
}

ReferenceQTBDynamics::~ReferenceQTBDynamics() {
}

void ReferenceQTBDynamics::calcSpectrum(ThreadPool& threads) {
    // Compute the standard spectrum.

    double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    double kT = BOLTZ*getTemperature();
    theta[0] = kT;
    for (int i = 1; i < numFreq; i++) {
        double w = M_PI*i/(numFreq*getDeltaT());
        theta[i] = hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
    }

    // Compute the deconvolved version.

    auto C = [&](double w0, double w) {
        double t = w*w-w0*w0;
        return (friction/M_PI)*w0*w0/(t*t+friction*friction*w*w);
    };
    vector<vector<double> > D(numFreq, vector<double>(numFreq));
    vector<double> h(numFreq);
    vector<double> fcurrent(numFreq), fnext(numFreq);
    for (int i = 0; i < numFreq; i++)
        fcurrent[i] = 0.5*theta[i];
    double dw = M_PI/(numFreq*getDeltaT());
    threads.execute([&] (ThreadPool& threads, int threadIndex) {
        for (int i = threadIndex; i < numFreq; i += threads.getNumThreads()) {
            double wi = M_PI*(i+0.5)/(numFreq*getDeltaT());
            h[i] = 0.0;
            for (int j = 0; j < numFreq; j++) {
                double wj = M_PI*(j+0.5)/(numFreq*getDeltaT());
                h[i] += dw*C(wj, wi)*fcurrent[j];
                D[i][j] = 0.0;
                for (int k = 0; k < numFreq; k++) {
                    double wk = M_PI*(k+0.5)/(numFreq*getDeltaT());
                    D[i][j] += dw*C(wk, wi)*C(wk, wj);
                }
            }
        }
    });
    threads.waitForThreads();
    for (int iteration = 0; iteration < 20; iteration++) {
        for (int i = 0; i < numFreq; i++) {
            double denom = 0.0;
            for (int j = 0; j < numFreq; j++)
                denom += dw*D[i][j]*fcurrent[j];
            fnext[i] = fcurrent[i]*h[i]/denom;
        }
        fcurrent = fnext;
    }
    thetad = fnext;
    lastTemperature = getTemperature();
}

void ReferenceQTBDynamics::updatePart1(int numParticles, vector<Vec3>& velocities, vector<Vec3>& forces) {
    for (int i = 0; i < numParticles; i++) {
        segmentVelocity[i*segmentLength+stepIndex] = velocities[i];
        if (inverseMasses[i] != 0.0)
            velocities[i] += (getDeltaT()*inverseMasses[i])*forces[i];
    }
}

void ReferenceQTBDynamics::updatePart2(int numParticles, vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& velocities, vector<Vec3>& xPrime) {
    const double dt = getDeltaT();
    const double halfdt = 0.5*dt;
    const double vscale = exp(-dt*friction);

    for (int i = 0; i < numParticles; i++) {
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            velocities[i] = vscale*velocities[i] + inverseMasses[i]*dt*randomForce[segmentLength*i+stepIndex];
            xPrime[i] = xPrime[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void ReferenceQTBDynamics::updatePart3(OpenMM::ContextImpl& context, int numParticles, vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& velocities, vector<Vec3>& xPrime) {
    for (int i = 0; i < numParticles; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/getDeltaT();
            atomCoordinates[i] = xPrime[i];
        }
    }
}

void ReferenceQTBDynamics::update(ContextImpl& context, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
            vector<double>& masses, double tolerance, ThreadPool& threads) {
    int numParticles = context.getSystem().getNumParticles();
    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (getTimeStep() == 0) {
        // Invert masses

        for (int ii = 0; ii < numParticles; ii++) {
            if (masses[ii] == 0.0)
                inverseMasses[ii] = 0.0;
            else
                inverseMasses[ii] = 1.0/masses[ii];
        }
    }
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    vector<Vec3>& forces = *data->forces;
    if (stepIndex%segmentLength == 0) {
        if (lastTemperature != getTemperature() || thetad.size() == 0)
            calcSpectrum(threads);
        adaptFriction();
        generateNoise(numParticles, masses, threads);
        stepIndex = 0;
    }

    // 1st update

    updatePart1(numParticles, velocities, forces);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);

    // 2nd update

    updatePart2(numParticles, atomCoordinates, velocities, xPrime);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

    // 3rd update

    updatePart3(context, numParticles, atomCoordinates, velocities, xPrime);

    getVirtualSites().computePositions(context.getSystem(), atomCoordinates);
    incrementTimeStep();
    stepIndex += 1;
}

void ReferenceQTBDynamics::generateNoise(int numParticles, vector<double>& masses, ThreadPool& threads) {
    // Update the buffer of white noise.

    for (int base = 0; base < noise.size(); base += 3*segmentLength) {
        for (int i = 0; i < 2*segmentLength; i++)
            noise[base+i] = noise[base+i+segmentLength];
        for (int i = 0; i < segmentLength; i++)
            noise[base+2*segmentLength+i] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
    }

    // Generate the random force for the next segment.

    double dt = getDeltaT();
    vector<ptrdiff_t> realStride = {(ptrdiff_t) sizeof(double)};
    vector<ptrdiff_t> complexStride = {(ptrdiff_t) sizeof(complex<double>)};
    vector<size_t> shape = {(size_t) 3*segmentLength}, axes = {0};
    threads.execute([&] (ThreadPool& threads, int threadIndex) {
        vector<complex<double> > recipData(numFreq);
        vector<double> force(3*segmentLength);
        for (int particle = threadIndex; particle < numParticles; particle += threads.getNumThreads()) {
            int type = particleType[particle];
            for (int axis = 0; axis < 3; axis++) {
                double* data = &noise[(3*particle+axis)*3*segmentLength];
                pocketfft::r2c(shape, realStride, complexStride, axes, true, data, recipData.data(), 1.0, 1);
                for (int i = 0; i < numFreq; i++) {
                    double w = M_PI*i/(numFreq*dt);
                    double gamma = adaptedFriction[type][i];
                    double cw = (1 - 2*exp(-dt*friction)*cos(w*dt) + exp(-2*friction*dt)) / ((friction*friction+w*w)*dt*dt);
                    recipData[i] *= sqrt(cutoffFunction[i]*thetad[i]*cw*gamma/friction);
                }
                pocketfft::c2r(shape, complexStride, realStride, axes, false, recipData.data(), force.data(), 1.0/(3*segmentLength), 1);
                for (int i = 0; i < segmentLength; i++)
                    randomForce[particle*segmentLength+i][axis] = sqrt(2*masses[particle]*friction/dt)*force[segmentLength+i];
            }
        }
    });
    threads.waitForThreads();
}

void ReferenceQTBDynamics::adaptFriction() {
    vector<double> vel(3*segmentLength, 0.0), force(3*segmentLength, 0.0);
    vector<complex<double> > recipVel(numFreq), recipForce(numFreq);
    vector<double> dfdt(numFreq);
    vector<ptrdiff_t> realStride = {(ptrdiff_t) sizeof(double)};
    vector<ptrdiff_t> complexStride = {(ptrdiff_t) sizeof(complex<double>)};
    vector<size_t> shape = {(size_t) 3*segmentLength}, axes = {0};
    for (int type = 0; type < typeParticles.size(); type++) {
        for (int i = 0; i < dfdt.size(); i++)
            dfdt[i] = 0;
        for (int particle : typeParticles[type]) {
            for (int axis = 0; axis < 3; axis++) {
                // Compute the Fourier transformed velocity and force.

                for (int i = 0; i < segmentLength; i++) {
                    vel[i] = segmentVelocity[particle*segmentLength+i][axis];
                    force[i] = randomForce[particle*segmentLength+i][axis];
                }
                pocketfft::r2c(shape, realStride, complexStride, axes, true, vel.data(), recipVel.data(), 1.0, 1);
                pocketfft::r2c(shape, realStride, complexStride, axes, true, force.data(), recipForce.data(), 1.0, 1);

                // Compute the error in the fluctuation dissipation theorem.

                double mass = typeMass[type];
                for (int i = 0; i < numFreq; i++) {
                    double cvv = norm(recipVel[i]);
                    complex<double> cvf = recipVel[i]*conj(recipForce[i]);
                    dfdt[i] += mass*adaptedFriction[type][i]*cvv - cvf.real();
                }
            }
        }

        // Average over particles and axes, and update the friction.

        double scale = getDeltaT()/(3*typeParticles[type].size()*segmentLength);
        for (int i = 0; i < adaptedFriction[type].size(); i++)
            adaptedFriction[type][i] = max(0.0, adaptedFriction[type][i]-scale*typeAdaptationRate[type]*dfdt[i]);
    }
}

void ReferenceQTBDynamics::getAdaptedFriction(int particle, vector<double>& friction) const {
    friction = adaptedFriction[particleType[particle]];
}

void ReferenceQTBDynamics::createCheckpoint(std::ostream& stream) const {
    stream.write((char*) &stepIndex, sizeof(int));
    stream.write((char*) noise.data(), sizeof(double)*noise.size());
    stream.write((char*) randomForce.data(), sizeof(Vec3)*randomForce.size());
    stream.write((char*) segmentVelocity.data(), sizeof(Vec3)*segmentVelocity.size());
    for (int i = 0; i < adaptedFriction.size(); i++)
        stream.write((char*) adaptedFriction[i].data(), sizeof(double)*adaptedFriction[i].size());
}

void ReferenceQTBDynamics::loadCheckpoint(std::istream& stream) {
    stream.read((char*) &stepIndex, sizeof(int));
    stream.read((char*) noise.data(), sizeof(double)*noise.size());
    stream.read((char*) randomForce.data(), sizeof(Vec3)*randomForce.size());
    stream.read((char*) segmentVelocity.data(), sizeof(Vec3)*segmentVelocity.size());
    for (int i = 0; i < adaptedFriction.size(); i++)
        stream.read((char*) adaptedFriction[i].data(), sizeof(double)*adaptedFriction[i].size());
}
