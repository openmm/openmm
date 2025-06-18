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
#include <cmath>
#include <map>

#ifdef _MSC_VER
  #define POCKETFFT_NO_VECTORS
#endif
#include "pocketfft_hdronly.h"

using std::complex;
using std::map;
using std::vector;
using namespace OpenMM;

ReferenceQTBDynamics::ReferenceQTBDynamics(const System& system, const QTBIntegrator& integrator) :
           ReferenceDynamics(system.getNumParticles(), integrator.getStepSize(), integrator.getTemperature()), friction(integrator.getFriction()), stepIndex(0) {
    defaultAdaptationRate = integrator.getDefaultAdaptationRate();
    segmentLength = (int) ceil(integrator.getSegmentLength()/integrator.getStepSize());
    int numParticles = system.getNumParticles();
    xPrime.resize(numParticles);
    oldx.resize(numParticles);
    inverseMasses.resize(numParticles);
    noise.resize(3*3*segmentLength*numParticles);
    randomForce.resize(segmentLength*numParticles);
    for (int i = 0; i < noise.size(); i++)
        noise[i] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();

    // Record information about groups defined by particle types.

    map<int, int> typeIndex;
    for (auto particle : integrator.getParticleTypes()) {
        int type = particle.second;
        if (typeIndex.find(type) == typeIndex.end()) {
            typeIndex[type] = typeIndex.size();
            double rate = defaultAdaptationRate;
            const auto& typeRates = integrator.getTypeAdaptationRates();
            if (typeRates.find(type) != typeRates.end())
                rate = typeRates.at(type);
            typeAdaptationRate.push_back(rate);
            typeParticles.push_back(vector<int>());
        }
    }

    // Calculate the target energy distribution.

    theta.resize((3*segmentLength+1)/2);
    double hbar = 1.054571628e-34*AVOGADRO/(1000*1e-12);
    double kT = BOLTZ*getTemperature();
    double cutoff = integrator.getCutoffFrequency();
    double cutoffWidth = cutoff/100;
    theta[0] = kT;
    for (int i = 1; i < theta.size(); i++) {
        double w = M_PI*i/(theta.size()*getDeltaT());
        double f = 1.0/(1.0+exp((w-cutoff)/cutoffWidth));
        theta[i] = f*hbar*w*(0.5+1/(exp(hbar*w/kT)-1));
    }
}

ReferenceQTBDynamics::~ReferenceQTBDynamics() {
}

void ReferenceQTBDynamics::updatePart1(int numParticles, vector<Vec3>& velocities, vector<Vec3>& forces) {
    for (int i = 0; i < numParticles; i++)
        if (inverseMasses[i] != 0.0)
            velocities[i] += (getDeltaT()*inverseMasses[i])*forces[i];
}

void ReferenceQTBDynamics::updatePart2(int numParticles, vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& velocities, vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    const double vscale = exp(-getDeltaT()*friction);
    const double noisescale = sqrt(1-vscale*vscale);

    for (int i = 0; i < numParticles; i++) {
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            velocities[i] = vscale*velocities[i] + noisescale*randomForce[segmentLength*i+stepIndex];
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

void ReferenceQTBDynamics::update(ContextImpl& context, vector<Vec3>& atomCoordinates,
                                    vector<Vec3>& velocities, vector<double>& masses, double tolerance) {
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
        generateNoise(numParticles, masses);
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

void ReferenceQTBDynamics::generateNoise(int numParticles, vector<double>& masses) {
    // Update the buffer of white noise.

    for (int base = 0; base < noise.size(); base += 3*segmentLength) {
        for (int i = 0; i < 2*segmentLength; i++)
            noise[base+i] = noise[base+i+segmentLength];
        for (int i = 0; i < segmentLength; i++)
            noise[base+2*segmentLength+i] = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
    }

    // Generate the random force for the next segment.

    vector<complex<float> > recipData(theta.size());
    vector<float> force(3*segmentLength);
    vector<ptrdiff_t> realStride = {(ptrdiff_t) sizeof(float)};
    vector<ptrdiff_t> complexStride = {(ptrdiff_t) sizeof(complex<float>)};
    vector<size_t> shape = {0};
    for (int particle = 0; particle < numParticles; particle++) {
        for (int axis = 0; axis < 3; axis++) {
            float* data = &noise[(3*particle+axis)*3*segmentLength];
            pocketfft::r2c({(unsigned long) (3*segmentLength)}, realStride, complexStride, shape, true, data, recipData.data(), 1.0f, 1);
            for (int i = 0; i < theta.size(); i++)
                recipData[i] *= sqrt(theta[i]);
            pocketfft::c2r({(unsigned long) (3*segmentLength)}, complexStride, realStride, shape, false, recipData.data(), force.data(), (float) (1.0/(3*segmentLength)), 1);
            for (int i = 0; i < segmentLength; i++)
                randomForce[particle*segmentLength+i][axis] = force[segmentLength+i];
        }
    }
}
