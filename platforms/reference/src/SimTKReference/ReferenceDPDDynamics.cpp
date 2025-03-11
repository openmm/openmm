
/* Portions copyright (c) 2025 Stanford University and Simbios.
 * Contributors: Pande Group
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
#include "ReferenceDPDDynamics.h"
#include "ReferencePlatform.h"
#include "ReferenceVirtualSites.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/DPDIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/DPDIntegratorUtilities.h"
#include "ReferenceForce.h"
#include <set>

using namespace std;
using namespace OpenMM;

ReferenceDPDDynamics::ReferenceDPDDynamics(const System& system, const DPDIntegrator& integrator) :
           ReferenceDynamics(system.getNumParticles(), integrator.getStepSize(), integrator.getTemperature()) {
    int numParticles = system.getNumParticles();
    xPrime.resize(numParticles);
    oldx.resize(numParticles);
    periodic = system.usesPeriodicBoundaryConditions();
    int numTypes;
    DPDIntegratorUtilities::createTypeTables(integrator, numParticles, numTypes, particleType, frictionTable, cutoffTable, maxCutoff);
}

ReferenceDPDDynamics::~ReferenceDPDDynamics() {
}

void ReferenceDPDDynamics::setPeriodicBoxVectors(OpenMM::Vec3* vectors) {
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];
}

double ReferenceDPDDynamics::getMaxCutoff() const {
    return maxCutoff;
}

void ReferenceDPDDynamics::updatePart1(int numParticles, vector<Vec3>& velocities, vector<Vec3>& forces) {
    for (int i = 0; i < numParticles; i++)
        if (inverseMasses[i] != 0.0)
            velocities[i] += (getDeltaT()*inverseMasses[i])*forces[i];
}

void ReferenceDPDDynamics::updatePart2(int numParticles, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
                                       vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    const double kT = BOLTZ*getTemperature();

    // First position update.

    for (int i = 0; i < numParticles; i++) {
        xPrime[i] = atomCoordinates[i];
        if (inverseMasses[i] != 0.0)
            xPrime[i] += velocities[i]*halfdt;
    }

    // Apply friction and noise to velocities.

    vector<set<int> > exclusions(numParticles);
    computeNeighborListVoxelHash(neighborList, numParticles, atomCoordinates, exclusions, periodicBoxVectors, periodic, maxCutoff, 0.0);
    for (auto& pair : neighborList) {
        int i = pair.first;
        int j = pair.second;
        if (masses[i] == 0.0 && masses[j] == 0.0)
            continue;
        int type1 = particleType[i];
        int type2 = particleType[j];
        double friction = frictionTable[type1][type2];
        double cutoff = cutoffTable[type1][type2];
        double deltaR[ReferenceForce::LastDeltaRIndex];
        if (periodic)
            ReferenceForce::getDeltaRPeriodic(xPrime[i], xPrime[j], periodicBoxVectors, deltaR);
        else
            ReferenceForce::getDeltaR(xPrime[i], xPrime[j], deltaR);
        double r = deltaR[ReferenceForce::RIndex];
        if (r >= cutoff)
            continue;
        double m = masses[i]*masses[j]/(masses[i]+masses[j]);
        double omega = 1.0-(r/cutoff);
        double vscale = exp(-getDeltaT()*2*friction*omega*omega);
        double noisescale = sqrt(1-vscale*vscale);
        Vec3 dir = Vec3(deltaR[ReferenceForce::XIndex], deltaR[ReferenceForce::YIndex], deltaR[ReferenceForce::ZIndex])/r;
        Vec3 v = velocities[j]-velocities[i];
        double dv = (1.0-vscale)*v.dot(dir) + noisescale*sqrt(kT/m)*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        if (masses[i] != 0.0)
            velocities[i] += (m/masses[i])*dv*dir;
        if (masses[j] != 0.0)
            velocities[j] -= (m/masses[j])*dv*dir;
    }

    // Second position update.

    for (int i = 0; i < numParticles; i++)
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = xPrime[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
}

void ReferenceDPDDynamics::updatePart3(OpenMM::ContextImpl& context, int numParticles, vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& velocities, vector<Vec3>& xPrime) {
    for (int i = 0; i < numParticles; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/getDeltaT();
            atomCoordinates[i] = xPrime[i];
        }
    }
}

void ReferenceDPDDynamics::update(ContextImpl& context, vector<Vec3>& atomCoordinates,
                                    vector<Vec3>& velocities, vector<double>& masses, double tolerance) {
    int numParticles = context.getSystem().getNumParticles();
    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (this->masses.size() == 0) {
        this->masses = masses;
        inverseMasses.resize(masses.size());
        for (int i = 0; i < masses.size(); i++) {
            if (masses[i] == 0.0)
                inverseMasses[i] = 0.0;
            else
                inverseMasses[i] = 1.0/masses[i];
        }
    }

    // 1st update

    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    updatePart1(numParticles, velocities, *data->forces);
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
}
