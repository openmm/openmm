
/* Portions copyright (c) 2006-2019 Stanford University and Simbios.
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

#include <cstring>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceBAOABDynamics.h"
#include "ReferencePlatform.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"

using std::vector;
using namespace OpenMM;

ReferenceBAOABDynamics::ReferenceBAOABDynamics(int numberOfAtoms,
                                               double deltaT, double friction,
                                               double temperature) : 
           ReferenceDynamics(numberOfAtoms, deltaT, temperature), friction(friction) {
   xPrime.resize(numberOfAtoms);
   oldx.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

ReferenceBAOABDynamics::~ReferenceBAOABDynamics() {
}

double ReferenceBAOABDynamics::getFriction() const {
   return friction;
}

void ReferenceBAOABDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    for (int i = 0; i < numberOfAtoms; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (halfdt*inverseMasses[i])*forces[i];
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void ReferenceBAOABDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    const double kT = BOLTZ*getTemperature();
    const double friction = getFriction();
    const double vscale = exp(-getDeltaT()*friction);
    const double noisescale = sqrt(1-vscale*vscale);

    for (int i = 0; i < numberOfAtoms; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/halfdt;
            velocities[i] = vscale*velocities[i] + noisescale*sqrt(kT*inverseMasses[i])*Vec3(
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            atomCoordinates[i] = xPrime[i];
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void ReferenceBAOABDynamics::updatePart3(OpenMM::ContextImpl& context, int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                              vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& inverseMasses,
                                              vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    for (int i = 0; i < numberOfAtoms; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/halfdt;
            atomCoordinates[i] = xPrime[i];
        }
    }
    context.calcForcesAndEnergy(true, false);
    for (int i = 0; i < numberOfAtoms; i++) {
        if (inverseMasses[i] != 0.0)
            velocities[i] += (halfdt*inverseMasses[i])*forces[i];
    }
}

void ReferenceBAOABDynamics::update(ContextImpl& context, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities, vector<double>& masses, bool& forcesAreValid, double tolerance) {
    int numberOfAtoms = context.getSystem().getNumParticles();
    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (getTimeStep() == 0) {
        // Invert masses

        for (int ii = 0; ii < numberOfAtoms; ii++) {
            if (masses[ii] == 0.0)
                inverseMasses[ii] = 0.0;
            else
                inverseMasses[ii] = 1.0/masses[ii];
        }
    }
    if (!forcesAreValid) {
        context.calcForcesAndEnergy(true, false);
        forcesAreValid = true;
    }
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    vector<Vec3>& forces = *data->forces;

    // 1st update

    updatePart1(numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

    // 2nd update

    updatePart2(numberOfAtoms, atomCoordinates, velocities, inverseMasses, xPrime);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

    // 3rd update

    updatePart3(context, numberOfAtoms, atomCoordinates, velocities, forces, inverseMasses, xPrime);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);

    ReferenceVirtualSites::computePositions(context.getSystem(), atomCoordinates);
    incrementTimeStep();
}
