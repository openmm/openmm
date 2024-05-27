
/* Portions copyright (c) 2006-2024 Stanford University and Simbios.
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
#include <algorithm>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceVariableStochasticDynamics.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param friction       friction coefficient
   @param temperature    temperature
   @param accuracy       required accuracy

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::ReferenceVariableStochasticDynamics(int numberOfAtoms,
                                                          double friction, double temperature,
                                                          double accuracy) :
           ReferenceDynamics(numberOfAtoms, 0.0f, temperature), friction(friction), _accuracy(accuracy) {
   xPrime.resize(numberOfAtoms);
   oldx.resize(numberOfAtoms);
   inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceVariableStochasticDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVariableStochasticDynamics::~ReferenceVariableStochasticDynamics() {
}

/**---------------------------------------------------------------------------------------

   Get the required accuracy

   @return accuracy

 --------------------------------------------------------------------------------------- */

double ReferenceVariableStochasticDynamics::getAccuracy() const {
    return _accuracy;
}

/**---------------------------------------------------------------------------------------

   Set the required accuracy

 --------------------------------------------------------------------------------------- */

void ReferenceVariableStochasticDynamics::setAccuracy(double accuracy) {
    _accuracy = accuracy;
}

/**---------------------------------------------------------------------------------------

   Get friction coefficient

   --------------------------------------------------------------------------------------- */

double ReferenceVariableStochasticDynamics::getFriction() const {
   return friction;
}

void ReferenceVariableStochasticDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& inverseMasses, double maxStepSize) {
   // Select the step size to use
    double error = 0;
    for (int i = 0; i < numberOfAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            double xerror = inverseMasses[i]*forces[i][j];
            error += xerror*xerror;
        }
    }
    error = sqrt(error/(numberOfAtoms*3));
    double dt = sqrt(getAccuracy()/error);
    if (getDeltaT() > 0.0f)
        dt = std::min(dt, getDeltaT()*2.0f); // For safety, limit how quickly dt can increase.
    if (dt > getDeltaT() && dt < 1.2f*getDeltaT())
        dt = getDeltaT(); // Keeping dt constant between steps improves the behavior of the integrator.
    if (dt > maxStepSize)
        dt = maxStepSize;
    setDeltaT(dt);
 
    // perform first update

    for (int i = 0; i < numberOfAtoms; i++)
        if (inverseMasses[i] != 0.0)
            velocities[i] += (dt*inverseMasses[i])*forces[i];
}

void ReferenceVariableStochasticDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& velocities, vector<double>& inverseMasses,
                                         vector<Vec3>& xPrime) {
    const double halfdt = 0.5*getDeltaT();
    const double kT = BOLTZ*getTemperature();
    const double friction = getFriction();
    const double vscale = exp(-getDeltaT()*friction);
    const double noisescale = sqrt(1-vscale*vscale);

    for (int i = 0; i < numberOfAtoms; i++) {
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            velocities[i] = vscale*velocities[i] + noisescale*sqrt(kT*inverseMasses[i])*Vec3(
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber(),
                    SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            xPrime[i] = xPrime[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing stochastic dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceVariableStochasticDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities,
                                          vector<Vec3>& forces, vector<double>& masses, double maxStepSize, double tolerance) {
    int numberOfAtoms = system.getNumParticles();
    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (getTimeStep() == 0) {
        // Invert masses

        for (int ii = 0; ii < numberOfAtoms; ii++) {
            if (masses[ii] == 0.0)
                inverseMasses[ii] = 0.0;
            else
                inverseMasses[ii] = 1/masses[ii];
        }
    }

    // 1st update

    updatePart1(numberOfAtoms, velocities, forces, inverseMasses, maxStepSize);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);

    // 2nd update

    updatePart2(numberOfAtoms, atomCoordinates, velocities, inverseMasses, xPrime);
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

    // copy xPrime -> atomCoordinates

    double invStepSize = 1.0/getDeltaT();
    for (int i = 0; i < numberOfAtoms; i++) {
        if (masses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])*invStepSize;
            atomCoordinates[i] = xPrime[i];
        }
    }

    getVirtualSites().computePositions(system, atomCoordinates);
    incrementTimeStep();
}
