
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
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

#include <string.h>
#include <sstream>
#include <algorithm>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceVariableVerletDynamics.h"
#include "ReferenceVirtualSites.h"

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceVariableVerletDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         initial delta t for dynamics
   @param accuracy       required accuracy

   --------------------------------------------------------------------------------------- */

ReferenceVariableVerletDynamics::ReferenceVariableVerletDynamics(int numberOfAtoms, double accuracy) :
           ReferenceDynamics(numberOfAtoms, 0.0f, 0.0f), _accuracy(accuracy) {
    xPrime.resize(numberOfAtoms);
    inverseMasses.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceVariableVerletDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceVariableVerletDynamics::~ReferenceVariableVerletDynamics() {
}

/**---------------------------------------------------------------------------------------

   Get the required accuracy

   @return accuracy

 --------------------------------------------------------------------------------------- */

double ReferenceVariableVerletDynamics::getAccuracy() const {
    return _accuracy;
}

/**---------------------------------------------------------------------------------------

   Set the required accuracy

 --------------------------------------------------------------------------------------- */

void ReferenceVariableVerletDynamics::setAccuracy(double accuracy) {
    _accuracy = accuracy;
}

/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Verlet dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses
   @param maxStepSize         maximum time step

   --------------------------------------------------------------------------------------- */

void ReferenceVariableVerletDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities,
                                          vector<Vec3>& forces, vector<double>& masses, double maxStepSize, double tolerance) {
    // first-time-through initialization

    int numberOfAtoms = system.getNumParticles();
    if (getTimeStep() == 0) {
       // invert masses

       for (int ii = 0; ii < numberOfAtoms; ii++) {
          if (masses[ii] == 0.0)
              inverseMasses[ii] = 0.0;
          else
              inverseMasses[ii] = 1.0/masses[ii];
       }
    }

    double error = 0.0;
    for (int i = 0; i < numberOfAtoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            double xerror = inverseMasses[i]*forces[i][j];
            error += xerror*xerror;
        }
    }
    error = sqrt(error/(numberOfAtoms*3));
    double newStepSize = sqrt(getAccuracy()/error);
    if (getDeltaT() > 0.0f)
        newStepSize = std::min(newStepSize, getDeltaT()*2.0f); // For safety, limit how quickly dt can increase.
    if (newStepSize > getDeltaT() && newStepSize < 1.2f*getDeltaT())
        newStepSize = getDeltaT(); // Keeping dt constant between steps improves the behavior of the integrator.
    if (newStepSize > maxStepSize)
        newStepSize = maxStepSize;
    double vstep = 0.5f*(newStepSize+getDeltaT()); // The time interval by which to advance the velocities
    setDeltaT(newStepSize);
    for (int i = 0; i < numberOfAtoms; ++i) {
        if (masses[i] != 0.0)
            for (int j = 0; j < 3; ++j) {
                double vPrime = velocities[i][j] + inverseMasses[i]*forces[i][j]*vstep;
                xPrime[i][j] = atomCoordinates[i][j] + vPrime*getDeltaT();
            }
    }
    ReferenceConstraintAlgorithm* referenceConstraintAlgorithm = getReferenceConstraintAlgorithm();
    if (referenceConstraintAlgorithm)
        referenceConstraintAlgorithm->apply(atomCoordinates, xPrime, inverseMasses, tolerance);

   // Update the positions and velocities.

   double velocityScale = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i) {
       if (masses[i] != 0.0)
           for (int j = 0; j < 3; ++j) {
               velocities[i][j] = velocityScale*(xPrime[i][j] - atomCoordinates[i][j]);
               atomCoordinates[i][j] = xPrime[i][j];
           }
   }
   ReferenceVirtualSites::computePositions(system, atomCoordinates);
   incrementTimeStep();
}


