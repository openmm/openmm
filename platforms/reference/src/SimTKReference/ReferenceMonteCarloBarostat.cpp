
/* Portions copyright (c) 2010 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceMonteCarloBarostat.h"

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

  Constructor

  --------------------------------------------------------------------------------------- */

ReferenceMonteCarloBarostat::ReferenceMonteCarloBarostat(int numAtoms, const vector<vector<int> >& molecules, const std::vector<double>& masses) :
            molecules(molecules), masses(masses) {
    savedAtomPositions[0].resize(numAtoms);
    savedAtomPositions[1].resize(numAtoms);
    savedAtomPositions[2].resize(numAtoms);
}

/**---------------------------------------------------------------------------------------

  Destructor

  --------------------------------------------------------------------------------------- */

ReferenceMonteCarloBarostat::~ReferenceMonteCarloBarostat() {
}

/**---------------------------------------------------------------------------------------

  Save the positions before applying the barostat.

  @param atomPositions      atom positions

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::savePositions(vector<Vec3>& atomPositions) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            savedAtomPositions[j][i] = atomPositions[i][j];
}

/**---------------------------------------------------------------------------------------

  Apply the barostat at the start of a time step.

  @param atomPositions      atom positions
  @param boxVectors         the periodic box vectors
  @param scaleX             the factor by which to scale atom x-coordinates
  @param scaleY             the factor by which to scale atom y-coordinates
  @param scaleZ             the factor by which to scale atom z-coordinates

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::applyBarostat(vector<Vec3>& atomPositions, const Vec3* boxVectors, double scaleX, double scaleY, double scaleZ) {
    // Loop over molecules.

    for (auto& molecule : molecules) {
        // Find the molecule center.

        Vec3 pos(0, 0, 0);
        for (int atom : molecule) {
            Vec3& atomPos = atomPositions[atom];
            pos += atomPos;
        }
        pos /= molecule.size();

        // Now scale the position of the molecule center.

        Vec3 newPos = pos;
        newPos[0] *= scaleX;
        newPos[1] *= scaleY;
        newPos[2] *= scaleZ;
        Vec3 offset = newPos-pos;
        for (int atom : molecule) {
            Vec3& atomPos = atomPositions[atom];
            atomPos += offset;
        }
    }
}

/**---------------------------------------------------------------------------------------

  Restore atom positions to what they were before savePositions() was called.

  @param atomPositions      atom positions

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::restorePositions(vector<Vec3>& atomPositions) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            atomPositions[i][j] = savedAtomPositions[j][i];
}

void ReferenceMonteCarloBarostat::computeMolecularKineticEnergy(const vector<Vec3>& velocities, vector<double>& ke, int components) {
    ke.resize(components);
    for (int i = 0; i < components; i++)
        ke[i] = 0.0;
    for (auto& molecule : molecules) {
        Vec3 molVel;
        double molMass = 0.0;
        for (int atom : molecule) {
            molVel += masses[atom]*velocities[atom];
            molMass += masses[atom];
        }
        molVel /= molecule.size();
        if (components == 1)
            ke[0] += 0.5*molMass*molVel.dot(molVel);
        else {
            ke[0] += 0.5*molMass*molVel[0]*molVel[0];
            ke[1] += 0.5*molMass*molVel[1]*molVel[1];
            ke[2] += 0.5*molMass*molVel[2]*molVel[2];
            if (components == 6) {
                ke[3] += 0.5*molMass*molVel[1]*molVel[0];
                ke[4] += 0.5*molMass*molVel[2]*molVel[0];
                ke[5] += 0.5*molMass*molVel[2]*molVel[1];
            }
        }
    }
}
