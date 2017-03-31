
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

ReferenceMonteCarloBarostat::ReferenceMonteCarloBarostat(int numAtoms, const vector<vector<int> >& molecules) : molecules(molecules) {
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

  Apply the barostat at the start of a time step.

  @param atomPositions      atom positions
  @param boxVectors         the periodic box vectors
  @param scaleX             the factor by which to scale atom x-coordinates
  @param scaleY             the factor by which to scale atom y-coordinates
  @param scaleZ             the factor by which to scale atom z-coordinates

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::applyBarostat(vector<Vec3>& atomPositions, const Vec3* boxVectors, double scaleX, double scaleY, double scaleZ) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            savedAtomPositions[j][i] = atomPositions[i][j];

    // Loop over molecules.

    for (auto& molecule : molecules) {
        // Find the molecule center.

        Vec3 pos(0, 0, 0);
        for (int atom : molecule) {
            Vec3& atomPos = atomPositions[atom];
            pos += atomPos;
        }
        pos /= molecule.size();

        // Move it into the first periodic box.

        Vec3 newPos = pos;
        newPos -= boxVectors[2]*floor(newPos[2]/boxVectors[2][2]);
        newPos -= boxVectors[1]*floor(newPos[1]/boxVectors[1][1]);
        newPos -= boxVectors[0]*floor(newPos[0]/boxVectors[0][0]);

        // Now scale the position of the molecule center.

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

  Restore atom positions to what they were before applyBarostat() was called.

  @param atomPositions      atom positions

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::restorePositions(vector<Vec3>& atomPositions) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            atomPositions[i][j] = savedAtomPositions[j][i];
}
