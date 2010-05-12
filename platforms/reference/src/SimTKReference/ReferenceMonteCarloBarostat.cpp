
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

#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceMonteCarloBarostat.h"

using namespace std;

/**---------------------------------------------------------------------------------------

  Constructor

  --------------------------------------------------------------------------------------- */

ReferenceMonteCarloBarostat::ReferenceMonteCarloBarostat(int numAtoms, const vector<pair<int, int> >& constraints) {
    savedAtomPositions[0].resize(numAtoms);
    savedAtomPositions[1].resize(numAtoms);
    savedAtomPositions[2].resize(numAtoms);

    // First make a list of every other atom to which each atom is connect by a constraint.

    vector<vector<int> > atomBonds(numAtoms);
    for (int i = 0; i < (int) constraints.size(); i++) {
        atomBonds[constraints[i].first].push_back(constraints[i].second);
        atomBonds[constraints[i].second].push_back(constraints[i].first);
    }

    // Now tag atoms by which cluster they belong to.

    vector<int> atomCluster(numAtoms, -1);
    int numClusters = 0;
    for (int i = 0; i < numAtoms; i++)
        if (atomCluster[i] == -1)
            tagAtomsInCluster(i, numClusters++, atomCluster, atomBonds);
    clusters.resize(numClusters);
    for (int i = 0; i < numAtoms; i++)
        clusters[atomCluster[i]].push_back(i);
}

void ReferenceMonteCarloBarostat::tagAtomsInCluster(int atom, int cluster, vector<int>& atomCluster, vector<vector<int> >& atomBonds) {
    // Recursively tag atoms as belonging to a particular cluster.

    atomCluster[atom] = cluster;
    for (int i = 0; i < (int) atomBonds[atom].size(); i++)
        if (atomCluster[atomBonds[atom][i]] == -1)
            tagAtomsInCluster(atomBonds[atom][i], cluster, atomCluster, atomBonds);
}

/**---------------------------------------------------------------------------------------

  Destructor

  --------------------------------------------------------------------------------------- */

ReferenceMonteCarloBarostat::~ReferenceMonteCarloBarostat( ) {
}

/**---------------------------------------------------------------------------------------

  Apply the barostat at the start of a time step.

  @param atomPositions      atom positions
  @param boxSize            the periodic box dimensions
  @param scale              the factor by which to scale atom positions

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::applyBarostat(RealOpenMM** atomPositions, RealOpenMM* boxSize, RealOpenMM scale) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            savedAtomPositions[j][i] = atomPositions[i][j];

    // Loop over clusters.

    for (int i = 0; i < (int) clusters.size(); i++) {
        // Find the cluster center.

        RealOpenMM pos[3] = {0, 0, 0};
        for (int j = 0; j < (int) clusters[i].size(); j++) {
            RealOpenMM* atomPos = atomPositions[clusters[i][j]];
            pos[0] += atomPos[0];
            pos[1] += atomPos[1];
            pos[2] += atomPos[2];
        }
        pos[0] /= clusters[i].size();
        pos[1] /= clusters[i].size();
        pos[2] /= clusters[i].size();

        // Move it into the first periodic box.

        int xcell = (int) floor(pos[0]/boxSize[0]);
        int ycell = (int) floor(pos[1]/boxSize[1]);
        int zcell = (int) floor(pos[2]/boxSize[2]);
        float dx = xcell*boxSize[0];
        float dy = ycell*boxSize[1];
        float dz = zcell*boxSize[2];
        pos[0] -= dx;
        pos[1] -= dy;
        pos[2] -= dz;
        for (int j = 0; j < (int) clusters[i].size(); j++) {
            RealOpenMM* atomPos = atomPositions[clusters[i][j]];
            atomPos[0] -= dx;
            atomPos[1] -= dy;
            atomPos[2] -= dz;
        }

        // Now scale the position of the cluster center.

        dx = pos[0]*(scale-1);
        dy = pos[1]*(scale-1);
        dz = pos[2]*(scale-1);
        for (int j = 0; j < (int) clusters[i].size(); j++) {
            RealOpenMM* atomPos = atomPositions[clusters[i][j]];
            atomPos[0] += dx;
            atomPos[1] += dy;
            atomPos[2] += dz;
        }
    }
}

/**---------------------------------------------------------------------------------------

  Restore atom positions to what they were before applyBarostat() was called.

  @param atomPositions      atom positions

  --------------------------------------------------------------------------------------- */

void ReferenceMonteCarloBarostat::restorePositions(RealOpenMM** atomPositions) {
    int numAtoms = savedAtomPositions[0].size();
    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < 3; j++)
            atomPositions[i][j] = savedAtomPositions[j][i];
}
