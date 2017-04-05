/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/**
 * This tests all the CPU implementation of neighbor list construction.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ThreadPool.h"
#include "AlignedArray.h"
#include "CpuNeighborList.h"
#include "CpuPlatform.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace std;

void testNeighborList(bool periodic, bool triclinic) {
    const int numParticles = 500;
    const float cutoff = 2.0f;
    Vec3 boxVectors[3];
    if (triclinic) {
        boxVectors[0] = Vec3(10, 0, 0);
        boxVectors[1] = Vec3(4, 9, 0);
        boxVectors[2] = Vec3(-3, -3.5, 11);
    }
    else {
        boxVectors[0] = Vec3(10, 0, 0);
        boxVectors[1] = Vec3(0, 9, 0);
        boxVectors[2] = Vec3(0, 0, 11);
    }
    const float boxSize[3] = {(float) boxVectors[0][0], (float) boxVectors[1][1], (float) boxVectors[2][2]};
    const int blockSize = 8;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    AlignedArray<float> positions(4*numParticles);
    for (int i = 0; i < 4*numParticles; i++)
        if (i%4 < 3)
            positions[i] = boxSize[i%4]*genrand_real2(sfmt);
    vector<set<int> > exclusions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        int num = min(i+1, 10);
        for (int j = 0; j < num; j++) {
            exclusions[i].insert(i-j);
            exclusions[i-j].insert(i);
        }
    }
    ThreadPool threads;
    CpuNeighborList neighborList(blockSize);
    neighborList.computeNeighborList(numParticles, positions, exclusions, boxVectors, periodic, cutoff, threads);
    
    // Convert the neighbor list to a set for faster lookup.
    
    set<pair<int, int> > neighbors;
    for (int i = 0; i < (int) neighborList.getSortedAtoms().size(); i++) {
        int blockIndex = i/blockSize;
        int indexInBlock = i-blockIndex*blockSize;
        char mask = 1<<indexInBlock;
        for (int j = 0; j < (int) neighborList.getBlockExclusions(blockIndex).size(); j++) {
            if ((neighborList.getBlockExclusions(blockIndex)[j] & mask) == 0) {
                int atom1 = neighborList.getSortedAtoms()[i];
                int atom2 = neighborList.getBlockNeighbors(blockIndex)[j];
                pair<int, int> entry = make_pair(min(atom1, atom2), max(atom1, atom2));
                ASSERT(neighbors.find(entry) == neighbors.end() && neighbors.find(make_pair(entry.second, entry.first)) == neighbors.end()); // No duplicates
                neighbors.insert(entry);
            }
        }
    }
    
    // Check each particle pair and figure out whether they should be in the neighbor list.

    for (int i = 0; i < numParticles; i++)
        for (int j = 0; j <= i; j++) {
            bool shouldInclude = (exclusions[i].find(j) == exclusions[i].end());
            Vec3 diff(positions[4*i]-positions[4*j], positions[4*i+1]-positions[4*j+1], positions[4*i+2]-positions[4*j+2]);
            if (periodic) {
                diff -= boxVectors[2]*floor(diff[2]/boxSize[2]+0.5);
                diff -= boxVectors[1]*floor(diff[1]/boxSize[1]+0.5);
                diff -= boxVectors[0]*floor(diff[0]/boxSize[0]+0.5);
            }
            if (diff.dot(diff) > cutoff*cutoff)
                shouldInclude = false;
            bool isIncluded = (neighbors.find(make_pair(i, j)) != neighbors.end() || neighbors.find(make_pair(j, i)) != neighbors.end());
            if (shouldInclude)
                ASSERT(isIncluded);
        }
}

int main() {
    try {
        if (!CpuPlatform::isProcessorSupported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 0;
        }
        testNeighborList(false, false);
        testNeighborList(true, false);
        testNeighborList(true, true);
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
