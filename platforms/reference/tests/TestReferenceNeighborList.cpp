/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/AssertionUtilities.h"
#include "ReferenceNeighborList.h"
#include "sfmt/SFMT.h"
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;
using namespace OpenMM;

void testNeighborList()
{
    vector<Vec3> particleList(2);
    particleList[0] = Vec3(13.6, 0, 0);
    particleList[1] = Vec3(0, 0, 0);
    vector<set<int> > exclusions(2);
    
    NeighborList neighborList;

    Vec3 boxVectors[3];
    computeNeighborListNaive(neighborList, 2, particleList, exclusions, boxVectors, false, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListNaive(neighborList, 2, particleList, exclusions, boxVectors, false, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    computeNeighborListVoxelHash(neighborList, 2, particleList, exclusions, boxVectors, false, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListVoxelHash(neighborList, 2, particleList, exclusions, boxVectors, false, 13.5, 0.01);
    assert(neighborList.size() == 0);
}

double distance2(Vec3& pos1, Vec3& pos2, const Vec3* periodicBoxVectors) {
    Vec3 diff = pos1-pos2;
    diff -= periodicBoxVectors[2]*floor(diff[2]/periodicBoxVectors[2][2]+0.5);
    diff -= periodicBoxVectors[1]*floor(diff[1]/periodicBoxVectors[1][1]+0.5);
    diff -= periodicBoxVectors[0]*floor(diff[0]/periodicBoxVectors[0][0]+0.5);
    return diff.dot(diff);
}

void verifyNeighborList(NeighborList& list, int numParticles, vector<Vec3>& positions, const Vec3* periodicBoxVectors, double cutoff) {
    for (int i = 0; i < (int) list.size(); i++) {
        int particle1 = list[i].first;
        int particle2 = list[i].second;
        ASSERT(distance2(positions[particle1], positions[particle2], periodicBoxVectors) <= cutoff*cutoff);
    }
    int count = 0;
    for (int i = 0; i < numParticles; i++)
        for (int j = i+1; j < numParticles; j++)
            if (distance2(positions[i], positions[j], periodicBoxVectors) <= cutoff*cutoff)
                count++;
    ASSERT_EQUAL(count, list.size());
}

void testPeriodic() {
    const int numParticles = 100;
    const double cutoff = 3.0;
    Vec3 periodicBoxVectors[3];
    periodicBoxVectors[0] = Vec3(20, 0, 0);
    periodicBoxVectors[1] = Vec3(0, 15, 0);
    periodicBoxVectors[2] = Vec3(0, 0, 22);
    vector<Vec3> particleList(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i <numParticles; i++) {
        particleList[i][0] = genrand_real2(sfmt)*periodicBoxVectors[0][0]*3;
        particleList[i][1] = genrand_real2(sfmt)*periodicBoxVectors[1][1]*3;
        particleList[i][2] = genrand_real2(sfmt)*periodicBoxVectors[2][2]*3;
    }
    vector<set<int> > exclusions(numParticles);
    NeighborList neighborList;
    computeNeighborListNaive(neighborList, numParticles, particleList, exclusions, periodicBoxVectors, true, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxVectors, cutoff);
    computeNeighborListVoxelHash(neighborList, numParticles, particleList, exclusions, periodicBoxVectors, true, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxVectors, cutoff);
}

void testTriclinic() {
    const int numParticles = 1000;
    const double cutoff = 3.0;
    Vec3 periodicBoxVectors[3];
    periodicBoxVectors[0] = Vec3(20, 0, 0);
    periodicBoxVectors[1] = Vec3(5, 15, 0);
    periodicBoxVectors[2] = Vec3(-3, -7, 22);
    vector<Vec3> particleList(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i <numParticles; i++) {
        particleList[i][0] = genrand_real2(sfmt)*periodicBoxVectors[0][0]*3;
        particleList[i][1] = genrand_real2(sfmt)*periodicBoxVectors[1][1]*3;
        particleList[i][2] = genrand_real2(sfmt)*periodicBoxVectors[2][2]*3;
    }
    vector<set<int> > exclusions(numParticles);
    NeighborList neighborList;
    computeNeighborListNaive(neighborList, numParticles, particleList, exclusions, periodicBoxVectors, true, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxVectors, cutoff);
    computeNeighborListVoxelHash(neighborList, numParticles, particleList, exclusions, periodicBoxVectors, true, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxVectors, cutoff);
}

int main() 
{
    try {
        testNeighborList();
        testPeriodic();
        testTriclinic();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

