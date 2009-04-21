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

#include "../../../tests/AssertionUtilities.h"
#include "../src/SimTKReference/ReferenceNeighborList.h"
#include "../src/sfmt/SFMT.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace OpenMM;

void testNeighborList()
{
    RealOpenMM* particleList[2];
    particleList[0] = new RealOpenMM[3];
    particleList[1] = new RealOpenMM[3];
    // particleList[2] = new RealOpenMM[3];
    particleList[0][0] = 13.6f;
    particleList[0][1] = 0;
    particleList[0][2] = 0;
    particleList[1][0] = 0;
    particleList[1][1] = 0;
    particleList[1][2] = 0;
    vector<set<int> > exclusions(2);
    
    NeighborList neighborList;
    
    computeNeighborListNaive(neighborList, 2, particleList, exclusions, NULL, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListNaive(neighborList, 2, particleList, exclusions, NULL, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    computeNeighborListVoxelHash(neighborList, 2, particleList, exclusions, NULL, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListVoxelHash(neighborList, 2, particleList, exclusions, NULL, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    delete[] particleList[0];
    delete[] particleList[1];
    // delete[] particleList[2];
}

double periodicDifference(double val1, double val2, double period) {
    double diff = val1-val2;
    double base = floor(diff/period+0.5)*period;
    return diff-base;
}

double distance2(RealOpenMM* pos1, RealOpenMM* pos2, const RealOpenMM* periodicBoxSize) {
    double dx = periodicDifference(pos1[0], pos2[0], periodicBoxSize[0]);
    double dy = periodicDifference(pos1[1], pos2[1], periodicBoxSize[1]);
    double dz = periodicDifference(pos1[2], pos2[2], periodicBoxSize[2]);
    return dx*dx+dy*dy+dz*dz;
}

void verifyNeighborList(NeighborList& list, int numParticles, RealOpenMM** positions, const RealOpenMM* periodicBoxSize, double cutoff) {
    for (int i = 0; i < (int) list.size(); i++) {
        int particle1 = list[i].first;
        int particle2 = list[i].second;
        ASSERT(distance2(positions[particle1], positions[particle2], periodicBoxSize) <= cutoff*cutoff);
    }
    int count = 0;
    for (int i = 0; i < numParticles; i++)
        for (int j = i+1; j < numParticles; j++)
            if (distance2(positions[i], positions[j], periodicBoxSize) <= cutoff*cutoff)
                count++;
    ASSERT(count == list.size());
}

void testPeriodic() {
    const int numParticles = 100;
    const double cutoff = 3.0;
    const RealOpenMM periodicBoxSize[3] = {20.0, 15.0, 22.0};
    RealOpenMM* particleList[numParticles];
    init_gen_rand(0);
    for (int i = 0; i <numParticles; i++) {
        particleList[i] = new RealOpenMM[3];
        particleList[i][0] = (RealOpenMM) (genrand_real2()*periodicBoxSize[0]*3);
        particleList[i][1] = (RealOpenMM) (genrand_real2()*periodicBoxSize[1]*3);
        particleList[i][2] = (RealOpenMM) (genrand_real2()*periodicBoxSize[2]*3);
    }
    vector<set<int> > exclusions(numParticles);
    NeighborList neighborList;
    computeNeighborListNaive(neighborList, numParticles, particleList, exclusions, periodicBoxSize, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxSize, cutoff);
    computeNeighborListVoxelHash(neighborList, numParticles, particleList, exclusions, periodicBoxSize, cutoff);
    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxSize, cutoff);
    for (int i = 0; i <numParticles; i++)
        delete[] particleList[i];
}

int main() 
{
try {
    testNeighborList();
    testPeriodic();
    
    cout << "Test Passed" << endl;
    return 0;
}
catch (...) {
    cerr << "*** ERROR: Test Failed ***" << endl;
    return 1;
}
}

