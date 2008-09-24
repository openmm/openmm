#include "../../../tests/AssertionUtilities.h"
#include "../src/SimTKReference/ReferenceNeighborList.h"
#include "../src/sfmt/SFMT.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace OpenMM;

void testNeighborList()
{
    RealOpenMM* atomList[2];
    atomList[0] = new RealOpenMM[3];
    atomList[1] = new RealOpenMM[3];
    atomList[2] = new RealOpenMM[3];
    atomList[0][0] = 13.6f;
    atomList[0][1] = 0;
    atomList[0][2] = 0;
    atomList[1][0] = 0;
    atomList[1][1] = 0;
    atomList[1][2] = 0;
    vector<set<int> > exclusions(2);
    
    NeighborList neighborList;
    
    computeNeighborListNaive(neighborList, 2, atomList, exclusions, NULL, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListNaive(neighborList, 2, atomList, exclusions, NULL, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    computeNeighborListVoxelHash(neighborList, 2, atomList, exclusions, NULL, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListVoxelHash(neighborList, 2, atomList, exclusions, NULL, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    delete[] atomList[0];
    delete[] atomList[1];
    delete[] atomList[2];
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

void verifyNeighborList(NeighborList& list, int numAtoms, RealOpenMM** positions, const RealOpenMM* periodicBoxSize, double cutoff) {
    for (int i = 0; i < (int) list.size(); i++) {
        int atom1 = list[i].first;
        int atom2 = list[i].second;
        ASSERT(distance2(positions[atom1], positions[atom2], periodicBoxSize) <= cutoff*cutoff);
    }
    int count = 0;
    for (int i = 0; i < numAtoms; i++)
        for (int j = i+1; j < numAtoms; j++)
            if (distance2(positions[i], positions[j], periodicBoxSize) <= cutoff*cutoff)
                count++;
    ASSERT(count == list.size());
}

void testPeriodic() {
    const int numAtoms = 100;
    const double cutoff = 3.0;
    const RealOpenMM periodicBoxSize[3] = {20.0, 15.0, 22.0};
    RealOpenMM* atomList[numAtoms];
    init_gen_rand(0);
    for (int i = 0; i <numAtoms; i++) {
        atomList[i] = new RealOpenMM[3];
        atomList[i][0] = (RealOpenMM) (genrand_real2()*periodicBoxSize[0]*3);
        atomList[i][1] = (RealOpenMM) (genrand_real2()*periodicBoxSize[1]*3);
        atomList[i][2] = (RealOpenMM) (genrand_real2()*periodicBoxSize[2]*3);
    }
    vector<set<int> > exclusions(numAtoms);
    NeighborList neighborList;
    computeNeighborListNaive(neighborList, numAtoms, atomList, exclusions, periodicBoxSize, cutoff);
    verifyNeighborList(neighborList, numAtoms, atomList, periodicBoxSize, cutoff);
    computeNeighborListVoxelHash(neighborList, numAtoms, atomList, exclusions, periodicBoxSize, cutoff);
    verifyNeighborList(neighborList, numAtoms, atomList, periodicBoxSize, cutoff);
    for (int i = 0; i <numAtoms; i++)
        delete[] atomList[i];
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

