#include "../src/SimTKReference/NeighborList.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace OpenMM;

void testNeighborList()
{
    AtomLocationList atomList;
    atomList.push_back(Vec3(13.6, 0, 0));
    atomList.push_back(Vec3(0, 0, 0));
    
    NeighborList neighborList;
    
    computeNeighborListNaive(neighborList, atomList, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListNaive(neighborList, atomList, 13.5, 0.01);
    assert(neighborList.size() == 0);
    
    computeNeighborListVoxelHash(neighborList, atomList, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListVoxelHash(neighborList, atomList, 13.5, 0.01);
    assert(neighborList.size() == 0);
}

int main() 
{
try {
    testNeighborList();
    
    cout << "Test Passed" << endl;
    return 0;
}
catch (...) {
    cerr << "*** ERROR: Test Failed ***" << endl;
    return 1;
}
}

