#ifndef OPENMM_REFERENCE_NEIGHBORLIST_H_
#define OPENMM_REFERENCE_NEIGHBORLIST_H_

#include "Vec3.h"
#include <vector>

namespace OpenMM {

typedef std::vector<Vec3> AtomLocationList;
typedef unsigned int AtomIndex;
typedef std::pair<AtomIndex, AtomIndex> AtomPair;
typedef std::vector<AtomPair>  NeighborList;

// Ridiculous O(n^2) version of neighbor list
// for pedagogical purposes and simplicity
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void computeNeighborListNaive(
                              NeighborList& neighborList,
                              const AtomLocationList& atomLocations, 
                              double maxDistance,
                              double minDistance = 0.0,
                              bool reportSymmetricPairs = false
                             );

// O(n) neighbor list method using voxel hash data structure
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              const AtomLocationList& atomLocations, 
                              double maxDistance,
                              double minDistance,
                              bool reportSymmetricPairs = false
                             );

} // namespace OpenMM

#endif // OPENMM_REFERENCE_NEIGHBORLIST_H_
