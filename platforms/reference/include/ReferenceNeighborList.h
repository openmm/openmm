#ifndef OPENMM_REFERENCE_NEIGHBORLIST_H_
#define OPENMM_REFERENCE_NEIGHBORLIST_H_

#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"
#include <set>
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
void OPENMM_EXPORT computeNeighborListNaive(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const std::vector<std::set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance = 0.0,
                              bool reportSymmetricPairs = false
                            );

// O(n) neighbor list method using voxel hash data structure
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void OPENMM_EXPORT computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations,
                              const std::vector<std::set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance = 0.0,
                              bool reportSymmetricPairs = false
                            );

} // namespace OpenMM

#endif // OPENMM_REFERENCE_NEIGHBORLIST_H_
