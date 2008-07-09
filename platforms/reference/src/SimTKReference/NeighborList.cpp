#include "NeighborList.h"
#include <set>
#include <map>
#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMM {

typedef std::vector<AtomIndex> AtomList;

// squared distance between two points
double compPairDistanceSquared(const Vec3& pos1, const Vec3& pos2) {
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];
    return dx*dx + dy*dy + dz*dz;
}

// Ridiculous O(n^2) version of neighbor list
// for pedagogical purposes and simplicity
void computeNeighborListNaive(
                              NeighborList& neighborList,
                              const AtomLocationList& atomLocations, 
                              double maxDistance,
                              double minDistance,
                              bool reportSymmetricPairs
                             )
{
    neighborList.clear();
    
    int nAtoms = atomLocations.size();

    double maxDistanceSquared = maxDistance * maxDistance;
    double minDistanceSquared = minDistance * minDistance;

    for (AtomIndex atomI = 0; atomI < (nAtoms - 1); ++atomI)
    {
        for (AtomIndex atomJ = atomI + 1; atomJ < nAtoms; ++atomJ)
        {
            double pairDistanceSquared = compPairDistanceSquared(atomLocations[atomI], atomLocations[atomJ]);
            if ( (pairDistanceSquared <= maxDistanceSquared)  && (pairDistanceSquared >= minDistanceSquared) )
            {
                neighborList.push_back( AtomPair(atomI, atomJ) );
                if (reportSymmetricPairs)
                    neighborList.push_back( AtomPair(atomI, atomJ) );
            }
        }
    }
}


class VoxelIndex 
{
public:
    VoxelIndex(int xx, int yy, int zz) : x(xx), y(yy), z(zz) {}

    // operator<() needed for map
    bool operator<(const VoxelIndex& other) const {
        if      (x < other.x) return true;
        else if (x > other.x) return false;
        else if (y < other.y) return true;
        else if (y > other.y) return false;
        else if (z < other.z) return true;
        else return false;
    }

    int x;
    int y;
    int z;
};


typedef std::pair<Vec3, AtomIndex> VoxelItem;
typedef std::vector< VoxelItem > Voxel;

class VoxelHash
{
public:
    VoxelHash(double vs) : voxelSize(vs) {}

    void insert(const AtomIndex& item, const Vec3& location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if ( voxelMap.find(voxelIndex) == voxelMap.end() ) voxelMap[voxelIndex] = Voxel(); 
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back( VoxelItem(location, item) );
    }


    VoxelIndex getVoxelIndex(const Vec3& location) const {
        int x = int(floor(location[0] / voxelSize));
        int y = int(floor(location[1] / voxelSize));
        int z = int(floor(location[2] / voxelSize));
        
        return VoxelIndex(x, y, z);
    }

    void getNeighbors(
            NeighborList& neighbors, 
            const VoxelItem& referencePoint, 
            bool reportSymmetricPairs,
            double maxDistance, 
            double minDistance) const 
    {

        // Loop over neighboring voxels
        // TODO use more clever selection of neighboring voxels
        assert(maxDistance > 0);
        assert(minDistance >= 0);
        assert(voxelSize > 0);

        const AtomIndex atomI = referencePoint.second;
        const Vec3& locationI = referencePoint.first;
        
        double maxDistanceSquared = maxDistance * maxDistance;
        double minDistanceSquared = minDistance * minDistance;

        int dIndex = int(maxDistance / voxelSize) + 1; // How may voxels away do we have to look?
        VoxelIndex centerVoxelIndex = getVoxelIndex(locationI);
        for (int x = centerVoxelIndex.x - dIndex; x <= centerVoxelIndex.x + dIndex; ++x)
        {
            for (int y = centerVoxelIndex.y - dIndex; y <= centerVoxelIndex.y + dIndex; ++y)
            {
                for (int z = centerVoxelIndex.z - dIndex; z <= centerVoxelIndex.z + dIndex; ++z)
                {
                    VoxelIndex voxelIndex(x, y, z);
                    if (voxelMap.find(voxelIndex) == voxelMap.end()) continue; // no such voxel; skip
                    const Voxel& voxel = voxelMap.find(voxelIndex)->second;
                    for (Voxel::const_iterator itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter)
                    {
                        const AtomIndex atomJ = itemIter->second;
                        const Vec3& locationJ = itemIter->first;
                        
                        // Ignore self hits
                        if (atomI == atomJ) continue;
                        
                        Vec3 dv = locationI - locationJ;
                        double dSquared = dv.dot(dv);
                        if (dSquared > maxDistanceSquared) continue;
                        if (dSquared < minDistanceSquared) continue;
                        
                        neighbors.push_back( AtomPair(atomI, atomJ) );
                        if (reportSymmetricPairs)
                            neighbors.push_back( AtomPair(atomJ, atomI) );
                    }
                }
            }
        }
    }

private:
    double voxelSize;
    std::map<VoxelIndex, Voxel> voxelMap;
};


// O(n) neighbor list method using voxel hash data structure
void computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              const AtomLocationList& atomLocations, 
                              double maxDistance,
                              double minDistance,
                              bool reportSymmetricPairs
                             )
{
    neighborList.clear();
    
    const int nAtoms = atomLocations.size();

    const double edgeSize = maxDistance; // TODO - adjust this as needed
    VoxelHash voxelHash(edgeSize);
    for (AtomIndex atomJ = 0; atomJ < nAtoms; ++atomJ) // use "j", because j > i for pairs
    {
        // 1) Find other atoms that are close to this one
        const Vec3& location = atomLocations[atomJ];
        voxelHash.getNeighbors( 
            neighborList, 
            VoxelItem(location, atomJ), 
            reportSymmetricPairs, 
            maxDistance, 
            minDistance);
            
        // 2) Add this atom to the voxelHash
        voxelHash.insert(atomJ, location);
    }
}

} // namespace OpenMM
