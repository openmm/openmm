#include "ReferenceNeighborList.h"
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace std;

namespace OpenMM {

typedef std::vector<AtomIndex> AtomList;

static double periodicDifference(double val1, double val2, double period) {
    double diff = val1-val2;
    double base = floor(diff/period+0.5)*period;
    return diff-base;
}

// squared distance between two points
static double compPairDistanceSquared(const RealVec& pos1, const RealVec& pos2, const RealVec& periodicBoxSize, bool usePeriodic) {
    double dx, dy, dz;
    if (!usePeriodic) {
        dx = pos2[0] - pos1[0];
        dy = pos2[1] - pos1[1];
        dz = pos2[2] - pos1[2];
    }
    else {
        dx = periodicDifference(pos2[0], pos1[0], periodicBoxSize[0]);
        dy = periodicDifference(pos2[1], pos1[1], periodicBoxSize[1]);
        dz = periodicDifference(pos2[2], pos1[2], periodicBoxSize[2]);
    }
    return dx*dx + dy*dy + dz*dz;
}

// Ridiculous O(n^2) version of neighbor list
// for pedagogical purposes and simplicity
void OPENMM_EXPORT computeNeighborListNaive(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const vector<set<int> >& exclusions,
                              const RealVec& periodicBoxSize,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance,
                              bool reportSymmetricPairs) {
    neighborList.clear();
    
    double maxDistanceSquared = maxDistance * maxDistance;
    double minDistanceSquared = minDistance * minDistance;

    for (AtomIndex atomI = 0; atomI < (AtomIndex) (nAtoms - 1); ++atomI)
    {
        for (AtomIndex atomJ = atomI + 1; atomJ < (AtomIndex) nAtoms; ++atomJ)
        {
            double pairDistanceSquared = compPairDistanceSquared(atomLocations[atomI], atomLocations[atomJ], periodicBoxSize, usePeriodic);
            if ( (pairDistanceSquared <= maxDistanceSquared)  && (pairDistanceSquared >= minDistanceSquared))
                if (exclusions[atomI].find(atomJ) == exclusions[atomI].end())
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


typedef std::pair<const RealVec*, AtomIndex> VoxelItem;
typedef std::vector< VoxelItem > Voxel;

class VoxelHash
{
public:
    VoxelHash(double vsx, double vsy, double vsz, const RealVec& periodicBoxSize, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), voxelSizeZ(vsz), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floor(periodicBoxSize[0]/voxelSizeX+0.5);
            ny = (int) floor(periodicBoxSize[1]/voxelSizeY+0.5);
            nz = (int) floor(periodicBoxSize[2]/voxelSizeZ+0.5);
            voxelSizeX = periodicBoxSize[0]/nx;
            voxelSizeY = periodicBoxSize[1]/ny;
            voxelSizeZ = periodicBoxSize[2]/nz;
        }
    }

    void insert(const AtomIndex& item, const RealVec& location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if (voxelMap.find(voxelIndex) == voxelMap.end())
            voxelMap[voxelIndex] = Voxel(); 
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back(VoxelItem(&location, item));
    }


    VoxelIndex getVoxelIndex(const RealVec& location) const {
        double xperiodic, yperiodic, zperiodic;
        if (!usePeriodic) {
            xperiodic = location[0];
            yperiodic = location[1];
            zperiodic = location[2];
        }
        else {
            xperiodic = location[0]-periodicBoxSize[0]*floor(location[0]/periodicBoxSize[0]);
            yperiodic = location[1]-periodicBoxSize[1]*floor(location[1]/periodicBoxSize[1]);
            zperiodic = location[2]-periodicBoxSize[2]*floor(location[2]/periodicBoxSize[2]);
        }
        int x = int(floor(xperiodic / voxelSizeX));
        int y = int(floor(yperiodic / voxelSizeY));
        int z = int(floor(zperiodic / voxelSizeZ));
        
        return VoxelIndex(x, y, z);
    }

    void getNeighbors(
            NeighborList& neighbors, 
            const VoxelItem& referencePoint, 
            const vector<set<int> >& exclusions,
            bool reportSymmetricPairs,
            double maxDistance, 
            double minDistance) const 
    {

        // Loop over neighboring voxels
        // TODO use more clever selection of neighboring voxels
        assert(maxDistance > 0);
        assert(minDistance >= 0);
        assert(voxelSizeX > 0);
        assert(voxelSizeY > 0);
        assert(voxelSizeZ > 0);

        const AtomIndex atomI = referencePoint.second;
        const RealVec& locationI = *referencePoint.first;
        
        double maxDistanceSquared = maxDistance * maxDistance;
        double minDistanceSquared = minDistance * minDistance;

        int dIndexX = int(maxDistance / voxelSizeX) + 1; // How may voxels away do we have to look?
        int dIndexY = int(maxDistance / voxelSizeY) + 1;
        int dIndexZ = int(maxDistance / voxelSizeZ) + 1;
        VoxelIndex centerVoxelIndex = getVoxelIndex(locationI);
        int lastx = centerVoxelIndex.x+dIndexX;
        int lasty = centerVoxelIndex.y+dIndexY;
        int lastz = centerVoxelIndex.z+dIndexZ;
        if (usePeriodic) {
            lastx = min(lastx, centerVoxelIndex.x-dIndexX+nx-1);
            lasty = min(lasty, centerVoxelIndex.y-dIndexY+ny-1);
            lastz = min(lastz, centerVoxelIndex.z-dIndexZ+nz-1);
        }
        for (int x = centerVoxelIndex.x - dIndexX; x <= lastx; ++x)
        {
            for (int y = centerVoxelIndex.y - dIndexY; y <= lasty; ++y)
            {
                for (int z = centerVoxelIndex.z - dIndexZ; z <= lastz; ++z)
                {
                    VoxelIndex voxelIndex(x, y, z);
                    if (usePeriodic) {
                        voxelIndex.x = (x+nx)%nx;
                        voxelIndex.y = (y+ny)%ny;
                        voxelIndex.z = (z+nz)%nz;
                    }
                    const map<VoxelIndex, Voxel>::const_iterator voxelEntry = voxelMap.find(voxelIndex);
                    if (voxelEntry == voxelMap.end()) continue; // no such voxel; skip
                    const Voxel& voxel = voxelEntry->second;
                    for (Voxel::const_iterator itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter)
                    {
                        const AtomIndex atomJ = itemIter->second;
                        const RealVec& locationJ = *itemIter->first;
                        
                        // Ignore self hits
                        if (atomI == atomJ) continue;
                        
                        double dSquared = compPairDistanceSquared(locationI, locationJ, periodicBoxSize, usePeriodic);
                        if (dSquared > maxDistanceSquared) continue;
                        if (dSquared < minDistanceSquared) continue;
                        
                        // Ignore exclusions.
                        if (exclusions[atomI].find(atomJ) != exclusions[atomI].end()) continue;
                        
                        neighbors.push_back( AtomPair(atomI, atomJ) );
                        if (reportSymmetricPairs)
                            neighbors.push_back( AtomPair(atomJ, atomI) );
                    }
                }
            }
        }
    }

private:
    double voxelSizeX, voxelSizeY, voxelSizeZ;
    int nx, ny, nz;
    const RealVec& periodicBoxSize;
    const bool usePeriodic;
    std::map<VoxelIndex, Voxel> voxelMap;
};


// O(n) neighbor list method using voxel hash data structure
void OPENMM_EXPORT computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const vector<set<int> >& exclusions,
                              const RealVec& periodicBoxSize,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance,
                              bool reportSymmetricPairs
                             )
{
    neighborList.clear();

    double edgeSizeX, edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = 0.5*periodicBoxSize[0]/floor(periodicBoxSize[0]/maxDistance);
        edgeSizeY = 0.5*periodicBoxSize[1]/floor(periodicBoxSize[1]/maxDistance);
        edgeSizeZ = 0.5*periodicBoxSize[2]/floor(periodicBoxSize[2]/maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, periodicBoxSize, usePeriodic);
    for (AtomIndex atomJ = 0; atomJ < (AtomIndex) nAtoms; ++atomJ) // use "j", because j > i for pairs
    {
        // 1) Find other atoms that are close to this one
        const RealVec& location = atomLocations[atomJ];
        voxelHash.getNeighbors( 
            neighborList, 
            VoxelItem(&location, atomJ),
            exclusions,
            reportSymmetricPairs, 
            maxDistance, 
            minDistance);
            
        // 2) Add this atom to the voxelHash
        voxelHash.insert(atomJ, location);
    }
}

} // namespace OpenMM
