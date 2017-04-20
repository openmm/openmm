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

// squared distance between two points
static double compPairDistanceSquared(const Vec3& pos1, const Vec3& pos2, const Vec3* periodicBoxVectors, bool usePeriodic) {
    Vec3 diff = pos2-pos1;
    if (usePeriodic) {
        diff -= periodicBoxVectors[2]*floor(diff[2]/periodicBoxVectors[2][2]+0.5);
        diff -= periodicBoxVectors[1]*floor(diff[1]/periodicBoxVectors[1][1]+0.5);
        diff -= periodicBoxVectors[0]*floor(diff[0]/periodicBoxVectors[0][0]+0.5);
    }
    return diff.dot(diff);
}

// Ridiculous O(n^2) version of neighbor list
// for pedagogical purposes and simplicity
void OPENMM_EXPORT computeNeighborListNaive(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const vector<set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
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
            double pairDistanceSquared = compPairDistanceSquared(atomLocations[atomI], atomLocations[atomJ], periodicBoxVectors, usePeriodic);
            if ((pairDistanceSquared <= maxDistanceSquared)  && (pairDistanceSquared >= minDistanceSquared))
                if (exclusions[atomI].find(atomJ) == exclusions[atomI].end())
                {
                    neighborList.push_back(AtomPair(atomI, atomJ));
                    if (reportSymmetricPairs)
                        neighborList.push_back(AtomPair(atomI, atomJ));
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


typedef std::pair<const Vec3*, AtomIndex> VoxelItem;
typedef std::vector< VoxelItem > Voxel;

class VoxelHash
{
public:
    VoxelHash(double vsx, double vsy, double vsz, const Vec3* periodicBoxVectors, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), voxelSizeZ(vsz), periodicBoxVectors(periodicBoxVectors), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floor(periodicBoxVectors[0][0]/voxelSizeX+0.5);
            ny = (int) floor(periodicBoxVectors[1][1]/voxelSizeY+0.5);
            nz = (int) floor(periodicBoxVectors[2][2]/voxelSizeZ+0.5);
            voxelSizeX = periodicBoxVectors[0][0]/nx;
            voxelSizeY = periodicBoxVectors[1][1]/ny;
            voxelSizeZ = periodicBoxVectors[2][2]/nz;
        }
    }

    void insert(const AtomIndex& item, const Vec3& location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if (voxelMap.find(voxelIndex) == voxelMap.end())
            voxelMap[voxelIndex] = Voxel(); 
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back(VoxelItem(&location, item));
    }


    VoxelIndex getVoxelIndex(const Vec3& location) const {
        Vec3 r = location;
        if (usePeriodic) {
            r -= periodicBoxVectors[2]*floor(r[2]/periodicBoxVectors[2][2]);
            r -= periodicBoxVectors[1]*floor(r[1]/periodicBoxVectors[1][1]);
            r -= periodicBoxVectors[0]*floor(r[0]/periodicBoxVectors[0][0]);
        }
        int x = int(floor(r[0]/voxelSizeX));
        int y = int(floor(r[1]/voxelSizeY));
        int z = int(floor(r[2]/voxelSizeZ));

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
        const Vec3& locationI = *referencePoint.first;
        
        double maxDistanceSquared = maxDistance * maxDistance;
        double minDistanceSquared = minDistance * minDistance;

        int dIndexX = int(maxDistance / voxelSizeX) + 1; // How may voxels away do we have to look?
        int dIndexY = int(maxDistance / voxelSizeY) + 1;
        int dIndexZ = int(maxDistance / voxelSizeZ) + 1;
        VoxelIndex centerVoxelIndex = getVoxelIndex(locationI);
        int minz = centerVoxelIndex.z-dIndexZ;
        int maxz = centerVoxelIndex.z+dIndexZ;
        if (usePeriodic)
            maxz = min(maxz, minz+nz-1);
        for (int z = minz; z <= maxz; ++z)
        {
            int boxz = (int) floor((float) z/nz);
            int miny = centerVoxelIndex.y-dIndexY;
            int maxy = centerVoxelIndex.y+dIndexY;
            if (usePeriodic) {
                double yoffset = boxz*periodicBoxVectors[2][1]/voxelSizeY;
                miny -= (int) ceil(yoffset);
                maxy -= (int) floor(yoffset);
                maxy = min(maxy, miny+ny-1);
            }
            for (int y = miny; y <= maxy; ++y)
            {
                int boxy = (int) floor((float) y/ny);
                int minx = centerVoxelIndex.x-dIndexX;
                int maxx = centerVoxelIndex.x+dIndexX;
                if (usePeriodic) {
                    double xoffset = (boxy*periodicBoxVectors[1][0]+boxz*periodicBoxVectors[2][0])/voxelSizeX;
                    minx -= (int) ceil(xoffset);
                    maxx -= (int) floor(xoffset);
                    maxx = min(maxx, minx+nx-1);
                }
                for (int x = minx; x <= maxx; ++x)
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
                    for (auto& item : voxel)
                    {
                        const AtomIndex atomJ = item.second;
                        const Vec3& locationJ = *item.first;
                        
                        // Ignore self hits
                        if (atomI == atomJ) continue;
                        
                        double dSquared = compPairDistanceSquared(locationI, locationJ, periodicBoxVectors, usePeriodic);
                        if (dSquared > maxDistanceSquared) continue;
                        if (dSquared < minDistanceSquared) continue;
                        
                        // Ignore exclusions.
                        if (exclusions[atomI].find(atomJ) != exclusions[atomI].end()) continue;
                        
                        neighbors.push_back(AtomPair(atomI, atomJ));
                        if (reportSymmetricPairs)
                            neighbors.push_back(AtomPair(atomJ, atomI));
                    }
                }
            }
        }
    }

private:
    double voxelSizeX, voxelSizeY, voxelSizeZ;
    int nx, ny, nz;
    const Vec3* periodicBoxVectors;
    const bool usePeriodic;
    std::map<VoxelIndex, Voxel> voxelMap;
};


// O(n) neighbor list method using voxel hash data structure
void OPENMM_EXPORT computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations,
                              const vector<set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
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
        edgeSizeX = 0.5*periodicBoxVectors[0][0]/floor(periodicBoxVectors[0][0]/maxDistance);
        edgeSizeY = 0.5*periodicBoxVectors[1][1]/floor(periodicBoxVectors[1][1]/maxDistance);
        edgeSizeZ = 0.5*periodicBoxVectors[2][2]/floor(periodicBoxVectors[2][2]/maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, periodicBoxVectors, usePeriodic);
    for (AtomIndex atomJ = 0; atomJ < (AtomIndex) nAtoms; ++atomJ) // use "j", because j > i for pairs
    {
        // 1) Find other atoms that are close to this one
        const Vec3& location = atomLocations[atomJ];
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
