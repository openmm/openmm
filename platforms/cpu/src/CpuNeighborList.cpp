#include "CpuNeighborList.h"
#include <set>
#include <map>
#include <cmath>

using namespace std;

namespace OpenMM {

static float periodicDifference(float val1, float val2, float period) {
    float diff = val1-val2;
    float base = floorf(diff/period+0.5f)*period;
    return diff-base;
}

// squared distance between two points
static float compPairDistanceSquared(const float* pos1, const float* pos2, const float* periodicBoxSize, bool usePeriodic) {
    float dx, dy, dz;
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

typedef std::pair<const float*, int> VoxelItem;
typedef std::vector< VoxelItem > Voxel;

class VoxelHash
{
public:
    VoxelHash(float vsx, float vsy, float vsz, const float* periodicBoxSize, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), voxelSizeZ(vsz), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floorf(periodicBoxSize[0]/voxelSizeX+0.5f);
            ny = (int) floorf(periodicBoxSize[1]/voxelSizeY+0.5f);
            nz = (int) floorf(periodicBoxSize[2]/voxelSizeZ+0.5f);
        }
    }

    void insert(const int& item, const float* location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if (voxelMap.find(voxelIndex) == voxelMap.end()) voxelMap[voxelIndex] = Voxel(); 
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back(VoxelItem(location, item));
    }


    VoxelIndex getVoxelIndex(const float* location) const {
        float xperiodic, yperiodic, zperiodic;
        if (!usePeriodic) {
            xperiodic = location[0];
            yperiodic = location[1];
            zperiodic = location[2];
        }
        else {
            xperiodic = location[0]-periodicBoxSize[0]*floorf(location[0]/periodicBoxSize[0]);
            yperiodic = location[1]-periodicBoxSize[1]*floorf(location[1]/periodicBoxSize[1]);
            zperiodic = location[2]-periodicBoxSize[2]*floorf(location[2]/periodicBoxSize[2]);
        }
        int x = int(floorf(xperiodic / voxelSizeX));
        int y = int(floorf(yperiodic / voxelSizeY));
        int z = int(floorf(zperiodic / voxelSizeZ));
        
        return VoxelIndex(x, y, z);
    }

    void getNeighbors(
            vector<pair<int, int> >& neighbors, 
            const VoxelItem& referencePoint, 
            const vector<set<int> >& exclusions,
            bool reportSymmetricPairs,
            float maxDistance, 
            float minDistance) const 
    {

        // Loop over neighboring voxels
        // TODO use more clever selection of neighboring voxels

        const int atomI = referencePoint.second;
        const float* locationI = referencePoint.first;
        
        float maxDistanceSquared = maxDistance * maxDistance;
        float minDistanceSquared = minDistance * minDistance;

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
                    if (voxelMap.find(voxelIndex) == voxelMap.end()) continue; // no such voxel; skip
                    const Voxel& voxel = voxelMap.find(voxelIndex)->second;
                    for (Voxel::const_iterator itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter)
                    {
                        const int atomJ = itemIter->second;
                        const float* locationJ = itemIter->first;
                        
                        // Ignore self hits
                        if (atomI == atomJ) continue;
                        
                        // Ignore exclusions.
                        if (exclusions[atomI].find(atomJ) != exclusions[atomI].end()) continue;
                        
                        float dSquared = compPairDistanceSquared(locationI, locationJ, periodicBoxSize, usePeriodic);
                        if (dSquared > maxDistanceSquared) continue;
                        if (dSquared < minDistanceSquared) continue;
                        
                        neighbors.push_back(make_pair(atomI, atomJ));
                        if (reportSymmetricPairs)
                            neighbors.push_back(make_pair(atomJ, atomI));
                    }
                }
            }
        }
    }

private:
    float voxelSizeX, voxelSizeY, voxelSizeZ;
    int nx, ny, nz;
    const float* periodicBoxSize;
    const bool usePeriodic;
    std::map<VoxelIndex, Voxel> voxelMap;
};


// O(n) neighbor list method using voxel hash data structure
void CpuNeighborList::computeNeighborList(
                              int nAtoms,
                              const vector<float>& atomLocations, 
                              const vector<set<int> >& exclusions,
                              const float* periodicBoxSize,
                              bool usePeriodic,
                              float maxDistance,
                              float minDistance,
                              bool reportSymmetricPairs)
{
    neighbors.clear();

    float edgeSizeX, edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = periodicBoxSize[0]/floorf(periodicBoxSize[0]/maxDistance);
        edgeSizeY = periodicBoxSize[1]/floorf(periodicBoxSize[1]/maxDistance);
        edgeSizeZ = periodicBoxSize[2]/floorf(periodicBoxSize[2]/maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, periodicBoxSize, usePeriodic);
    for (int atomJ = 0; atomJ < (int) nAtoms; ++atomJ) // use "j", because j > i for pairs
    {
        // 1) Find other atoms that are close to this one
        const float location[3] = {atomLocations[4*atomJ], atomLocations[4*atomJ+1], atomLocations[4*atomJ+2]};
        voxelHash.getNeighbors(
            neighbors, 
            VoxelItem(location, atomJ),
            exclusions,
            reportSymmetricPairs, 
            maxDistance, 
            minDistance);
            
        // 2) Add this atom to the voxelHash
        voxelHash.insert(atomJ, location);
    }
}

const vector<pair<int, int> >& CpuNeighborList::getNeighbors() {
    return neighbors;
}

} // namespace OpenMM
