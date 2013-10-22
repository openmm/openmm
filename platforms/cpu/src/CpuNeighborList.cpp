#include "CpuNeighborList.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorize.h"
#include "hilbert.h"
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <smmintrin.h>

using namespace std;

namespace OpenMM {

const int CpuNeighborList::BlockSize = 4;

class VoxelIndex 
{
public:
    VoxelIndex(int xx, int yy, int zz) : x(xx), y(yy), z(zz) {
    }

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

typedef pair<const float*, int> VoxelItem;
typedef vector< VoxelItem > Voxel;

class CpuNeighborList::VoxelHash {
public:
    VoxelHash(float vsx, float vsy, float vsz, const float* periodicBoxSize, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), voxelSizeZ(vsz), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floorf(periodicBoxSize[0]/voxelSizeX+0.5f);
            ny = (int) floorf(periodicBoxSize[1]/voxelSizeY+0.5f);
            nz = (int) floorf(periodicBoxSize[2]/voxelSizeZ+0.5f);
            voxelSizeX = periodicBoxSize[0]/nx;
            voxelSizeY = periodicBoxSize[1]/ny;
            voxelSizeZ = periodicBoxSize[2]/nz;
        }
    }

    void insert(const int& item, const float* location) {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if (voxelMap.find(voxelIndex) == voxelMap.end())
            voxelMap[voxelIndex] = Voxel(); 
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

    void getNeighbors(vector<int>& neighbors, int blockIndex, fvec4 blockCenter, fvec4 blockWidth, const vector<int>& sortedAtoms, vector<char>& exclusions, float maxDistance, const vector<int> blockAtoms, const float* atomLocations) const {
        neighbors.resize(0);
        exclusions.resize(0);
        fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
        fvec4 invBoxSize(1/periodicBoxSize[0], 1/periodicBoxSize[1], 1/periodicBoxSize[2], 0);
        
        float maxDistanceSquared = maxDistance * maxDistance;
        float refineCutoff = maxDistance-max(max(blockWidth[0], blockWidth[1]), blockWidth[2]);
        float refineCutoffSquared = refineCutoff*refineCutoff;

        int dIndexX = int((maxDistance+blockWidth[0])/voxelSizeX) + 1; // How may voxels away do we have to look?
        int dIndexY = int((maxDistance+blockWidth[1])/voxelSizeY) + 1;
        int dIndexZ = int((maxDistance+blockWidth[2])/voxelSizeZ) + 1;
        float centerPos[4];
        blockCenter.store(centerPos);
        VoxelIndex centerVoxelIndex = getVoxelIndex(centerPos);
        int lastx = centerVoxelIndex.x+dIndexX;
        int lasty = centerVoxelIndex.y+dIndexY;
        int lastz = centerVoxelIndex.z+dIndexZ;
        if (usePeriodic) {
            lastx = min(lastx, centerVoxelIndex.x-dIndexX+nx-1);
            lasty = min(lasty, centerVoxelIndex.y-dIndexY+ny-1);
            lastz = min(lastz, centerVoxelIndex.z-dIndexZ+nz-1);
        }
        int lastSortedIndex = BlockSize*(blockIndex+1);
        VoxelIndex voxelIndex(0, 0, 0);
        for (int x = centerVoxelIndex.x - dIndexX; x <= lastx; ++x) {
            voxelIndex.x = x;
            if (usePeriodic)
                voxelIndex.x = (x < 0 ? x+nx : (x >= nx ? x-nx : x));
            for (int y = centerVoxelIndex.y - dIndexY; y <= lasty; ++y) {
                voxelIndex.y = y;
                if (usePeriodic)
                    voxelIndex.y = (y < 0 ? y+ny : (y >= ny ? y-ny : y));
                for (int z = centerVoxelIndex.z - dIndexZ; z <= lastz; ++z) {
                    voxelIndex.z = z;
                    if (usePeriodic)
                        voxelIndex.z = (z < 0 ? z+nz : (z >= nz ? z-nz : z));
                    const map<VoxelIndex, Voxel>::const_iterator voxelEntry = voxelMap.find(voxelIndex);
                    if (voxelEntry == voxelMap.end())
                        continue; // no such voxel; skip
                    const Voxel& voxel = voxelEntry->second;
                    for (Voxel::const_iterator itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter) {
                        const int sortedIndex = itemIter->second;

                        // Avoid duplicate entries.
                        if (sortedIndex >= lastSortedIndex)
                            break;
                        
                        fvec4 atomPos(itemIter->first);
                        fvec4 delta = atomPos-centerPos;
                        if (usePeriodic) {
                            fvec4 base = round(delta*invBoxSize)*boxSize;
                            delta = delta-base;
                        }
                        delta = max(0.0f, abs(delta)-blockWidth);
                        float dSquared = dot3(delta, delta);
                        if (dSquared > maxDistanceSquared)
                            continue;
                        
                        if (dSquared > refineCutoffSquared) {
                            // The distance is large enough that there might not be any actual interactions.
                            // Check individual atom pairs to be sure.
                            
                            bool any = false;
                            for (int k = 0; k < (int) blockAtoms.size(); k++) {
                                fvec4 pos1(&atomLocations[4*blockAtoms[k]]);
                                delta = atomPos-pos1;
                                if (usePeriodic) {
                                    fvec4 base = round(delta*invBoxSize)*boxSize;
                                    delta = delta-base;
                                }
                                float r2 = dot3(delta, delta);
                                if (r2 < maxDistanceSquared) {
                                    any = true;
                                    break;
                                }
                            }
                            if (!any)
                                continue;
                        }
                        
                        // Add this atom to the list of neighbors.
                        
                        neighbors.push_back(sortedAtoms[sortedIndex]);
                        if (sortedIndex < BlockSize*blockIndex)
                            exclusions.push_back(0);
                        else
                            exclusions.push_back(0xF & (0xF<<(sortedIndex-BlockSize*blockIndex)));
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
    map<VoxelIndex, Voxel> voxelMap;
};

class CpuNeighborList::ThreadData {
public:
    ThreadData(int index, CpuNeighborList& owner) : index(index), owner(owner) {
    }
    int index;
    CpuNeighborList& owner;
};

static void* threadBody(void* args) {
    CpuNeighborList::ThreadData& data = *reinterpret_cast<CpuNeighborList::ThreadData*>(args);
    data.owner.runThread(data.index);
    delete &data;
    return 0;
}

CpuNeighborList::CpuNeighborList() {
    isDeleted = false;
    numThreads = getNumProcessors();
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    thread.resize(numThreads);
    pthread_mutex_lock(&lock);
    waitCount = 0;
    for (int i = 0; i < numThreads; i++) {
        ThreadData* data = new ThreadData(i, *this);
        threadData.push_back(data);
        pthread_create(&thread[i], NULL, threadBody, data);
    }
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
}

CpuNeighborList::~CpuNeighborList() {
    isDeleted = true;
    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
    for (int i = 0; i < (int) thread.size(); i++)
        pthread_join(thread[i], NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
}

void CpuNeighborList::computeNeighborList(int numAtoms, const vector<float>& atomLocations, const vector<set<int> >& exclusions,
            const float* periodicBoxSize, bool usePeriodic, float maxDistance) {
    int numBlocks = (numAtoms+BlockSize-1)/BlockSize;
    blockNeighbors.resize(numBlocks);
    blockExclusions.resize(numBlocks);
    sortedAtoms.resize(numAtoms);
    
    // Sort the atoms based on a Hilbert curve.
    
    float minx = atomLocations[0], maxx = atomLocations[0];
    float miny = atomLocations[1], maxy = atomLocations[1];
    float minz = atomLocations[2], maxz = atomLocations[2];
    for (int i = 0; i < numAtoms; i++) {
        const float* pos = &atomLocations[4*i];
        if (pos[0] < minx)
            minx = pos[0];
        if (pos[1] < miny)
            miny = pos[1];
        if (pos[2] < minz)
            minz = pos[2];
        if (pos[0] > maxx)
            maxx = pos[0];
        if (pos[1] > maxy)
            maxy = pos[1];
        if (pos[2] > maxz)
            maxz = pos[2];
    }
    float binWidth = max(max(maxx-minx, maxy-miny), maxz-minz)/255.0f;
    float invBinWidth = 1.0f/binWidth;
    vector<pair<int, int> > atomBins(numAtoms);
    bitmask_t coords[3];
    for (int i = 0; i < numAtoms; i++) {
        const float* pos = &atomLocations[4*i];
        coords[0] = (bitmask_t) ((pos[0]-minx)*invBinWidth);
        coords[1] = (bitmask_t) ((pos[1]-miny)*invBinWidth);
        coords[2] = (bitmask_t) ((pos[2]-minz)*invBinWidth);
        int bin = (int) hilbert_c2i(3, 8, coords);
        atomBins[i] = pair<int, int>(bin, i);
    }
    sort(atomBins.begin(), atomBins.end());

    // Build the voxel hash.
    
    float edgeSizeX, edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = 0.5f*periodicBoxSize[0]/floorf(periodicBoxSize[0]/maxDistance);
        edgeSizeY = 0.5f*periodicBoxSize[1]/floorf(periodicBoxSize[1]/maxDistance);
        edgeSizeZ = 0.5f*periodicBoxSize[2]/floorf(periodicBoxSize[2]/maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, periodicBoxSize, usePeriodic);
    for (int i = 0; i < numAtoms; i++) {
        int atomIndex = atomBins[i].second;
        sortedAtoms[i] = atomIndex;
        voxelHash.insert(i, &atomLocations[4*atomIndex]);
    }
    
    // Record the parameters for the threads.
    
    this->voxelHash = &voxelHash;
    this->exclusions = &exclusions;
    this->atomLocations = &atomLocations[0];
    this->periodicBoxSize = periodicBoxSize;
    this->numAtoms = numAtoms;
    this->usePeriodic = usePeriodic;
    this->maxDistance = maxDistance;
    
    // Signal the threads to start running and wait for them to finish.
    
    pthread_mutex_lock(&lock);
    waitCount = 0;
    pthread_cond_broadcast(&startCondition);
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
    
    // Add padding atoms to fill up the last block.
    
    int numPadding = numBlocks*BlockSize-numAtoms;
    if (numPadding > 0) {
        char mask = (0xF0 >> numPadding) & 0xF;
        for (int i = 0; i < numPadding; i++)
            sortedAtoms.push_back(0);
        vector<char>& exc = blockExclusions[blockExclusions.size()-1];
        for (int i = 0; i < (int) exc.size(); i++)
            exc[i] |= mask;
    }
}

int CpuNeighborList::getNumBlocks() const {
    return sortedAtoms.size()/BlockSize;
}

const std::vector<int>& CpuNeighborList::getSortedAtoms() const {
    return sortedAtoms;
}

const std::vector<int>& CpuNeighborList::getBlockNeighbors(int blockIndex) const {
    return blockNeighbors[blockIndex];
}

const std::vector<char>& CpuNeighborList::getBlockExclusions(int blockIndex) const {
    return blockExclusions[blockIndex];
    
}

void CpuNeighborList::runThread(int index) {
    while (true) {
        // Wait for the signal to start running.
        
        pthread_mutex_lock(&lock);
        waitCount++;
        pthread_cond_signal(&endCondition);
        pthread_cond_wait(&startCondition, &lock);
        pthread_mutex_unlock(&lock);
        if (isDeleted)
            break;
        
        // Compute this thread's subset of neighbors.
        
        int numBlocks = blockNeighbors.size();
        vector<int> blockAtoms;
        for (int i = index; i < numBlocks; i += numThreads) {
            {
            int firstIndex = BlockSize*i;
            int atomsInBlock = min(BlockSize, numAtoms-firstIndex);
            blockAtoms.resize(atomsInBlock);
            for (int j = 0; j < atomsInBlock; j++)
                blockAtoms[j] = sortedAtoms[firstIndex+j];
            }

                        
            int firstIndex = BlockSize*i;
            fvec4 minPos(&atomLocations[4*sortedAtoms[firstIndex]]);
            fvec4 maxPos = minPos;
            int atomsInBlock = min(BlockSize, numAtoms-firstIndex);
            for (int j = 1; j < atomsInBlock; j++) {
                fvec4 pos(&atomLocations[4*sortedAtoms[firstIndex+j]]);
                minPos = min(minPos, pos);
                maxPos = max(maxPos, pos);
            }
            voxelHash->getNeighbors(blockNeighbors[i], i, (maxPos+minPos)*0.5f, (maxPos-minPos)*0.5f, sortedAtoms, blockExclusions[i], maxDistance, blockAtoms, atomLocations);
            
            // Record the exclusions for this block.
            
            for (int j = 0; j < atomsInBlock; j++) {
                const set<int>& atomExclusions = (*exclusions)[sortedAtoms[firstIndex+j]];
                char mask = 1<<j;
                for (int k = 0; k < (int) blockNeighbors[i].size(); k++) {
                    int atomIndex = blockNeighbors[i][k];
                    if (atomExclusions.find(atomIndex) != atomExclusions.end())
                        blockExclusions[i][k] |= mask;
                }
            }
        }
    }
}

} // namespace OpenMM
