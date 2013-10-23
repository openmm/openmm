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
    VoxelIndex(int x, int y) : x(x), y(y) {
    }
    int x;
    int y;
};

typedef pair<const float*, int> VoxelItem;

/**
 * This data structure organizes the particles spatially.  It divides them into bins along the x and y axes,
 * then sorts each bin along the z axis so ranges can be identified quickly with a binary search.
 */
class CpuNeighborList::Voxels {
public:
    Voxels(float vsx, float vsy, float minx, float maxx, float miny, float maxy, const float* periodicBoxSize, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), minx(minx), maxx(maxx), miny(miny), maxy(maxy), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floorf(periodicBoxSize[0]/voxelSizeX+0.5f);
            ny = (int) floorf(periodicBoxSize[1]/voxelSizeY+0.5f);
            voxelSizeX = periodicBoxSize[0]/nx;
            voxelSizeY = periodicBoxSize[1]/ny;
        }
        else {
            nx = max(1, (int) floorf((maxx-minx)/voxelSizeX+0.5f));
            ny = max(1, (int) floorf((maxy-miny)/voxelSizeY+0.5f));
            if (maxx > minx)
                voxelSizeX = (maxx-minx)/nx;
            if (maxy > miny)
                voxelSizeY = (maxy-miny)/ny;
        }
        bins.resize(nx);
        for (int i = 0; i < nx; i++) {
            bins[i].resize(ny);
            for (int j = 0; j < ny; j++)
                bins[i][j].resize(0);
        }
    }

    /**
     * Insert a particle into the voxel data structure.
     */
    void insert(const int& atom, const float* location) {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        bins[voxelIndex.x][voxelIndex.y].push_back(make_pair(location[2], VoxelItem(location, atom)));
    }
    
    /**
     * Sort the particles in each voxel by z coordinate.
     */
    void sortItems() {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                sort(bins[i][j].begin(), bins[i][j].end());
    }
    
    /**
     * Find the index of the first particle in voxel (x,y) whose z coordinate in >= the specified value.
     */
    int findLowerBound(int x, int y, double z) const {
        const vector<pair<float, VoxelItem> >& bin = bins[x][y];
        int lower = 0;
        int upper = bin.size();
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first < z)
                lower = middle+1;
            else
                upper = middle;
        }
        return lower;
    }
    
    /**
     * Find the index of the first particle in voxel (x,y) whose z coordinate in greater than the specified value.
     */
    int findUpperBound(int x, int y, double z) const {
        const vector<pair<float, VoxelItem> >& bin = bins[x][y];
        int lower = 0;
        int upper = bin.size();
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first > z)
                upper = middle;
            else
                lower = middle+1;
        }
        return upper;
    }

    /**
     * Get the voxel index containing a particular location.
     */
    VoxelIndex getVoxelIndex(const float* location) const {
        float xperiodic, yperiodic;
        if (!usePeriodic) {
            xperiodic = location[0]-minx;
            yperiodic = location[1]-miny;
        }
        else {
            xperiodic = location[0]-periodicBoxSize[0]*floorf(location[0]/periodicBoxSize[0]);
            yperiodic = location[1]-periodicBoxSize[1]*floorf(location[1]/periodicBoxSize[1]);
        }
        int x = min(nx-1, int(floorf(xperiodic / voxelSizeX)));
        int y = min(ny-1, int(floorf(yperiodic / voxelSizeY)));
        
        return VoxelIndex(x, y);
    }

    void getNeighbors(vector<int>& neighbors, int blockIndex, fvec4 blockCenter, fvec4 blockWidth, const vector<int>& sortedAtoms, vector<char>& exclusions, float maxDistance, const vector<int> blockAtoms, const float* atomLocations) const {
        neighbors.resize(0);
        exclusions.resize(0);
        fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
        fvec4 invBoxSize(1/periodicBoxSize[0], 1/periodicBoxSize[1], 1/periodicBoxSize[2], 0);
        
        float maxDistanceSquared = maxDistance * maxDistance;
        float refineCutoff = maxDistance-max(max(blockWidth[0], blockWidth[1]), blockWidth[2]);
        float refineCutoffSquared = refineCutoff*refineCutoff;

        int dIndexX = min(nx/2, int((maxDistance+blockWidth[0])/voxelSizeX)+1); // How may voxels away do we have to look?
        int dIndexY = min(ny/2, int((maxDistance+blockWidth[1])/voxelSizeY)+1);
        float centerPos[4];
        blockCenter.store(centerPos);
        VoxelIndex centerVoxelIndex = getVoxelIndex(centerPos);
        int startx = centerVoxelIndex.x-dIndexX;
        int starty = centerVoxelIndex.y-dIndexY;
        int endx = centerVoxelIndex.x+dIndexX;
        int endy = centerVoxelIndex.y+dIndexY;
        int numRanges;
        if (usePeriodic) {
            endx = min(endx, centerVoxelIndex.x-dIndexX+nx-1);
            endy = min(endy, centerVoxelIndex.y-dIndexY+ny-1);
            numRanges = 2;
        }
        else {
            startx = max(startx, 0);
            starty = max(starty, 0);
            endx = min(endx, nx-1);
            endy = min(endy, ny-1);
            numRanges = 1;
        }
        int lastSortedIndex = BlockSize*(blockIndex+1);
        VoxelIndex voxelIndex(0, 0);
        for (int x = startx; x <= endx; ++x) {
            voxelIndex.x = x;
            if (usePeriodic)
                voxelIndex.x = (x < 0 ? x+nx : (x >= nx ? x-nx : x));
            float dx = max(0.0f, voxelSizeX*max(0, abs(centerVoxelIndex.x-x)-1)-blockWidth[0]);
            for (int y = starty; y <= endy; ++y) {
                voxelIndex.y = y;
                if (usePeriodic)
                    voxelIndex.y = (y < 0 ? y+ny : (y >= ny ? y-ny : y));
                float dy = max(0.0f, voxelSizeY*max(0, abs(centerVoxelIndex.y-y)-1)-blockWidth[1]);
                
                // Identify the range of atoms within this bin we need to search.  When using periodic boundary
                // conditions, there may be two separate ranges.
                
                float dz = maxDistance+blockWidth[2];
                dz = sqrtf(max(0.0f, dz*dz-dx*dx-dy*dy));
                int rangeStart[2];
                int rangeEnd[2];
                rangeStart[0] = findLowerBound(voxelIndex.x, voxelIndex.y, centerPos[2]-dz);
                if (usePeriodic) {
                    rangeEnd[0] = findUpperBound(voxelIndex.x, voxelIndex.y, centerPos[2]+dz);
                    if (rangeStart[0] > 0) {
                        rangeStart[1] = 0;
                        rangeEnd[1] = min(findUpperBound(voxelIndex.x, voxelIndex.y, centerPos[2]+dz-periodicBoxSize[2]), rangeStart[0]);
                    }
                    else {
                        rangeStart[1] = max(findLowerBound(voxelIndex.x, voxelIndex.y, centerPos[2]-dz+periodicBoxSize[2]), rangeEnd[0]);
                        rangeEnd[1] = bins[voxelIndex.x][voxelIndex.y].size();
                    }
                }
                else
                    rangeEnd[0] = findUpperBound(voxelIndex.x, voxelIndex.y, centerPos[2]+dz);
                
                // Loop over atoms and check to see if they are neighbors of this block.
                
                for (int range = 0; range < numRanges; range++) {
                    for (int item = rangeStart[range]; item < rangeEnd[range]; item++) {
                        const int sortedIndex = bins[voxelIndex.x][voxelIndex.y][item].second.second;

                        // Avoid duplicate entries.
                        if (sortedIndex >= lastSortedIndex)
                            continue;
                        
                        fvec4 atomPos(bins[voxelIndex.x][voxelIndex.y][item].second.first);
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
    float voxelSizeX, voxelSizeY;
    float minx, maxx, miny, maxy;
    int nx, ny;
    const float* periodicBoxSize;
    const bool usePeriodic;
    vector<vector<vector<pair<float, VoxelItem> > > > bins;
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
    
    float edgeSizeX, edgeSizeY;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = 0.5f*periodicBoxSize[0]/floorf(periodicBoxSize[0]/maxDistance);
        edgeSizeY = 0.5f*periodicBoxSize[1]/floorf(periodicBoxSize[1]/maxDistance);
    }
    Voxels voxels(edgeSizeX, edgeSizeY, minx, maxx, miny, maxy, periodicBoxSize, usePeriodic);
    for (int i = 0; i < numAtoms; i++) {
        int atomIndex = atomBins[i].second;
        sortedAtoms[i] = atomIndex;
        voxels.insert(i, &atomLocations[4*atomIndex]);
    }
    voxels.sortItems();
    
    // Record the parameters for the threads.
    
    this->voxels = &voxels;
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
            voxels->getNeighbors(blockNeighbors[i], i, (maxPos+minPos)*0.5f, (maxPos-minPos)*0.5f, sortedAtoms, blockExclusions[i], maxDistance, blockAtoms, atomLocations);
            
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
