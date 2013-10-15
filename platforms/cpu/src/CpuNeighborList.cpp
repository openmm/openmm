#include "CpuNeighborList.h"
#include "openmm/internal/hardware.h"
#include <set>
#include <map>
#include <cmath>
#include <smmintrin.h>

using namespace std;

namespace OpenMM {

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
            voxelSizeY = periodicBoxSize[0]/ny;
            voxelSizeZ = periodicBoxSize[0]/nz;
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

    void getNeighbors(vector<pair<int, int> >& neighbors, const VoxelItem& referencePoint, const vector<set<int> >& exclusions, float maxDistance) const {

        // Loop over neighboring voxels
        // TODO use more clever selection of neighboring voxels

        const int atomI = referencePoint.second;
        const float* locationI = referencePoint.first;
        __m128 posI = _mm_loadu_ps(locationI);
        __m128 boxSize = _mm_set_ps(0, periodicBoxSize[2], periodicBoxSize[1], periodicBoxSize[0]);
        __m128 invBoxSize = _mm_set_ps(0, (1/periodicBoxSize[2]), (1/periodicBoxSize[1]), (1/periodicBoxSize[0]));
        __m128 half = _mm_set1_ps(0.5);
        
        float maxDistanceSquared = maxDistance * maxDistance;

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
                        const int atomJ = itemIter->second;

                        // Avoid duplicate entries.
                        if (atomJ >= atomI)
                            break;
                        
                        __m128 posJ = _mm_loadu_ps(itemIter->first);
                        __m128 delta = _mm_sub_ps(posJ, posI);
                        if (usePeriodic) {
                            __m128 base = _mm_mul_ps(_mm_floor_ps(_mm_add_ps(_mm_mul_ps(delta, invBoxSize), half)), boxSize);
                            delta = _mm_sub_ps(delta, base);
                        }
                        float dSquared = _mm_cvtss_f32(_mm_dp_ps(delta, delta, 0x71));
                        if (dSquared > maxDistanceSquared)
                            continue;
                        
                        // Ignore exclusions.
                        if (exclusions[atomI].find(atomJ) != exclusions[atomI].end())
                            continue;

                        neighbors.push_back(make_pair(atomI, atomJ));
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
    vector<pair<int, int> > threadNeighbors;
};

static void* threadBody(void* args) {
    CpuNeighborList::ThreadData& data = *reinterpret_cast<CpuNeighborList::ThreadData*>(args);
    data.owner.runThread(data.index, data.threadNeighbors);
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
    for (int i = 0; i < numAtoms; i++)
        voxelHash.insert(i, &atomLocations[4*i]);
    
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
    
    // Combine the results from all the threads.
    
    neighbors.clear();
    for (int i = 0; i < numThreads; i++)
        neighbors.insert(neighbors.end(), threadData[i]->threadNeighbors.begin(), threadData[i]->threadNeighbors.end());
}

const vector<pair<int, int> >& CpuNeighborList::getNeighbors() {
    return neighbors;
}

void CpuNeighborList::runThread(int index, vector<pair<int, int> >& threadNeighbors) {
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
        
        threadNeighbors.clear();
        for (int i = index; i < numAtoms; i += numThreads)
            voxelHash->getNeighbors(threadNeighbors, VoxelItem(&atomLocations[4*i], i), *exclusions, maxDistance);
    }
}

} // namespace OpenMM
