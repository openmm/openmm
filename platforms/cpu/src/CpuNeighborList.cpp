/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CpuNeighborList.h"
#include "openmm/internal/hardware.h"
#include "openmm/internal/vectorize.h"
#include "hilbert.h"
#include <algorithm>
#include <set>
#include <map>
#include <cmath>

using namespace std;

namespace OpenMM {

class VoxelIndex 
{
public:
    VoxelIndex() : x(0), y(0) {
    }
    VoxelIndex(int x, int y) : x(x), y(y) {
    }
    int x;
    int y;
};

/**
 * This data structure organizes the particles spatially.  It divides them into bins along the x and y axes,
 * then sorts each bin along the z axis so ranges can be identified quickly with a binary search.
 */
class CpuNeighborList::Voxels {
public:
    Voxels(int blockSize, float vsx, float vsy, float minx, float maxx, float miny, float maxy, const float* periodicBoxSize, bool usePeriodic) :
            blockSize(blockSize), voxelSizeX(vsx), voxelSizeY(vsy), minx(minx), maxx(maxx), miny(miny), maxy(maxy), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
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
        bins[voxelIndex.x][voxelIndex.y].push_back(make_pair(location[2], atom));
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
        const vector<pair<float, int> >& bin = bins[x][y];
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
        const vector<pair<float, int> >& bin = bins[x][y];
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

    void getNeighbors(vector<int>& neighbors, int blockIndex, const fvec4& blockCenter, const fvec4& blockWidth, const vector<int>& sortedAtoms, vector<char>& exclusions, float maxDistance, const vector<int>& blockAtoms, const float* atomLocations, const vector<VoxelIndex>& atomVoxelIndex) const {
        neighbors.resize(0);
        exclusions.resize(0);
        fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
        fvec4 invBoxSize(1/periodicBoxSize[0], 1/periodicBoxSize[1], 1/periodicBoxSize[2], 0);
        
        float maxDistanceSquared = maxDistance * maxDistance;
        float refineCutoff = maxDistance-max(max(blockWidth[0], blockWidth[1]), blockWidth[2]);
        float refineCutoffSquared = refineCutoff*refineCutoff;

        int dIndexX = int((maxDistance+blockWidth[0])/voxelSizeX)+1; // How may voxels away do we have to look?
        int dIndexY = int((maxDistance+blockWidth[1])/voxelSizeY)+1;
        if (usePeriodic) {
            dIndexX = min(nx/2, dIndexX);
            dIndexY = min(ny/2, dIndexY);
        }
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
        }
        else {
            startx = max(startx, 0);
            starty = max(starty, 0);
            endx = min(endx, nx-1);
            endy = min(endy, ny-1);
        }
        int lastSortedIndex = blockSize*(blockIndex+1);
        VoxelIndex voxelIndex(0, 0);
        for (int x = startx; x <= endx; ++x) {
            voxelIndex.x = x;
            if (usePeriodic)
                voxelIndex.x = (x < 0 ? x+nx : (x >= nx ? x-nx : x));
            for (int y = starty; y <= endy; ++y) {
                voxelIndex.y = y;
                if (usePeriodic)
                    voxelIndex.y = (y < 0 ? y+ny : (y >= ny ? y-ny : y));
                
                // Identify the range of atoms within this bin we need to search.  When using periodic boundary
                // conditions, there may be two separate ranges.
                
                float minz = centerPos[2];
                float maxz = centerPos[2];
                fvec4 offset(voxelSizeX*x+(usePeriodic ? 0.0f : minx), voxelSizeY*y+(usePeriodic ? 0.0f : miny), 0, 0);
                for (int k = 0; k < (int) blockAtoms.size(); k++) {
                    const float* atomPos = &atomLocations[4*blockAtoms[k]];
                    fvec4 posVec(atomPos);
                    fvec4 delta1 = offset-posVec;
                    fvec4 delta2 = delta1+fvec4(voxelSizeX, voxelSizeY, 0, 0);
                    if (usePeriodic) {
                        delta1 -= round(delta1*invBoxSize)*boxSize;
                        delta2 -= round(delta2*invBoxSize)*boxSize;
                    }
                    fvec4 delta = min(abs(delta1), abs(delta2));
                    float dx = (x == atomVoxelIndex[k].x ? 0.0f : delta[0]);
                    float dy = (y == atomVoxelIndex[k].y ? 0.0f : delta[1]);
                    float dist2 = maxDistanceSquared-dx*dx-dy*dy;
                    if (dist2 > 0) {
                        float dist = sqrtf(dist2);
                        minz = min(minz, atomPos[2]-dist);
                        maxz = max(maxz, atomPos[2]+dist);
                    }
                }
                if (minz == maxz)
                    continue;
                bool needPeriodic = (centerPos[0]-blockWidth[0] < maxDistance || centerPos[0]+blockWidth[0] > periodicBoxSize[0]-maxDistance ||
                                     centerPos[1]-blockWidth[1] < maxDistance || centerPos[1]+blockWidth[1] > periodicBoxSize[1]-maxDistance ||
                                     minz < 0.0f || maxz > periodicBoxSize[2]);
                int rangeStart[2];
                int rangeEnd[2];
                rangeStart[0] = findLowerBound(voxelIndex.x, voxelIndex.y, minz);
                if (needPeriodic) {
                    numRanges = 2;
                    rangeEnd[0] = findUpperBound(voxelIndex.x, voxelIndex.y, maxz);
                    if (rangeStart[0] > 0) {
                        rangeStart[1] = 0;
                        rangeEnd[1] = min(findUpperBound(voxelIndex.x, voxelIndex.y, maxz-periodicBoxSize[2]), rangeStart[0]);
                    }
                    else {
                        rangeStart[1] = max(findLowerBound(voxelIndex.x, voxelIndex.y, minz+periodicBoxSize[2]), rangeEnd[0]);
                        rangeEnd[1] = bins[voxelIndex.x][voxelIndex.y].size();
                    }
                }
                else {
                    numRanges = 1;
                    rangeEnd[0] = findUpperBound(voxelIndex.x, voxelIndex.y, maxz);
                }
                
                // Loop over atoms and check to see if they are neighbors of this block.
                
                for (int range = 0; range < numRanges; range++) {
                    for (int item = rangeStart[range]; item < rangeEnd[range]; item++) {
                        const int sortedIndex = bins[voxelIndex.x][voxelIndex.y][item].second;

                        // Avoid duplicate entries.
                        if (sortedIndex >= lastSortedIndex)
                            continue;
                        
                        fvec4 atomPos(atomLocations+4*sortedAtoms[sortedIndex]);
                        fvec4 delta = atomPos-centerPos;
                        if (needPeriodic) {
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
                                if (needPeriodic) {
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
                        if (sortedIndex < blockSize*blockIndex)
                            exclusions.push_back(0);
                        else {
                            int mask = (1<<blockSize)-1;
                            exclusions.push_back(mask & (mask<<(sortedIndex-blockSize*blockIndex)));
                        }
                    }
                }
            }
        }
    }

private:
    int blockSize;
    float voxelSizeX, voxelSizeY;
    float minx, maxx, miny, maxy;
    int nx, ny;
    const float* periodicBoxSize;
    const bool usePeriodic;
    vector<vector<vector<pair<float, int> > > > bins;
};

class CpuNeighborList::ThreadTask : public ThreadPool::Task {
public:
    ThreadTask(CpuNeighborList& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadComputeNeighborList(threads, threadIndex);
    }
    CpuNeighborList& owner;
};

CpuNeighborList::CpuNeighborList(int blockSize) : blockSize(blockSize) {
}

void CpuNeighborList::computeNeighborList(int numAtoms, const AlignedArray<float>& atomLocations, const vector<set<int> >& exclusions,
            const float* periodicBoxSize, bool usePeriodic, float maxDistance, ThreadPool& threads) {
    int numBlocks = (numAtoms+blockSize-1)/blockSize;
    blockNeighbors.resize(numBlocks);
    blockExclusions.resize(numBlocks);
    sortedAtoms.resize(numAtoms);
    
    // Record the parameters for the threads.
    
    this->exclusions = &exclusions;
    this->atomLocations = &atomLocations[0];
    this->periodicBoxSize = periodicBoxSize;
    this->numAtoms = numAtoms;
    this->usePeriodic = usePeriodic;
    this->maxDistance = maxDistance;
    
    // Identify the range of atom positions along each axis.
    
    fvec4 minPos(&atomLocations[0]);
    fvec4 maxPos = minPos;
    for (int i = 0; i < numAtoms; i++) {
        fvec4 pos(&atomLocations[4*i]);
        minPos = min(minPos, pos);
        maxPos = max(maxPos, pos);
    }
    minx = minPos[0];
    maxx = maxPos[0];
    miny = minPos[1];
    maxy = maxPos[1];
    minz = minPos[2];
    maxz = maxPos[2];
    
    // Sort the atoms based on a Hilbert curve.
    
    atomBins.resize(numAtoms);
    ThreadTask task(*this);
    threads.execute(task);
    threads.waitForThreads();
    sort(atomBins.begin(), atomBins.end());

    // Build the voxel hash.
    
    float edgeSizeX, edgeSizeY;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = 0.6f*periodicBoxSize[0]/floorf(periodicBoxSize[0]/maxDistance);
        edgeSizeY = 0.6f*periodicBoxSize[1]/floorf(periodicBoxSize[1]/maxDistance);
    }
    Voxels voxels(blockSize, edgeSizeX, edgeSizeY, minx, maxx, miny, maxy, periodicBoxSize, usePeriodic);
    for (int i = 0; i < numAtoms; i++) {
        int atomIndex = atomBins[i].second;
        sortedAtoms[i] = atomIndex;
        voxels.insert(i, &atomLocations[4*atomIndex]);
    }
    voxels.sortItems();
    this->voxels = &voxels;
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.resumeThreads();
    threads.waitForThreads();
    
    // Add padding atoms to fill up the last block.
    
    int numPadding = numBlocks*blockSize-numAtoms;
    if (numPadding > 0) {
        char mask = ((0xFFFF-(1<<blockSize)+1) >> numPadding);
        for (int i = 0; i < numPadding; i++)
            sortedAtoms.push_back(0);
        vector<char>& exc = blockExclusions[blockExclusions.size()-1];
        for (int i = 0; i < (int) exc.size(); i++)
            exc[i] |= mask;
    }
}

int CpuNeighborList::getNumBlocks() const {
    return sortedAtoms.size()/blockSize;
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

void CpuNeighborList::threadComputeNeighborList(ThreadPool& threads, int threadIndex) {
    // Compute the positions of atoms along the Hilbert curve.

    float binWidth = max(max(maxx-minx, maxy-miny), maxz-minz)/255.0f;
    float invBinWidth = 1.0f/binWidth;
    bitmask_t coords[3];
    int numThreads = threads.getNumThreads();
    for (int i = threadIndex; i < numAtoms; i += numThreads) {
        const float* pos = &atomLocations[4*i];
        coords[0] = (bitmask_t) ((pos[0]-minx)*invBinWidth);
        coords[1] = (bitmask_t) ((pos[1]-miny)*invBinWidth);
        coords[2] = (bitmask_t) ((pos[2]-minz)*invBinWidth);
        int bin = (int) hilbert_c2i(3, 8, coords);
        atomBins[i] = pair<int, int>(bin, i);
    }
    threads.syncThreads();

    // Compute this thread's subset of neighbors.

    int numBlocks = blockNeighbors.size();
    vector<int> blockAtoms;
    vector<VoxelIndex> atomVoxelIndex;
    for (int i = threadIndex; i < numBlocks; i += numThreads) {
        // Find the atoms in this block and compute their bounding box.
        
        int firstIndex = blockSize*i;
        int atomsInBlock = min(blockSize, numAtoms-firstIndex);
        blockAtoms.resize(atomsInBlock);
        atomVoxelIndex.resize(atomsInBlock);
        for (int j = 0; j < atomsInBlock; j++) {
            blockAtoms[j] = sortedAtoms[firstIndex+j];
            atomVoxelIndex[j] = voxels->getVoxelIndex(&atomLocations[4*blockAtoms[j]]);
        }
        fvec4 minPos(&atomLocations[4*sortedAtoms[firstIndex]]);
        fvec4 maxPos = minPos;
        for (int j = 1; j < atomsInBlock; j++) {
            fvec4 pos(&atomLocations[4*sortedAtoms[firstIndex+j]]);
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        voxels->getNeighbors(blockNeighbors[i], i, (maxPos+minPos)*0.5f, (maxPos-minPos)*0.5f, sortedAtoms, blockExclusions[i], maxDistance, blockAtoms, atomLocations, atomVoxelIndex);

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

} // namespace OpenMM
