/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2018 Stanford University and the Authors.      *
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
    VoxelIndex() : y(0), z(0) {
    }
    VoxelIndex(int y, int z) : y(y), z(z) {
    }
    int y;
    int z;
};

/**
 * This data structure organizes the particles spatially.  It divides them into bins along the x and y axes,
 * then sorts each bin along the z axis so ranges can be identified quickly with a binary search.
 */
class CpuNeighborList::Voxels {
public:
    Voxels(int blockSize, float vsy, float vsz, float miny, float maxy, float minz, float maxz, const Vec3* boxVectors, bool usePeriodic) :
            blockSize(blockSize), voxelSizeY(vsy), voxelSizeZ(vsz), miny(miny), maxy(maxy), minz(minz), maxz(maxz), usePeriodic(usePeriodic) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                periodicBoxVectors[i][j] = (float) boxVectors[i][j];
        periodicBoxSize[0] = (float) boxVectors[0][0];
        periodicBoxSize[1] = (float) boxVectors[1][1];
        periodicBoxSize[2] = (float) boxVectors[2][2];
        recipBoxSize[0] = (float) (1/boxVectors[0][0]);
        recipBoxSize[1] = (float) (1/boxVectors[1][1]);
        recipBoxSize[2] = (float) (1/boxVectors[2][2]);
        triclinic = (boxVectors[0][1] != 0.0 || boxVectors[0][2] != 0.0 ||
                     boxVectors[1][0] != 0.0 || boxVectors[1][2] != 0.0 ||
                     boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0.0);
        if (usePeriodic) {
            ny = (int) floorf(boxVectors[1][1]/voxelSizeY+0.5f);
            nz = (int) floorf(boxVectors[2][2]/voxelSizeZ+0.5f);
            voxelSizeY = boxVectors[1][1]/ny;
            voxelSizeZ = boxVectors[2][2]/nz;
        }
        else {
            ny = max(1, min(500, (int) floorf((maxy-miny)/voxelSizeY+0.5f)));
            nz = max(1, min(500, (int) floorf((maxz-minz)/voxelSizeZ+0.5f)));
            if (maxy > miny)
                voxelSizeY = (maxy-miny)/ny;
            if (maxz > minz)
                voxelSizeZ = (maxz-minz)/nz;
        }
        bins.resize(ny);
        for (int i = 0; i < ny; i++) {
            bins[i].resize(nz);
            for (int j = 0; j < nz; j++)
                bins[i][j].resize(0);
        }
    }

    /**
     * Insert a particle into the voxel data structure.
     */
    void insert(const int& atom, const float* location) {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        bins[voxelIndex.y][voxelIndex.z].push_back(make_pair(location[0], atom));
    }
    
    /**
     * Sort the particles in each voxel by x coordinate.
     */
    void sortItems() {
        for (int i = 0; i < ny; i++)
            for (int j = 0; j < nz; j++)
                sort(bins[i][j].begin(), bins[i][j].end());
    }
    
    /**
     * Find the index of the first particle in voxel (y,z) whose x coordinate is >= the specified value.
     */
    int findLowerBound(int y, int z, double x, int lower, int upper) const {
        const vector<pair<float, int> >& bin = bins[y][z];
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first < x)
                lower = middle+1;
            else
                upper = middle;
        }
        return lower;
    }
    
    /**
     * Find the index of the first particle in voxel (y,z) whose x coordinate is greater than the specified value.
     */
    int findUpperBound(int y, int z, double x, int lower, int upper) const {
        const vector<pair<float, int> >& bin = bins[y][z];
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first > x)
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
        float yperiodic, zperiodic;
        if (!usePeriodic) {
            yperiodic = location[1]-miny;
            zperiodic = location[2]-minz;
        }
        else {
            float scale2 = floorf(location[2]*recipBoxSize[2]);
            yperiodic = location[1]-periodicBoxVectors[2][1]*scale2;
            zperiodic = location[2]-periodicBoxVectors[2][2]*scale2;
            float scale1 = floorf(yperiodic*recipBoxSize[1]);
            yperiodic -= periodicBoxVectors[1][0]*scale1;
        }
        int y = max(0, min(ny-1, int(floorf(yperiodic / voxelSizeY))));
        int z = max(0, min(nz-1, int(floorf(zperiodic / voxelSizeZ))));
        
        return VoxelIndex(y, z);
    }
        
    void getNeighbors(vector<int>& neighbors, int blockIndex, const fvec4& blockCenter, const fvec4& blockWidth, const vector<int>& sortedAtoms, vector<char>& exclusions, float maxDistance, const vector<int>& blockAtoms, const vector<float>& blockAtomX, const vector<float>& blockAtomY, const vector<float>& blockAtomZ, const vector<float>& sortedPositions, const vector<VoxelIndex>& atomVoxelIndex) const {
        neighbors.resize(0);
        exclusions.resize(0);
        fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
        fvec4 invBoxSize(recipBoxSize[0], recipBoxSize[1], recipBoxSize[2], 0);
        fvec4 periodicBoxVec4[3];
        periodicBoxVec4[0] = fvec4(periodicBoxVectors[0][0], periodicBoxVectors[0][1], periodicBoxVectors[0][2], 0);
        periodicBoxVec4[1] = fvec4(periodicBoxVectors[1][0], periodicBoxVectors[1][1], periodicBoxVectors[1][2], 0);
        periodicBoxVec4[2] = fvec4(periodicBoxVectors[2][0], periodicBoxVectors[2][1], periodicBoxVectors[2][2], 0);

        float maxDistanceSquared = maxDistance * maxDistance;
        float refineCutoff = maxDistance-max(max(blockWidth[0], blockWidth[1]), blockWidth[2]);
        float refineCutoffSquared = refineCutoff*refineCutoff;

        int dIndexY = int((maxDistance+blockWidth[1])/voxelSizeY)+1; // How may voxels away do we have to look?
        int dIndexZ = int((maxDistance+blockWidth[2])/voxelSizeZ)+1;
        if (usePeriodic) {
            dIndexY = min(ny/2, dIndexY);
            dIndexZ = min(nz/2, dIndexZ);
        }
        float centerPos[4];
        blockCenter.store(centerPos);
        VoxelIndex centerVoxelIndex = getVoxelIndex(centerPos);

        // Loop over voxels along the z axis.

        int startz = centerVoxelIndex.z-dIndexZ;
        int endz = centerVoxelIndex.z+dIndexZ;
        if (usePeriodic)
            endz = min(endz, startz+nz-1);
        else {
            startz = max(startz, 0);
            endz = min(endz, nz-1);
        }
        int lastSortedIndex = blockSize*(blockIndex+1);
        VoxelIndex voxelIndex(0, 0);
        for (int z = startz; z <= endz; ++z) {
            voxelIndex.z = z;
            if (usePeriodic)
                voxelIndex.z = (z < 0 ? z+nz : (z >= nz ? z-nz : z));

            // Loop over voxels along the y axis.

            float boxz = floor((float) z/nz);
            int starty = centerVoxelIndex.y-dIndexY;
            int endy = centerVoxelIndex.y+dIndexY;
            float yoffset = (float) (usePeriodic ? boxz*periodicBoxVectors[2][1] : 0);
            if (usePeriodic) {
                starty -= (int) ceil(yoffset/voxelSizeY);
                endy -= (int) floor(yoffset/voxelSizeY);
                endy = min(endy, starty+ny-1);
            }
            else {
                starty = max(starty, 0);
                endy = min(endy, ny-1);
            }
            for (int y = starty; y <= endy; ++y) {
                voxelIndex.y = y;
                if (usePeriodic)
                    voxelIndex.y = (y < 0 ? y+ny : (y >= ny ? y-ny : y));
                float boxy = floor((float) y/ny);
                
                // Identify the range of atoms within this bin we need to search.  When using periodic boundary
                // conditions, there may be two separate ranges.
                
                float minx = centerPos[0];
                float maxx = centerPos[0];
                if (usePeriodic && triclinic) {
                    for (int k = 0; k < (int) blockAtoms.size(); k++) {
                        const float* atomPos = &sortedPositions[4*(blockSize*blockIndex+k)];
                        fvec4 delta1(0, voxelSizeY*voxelIndex.y-atomPos[1], voxelSizeZ*voxelIndex.z-atomPos[2], 0);
                        fvec4 delta2 = delta1+fvec4(0, voxelSizeY, 0, 0);
                        fvec4 delta3 = delta1+fvec4(0, 0, voxelSizeZ, 0);
                        fvec4 delta4 = delta1+fvec4(0, voxelSizeY, voxelSizeZ, 0);
                        delta1 -= periodicBoxVec4[2]*floorf(delta1[2]*recipBoxSize[2]+0.5f);
                        delta1 -= periodicBoxVec4[1]*floorf(delta1[1]*recipBoxSize[1]+0.5f);
                        delta1 -= periodicBoxVec4[0]*floorf(delta1[0]*recipBoxSize[0]+0.5f);
                        delta2 -= periodicBoxVec4[2]*floorf(delta2[2]*recipBoxSize[2]+0.5f);
                        delta2 -= periodicBoxVec4[1]*floorf(delta2[1]*recipBoxSize[1]+0.5f);
                        delta2 -= periodicBoxVec4[0]*floorf(delta2[0]*recipBoxSize[0]+0.5f);
                        delta3 -= periodicBoxVec4[2]*floorf(delta3[2]*recipBoxSize[2]+0.5f);
                        delta3 -= periodicBoxVec4[1]*floorf(delta3[1]*recipBoxSize[1]+0.5f);
                        delta3 -= periodicBoxVec4[0]*floorf(delta3[0]*recipBoxSize[0]+0.5f);
                        delta4 -= periodicBoxVec4[2]*floorf(delta4[2]*recipBoxSize[2]+0.5f);
                        delta4 -= periodicBoxVec4[1]*floorf(delta4[1]*recipBoxSize[1]+0.5f);
                        delta4 -= periodicBoxVec4[0]*floorf(delta4[0]*recipBoxSize[0]+0.5f);
                        if (delta1[1] < 0 && delta1[1]+voxelSizeY > 0)
                            delta1 = fvec4(delta1[0], 0, delta1[2], 0);
                        if (delta1[2] < 0 && delta1[2]+voxelSizeZ > 0)
                            delta1 = fvec4(delta1[0], delta1[1], 0, 0);
                        if (delta3[1] < 0 && delta3[1]+voxelSizeY > 0)
                            delta3 = fvec4(delta3[0], 0, delta3[2], 0);
                        if (delta2[2] < 0 && delta2[2]+voxelSizeZ > 0)
                            delta2 = fvec4(delta2[0], delta2[1], 0, 0);
                        fvec4 delta = min(min(min(abs(delta1), abs(delta2)), abs(delta3)), abs(delta4));
                        float dy = (voxelIndex.y == atomVoxelIndex[k].y ? 0.0f : delta[1]);
                        float dz = (voxelIndex.z == atomVoxelIndex[k].z ? 0.0f : delta[2]);
                        float dist2 = maxDistanceSquared-dy*dy-dz*dz;
                        if (dist2 > 0) {
                            float dist = sqrtf(dist2);
                            minx = min(minx, atomPos[0]-dist-max(max(max(delta1[0], delta2[0]), delta3[0]), delta4[0]));
                            maxx = max(maxx, atomPos[0]+dist-min(min(min(delta1[0], delta2[0]), delta3[0]), delta4[0]));
                        }
                    }
                }
                else {
                    float xoffset = (float) (usePeriodic ? boxy*periodicBoxVectors[1][0]+boxz*periodicBoxVectors[2][0] : 0);
                    fvec4 offset(-xoffset, -yoffset+voxelSizeY*y+(usePeriodic ? 0.0f : miny), voxelSizeZ*z+(usePeriodic ? 0.0f : minz), 0);
                    for (int k = 0; k < (int) blockAtoms.size(); k++) {
                        const float* atomPos = &sortedPositions[4*(blockSize*blockIndex+k)];
                        fvec4 posVec(atomPos);
                        fvec4 delta1 = offset-posVec;
                        fvec4 delta2 = delta1+fvec4(0, voxelSizeY, voxelSizeZ, 0);
                        if (usePeriodic) {
                            delta1 -= round(delta1*invBoxSize)*boxSize;
                            delta2 -= round(delta2*invBoxSize)*boxSize;
                        }
                        fvec4 delta = min(abs(delta1), abs(delta2));
                        float dy = (y == atomVoxelIndex[k].y ? 0.0f : delta[1]);
                        float dz = (z == atomVoxelIndex[k].z ? 0.0f : delta[2]);
                        float dist2 = maxDistanceSquared-dy*dy-dz*dz;
                        if (dist2 > 0) {
                            float dist = sqrtf(dist2);
                            minx = min(minx, atomPos[0]-dist-xoffset);
                            maxx = max(maxx, atomPos[0]+dist-xoffset);
                        }
                    }
                }
                if (minx == maxx)
                    continue;
                bool needPeriodic = usePeriodic && (centerPos[1]-blockWidth[1] < maxDistance || centerPos[1]+blockWidth[1] > periodicBoxSize[1]-maxDistance ||
                                                    centerPos[2]-blockWidth[2] < maxDistance || centerPos[2]+blockWidth[2] > periodicBoxSize[2]-maxDistance ||
                                                    minx < 0.0f || maxx > periodicBoxVectors[0][0]);
                int numRanges;
                int rangeStart[2];
                int rangeEnd[2];
                int binSize = bins[voxelIndex.y][voxelIndex.z].size();
                rangeStart[0] = findLowerBound(voxelIndex.y, voxelIndex.z, minx, 0, binSize);
                if (needPeriodic) {
                    numRanges = 2;
                    rangeEnd[0] = findUpperBound(voxelIndex.y, voxelIndex.z, maxx, rangeStart[0], binSize);
                    if (rangeStart[0] > 0 && rangeEnd[0] < binSize)
                        numRanges = 1;
                    else if (rangeStart[0] > 0) {
                        rangeStart[1] = 0;
                        rangeEnd[1] = min(findUpperBound(voxelIndex.y, voxelIndex.z, maxx-periodicBoxSize[0], 0, rangeStart[0]), rangeStart[0]);
                    }
                    else {
                        rangeStart[1] = max(findLowerBound(voxelIndex.y, voxelIndex.z, minx+periodicBoxSize[0], rangeEnd[0], binSize), rangeEnd[0]);
                        rangeEnd[1] = bins[voxelIndex.y][voxelIndex.z].size();
                    }
                }
                else {
                    numRanges = 1;
                    rangeEnd[0] = findUpperBound(voxelIndex.y, voxelIndex.z, maxx, rangeStart[0], binSize);
                }
                bool periodicRectangular = (needPeriodic && !triclinic);
                
                // Loop over atoms and check to see if they are neighbors of this block.
                
                const vector<pair<float, int> >& voxelBins = bins[voxelIndex.y][voxelIndex.z];
                for (int range = 0; range < numRanges; range++) {
                    for (int item = rangeStart[range]; item < rangeEnd[range]; item++) {
                        const int sortedIndex = voxelBins[item].second;

                        // Avoid duplicate entries.
                        if (sortedIndex >= lastSortedIndex)
                            continue;
                        
                        fvec4 atomPos(&sortedPositions[4*sortedIndex]);
                        fvec4 delta = atomPos-blockCenter;
                        if (periodicRectangular)
                            delta -= round(delta*invBoxSize)*boxSize;
                        else if (needPeriodic) {
                            delta -= periodicBoxVec4[2]*floorf(delta[2]*recipBoxSize[2]+0.5f);
                            delta -= periodicBoxVec4[1]*floorf(delta[1]*recipBoxSize[1]+0.5f);
                            delta -= periodicBoxVec4[0]*floorf(delta[0]*recipBoxSize[0]+0.5f);
                        }
                        delta = max(0.0f, abs(delta)-blockWidth);
                        float dSquared = dot3(delta, delta);
                        if (dSquared > maxDistanceSquared)
                            continue;
                        
                        if (dSquared > refineCutoffSquared) {
                            // The distance is large enough that there might not be any actual interactions.
                            // Check individual atom pairs to be sure.
                            
                            bool anyInteraction = false;
                            for (int k = 0; k < (int) blockAtoms.size(); k += 4) {
                                fvec4 dx = fvec4(&blockAtomX[k])-atomPos[0];
                                fvec4 dy = fvec4(&blockAtomY[k])-atomPos[1];
                                fvec4 dz = fvec4(&blockAtomZ[k])-atomPos[2];
                                if (periodicRectangular) {
                                    dx -= round(dx*invBoxSize[0])*boxSize[0];
                                    dy -= round(dy*invBoxSize[1])*boxSize[1];
                                    dz -= round(dz*invBoxSize[2])*boxSize[2];
                                }
                                else if (needPeriodic) {
                                    fvec4 scale3 = floor(dz*recipBoxSize[2]+0.5f);
                                    dx -= scale3*periodicBoxVectors[2][0];
                                    dy -= scale3*periodicBoxVectors[2][1];
                                    dz -= scale3*periodicBoxVectors[2][2];
                                    fvec4 scale2 = floor(dy*recipBoxSize[1]+0.5f);
                                    dx -= scale2*periodicBoxVectors[1][0];
                                    dy -= scale2*periodicBoxVectors[1][1];
                                    fvec4 scale1 = floor(dx*recipBoxSize[0]+0.5f);
                                    dx -= scale1*periodicBoxVectors[0][0];
                                }
                                fvec4 r2 = dx*dx + dy*dy + dz*dz;
                                if (any(r2 < maxDistanceSquared)) {
                                    anyInteraction = true;
                                    break;
                                }
                            }
                            if (!anyInteraction)
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
    float voxelSizeY, voxelSizeZ;
    float miny, maxy, minz, maxz;
    int ny, nz;
    float periodicBoxSize[3], recipBoxSize[3];
    bool triclinic;
    float periodicBoxVectors[3][3];
    const bool usePeriodic;
    vector<vector<vector<pair<float, int> > > > bins;
};

CpuNeighborList::CpuNeighborList(int blockSize) : blockSize(blockSize) {
}

void CpuNeighborList::computeNeighborList(int numAtoms, const AlignedArray<float>& atomLocations, const vector<set<int> >& exclusions,
            const Vec3* periodicBoxVectors, bool usePeriodic, float maxDistance, ThreadPool& threads) {
    int numBlocks = (numAtoms+blockSize-1)/blockSize;
    blockNeighbors.resize(numBlocks);
    blockExclusions.resize(numBlocks);
    sortedAtoms.resize(numAtoms);
    sortedPositions.resize(4*numAtoms);
    
    // Record the parameters for the threads.
    
    this->exclusions = &exclusions;
    this->atomLocations = &atomLocations[0];
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
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
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeNeighborList(threads, threadIndex); });
    threads.waitForThreads();
    sort(atomBins.begin(), atomBins.end());

    // Build the voxel hash.

    float edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeY = 0.6f*periodicBoxVectors[1][1]/floorf(periodicBoxVectors[1][1]/maxDistance);
        edgeSizeZ = 0.6f*periodicBoxVectors[2][2]/floorf(periodicBoxVectors[2][2]/maxDistance);
    }
    Voxels voxels(blockSize, edgeSizeY, edgeSizeZ, miny, maxy, minz, maxz, periodicBoxVectors, usePeriodic);
    for (int i = 0; i < numAtoms; i++) {
        int atomIndex = atomBins[i].second;
        sortedAtoms[i] = atomIndex;
        fvec4 atomPos(&atomLocations[4*atomIndex]);
        atomPos.store(&sortedPositions[4*i]);
        voxels.insert(i, &atomLocations[4*atomIndex]);
    }
    voxels.sortItems();
    this->voxels = &voxels;

    // Signal the threads to start running and wait for them to finish.
    
    atomicCounter = 0;
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

int CpuNeighborList::getBlockSize() const {
    return blockSize;
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
    vector<float> blockAtomX(blockSize), blockAtomY(blockSize), blockAtomZ(blockSize);
    vector<VoxelIndex> atomVoxelIndex;
    while (true) {
        int i = atomicCounter++;
        if (i >= numBlocks)
            break;

        // Find the atoms in this block and compute their bounding box.
        
        int firstIndex = blockSize*i;
        int atomsInBlock = min(blockSize, numAtoms-firstIndex);
        blockAtoms.resize(atomsInBlock);
        atomVoxelIndex.resize(atomsInBlock);
        for (int j = 0; j < atomsInBlock; j++) {
            blockAtoms[j] = sortedAtoms[firstIndex+j];
            atomVoxelIndex[j] = voxels->getVoxelIndex(&atomLocations[4*blockAtoms[j]]);
        }
        fvec4 minPos(&sortedPositions[4*firstIndex]);
        fvec4 maxPos = minPos;
        for (int j = 1; j < atomsInBlock; j++) {
            fvec4 pos(&sortedPositions[4*(firstIndex+j)]);
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        for (int j = 0; j < atomsInBlock; j++) {
            blockAtomX[j] = sortedPositions[4*(firstIndex+j)];
            blockAtomY[j] = sortedPositions[4*(firstIndex+j)+1];
            blockAtomZ[j] = sortedPositions[4*(firstIndex+j)+2];
        }
        for (int j = atomsInBlock; j < blockSize; j++) {
            blockAtomX[j] = 1e10;
            blockAtomY[j] = 1e10;
            blockAtomZ[j] = 1e10;
        }
        voxels->getNeighbors(blockNeighbors[i], i, (maxPos+minPos)*0.5f, (maxPos-minPos)*0.5f, sortedAtoms, blockExclusions[i], maxDistance, blockAtoms, blockAtomX, blockAtomY, blockAtomZ, sortedPositions, atomVoxelIndex);

        // Record the exclusions for this block.

        map<int, char> atomFlags;
        for (int j = 0; j < atomsInBlock; j++) {
            const set<int>& atomExclusions = (*exclusions)[sortedAtoms[firstIndex+j]];
            char mask = 1<<j;
            for (int exclusion : atomExclusions) {
                map<int, char>::iterator thisAtomFlags = atomFlags.find(exclusion);
                if (thisAtomFlags == atomFlags.end())
                    atomFlags[exclusion] = mask;
                else
                    thisAtomFlags->second |= mask;
            }
        }
        int numNeighbors = blockNeighbors[i].size();
        for (int k = 0; k < numNeighbors; k++) {
            int atomIndex = blockNeighbors[i][k];
            map<int, char>::iterator thisAtomFlags = atomFlags.find(atomIndex);
            if (thisAtomFlags != atomFlags.end())
                blockExclusions[i][k] |= thisAtomFlags->second;
        }
    }
}

} // namespace OpenMM
