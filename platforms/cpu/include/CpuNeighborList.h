#ifndef OPENMM_CPU_NEIGHBORLIST_H_
#define OPENMM_CPU_NEIGHBORLIST_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2022 Stanford University and the Authors.      *
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

#include "AlignedArray.h"
#include "openmm/Vec3.h"
#include "windowsExportCpu.h"
#include "openmm/internal/ThreadPool.h"
#include <atomic>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT_CPU CpuNeighborList {
public:
    class Voxels;
    class NeighborIterator;
    CpuNeighborList(int blockSize);
    /**
     * Compute the neighbor list based on the current positions of atoms.
     * 
     * @param numAtoms            the number of atoms in the system
     * @param atomLocations       the positions of the atoms
     * @param exclusions          exclusions[i] contains the indices of all atoms with which atom i should not interact
     * @param periodicBoxVectors  the current periodic box vectors
     * @param usePeriodic         whether to apply periodic boundary conditions
     * @param maxDistance         the neighbor list will contain all pairs that are within this distance of each other
     * @param threads             used for parallelization
     */
    void computeNeighborList(int numAtoms, const AlignedArray<float>& atomLocations, const std::vector<std::set<int> >& exclusions,
            const Vec3* periodicBoxVectors, bool usePeriodic, float maxDistance, ThreadPool& threads);
    /**
     * Build a dense neighbor list, in which every atom interacts with every other (except exclusions), regardless of distance.
     * 
     * @param numAtoms            the number of atoms in the system
     * @param exclusions          exclusions[i] contains the indices of all atoms with which atom i should not interact
     */
    void createDenseNeighborList(int numAtoms, const std::vector<std::set<int> >& exclusions);
    int getNumBlocks() const;
    int getBlockSize() const;
    /**
     * Get an object for iterating over the neighbors of an atom block.
     */
    NeighborIterator getNeighborIterator(int blockIndex) const;
    const std::vector<int32_t>& getSortedAtoms() const;
    const std::vector<int>& getBlockNeighbors(int blockIndex) const;

    /**
     * Bitset for a single block, marking which indexes should be excluded. This data type needs to be big
     * enough to store all the bits for any possible block size.
     */
    using BlockExclusionMask = int16_t;

    const std::vector<BlockExclusionMask>& getBlockExclusions(int blockIndex) const;

    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeNeighborList(ThreadPool& threads, int threadIndex);
    void runThread(int index);
private:
    int blockSize;
    std::vector<int> sortedAtoms;
    std::vector<float> sortedPositions;
    std::vector<std::vector<int> > blockNeighbors, blockExclusionIndices;
    std::vector<std::vector<BlockExclusionMask> > blockExclusions;
    // The following variables are used to make information accessible to the individual threads.
    float minx, maxx, miny, maxy, minz, maxz;
    std::vector<std::pair<int, int> > atomBins;
    Voxels* voxels;
    const std::vector<std::set<int> >* exclusions;
    const float* atomLocations;
    Vec3 periodicBoxVectors[3];
    int numAtoms;
    bool usePeriodic, dense;
    float maxDistance;
    std::atomic<int> atomicCounter;
};

class OPENMM_EXPORT_CPU CpuNeighborList::NeighborIterator {
public:
    /**
     * This constructor is used for standard neighbor lists.  Do not call it directly.  Obtain a
     * NeighborIterator by calling getNeighborIterator() on a neighbor list.
     */
    NeighborIterator(const std::vector<int>& neighbors, const std::vector<BlockExclusionMask>& exclusions);
    /**
     * This constructor is used for dense neighbor lists.  Do not call it directly.  Obtain a
     * NeighborIterator by calling getNeighborIterator() on a neighbor list.
     */
    NeighborIterator(int firstAtom, int lastAtom, const std::vector<int>& exclusionIndices, const std::vector<BlockExclusionMask>& exclusions);
    /**
     * Advance the iterator to the next neighbor.
     * 
     * @return false if there are no more neighbors, true otherwise.
     */
    bool next();
    /**
     * Get the index of the current neighbor.
     */
    int getNeighbor() const;
    /**
     * Get bit flags marking which atoms in the block the current atom is excluded from interacting with.
     */
    BlockExclusionMask getExclusions() const;
private:
    bool dense;
    int currentAtom, currentIndex, lastAtom;
    BlockExclusionMask currentExclusions;
    const std::vector<int>* neighbors;
    const std::vector<int>* exclusionIndices;
    const std::vector<BlockExclusionMask>* exclusions;
};

} // namespace OpenMM

#endif // OPENMM_CPU_NEIGHBORLIST_H_
