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
    CpuNeighborList(int blockSize);
    void computeNeighborList(int numAtoms, const AlignedArray<float>& atomLocations, const std::vector<std::set<int> >& exclusions,
            const Vec3* periodicBoxVectors, bool usePeriodic, float maxDistance, ThreadPool& threads);
    int getNumBlocks() const;
    int getBlockSize() const;
    const std::vector<int>& getSortedAtoms() const;
    const std::vector<int>& getBlockNeighbors(int blockIndex) const;
    const std::vector<char>& getBlockExclusions(int blockIndex) const;
    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeNeighborList(ThreadPool& threads, int threadIndex);
    void runThread(int index);
private:
    int blockSize;
    std::vector<int> sortedAtoms;
    std::vector<float> sortedPositions;
    std::vector<std::vector<int> > blockNeighbors;
    std::vector<std::vector<char> > blockExclusions;
    // The following variables are used to make information accessible to the individual threads.
    float minx, maxx, miny, maxy, minz, maxz;
    std::vector<std::pair<int, int> > atomBins;
    Voxels* voxels;
    const std::vector<std::set<int> >* exclusions;
    const float* atomLocations;
    Vec3 periodicBoxVectors[3];
    int numAtoms;
    bool usePeriodic;
    float maxDistance;
    std::atomic<int> atomicCounter;
};

} // namespace OpenMM

#endif // OPENMM_CPU_NEIGHBORLIST_H_
