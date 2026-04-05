/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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

#ifndef OPENMM_CPU_CLUSTER_PAIR_LIST_H_
#define OPENMM_CPU_CLUSTER_PAIR_LIST_H_

#include "windowsExportCpu.h"
#include "AlignedArray.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include <atomic>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {

/**
 * Cluster size for the pair list. Must be a power of 2 and <= SIMD width.
 * Using 4 gives a 4x4=16 pair interaction kernel that maps well to SSE/AVX.
 */
static const int CLUSTER_SIZE = 8;

/**
 * A cluster of CLUSTER_SIZE atoms stored in SoA (Structure of Arrays) format
 * for efficient SIMD access.
 */
struct OPENMM_EXPORT_CPU AtomCluster {
    float x[CLUSTER_SIZE];      // x positions
    float y[CLUSTER_SIZE];      // y positions
    float z[CLUSTER_SIZE];      // z positions
    float q[CLUSTER_SIZE];      // charges (scaled by ONE_4PI_EPS0)
    float sigma[CLUSTER_SIZE];  // LJ sigma (half-sigma form)
    float epsilon[CLUSTER_SIZE];// LJ epsilon (sqrt form)
    int atomIndex[CLUSTER_SIZE];// original atom indices
    int size;                   // actual number of atoms (may be < CLUSTER_SIZE for last cluster)

    /**
     * Compute the bounding box center and half-width of this cluster.
     */
    void computeBoundingBox(float& cx, float& cy, float& cz,
                            float& wx, float& wy, float& wz) const {
        float minx = x[0], maxx = x[0];
        float miny = y[0], maxy = y[0];
        float minz = z[0], maxz = z[0];
        for (int i = 1; i < size; i++) {
            if (x[i] < minx) minx = x[i]; else if (x[i] > maxx) maxx = x[i];
            if (y[i] < miny) miny = y[i]; else if (y[i] > maxy) maxy = y[i];
            if (z[i] < minz) minz = z[i]; else if (z[i] > maxz) maxz = z[i];
        }
        cx = 0.5f*(minx+maxx); cy = 0.5f*(miny+maxy); cz = 0.5f*(minz+maxz);
        wx = 0.5f*(maxx-minx); wy = 0.5f*(maxy-miny); wz = 0.5f*(maxz-minz);
    }
};

/**
 * A pair of clusters that should interact. The exclusion mask has one bit per
 * atom pair (i,j) in the 4x4 interaction matrix, where bit (i*4+j) indicates
 * that the pair should be excluded.
 */
/**
 * Compact layout: 16 bytes instead of 24 (supports up to 32767 clusters = 262K atoms).
 * Tighter packing reduces cache misses for large pair lists.
 */
struct ClusterPair {
    uint64_t exclusionMask; // 64-bit mask: bit (i*CLUSTER_SIZE+j) = exclude pair (i,j)
    uint16_t clusterI;      // index of the i-cluster (max 32767 = 262K atoms)
    uint16_t clusterJ;      // index of the j-cluster
    uint8_t jAtomMask;      // which j-atoms have at least one i-atom within cutoff
    uint8_t _pad[3];        // padding to 16-byte alignment
};

/**
 * A cluster-based pair list for nonbonded interactions on the CPU platform.
 *
 * Atoms are organized into spatial clusters of CLUSTER_SIZE atoms. A pair list
 * is built between clusters that are within the cutoff distance. The force
 * kernel then computes all CLUSTER_SIZE^2 pair interactions for each cluster
 * pair simultaneously using SIMD instructions.
 *
 * This approach follows the design described in:
 *   Páll & Hess, Comput. Phys. Commun. 184, 2641–2650 (2013)
 */
class OPENMM_EXPORT_CPU CpuClusterPairList {
public:
    CpuClusterPairList();

    /**
     * Build the cluster pair list from atom positions.
     *
     * @param numAtoms          number of atoms in the system
     * @param posq              interleaved positions+charges [x,y,z,q, x,y,z,q, ...]
     * @param atomParameters    per-atom LJ parameters (sigma, epsilon)
     * @param exclusions        exclusion lists per atom
     * @param periodicBoxVectors the periodic box vectors
     * @param usePeriodic       whether to use periodic boundary conditions
     * @param cutoff            the cutoff distance (including padding)
     * @param threads           thread pool for parallelization
     */
    void build(int numAtoms, const float* posq,
               const std::vector<std::pair<float,float>>& atomParameters,
               const std::vector<std::set<int>>& exclusions,
               const Vec3* periodicBoxVectors, bool usePeriodic,
               float cutoff, ThreadPool& threads);

    /**
     * Get the number of clusters.
     */
    int getNumClusters() const { return clusters.size(); }

    /**
     * Get the cluster data.
     */
    const std::vector<AtomCluster>& getClusters() const { return clusters; }
    std::vector<AtomCluster>& getMutableClusters() { return clusters; }

    /**
     * Get the cluster pair list for a given i-cluster.
     * Returns pairs where clusterI == iCluster.
     */
    const std::vector<ClusterPair>& getClusterPairs() const { return clusterPairs; }

    /**
     * Get the number of cluster pairs.
     */
    int getNumClusterPairs() const { return clusterPairs.size(); }

    /**
     * Get the per-i-cluster work index. Entry i is the start index in clusterPairs
     * for pairs where clusterI has changed from the previous entry. The pairs are
     * grouped by clusterI (sorted) to enable batched i-cluster processing.
     * The last entry equals numPairs (sentinel).
     */
    const std::vector<int>& getIClusterStarts() const { return iClusterStarts; }

    /**
     * Build cluster pair list from an existing CpuNeighborList.
     */
    void buildFromNeighborList(const class CpuNeighborList& neighborList,
                                int numAtoms, const float* posq,
                                const std::vector<std::pair<float,float>>& atomParameters,
                                const std::vector<std::set<int>>& exclusions,
                                float cutoff = 0, const Vec3* periodicBoxVectors = nullptr,
                                bool usePeriodic = false);

    /**
     * Fast direct cluster pair build using Hilbert-order window search.
     * Uses the NL's sorted atoms for cluster assignment, then finds pairs
     * via a window search (exploiting Hilbert spatial locality) + full BB sweep.
     * Skips the expensive NL neighbor traversal (~100ms → ~3ms).
     */
    void buildDirect(const class CpuNeighborList& neighborList,
                     int numAtoms, const float* posq,
                     const std::vector<std::pair<float,float>>& atomParameters,
                     const std::vector<std::set<int>>& exclusions,
                     float cutoff, const Vec3* periodicBoxVectors,
                     bool usePeriodic);

private:
    /**
     * Assign atoms to clusters using spatial binning.
     */
    void assignAtomsToClusters(int numAtoms, const float* posq,
                                const std::vector<std::pair<float,float>>& atomParameters,
                                const Vec3* periodicBoxVectors, bool usePeriodic);

    /**
     * Build the pair list between clusters (parallel).
     */
    void buildPairList(const std::vector<std::set<int>>& exclusions,
                       const Vec3* periodicBoxVectors, bool usePeriodic,
                       float cutoff, ThreadPool& threads);

    std::vector<AtomCluster> clusters;
    std::vector<ClusterPair> clusterPairs;
    std::vector<int> iClusterStarts; // per-i-cluster start indices after sorting by clusterI
    std::atomic<int> atomicCounter;
};

} // namespace OpenMM

#endif // OPENMM_CPU_CLUSTER_PAIR_LIST_H_
