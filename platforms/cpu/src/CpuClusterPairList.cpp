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

// This file requires AVX2. It is only compiled on x86 with AVX2 flags.
#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))

#include "CpuClusterPairList.h"
#include "CpuNeighborList.h"
#include "hilbert.h"
#include "openmm/internal/vectorizeAvx.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <immintrin.h>
#include <thread>

using namespace std;
using namespace OpenMM;

CpuClusterPairList::CpuClusterPairList() {
}

void CpuClusterPairList::assignAtomsToClusters(int numAtoms, const float* posq,
                                                 const vector<pair<float,float>>& atomParameters,
                                                 const Vec3* periodicBoxVectors, bool usePeriodic) {
    // Find the bounding box of all atoms.
    float minx = posq[0], maxx = posq[0];
    float miny = posq[1], maxy = posq[1];
    float minz = posq[2], maxz = posq[2];
    for (int i = 1; i < numAtoms; i++) {
        float x = posq[4*i], y = posq[4*i+1], z = posq[4*i+2];
        if (x < minx) minx = x; else if (x > maxx) maxx = x;
        if (y < miny) miny = y; else if (y > maxy) maxy = y;
        if (z < minz) minz = z; else if (z > maxz) maxz = z;
    }

    // Grid-based clustering: first assign atoms to spatial grid cells, then
    // Hilbert-sort within each cell. This produces tighter clusters than global Hilbert sort.
    float gridCell = 0.8f;  // tuned: larger cells reduce sort overhead while Hilbert sub-sort keeps clusters tight
    if (usePeriodic) {
        // Fit grid to box
        int ngx = max(1, (int)(periodicBoxVectors[0][0] / gridCell));
        int ngy = max(1, (int)(periodicBoxVectors[1][1] / gridCell));
        int ngz = max(1, (int)(periodicBoxVectors[2][2] / gridCell));
        float gcx = (float)periodicBoxVectors[0][0] / ngx;
        float gcy = (float)periodicBoxVectors[1][1] / ngy;
        float gcz = (float)periodicBoxVectors[2][2] / ngz;
        float invGcx = 1.0f/gcx, invGcy = 1.0f/gcy, invGcz = 1.0f/gcz;

        // Assign atoms to grid cells, then Hilbert-sort within each cell.
        // Use cell index as primary sort key, Hilbert index as secondary.
        float hilbertWidth = max(max(gcx, gcy), gcz) / 15.0f;
        if (hilbertWidth == 0.0f) hilbertWidth = 1.0f;
        float invHW = 1.0f / hilbertWidth;

        vector<pair<int64_t, int>> atomBins(numAtoms);
        for (int i = 0; i < numAtoms; i++) {
            float x = posq[4*i], y = posq[4*i+1], z = posq[4*i+2];
            int cx = min(ngx-1, max(0, (int)(x * invGcx)));
            int cy = min(ngy-1, max(0, (int)(y * invGcy)));
            int cz = min(ngz-1, max(0, (int)(z * invGcz)));
            int cellIdx = cx * ngy * ngz + cy * ngz + cz;

            // Sub-cell Hilbert index (4-bit resolution within cell)
            bitmask_t coords[3];
            coords[0] = (bitmask_t)min(15.0f, max(0.0f, (x - cx*gcx) * invHW));
            coords[1] = (bitmask_t)min(15.0f, max(0.0f, (y - cy*gcy) * invHW));
            coords[2] = (bitmask_t)min(15.0f, max(0.0f, (z - cz*gcz) * invHW));
            int subBin = (int)hilbert_c2i(3, 4, coords);

            atomBins[i] = make_pair(((int64_t)cellIdx << 20) | subBin, i);
        }
        // Two-level counting sort: by cell (O(N+cells)), then by subBin within cell.
        {
            int totalCells = ngx * ngy * ngz;
            // First pass: count atoms per cell.
            vector<int> cc(totalCells + 1, 0);
            for (int i = 0; i < numAtoms; i++) cc[(int)(atomBins[i].first >> 20) + 1]++;
            for (int i = 1; i <= totalCells; i++) cc[i] += cc[i-1];
            // Second pass: place atoms by cell.
            vector<pair<int64_t, int>> temp(numAtoms);
            vector<int> cp(cc.begin(), cc.end());
            for (int i = 0; i < numAtoms; i++) temp[cp[(int)(atomBins[i].first >> 20)]++] = atomBins[i];
            // Third pass: sort within each cell by subBin (cells have ~10 atoms, insertion sort is fast).
            for (int c = 0; c < totalCells; c++) {
                int s = cc[c], e = cc[c+1];
                for (int i = s + 1; i < e; i++) {
                    auto key = temp[i];
                    int j = i - 1;
                    while (j >= s && temp[j].first > key.first) { temp[j+1] = temp[j]; j--; }
                    temp[j+1] = key;
                }
            }
            atomBins = std::move(temp);
        }

        // Build clusters from sorted atoms.
        int numCl = (numAtoms + CLUSTER_SIZE - 1) / CLUSTER_SIZE;
        clusters.resize(numCl);
        for (int c = 0; c < numCl; c++) {
            AtomCluster& cluster = clusters[c];
            int start = c * CLUSTER_SIZE;
            cluster.size = min(CLUSTER_SIZE, numAtoms - start);
            for (int k = 0; k < cluster.size; k++) {
                int atomIdx = atomBins[start + k].second;
                cluster.atomIndex[k] = atomIdx;
                cluster.x[k] = posq[4*atomIdx];
                cluster.y[k] = posq[4*atomIdx+1];
                cluster.z[k] = posq[4*atomIdx+2];
                cluster.q[k] = posq[4*atomIdx+3];
                cluster.sigma[k] = atomParameters[atomIdx].first;
                cluster.epsilon[k] = atomParameters[atomIdx].second;
            }
            for (int k = cluster.size; k < CLUSTER_SIZE; k++) {
                cluster.atomIndex[k] = cluster.atomIndex[0]; // safe: points to real atom
                cluster.x[k] = 1e10f; cluster.y[k] = 1e10f; cluster.z[k] = 1e10f;
                cluster.q[k] = 0.0f; cluster.sigma[k] = 0.0f; cluster.epsilon[k] = 0.0f;
            }
        }
        return;  // skip the global Hilbert path below
    }

    // Non-periodic: global Hilbert sort (original approach).
    float binWidth = max(max(maxx-minx, maxy-miny), maxz-minz)/255.0f;
    if (binWidth == 0.0f) binWidth = 1.0f;
    float invBinWidth = 1.0f/binWidth;

    vector<pair<int, int>> atomBins(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        bitmask_t coords[3];
        coords[0] = (bitmask_t) ((posq[4*i]-minx)*invBinWidth);
        coords[1] = (bitmask_t) ((posq[4*i+1]-miny)*invBinWidth);
        coords[2] = (bitmask_t) ((posq[4*i+2]-minz)*invBinWidth);
        int bin = (int) hilbert_c2i(3, 8, coords);
        atomBins[i] = make_pair(bin, i);
    }
    sort(atomBins.begin(), atomBins.end());

    // Group sorted atoms into clusters of CLUSTER_SIZE.
    int numClusters = (numAtoms + CLUSTER_SIZE - 1) / CLUSTER_SIZE;
    clusters.resize(numClusters);

    for (int c = 0; c < numClusters; c++) {
        AtomCluster& cluster = clusters[c];
        int start = c * CLUSTER_SIZE;
        cluster.size = min(CLUSTER_SIZE, numAtoms - start);

        for (int k = 0; k < cluster.size; k++) {
            int atomIdx = atomBins[start + k].second;
            cluster.atomIndex[k] = atomIdx;
            cluster.x[k] = posq[4*atomIdx];
            cluster.y[k] = posq[4*atomIdx+1];
            cluster.z[k] = posq[4*atomIdx+2];
            cluster.q[k] = posq[4*atomIdx+3];
            cluster.sigma[k] = atomParameters[atomIdx].first;
            cluster.epsilon[k] = atomParameters[atomIdx].second;
        }

        // Pad incomplete clusters with dummy data (far away, zero charge/epsilon).
        for (int k = cluster.size; k < CLUSTER_SIZE; k++) {
            cluster.atomIndex[k] = 0;
            cluster.x[k] = 1e10f;
            cluster.y[k] = 1e10f;
            cluster.z[k] = 1e10f;
            cluster.q[k] = 0.0f;
            cluster.sigma[k] = 0.0f;
            cluster.epsilon[k] = 0.0f;
        }
    }
}

void CpuClusterPairList::buildPairList(const vector<set<int>>& exclusions,
                                         const Vec3* periodicBoxVectors, bool usePeriodic,
                                         float cutoff, ThreadPool& threads) {
    int numClusters = clusters.size();
    float cutoff2 = cutoff * cutoff;

    float boxSizeArr[3] = {0, 0, 0};
    float invBoxSizeArr[3] = {0, 0, 0};
    if (usePeriodic) {
        boxSizeArr[0] = (float) periodicBoxVectors[0][0];
        boxSizeArr[1] = (float) periodicBoxVectors[1][1];
        boxSizeArr[2] = (float) periodicBoxVectors[2][2];
        invBoxSizeArr[0] = 1.0f/boxSizeArr[0];
        invBoxSizeArr[1] = 1.0f/boxSizeArr[1];
        invBoxSizeArr[2] = 1.0f/boxSizeArr[2];
    }

    // Compute bounding box for each cluster.
    vector<float> centerX(numClusters), centerY(numClusters), centerZ(numClusters);
    vector<float> halfWidthX(numClusters), halfWidthY(numClusters), halfWidthZ(numClusters);
    for (int c = 0; c < numClusters; c++) {
        clusters[c].computeBoundingBox(centerX[c], centerY[c], centerZ[c],
                                        halfWidthX[c], halfWidthY[c], halfWidthZ[c]);
    }

    // Compute max cluster half-width for cell sizing.
    float maxHalfW = 0;
    for (int c = 0; c < numClusters; c++)
        maxHalfW = max(maxHalfW, max({halfWidthX[c], halfWidthY[c], halfWidthZ[c]}));

    // Cell size ensures cellRange=2 covers cutoff + cluster diameter.
    float cellSize = max(cutoff * 0.5f, (cutoff + 2.0f*maxHalfW) * 0.5f);
    float minX = centerX[0], maxX = centerX[0];
    float minY = centerY[0], maxY = centerY[0];
    float minZ = centerZ[0], maxZ = centerZ[0];
    for (int c = 1; c < numClusters; c++) {
        if (centerX[c] < minX) minX = centerX[c]; else if (centerX[c] > maxX) maxX = centerX[c];
        if (centerY[c] < minY) minY = centerY[c]; else if (centerY[c] > maxY) maxY = centerY[c];
        if (centerZ[c] < minZ) minZ = centerZ[c]; else if (centerZ[c] > maxZ) maxZ = centerZ[c];
    }

    int ncx, ncy, ncz;
    if (usePeriodic) {
        ncx = max(1, (int)(boxSizeArr[0]/cellSize));
        ncy = max(1, (int)(boxSizeArr[1]/cellSize));
        ncz = max(1, (int)(boxSizeArr[2]/cellSize));
        cellSize = max(max(boxSizeArr[0]/ncx, boxSizeArr[1]/ncy), boxSizeArr[2]/ncz);
    } else {
        ncx = max(1, (int)((maxX-minX)/cellSize) + 1);
        ncy = max(1, (int)((maxY-minY)/cellSize) + 1);
        ncz = max(1, (int)((maxZ-minZ)/cellSize) + 1);
    }

    // Assign clusters to cells, sorted by x-coordinate within each cell.
    // Sorting enables binary search for x-range (like the NL's voxel approach),
    // dramatically reducing the number of candidates checked per cluster.
    vector<vector<pair<float, int>>> cells(ncx*ncy*ncz);
    for (int c = 0; c < numClusters; c++) {
        int cx, cy, cz;
        if (usePeriodic) {
            cx = min(ncx-1, max(0, (int)(centerX[c]*invBoxSizeArr[0]*ncx)));
            cy = min(ncy-1, max(0, (int)(centerY[c]*invBoxSizeArr[1]*ncy)));
            cz = min(ncz-1, max(0, (int)(centerZ[c]*invBoxSizeArr[2]*ncz)));
        } else {
            cx = min(ncx-1, max(0, (int)((centerX[c]-minX)/cellSize)));
            cy = min(ncy-1, max(0, (int)((centerY[c]-minY)/cellSize)));
            cz = min(ncz-1, max(0, (int)((centerZ[c]-minZ)/cellSize)));
        }
        cells[cx*ncy*ncz + cy*ncz + cz].push_back(make_pair(centerX[c], c));
    }
    // Sort each cell by x-coordinate for binary search.
    for (auto& cell : cells)
        sort(cell.begin(), cell.end());
    float searchRadius = cutoff + 2.0f*maxHalfW;

    // Build cluster pairs by iterating over cell pairs.
    int numThreads = threads.getNumThreads();
    vector<vector<ClusterPair>> threadPairs(numThreads);
    atomicCounter = 0;

    int cellRange = (int)ceil(cutoff / cellSize);

    // Check if pre-shift is safe: box must be large enough that cluster-level
    // periodic shift is valid for all atoms within the clusters.
    float maxDiameter = 0;
    for (int c = 0; c < numClusters; c++) {
        float diam = 2.0f * max({halfWidthX[c], halfWidthY[c], halfWidthZ[c]});
        if (diam > maxDiameter) maxDiameter = diam;
    }
    bool safePreShift = !usePeriodic;
    if (usePeriodic) {
        float minHalfBox = 0.5f * min({boxSizeArr[0], boxSizeArr[1], boxSizeArr[2]});
        safePreShift = (minHalfBox > cutoff + 2.0f * maxDiameter);
    }

    int numAtoms = (int)exclusions.size();

    threads.execute([&](ThreadPool& threads, int threadIndex) {
        vector<ClusterPair>& myPairs = threadPairs[threadIndex];
        // Sparse seen-set: track which cj indices were set, clear only those.
        vector<uint8_t> seen(numClusters, 0);
        vector<int> seenList;
        seenList.reserve(512);
        // Flat exclusion lookup per-thread: exclFlat[atomJ] = bitmask of which i-slots exclude atomJ.
        vector<uint8_t> exclFlat(numAtoms, 0);
        vector<int> exclFlaggedAtoms;
        exclFlaggedAtoms.reserve(256);

        while (true) {
            int ci = atomicCounter.fetch_add(1, std::memory_order_relaxed);
            if (ci >= numClusters)
                break;

            // Sparse clear seen flags.
            for (int idx : seenList) seen[idx] = 0;
            seenList.clear();

            // Pre-populate flat exclusion array for atoms in ci.
            for (int idx : exclFlaggedAtoms) exclFlat[idx] = 0;
            exclFlaggedAtoms.clear();
            for (int k = 0; k < clusters[ci].size; k++) {
                int ai = clusters[ci].atomIndex[k];
                for (int excl : exclusions[ai]) {
                    if (exclFlat[excl] == 0)
                        exclFlaggedAtoms.push_back(excl);
                    exclFlat[excl] |= (1 << k);
                }
            }

            float cix = centerX[ci], ciy = centerY[ci], ciz = centerZ[ci];
            float wix = halfWidthX[ci], wiy = halfWidthY[ci], wiz = halfWidthZ[ci];

            // Pre-load i-cluster positions for SIMD distance checks.
            fvec8 ciX(clusters[ci].x), ciY(clusters[ci].y), ciZ(clusters[ci].z);

            // Pre-compute i-cluster padding mask (rows for padding i-atoms).
            uint64_t iPadMask = 0;
            for (int k = clusters[ci].size; k < CLUSTER_SIZE; k++)
                for (int l = 0; l < CLUSTER_SIZE; l++)
                    iPadMask |= (1ULL << (k * CLUSTER_SIZE + l));
            fvec8 cut2v(cutoff2);

            // Add self-interaction pair.
            {
                uint64_t selfExclMask = 0;
                for (int ki = 0; ki < CLUSTER_SIZE; ki++) {
                    for (int kj = 0; kj < CLUSTER_SIZE; kj++) {
                        bool exclude = (ki >= kj) ||
                                       (ki >= clusters[ci].size || kj >= clusters[ci].size);
                        if (!exclude && (exclFlat[clusters[ci].atomIndex[kj]] & (1 << ki)))
                            exclude = true;
                        if (exclude)
                            selfExclMask |= (1ULL << (ki * CLUSTER_SIZE + kj));
                    }
                }
                ClusterPair selfPair;
                selfPair.clusterI = ci;
                selfPair.clusterJ = ci;
                selfPair.exclusionMask = selfExclMask;
                selfPair.jAtomMask = (uint8_t)((1 << clusters[ci].size) - 1);
                myPairs.push_back(selfPair);
            }

            // Find which cell this cluster belongs to.
            int cellX, cellY, cellZ;
            if (usePeriodic) {
                cellX = min(ncx-1, max(0, (int)(cix*invBoxSizeArr[0]*ncx)));
                cellY = min(ncy-1, max(0, (int)(ciy*invBoxSizeArr[1]*ncy)));
                cellZ = min(ncz-1, max(0, (int)(ciz*invBoxSizeArr[2]*ncz)));
            } else {
                cellX = min(ncx-1, max(0, (int)((cix-minX)/cellSize)));
                cellY = min(ncy-1, max(0, (int)((ciy-minY)/cellSize)));
                cellZ = min(ncz-1, max(0, (int)((ciz-minZ)/cellSize)));
            }

            // Check neighboring cells.
            for (int dx = -cellRange; dx <= cellRange; dx++) {
                int nx = cellX + dx;
                if (usePeriodic) { nx = ((nx % ncx) + ncx) % ncx; }
                else if (nx < 0 || nx >= ncx) continue;

                // Pre-compute cell-level periodic shift for x.
                float cellShX = 0;
                if (usePeriodic && safePreShift)
                    cellShX = roundf(((cellX+0.5f)/ncx - (nx+0.5f)/ncx)) * boxSizeArr[0];

                for (int dy = -cellRange; dy <= cellRange; dy++) {
                    int ny = cellY + dy;
                    if (usePeriodic) { ny = ((ny % ncy) + ncy) % ncy; }
                    else if (ny < 0 || ny >= ncy) continue;

                    float cellShY = 0;
                    if (usePeriodic && safePreShift)
                        cellShY = roundf(((cellY+0.5f)/ncy - (ny+0.5f)/ncy)) * boxSizeArr[1];

                    for (int dz = -cellRange; dz <= cellRange; dz++) {
                        int nz = cellZ + dz;
                        if (usePeriodic) { nz = ((nz % ncz) + ncz) % ncz; }
                        else if (nz < 0 || nz >= ncz) continue;

                        float cellShZ = 0;
                        if (usePeriodic && safePreShift)
                            cellShZ = roundf(((cellZ+0.5f)/ncz - (nz+0.5f)/ncz)) * boxSizeArr[2];

                        const auto& cell = cells[nx*ncy*ncz + ny*ncz + nz];
                        // Binary search for x-range within cell (like NL voxel approach).
                        float xshift = usePeriodic ? cellShX : 0;
                        float xlo = cix - searchRadius + xshift;
                        float xhi = cix + searchRadius + xshift;
                        // lower_bound on the x-sorted cell
                        int idxStart = (int)(lower_bound(cell.begin(), cell.end(), make_pair(xlo, 0)) - cell.begin());
                        int idxEnd = (int)(upper_bound(cell.begin(), cell.end(), make_pair(xhi, INT_MAX)) - cell.begin());
                        for (int idx = idxStart; idx < idxEnd; idx++) {
                            int cj = cell[idx].second;
                            if (cj <= ci || seen[cj])
                                continue; // avoid duplicates

                            // Bounding box check using cell-level pre-shift.
                            float shX, shY, shZ;
                            float ddx, ddy, ddz;
                            if (usePeriodic && safePreShift) {
                                // Large box: cell-level shift is exact for all clusters.
                                shX = cellShX; shY = cellShY; shZ = cellShZ;
                                ddx = centerX[cj] - cix + shX;
                                ddy = centerY[cj] - ciy + shY;
                                ddz = centerZ[cj] - ciz + shZ;
                            } else if (usePeriodic) {
                                // Small box: per-cluster shift needed.
                                ddx = centerX[cj] - cix;
                                ddy = centerY[cj] - ciy;
                                ddz = centerZ[cj] - ciz;
                                shX = roundf(ddx*invBoxSizeArr[0])*boxSizeArr[0];
                                shY = roundf(ddy*invBoxSizeArr[1])*boxSizeArr[1];
                                shZ = roundf(ddz*invBoxSizeArr[2])*boxSizeArr[2];
                                ddx -= shX; ddy -= shY; ddz -= shZ;
                            } else {
                                shX = shY = shZ = 0;
                                ddx = centerX[cj] - cix;
                                ddy = centerY[cj] - ciy;
                                ddz = centerZ[cj] - ciz;
                            }
                            float bbDx = max(0.0f, fabsf(ddx) - halfWidthX[cj] - wix);
                            float bbDy = max(0.0f, fabsf(ddy) - halfWidthY[cj] - wiy);
                            float bbDz = max(0.0f, fabsf(ddz) - halfWidthZ[cj] - wiz);
                            if (bbDx*bbDx + bbDy*bbDy + bbDz*bbDz > cutoff2)
                                continue;

                            // Fine-grained AVX check.
                            uint8_t jActive = 0;
                            if (usePeriodic && !safePreShift) {
                                // Small box: need per-lane periodic wrapping.
                                for (int kj = 0; kj < clusters[cj].size; kj++) {
                                    fvec8 dx8 = ciX - fvec8(clusters[cj].x[kj]);
                                    fvec8 dy8 = ciY - fvec8(clusters[cj].y[kj]);
                                    fvec8 dz8 = ciZ - fvec8(clusters[cj].z[kj]);
                                    dx8 -= round(dx8*invBoxSizeArr[0])*boxSizeArr[0];
                                    dy8 -= round(dy8*invBoxSizeArr[1])*boxSizeArr[1];
                                    dz8 -= round(dz8*invBoxSizeArr[2])*boxSizeArr[2];
                                    fvec8 r2 = dx8*dx8 + dy8*dy8 + dz8*dz8;
                                    if (any(ivec8(r2 < cut2v))) jActive |= (1 << kj);
                                }
                            } else {
                                // Large box or non-periodic: pre-shift eliminates per-lane wrapping.
                                for (int kj = 0; kj < clusters[cj].size; kj++) {
                                    fvec8 dx8 = ciX - fvec8(clusters[cj].x[kj] + shX);
                                    fvec8 dy8 = ciY - fvec8(clusters[cj].y[kj] + shY);
                                    fvec8 dz8 = ciZ - fvec8(clusters[cj].z[kj] + shZ);
                                    fvec8 r2 = dx8*dx8 + dy8*dy8 + dz8*dz8;
                                    if (any(ivec8(r2 < cut2v))) jActive |= (1 << kj);
                                }
                            }
                            if (jActive == 0)
                                continue;

                            // Build exclusion mask using flat exclusion array (O(1) per j-atom).
                            uint64_t exclMask = 0;
                            for (int kj = 0; kj < clusters[cj].size; kj++) {
                                uint8_t bits = exclFlat[clusters[cj].atomIndex[kj]];
                                if (bits) {
                                    for (int ki = 0; ki < 8; ki++)
                                        if (bits & (1 << ki))
                                            exclMask |= (1ULL << (ki * CLUSTER_SIZE + kj));
                                }
                            }
                            // Mask out padding atoms (i-padding precomputed, j-padding inline).
                            exclMask |= iPadMask;
                            if (clusters[cj].size < CLUSTER_SIZE) {
                                for (int k = clusters[cj].size; k < CLUSTER_SIZE; k++)
                                    for (int l = 0; l < CLUSTER_SIZE; l++)
                                        exclMask |= (1ULL << (l * CLUSTER_SIZE + k));
                            }

                            ClusterPair cpair;
                            cpair.clusterI = ci;
                            cpair.clusterJ = cj;
                            cpair.exclusionMask = exclMask;
                            cpair.jAtomMask = jActive;
                            myPairs.push_back(cpair);
                            seen[cj] = 1;
                            seenList.push_back(cj);
                        }
                    }
                }
            }
        }
    });
    threads.waitForThreads();

    // Merge thread-local pair lists.
    clusterPairs.clear();
    for (int t = 0; t < numThreads; t++)
        clusterPairs.insert(clusterPairs.end(), threadPairs[t].begin(), threadPairs[t].end());

    // Counting sort by clusterI for grouped kernel processing.
    {
        vector<int> counts(numClusters + 1, 0);
        for (const auto& cp : clusterPairs) counts[cp.clusterI + 1]++;
        for (int i = 1; i <= numClusters; i++) counts[i] += counts[i-1];
        vector<ClusterPair> sorted(clusterPairs.size());
        for (const auto& cp : clusterPairs) sorted[counts[cp.clusterI]++] = cp;
        clusterPairs.swap(sorted);
    }

    // Build i-cluster start index.
    iClusterStarts.clear();
    iClusterStarts.push_back(0);
    for (int i = 1; i < (int)clusterPairs.size(); i++)
        if (clusterPairs[i].clusterI != clusterPairs[i-1].clusterI)
            iClusterStarts.push_back(i);
    iClusterStarts.push_back((int)clusterPairs.size());
}

void CpuClusterPairList::buildFromNeighborList(const CpuNeighborList& neighborList,
                                                int numAtoms, const float* posq,
                                                const vector<pair<float,float>>& atomParameters,
                                                const vector<set<int>>& exclusions,
                                                float cutoff, const Vec3* periodicBoxVectors,
                                                bool usePeriodic) {
    // The existing neighbor list already has blocks of 8 atoms (Hilbert-sorted)
    // and per-block neighbor atom lists. Convert to cluster pairs.

    int nlBlockSize = neighborList.getBlockSize();
    int nlNumBlocks = neighborList.getNumBlocks();
    const vector<int32_t>& sortedAtoms = neighborList.getSortedAtoms();

    // Build CLUSTER_SIZE-atom clusters by merging consecutive NL blocks if needed.
    // NL blockSize may be 4 (SSE) while CLUSTER_SIZE is 8 (AVX).
    int numBlocks = (numAtoms + CLUSTER_SIZE - 1) / CLUSTER_SIZE;
    clusters.resize(numBlocks);
    for (int b = 0; b < numBlocks; b++) {
        AtomCluster& cluster = clusters[b];
        int start = b * CLUSTER_SIZE;
        cluster.size = min(CLUSTER_SIZE, numAtoms - start);
        for (int k = 0; k < cluster.size; k++) {
            int sortIdx = start + k;
            int atomIdx = (sortIdx < (int)sortedAtoms.size()) ? sortedAtoms[sortIdx] : 0;
            cluster.atomIndex[k] = atomIdx;
            cluster.x[k] = posq[4*atomIdx];
            cluster.y[k] = posq[4*atomIdx+1];
            cluster.z[k] = posq[4*atomIdx+2];
            cluster.q[k] = posq[4*atomIdx+3];
            cluster.sigma[k] = atomParameters[atomIdx].first;
            cluster.epsilon[k] = atomParameters[atomIdx].second;
        }
        for (int k = cluster.size; k < CLUSTER_SIZE; k++) {
            cluster.atomIndex[k] = 0;
            cluster.x[k] = 1e10f; cluster.y[k] = 1e10f; cluster.z[k] = 1e10f;
            cluster.q[k] = 0; cluster.sigma[k] = 0; cluster.epsilon[k] = 0;
        }
    }

    // Build atom-to-block mapping.
    vector<int> atomToBlock(numAtoms, -1);
    vector<int> atomToSlot(numAtoms, -1);
    for (int b = 0; b < numBlocks; b++)
        for (int k = 0; k < clusters[b].size; k++) {
            atomToBlock[clusters[b].atomIndex[k]] = b;
            atomToSlot[clusters[b].atomIndex[k]] = k;
        }

    // Convert block-atom neighbor lists to cluster-cluster pairs.
    clusterPairs.clear();

    // Pre-compute bounding boxes for distance pruning.
    vector<float> centerX(numBlocks), centerY(numBlocks), centerZ(numBlocks);
    vector<float> halfWX(numBlocks), halfWY(numBlocks), halfWZ(numBlocks);
    for (int b = 0; b < numBlocks; b++)
        clusters[b].computeBoundingBox(centerX[b], centerY[b], centerZ[b],
                                       halfWX[b], halfWY[b], halfWZ[b]);
    float bbCutoff2 = (cutoff > 0) ? cutoff*cutoff : 0;
    float bsx = 0, bsy = 0, bsz = 0, ibsx = 0, ibsy = 0, ibsz = 0;
    bool needWrap = false;
    if (usePeriodic && periodicBoxVectors) {
        bsx = (float)periodicBoxVectors[0][0]; bsy = (float)periodicBoxVectors[1][1];
        bsz = (float)periodicBoxVectors[2][2];
        ibsx = 1.0f/bsx; ibsy = 1.0f/bsy; ibsz = 1.0f/bsz;
        // Check if per-atom wrapping is needed for SIMD distance checks.
        // If any cluster has spread comparable to box size, the simple per-cluster
        // shift becomes inaccurate and we must wrap per-atom.
        float maxD = 0;
        int nS = min((int)clusters.size(), 200);
        for (int c = 0; c < nS; c++)
            for (int k = 1; k < clusters[c].size; k++) {
                float d2 = (clusters[c].x[k]-clusters[c].x[0])*(clusters[c].x[k]-clusters[c].x[0])
                          +(clusters[c].y[k]-clusters[c].y[0])*(clusters[c].y[k]-clusters[c].y[0])
                          +(clusters[c].z[k]-clusters[c].z[0])*(clusters[c].z[k]-clusters[c].z[0]);
                if (d2 > maxD) maxD = d2;
            }
        needWrap = (0.5f*min({bsx,bsy,bsz}) <= cutoff + 2*sqrtf(maxD));
    }

    // Map NL blocks to our clusters.
    vector<vector<int>> clusterNLBlocks(numBlocks);
    for (int nb = 0; nb < nlNumBlocks; nb++) {
        int firstSortIdx = nb * nlBlockSize;
        if (firstSortIdx < (int)sortedAtoms.size() && firstSortIdx < numAtoms) {
            int ci = atomToBlock[sortedAtoms[firstSortIdx]];
            if (ci >= 0 && ci < numBlocks)
                clusterNLBlocks[ci].push_back(nb);
        }
    }

    // Parallel pair generation: each thread processes a range of i-clusters.
    // Thread-local dedup arrays and exclusion tracking avoid synchronization.
    int nBuildThreads = min(4, max(1, (int)thread::hardware_concurrency()));
    int biChunk = (numBlocks + nBuildThreads - 1) / nBuildThreads;
    vector<vector<ClusterPair>> threadPairs(nBuildThreads);

    auto buildWork = [&](int startBI, int endBI, int threadIdx) {
        auto& myPairs = threadPairs[threadIdx];
        myPairs.reserve((endBI - startBI) * 100); // rough estimate

        // Thread-local dedup: generation counter + per-block gen stamp.
        vector<int> myPairGen(numBlocks, -1);
        int myGen = 0;

        // Thread-local exclusion tracking.
        vector<uint8_t> myExclFlat(numAtoms, 0);
        vector<int> myExclFlagged;

        for (int bi = startBI; bi < endBI; bi++) {
            myGen++;
            // Pre-populate exclusion flags for atoms in this cluster.
            myExclFlagged.clear();
            for (int k = 0; k < clusters[bi].size; k++) {
                int ai = clusters[bi].atomIndex[k];
                for (int excl : exclusions[ai]) {
                    if (myExclFlat[excl] == 0)
                        myExclFlagged.push_back(excl);
                    myExclFlat[excl] |= (1 << k);
                }
            }

            // Self-pair.
            {
                uint64_t selfExclMask = 0;
                for (int ki = 0; ki < CLUSTER_SIZE; ki++)
                    for (int kj = 0; kj < CLUSTER_SIZE; kj++) {
                        bool exclude = (ki >= kj) ||
                                       (ki >= clusters[bi].size || kj >= clusters[bi].size);
                        if (!exclude) {
                            int ai = clusters[bi].atomIndex[ki];
                            int aj = clusters[bi].atomIndex[kj];
                            if (exclusions[ai].count(aj)) exclude = true;
                        }
                        if (exclude) selfExclMask |= (1ULL << (ki*CLUSTER_SIZE+kj));
                    }
                ClusterPair sp;
                sp.clusterI = bi;
                sp.clusterJ = bi;
                sp.jAtomMask = (uint8_t)((1 << clusters[bi].size) - 1);
                sp.exclusionMask = selfExclMask;
                myPairs.push_back(sp);
            }

            // Inter-cluster pairs from NL.
            for (int nbIdx = 0; nbIdx < (int)clusterNLBlocks[bi].size(); nbIdx++) {
                int nb = clusterNLBlocks[bi][nbIdx];
                const vector<int>& neighbors = neighborList.getBlockNeighbors(nb);
                for (int idx = 0; idx < (int)neighbors.size(); idx++) {
                    int neighborAtom = neighbors[idx];
                    int bj = atomToBlock[neighborAtom];
                    if (bj < 0 || bj >= bi) continue;

                    if (myPairGen[bj] == myGen)
                        continue; // dedup

                    // BB distance pruning.
                    if (bbCutoff2 > 0) {
                        float ddx = centerX[bi] - centerX[bj];
                        float ddy = centerY[bi] - centerY[bj];
                        float ddz = centerZ[bi] - centerZ[bj];
                        if (usePeriodic) {
                            ddx -= roundf(ddx*ibsx)*bsx;
                            ddy -= roundf(ddy*ibsy)*bsy;
                            ddz -= roundf(ddz*ibsz)*bsz;
                        }
                        float bbDx = max(0.0f, fabsf(ddx) - halfWX[bj] - halfWX[bi]);
                        float bbDy = max(0.0f, fabsf(ddy) - halfWY[bj] - halfWY[bi]);
                        float bbDz = max(0.0f, fabsf(ddz) - halfWZ[bj] - halfWZ[bi]);
                        if (bbDx*bbDx + bbDy*bbDy + bbDz*bbDz > bbCutoff2) {
                            myPairGen[bj] = myGen;
                            continue;
                        }
                    }

                    // Build exclusion mask.
                    uint64_t exclMask = 0;
                    if (!myExclFlagged.empty()) {
                        for (int ki = 0; ki < clusters[bj].size; ki++) {
                            int ai = clusters[bj].atomIndex[ki];
                            uint8_t bits = myExclFlat[ai];
                            if (bits) {
                                for (int kj = 0; kj < 8; kj++)
                                    if (bits & (1 << kj))
                                        exclMask |= (1ULL << (ki*CLUSTER_SIZE+kj));
                            }
                        }
                    }
                    for (int k = clusters[bj].size; k < CLUSTER_SIZE; k++)
                        for (int l = 0; l < CLUSTER_SIZE; l++)
                            exclMask |= (1ULL << (k*CLUSTER_SIZE+l));
                    for (int k = clusters[bi].size; k < CLUSTER_SIZE; k++)
                        for (int l = 0; l < CLUSTER_SIZE; l++)
                            exclMask |= (1ULL << (l*CLUSTER_SIZE+k));

                    // Inline j-mask refinement: SIMD distance check during pair creation.
                    // This eliminates the separate 12ms refinement pass.
                    uint8_t jActive = 0;
                    if (bbCutoff2 > 0) {
                        const auto& ciRef = clusters[bj];  // i-cluster in pair
                        const auto& cjRef = clusters[bi];   // j-cluster in pair
                        fvec8 ciX(ciRef.x), ciY(ciRef.y), ciZ(ciRef.z);
                        fvec8 cut2v(bbCutoff2);
                        float shX = 0, shY = 0, shZ = 0;
                        if (usePeriodic && !needWrap) {
                            shX = roundf((ciRef.x[0]-cjRef.x[0])*ibsx)*bsx;
                            shY = roundf((ciRef.y[0]-cjRef.y[0])*ibsy)*bsy;
                            shZ = roundf((ciRef.z[0]-cjRef.z[0])*ibsz)*bsz;
                        }
                        for (int kj = 0; kj < cjRef.size; kj++) {
                            fvec8 dx8, dy8, dz8;
                            if (usePeriodic && needWrap) {
                                dx8 = ciX - fvec8(cjRef.x[kj]); dy8 = ciY - fvec8(cjRef.y[kj]); dz8 = ciZ - fvec8(cjRef.z[kj]);
                                dx8 -= round(dx8*ibsx)*bsx; dy8 -= round(dy8*ibsy)*bsy; dz8 -= round(dz8*ibsz)*bsz;
                            } else {
                                dx8 = ciX - fvec8(cjRef.x[kj]+shX); dy8 = ciY - fvec8(cjRef.y[kj]+shY); dz8 = ciZ - fvec8(cjRef.z[kj]+shZ);
                            }
                            if (any(ivec8(dx8*dx8 + dy8*dy8 + dz8*dz8 < cut2v)))
                                jActive |= (1 << kj);
                        }
                    } else {
                        jActive = (uint8_t)((1 << clusters[bi].size) - 1);
                    }

                    myPairGen[bj] = myGen;
                    if (jActive == 0)
                        continue; // no j-atoms within cutoff - skip entirely

                    ClusterPair cp;
                    cp.clusterI = bj;
                    cp.clusterJ = bi;
                    cp.exclusionMask = exclMask;
                    cp.jAtomMask = jActive;
                    myPairs.push_back(cp);
                }
            }

            // Clear exclusion flags for next cluster.
            for (int idx : myExclFlagged)
                myExclFlat[idx] = 0;
        }
    };

    // Launch threads.
    {
        vector<thread> buildThreads;
        for (int t = 0; t < nBuildThreads; t++) {
            int s = t * biChunk, e = min(s + biChunk, numBlocks);
            if (s < numBlocks)
                buildThreads.emplace_back(buildWork, s, e, t);
        }
        for (auto& th : buildThreads) th.join();
    }

    // Merge per-thread pairs.
    {
        int totalPairs = 0;
        for (auto& tp : threadPairs) totalPairs += (int)tp.size();
        clusterPairs.reserve(totalPairs);
        for (auto& tp : threadPairs) {
            clusterPairs.insert(clusterPairs.end(), tp.begin(), tp.end());
            tp.clear(); tp.shrink_to_fit(); // free memory
        }
    }

    // Sort pairs by clusterI for grouped processing in the force kernel.
    // j-mask refinement is now done inline during pair generation above.
    {
        int numCl = (int)clusters.size();
        vector<int> counts(numCl + 1, 0);
        for (const auto& cp : clusterPairs)
            counts[cp.clusterI + 1]++;
        for (int i = 1; i <= numCl; i++)
            counts[i] += counts[i-1];
        vector<ClusterPair> sorted(clusterPairs.size());
        for (const auto& cp : clusterPairs)
            sorted[counts[cp.clusterI]++] = cp;
        clusterPairs.swap(sorted);
    }

    // Build i-cluster start index.
    iClusterStarts.clear();
    iClusterStarts.push_back(0);
    for (int i = 1; i < (int)clusterPairs.size(); i++) {
        if (clusterPairs[i].clusterI != clusterPairs[i-1].clusterI)
            iClusterStarts.push_back(i);
    }
    iClusterStarts.push_back((int)clusterPairs.size()); // sentinel
}

void CpuClusterPairList::build(int numAtoms, const float* posq,
                                const vector<pair<float,float>>& atomParameters,
                                const vector<set<int>>& exclusions,
                                const Vec3* periodicBoxVectors, bool usePeriodic,
                                float cutoff, ThreadPool& threads) {
    assignAtomsToClusters(numAtoms, posq, atomParameters, periodicBoxVectors, usePeriodic);
    buildPairList(exclusions, periodicBoxVectors, usePeriodic, cutoff, threads);
}

void CpuClusterPairList::buildDirect(const CpuNeighborList& neighborList,
                                      int numAtoms, const float* posq,
                                      const vector<pair<float,float>>& atomParameters,
                                      const vector<set<int>>& exclusions,
                                      float cutoff, const Vec3* periodicBoxVectors,
                                      bool usePeriodic) {
    // Reuse NL's Hilbert-sorted atoms for cluster assignment (same as buildFromNeighborList).
    const vector<int32_t>& sortedAtoms = neighborList.getSortedAtoms();
    int numBlocks = (numAtoms + CLUSTER_SIZE - 1) / CLUSTER_SIZE;
    clusters.resize(numBlocks);
    for (int b = 0; b < numBlocks; b++) {
        AtomCluster& cluster = clusters[b];
        int start = b * CLUSTER_SIZE;
        cluster.size = min(CLUSTER_SIZE, numAtoms - start);
        for (int k = 0; k < cluster.size; k++) {
            int sortIdx = start + k;
            int atomIdx = (sortIdx < (int)sortedAtoms.size()) ? sortedAtoms[sortIdx] : 0;
            cluster.atomIndex[k] = atomIdx;
            cluster.x[k] = posq[4*atomIdx];
            cluster.y[k] = posq[4*atomIdx+1];
            cluster.z[k] = posq[4*atomIdx+2];
            cluster.q[k] = posq[4*atomIdx+3];
            cluster.sigma[k] = atomParameters[atomIdx].first;
            cluster.epsilon[k] = atomParameters[atomIdx].second;
        }
        for (int k = cluster.size; k < CLUSTER_SIZE; k++) {
            cluster.atomIndex[k] = 0;
            cluster.x[k] = 1e10f; cluster.y[k] = 1e10f; cluster.z[k] = 1e10f;
            cluster.q[k] = 0; cluster.sigma[k] = 0; cluster.epsilon[k] = 0;
        }
    }

    // Precompute bounding boxes.
    vector<float> centerX(numBlocks), centerY(numBlocks), centerZ(numBlocks);
    vector<float> halfWX(numBlocks), halfWY(numBlocks), halfWZ(numBlocks);
    for (int b = 0; b < numBlocks; b++)
        clusters[b].computeBoundingBox(centerX[b], centerY[b], centerZ[b],
                                       halfWX[b], halfWY[b], halfWZ[b]);

    float bbCutoff2 = cutoff * cutoff;
    float bsx = 0, bsy = 0, bsz = 0, ibsx = 0, ibsy = 0, ibsz = 0;
    bool needWrap = false;
    if (usePeriodic && periodicBoxVectors) {
        bsx = (float)periodicBoxVectors[0][0]; bsy = (float)periodicBoxVectors[1][1];
        bsz = (float)periodicBoxVectors[2][2];
        ibsx = 1.0f/bsx; ibsy = 1.0f/bsy; ibsz = 1.0f/bsz;
        float maxD = 0;
        int nS = min(numBlocks, 200);
        for (int c = 0; c < nS; c++)
            for (int k = 1; k < clusters[c].size; k++) {
                float d2 = (clusters[c].x[k]-clusters[c].x[0])*(clusters[c].x[k]-clusters[c].x[0])
                          +(clusters[c].y[k]-clusters[c].y[0])*(clusters[c].y[k]-clusters[c].y[0])
                          +(clusters[c].z[k]-clusters[c].z[0])*(clusters[c].z[k]-clusters[c].z[0]);
                if (d2 > maxD) maxD = d2;
            }
        needWrap = (0.5f*min({bsx,bsy,bsz}) <= cutoff + 2*sqrtf(maxD));
    }

    // Build atom-to-block mapping for exclusion computation.
    vector<int> atomToBlock(numAtoms, -1);
    for (int b = 0; b < numBlocks; b++)
        for (int k = 0; k < clusters[b].size; k++)
            atomToBlock[clusters[b].atomIndex[k]] = b;

    // Direct pair build: parallel over i-clusters.
    // For each i-cluster, check ALL other clusters via BB distance.
    // The O(N²) BB checks are cheap (~20 ops each, ~1106² total = 1.2M ops).
    int nBuildThreads = min(4, max(1, (int)thread::hardware_concurrency()));
    int biChunk = (numBlocks + nBuildThreads - 1) / nBuildThreads;
    vector<vector<ClusterPair>> threadPairs(nBuildThreads);

    auto directWork = [&](int startBI, int endBI, int threadIdx) {
        auto& myPairs = threadPairs[threadIdx];
        myPairs.reserve((endBI - startBI) * 200);
        vector<uint8_t> exclFlat(numAtoms, 0);
        vector<int> exclFlagged;

        for (int bi = startBI; bi < endBI; bi++) {
            // Build exclusion flat for this cluster.
            exclFlagged.clear();
            for (int k = 0; k < clusters[bi].size; k++) {
                int ai = clusters[bi].atomIndex[k];
                for (int excl : exclusions[ai]) {
                    if (exclFlat[excl] == 0) exclFlagged.push_back(excl);
                    exclFlat[excl] |= (1 << k);
                }
            }

            // Self-pair.
            {
                uint64_t selfExclMask = 0;
                for (int ki = 0; ki < CLUSTER_SIZE; ki++)
                    for (int kj = 0; kj < CLUSTER_SIZE; kj++) {
                        bool exclude = (ki >= kj) || (ki >= clusters[bi].size || kj >= clusters[bi].size);
                        if (!exclude && exclusions[clusters[bi].atomIndex[ki]].count(clusters[bi].atomIndex[kj]))
                            exclude = true;
                        if (exclude) selfExclMask |= (1ULL << (ki*CLUSTER_SIZE+kj));
                    }
                ClusterPair sp;
                sp.clusterI = bi; sp.clusterJ = bi;
                sp.exclusionMask = selfExclMask;
                sp.jAtomMask = (uint8_t)((1 << clusters[bi].size) - 1);
                myPairs.push_back(sp);
            }

            // Check all other clusters (bj < bi for N3L).
            for (int bj = 0; bj < bi; bj++) {
                // BB distance check.
                float ddx = centerX[bi] - centerX[bj];
                float ddy = centerY[bi] - centerY[bj];
                float ddz = centerZ[bi] - centerZ[bj];
                if (usePeriodic) {
                    ddx -= roundf(ddx*ibsx)*bsx;
                    ddy -= roundf(ddy*ibsy)*bsy;
                    ddz -= roundf(ddz*ibsz)*bsz;
                }
                float bbDx = max(0.0f, fabsf(ddx) - halfWX[bj] - halfWX[bi]);
                float bbDy = max(0.0f, fabsf(ddy) - halfWY[bj] - halfWY[bi]);
                float bbDz = max(0.0f, fabsf(ddz) - halfWZ[bj] - halfWZ[bi]);
                if (bbDx*bbDx + bbDy*bbDy + bbDz*bbDz > bbCutoff2)
                    continue;

                // Build exclusion mask.
                uint64_t exclMask = 0;
                for (int ki = 0; ki < clusters[bj].size; ki++) {
                    uint8_t bits = exclFlat[clusters[bj].atomIndex[ki]];
                    if (bits) {
                        for (int kj = 0; kj < 8; kj++)
                            if (bits & (1 << kj))
                                exclMask |= (1ULL << (ki*CLUSTER_SIZE+kj));
                    }
                }
                for (int k = clusters[bj].size; k < CLUSTER_SIZE; k++)
                    for (int l = 0; l < CLUSTER_SIZE; l++)
                        exclMask |= (1ULL << (k*CLUSTER_SIZE+l));
                for (int k = clusters[bi].size; k < CLUSTER_SIZE; k++)
                    for (int l = 0; l < CLUSTER_SIZE; l++)
                        exclMask |= (1ULL << (l*CLUSTER_SIZE+k));

                // J-mask SIMD refinement.
                uint8_t jActive = 0;
                const auto& ciRef = clusters[bj];
                const auto& cjRef = clusters[bi];
                fvec8 ciX(ciRef.x), ciY(ciRef.y), ciZ(ciRef.z);
                fvec8 cut2v(bbCutoff2);
                float shX = 0, shY = 0, shZ = 0;
                if (usePeriodic && !needWrap) {
                    shX = roundf((ciRef.x[0]-cjRef.x[0])*ibsx)*bsx;
                    shY = roundf((ciRef.y[0]-cjRef.y[0])*ibsy)*bsy;
                    shZ = roundf((ciRef.z[0]-cjRef.z[0])*ibsz)*bsz;
                }
                for (int kj = 0; kj < cjRef.size; kj++) {
                    fvec8 dx8, dy8, dz8;
                    if (usePeriodic && needWrap) {
                        dx8 = ciX - fvec8(cjRef.x[kj]); dy8 = ciY - fvec8(cjRef.y[kj]); dz8 = ciZ - fvec8(cjRef.z[kj]);
                        dx8 -= round(dx8*ibsx)*bsx; dy8 -= round(dy8*ibsy)*bsy; dz8 -= round(dz8*ibsz)*bsz;
                    } else {
                        dx8 = ciX - fvec8(cjRef.x[kj]+shX); dy8 = ciY - fvec8(cjRef.y[kj]+shY); dz8 = ciZ - fvec8(cjRef.z[kj]+shZ);
                    }
                    if (any(ivec8(dx8*dx8 + dy8*dy8 + dz8*dz8 < cut2v)))
                        jActive |= (1 << kj);
                }
                if (jActive == 0) continue;

                ClusterPair cp;
                cp.clusterI = bj; cp.clusterJ = bi;
                cp.exclusionMask = exclMask;
                cp.jAtomMask = jActive;
                myPairs.push_back(cp);
            }

            for (int idx : exclFlagged) exclFlat[idx] = 0;
        }
    };

    // Launch parallel threads.
    {
        vector<thread> buildThreads;
        for (int t = 0; t < nBuildThreads; t++) {
            int s = t * biChunk, e = min(s + biChunk, numBlocks);
            if (s < numBlocks)
                buildThreads.emplace_back(directWork, s, e, t);
        }
        for (auto& th : buildThreads) th.join();
    }

    // Merge + sort.
    clusterPairs.clear();
    {
        int totalPairs = 0;
        for (auto& tp : threadPairs) totalPairs += (int)tp.size();
        clusterPairs.reserve(totalPairs);
        for (auto& tp : threadPairs) {
            clusterPairs.insert(clusterPairs.end(), tp.begin(), tp.end());
            tp.clear(); tp.shrink_to_fit();
        }
    }
    {
        vector<int> counts(numBlocks + 1, 0);
        for (const auto& cp : clusterPairs) counts[cp.clusterI + 1]++;
        for (int i = 1; i <= numBlocks; i++) counts[i] += counts[i-1];
        vector<ClusterPair> sorted(clusterPairs.size());
        for (const auto& cp : clusterPairs) sorted[counts[cp.clusterI]++] = cp;
        clusterPairs.swap(sorted);
    }
    iClusterStarts.clear();
    iClusterStarts.push_back(0);
    for (int i = 1; i < (int)clusterPairs.size(); i++)
        if (clusterPairs[i].clusterI != clusterPairs[i-1].clusterI)
            iClusterStarts.push_back(i);
    iClusterStarts.push_back((int)clusterPairs.size());
}

#endif // __AVX2__
