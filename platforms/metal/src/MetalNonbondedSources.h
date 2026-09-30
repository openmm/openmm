/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Ported from OpenMM's OpenCL/CUDA/HIP nonbonded kernels and utilities.
 * Original code: Portions copyright (c) 2009-2025 Stanford University and
 * the Authors. Authors: Peter Eastman and the OpenMM contributors.
 * Metal code: Portions copyright (c) 2026 Chun-Chi Hung.
 * Author: Chun-Chi Hung.
 * This program is free software under the GNU Lesser General Public License,
 * version 3 or (at your option) any later version. It is distributed without
 * any warranty; see <http://www.gnu.org/licenses/> for the license.
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALNONBONDEDSOURCES_H_
#define OPENMM_METALNONBONDEDSOURCES_H_

#include "MetalCudaKernelSources.h"
#include "MetalKernelSources.h"
#include "openmm/OpenMMException.h"
#include <string>

namespace OpenMM {
namespace MetalNonbondedSources {

/** @brief Apply a checked edit to a known upstream template, never arbitrary user code. */
inline void replace(std::string& source, const std::string& before, const std::string& after) {
    size_t offset = source.find(before);
    if (offset == std::string::npos)
        throw OpenMMException("OpenMM nonbonded template changed: missing Metal adaptation for "+before);
    do {
        source.replace(offset, before.size(), after);
        offset = source.find(before, offset+after.size());
    } while (offset != std::string::npos);
}

/** @brief Adapt only CUDA's sparse-pair fragments; tiled interactions remain OpenCL. */
inline void addOpenCLPairs(std::string& source) {
    const std::string& cuda = MetalCudaKernelSources::nonbonded;
    size_t helperBegin = cuda.find("__device__ void saveSingleForce");
    size_t helperEnd = cuda.find("/**", helperBegin);
    size_t pairsBegin = cuda.find("    // Third loop: single pairs");
    size_t pairsEnd = cuda.find("#ifdef INCLUDE_ENERGY", pairsBegin);
    if (helperBegin == std::string::npos || helperEnd == std::string::npos ||
            pairsBegin == std::string::npos || pairsEnd == std::string::npos)
        throw OpenMMException("CUDA sparse-pair template changed");
    std::string helper = cuda.substr(helperBegin, helperEnd-helperBegin);
    replace(helper, "static_cast<unsigned long long>", "(mm_ulong)");
    replace(helper, "unsigned long long", "mm_ulong");
    replace(helper, "__device__", "DEVICE");
    replace(helper, "atomicAdd(", "ATOMIC_ADD(");
    for (const std::string offset : {"atom", "atom+PADDED_NUM_ATOMS", "atom+2*PADDED_NUM_ATOMS"})
        replace(helper, "ATOMIC_ADD(&forceBuffers["+offset+"],",
            "METAL_ACCUMULATE_SPARSE_FORCE(forceBuffers, "+offset+",");
    replace(helper, "mm_ulong* forceBuffers", "GLOBAL mm_ulong* forceBuffers");
    std::string pairs = cuda.substr(pairsBegin, pairsEnd-pairsBegin);
    replace(pairs, "blockIdx.x", "get_group_id(0)");
    replace(pairs, "blockDim.x", "get_local_size(0)");
    replace(pairs, "gridDim.x", "get_num_groups(0)");
    replace(pairs, "threadIdx.x", "get_local_id(0)");
    replace(pairs, "#if USE_NEIGHBOR_LIST", "#ifdef USE_SPARSE_PAIRS");
    const std::string anchor = "#ifdef INCLUDE_ENERGY\n    energyBuffer[get_global_id(0)] += energy;";
    replace(source, anchor, pairs+anchor);
    replace(source, "__global const int* restrict interactingAtoms\n#endif",
        "__global const int* restrict interactingAtoms\n"
        "#ifdef USE_SPARSE_PAIRS\n, unsigned int maxSinglePairs, GLOBAL const int2* singlePairs\n#endif\n#endif");
    source = helper+source;
}

/** @brief Shared by the production compression and its boundary-value regression. */
inline std::string halfBoundsSource() {
    return R"(
/** Round nonnegative bounds upward, preserving CUDA/HIP's conservative rule. */
DEVICE half metalBoundHalf(real value) {
    half result = half(value);
    ushort bits = as_type<ushort>(result);
    if (float(result) < value) result = as_type<half>(ushort(bits+1));
    return result;
}
DEVICE half4 metalBoundsHalf(real4 value) {
    return half4(metalBoundHalf(value.x), metalBoundHalf(value.y), metalBoundHalf(value.z), half(0));
}
)";
}

/** @brief Apply independent bounds, ballot/compaction, and sparse-pair options to OpenCL. */
inline std::string neighbors(std::string source, bool sparsePairs) {
#if OPENMM_METAL_FAST_BLOCK_BOUNDS
    {
        size_t first = source.find("__kernel void findBlockBounds(");
        size_t last = source.find("__kernel void computeSortKeys(", first);
        if (first == std::string::npos || last == std::string::npos)
            throw OpenMMException("OpenCL block bounds template changed");
        source.replace(first, last-first, MetalKernelSources::neighborBounds+"\n");
    }
#endif
#if OPENMM_METAL_FAST_FP16_BOUNDS
    // Only sorted and large boxes are compressed. Public Common block bounds
    // keep float4 storage; every consumer expands half4 explicitly to float4.
    replace(source, "real4* restrict sortedBlockBoundingBox", "half4* restrict sortedBlockBoundingBox");
    replace(source, "real4* restrict largeBlockBoundingBox", "half4* restrict largeBlockBoundingBox");
    replace(source, "sortedBlockBoundingBox[i] = blockBoundingBox[index];",
        "sortedBlockBoundingBox[i] = metalBoundsHalf(blockBoundingBox[index]);");
    replace(source, "largeBlockBoundingBox[i] = 0.5f*(maxPos-minPos);",
        "largeBlockBoundingBox[i] = metalBoundsHalf(0.5f*(maxPos-minPos));");
    replace(source, "= sortedBlockBoundingBox[block1];", "= float4(sortedBlockBoundingBox[block1]);");
    replace(source, "= sortedBlockBoundingBox[block2];", "= float4(sortedBlockBoundingBox[block2]);");
    replace(source, "= largeBlockBoundingBox[largeBlockIndex];", "= float4(largeBlockBoundingBox[largeBlockIndex]);");
    source = halfBoundsSource()+source;
#endif
    if (sparsePairs) {
        // Count sparse atom-to-block interactions before OpenCL's compaction.
        // This deliberately does not depend on the optional ballot algorithm.
        replace(source, "interactionCount[0] = 0;", "interactionCount[0] = 0;\n        interactionCount[1] = 0;");
        const std::string argumentEnd = "__global const int* restrict rebuildNeighborList\n#ifdef USE_LARGE_BLOCKS";
        replace(source, argumentEnd,
            "__global const int* restrict rebuildNeighborList, unsigned int maxSinglePairs, GLOBAL int2* singlePairs\n#ifdef USE_LARGE_BLOCKS");
        replace(source, "bool interacts = false;", "bool interacts = false;\n                    uint pairMask = 0;");
        replace(source, "interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED);",
            "if (x*TILE_SIZE+j < NUM_ATOMS && dot(delta, delta) < PADDED_CUTOFF_SQUARED) pairMask |= 1u<<j;");
        const std::string compact = "                    // Do a prefix sum to compact the list of atoms.";
        replace(source, compact, R"(
                    interacts = pairMask != 0;
                    uint sparseCount = popcount(pairMask);
                    sparseCount = sparseCount <= 3 ? sparseCount : 0;
                    uint sparseTotal = simdReduceAdd(sparseCount);
                    uint sparseRank = simdPrefixExclusiveAdd(sparseCount);
                    uint sparseFirst = 0;
                    if (sparseTotal != 0 && indexInWarp == 0)
                        sparseFirst = ATOMIC_ADD(interactionCount+1, sparseTotal);
                    sparseFirst = simdShuffle(sparseFirst, 0);
                    if (sparseCount != 0) {
                        uint first = sparseFirst+sparseRank;
                        while (pairMask != 0) {
                            uint lane = ctz(pairMask);
                            pairMask &= pairMask-1;
                            if (first < maxSinglePairs)
                                singlePairs[first] = make_int2(x*TILE_SIZE+lane, atom2);
                            first++;
                        }
                        interacts = false;
                    }
)"+compact);
    }
#if OPENMM_METAL_FAST_NEIGHBOR_BALLOT
    replace(source, "    __local bool includeBlockFlags[GROUP_SIZE];", "");
    replace(source, "    __local volatile short2 atomCountBuffer[GROUP_SIZE];", "");
    const std::string loopStart =
        "            includeBlockFlags[get_local_id(0)] = includeBlock2;\n"
        "            SYNC_WARPS;\n"
        "            for (int i = 0; i < TILE_SIZE; i++) {\n"
        "                while (i < TILE_SIZE && !includeBlockFlags[warpStart+i])\n"
        "                    i++;\n"
        "                if (i < TILE_SIZE) {";
    replace(source, loopStart,
        "            uint includeBlockMask = simdBallot(includeBlock2);\n"
        "            while (includeBlockMask != 0) {\n"
        "                int i = ctz(includeBlockMask);\n"
        "                includeBlockMask &= includeBlockMask-1;\n                {");
    replace(source, "                else {\n                    SYNC_WARPS;\n                }", "");
    size_t first = source.find("                    atomCountBuffer[get_local_id(0)].x = (interacts ? 1 : 0);");
    size_t last = source.find("                    if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE)", first);
    if (first == std::string::npos || last == std::string::npos)
        throw OpenMMException("OpenCL neighbor compaction template changed");
    source.replace(first, last-first, R"(
                    uint atomMask = simdBallot(interacts);
                    uint lowerLanes = (1u<<uint(indexInWarp))-1u;
                    if (interacts)
                        buffer[neighborsInBuffer+popcount(atomMask&lowerLanes)] = atom2;
                    neighborsInBuffer += popcount(atomMask);
                    SYNC_WARPS;
)");
#endif
    return source;
}

} // namespace MetalNonbondedSources
} // namespace OpenMM
#endif
