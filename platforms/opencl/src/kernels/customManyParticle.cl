#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

/**
 * Record the force on an atom to global memory.
 */
inline void storeForce(int atom, real4 force, __global long* restrict forceBuffers) {
    atom_add(&forceBuffers[atom], (long) (force.x*0x100000000));
    atom_add(&forceBuffers[atom+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
    atom_add(&forceBuffers[atom+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
}

/**
 * Compute the difference between two vectors, taking periodic boundary conditions into account
 * and setting the fourth component to the squared magnitude.
 */
inline real4 delta(real4 vec1, real4 vec2, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = (real4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
#ifdef USE_PERIODIC
    APPLY_PERIODIC_TO_DELTA(result)
#endif
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
real computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real4 crossProduct = cross(vec1, vec2);
        real scale = vec1.w*vec2.w;
        angle = asin(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0.0f)
            angle = M_PI-angle;
    }
    else
       angle = acos(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
inline real4 computeCross(real4 vec1, real4 vec2) {
    real4 cp = cross(vec1, vec2);
    return (real4) (cp.x, cp.y, cp.z, cp.x*cp.x+cp.y*cp.y+cp.z*cp.z);
}

/**
 * Determine whether a particular interaction is in the list of exclusions.
 */
inline bool isInteractionExcluded(int atom1, int atom2, __global const int* restrict exclusions, __global const int* restrict exclusionStartIndex) {
    if (atom1 > atom2) {
        int temp = atom1;
        atom1 = atom2;
        atom2 = temp;
    }
    int first = exclusionStartIndex[atom1];
    int last = exclusionStartIndex[atom1+1];
    for (int i = last-1; i >= first; i--) {
        int excluded = exclusions[i];
        if (excluded == atom2)
            return true;
        if (excluded <= atom1)
            return false;
    }
    return false;
}

/**
 * Compute the interaction.
 */
__kernel void computeInteraction(
        __global long* restrict forceBuffers, __global mixed* restrict energyBuffer, __global const real4* restrict posq,
        real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ
#ifdef USE_CUTOFF
        , __global const int* restrict neighbors, __global const int* restrict neighborStartIndex
#endif
#ifdef USE_FILTERS
        , __global int* restrict particleTypes, __global int* restrict orderIndex, __global int* restrict particleOrder
#endif
#ifdef USE_EXCLUSIONS
        , __global int* restrict exclusions, __global int* restrict exclusionStartIndex
#endif
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    
    // Loop over particles to be the first one in the set.
    
    for (int p1 = get_group_id(0); p1 < NUM_ATOMS; p1 += get_num_groups(0)) {
#ifdef USE_CENTRAL_PARTICLE
        const int a1 = p1;
#else
        const int a1 = 0;
#endif
#ifdef USE_CUTOFF
        int firstNeighbor = neighborStartIndex[p1];
        int numNeighbors = neighborStartIndex[p1+1]-firstNeighbor;
#else
  #ifdef USE_CENTRAL_PARTICLE
        int numNeighbors = NUM_ATOMS;
  #else
        int numNeighbors = NUM_ATOMS-p1-1;
  #endif
#endif
        int numCombinations = NUM_CANDIDATE_COMBINATIONS;
        for (int index = get_local_id(0); index < numCombinations; index += get_local_size(0)) {
            FIND_ATOMS_FOR_COMBINATION_INDEX;
            bool includeInteraction = IS_VALID_COMBINATION;
#ifdef USE_CUTOFF
            if (includeInteraction) {
                VERIFY_CUTOFF;
            }
#endif
#ifdef USE_FILTERS
            int order = orderIndex[COMPUTE_TYPE_INDEX];
            if (order == -1)
                includeInteraction = false;
#endif
#ifdef USE_EXCLUSIONS
            if (includeInteraction) {
                VERIFY_EXCLUSIONS;
            }
#endif
            if (includeInteraction) {
                PERMUTE_ATOMS;
                LOAD_PARTICLE_DATA;
                COMPUTE_INTERACTION;
            }
        }
    }
    energyBuffer[get_global_id(0)] += energy;
}

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global const real4* restrict posq, __global real4* restrict blockCenter, __global real4* restrict blockBoundingBox, __global int* restrict numNeighborPairs) {
    int index = get_global_id(0);
    int base = index*TILE_SIZE;
    while (base < NUM_ATOMS) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, NUM_ATOMS);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = (real4) (min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = (real4) (max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += get_global_size(0);
        base = index*TILE_SIZE;
    }
    if (get_group_id(0) == 0 && get_local_id(0) == 0)
        *numNeighborPairs = 0;
}

/**
 * Find a list of neighbors for each atom.
 */
__kernel void findNeighbors(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global const real4* restrict posq, __global const real4* restrict blockCenter, __global const real4* restrict blockBoundingBox, __global int2* restrict neighborPairs,
        __global int* restrict numNeighborPairs, __global int* restrict numNeighborsForAtom, int maxNeighborPairs
#ifdef USE_EXCLUSIONS
        , __global const int* restrict exclusions, __global const int* restrict exclusionStartIndex
#endif
        ) {
    __local real4 positionCache[FIND_NEIGHBORS_WORKGROUP_SIZE];
    __local bool includeBlockFlags[FIND_NEIGHBORS_WORKGROUP_SIZE];
    int indexInWarp = get_local_id(0)%32;
    int warpStart = get_local_id(0)-indexInWarp;
    for (int atom1 = get_global_id(0); atom1 < PADDED_NUM_ATOMS; atom1 += get_global_size(0)) {
        // Load data for this atom.  Note that all threads in a warp are processing atoms from the same block.
        
        real4 pos1 = posq[atom1];
        int block1 = atom1/TILE_SIZE;
        real4 blockCenter1 = blockCenter[block1];
        real4 blockSize1 = blockBoundingBox[block1];
        int totalNeighborsForAtom1 = 0;
        
        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

#ifdef USE_CENTRAL_PARTICLE
        int startBlock = 0;
#else
        int startBlock = block1;
#endif
        for (int block2Base = startBlock; block2Base < NUM_BLOCKS; block2Base += 32) {
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            if (includeBlock2) {
                real4 blockCenter2 = blockCenter[block2];
                real4 blockSize2 = blockBoundingBox[block2];
                real4 blockDelta = blockCenter1-blockCenter2;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                blockDelta.x = max((real) 0, fabs(blockDelta.x)-blockSize1.x-blockSize2.x);
                blockDelta.y = max((real) 0, fabs(blockDelta.y)-blockSize1.y-blockSize2.y);
                blockDelta.z = max((real) 0, fabs(blockDelta.z)-blockSize1.z-blockSize2.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < CUTOFF_SQUARED);
            }
            
            // Loop over any blocks we identified as potentially containing neighbors.
            
            includeBlockFlags[get_local_id(0)] = includeBlock2;
            SYNC_WARPS;
            for (int i = 0; i < TILE_SIZE; i++) {
                if (includeBlockFlags[warpStart+i]) {
                    int block2 = block2Base+i;

                    // Loop over atoms in this block.

                    int start = block2*TILE_SIZE;
                    int included[TILE_SIZE];
                    int numIncluded = 0;
                    SYNC_WARPS;
                    positionCache[get_local_id(0)] = posq[start+indexInWarp];
                    SYNC_WARPS;
                    if (atom1 < NUM_ATOMS) {
                        for (int j = 0; j < 32; j++) {
                            int atom2 = start+j;
                            real4 pos2 = positionCache[get_local_id(0)-indexInWarp+j];

                            // Decide whether to include this atom pair in the neighbor list.

                            real4 atomDelta = delta(pos1, pos2, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
#ifdef USE_CENTRAL_PARTICLE
                            bool includeAtom = (atom2 != atom1 && atom2 < NUM_ATOMS && atomDelta.w < CUTOFF_SQUARED);
#else
                            bool includeAtom = (atom2 > atom1 && atom2 < NUM_ATOMS && atomDelta.w < CUTOFF_SQUARED);
#endif
#ifdef USE_EXCLUSIONS
                            if (includeAtom)
                                includeAtom &= !isInteractionExcluded(atom1, atom2, exclusions, exclusionStartIndex);
#endif
                            if (includeAtom)
                                included[numIncluded++] = atom2;
                        }
                    }

                    // If we found any neighbors, store them to the neighbor list.

                    if (numIncluded > 0) {
                        int baseIndex = atom_add(numNeighborPairs, numIncluded);
                        if (baseIndex+numIncluded <= maxNeighborPairs)
                            for (int j = 0; j < numIncluded; j++)
                                neighborPairs[baseIndex+j] = (int2) (atom1, included[j]);
                        totalNeighborsForAtom1 += numIncluded;
                    }
                }
            }
        }
        if (atom1 < NUM_ATOMS)
            numNeighborsForAtom[atom1] = totalNeighborsForAtom1;
        SYNC_WARPS;
    }
}

/**
 * Sum the neighbor counts to compute the start position of each atom.  This kernel
 * is executed as a single work group.
 */
__kernel void computeNeighborStartIndices(__global int* restrict numNeighborsForAtom, __global int* restrict neighborStartIndex,
            __global int* restrict numNeighborPairs, int maxNeighborPairs) {
    __local unsigned int posBuffer[256];
    if (*numNeighborPairs > maxNeighborPairs) {
        // There wasn't enough memory for the neighbor list, so we'll need to rebuild it.  Set the neighbor start
        // indices to indicate no neighbors for any atom.
        
        for (int i = get_local_id(0); i <= NUM_ATOMS; i += get_local_size(0))
            neighborStartIndex[i] = 0;
        return;
    }
    unsigned int globalOffset = 0;
    for (unsigned int startAtom = 0; startAtom < NUM_ATOMS; startAtom += get_local_size(0)) {
        // Load the neighbor counts into local memory.

        unsigned int globalIndex = startAtom+get_local_id(0);
        posBuffer[get_local_id(0)] = (globalIndex < NUM_ATOMS ? numNeighborsForAtom[globalIndex] : 0);
        barrier(CLK_LOCAL_MEM_FENCE);

        // Perform a parallel prefix sum.

        for (unsigned int step = 1; step < get_local_size(0); step *= 2) {
            unsigned int add = (get_local_id(0) >= step ? posBuffer[get_local_id(0)-step] : 0);
            barrier(CLK_LOCAL_MEM_FENCE);
            posBuffer[get_local_id(0)] += add;
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // Write the results back to global memory.

        if (globalIndex < NUM_ATOMS) {
            neighborStartIndex[globalIndex+1] = posBuffer[get_local_id(0)]+globalOffset;
            numNeighborsForAtom[globalIndex] = 0; // Clear this so the next kernel can use it as a counter
        }
        globalOffset += posBuffer[get_local_size(0)-1];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0)
        neighborStartIndex[0] = 0;
}

/**
 * Assemble the final neighbor list.
 */
__kernel void copyPairsToNeighborList(__global const int2* restrict neighborPairs, __global int* restrict neighbors, __global int* restrict numNeighborPairs,
            int maxNeighborPairs, __global int* restrict numNeighborsForAtom, __global const int* restrict neighborStartIndex) {
    int actualPairs = *numNeighborPairs;
    if (actualPairs > maxNeighborPairs)
        return; // There wasn't enough memory for the neighbor list, so we'll need to rebuild it.
    for (unsigned int index = get_global_id(0); index < actualPairs; index += get_global_size(0)) {
        int2 pair = neighborPairs[index];
        int startIndex = neighborStartIndex[pair.x];
        int offset = atom_add(numNeighborsForAtom+pair.x, 1);
        neighbors[startIndex+offset] = pair.y;
    }
}
