#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif

#define TILE_SIZE 32

// Cannot use float3 as OpenCL defines it to be 4 DWORD aligned. This would
// cause every element of array to have DWORD of padding to make it 4 DWORD
// aligned which wastes space and causes LDS bank conflicts as stride is no
// longer odd DWORDS.
typedef struct {
    real x, y, z;
} UnalignedReal3;

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;

/**
 * Compute nonbonded interactions.
 */

__kernel __attribute__((reqd_work_group_size(FORCE_WORK_GROUP_SIZE, 1, 1)))
void computeNonbonded(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global real* restrict energyBuffer, __global const real4* restrict posq, __global const unsigned int* restrict exclusions,
        __global const unsigned int* restrict exclusionIndices, __global const unsigned int* restrict exclusionRowIndices,
        unsigned int startTileIndex, unsigned int endTileIndex,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, __global const unsigned int* restrict interactionFlags
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = (numTiles > maxTiles ? startTileIndex+get_group_id(0)*(endTileIndex-startTileIndex)/get_num_groups(0) : get_group_id(0)*numTiles/get_num_groups(0));
    unsigned int end = (numTiles > maxTiles ? startTileIndex+(get_group_id(0)+1)*(endTileIndex-startTileIndex)/get_num_groups(0) : (get_group_id(0)+1)*numTiles/get_num_groups(0));
#else
    unsigned int pos = startTileIndex+get_group_id(0)*numTiles/get_num_groups(0);
    unsigned int end = startTileIndex+(get_group_id(0)+1)*numTiles/get_num_groups(0);
#endif
    real energy = 0;
    unsigned int lasty = 0xFFFFFFFF;
    __local AtomData localData[TILE_SIZE];
    __local UnalignedReal3 localForce[FORCE_WORK_GROUP_SIZE];
#ifdef USE_EXCLUSIONS
    __local unsigned int exclusionRange[2];
    __local int exclusionIndex[1];
#endif

    while (pos < end) {
        // Extract the coordinates of this tile
        unsigned int x, y;
#ifdef USE_CUTOFF
        if (numTiles <= maxTiles) {
            ushort2 tileIndices = tiles[pos];
            x = tileIndices.x;
            y = tileIndices.y;
        }
        else
#endif
        {
            y = (unsigned int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
        }
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int localForceOffset = get_local_id(0) & ~(TILE_SIZE-1);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 force = 0;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS

        // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
        if (get_local_id(0) < 2)
            exclusionRange[get_local_id(0)] = exclusionRowIndices[x+get_local_id(0)];
        if (get_local_id(0) == 0)
            exclusionIndex[0] = -1;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i = exclusionRange[0]+get_local_id(0); i < exclusionRange[1]; i += FORCE_WORK_GROUP_SIZE)
            if (exclusionIndices[i] == y)
                exclusionIndex[0] = i*TILE_SIZE;
        barrier(CLK_LOCAL_MEM_FENCE);
        bool hasExclusions = (exclusionIndex[0] > -1);
#endif
        if (x == y) {
            // This tile is on the diagonal.

            if (get_local_id(0) < TILE_SIZE) {
                const unsigned int localAtomIndex = tgx;
                localData[localAtomIndex].x = posq1.x;
                localData[localAtomIndex].y = posq1.y;
                localData[localAtomIndex].z = posq1.z;
                localData[localAtomIndex].q = posq1.w;
                LOAD_LOCAL_PARAMETERS_FROM_1
            }
            barrier(CLK_LOCAL_MEM_FENCE);
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndex[0]+tgx] >> baseLocalAtom;
#endif
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                unsigned int atom2 = baseLocalAtom+j;
                real4 posq2 = (real4) (localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+baseLocalAtom+j;
#ifdef USE_SYMMETRIC
                real dEdR = 0;
#else
                real4 dEdR1 = (real4) 0;
                real4 dEdR2 = (real4) 0;
#endif
                real tempEnergy = 0;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
#ifdef USE_SYMMETRIC
                force.xyz -= delta.xyz*dEdR;
#else
                force.xyz -= dEdR1.xyz;
#endif
                excl >>= 1;
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE) {
                localData[tgx].fx = force.x;
                localData[tgx].fy = force.y;
                localData[tgx].fz = force.z;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                const unsigned int offset = x*TILE_SIZE + tgx;
                atom_add(&forceBuffers[offset], (long) ((force.x + localData[tgx].fx)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) ((force.y + localData[tgx].fy)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) ((force.z + localData[tgx].fz)*0xFFFFFFFF));
#else
                force.x += localData[tgx].fx;
                force.y += localData[tgx].fy;
                force.z += localData[tgx].fz;
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                // Cheaper to load/store real4 than real3.
                real4 sum = forceBuffers[offset];
                sum.xyz += force.xyz;
                forceBuffers[offset] = sum;
#endif
            }
            // barrier not required here as localData[*].temp is not accessed before encountering another barrier.
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
                const unsigned int localAtomIndex = tgx;
                unsigned int j = y*TILE_SIZE + tgx;
                real4 tempPosq = posq[j];
                localData[localAtomIndex].x = tempPosq.x;
                localData[localAtomIndex].y = tempPosq.y;
                localData[localAtomIndex].z = tempPosq.z;
                localData[localAtomIndex].q = tempPosq.w;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            localForce[get_local_id(0)].x = 0;
            localForce[get_local_id(0)].y = 0;
            localForce[get_local_id(0)].z = 0;
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

            unsigned int tj = (tgx+baseLocalAtom) & (TILE_SIZE-1);
#ifdef USE_EXCLUSIONS
            unsigned int excl = (hasExclusions ? exclusions[exclusionIndex[0]+tgx] : 0xFFFFFFFF);
            excl = (excl >> tj) | (excl << (TILE_SIZE - tj));
#endif
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                real4 posq2 = (real4) (localData[tj].x, localData[tj].y, localData[tj].z, localData[tj].q);
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                int atom2 = tj;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+tj;
#ifdef USE_SYMMETRIC
                real dEdR = 0;
#else
                real4 dEdR1 = (real4) 0;
                real4 dEdR2 = (real4) 0;
#endif
                real tempEnergy = 0;
                COMPUTE_INTERACTION
                energy += tempEnergy;
#ifdef USE_SYMMETRIC
                delta.xyz *= dEdR;
                force.xyz -= delta.xyz;
                localForce[tj+localForceOffset].x += delta.x;
                localForce[tj+localForceOffset].y += delta.y;
                localForce[tj+localForceOffset].z += delta.z;
#else
                force.xyz -= dEdR1.xyz;
                localForce[tj+localForceOffset].x += dEdR2.x;
                localForce[tj+localForceOffset].y += dEdR2.y;
                localForce[tj+localForceOffset].z += dEdR2.z;
#endif
                barrier(CLK_LOCAL_MEM_FENCE);
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                tj = (tj+1) & (TILE_SIZE-1);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE) {
                localData[tgx].fx = force.x;
                localData[tgx].fy = force.y;
                localData[tgx].fz = force.z;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                const unsigned int offset1 = x*TILE_SIZE + tgx;
                const unsigned int offset2 = y*TILE_SIZE + tgx;
                atom_add(&forceBuffers[offset1], (long) ((force.x + localData[tgx].fx)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset1+PADDED_NUM_ATOMS], (long) ((force.y + localData[tgx].fy)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset1+2*PADDED_NUM_ATOMS], (long) ((force.z + localData[tgx].fz)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2], (long) ((localForce[tgx].x + localForce[tgx+TILE_SIZE].x)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2+PADDED_NUM_ATOMS], (long) ((localForce[tgx].y + localForce[tgx+TILE_SIZE].y)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2+2*PADDED_NUM_ATOMS], (long) ((localForce[tgx].z + localForce[tgx+TILE_SIZE].z)*0xFFFFFFFF));
#else
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                const unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                const unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                // Cheaper to load/store real4 than real3. Do all loads before all stores to minimize store-load waits.
                real4 sum1 = forceBuffers[offset1];
                real4 sum2 = forceBuffers[offset2];
                sum1.x += localData[tgx].fx + force.x;
                sum1.y += localData[tgx].fy + force.y;
                sum1.z += localData[tgx].fz + force.z;
                sum2.x += localForce[tgx].x + localForce[tgx+TILE_SIZE].x;
                sum2.y += localForce[tgx].y + localForce[tgx+TILE_SIZE].y;
                sum2.z += localForce[tgx].z + localForce[tgx+TILE_SIZE].z;
                forceBuffers[offset1] = sum1;
                forceBuffers[offset2] = sum2;
#endif
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        lasty = y;
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
