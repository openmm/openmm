#define TILE_SIZE 32

/**
 * Compute a value based on pair interactions.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeN2Value(__global float4* posq, __local float4* local_posq, __global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, __global unsigned int* exclusionRowIndices, __global float* global_value, __local float* local_value,
        __local float* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionFlags, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
#endif
    unsigned int pos = get_group_id(0)*numTiles/get_num_groups(0);
    unsigned int end = (get_group_id(0)+1)*numTiles/get_num_groups(0);
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;
    __local unsigned int exclusionRange[2];
    __local int exclusionIndex[1];

    while (pos < end) {
        // Extract the coordinates of this tile
#ifdef USE_CUTOFF
        ushort2 tileIndices = tiles[pos];
        unsigned int x = tileIndices.x;
        unsigned int y = tileIndices.y;
#else
        unsigned int y = (unsigned int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        unsigned int x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y++;
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }
#endif
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int valueBufferOffset = (tgx < TILE_SIZE/2 ? 0 : TILE_SIZE);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        float value = 0.0f;
        float4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS

        // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
        if (get_local_id(0) < 2)
            exclusionRange[get_local_id(0)] = exclusionRowIndices[x+get_local_id(0)];
        if (tgx == 0)
            exclusionIndex[0] = -1;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i = exclusionRange[0]+tgx; i < exclusionRange[1]; i += TILE_SIZE)
            if (exclusionIndices[i] == y)
                exclusionIndex[0] = i*TILE_SIZE;
        barrier(CLK_LOCAL_MEM_FENCE);
        bool hasExclusions = (exclusionIndex[0] > -1);
#endif
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
            barrier(CLK_LOCAL_MEM_FENCE);
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndex[0]+tgx] >> baseLocalAtom;
#endif
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                int atom2 = baseLocalAtom+j;
                float4 posq2 = local_posq[atom2];
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                float r = SQRT(r2);
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+baseLocalAtom+j;
                float tempValue1 = 0.0f;
                float tempValue2 = 0.0f;
#ifdef USE_EXCLUSIONS
                if (!isExcluded && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#else
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
                    COMPUTE_VALUE
                }
                value += tempValue1;
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
            }

            // Sum the values and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = value;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                global_value[offset] += value+tempBuffer[get_local_id(0)+TILE_SIZE];
            }
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
                unsigned int j = y*TILE_SIZE + tgx;
                local_posq[get_local_id(0)] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_value[get_local_id(0)] = 0.0f;
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

#ifdef USE_EXCLUSIONS
            unsigned int excl = (hasExclusions ? exclusions[exclusionIndex[0]+tgx] : 0xFFFFFFFF);
            excl = (excl >> baseLocalAtom) & 0xFFFF;
            excl += excl << 16;
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
            unsigned int tj = tgx%(TILE_SIZE/2);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                int atom2 = baseLocalAtom+tj;
                float4 posq2 = local_posq[atom2];
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                float r = SQRT(r2);
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+baseLocalAtom+tj;
                float tempValue1 = 0.0f;
                float tempValue2 = 0.0f;
#ifdef USE_EXCLUSIONS
                if (!isExcluded && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#else
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                    COMPUTE_VALUE
                }
                value += tempValue1;
                local_value[baseLocalAtom+tj+valueBufferOffset] += tempValue2;
#ifdef USE_CUTOFF
                }
#endif
                barrier(CLK_LOCAL_MEM_FENCE);
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                tj = (tj+1)%(TILE_SIZE/2);
            }

            // Sum the values and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = value;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                global_value[offset1] += value+tempBuffer[get_local_id(0)+TILE_SIZE];
                global_value[offset2] += local_value[get_local_id(0)]+local_value[get_local_id(0)+TILE_SIZE];
            }
        }
        lasty = y;
        pos++;
    }
}
