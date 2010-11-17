#define TILE_SIZE 32

/**
 * Compute a value based on pair interactions.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeN2Value(__global float4* posq, __local float4* local_posq, __global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, __global unsigned int* exclusionRowIndices, __global float* global_value, __local float* local_value,
        __local float* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles, __global unsigned int* interactionFlags
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;
    __local unsigned int exclusionRange[4];
    __local int exclusionIndex[2];

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
            y = (unsigned int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y++;
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
        }
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        float value = 0.0f;
        float4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS

        // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
        int localGroupIndex = get_local_id(0)/TILE_SIZE;
        if (tgx < 2)
            exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
        if (tgx == 0)
            exclusionIndex[localGroupIndex] = -1;
        for (int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
            if (exclusionIndices[i] == y)
                exclusionIndex[localGroupIndex] = i*TILE_SIZE;
        bool hasExclusions = (exclusionIndex[localGroupIndex] > -1);
#else
        bool hasExclusions = false;
#endif
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = get_local_id(0);
            local_posq[localAtomIndex] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndex[localGroupIndex]+tgx];
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                int atom2 = tbx+j;
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
                atom2 = y*TILE_SIZE+j;
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

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_value[offset] += value;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y*TILE_SIZE + tgx;
                local_posq[get_local_id(0)] = posq[j];
                const unsigned int localAtomIndex = get_local_id(0);
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_value[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
            unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
            if (!hasExclusions && flags != 0xFFFFFFFF) {
                if (flags == 0) {
                    // No interactions in this tile.
                }
                else {
                    // Compute only a subset of the interactions in this tile.

                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        if ((flags&(1<<j)) != 0) {
                            int atom2 = tbx+j;
                            float4 posq2 = local_posq[atom2];
                            float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float tempValue1 = 0.0f;
                            float tempValue2 = 0.0f;
                            if (r2 < CUTOFF_SQUARED) {
                                float r = SQRT(r2);
                                LOAD_ATOM2_PARAMETERS
                                atom2 = y*TILE_SIZE+j;
                                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                                    COMPUTE_VALUE
                                }
                                value += tempValue1;
                            }
                            tempBuffer[get_local_id(0)] = tempValue2;

                            // Sum the forces on atom2.

                            if (tgx % 2 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+1];
                            if (tgx % 4 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+2];
                            if (tgx % 8 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+4];
                            if (tgx % 16 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+8];
                            if (tgx == 0)
                                local_value[tbx+j] += tempBuffer[get_local_id(0)] + tempBuffer[get_local_id(0)+16];
                        }
                    }
                }
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

#ifdef USE_EXCLUSIONS
                unsigned int excl = (hasExclusions ? exclusions[exclusionIndex[localGroupIndex]+tgx] : 0xFFFFFFFF);
                excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    int atom2 = tbx+tj;
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
                    atom2 = y*TILE_SIZE+tj;
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
                    local_value[tbx+tj] += tempValue2;
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
            unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_value[offset1] += value;
            global_value[offset2] += local_value[get_local_id(0)];
        }
        lasty = y;
        pos++;
    }
}
