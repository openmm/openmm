#define TILE_SIZE 32

/**
 * Compute a value based on pair interactions.
 */

__kernel void computeN2Value(__global float4* posq, __local float4* local_posq, __global float* global_value,
        __local float* local_value, __local float* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
#endif
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;

    while (pos < end) {
        // Extract the coordinates of this tile
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff)*TILE_SIZE;
        bool hasExclusions = (x & 0x1);
        x = (x>>17)*TILE_SIZE;
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int atom1 = x + tgx;
        float value = 0.0f;
        float4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
            unsigned int xi = x/TILE_SIZE;
            unsigned int tile = xi+xi*PADDED_NUM_ATOMS/TILE_SIZE-xi*(xi+1)/2;
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndices[tile]+tgx];
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                int atom2 = tbx+j;
                float4 posq2 = local_posq[atom2];
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float r = sqrt(r2);
                LOAD_ATOM2_PARAMETERS
                atom2 = y+j;
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
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_value[offset] += value;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y + tgx;
                local_posq[get_local_id(0)] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_value[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
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
                            delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                            delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                            delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float r = sqrt(r2);
                            LOAD_ATOM2_PARAMETERS
                            atom2 = y+j;
                            float tempValue1 = 0.0f;
                            float tempValue2 = 0.0f;
                            if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                                COMPUTE_VALUE
                            }
                            value += tempValue1;
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

                unsigned int xi = x/TILE_SIZE;
                unsigned int yi = y/TILE_SIZE;
                unsigned int tile = xi+yi*PADDED_NUM_ATOMS/TILE_SIZE-yi*(yi+1)/2;
#ifdef USE_EXCLUSIONS
                unsigned int excl = (hasExclusions ? exclusions[exclusionIndices[tile]+tgx] : 0xFFFFFFFF);
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
                    delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                    delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                    delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float r = sqrt(r2);
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y+tj;
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
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x + tgx + (y/TILE_SIZE)*PADDED_NUM_ATOMS;
            unsigned int offset2 = y + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset1 = x + tgx + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = y + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_value[offset1] += value;
            global_value[offset2] += local_value[get_local_id(0)];
            lasty = y;
        }
        pos++;
    }
}
