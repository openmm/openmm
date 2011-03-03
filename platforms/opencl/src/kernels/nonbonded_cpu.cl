#define TILE_SIZE 32

typedef struct {
    float x, y, z;
    float q;
    float fx, fy, fz;
    ATOM_PARAMETER_DATA
} AtomData;

/**
 * Compute nonbonded interactions.
 */

__kernel void computeNonbonded(__global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, __global unsigned int* exclusionRowIndices, __local AtomData* localData, __local float4* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles, __global unsigned int* interactionFlags
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = get_group_id(0)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0);
    unsigned int end = (get_group_id(0)+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0);
#else
    unsigned int pos = get_group_id(0)*numTiles/get_num_groups(0);
    unsigned int end = (get_group_id(0)+1)*numTiles/get_num_groups(0);
#endif
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;

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
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
        }

        // Locate the exclusion data for this tile.

#ifdef USE_EXCLUSIONS
        unsigned int exclusionStart = exclusionRowIndices[x];
        unsigned int exclusionEnd = exclusionRowIndices[x+1];
        int exclusionIndex = -1;
        for (int i = exclusionStart; i < exclusionEnd; i++)
            if (exclusionIndices[i] == y) {
                exclusionIndex = i*TILE_SIZE;
                break;
            }
        bool hasExclusions = (exclusionIndex > -1);
#endif

        // Load the data for this tile if we don't already have it cached.

        if (lasty != y) {
            for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
                unsigned int j = y*TILE_SIZE + localAtomIndex;
                float4 tempPosq = posq[j];
                localData[localAtomIndex].x = tempPosq.x;
                localData[localAtomIndex].y = tempPosq.y;
                localData[localAtomIndex].z = tempPosq.z;
                localData[localAtomIndex].q = tempPosq.w;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
        }
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[exclusionIndex+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                float4 force = 0.0f;
                float4 posq1 = posq[atom1];
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    float4 posq2 = (float4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
                    float r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                    float invR = RSQRT(r2);
                    float r = RECIP(invR);
                    unsigned int atom2 = j;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                    float dEdR = 0.0f;
#else
                    float4 dEdR1 = (float4) 0.0f;
                    float4 dEdR2 = (float4) 0.0f;
#endif
                    float tempEnergy = 0.0f;
                    COMPUTE_INTERACTION
                    energy += 0.5f*tempEnergy;
#ifdef USE_SYMMETRIC
                    force.xyz -= delta.xyz*dEdR;
#else
                    force.xyz -= dEdR1.xyz;
#endif
#ifdef USE_CUTOFF
                    }
#endif
                    excl >>= 1;
                }

                // Write results.

                unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
                localData[tgx].fx = 0.0f;
                localData[tgx].fy = 0.0f;
                localData[tgx].fz = 0.0f;
            }
#ifdef USE_CUTOFF
            unsigned int flags1 = (numTiles <= maxTiles ? interactionFlags[2*pos] : 0xFFFFFFFF);
            unsigned int flags2 = (numTiles <= maxTiles ? interactionFlags[2*pos+1] : 0xFFFFFFFF);
            if (!hasExclusions && (flags1 != 0xFFFFFFFF || flags2 != 0xFFFFFFFF)) {
                // Compute only a subset of the interactions in this tile.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    if ((flags2&(1<<tgx)) != 0) {
                        unsigned int atom1 = x*TILE_SIZE+tgx;
                        float4 force = 0.0f;
                        float4 posq1 = posq[atom1];
                        LOAD_ATOM1_PARAMETERS
                        for (unsigned int j = 0; j < TILE_SIZE; j++) {
                            if ((flags1&(1<<j)) != 0) {
                                bool isExcluded = false;
                                float4 posq2 = (float4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                                delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
                                float r2 = dot(delta.xyz, delta.xyz);
                                if (r2 < CUTOFF_SQUARED) {
                                    float invR = RSQRT(r2);
                                    float r = RECIP(invR);
                                    unsigned int atom2 = j;
                                    LOAD_ATOM2_PARAMETERS
                                    atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                                    float dEdR = 0.0f;
#else
                                    float4 dEdR1 = (float4) 0.0f;
                                    float4 dEdR2 = (float4) 0.0f;
#endif
                                    float tempEnergy = 0.0f;
                                    COMPUTE_INTERACTION
                                    energy += tempEnergy;
#ifdef USE_SYMMETRIC
                                    delta.xyz *= dEdR;
                                    force.xyz -= delta.xyz;
                                    localData[j].fx += delta.x;
                                    localData[j].fy += delta.y;
                                    localData[j].fz += delta.z;
#else
                                    force.xyz -= dEdR1.xyz;
                                    localData[j].fx += dEdR2.x;
                                    localData[j].fy += dEdR2.y;
                                    localData[j].fz += dEdR2.z;
#endif
                                }
                            }
                        }

                        // Write results for atom1.

                        unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
                    }
                }
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    float4 force = 0.0f;
                    float4 posq1 = posq[atom1];
                    LOAD_ATOM1_PARAMETERS
#ifdef USE_EXCLUSIONS
                    unsigned int excl = (hasExclusions ? exclusions[exclusionIndex+tgx] : 0xFFFFFFFF);
#endif
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
#ifdef USE_EXCLUSIONS
                        bool isExcluded = !(excl & 0x1);
#endif
                        float4 posq2 = (float4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                        float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                        delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;
#endif
                        float r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
                        float invR = RSQRT(r2);
                        float r = RECIP(invR);
                        unsigned int atom2 = j;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                        float dEdR = 0.0f;
#else
                        float4 dEdR1 = (float4) 0.0f;
                        float4 dEdR2 = (float4) 0.0f;
#endif
                        float tempEnergy = 0.0f;
                        COMPUTE_INTERACTION
                        energy += tempEnergy;
#ifdef USE_SYMMETRIC
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        localData[j].fx += delta.x;
                        localData[j].fy += delta.y;
                        localData[j].fz += delta.z;
#else
                        force.xyz -= dEdR1.xyz;
                        localData[j].fx += dEdR2.x;
                        localData[j].fy += dEdR2.y;
                        localData[j].fz += dEdR2.z;
#endif
#ifdef USE_CUTOFF
                        }
#endif
#ifdef USE_EXCLUSIONS
                        excl >>= 1;
#endif
                    }

                   // Write results for atom1.

                    unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
                }
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int offset = y*TILE_SIZE+tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                float4 f = forceBuffers[offset];
                f.x += localData[tgx].fx;
                f.y += localData[tgx].fy;
                f.z += localData[tgx].fz;
                forceBuffers[offset] = f;
            }
        }
        lasty = y;
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
