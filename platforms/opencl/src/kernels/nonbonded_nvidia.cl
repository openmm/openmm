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

__attribute__((reqd_work_group_size(64, 1, 1)))
__kernel void computeNonbonded(__global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, __local AtomData* localData, __local float4* tempBuffer, __global unsigned int* tiles,
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
        float4 force = 0.0f;
        float4 posq1 = posq[atom1];
        __local AtomData* atom1Data = &localData[get_local_id(0)];
        LOAD_ATOM1_PARAMETERS
        if (x == y) {
            // This tile is on the diagonal.

            atom1Data->x = posq1.x;
            atom1Data->y = posq1.y;
            atom1Data->z = posq1.z;
            atom1Data->q = posq1.w;
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
                __local AtomData* atom2Data = &localData[tbx+j];
                float4 posq2 = (float4) (atom2Data->x, atom2Data->y, atom2Data->z, atom2Data->q);
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*INV_PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                delta.y -= floor(delta.y*INV_PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                delta.z -= floor(delta.z*INV_PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float r = sqrt(r2);
                float invR = 1.0f/r;
                LOAD_ATOM2_PARAMETERS
                int atom2 = y+j;
                float dEdR = 0.0f;
                float tempEnergy = 0.0f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
                delta.xyz *= dEdR;
                force.xyz -= delta.xyz;
                excl >>= 1;
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x + tgx + warp*PADDED_NUM_ATOMS;
#endif
            forceBuffers[offset].xyz += force.xyz;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y + tgx;
                float4 tempPosq = posq[j];
                atom1Data->x = tempPosq.x;
                atom1Data->y = tempPosq.y;
                atom1Data->z = tempPosq.z;
                atom1Data->q = tempPosq.w;
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            atom1Data->fx = 0.0f;
            atom1Data->fy = 0.0f;
            atom1Data->fz = 0.0f;
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
                            bool isExcluded = false;
                            __local AtomData* atom2Data = &localData[tbx+j];
                            float4 posq2 = (float4) (atom2Data->x, atom2Data->y, atom2Data->z, atom2Data->q);
                            float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x*INV_PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                            delta.y -= floor(delta.y*INV_PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                            delta.z -= floor(delta.z*INV_PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float invR = RSQRT(r2);
                            float r = RECIP(invR);
                            LOAD_ATOM2_PARAMETERS
                            int atom2 = y+j;
                            float dEdR = 0.0f;
                            float tempEnergy = 0.0f;
                            COMPUTE_INTERACTION
			    energy += tempEnergy;
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            tempBuffer[get_local_id(0)] = delta;

                            // Sum the forces on atom2.

                            if (tgx % 2 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+1].xyz;
                            if (tgx % 4 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+2].xyz;
                            if (tgx % 8 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+4].xyz;
                            if (tgx % 16 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+8].xyz;
                            if (tgx == 0) {
                                atom2Data->fx += tempBuffer[get_local_id(0)].x + tempBuffer[get_local_id(0)+16].x;
                                atom2Data->fy += tempBuffer[get_local_id(0)].y + tempBuffer[get_local_id(0)+16].y;
                                atom2Data->fz += tempBuffer[get_local_id(0)].z + tempBuffer[get_local_id(0)+16].z;
                            }
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
                    __local AtomData* atom2Data = &localData[tbx+tj];
                    float4 posq2 = (float4) (atom2Data->x, atom2Data->y, atom2Data->z, atom2Data->q);
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*INV_PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                    delta.y -= floor(delta.y*INV_PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                    delta.z -= floor(delta.z*INV_PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float invR = RSQRT(r2);
                    float r = RECIP(invR);
                    LOAD_ATOM2_PARAMETERS
                    int atom2 = y+tj;
                    float dEdR = 0.0f;
                    float tempEnergy = 0.0f;
                    COMPUTE_INTERACTION
		    energy += tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                    atom2Data->fx += delta.x;
                    atom2Data->fy += delta.y;
                    atom2Data->fz += delta.z;
                    excl >>= 1;
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
            forceBuffers[offset1].xyz += force.xyz;
            forceBuffers[offset2] += (float4) (atom1Data->fx, atom1Data->fy, atom1Data->fz, 0.0f);
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
