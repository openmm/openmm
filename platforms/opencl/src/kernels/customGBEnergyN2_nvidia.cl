#define TILE_SIZE 32
#define STORE_DERIVATIVE_1(INDEX) derivBuffers##INDEX[offset1] += deriv##INDEX##_1;
#define STORE_DERIVATIVE_2(INDEX) derivBuffers##INDEX[offset2] += local_deriv##INDEX[get_local_id(0)];

/**
 * Compute a force based on pair interactions.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeN2Energy(__global float4* forceBuffers, __global float* energyBuffer, __local float4* local_force,
	__global float4* posq, __local float4* local_posq, __global unsigned int* exclusions, __global unsigned int* exclusionIndices,
        __global unsigned int* exclusionRowIndices, __local float4* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize
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
    __local unsigned int exclusionRange[4];
    __local int exclusionIndex[2];

    while (pos < end) {
        // Extract the coordinates of this tile
#ifdef USE_CUTOFF
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff);
        x = (x>>17);
#else
        unsigned int y = (unsigned int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        unsigned int x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y++;
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }
#endif
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        float4 force = 0.0f;
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

            local_posq[get_local_id(0)] = posq1;
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
                float dEdR = 0.0f;
                float tempEnergy = 0.0f;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
                    COMPUTE_INTERACTION
                    dEdR /= -r;
                }
                energy += 0.5f*tempEnergy;
                delta.xyz *= dEdR;
                force.xyz -= delta.xyz;
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            forceBuffers[offset1].xyz += force.xyz;
            STORE_DERIVATIVES_1
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y*TILE_SIZE + tgx;
                local_posq[get_local_id(0)] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_force[get_local_id(0)] = 0.0f;
            CLEAR_LOCAL_DERIVATIVES
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
            if (!hasExclusions && flags == 0) {
                // No interactions in this tile.
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
                    float dEdR = 0.0f;
                    float tempEnergy = 0.0f;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        COMPUTE_INTERACTION
                        dEdR /= -r;
                    }
		    energy += tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                    atom2 = tbx+tj;
                    local_force[atom2].xyz += delta.xyz;
                    RECORD_DERIVATIVE_2
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
            forceBuffers[offset1].xyz += force.xyz;
            forceBuffers[offset2].xyz += local_force[get_local_id(0)].xyz;
            STORE_DERIVATIVES_1
            STORE_DERIVATIVES_2
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
