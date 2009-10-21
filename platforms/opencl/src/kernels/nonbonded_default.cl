const unsigned int TileSize = 32;

/**
 * Compute nonbonded interactions.
 */

__kernel void computeNonbonded(__global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global unsigned int* exclusions,
        __global unsigned int* exclusionIndices, __local float4* local_posq, __local float4* local_force, __local float4* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount
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

    while (pos < end) {
        // Extract the coordinates of this tile
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff)*TileSize;
        bool hasExclusions = (x & 0x1);
        x = (x>>17)*TileSize;
        unsigned int baseLocalAtom = (get_local_id(0) < TileSize ? 0 : TileSize/2);
        unsigned int tgx = get_local_id(0) & (TileSize-1);
        unsigned int forceBufferOffset = (tgx < TileSize/2 ? 0 : TileSize);
        unsigned int atom1 = x + tgx;
        float4 force = 0.0f;
        float4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
            barrier(CLK_LOCAL_MEM_FENCE);
            unsigned int xi = x/TileSize;
            unsigned int tile = xi+xi*PADDED_NUM_ATOMS/TileSize-xi*(xi+1)/2;
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndices[tile]+tgx] >> baseLocalAtom;
#endif
            for (unsigned int j = 0; j < TileSize/2; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                int atom2 = baseLocalAtom+j;
                float4 posq2 = local_posq[atom2];
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float r = sqrt(r2);
                float invR = 1.0f/r;
                LOAD_ATOM2_PARAMETERS
                atom2 = y+baseLocalAtom+j;
                float dEdR = 0.0f;
                float tempEnergy = 0.0f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
                delta.xyz *= dEdR;
                force.xyz -= delta.xyz;
                excl >>= 1;
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TileSize)
                tempBuffer[get_local_id(0)] = force;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TileSize) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset = x + tgx + (x/TileSize)*PADDED_NUM_ATOMS;
#else
                unsigned int offset = x + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                forceBuffers[offset].xyz += force.xyz+tempBuffer[get_local_id(0)+TileSize].xyz;
            }
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TileSize) {
                unsigned int j = y + tgx;
                local_posq[get_local_id(0)] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_force[get_local_id(0)] = 0.0f;
            barrier(CLK_LOCAL_MEM_FENCE);
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
            if (!hasExclusions && flags == 0) {
                // No interactions in this tile.
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

                unsigned int xi = x/TileSize;
                unsigned int yi = y/TileSize;
                unsigned int tile = xi+yi*PADDED_NUM_ATOMS/TileSize-yi*(yi+1)/2;
#ifdef USE_EXCLUSIONS
                unsigned int excl = (hasExclusions ? exclusions[exclusionIndices[tile]+tgx] : 0xFFFFFFFF);
                excl = (excl >> tgx) | (excl << (TileSize - tgx));
                excl >>= baseLocalAtom;
#endif
                unsigned int tj = tgx%(TileSize/2);
                for (unsigned int j = 0; j < TileSize/2; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    int atom2 = baseLocalAtom+tj;
                    float4 posq2 = local_posq[atom2];
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                    delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                    delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float r = sqrt(r2);
                    float invR = 1.0f/r;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y+baseLocalAtom+tj;
                    float dEdR = 0.0f;
                    float tempEnergy = 0.0f;
                    COMPUTE_INTERACTION
		    energy += tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                    local_force[baseLocalAtom+tj+forceBufferOffset].xyz += delta.xyz;
                    barrier(CLK_LOCAL_MEM_FENCE);
                    excl >>= 1;
                    tj = (tj+1)%(TileSize/2);
                }
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TileSize)
                tempBuffer[get_local_id(0)] = force;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TileSize) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset1 = x + tgx + (y/TileSize)*PADDED_NUM_ATOMS;
                unsigned int offset2 = y + tgx + (x/TileSize)*PADDED_NUM_ATOMS;
#else
                unsigned int offset1 = x + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                unsigned int offset2 = y + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                forceBuffers[offset1].xyz += force.xyz+tempBuffer[get_local_id(0)+TileSize].xyz;
                forceBuffers[offset2].xyz += local_force[get_local_id(0)].xyz+local_force[get_local_id(0)+TileSize].xyz;
            }
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
