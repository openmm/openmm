const unsigned int TileSize = 32;

/**
 * Compute nonbonded interactions.
 */

__kernel void computeNonbonded(int paddedNumAtoms, __global float4* forceBuffers, __global float* energyBuffer, __global float4* posq,
        __global unsigned int* exclusions,  __global unsigned int* exclusionIndices, __local float4* local_posq, __local float4* local_force, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        float cutoffSquared, float4 periodicBoxSize, __global unsigned int* interactionFlags, __global unsigned int* interactionCount, __local float4* tempBuffer
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
#endif
    unsigned int totalWarps = get_global_size(0)/TileSize;
    unsigned int warp = get_global_id(0)/TileSize;
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;

    while (pos < end) {
        // Extract the coordinates of this tile
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff)*TileSize;
        bool hasExclusions = (x & 0x1);
        x = (x>>17)*TileSize;
        unsigned int tgx = get_local_id(0) & (TileSize-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int tj = tgx;
        unsigned int i = x + tgx;
        float4 force = 0.0f;
        float4 posq1 = posq[i];
        LOAD_ATOM1_PARAMETERS
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
            unsigned int xi = x/TileSize;
            unsigned int tile = xi+xi*paddedNumAtoms/TileSize-xi*(xi+1)/2;
#ifdef USE_EXCLUSIONS
            unsigned int excl = exclusions[exclusionIndices[tile]+tgx];
#endif
            for (unsigned int j = 0; j < TileSize; j++) {
#ifdef USE_EXCLUSIONS
                bool isExcluded = !(excl & 0x1);
#endif
                float4 delta = (float4) (local_posq[tbx+j].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float r = sqrt(r2);
                float invR = 1.0f / r;
                float4 posq2 = local_posq[tbx+j];
                LOAD_ATOM2_PARAMETERS_J
                float dEdR = 0.0f;
                float tempEnergy = 0.0f;
                COMPUTE_INTERACTION
#if defined USE_CUTOFF || defined USE_EXCLUSIONS
    #if defined USE_CUTOFF && defined USE_EXCLUSIONS
                excl >>= 1;
                if (isExcluded || r2 > cutoffSquared) {
    #elif defined USE_CUTOFF
                if (r2 > cutoffSquared) {
    #elif defined USE_EXCLUSIONS
                excl >>= 1;
                if (isExcluded) {
    #endif
                    dEdR = 0.0f;
                    tempEnergy  = 0.0f;
                }
#endif
                energy += 0.5f*tempEnergy;
                delta.xyz *= dEdR;
                force.xyz -= delta.xyz;
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x + tgx + (x/TileSize)*paddedNumAtoms;
#else
            unsigned int offset = x + tgx + warp*paddedNumAtoms;
#endif
            of = forceBuffers[offset];
            of.xyz += force.xyz;
            forceBuffers[offset] = of;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y + tgx;
                local_posq[get_local_id(0)] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            local_force[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
            if (!hasExclusions && flags != 0xFFFFFFFF) {
                if (flags == 0) {
                    // No interactions in this tile.
                }
                else {
                    // Compute only a subset of the interactions in this tile.

                    for (unsigned int j = 0; j < TileSize; j++) {
                        if ((flags&(1<<j)) != 0) {
                            bool isExcluded = false;
                            float4 delta = (float4) (local_posq[tbx+j].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float r = sqrt(r2);
                            float invR = 1.0f / r;
                            float4 posq2 = local_posq[tbx+j];
                            LOAD_ATOM2_PARAMETERS_J
                            float dEdR = 0.0f;
                            float tempEnergy = 0.0f;
                            COMPUTE_INTERACTION
#ifdef USE_CUTOFF
                            if (r2 > cutoffSquared) {
                                dEdR = 0.0f;
				tempEnergy = 0.0f;
                            }
#endif
			    energy += tempEnergy;
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            tempBuffer[get_local_id(0)] = delta;

                            // Sum the forces on atom j.

                            if (tgx % 2 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+1].xyz;
                            if (tgx % 4 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+2].xyz;
                            if (tgx % 8 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+4].xyz;
                            if (tgx % 16 == 0)
                                tempBuffer[get_local_id(0)].xyz += tempBuffer[get_local_id(0)+8].xyz;
                            if (tgx == 0)
                                local_force[tbx+j].xyz += tempBuffer[get_local_id(0)].xyz + tempBuffer[get_local_id(0)+16].xyz;
                        }
                    }
                }
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

                unsigned int xi = x/TileSize;
                unsigned int yi = y/TileSize;
                unsigned int tile = xi+yi*paddedNumAtoms/TileSize-yi*(yi+1)/2;
#ifdef USE_EXCLUSIONS
                unsigned int excl = (hasExclusions ? exclusions[exclusionIndices[tile]+tgx] : 0xFFFFFFFF);
                excl = (excl >> tgx) | (excl << (TileSize - tgx));
#endif
                for (unsigned int j = 0; j < TileSize; j++) {
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    float4 delta = (float4) (local_posq[tbx+tj].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float r = sqrt(r2);
                    float invR = 1.0f / r;
                    float4 posq2 = local_posq[tbx+tj];
                    LOAD_ATOM2_PARAMETERS_TJ
                    float dEdR = 0.0f;
                    float tempEnergy = 0.0f;
                    COMPUTE_INTERACTION
#if defined USE_CUTOFF || defined USE_EXCLUSIONS
    #if defined USE_CUTOFF && defined USE_EXCLUSIONS
                    excl >>= 1;
                    if (isExcluded || r2 > cutoffSquared) {
    #elif defined USE_CUTOFF
                    if (r2 > cutoffSquared) {
    #elif defined USE_EXCLUSIONS
                    excl >>= 1;
                    if (isExcluded) {
    #endif
                        dEdR = 0.0f;
			tempEnergy  = 0.0f;
                    }
#endif
		    energy += tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                    local_force[tbx+tj].xyz += delta.xyz;
                    tj = (tj + 1) & (TileSize - 1);
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x + tgx + (y/TileSize)*paddedNumAtoms;
            unsigned int offset2 = y + tgx + (x/TileSize)*paddedNumAtoms;
#else
            unsigned int offset1 = x + tgx + warp*paddedNumAtoms;
            unsigned int offset2 = y + tgx + warp*paddedNumAtoms;
#endif
            of = forceBuffers[offset1];
            of.xyz += force.xyz;
            forceBuffers[offset1] = of;
            of = forceBuffers[offset2];
            of.xyz += local_force[get_local_id(0)].xyz;
            forceBuffers[offset2] = of;
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
