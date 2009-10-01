const unsigned int TileSize = 32;
const float EpsilonFactor = 138.935485f;

/**
 * Compute nonbonded interactions.
 */

__kernel void computeNonbonded(int numTiles, int paddedNumAtoms, float cutoffSquared, float4 periodicBoxSize,
        __global float4* forceBuffers, __global float* energyBuffer, __global float4* posq, __global unsigned int* tiles,
        __global unsigned int* exclusions,  __global unsigned int* exclusionIndices, __local float4* local_posq, __local float4* local_force,
        __global float2* sigmaEpsilon, __local float2* local_sigmaEpsilon) {
    unsigned int totalWarps = get_global_size(0)/TileSize;
    unsigned int warp = get_global_id(0)/TileSize;
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
    float energy = 0.0f;
#ifdef USE_CUTOFF
    float3* tempBuffer = (float3*) &sA[cSim.nonbond_threads_per_block];
#endif
    unsigned int lasty = 0xFFFFFFFF;

    while (pos < end) {
        // Extract tile coordinates from appropriate work unit
        unsigned int x = tiles[pos];
        unsigned int y = ((x >> 2) & 0x7fff)*TileSize;
        bool hasExclusions = (x & 0x1);
        x = (x>>17)*TileSize;
        float4 apos;   // Local atom x, y, z, q
        float4 af = 0.0f;     // Local atom fx, fy, fz
        unsigned int tgx = get_local_id(0) & (TileSize-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int tj = tgx;
        unsigned int i = x + tgx;
        apos = posq[i];
        float2 a = sigmaEpsilon[i];
        if (x == y) {
            // Handle diagonals uniquely at 50% efficiency
            // Read fixed atom data into registers and GRF

            local_posq[get_local_id(0)] = apos;
            local_sigmaEpsilon[get_local_id(0)] = a;
            apos.w *= EpsilonFactor;
            unsigned int xi = x/TileSize;
            unsigned int tile = xi+xi*paddedNumAtoms/TileSize-xi*(xi+1)/2;
            unsigned int excl = exclusions[exclusionIndices[tile]+tgx];
            for (unsigned int j = 0; j < TileSize; j++) {
                bool isExcluded = !(excl & 0x1);
                float4 delta = (float4) (local_posq[tbx+j].xyz - apos.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                float invR = 1.0f / sqrt(r2);
                float sig = a.x + local_sigmaEpsilon[tbx+j].x;
                float sig2 = invR * sig;
                sig2 *= sig2;
                float sig6 = sig2 * sig2 * sig2;
                float eps = a.y * local_sigmaEpsilon[tbx+j].y;
                float dEdR = eps * (12.0f * sig6 - 6.0f) * sig6;
                float tempEnergy = eps * (sig6 - 1.0f) * sig6;
#ifdef USE_CUTOFF
                dEdR += apos.w * local_posq[tbx+j].w * (invR - 2.0f * cSim.reactionFieldK * r2);
                tempEnergy += apos.w * local_posq[tbx+j].w * (invR + cSim.reactionFieldK * r2 - cSim.reactionFieldC);
#else
                dEdR += apos.w * local_posq[tbx+j].w * invR;
                tempEnergy += apos.w * local_posq[tbx+j].w * invR;
#endif
                dEdR *= invR * invR;
#ifdef USE_CUTOFF
                if (isExcluded || r2 > cutoffSquared) {
#else
                if (isExcluded) {
#endif
                    dEdR = 0.0f;
                    tempEnergy  = 0.0f;
                }
                energy += 0.5f*tempEnergy;
                delta.xyz *= dEdR;
                af.xyz -= delta.xyz;
                excl >>= 1;
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset = x + tgx + warp*paddedNumAtoms;
            of = forceBuffers[offset];
            of.xyz += af.xyz;
            forceBuffers[offset] = of;
#else
            of.xyz = af.xyz;
            of.w = 0.0f;
            unsigned int offset = x + tgx + (x/TileSize) * paddedNumAtoms;
            forceBuffers[offset] = of;
#endif
        }
        else {
            // 100% utilization
            // Read fixed atom data into registers and GRF
            if (lasty != y) {
                unsigned int j = y + tgx;
                float2 temp1 = sigmaEpsilon[j];
                local_posq[get_local_id(0)] = posq[j];
                local_sigmaEpsilon[get_local_id(0)] = sigmaEpsilon[j];
            }
            local_force[get_local_id(0)] = 0.0f;
            apos.w *= EpsilonFactor;
#ifdef USE_CUTOFF
            unsigned int flags = cSim.pInteractionFlag[pos];
            if (!hasExclusions && flags != 0xFFFFFFFF) {
                if (flags == 0) {
                    // No interactions in this tile.
                }
                else {
                    // Compute only a subset of the interactions in this tile.

                    for (unsigned int j = 0; j < TileSize; j++) {
                        if ((flags&(1<<j)) != 0) {
                            bool isExcluded = false;
                            float4 delta = (float4) (local_posq[tbx+j].xyz - apos.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            float r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                            float invR = 1.0f / sqrt(r2);
                            float sig = a.x + local_sigmaEpsilon[tbx+j].x;
                            float sig2 = invR * sig;
                            sig2 *= sig2;
                            float sig6 = sig2 * sig2 * sig2;
                            float eps = a.y * local_sigmaEpsilon[tbx+j].Y;
                            float dEdR = eps * (12.0f * sig6 - 6.0f) * sig6;
			    float tempEnergy = eps * (sig6 - 1.0f) * sig6;
#ifdef USE_CUTOFF
                            dEdR += apos.w * local_posq[tbx+j].w * (invR - 2.0f * cSim.reactionFieldK * r2);
                            tempEnergy += apos.w * local_posq[tbx+j].w * (invR + cSim.reactionFieldK * r2 - cSim.reactionFieldC);
#else
                            dEdR += apos.w * local_posq[tbx+j].w * invR;
                            tempEnergy += apos.w * local_posq[tbx+j].w * invR;
#endif
                            dEdR *= invR * invR;
#ifdef USE_CUTOFF
                            if (r2 > cutoffSquared) {
                                dEdR = 0.0f;
				tempEnergy = 0.0f;
                            }
#endif
			    energy += tempEnergy;
                            delta.xyz *= dEdR;
                            af.xyz -= delta.xyz;
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
            else  // bExclusion
#endif
            {
                // Read fixed atom data into registers and GRF
                unsigned int xi = x/TileSize;
                unsigned int yi = y/TileSize;
                unsigned int tile = xi+yi*paddedNumAtoms/TileSize-yi*(yi+1)/2;
                unsigned int excl = exclusions[exclusionIndices[tile]+tgx];
                excl = (excl >> tgx) | (excl << (TileSize - tgx));
                for (unsigned int j = 0; j < TileSize; j++) {
                    bool isExcluded = !(excl & 0x1);
                    float4 delta = (float4) (local_posq[tbx+tj].xyz - apos.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x/periodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y/periodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z/periodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    float r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                    float invR = 1.0f / sqrt(r2);
                    float sig = a.x + local_sigmaEpsilon[tbx+j].x;
                    float sig2 = invR * sig;
                    sig2 *= sig2;
                    float sig6 = sig2 * sig2 * sig2;
                    float eps = a.y * local_sigmaEpsilon[tbx+j].y;
                    float dEdR = eps * (12.0f * sig6 - 6.0f) * sig6;
		    float tempEnergy = eps * (sig6 - 1.0f) * sig6;
#ifdef USE_CUTOFF
                    dEdR += apos.w * local_posq[tbx+j].w * (invR - 2.0f * cSim.reactionFieldK * r2);
		    tempEnergy += apos.w * local_posq[tbx+j].w * (invR + cSim.reactionFieldK * r2 - cSim.reactionFieldC);
#else
                    dEdR += apos.w * local_posq[tbx+j].w * invR;
                    tempEnergy += apos.w * local_posq[tbx+j].w * invR;
#endif
                    dEdR *= invR * invR;
#ifdef USE_CUTOFF
                    if (isExcluded || r2 > cutoffSquared) {
#else
                    if (isExcluded) {
#endif
                        dEdR = 0.0f;
			tempEnergy  = 0.0f;
                    }
		    energy += tempEnergy;
                    delta.xyz *= dEdR;
                    af.xyz -= delta.xyz;
                    local_force[tbx+tj].xyz += delta.xyz;
                    excl >>= 1;
                    tj = (tj + 1) & (TileSize - 1);
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_WARP
            unsigned int offset = x + tgx + warp*paddedNumAtoms;
            of = forceBuffers[offset];
            of.xyz += af.xyz;
            forceBuffers[offset] = of;
            offset = y + tgx + warp*paddedNumAtoms;
            of = forceBuffers[offset];
            of.xyz += local_foce[get_local_id(0)].xyz;
            forceBuffers[offset] = of;
#else
            of.xyz = af.xyz;
            of.w = 0.0f;
            unsigned int offset = x + tgx + (y/TileSize) * paddedNumAtoms;
            forceBuffers[offset] = of;
            of = local_force[get_local_id(0)];
            offset = y + tgx + (x/TileSize) * paddedNumAtoms;
            forceBuffers[offset] = of;
#endif
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
