#define TILE_SIZE 32

/**
 * Compute the Born sum.
 */

__kernel void computeBornSum(__global float* global_bornSum, __local float* local_bornSum, __global float4* posq,
        __local float4* local_posq, __global float2* global_params, __local float2* local_params, __local float* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount) {
#else
        unsigned int numTiles) {
#endif
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
        x = (x>>17)*TILE_SIZE;
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int atom1 = x + tgx;
        float bornSum = 0.0f;
        float4 posq1 = posq[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            local_params[get_local_id(0)] = params1;
            unsigned int xi = x/TILE_SIZE;
            unsigned int tile = xi+xi*PADDED_NUM_ATOMS/TILE_SIZE-xi*(xi+1)/2;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                float4 delta = (float4) (local_posq[tbx+j].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (atom1 < NUM_ATOMS && y+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                if (atom1 < NUM_ATOMS && y+j < NUM_ATOMS) {
#endif
                    float r = sqrt(r2);
                    float invR = 1.0f/r;
                    float2 params2 = local_params[tbx+j];
                    float rScaledRadiusJ = r+params2.y;
                    if ((j != tgx) && (params1.x < rScaledRadiusJ)) {
                        float l_ij = 1.0f/max(params1.x, fabs(r-params2.y));
                        float u_ij = 1.0f/rScaledRadiusJ;
                        float l_ij2 = l_ij*l_ij;
                        float u_ij2 = u_ij*u_ij;
                        float ratio = log(u_ij / l_ij);
                        bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                         (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                        if (params1.x < params2.x-r)
                            bornSum += 2.0f*(1.0f/params1.x-l_ij);
                    }
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_bornSum[offset] += bornSum;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y + tgx;
                local_posq[get_local_id(0)] = posq[j];
                local_params[get_local_id(0)] = global_params[j];
            }
            local_bornSum[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
            if (flags != 0xFFFFFFFF) {
                if (flags == 0) {
                    // No interactions in this tile.
                }
                else {
                    // Compute only a subset of the interactions in this tile.

                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        if ((flags&(1<<j)) != 0) {
                            float4 delta = (float4) (local_posq[tbx+j].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                            delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                            delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            tempBuffer[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
                            if (atom1 < NUM_ATOMS && y+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                            if (atom1 < NUM_ATOMS && y+j < NUM_ATOMS) {
#endif
                                float r = sqrt(r2);
                                float invR = 1.0f/r;
                                float2 params2 = local_params[tbx+j];
                                float rScaledRadiusJ = r+params2.y;
                                if (params1.x < rScaledRadiusJ) {
                                    float l_ij = 1.0f/max(params1.x, fabs(r-params2.y));
                                    float u_ij = 1.0f/rScaledRadiusJ;
                                    float l_ij2 = l_ij*l_ij;
                                    float u_ij2 = u_ij*u_ij;
                                    float ratio = log(u_ij / l_ij);
                                    bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                                     (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                                    if (params1.x < params2.x-r)
                                        bornSum += 2.0f*(1.0f/params1.x-l_ij);
                                }
                                float rScaledRadiusI = r+params1.y;
                                if (params2.x < rScaledRadiusI) {
                                    float l_ij = 1.0f/max(params2.x, fabs(r-params1.y));
                                    float u_ij = 1.0f/rScaledRadiusI;
                                    float l_ij2 = l_ij*l_ij;
                                    float u_ij2 = u_ij*u_ij;
                                    float ratio = log(u_ij / l_ij);
                                    float term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                                     (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                                    if (params2.x < params1.x-r)
                                        term += 2.0f*(1.0f/params2.x-l_ij);
                                    tempBuffer[get_local_id(0)] = term;
                                }
                            }

                            // Sum the forces on atom j.

                            if (tgx % 2 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+1];
                            if (tgx % 4 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+2];
                            if (tgx % 8 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+4];
                            if (tgx % 16 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+8];
                            if (tgx == 0)
                                local_bornSum[tbx+j] += tempBuffer[get_local_id(0)] + tempBuffer[get_local_id(0)+16];
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
                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    float4 delta = (float4) (local_posq[tbx+tj].xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                    delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                    delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && y+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && y+tj < NUM_ATOMS) {
#endif
                        float r = sqrt(r2);
                        float invR = 1.0f/r;
                        float2 params2 = local_params[tbx+tj];
                        float rScaledRadiusJ = r+params2.y;
                        if (params1.x < rScaledRadiusJ) {
                            float l_ij = 1.0f/max(params1.x, fabs(r-params2.y));
                            float u_ij = 1.0f/rScaledRadiusJ;
                            float l_ij2 = l_ij*l_ij;
                            float u_ij2 = u_ij*u_ij;
                            float ratio = log(u_ij / l_ij);
                            bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                             (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                            if (params1.x < params2.x-r)
                                bornSum += 2.0f*(1.0f/params1.x-l_ij);
                        }
                        float rScaledRadiusI = r+params1.y;
                        if (params2.x < rScaledRadiusI) {
                            float l_ij = 1.0f/max(params2.x, fabs(r-params1.y));
                            float u_ij = 1.0f/rScaledRadiusI;
                            float l_ij2 = l_ij*l_ij;
                            float u_ij2 = u_ij*u_ij;
                            float ratio = log(u_ij / l_ij);
                            float term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                             (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                            if (params2.x < params1.x-r)
                                term += 2.0f*(1.0f/params2.x-l_ij);
                            local_bornSum[tbx+tj] += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x + tgx + (y/TILE_SIZE)*PADDED_NUM_ATOMS;
            unsigned int offset2 = y + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset1 = x + tgx + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = y + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_bornSum[offset1] += bornSum;
            global_bornSum[offset2] += local_bornSum[get_local_id(0)];
            lasty = y;
        }
        pos++;
    }
}

/**
 * First part of computing the GBSA interaction.
 */

__kernel void computeGBSAForce1(__global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __local float4* local_posq, __local float4* local_force, __global float* global_bornRadii, __local float* local_bornRadii,
        __global float* global_bornForce, __local float* local_bornForce, __local float4* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount) {
#else
        unsigned int numTiles) {
#endif
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
        x = (x>>17)*TILE_SIZE;
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int tbx = get_local_id(0) - tgx;
        unsigned int atom1 = x + tgx;
        float4 force = 0.0f;
        float4 posq1 = posq[atom1];
        float bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            local_posq[get_local_id(0)] = posq1;
            local_bornRadii[get_local_id(0)] = bornRadius1;
            unsigned int xi = x/TILE_SIZE;
            unsigned int tile = xi+xi*PADDED_NUM_ATOMS/TILE_SIZE-xi*(xi+1)/2;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y+j < NUM_ATOMS) {
                    float4 posq2 = local_posq[tbx+j];
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                    delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                    delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float r = sqrt(r2);
                    float invR = 1.0f/r;
                    float bornRadius2 = local_bornRadii[tbx+j];
                    float alpha2_ij = bornRadius1*bornRadius2;
                    float D_ij = r2/(4.0f*alpha2_ij);
                    float expTerm = exp(-D_ij);
                    float denominator2 = r2 + alpha2_ij*expTerm;
                    float denominator = sqrt(denominator2);
                    float tempEnergy = (PREFACTOR*posq1.w*posq2.w)/denominator;
                    float Gpol = tempEnergy/denominator2;
                    float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                    force.w += dGpol_dalpha2_ij*bornRadius2;
                    float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                    if (r2 > CUTOFF_SQUARED) {
                        dEdR = 0.0f;
                        tempEnergy  = 0.0f;
                    }
#endif
                    energy += 0.5f*tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x + tgx + (x/TILE_SIZE)*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x + tgx + warp*PADDED_NUM_ATOMS;
#endif
            forceBuffers[offset].xyz += force.xyz;
            global_bornForce[offset] += force.w;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y + tgx;
                local_posq[get_local_id(0)] = posq[j];
                local_bornRadii[get_local_id(0)] = global_bornRadii[j];
            }
            local_force[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
            unsigned int flags = interactionFlags[pos];
            if (flags != 0xFFFFFFFF) {
                if (flags == 0) {
                    // No interactions in this tile.
                }
                else {
                    // Compute only a subset of the interactions in this tile.

                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        if ((flags&(1<<j)) != 0) {
                            float4 posq2 = local_posq[tbx+j];
                            float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                            delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                            delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float r = sqrt(r2);
                            float invR = 1.0f/r;
                            float bornRadius2 = local_bornRadii[tbx+j];
                            float alpha2_ij = bornRadius1*bornRadius2;
                            float D_ij = r2/(4.0f*alpha2_ij);
                            float expTerm = exp(-D_ij);
                            float denominator2 = r2 + alpha2_ij*expTerm;
                            float denominator = sqrt(denominator2);
                            float tempEnergy = (PREFACTOR*posq1.w*posq2.w)/denominator;
                            float Gpol = tempEnergy/denominator2;
                            float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
                            float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                            if (atom1 >= NUM_ATOMS || y+j >= NUM_ATOMS || r2 > CUTOFF_SQUARED) {
#else
                            if (atom1 >= NUM_ATOMS || y+j >= NUM_ATOMS) {
#endif
                                dEdR = 0.0f;
				tempEnergy = 0.0f;
                            }
			    energy += tempEnergy;
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            tempBuffer[get_local_id(0)] = (float4) (delta.xyz, dGpol_dalpha2_ij*bornRadius1);

                            // Sum the forces on atom j.

                            if (tgx % 2 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+1];
                            if (tgx % 4 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+2];
                            if (tgx % 8 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+4];
                            if (tgx % 16 == 0)
                                tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+8];
                            if (tgx == 0)
                                local_force[tbx+j] += tempBuffer[get_local_id(0)] + tempBuffer[get_local_id(0)+16];
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
                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    if (atom1 < NUM_ATOMS && y+tj < NUM_ATOMS) {
                        float4 posq2 = local_posq[tbx+tj];
                        float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x/PERIODIC_BOX_SIZE_X+0.5f)*PERIODIC_BOX_SIZE_X;
                        delta.y -= floor(delta.y/PERIODIC_BOX_SIZE_Y+0.5f)*PERIODIC_BOX_SIZE_Y;
                        delta.z -= floor(delta.z/PERIODIC_BOX_SIZE_Z+0.5f)*PERIODIC_BOX_SIZE_Z;
#endif
                        float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        float r = sqrt(r2);
                        float invR = 1.0f/r;
                        float bornRadius2 = local_bornRadii[tbx+tj];
                        float alpha2_ij = bornRadius1*bornRadius2;
                        float D_ij = r2/(4.0f*alpha2_ij);
                        float expTerm = exp(-D_ij);
                        float denominator2 = r2 + alpha2_ij*expTerm;
                        float denominator = sqrt(denominator2);
                        float tempEnergy = (PREFACTOR*posq1.w*posq2.w)/denominator;
                        float Gpol = tempEnergy/denominator2;
                        float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                        force.w += dGpol_dalpha2_ij*bornRadius2;
                        float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                        if (r2 > CUTOFF_SQUARED) {
                            dEdR = 0.0f;
                            tempEnergy  = 0.0f;
                        }
#endif
                        energy += tempEnergy;
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        local_force[tbx+tj] += (float4) (delta.xyz, dGpol_dalpha2_ij*bornRadius1);
                    }
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
            forceBuffers[offset2].xyz += local_force[get_local_id(0)].xyz;
            global_bornForce[offset1] += force.w;
            global_bornForce[offset2] += local_force[get_local_id(0)].w;
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
