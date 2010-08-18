#define TILE_SIZE 32

typedef struct {
    float x, y, z;
    float q;
    float fx, fy, fz, fw;
    float radius, scaledRadius;
    float bornSum;
    float bornRadius;
    float bornForce;
} AtomData;

/**
 * Compute the Born sum.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeBornSum(__global float* global_bornSum, __global float4* posq, __global float2* global_params, __local AtomData* localData, __local float* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
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
        float bornSum = 0.0f;
        float4 posq1 = posq[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            localData[get_local_id(0)].x = posq1.x;
            localData[get_local_id(0)].y = posq1.y;
            localData[get_local_id(0)].z = posq1.z;
            localData[get_local_id(0)].q = posq1.w;
            localData[get_local_id(0)].radius = params1.x;
            localData[get_local_id(0)].scaledRadius = params1.y;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                float4 delta = (float4) (localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                    float invR = RSQRT(r2);
                    float r = RECIP(invR);
                    float2 params2 = (float2) (localData[tbx+j].radius, localData[tbx+j].scaledRadius);
                    float rScaledRadiusJ = r+params2.y;
                    if ((j != tgx) && (params1.x < rScaledRadiusJ)) {
                        float l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                        float u_ij = RECIP(rScaledRadiusJ);
                        float l_ij2 = l_ij*l_ij;
                        float u_ij2 = u_ij*u_ij;
                        float ratio = LOG(u_ij * RECIP(l_ij));
                        bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                         (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                        if (params1.x < params2.x-r)
                            bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                    }
                }
            }

            // Write results
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_bornSum[offset] += bornSum;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y*TILE_SIZE + tgx;
                float4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                localData[get_local_id(0)].q = tempPosq.w;
                float2 tempParams = global_params[j];
                localData[get_local_id(0)].radius = tempParams.x;
                localData[get_local_id(0)].scaledRadius = tempParams.y;
            }
            localData[get_local_id(0)].bornSum = 0.0f;
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
                            float4 delta = (float4) (localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            tempBuffer[get_local_id(0)] = 0.0f;
#ifdef USE_CUTOFF
                            if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                            if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                                float invR = RSQRT(r2);
                                float r = RECIP(invR);
                                float2 params2 = (float2) (localData[tbx+j].radius, localData[tbx+j].scaledRadius);
                                float rScaledRadiusJ = r+params2.y;
                                if (params1.x < rScaledRadiusJ) {
                                    float l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                                    float u_ij = RECIP(rScaledRadiusJ);
                                    float l_ij2 = l_ij*l_ij;
                                    float u_ij2 = u_ij*u_ij;
                                    float ratio = LOG(u_ij * RECIP(l_ij));
                                    bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                                     (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                                    if (params1.x < params2.x-r)
                                        bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                                }
                                float rScaledRadiusI = r+params1.y;
                                if (params2.x < rScaledRadiusI) {
                                    float l_ij = RECIP(max(params2.x, fabs(r-params1.y)));
                                    float u_ij = RECIP(rScaledRadiusI);
                                    float l_ij2 = l_ij*l_ij;
                                    float u_ij2 = u_ij*u_ij;
                                    float ratio = LOG(u_ij * RECIP(l_ij));
                                    float term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                                     (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                                    if (params2.x < params1.x-r)
                                        term += 2.0f*(RECIP(params2.x)-l_ij);
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
                                localData[tbx+j].bornSum += tempBuffer[get_local_id(0)] + tempBuffer[get_local_id(0)+16];
                        }
                    }
                }
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    float4 delta = (float4) (localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
#endif
                        float invR = RSQRT(r2);
                        float r = RECIP(invR);
                        float2 params2 = (float2) (localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
                        float rScaledRadiusJ = r+params2.y;
                        if (params1.x < rScaledRadiusJ) {
                            float l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                            float u_ij = RECIP(rScaledRadiusJ);
                            float l_ij2 = l_ij*l_ij;
                            float u_ij2 = u_ij*u_ij;
                            float ratio = LOG(u_ij * RECIP(l_ij));
                            bornSum += l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                             (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2);
                            if (params1.x < params2.x-r)
                                bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                        }
                        float rScaledRadiusI = r+params1.y;
                        if (params2.x < rScaledRadiusI) {
                            float l_ij = RECIP(max(params2.x, fabs(r-params1.y)));
                            float u_ij = RECIP(rScaledRadiusI);
                            float l_ij2 = l_ij*l_ij;
                            float u_ij2 = u_ij*u_ij;
                            float ratio = LOG(u_ij * RECIP(l_ij));
                            float term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                             (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                            if (params2.x < params1.x-r)
                                term += 2.0f*(RECIP(params2.x)-l_ij);
                            localData[tbx+tj].bornSum += term;
                        }
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results
            float4 of;
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
            unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
            unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            global_bornSum[offset1] += bornSum;
            global_bornSum[offset2] += localData[get_local_id(0)].bornSum;
            lasty = y;
        }
        pos++;
    }
}

/**
 * First part of computing the GBSA interaction.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeGBSAForce1(__global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float* global_bornRadii,
        __global float* global_bornForce, __local AtomData* localData, __local float4* tempBuffer, __global unsigned int* tiles,
#ifdef USE_CUTOFF
        __global unsigned int* interactionFlags, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize) {
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
        float bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            localData[get_local_id(0)].x = posq1.x;
            localData[get_local_id(0)].y = posq1.y;
            localData[get_local_id(0)].z = posq1.z;
            localData[get_local_id(0)].q = posq1.w;
            localData[get_local_id(0)].bornRadius = bornRadius1;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
                    float4 posq2 = (float4) (localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z, localData[tbx+j].q);
                    float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    float invR = RSQRT(r2);
                    float r = RECIP(invR);
                    float bornRadius2 = localData[tbx+j].bornRadius;
                    float alpha2_ij = bornRadius1*bornRadius2;
                    float D_ij = r2*RECIP(4.0f*alpha2_ij);
                    float expTerm = EXP(-D_ij);
                    float denominator2 = r2 + alpha2_ij*expTerm;
                    float denominator = SQRT(denominator2);
                    float tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                    float Gpol = tempEnergy*RECIP(denominator2);
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
            unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
            unsigned int offset = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
#endif
            forceBuffers[offset].xyz += force.xyz;
            global_bornForce[offset] += force.w;
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y) {
                unsigned int j = y*TILE_SIZE + tgx;
                float4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                localData[get_local_id(0)].q = tempPosq.w;
                localData[get_local_id(0)].bornRadius = global_bornRadii[j];
            }
            localData[get_local_id(0)].fx = 0.0f;
            localData[get_local_id(0)].fy = 0.0f;
            localData[get_local_id(0)].fz = 0.0f;
            localData[get_local_id(0)].fw = 0.0f;
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
                            float4 posq2 = (float4) (localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z, localData[tbx+j].q);
                            float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                            float invR = RSQRT(r2);
                            float r = RECIP(invR);
                            float bornRadius2 = localData[tbx+j].bornRadius;
                            float alpha2_ij = bornRadius1*bornRadius2;
                            float D_ij = r2*RECIP(4.0f*alpha2_ij);
                            float expTerm = EXP(-D_ij);
                            float denominator2 = r2 + alpha2_ij*expTerm;
                            float denominator = SQRT(denominator2);
                            float tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                            float Gpol = tempEnergy*RECIP(denominator2);
                            float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
                            float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                            if (atom1 >= NUM_ATOMS || y*TILE_SIZE+j >= NUM_ATOMS || r2 > CUTOFF_SQUARED) {
#else
                            if (atom1 >= NUM_ATOMS || y*TILE_SIZE+j >= NUM_ATOMS) {
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
                            if (tgx == 0) {
                                float4 sum = tempBuffer[get_local_id(0)] + tempBuffer[get_local_id(0)+16];
                                localData[tbx+j].fx += sum.x;
                                localData[tbx+j].fy += sum.y;
                                localData[tbx+j].fz += sum.z;
                                localData[tbx+j].fw += sum.w;
                            }
                        }
                    }
                }
            }
            else
#endif
            {
                // Compute the full set of interactions in this tile.

                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
                        float4 posq2 = (float4) (localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z, localData[tbx+tj].q);
                        float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                        float invR = RSQRT(r2);
                        float r = RECIP(invR);
                        float bornRadius2 = localData[tbx+tj].bornRadius;
                        float alpha2_ij = bornRadius1*bornRadius2;
                        float D_ij = r2*RECIP(4.0f*alpha2_ij);
                        float expTerm = EXP(-D_ij);
                        float denominator2 = r2 + alpha2_ij*expTerm;
                        float denominator = SQRT(denominator2);
                        float tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                        float Gpol = tempEnergy*RECIP(denominator2);
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
                        localData[tbx+tj].fx += delta.x;
                        localData[tbx+tj].fy += delta.y;
                        localData[tbx+tj].fz += delta.z;
                        localData[tbx+tj].fw += dGpol_dalpha2_ij*bornRadius1;
                    }
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
            forceBuffers[offset2] += (float4) (localData[get_local_id(0)].fx, localData[get_local_id(0)].fy, localData[get_local_id(0)].fz, 0);
            global_bornForce[offset1] += force.w;
            global_bornForce[offset2] += localData[get_local_id(0)].fw;
            lasty = y;
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
