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
void computeBornSum(__global float* global_bornSum, __global float4* posq, __global float2* global_params, __local AtomData* localData, __local float* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles) {
#else
        unsigned int numTiles) {
#endif
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = get_group_id(0)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0);
    unsigned int end = (get_group_id(0)+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0);
#else
    unsigned int pos = get_group_id(0)*numTiles/get_num_groups(0);
    unsigned int end = (get_group_id(0)+1)*numTiles/get_num_groups(0);
#endif
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
            if (x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y++;
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
        }
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int forceBufferOffset = (tgx < TILE_SIZE/2 ? 0 : TILE_SIZE);
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
            barrier(CLK_LOCAL_MEM_FENCE);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                float4 delta = (float4) (localData[baseLocalAtom+j].x-posq1.x, localData[baseLocalAtom+j].y-posq1.y, localData[baseLocalAtom+j].z-posq1.z, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float invR = RSQRT(r2);
                float r = RECIP(invR);
                float2 params2 = (float2) (localData[baseLocalAtom+j].radius, localData[baseLocalAtom+j].scaledRadius);
                float rScaledRadiusJ = r+params2.y;
#ifdef USE_CUTOFF
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS && r2 < CUTOFF_SQUARED && (j+baseLocalAtom != tgx) && (params1.x < rScaledRadiusJ));
#else
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS && (j+baseLocalAtom != tgx) && (params1.x < rScaledRadiusJ));
#endif
                float l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                float u_ij = RECIP(rScaledRadiusJ);
                float l_ij2 = l_ij*l_ij;
                float u_ij2 = u_ij*u_ij;
                float ratio = LOG(u_ij * RECIP(l_ij));
                bornSum += select(0.0f, l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                 (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2), includeInteraction);
                bornSum += select(0.0f, 2.0f*(RECIP(params1.x)-l_ij), includeInteraction && params1.x < params2.x-r);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = bornSum;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                global_bornSum[offset] += bornSum+tempBuffer[get_local_id(0)+TILE_SIZE];
            }
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
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
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx%(TILE_SIZE/2);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                float4 delta = (float4) (localData[baseLocalAtom+tj].x-posq1.x, localData[baseLocalAtom+tj].y-posq1.y, localData[baseLocalAtom+tj].z-posq1.z, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED);
#else
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+tj < NUM_ATOMS);
#endif
                float invR = RSQRT(r2);
                float r = RECIP(invR);
                float2 params2 = (float2) (localData[baseLocalAtom+tj].radius, localData[baseLocalAtom+tj].scaledRadius);
                float rScaledRadiusJ = r+params2.y;
                {
                    float l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                    float u_ij = RECIP(rScaledRadiusJ);
                    float l_ij2 = l_ij*l_ij;
                    float u_ij2 = u_ij*u_ij;
                    float ratio = LOG(u_ij * RECIP(l_ij));
                    unsigned int includeTerm = (includeInteraction && params1.x < rScaledRadiusJ);
                    bornSum += select(0.0f, l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                     (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2), includeTerm);
                    bornSum += select(0.0f, 2.0f*(RECIP(params1.x)-l_ij), includeTerm && params1.x < params2.x-r);
                }
                float rScaledRadiusI = r+params1.y;
                {
                    float l_ij = RECIP(max(params2.x, fabs(r-params1.y)));
                    float u_ij = RECIP(rScaledRadiusI);
                    float l_ij2 = l_ij*l_ij;
                    float u_ij2 = u_ij*u_ij;
                    float ratio = LOG(u_ij * RECIP(l_ij));
                    float term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                     (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                    term += select(0.0f, 2.0f*(RECIP(params2.x)-l_ij), params2.x < params1.x-r);
                    localData[baseLocalAtom+tj+forceBufferOffset].bornSum += select(0.0f, term, includeInteraction && params2.x < rScaledRadiusI);
                }
                barrier(CLK_LOCAL_MEM_FENCE);
                tj = (tj+1)%(TILE_SIZE/2);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = bornSum;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                global_bornSum[offset1] += bornSum+tempBuffer[get_local_id(0)+TILE_SIZE];
                global_bornSum[offset2] += localData[get_local_id(0)].bornSum+localData[get_local_id(0)+TILE_SIZE].bornSum;
            }
        }
        lasty = y;
        pos++;
    }
}

/**
 * First part of computing the GBSA interaction.
 */

__kernel __attribute__((reqd_work_group_size(WORK_GROUP_SIZE, 1, 1)))
void computeGBSAForce1(__global float4* forceBuffers, __global float* energyBuffer,
        __global float4* posq, __global float* global_bornRadii,
        __global float* global_bornForce, __local AtomData* localData, __local float4* tempBuffer,
#ifdef USE_CUTOFF
        __global ushort2* tiles, __global unsigned int* interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles) {
#else
        unsigned int numTiles) {
#endif
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
            if (x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y++;
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
        }
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int forceBufferOffset = (tgx < TILE_SIZE/2 ? 0 : TILE_SIZE);
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
            barrier(CLK_LOCAL_MEM_FENCE);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS);
                float4 posq2 = (float4) (localData[baseLocalAtom+j].x, localData[baseLocalAtom+j].y, localData[baseLocalAtom+j].z, localData[baseLocalAtom+j].q);
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float invR = RSQRT(r2);
                float r = RECIP(invR);
                float bornRadius2 = localData[baseLocalAtom+j].bornRadius;
                float alpha2_ij = bornRadius1*bornRadius2;
                float D_ij = r2*RECIP(4.0f*alpha2_ij);
                float expTerm = EXP(-D_ij);
                float denominator2 = r2 + alpha2_ij*expTerm;
                float denominator = SQRT(denominator2);
                float tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                float Gpol = tempEnergy*RECIP(denominator2);
                float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                dEdR = select(dEdR, 0.0f, r2 > CUTOFF_SQUARED);
                tempEnergy = select(tempEnergy, 0.0f, r2 > CUTOFF_SQUARED);
                dGpol_dalpha2_ij = select(dGpol_dalpha2_ij, 0.0f, r2 > CUTOFF_SQUARED);
#endif
                force.w += select(0.0f, dGpol_dalpha2_ij*bornRadius2, includeInteraction);
                energy += select(0.0f, 0.5f*tempEnergy, includeInteraction);
                delta.xyz *= select(0.0f, dEdR, includeInteraction);
                force.xyz -= delta.xyz;
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = force;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz+tempBuffer[get_local_id(0)+TILE_SIZE].xyz;
                global_bornForce[offset] += force.w+tempBuffer[get_local_id(0)+TILE_SIZE].w;
            }
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
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
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx%(TILE_SIZE/2);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+tj < NUM_ATOMS);
                float4 posq2 = (float4) (localData[baseLocalAtom+tj].x, localData[baseLocalAtom+tj].y, localData[baseLocalAtom+tj].z, localData[baseLocalAtom+tj].q);
                float4 delta = (float4) (posq2.xyz - posq1.xyz, 0.0f);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                float r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                float invR = RSQRT(r2);
                float r = RECIP(invR);
                float bornRadius2 = localData[baseLocalAtom+tj].bornRadius;
                float alpha2_ij = bornRadius1*bornRadius2;
                float D_ij = r2*RECIP(4.0f*alpha2_ij);
                float expTerm = EXP(-D_ij);
                float denominator2 = r2 + alpha2_ij*expTerm;
                float denominator = SQRT(denominator2);
                float tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                float Gpol = tempEnergy*RECIP(denominator2);
                float dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                dEdR = select(dEdR, 0.0f, r2 > CUTOFF_SQUARED);
                tempEnergy = select(tempEnergy, 0.0f, r2 > CUTOFF_SQUARED);
                dGpol_dalpha2_ij = select(dGpol_dalpha2_ij, 0.0f, r2 > CUTOFF_SQUARED);
#endif
                force.w += select(0.0f, dGpol_dalpha2_ij*bornRadius2, includeInteraction);
                energy += select(0.0f, tempEnergy, includeInteraction);
                delta.xyz *= select(0.0f, dEdR, includeInteraction);
                force.xyz -= delta.xyz;
                localData[baseLocalAtom+tj+forceBufferOffset].fx += delta.x;
                localData[baseLocalAtom+tj+forceBufferOffset].fy += delta.y;
                localData[baseLocalAtom+tj+forceBufferOffset].fz += delta.z;
                localData[baseLocalAtom+tj+forceBufferOffset].fw += select(0.0f, dGpol_dalpha2_ij*bornRadius1, includeInteraction);
                barrier(CLK_LOCAL_MEM_FENCE);
                tj = (tj+1)%(TILE_SIZE/2);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                tempBuffer[get_local_id(0)] = force;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                forceBuffers[offset1].xyz = forceBuffers[offset1].xyz+force.xyz+tempBuffer[get_local_id(0)+TILE_SIZE].xyz;
                float4 sum = (float4) (localData[get_local_id(0)].fx+localData[get_local_id(0)+TILE_SIZE].fx,
                                       localData[get_local_id(0)].fy+localData[get_local_id(0)+TILE_SIZE].fy,
                                       localData[get_local_id(0)].fz+localData[get_local_id(0)+TILE_SIZE].fz,
                                       localData[get_local_id(0)].fw+localData[get_local_id(0)+TILE_SIZE].fw);
                forceBuffers[offset2].xyz = forceBuffers[offset2].xyz+sum.xyz;
                global_bornForce[offset1] += force.w+tempBuffer[get_local_id(0)+TILE_SIZE].w;
                global_bornForce[offset2] += sum.w;
            }
        }
        lasty = y;
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
