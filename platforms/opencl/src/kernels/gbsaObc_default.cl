#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif
#define TILE_SIZE 32

typedef struct {
    real x, y, z;
    float radius, scaledRadius;
} AtomData1;

/**
 * Compute the Born sum.
 */

__kernel __attribute__((reqd_work_group_size(FORCE_WORK_GROUP_SIZE, 1, 1)))
void computeBornSum(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict global_bornSum,
#else
        __global real* restrict global_bornSum,
#endif
        __global const real4* restrict posq, __global const float2* restrict global_params,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles) {
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
    __local AtomData1 localData[TILE_SIZE];
    __local real localBornSum[FORCE_WORK_GROUP_SIZE];
    __local real localTemp[TILE_SIZE];

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
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int localForceOffset = get_local_id(0) & ~(TILE_SIZE-1);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real bornSum = 0.0f;
        real4 posq1 = posq[atom1];
        float2 params1 = global_params[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            if (get_local_id(0) < TILE_SIZE) {
                localData[get_local_id(0)].x = posq1.x;
                localData[get_local_id(0)].y = posq1.y;
                localData[get_local_id(0)].z = posq1.z;
                localData[get_local_id(0)].radius = params1.x;
                localData[get_local_id(0)].scaledRadius = params1.y;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                real4 delta = (real4) (localData[baseLocalAtom+j].x-posq1.x, localData[baseLocalAtom+j].y-posq1.y, localData[baseLocalAtom+j].z-posq1.z, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                float2 params2 = (float2) (localData[baseLocalAtom+j].radius, localData[baseLocalAtom+j].scaledRadius);
                real rScaledRadiusJ = r+params2.y;
#ifdef USE_CUTOFF
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS && r2 < CUTOFF_SQUARED && (j+baseLocalAtom != tgx) && (params1.x < rScaledRadiusJ));
#else
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS && (j+baseLocalAtom != tgx) && (params1.x < rScaledRadiusJ));
#endif
                real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                real u_ij = RECIP(rScaledRadiusJ);
                real l_ij2 = l_ij*l_ij;
                real u_ij2 = u_ij*u_ij;
                real ratio = LOG(u_ij * RECIP(l_ij));
                bornSum += (includeInteraction ? l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                 (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2) : (real) 0);
                bornSum += (includeInteraction && params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : (real) 0);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                localTemp[tgx] = bornSum;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                const unsigned int offset = x*TILE_SIZE + tgx;
                atom_add(&global_bornSum[offset], (long) ((bornSum + localTemp[tgx])*0xFFFFFFFF));
#else
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                const unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                const unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                global_bornSum[offset] += bornSum + localTemp[tgx];
#endif
            }
            // barrier not required here as localTemp is not accessed before encountering another barrier.
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
                unsigned int j = y*TILE_SIZE + tgx;
                real4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                float2 tempParams = global_params[j];
                localData[get_local_id(0)].radius = tempParams.x;
                localData[get_local_id(0)].scaledRadius = tempParams.y;
            }
            localBornSum[get_local_id(0)] = 0.0f;
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

            unsigned int tj = (tgx+baseLocalAtom) & (TILE_SIZE-1);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                real4 delta = (real4) (localData[tj].x-posq1.x, localData[tj].y-posq1.y, localData[tj].z-posq1.z, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED);
#else
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS);
#endif
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                float2 params2 = (float2) (localData[tj].radius, localData[tj].scaledRadius);
                real rScaledRadiusJ = r+params2.y;
                {
                    real l_ij = RECIP(max((real) params1.x, fabs(r-params2.y)));
                    real u_ij = RECIP(rScaledRadiusJ);
                    real l_ij2 = l_ij*l_ij;
                    real u_ij2 = u_ij*u_ij;
                    real ratio = LOG(u_ij * RECIP(l_ij));
                    unsigned int includeTerm = (includeInteraction && params1.x < rScaledRadiusJ);
                    bornSum += (includeTerm ? l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                     (0.25f*params2.y*params2.y*invR)*(l_ij2-u_ij2) : (real) 0);
                    bornSum += (includeTerm && params1.x < params2.y-r ? 2.0f*(RECIP(params1.x)-l_ij) : (real) 0);
                }
                real rScaledRadiusI = r+params1.y;
                {
                    real l_ij = RECIP(max((real) params2.x, fabs(r-params1.y)));
                    real u_ij = RECIP(rScaledRadiusI);
                    real l_ij2 = l_ij*l_ij;
                    real u_ij2 = u_ij*u_ij;
                    real ratio = LOG(u_ij * RECIP(l_ij));
                    real term = l_ij - u_ij + 0.25f*r*(u_ij2-l_ij2) + (0.50f*invR*ratio) +
                                     (0.25f*params1.y*params1.y*invR)*(l_ij2-u_ij2);
                    term += (params2.x < params1.y-r ? 2.0f*(RECIP(params2.x)-l_ij) : (real) 0);
                    localBornSum[tj+localForceOffset] += (includeInteraction && params2.x < rScaledRadiusI ? term : (real) 0);
                }
                barrier(CLK_LOCAL_MEM_FENCE);
                tj = (tj+1) & (TILE_SIZE-1);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE)
                localTemp[tgx] = bornSum;
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                const unsigned int offset1 = x*TILE_SIZE + tgx;
                const unsigned int offset2 = y*TILE_SIZE + tgx;
                atom_add(&global_bornSum[offset1], (long) ((bornSum + localTemp[tgx])*0xFFFFFFFF));
                atom_add(&global_bornSum[offset2], (long) ((localBornSum[get_local_id(0)] + localBornSum[get_local_id(0)+TILE_SIZE])*0xFFFFFFFF));
#else
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                const unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                const unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                // Do both loads before both stores to minimize store-load waits.
                real sum1 = global_bornSum[offset1];
                real sum2 = global_bornSum[offset2];
                sum1 += bornSum + localTemp[tgx];
                sum2 += localBornSum[get_local_id(0)] + localBornSum[get_local_id(0)+TILE_SIZE];
                global_bornSum[offset1] = sum1;
                global_bornSum[offset2] = sum2;
#endif
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        lasty = y;
        pos++;
    }
}

typedef struct {
    real x, y, z, w;
    real padding;
} PaddedUnalignedFloat4;

typedef struct {
    real x, y, z;
    real q;
    real bornRadius;
    real temp_x, temp_y, temp_z, temp_w;
} AtomData2;

/**
 * First part of computing the GBSA interaction.
 */

__kernel __attribute__((reqd_work_group_size(FORCE_WORK_GROUP_SIZE, 1, 1)))
void computeGBSAForce1(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers, __global long* restrict global_bornForce,
#else
        __global real4* restrict forceBuffers, __global real* restrict global_bornForce,
#endif
        __global real* restrict energyBuffer, __global const real4* restrict posq, __global const real* restrict global_bornRadii,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles) {
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
    real energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;
    __local AtomData2 localData[TILE_SIZE];
    __local PaddedUnalignedFloat4 localForce[FORCE_WORK_GROUP_SIZE];

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
        unsigned int baseLocalAtom = (get_local_id(0) < TILE_SIZE ? 0 : TILE_SIZE/2);
        unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        unsigned int localForceOffset = get_local_id(0) & ~(TILE_SIZE-1);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 force = 0.0f;
        real4 posq1 = posq[atom1];
        real bornRadius1 = global_bornRadii[atom1];
        if (x == y) {
            // This tile is on the diagonal.

            if (get_local_id(0) < TILE_SIZE) {
                localData[get_local_id(0)].x = posq1.x;
                localData[get_local_id(0)].y = posq1.y;
                localData[get_local_id(0)].z = posq1.z;
                localData[get_local_id(0)].q = posq1.w;
                localData[get_local_id(0)].bornRadius = bornRadius1;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+baseLocalAtom+j < NUM_ATOMS);
                real4 posq2 = (real4) (localData[baseLocalAtom+j].x, localData[baseLocalAtom+j].y, localData[baseLocalAtom+j].z, localData[baseLocalAtom+j].q);
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                real bornRadius2 = localData[baseLocalAtom+j].bornRadius;
                real alpha2_ij = bornRadius1*bornRadius2;
                real D_ij = r2*RECIP(4.0f*alpha2_ij);
                real expTerm = EXP(-D_ij);
                real denominator2 = r2 + alpha2_ij*expTerm;
                real denominator = SQRT(denominator2);
                real tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                real Gpol = tempEnergy*RECIP(denominator2);
                real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                real dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                dEdR = (r2 > CUTOFF_SQUARED ? (real) 0 : dEdR);
                tempEnergy = (r2 > CUTOFF_SQUARED ? (real) 0 : tempEnergy);
                dGpol_dalpha2_ij = (r2 > CUTOFF_SQUARED ? (real) 0 : dGpol_dalpha2_ij);
#endif
                force.w += (includeInteraction ? dGpol_dalpha2_ij*bornRadius2 : (real) 0);
                energy += (includeInteraction ? 0.5f*tempEnergy : (real) 0);
                delta.xyz *= (includeInteraction ? dEdR : (real) 0);
                force.xyz -= delta.xyz;
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE) {
                localData[tgx].temp_x = force.x;
                localData[tgx].temp_y = force.y;
                localData[tgx].temp_z = force.z;
                localData[tgx].temp_w = force.w;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                const unsigned int offset = x*TILE_SIZE + tgx;
                atom_add(&forceBuffers[offset], (long) ((force.x + localData[tgx].temp_x)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) ((force.y + localData[tgx].temp_y)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) ((force.z + localData[tgx].temp_z)*0xFFFFFFFF));
                atom_add(&global_bornForce[offset], (long) ((force.w + localData[tgx].temp_w)*0xFFFFFFFF));
#else
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                const unsigned int offset = x*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                const unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                // Cheaper to load/store real4 than real3. Do all loads before all stores to minimize store-load waits.
                real4 sum = forceBuffers[offset];
                real global_sum = global_bornForce[offset];
                sum.x += force.x + localData[tgx].temp_x;
                sum.y += force.y + localData[tgx].temp_y;
                sum.z += force.z + localData[tgx].temp_z;
                global_sum += force.w + localData[tgx].temp_w;
                forceBuffers[offset] = sum;
                global_bornForce[offset] = global_sum;
#endif
            }
            // barrier not required here as localData[*]/temp_* is not accessed before encountering another barrier.
        }
        else {
            // This is an off-diagonal tile.

            if (lasty != y && get_local_id(0) < TILE_SIZE) {
                unsigned int j = y*TILE_SIZE + tgx;
                real4 tempPosq = posq[j];
                localData[get_local_id(0)].x = tempPosq.x;
                localData[get_local_id(0)].y = tempPosq.y;
                localData[get_local_id(0)].z = tempPosq.z;
                localData[get_local_id(0)].q = tempPosq.w;
                localData[get_local_id(0)].bornRadius = global_bornRadii[j];
            }
            localForce[get_local_id(0)].x = 0.0f;
            localForce[get_local_id(0)].y = 0.0f;
            localForce[get_local_id(0)].z = 0.0f;
            localForce[get_local_id(0)].w = 0.0f;
            barrier(CLK_LOCAL_MEM_FENCE);

            // Compute the full set of interactions in this tile.

            unsigned int tj = (tgx+baseLocalAtom) & (TILE_SIZE-1);
            for (unsigned int j = 0; j < TILE_SIZE/2; j++) {
                unsigned int includeInteraction = (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS);
                real4 posq2 = (real4) (localData[tj].x, localData[tj].y, localData[tj].z, localData[tj].q);
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = RECIP(invR);
                real bornRadius2 = localData[tj].bornRadius;
                real alpha2_ij = bornRadius1*bornRadius2;
                real D_ij = r2*RECIP(4.0f*alpha2_ij);
                real expTerm = EXP(-D_ij);
                real denominator2 = r2 + alpha2_ij*expTerm;
                real denominator = SQRT(denominator2);
                real tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                real Gpol = tempEnergy*RECIP(denominator2);
                real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                real dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                dEdR = (r2 > CUTOFF_SQUARED ? (real) 0 : dEdR);
                tempEnergy = (r2 > CUTOFF_SQUARED ? (real) 0 : tempEnergy);
                dGpol_dalpha2_ij = (r2 > CUTOFF_SQUARED ? (real) 0 : dGpol_dalpha2_ij);
#endif
                force.w += (includeInteraction ? dGpol_dalpha2_ij*bornRadius2 : (real) 0);
                energy += (includeInteraction ? tempEnergy : (real) 0);
                delta.xyz *= (includeInteraction ? dEdR : (real) 0);
                force.xyz -= delta.xyz;
                localForce[tj+localForceOffset].x += delta.x;
                localForce[tj+localForceOffset].y += delta.y;
                localForce[tj+localForceOffset].z += delta.z;
                localForce[tj+localForceOffset].w += (includeInteraction ? dGpol_dalpha2_ij*bornRadius1 : (real) 0);
                barrier(CLK_LOCAL_MEM_FENCE);
                tj = (tj+1) & (TILE_SIZE-1);
            }

            // Sum the forces and write results.

            if (get_local_id(0) >= TILE_SIZE) {
                localData[tgx].temp_x = force.x;
                localData[tgx].temp_y = force.y;
                localData[tgx].temp_z = force.z;
                localData[tgx].temp_w = force.w;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (get_local_id(0) < TILE_SIZE) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                barrier(CLK_LOCAL_MEM_FENCE);
                const unsigned int offset1 = x*TILE_SIZE + tgx;
                const unsigned int offset2 = y*TILE_SIZE + tgx;
                atom_add(&forceBuffers[offset1], (long) ((force.x + localData[tgx].temp_x)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset1+PADDED_NUM_ATOMS], (long) ((force.y + localData[tgx].temp_y)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset1+2*PADDED_NUM_ATOMS], (long) ((force.z + localData[tgx].temp_z)*0xFFFFFFFF));
                atom_add(&global_bornForce[offset1], (long) ((force.w + localData[tgx].temp_w)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2], (long) ((localForce[get_local_id(0)].x + localForce[get_local_id(0)+TILE_SIZE].x)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2+PADDED_NUM_ATOMS], (long) ((localForce[get_local_id(0)].y + localForce[get_local_id(0)+TILE_SIZE].y)*0xFFFFFFFF));
                atom_add(&forceBuffers[offset2+2*PADDED_NUM_ATOMS], (long) ((localForce[get_local_id(0)].z + localForce[get_local_id(0)+TILE_SIZE].z)*0xFFFFFFFF));
                atom_add(&global_bornForce[offset2], (long) ((localForce[get_local_id(0)].w + localForce[get_local_id(0)+TILE_SIZE].w)*0xFFFFFFFF));
#else
#ifdef USE_OUTPUT_BUFFER_PER_BLOCK
                const unsigned int offset1 = x*TILE_SIZE + tgx + y*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + x*PADDED_NUM_ATOMS;
#else
                const unsigned int offset1 = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                const unsigned int offset2 = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
#endif
                // Cheaper to load/store real4 than real3. Do all loads before all stores to minimize store-load waits.
                real4 sum1 = forceBuffers[offset1];
                real4 sum2 = forceBuffers[offset2];
                real global_sum1 = global_bornForce[offset1];
                real global_sum2 = global_bornForce[offset2];
                sum1.x += force.x + localData[tgx].temp_x;
                sum1.y += force.y + localData[tgx].temp_y;
                sum1.z += force.z + localData[tgx].temp_z;
                global_sum1 += force.w + localData[tgx].temp_w;
                sum2.x += localForce[get_local_id(0)].x + localForce[get_local_id(0)+TILE_SIZE].x;
                sum2.y += localForce[get_local_id(0)].y + localForce[get_local_id(0)+TILE_SIZE].y;
                sum2.z += localForce[get_local_id(0)].z + localForce[get_local_id(0)+TILE_SIZE].z;
                global_sum2 += localForce[get_local_id(0)].w + localForce[get_local_id(0)+TILE_SIZE].w;
                forceBuffers[offset1] = sum1;
                forceBuffers[offset2] = sum2;
                global_bornForce[offset1] = global_sum1;
                global_bornForce[offset2] = global_sum2;
#endif
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        lasty = y;
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
}
