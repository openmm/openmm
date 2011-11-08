#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif
#define TILE_SIZE 32
#define WARPS_PER_GROUP (FORCE_WORK_GROUP_SIZE/TILE_SIZE)

typedef struct {
    float x, y, z;
    float q;
    float radius, scaledRadius;
    float bornSum;
} AtomData1;

/**
 * Compute the Born sum.
 */
__kernel void computeBornSum(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict global_bornSum,
#else
        __global float* restrict global_bornSum,
#endif
        __global const float4* restrict posq, __global const float2* restrict global_params,
        __local AtomData1* restrict localData,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles, __global const unsigned int* restrict interactionFlags,
#else
        unsigned int numTiles,
#endif
        __global unsigned int* exclusionIndices, __global unsigned int* exclusionRowIndices) {
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    unsigned int lasty = 0xFFFFFFFF;
    __local float tempBuffer[FORCE_WORK_GROUP_SIZE];
    __local int2 reservedBlocks[WARPS_PER_GROUP];
    __local unsigned int* exclusionRange = (__local unsigned int*) reservedBlocks;
    __local int exclusionIndex[WARPS_PER_GROUP];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        const unsigned int tbx = get_local_id(0) - tgx;
        const unsigned int localGroupIndex = get_local_id(0)/TILE_SIZE;
        unsigned int x, y;
        float bornSum = 0.0f;
        if (pos < end) {
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
            unsigned int atom1 = x*TILE_SIZE + tgx;
            float4 posq1 = posq[atom1];
            float2 params1 = global_params[atom1];
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
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
                            bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                             (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                            if (params1.x < params2.y-r)
                                bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                        }
                    }
                }
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
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                bool computeSubset = false;
                if (flags != 0xFFFFFFFF) {
                    if (tgx < 2)
                        exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
                    if (tgx == 0)
                        exclusionIndex[localGroupIndex] = -1;
                    for (int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
                        if (exclusionIndices[i] == y)
                            exclusionIndex[localGroupIndex] = i*TILE_SIZE;
                    computeSubset = (exclusionIndex[localGroupIndex] == -1);
                }
                if (computeSubset) {
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
                                        bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                         (params2.y*params2.y*invR)*(l_ij2-u_ij2));
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
                                        float term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                         (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                                        if (params2.x < params1.y-r)
                                            term += 2.0f*(RECIP(params2.x)-l_ij);
                                        tempBuffer[get_local_id(0)] = term;
                                    }
                                }

                                // Sum the forces on atom j.

                                if (tgx % 4 == 0)
                                    tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+1]+tempBuffer[get_local_id(0)+2]+tempBuffer[get_local_id(0)+3];
                                if (tgx == 0)
                                    localData[tbx+j].bornSum += tempBuffer[get_local_id(0)]+tempBuffer[get_local_id(0)+4]+tempBuffer[get_local_id(0)+8]+tempBuffer[get_local_id(0)+12]+tempBuffer[get_local_id(0)+16]+tempBuffer[get_local_id(0)+20]+tempBuffer[get_local_id(0)+24]+tempBuffer[get_local_id(0)+28];
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
                                bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                 (params2.y*params2.y*invR)*(l_ij2-u_ij2));
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
                                float term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                 (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                                if (params2.x < params1.y-r)
                                    term += 2.0f*(RECIP(params2.x)-l_ij);
                                localData[tbx+tj].bornSum += term;
                            }
                        }
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                }
            }
        }
        
        // Write results.  We need to coordinate between warps to make sure no two of them
        // ever try to write to the same piece of memory at the same time.
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atom_add(&global_bornSum[offset], (long) (bornSum*0xFFFFFFFF));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atom_add(&global_bornSum[offset], (long) (localData[get_local_id(0)].bornSum*0xFFFFFFFF));
        }
#else
        int writeX = (pos < end ? x : -1);
        int writeY = (pos < end && x != y ? y : -1);
        if (tgx == 0)
            reservedBlocks[localGroupIndex] = (int2)(writeX, writeY);
        bool done = false;
        int doneIndex = 0;
        int checkIndex = 0;
        while (true) {
            // See if any warp still needs to write its data.

            bool allDone = true;
            barrier(CLK_LOCAL_MEM_FENCE);
            while (doneIndex < WARPS_PER_GROUP && allDone) {
                if (reservedBlocks[doneIndex].x != -1)
                    allDone = false;
                else
                    doneIndex++;
            }
            if (allDone)
                break;
            if (!done) {
                // See whether this warp can write its data.  This requires that no previous warp
                // is trying to write to the same block of the buffer.

                bool canWrite = (writeX != -1);
                while (checkIndex < localGroupIndex && canWrite) {
                    if ((reservedBlocks[checkIndex].x == x || reservedBlocks[checkIndex].y == x) ||
                            (writeY != -1 && (reservedBlocks[checkIndex].x == y || reservedBlocks[checkIndex].y == y)))
                        canWrite = false;
                    else
                        checkIndex++;
                }
                if (canWrite) {
                    // Write the data to global memory, then mark this warp as done.

                    if (writeX > -1) {
                        const unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        global_bornSum[offset] += bornSum;
                    }
                    if (writeY > -1) {
                        const unsigned int offset = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        global_bornSum[offset] += localData[get_local_id(0)].bornSum;
                    }
                    done = true;
                    if (tgx == 0)
                        reservedBlocks[localGroupIndex] = (int2)(-1, -1);
                }
            }
        }
#endif
        lasty = y;
        pos++;
    } while (pos < end);
}

typedef struct {
    float x, y, z;
    float q;
    float fx, fy, fz, fw;
    float bornRadius;
} AtomData2;

/**
 * First part of computing the GBSA interaction.
 */

__kernel void computeGBSAForce1(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers, __global long* restrict global_bornForce,
#else
        __global float4* restrict forceBuffers, __global float* restrict global_bornForce,
#endif
        __global float* restrict energyBuffer, __global const float4* restrict posq, __global const float* restrict global_bornRadii,
        __local AtomData2* restrict localData,
#ifdef USE_CUTOFF
        __global const ushort2* restrict tiles, __global const unsigned int* restrict interactionCount, float4 periodicBoxSize, float4 invPeriodicBoxSize, unsigned int maxTiles, __global const unsigned int* restrict interactionFlags,
#else
        unsigned int numTiles,
#endif
        __global unsigned int* exclusionIndices, __global unsigned int* exclusionRowIndices) {
    unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    unsigned int warp = get_global_id(0)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    float energy = 0.0f;
    unsigned int lasty = 0xFFFFFFFF;
    __local float4 tempBuffer[FORCE_WORK_GROUP_SIZE];
    __local int2 reservedBlocks[WARPS_PER_GROUP];
    __local unsigned int* exclusionRange = (__local unsigned int*) reservedBlocks;
    __local int exclusionIndex[WARPS_PER_GROUP];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
        const unsigned int tbx = get_local_id(0) - tgx;
        const unsigned int localGroupIndex = get_local_id(0)/TILE_SIZE;
        unsigned int x, y;
        float4 force = 0.0f;
        if (pos < end) {
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
            unsigned int atom1 = x*TILE_SIZE + tgx;
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
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
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
                        float dEdR = Gpol*(1.0f - 0.25f*expTerm);
                        force.w += dGpol_dalpha2_ij*bornRadius2;
                        energy += 0.5f*tempEnergy;
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
#ifdef USE_CUTOFF
                        }
#endif
                    }
                }
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
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                bool computeSubset = false;
                if (flags != 0xFFFFFFFF) {
                    if (tgx < 2)
                        exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
                    if (tgx == 0)
                        exclusionIndex[localGroupIndex] = -1;
                    for (int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
                        if (exclusionIndices[i] == y)
                            exclusionIndex[localGroupIndex] = i*TILE_SIZE;
                    computeSubset = (exclusionIndex[localGroupIndex] == -1);
                }
                if (computeSubset) {
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
#ifdef USE_CUTOFF
                                if (r2 < CUTOFF_SQUARED) {
#endif
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
                                float dEdR = Gpol*(1.0f - 0.25f*expTerm);
#ifdef USE_CUTOFF
                                if (atom1 >= NUM_ATOMS || y*TILE_SIZE+j >= NUM_ATOMS || r2 > CUTOFF_SQUARED) {
#else
                                if (atom1 >= NUM_ATOMS || y*TILE_SIZE+j >= NUM_ATOMS) {
#endif
                                    dEdR = 0.0f;
                                    dGpol_dalpha2_ij = 0.0f;
                                    tempEnergy = 0.0f;
                                }
                                energy += tempEnergy;
                                force.w += dGpol_dalpha2_ij*bornRadius2;
                                delta.xyz *= dEdR;
                                force.xyz -= delta.xyz;
                                tempBuffer[get_local_id(0)] = (float4) (delta.xyz, dGpol_dalpha2_ij*bornRadius1);
#ifdef USE_CUTOFF
                                }
                                else
                                    tempBuffer[get_local_id(0)] = (float4) 0.0f;
#endif

                                // Sum the forces on atom j.

                                if (tgx % 4 == 0)
                                    tempBuffer[get_local_id(0)] += tempBuffer[get_local_id(0)+1]+tempBuffer[get_local_id(0)+2]+tempBuffer[get_local_id(0)+3];
                                if (tgx == 0) {
                                    float4 sum = tempBuffer[get_local_id(0)]+tempBuffer[get_local_id(0)+4]+tempBuffer[get_local_id(0)+8]+tempBuffer[get_local_id(0)+12]+tempBuffer[get_local_id(0)+16]+tempBuffer[get_local_id(0)+20]+tempBuffer[get_local_id(0)+24]+tempBuffer[get_local_id(0)+28];
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
#ifdef USE_CUTOFF
                            if (r2 < CUTOFF_SQUARED) {
#endif
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
                            float dEdR = Gpol*(1.0f - 0.25f*expTerm);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
                            energy += tempEnergy;
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            localData[tbx+tj].fx += delta.x;
                            localData[tbx+tj].fy += delta.y;
                            localData[tbx+tj].fz += delta.z;
                            localData[tbx+tj].fw += dGpol_dalpha2_ij*bornRadius1;
#ifdef USE_CUTOFF
                            }
#endif
                        }
                        tj = (tj + 1) & (TILE_SIZE - 1);
                    }
                }
            }
        }
        
        // Write results.  We need to coordinate between warps to make sure no two of them
        // ever try to write to the same piece of memory at the same time.
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (force.x*0xFFFFFFFF));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (force.y*0xFFFFFFFF));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (force.z*0xFFFFFFFF));
            atom_add(&global_bornForce[offset], (long) (force.w*0xFFFFFFFF));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (localData[get_local_id(0)].fx*0xFFFFFFFF));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fy*0xFFFFFFFF));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (localData[get_local_id(0)].fz*0xFFFFFFFF));
            atom_add(&global_bornForce[offset], (long) (localData[get_local_id(0)].fw*0xFFFFFFFF));
        }
#else
        int writeX = (pos < end ? x : -1);
        int writeY = (pos < end && x != y ? y : -1);
        if (tgx == 0)
            reservedBlocks[localGroupIndex] = (int2)(writeX, writeY);
        bool done = false;
        int doneIndex = 0;
        int checkIndex = 0;
        while (true) {
            // See if any warp still needs to write its data.

            bool allDone = true;
            barrier(CLK_LOCAL_MEM_FENCE);
            while (doneIndex < WARPS_PER_GROUP && allDone) {
                if (reservedBlocks[doneIndex].x != -1)
                    allDone = false;
                else
                    doneIndex++;
            }
            if (allDone)
                break;
            if (!done) {
                // See whether this warp can write its data.  This requires that no previous warp
                // is trying to write to the same block of the buffer.

                bool canWrite = (writeX != -1);
                while (checkIndex < localGroupIndex && canWrite) {
                    if ((reservedBlocks[checkIndex].x == x || reservedBlocks[checkIndex].y == x) ||
                            (writeY != -1 && (reservedBlocks[checkIndex].x == y || reservedBlocks[checkIndex].y == y)))
                        canWrite = false;
                    else
                        checkIndex++;
                }
                if (canWrite) {
                    // Write the data to global memory, then mark this warp as done.

                    if (writeX > -1) {
                        const unsigned int offset = x*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset].xyz += force.xyz;
                        global_bornForce[offset] += force.w;
                    }
                    if (writeY > -1) {
                        const unsigned int offset = y*TILE_SIZE + tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                        forceBuffers[offset] += (float4) (localData[get_local_id(0)].fx, localData[get_local_id(0)].fy, localData[get_local_id(0)].fz, 0.0f);
                        global_bornForce[offset] += localData[get_local_id(0)].fw;
                    }
                    done = true;
                    if (tgx == 0)
                        reservedBlocks[localGroupIndex] = (int2)(-1, -1);
                }
            }
        }
#endif
        lasty = y;
        pos++;
    } while (pos < end);
    energyBuffer[get_global_id(0)] += energy;
}
