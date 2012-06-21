#define TILE_SIZE 32
#define WARPS_PER_GROUP (FORCE_WORK_GROUP_SIZE/TILE_SIZE)

typedef struct {
    real x, y, z;
    real q;
    float radius, scaledRadius;
    real bornSum;
} AtomData1;

/**
 * Compute the Born sum.
 */
extern "C" __global__ void computeBornSum(unsigned long long* __restrict__ global_bornSum, const real4* __restrict__ posq, const float2* __restrict__ global_params,
#ifdef USE_CUTOFF
        const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags,
#else
        unsigned int numTiles,
#endif
        unsigned int* exclusionIndices, unsigned int* exclusionRowIndices) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    unsigned int lasty = 0xFFFFFFFF;
    __shared__ AtomData1 localData[FORCE_WORK_GROUP_SIZE];
    __shared__ real tempBuffer[FORCE_WORK_GROUP_SIZE];
    __shared__ unsigned int exclusionRange[2*WARPS_PER_GROUP];
    __shared__ int exclusionIndex[WARPS_PER_GROUP];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        const unsigned int localGroupIndex = threadIdx.x/TILE_SIZE;
        unsigned int x, y;
        real bornSum = 0;
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
            real4 posq1 = posq[atom1];
            float2 params1 = global_params[atom1];
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].x = posq1.x;
                localData[threadIdx.x].y = posq1.y;
                localData[threadIdx.x].z = posq1.z;
                localData[threadIdx.x].q = posq1.w;
                localData[threadIdx.x].radius = params1.x;
                localData[threadIdx.x].scaledRadius = params1.y;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 delta = make_real3(localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z);
#ifdef USE_PERIODIC
                    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                        real invR = RSQRT(r2);
                        real r = RECIP(invR);
                        float2 params2 = make_float2(localData[tbx+j].radius, localData[tbx+j].scaledRadius);
                        real rScaledRadiusJ = r+params2.y;
                        if ((j != tgx) && (params1.x < rScaledRadiusJ)) {
                            real l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                            real u_ij = RECIP(rScaledRadiusJ);
                            real l_ij2 = l_ij*l_ij;
                            real u_ij2 = u_ij*u_ij;
                            real ratio = LOG(u_ij * RECIP(l_ij));
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
                    real4 tempPosq = posq[j];
                    localData[threadIdx.x].x = tempPosq.x;
                    localData[threadIdx.x].y = tempPosq.y;
                    localData[threadIdx.x].z = tempPosq.z;
                    localData[threadIdx.x].q = tempPosq.w;
                    float2 tempParams = global_params[j];
                    localData[threadIdx.x].radius = tempParams.x;
                    localData[threadIdx.x].scaledRadius = tempParams.y;
                }
                localData[threadIdx.x].bornSum = 0.0f;
#ifdef USE_CUTOFF
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                bool computeSubset = false;
                if (flags != 0xFFFFFFFF) {
                    if (tgx < 2)
                        exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
                    if (tgx == 0)
                        exclusionIndex[localGroupIndex] = -1;
                    for (unsigned int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
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
                                real3 delta = make_real3(localData[tbx+j].x-posq1.x, localData[tbx+j].y-posq1.y, localData[tbx+j].z-posq1.z);
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                                tempBuffer[threadIdx.x] = 0.0f;
#ifdef USE_CUTOFF
                                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                                if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
#endif
                                    real invR = RSQRT(r2);
                                    real r = RECIP(invR);
                                    float2 params2 = make_float2(localData[tbx+j].radius, localData[tbx+j].scaledRadius);
                                    real rScaledRadiusJ = r+params2.y;
                                    if (params1.x < rScaledRadiusJ) {
                                        real l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                                        real u_ij = RECIP(rScaledRadiusJ);
                                        real l_ij2 = l_ij*l_ij;
                                        real u_ij2 = u_ij*u_ij;
                                        real ratio = LOG(u_ij * RECIP(l_ij));
                                        bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                         (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                                        if (params1.x < params2.y-r)
                                            bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                                    }
                                    real rScaledRadiusI = r+params1.y;
                                    if (params2.x < rScaledRadiusI) {
                                        real l_ij = RECIP(max(params2.x, fabs(r-params1.y)));
                                        real u_ij = RECIP(rScaledRadiusI);
                                        real l_ij2 = l_ij*l_ij;
                                        real u_ij2 = u_ij*u_ij;
                                        real ratio = LOG(u_ij * RECIP(l_ij));
                                        real term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                         (params1.y*params1.y*invR)*(l_ij2-u_ij2));
                                        if (params2.x < params1.y-r)
                                            term += 2.0f*(RECIP(params2.x)-l_ij);
                                        tempBuffer[threadIdx.x] = term;
                                    }
                                }

                                // Sum the forces on atom j.

                                if (tgx % 4 == 0)
                                    tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1]+tempBuffer[threadIdx.x+2]+tempBuffer[threadIdx.x+3];
                                if (tgx == 0)
                                    localData[tbx+j].bornSum += tempBuffer[threadIdx.x]+tempBuffer[threadIdx.x+4]+tempBuffer[threadIdx.x+8]+tempBuffer[threadIdx.x+12]+tempBuffer[threadIdx.x+16]+tempBuffer[threadIdx.x+20]+tempBuffer[threadIdx.x+24]+tempBuffer[threadIdx.x+28];
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
                        real3 delta = make_real3(localData[tbx+tj].x-posq1.x, localData[tbx+tj].y-posq1.y, localData[tbx+tj].z-posq1.z);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                        if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                        if (atom1 < NUM_ATOMS && y*TILE_SIZE+tj < NUM_ATOMS) {
#endif
                            real invR = RSQRT(r2);
                            real r = RECIP(invR);
                            float2 params2 = make_float2(localData[tbx+tj].radius, localData[tbx+tj].scaledRadius);
                            real rScaledRadiusJ = r+params2.y;
                            if (params1.x < rScaledRadiusJ) {
                                real l_ij = RECIP(max(params1.x, fabs(r-params2.y)));
                                real u_ij = RECIP(rScaledRadiusJ);
                                real l_ij2 = l_ij*l_ij;
                                real u_ij2 = u_ij*u_ij;
                                real ratio = LOG(u_ij * RECIP(l_ij));
                                bornSum += l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
                                                 (params2.y*params2.y*invR)*(l_ij2-u_ij2));
                                if (params1.x < params2.y-r)
                                    bornSum += 2.0f*(RECIP(params1.x)-l_ij);
                            }
                            real rScaledRadiusI = r+params1.y;
                            if (params2.x < rScaledRadiusI) {
                                real l_ij = RECIP(max(params2.x, fabs(r-params1.y)));
                                real u_ij = RECIP(rScaledRadiusI);
                                real l_ij2 = l_ij*l_ij;
                                real u_ij2 = u_ij*u_ij;
                                real ratio = LOG(u_ij * RECIP(l_ij));
                                real term = l_ij - u_ij + (0.50f*invR*ratio) + 0.25f*(r*(u_ij2-l_ij2) +
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
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&global_bornSum[offset], static_cast<unsigned long long>((long long) (bornSum*0xFFFFFFFF)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atomicAdd(&global_bornSum[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].bornSum*0xFFFFFFFF)));
        }
        lasty = y;
        pos++;
    } while (pos < end);
}

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz, fw;
    real bornRadius;
} AtomData2;

/**
 * First part of computing the GBSA interaction.
 */

extern "C" __global__ void computeGBSAForce1(unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ global_bornForce,
        real* __restrict__ energyBuffer, const real4* __restrict__ posq, const real* __restrict__ global_bornRadii,
#ifdef USE_CUTOFF
        const ushort2* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const unsigned int* __restrict__ interactionFlags,
#else
        unsigned int numTiles,
#endif
        unsigned int* exclusionIndices, unsigned int* exclusionRowIndices) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    unsigned int pos = warp*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
    unsigned int end = (warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*(NUM_BLOCKS+1)/2 : numTiles)/totalWarps;
#else
    unsigned int pos = warp*numTiles/totalWarps;
    unsigned int end = (warp+1)*numTiles/totalWarps;
#endif
    real energy = 0;
    unsigned int lasty = 0xFFFFFFFF;
    __shared__ AtomData2 localData[FORCE_WORK_GROUP_SIZE];
    __shared__ real4 tempBuffer[FORCE_WORK_GROUP_SIZE];
    __shared__ unsigned int exclusionRange[2*WARPS_PER_GROUP];
    __shared__ int exclusionIndex[WARPS_PER_GROUP];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        const unsigned int localGroupIndex = threadIdx.x/TILE_SIZE;
        unsigned int x, y;
        real4 force = make_real4(0);
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
            real4 posq1 = posq[atom1];
            real bornRadius1 = global_bornRadii[atom1];
            if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].x = posq1.x;
                localData[threadIdx.x].y = posq1.y;
                localData[threadIdx.x].z = posq1.z;
                localData[threadIdx.x].q = posq1.w;
                localData[threadIdx.x].bornRadius = bornRadius1;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    if (atom1 < NUM_ATOMS && y*TILE_SIZE+j < NUM_ATOMS) {
                        real4 posq2 = make_real4(localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z, localData[tbx+j].q);
                        real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                        delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                        delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                        delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                        if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = RECIP(invR);
                        real bornRadius2 = localData[tbx+j].bornRadius;
                        real alpha2_ij = bornRadius1*bornRadius2;
                        real D_ij = r2*RECIP(4.0f*alpha2_ij);
                        real expTerm = EXP(-D_ij);
                        real denominator2 = r2 + alpha2_ij*expTerm;
                        real denominator = SQRT(denominator2);
                        real tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                        real Gpol = tempEnergy*RECIP(denominator2);
                        real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                        real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                        force.w += dGpol_dalpha2_ij*bornRadius2;
                        energy += 0.5f*tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
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
                    real4 tempPosq = posq[j];
                    localData[threadIdx.x].x = tempPosq.x;
                    localData[threadIdx.x].y = tempPosq.y;
                    localData[threadIdx.x].z = tempPosq.z;
                    localData[threadIdx.x].q = tempPosq.w;
                    localData[threadIdx.x].bornRadius = global_bornRadii[j];
                }
                localData[threadIdx.x].fx = 0.0f;
                localData[threadIdx.x].fy = 0.0f;
                localData[threadIdx.x].fz = 0.0f;
                localData[threadIdx.x].fw = 0.0f;
#ifdef USE_CUTOFF
                unsigned int flags = (numTiles <= maxTiles ? interactionFlags[pos] : 0xFFFFFFFF);
                bool computeSubset = false;
                if (flags != 0xFFFFFFFF) {
                    if (tgx < 2)
                        exclusionRange[2*localGroupIndex+tgx] = exclusionRowIndices[x+tgx];
                    if (tgx == 0)
                        exclusionIndex[localGroupIndex] = -1;
                    for (unsigned int i = exclusionRange[2*localGroupIndex]+tgx; i < exclusionRange[2*localGroupIndex+1]; i += TILE_SIZE)
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
                                real4 posq2 = make_real4(localData[tbx+j].x, localData[tbx+j].y, localData[tbx+j].z, localData[tbx+j].q);
                                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                                delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                                delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                                delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                                if (r2 < CUTOFF_SQUARED) {
#endif
                                real invR = RSQRT(r2);
                                real r = RECIP(invR);
                                real bornRadius2 = localData[tbx+j].bornRadius;
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
                                delta *= dEdR;
                                force.x -= delta.x;
                                force.y -= delta.y;
                                force.z -= delta.z;
                                tempBuffer[threadIdx.x] = make_real4(delta.x, delta.y, delta.z, dGpol_dalpha2_ij*bornRadius1);
#ifdef USE_CUTOFF
                                }
                                else
                                    tempBuffer[threadIdx.x] = make_real4(0);
#endif

                                // Sum the forces on atom j.

                                if (tgx % 4 == 0)
                                    tempBuffer[threadIdx.x] += tempBuffer[threadIdx.x+1]+tempBuffer[threadIdx.x+2]+tempBuffer[threadIdx.x+3];
                                if (tgx == 0) {
                                    real4 sum = tempBuffer[threadIdx.x]+tempBuffer[threadIdx.x+4]+tempBuffer[threadIdx.x+8]+tempBuffer[threadIdx.x+12]+tempBuffer[threadIdx.x+16]+tempBuffer[threadIdx.x+20]+tempBuffer[threadIdx.x+24]+tempBuffer[threadIdx.x+28];
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
                            real4 posq2 = make_real4(localData[tbx+tj].x, localData[tbx+tj].y, localData[tbx+tj].z, localData[tbx+tj].q);
                            real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                            delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
                            delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
                            delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;
#endif
                            real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                            if (r2 < CUTOFF_SQUARED) {
#endif
                            real invR = RSQRT(r2);
                            real r = RECIP(invR);
                            real bornRadius2 = localData[tbx+tj].bornRadius;
                            real alpha2_ij = bornRadius1*bornRadius2;
                            real D_ij = r2*RECIP(4.0f*alpha2_ij);
                            real expTerm = EXP(-D_ij);
                            real denominator2 = r2 + alpha2_ij*expTerm;
                            real denominator = SQRT(denominator2);
                            real tempEnergy = (PREFACTOR*posq1.w*posq2.w)*RECIP(denominator);
                            real Gpol = tempEnergy*RECIP(denominator2);
                            real dGpol_dalpha2_ij = -0.5f*Gpol*expTerm*(1.0f+D_ij);
                            real dEdR = Gpol*(1.0f - 0.25f*expTerm);
                            force.w += dGpol_dalpha2_ij*bornRadius2;
                            energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
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
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (force.x*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0xFFFFFFFF)));
            atomicAdd(&global_bornForce[offset], static_cast<unsigned long long>((long long) (force.w*0xFFFFFFFF)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0xFFFFFFFF)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0xFFFFFFFF)));
            atomicAdd(&global_bornForce[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fw*0xFFFFFFFF)));
        }
        lasty = y;
        pos++;
    } while (pos < end);
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
