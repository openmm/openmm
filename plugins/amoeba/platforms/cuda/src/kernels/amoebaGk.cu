#define TILE_SIZE 32

/**
 * Reduce the Born sums to compute the Born radii.
 */
extern "C" __global__ void reduceBornSum(const long long* __restrict__ bornSum, const float2* __restrict__ params, real* __restrict__ bornRadii) {
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Get summed Born data

        real sum = RECIP(0x100000000)*bornSum[index];

        // Now calculate Born radius.

        float radius = params[index].x;
        radius = RECIP(radius*radius*radius);
        sum = radius-sum;
        sum = (sum <= 0 ? (real) 1000 : POW(sum, -1/(real) 3));
        bornRadii[index] = sum;
    }
}

#ifdef SURFACE_AREA_FACTOR
/**
 * Apply the surface area term to the force and energy.
 */
extern "C" __global__ void computeSurfaceAreaForce(long long* __restrict__ bornForce, mixed* __restrict__ energyBuffer, const float2* __restrict__ params, const real* __restrict__ bornRadii) {
    mixed energy = 0;
    for (unsigned int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real bornRadius = bornRadii[index];
        float radius = params[index].x;
        real r = radius + DIELECTRIC_OFFSET + PROBE_RADIUS;
        real ratio6 = (radius+DIELECTRIC_OFFSET)/bornRadius;
        ratio6 = ratio6*ratio6*ratio6;
        ratio6 = ratio6*ratio6;
        real saTerm = SURFACE_AREA_FACTOR * r * r * ratio6;
        bornForce[index] += (long long) (saTerm*0x100000000/bornRadius);
        energy += saTerm;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] -= energy/6;
}
#endif

/**
 * Data structure used by computeBornSum().
 */
typedef struct {
    real3 pos;
    real bornSum;
    float radius, scaledRadius, padding;
} AtomData1;

__device__ void computeBornSumOneInteraction(AtomData1& atom1, AtomData1& atom2) {
    if (atom1.radius <= 0)
        return; // Ignore this interaction
    real3 delta = atom2.pos - atom1.pos;
    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    float sk = atom2.scaledRadius;

    if (atom1.radius > r + sk)
        return; // No descreening due to atom1 engulfing atom2.

    real sk2 = sk*sk;
    if (atom1.radius+r < sk) {
        real lik = atom1.radius;
        real uik = sk - r; 
        atom1.bornSum -= RECIP(uik*uik*uik) - RECIP(lik*lik*lik);
    }
    real uik = r+sk;
    real lik;
    if (atom1.radius+r < sk)
        lik = sk-r;
    else if (r < atom1.radius+sk)
        lik = atom1.radius;
    else
        lik = r-sk;
    real l2 = lik*lik; 
    real l4 = l2*l2;
    real lr = lik*r;
    real l4r = l4*r; 
    real u2 = uik*uik;
    real u4 = u2*u2;
    real ur = uik*r; 
    real u4r = u4*r;
    real term = (3*(r2-sk2)+6*u2-8*ur)/u4r - (3*(r2-sk2)+6*l2-8*lr)/l4r;
    atom1.bornSum += term/16;
}

/**
 * Compute the Born sum.
 */
extern "C" __global__ void computeBornSum(unsigned long long* __restrict__ bornSum, const real4* __restrict__ posq,
        const float2* __restrict__ params, unsigned int numTiles) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    unsigned int pos = (unsigned int) (warp*(long long)numTiles/totalWarps);
    unsigned int end = (unsigned int) ((warp+1)*(long long)numTiles/totalWarps);
    unsigned int lasty = 0xFFFFFFFF;
    __shared__ AtomData1 localData[BORN_SUM_THREAD_BLOCK_SIZE];
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        int x, y;
        AtomData1 data;
        data.bornSum = 0;
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-sqrt((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            data.pos = trimTo3(posq[atom1]);
            float2 params1 = params[atom1];
            data.radius = params1.x;
            data.scaledRadius = params1.y;
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].pos = data.pos;
                localData[threadIdx.x].radius = params1.x;
                localData[threadIdx.x].scaledRadius = params1.y;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2)
                        computeBornSumOneInteraction(data, localData[tbx+j]);
                }
            }
            else {
                // This is an off-diagonal tile.

                if (lasty != y) {
                    unsigned int j = y*TILE_SIZE + tgx;
                    real4 tempPosq = posq[j];
                    localData[threadIdx.x].pos = trimTo3(tempPosq);
                    float2 tempParams = params[j];
                    localData[threadIdx.x].radius = tempParams.x;
                    localData[threadIdx.x].scaledRadius = tempParams.y;
                }
                localData[threadIdx.x].bornSum = 0;
                
                // Compute the full set of interactions in this tile.

                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        computeBornSumOneInteraction(data, localData[tbx+tj]);
                        computeBornSumOneInteraction(localData[tbx+tj], data);
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }
        }
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&bornSum[offset], static_cast<unsigned long long>((long long) (data.bornSum*0x100000000)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            atomicAdd(&bornSum[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].bornSum*0x100000000)));
        }
        lasty = y;
        pos++;
    } while (pos < end);
}

/**
 * Data structure used by computeGKForces().
 */
typedef struct {
    real3 pos, force, dipole, inducedDipole, inducedDipolePolar;
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ, quadrupoleZZ;
    real q, bornRadius, bornForce;
} AtomData2;

__device__ void computeOneInteractionF1(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force);
__device__ void computeOneInteractionF2(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force);
__device__ void computeOneInteractionT1(AtomData2& atom1, volatile AtomData2& atom2, real3& torque);
__device__ void computeOneInteractionT2(AtomData2& atom1, volatile AtomData2& atom2, real3& torque);
__device__ void computeOneInteractionB1B2(AtomData2& atom1, volatile AtomData2& atom2);

inline __device__ void loadAtomData2(AtomData2& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const real* __restrict__ bornRadius) {
    real4 atomPosq = posq[atom];
    data.pos = trimTo3(atomPosq);
    data.q = atomPosq.w;
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.quadrupoleZZ = -(data.quadrupoleXX+data.quadrupoleYY);
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    data.bornRadius = bornRadius[atom];
}

inline __device__ void zeroAtomData(AtomData2& data) {
    data.force = make_real3(0);
    data.bornForce = 0;
}

/**
 * Compute electrostatic interactions.
 */
extern "C" __global__ void computeGKForces(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ torqueBuffers, mixed* __restrict__ energyBuffer,
        const real4* __restrict__ posq, unsigned int startTileIndex, unsigned int numTileIndices, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
        const real* __restrict__ bornRadii, unsigned long long* __restrict__ bornForce) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = (unsigned int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    unsigned int end = (unsigned int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
    mixed energy = 0;
    __shared__ AtomData2 localData[GK_FORCE_THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        int x, y;
        AtomData2 data;
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            loadAtomData2(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, bornRadii);
            zeroAtomData(data);
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].pos = data.pos;
                localData[threadIdx.x].q = data.q;
                localData[threadIdx.x].dipole = data.dipole;
                localData[threadIdx.x].quadrupoleXX = data.quadrupoleXX;
                localData[threadIdx.x].quadrupoleXY = data.quadrupoleXY;
                localData[threadIdx.x].quadrupoleXZ = data.quadrupoleXZ;
                localData[threadIdx.x].quadrupoleYY = data.quadrupoleYY;
                localData[threadIdx.x].quadrupoleYZ = data.quadrupoleYZ;
                localData[threadIdx.x].quadrupoleZZ = data.quadrupoleZZ;
                localData[threadIdx.x].inducedDipole = data.inducedDipole;
                localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
                localData[threadIdx.x].bornRadius = data.bornRadius;
                
                // Compute forces.
                
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        real tempEnergy;
                        computeOneInteractionF1(data, localData[tbx+j], tempEnergy, tempForce);
                        computeOneInteractionF2(data, localData[tbx+j], tempEnergy, tempForce);
                        data.force += tempForce;
                        energy += 0.5f*tempEnergy;
                    }
                }
                data.force *= 0.5f;
                atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                
                // Compute torques.
                
                zeroAtomData(data);
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempTorque;
                        computeOneInteractionT1(data, localData[tbx+j], tempTorque);
                        computeOneInteractionT2(data, localData[tbx+j], tempTorque);
                        data.force += tempTorque;
                    }
                }
                atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                
                // Compute chain rule terms.
                
                zeroAtomData(data);
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                        computeOneInteractionB1B2(data, localData[tbx+j]);
                }
                atomicAdd(&bornForce[atom1], static_cast<unsigned long long>((long long) (data.bornForce*0x100000000)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData2(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, bornRadii);
                zeroAtomData(localData[threadIdx.x]);
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        real tempEnergy;
                        computeOneInteractionF1(data, localData[tbx+tj], tempEnergy, tempForce);
                        computeOneInteractionF2(data, localData[tbx+tj], tempEnergy, tempForce);
                        data.force += tempForce;
                        localData[tbx+tj].force -= tempForce;
                        energy += tempEnergy;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
                data.force *= 0.5f;
                localData[threadIdx.x].force *= 0.5f;
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                    atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                    atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                    atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                    atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                }

                // Compute torques.

                zeroAtomData(data);
                zeroAtomData(localData[threadIdx.x]);
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempTorque;
                        computeOneInteractionT1(data, localData[tbx+tj], tempTorque);
                        computeOneInteractionT2(data, localData[tbx+tj], tempTorque);
                        data.force += tempTorque;
                        computeOneInteractionT1(localData[tbx+tj], data, tempTorque);
                        computeOneInteractionT2(localData[tbx+tj], data, tempTorque);
                        localData[tbx+tj].force += tempTorque;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                    atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                    atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                    atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                    atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                }

                // Compute chain rule terms.

                zeroAtomData(data);
                zeroAtomData(localData[threadIdx.x]);
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                        computeOneInteractionB1B2(data, localData[tbx+tj]);
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    atomicAdd(&bornForce[offset], static_cast<unsigned long long>((long long) (data.bornForce*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    atomicAdd(&bornForce[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].bornForce*0x100000000)));
                }
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*0.5f;
}


/**
 * Data structure used by computeChainRuleForce().
 */
typedef struct {
    real3 pos, force;
    real radius, scaledRadius, bornRadius, bornForce;
} AtomData3;

inline __device__ void loadAtomData3(AtomData3& data, int atom, const real4* __restrict__ posq, const float2* __restrict__ params, const real* __restrict__ bornRadius, const long long* __restrict__ bornForce) {
    data.pos = trimTo3(posq[atom]);
    data.bornRadius = bornRadius[atom];
    float2 params1 = params[atom];
    data.radius = params1.x;
    data.scaledRadius = params1.y;
    data.bornForce = bornForce[atom]/(real) 0x100000000;
}

__device__ void computeBornChainRuleInteraction(AtomData3& atom1, AtomData3& atom2, real3& force) {
    real third = 1/(real) 3;
    real pi43 = 4*third*M_PI;
    real factor = -POW(M_PI, third)*POW((real) 6, 2/(real) 3)/9;
    real term = pi43/(atom1.bornRadius*atom1.bornRadius*atom1.bornRadius);
    term = factor/POW(term, 4/(real) 3);

    real3 delta = atom2.pos-atom1.pos;

    float sk = atom2.scaledRadius;
    real sk2 = sk*sk;
    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    real de = 0;

    if (atom1.radius > r + sk)
        return; // No descreening due to atom1 engulfing atom2.

    if (atom1.radius+r < sk) {
        real uik = sk-r;
        real uik4 = uik*uik;
        uik4 = uik4*uik4;
        de = -4*M_PI/uik4;
        real lik = sk - r;
        real lik4 = lik*lik;
        lik4 = lik4*lik4;
        de += 0.25f*M_PI*(sk2-4*sk*r+17*r2)/(r2*lik4);
    }
    else if (r < atom1.radius+sk) {
        real lik = atom1.radius;
        real lik4 = lik*lik;
        lik4 = lik4*lik4;
        de += 0.25f*M_PI*(2*atom1.radius*atom1.radius-sk2-r2)/(r2*lik4);
    }
    else {
        real lik = r-sk;
        real lik4 = lik*lik;
        lik4 = lik4*lik4;
        de += 0.25f*M_PI*(sk2-4*sk*r+r2)/(r2*lik4);
    }
    real uik = r+sk;
    real uik4 = uik*uik;
    uik4 = uik4*uik4;
    de -= 0.25f*M_PI*(sk2+4*sk*r+r2)/(r2*uik4);
    real dbr = term*de/r;
    de = dbr*atom1.bornForce;
    force = delta*de;
}

/**
 * Compute chain rule terms.
 */
extern "C" __global__ void computeChainRuleForce(
        unsigned long long* __restrict__ forceBuffers, const real4* __restrict__ posq, unsigned int startTileIndex, unsigned int numTileIndices,
        const float2* __restrict__ params, const real* __restrict__ bornRadii, const long long* __restrict__ bornForce) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = startTileIndex+warp*numTiles/totalWarps;
    unsigned int end = startTileIndex+(warp+1)*numTiles/totalWarps;
    __shared__ AtomData3 localData[CHAIN_RULE_THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        int x, y;
        AtomData3 data;
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            loadAtomData3(data, atom1, posq, params, bornRadii, bornForce);
            data.force = make_real3(0);
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.

                localData[threadIdx.x].pos = data.pos;
                localData[threadIdx.x].radius = data.radius;
                localData[threadIdx.x].scaledRadius = data.scaledRadius;
                localData[threadIdx.x].bornRadius = data.bornRadius;
                localData[threadIdx.x].bornForce = data.bornForce;
                localData[threadIdx.x].force = make_real3(0);
                
                // Compute forces.
                
                for (unsigned int j = (tgx+1)&(TILE_SIZE-1); j != tgx; j = (j+1)&(TILE_SIZE-1)) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        computeBornChainRuleInteraction(data, localData[tbx+j], tempForce);
                        data.force -= tempForce;
                        localData[tbx+j].force += tempForce;
                    }
                }
                atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) ((data.force.x+localData[threadIdx.x].force.x)*0x100000000)));
                atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((data.force.y+localData[threadIdx.x].force.y)*0x100000000)));
                atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) ((data.force.z+localData[threadIdx.x].force.z)*0x100000000)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData3(localData[threadIdx.x], j, posq, params, bornRadii, bornForce);
                localData[threadIdx.x].force = make_real3(0);
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        computeBornChainRuleInteraction(data, localData[tbx+tj], tempForce);
                        data.force -= tempForce;
                        localData[tbx+tj].force += tempForce;
                        computeBornChainRuleInteraction(localData[tbx+tj], data, tempForce);
                        data.force += tempForce;
                        localData[tbx+tj].force -= tempForce;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
                    atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
                    atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                    atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                    atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
                }
            }
        }
        pos++;
    } while (pos < end);
}

typedef struct {
    real3 pos, force, dipole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS;
    real q, quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ, quadrupoleZZ;
    float thole, damp;
} AtomData4;

__device__ void computeOneEDiffInteractionF1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real& outputEnergy, real3& outputForce);
__device__ void computeOneEDiffInteractionT1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce);
__device__ void computeOneEDiffInteractionT3(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce);

inline __device__ void loadAtomData4(AtomData4& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
        const real* __restrict__ inducedDipoleS, const real* __restrict__ inducedDipolePolarS, const float2* __restrict__ dampingAndThole) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.q = atomPosq.w;
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.quadrupoleZZ = -(data.quadrupoleXX+data.quadrupoleYY);
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    data.inducedDipoleS.x = inducedDipoleS[atom*3];
    data.inducedDipoleS.y = inducedDipoleS[atom*3+1];
    data.inducedDipoleS.z = inducedDipoleS[atom*3+2];
    data.inducedDipolePolarS.x = inducedDipolePolarS[atom*3];
    data.inducedDipolePolarS.y = inducedDipolePolarS[atom*3+1];
    data.inducedDipolePolarS.z = inducedDipolePolarS[atom*3+2];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

__device__ real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

/**
 * Compute electrostatic interactions.
 */
extern "C" __global__ void computeEDiffForce(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ torqueBuffers, mixed* __restrict__ energyBuffer,
        const real4* __restrict__ posq, const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags,
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const real* __restrict__ inducedDipoleS, const real* __restrict__ inducedDipolePolarS,
        const float2* __restrict__ dampingAndThole) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    mixed energy = 0;
    __shared__ AtomData4 localData[EDIFF_THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData4 data;
        data.force = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData4(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
        uint2 covalent = covalentFlags[pos*TILE_SIZE+tgx];
        unsigned int polarizationGroup = polarizationGroupFlags[pos*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].q = data.q;
            localData[threadIdx.x].dipole = data.dipole;
            localData[threadIdx.x].quadrupoleXX = data.quadrupoleXX;
            localData[threadIdx.x].quadrupoleXY = data.quadrupoleXY;
            localData[threadIdx.x].quadrupoleXZ = data.quadrupoleXZ;
            localData[threadIdx.x].quadrupoleYY = data.quadrupoleYY;
            localData[threadIdx.x].quadrupoleYZ = data.quadrupoleYZ;
            localData[threadIdx.x].quadrupoleZZ = data.quadrupoleZZ;
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].inducedDipoleS = data.inducedDipoleS;
            localData[threadIdx.x].inducedDipolePolarS = data.inducedDipolePolarS;
            localData[threadIdx.x].thole = data.thole;
            localData[threadIdx.x].damp = data.damp;

            // Compute forces.

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneEDiffInteractionF1(data, localData[tbx+j], d, p, tempEnergy, tempForce);
                    energy += 0.25f*tempEnergy;
                    data.force += tempForce;
                }
            }
            data.force *= ENERGY_SCALE_FACTOR;
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));

            // Compute torques.

            data.force = make_real3(0);
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneEDiffInteractionT1(data, localData[tbx+j], d, p, tempTorque);
                    data.force += tempTorque;
                }
            }
            data.force *= ENERGY_SCALE_FACTOR;
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData4(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[threadIdx.x].force = make_real3(0);

            // Compute forces.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneEDiffInteractionF1(data, localData[tbx+tj], d, p, tempEnergy, tempForce);
                    energy += 0.5f*tempEnergy;
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));

            // Compute torques.

            data.force = make_real3(0);
            localData[threadIdx.x].force = make_real3(0);
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneEDiffInteractionT1(data, localData[tbx+tj], d, p, tempTorque);
                    data.force += tempTorque;
                    computeOneEDiffInteractionT3(data, localData[tbx+tj], d, p, tempTorque);
                    localData[tbx+tj].force += tempTorque;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions (by enumerating all of them, since there's no cutoff).

    const unsigned int numTiles = numTileIndices;
    int pos = startTileIndex+warp*numTiles/totalWarps;
    int end = startTileIndex+(warp+1)*numTiles/totalWarps;
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ volatile int skipTiles[EDIFF_THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;

    while (pos < end) {
        // Extract the coordinates of this tile.

        int x, y;
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[threadIdx.x] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        bool includeTile = (skipTiles[currentSkipIndex] != pos);
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData4 data;
            data.force = make_real3(0);
            loadAtomData4(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            loadAtomData4(localData[threadIdx.x], atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData4(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[threadIdx.x].force = make_real3(0);

            // Compute forces.

            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    computeOneEDiffInteractionF1(data, localData[tbx+tj], 1, 1, tempEnergy, tempForce);
                    energy += 0.5f*tempEnergy;
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));

            // Compute torques.

            data.force = make_real3(0);
            localData[threadIdx.x].force = make_real3(0);
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    computeOneEDiffInteractionT1(data, localData[tbx+tj], 1, 1, tempTorque);
                    data.force += tempTorque;
                    computeOneEDiffInteractionT3(data, localData[tbx+tj], 1, 1, tempTorque);
                    localData[tbx+tj].force += tempTorque;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
        }
        pos++;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*ENERGY_SCALE_FACTOR;
}
