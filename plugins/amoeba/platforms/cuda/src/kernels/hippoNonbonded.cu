// This is a modified version of the standard nonbonded kernel for computing HippoNonbondedForce.
// This is needed because of two ways in which it differs from most nonbonded interactions:
// the force between two atoms doesn't always point along the line between them, and we need
// to accumulate torques as well as forces.

#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

#ifndef ENABLE_SHUFFLE
typedef struct {
    real x, y, z;
    real q;
#ifndef ENABLE_SHUFFLE
    real fx, fy, fz;
    real tx, ty, tz;
#endif
    ATOM_PARAMETER_DATA
#ifndef PARAMETER_SIZE_IS_EVEN
    real padding;
#endif
} AtomData;
#endif

#ifdef ENABLE_SHUFFLE
//support for 64 bit shuffles
static __inline__ __device__ float real_shfl(float var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ double real_shfl(double var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "d"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    return __hiloint2double( hi, lo );
}

static __inline__ __device__ long long real_shfl(long long var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "l"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    // unforunately there isn't an __nv_hiloint2long(hi,lo) intrinsic cast
    int2 fuse; fuse.x = lo; fuse.y = hi;
    return *reinterpret_cast<long long*>(&fuse);
}
#endif

__device__ void formQIRotationMatrix(real3 deltaR, real rInv, real (&rotationMatrix)[3][3]) {
    real3 vectorZ = deltaR*rInv;
    real3 vectorX = make_real3(0);
    if (fabs(vectorZ.y) > fabs(vectorZ.x))
        vectorX.x = 1;
    else
        vectorX.y = 1;

    vectorX -= vectorZ*dot(vectorZ, vectorX);
    vectorX = normalize(vectorX);
    real3 vectorY = cross(vectorZ, vectorX);

    rotationMatrix[0][0] = vectorX.x;
    rotationMatrix[0][1] = vectorX.y;
    rotationMatrix[0][2] = vectorX.z;
    rotationMatrix[1][0] = vectorY.x;
    rotationMatrix[1][1] = vectorY.y;
    rotationMatrix[1][2] = vectorY.z;
    rotationMatrix[2][0] = vectorZ.x;
    rotationMatrix[2][1] = vectorZ.y;
    rotationMatrix[2][2] = vectorZ.z;
}

__device__ real3 rotateVectorToQI(real3 v, const real (&mat)[3][3]) {
    return make_real3(mat[0][0]*v.x + mat[0][1]*v.y + mat[0][2]*v.z,
                      mat[1][0]*v.x + mat[1][1]*v.y + mat[1][2]*v.z,
                      mat[2][0]*v.x + mat[2][1]*v.y + mat[2][2]*v.z);
}

__device__ real3 rotateVectorFromQI(real3 v, const real (&mat)[3][3]) {
    return make_real3(mat[0][0]*v.x + mat[1][0]*v.y + mat[2][0]*v.z,
                      mat[0][1]*v.x + mat[1][1]*v.y + mat[2][1]*v.z,
                      mat[0][2]*v.x + mat[1][2]*v.y + mat[2][2]*v.z);
}

__device__ void rotateQuadrupoleToQI(real qXX, real qXY, real qXZ, real qYY, real qYZ, 
            real &qiQXX, real &qiQXY, real &qiQXZ, real &qiQYY, real &qiQYZ, real &qiQZZ, const real (&mat)[3][3]) {
    real qZZ = -qXX-qYY;
    qiQXX = mat[0][0]*(mat[0][0]*qXX + 2*(mat[0][1]*qXY + mat[0][2]*qXZ)) + mat[0][1]*(mat[0][1]*qYY + 2*mat[0][2]*qYZ) + mat[0][2]*mat[0][2]*qZZ;
    qiQYY = mat[1][0]*(mat[1][0]*qXX + 2*(mat[1][1]*qXY + mat[1][2]*qXZ)) + mat[1][1]*(mat[1][1]*qYY + 2*mat[1][2]*qYZ) + mat[1][2]*mat[1][2]*qZZ;
    qiQXY = mat[0][0]*mat[1][0]*qXX + mat[0][1]*mat[1][1]*qYY + mat[0][2]*mat[1][2]*qZZ + (mat[0][0]*mat[1][1] + mat[0][1]*mat[1][0])*qXY + (mat[0][0]*mat[1][2] + mat[0][2]*mat[1][0])*qXZ + (mat[0][1]*mat[1][2] + mat[0][2]*mat[1][1])*qYZ;
    qiQXZ = mat[0][0]*mat[2][0]*qXX + mat[0][1]*mat[2][1]*qYY + mat[0][2]*mat[2][2]*qZZ + (mat[0][0]*mat[2][1] + mat[0][1]*mat[2][0])*qXY + (mat[0][0]*mat[2][2] + mat[0][2]*mat[2][0])*qXZ + (mat[0][1]*mat[2][2] + mat[0][2]*mat[2][1])*qYZ;
    qiQYZ = mat[1][0]*mat[2][0]*qXX + mat[1][1]*mat[2][1]*qYY + mat[1][2]*mat[2][2]*qZZ + (mat[1][0]*mat[2][1] + mat[1][1]*mat[2][0])*qXY + (mat[1][0]*mat[2][2] + mat[1][2]*mat[2][0])*qXZ + (mat[1][1]*mat[2][2] + mat[1][2]*mat[2][1])*qYZ;
    qiQZZ = -qiQXX-qiQYY;
}

__device__ void computeDispersionDampingFactors(float alphaI, float alphaJ, real r, real& fdamp, real& ddamp) {
    real arI = alphaI*r;
    real arI2 = arI*arI;
    real arI3 = arI2*arI;
    real expARI = EXP(-arI);
    real fdamp3, fdamp5;
    real one = 1;
    real seven = 7;
    if (alphaI == alphaJ) {
        real arI4 = arI3*arI;
        real arI5 = arI4*arI;
        fdamp3 = 1 - (1 + arI + arI2*(one/2) + arI3*(seven/48) + arI4*(one/48))*expARI;
        fdamp5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/144))*expARI;
        ddamp = alphaI*(arI5 - 3*arI3 - 3*arI2)*expARI*(one/96);
    }
    else {
        real arJ = alphaJ*r;
        real arJ2 = arJ*arJ;
        real arJ3 = arJ2*arJ;
        real expARJ = EXP(-arJ);
        real aI2 = alphaI*alphaI;
        real aJ2 = alphaJ*alphaJ;
        real A = aJ2/(aJ2-aI2);
        real B = aI2/(aI2-aJ2);
        real A2 = A*A;
        real B2 = B*B;
        fdamp3 = 1 - A2*(1 + arI + arI2*(one/2))*expARI -
                     B2*(1 + arJ + arJ2*(one/2))*expARJ -
                     2*A2*B*(1 + arI)*expARI -
                     2*B2*A*(1 + arJ)*expARJ;
        fdamp5 = 1 - A2*(1 + arI + arI2*(one/2) + arI3*(one/6))*expARI -
                     B2*(1 + arJ + arJ2*(one/2) + arJ3*(one/6))*expARJ -
                     2*A2*B*(1 + arI + arI2*(one/3))*expARI -
                     2*B2*A*(1 + arJ + arJ2*(one/3))*expARJ;
        ddamp = (A2*arI2*alphaI*expARI*(r*alphaI + 4*B - 1) +
                (B2*arJ2*alphaJ*expARJ*(r*alphaJ + 4*A - 1)))*(one/4);
    }
    fdamp = 1.5f*fdamp5 - 0.5f*fdamp3;
}

__device__ void computeRepulsionDampingFactors(real pauliAlphaI, real pauliAlphaJ, real r,
            real& fdamp1, real& fdamp3, real& fdamp5, real& fdamp7, real& fdamp9, real& fdamp11) {
    real r2 = r*r;
    real r3 = r2*r;
    real r4 = r2*r2;
    real r5 = r3*r2;
    real r6 = r3*r3;
    real aI2 = 0.5f*pauliAlphaI;
    real arI2 = aI2*r;
    real expI = EXP(-arI2);
    real aI2_2 = aI2*aI2;
    real aI2_3 = aI2_2*aI2;
    real aI2_4 = aI2_2*aI2_2;
    real aI2_5 = aI2_3*aI2_2;
    real aI2_6 = aI2_3*aI2_3;
    real fexp, fexp1, fexp2, fexp3, fexp4, fexp5, pre;
    real one = 1;
    real two = 2;
    real four = 4;
    real eight = 8;
    real twelve = 12;
    real sixteen = 16;
    if (pauliAlphaI == pauliAlphaJ) {
        real r7 = r4*r3;
        real r8 = r4*r4;
        real aI2_7 = aI2_4*aI2_3;
        pre = 128;
        fexp = (r + aI2*r2 + aI2_2*r3*(one/3))*expI;
        fexp1 = (aI2_2*r3 + aI2_3*r4)*expI*(one/3);
        fexp2 = aI2_4*expI*r5*(one/9);
        fexp3 = aI2_5*expI*r6*(one/45);
        fexp4 = (aI2_5*r6 + aI2_6*r7)*expI*(one/315);
        fexp5 = (aI2_5*r6 + aI2_6*r7 + aI2_7*r8*(one/3))*expI*(one/945);
    }
    else {
        real aJ2 = 0.5f*pauliAlphaJ;
        real arJ2 = aJ2*r;
        real expJ = EXP(-arJ2);
        real aJ2_2 = aJ2*aJ2;
        real aJ2_3 = aJ2_2*aJ2;
        real aJ2_4 = aJ2_2*aJ2_2;
        real aJ2_5 = aJ2_3*aJ2_2;
        real aJ2_6 = aJ2_3*aJ2_3;
        real scale = 1/(aI2_2-aJ2_2);
        pre = 8192*aI2_3*aJ2_3*(scale*scale*scale*scale);
        real tmp = 4*aI2*aJ2*scale;
        fexp = (arI2-tmp)*expJ + (arJ2+tmp)*expI;
        fexp1 = (aI2*aJ2*r2 - 4*aI2*aJ2_2*r*scale - 4*aI2*aJ2*scale)*expJ +
             (aI2*aJ2*r2 + 4*aI2_2*aJ2*r*scale + 4*aI2*aJ2*scale)*expI;
        fexp2 = (aI2*aJ2*r2*(one/3) + aI2*aJ2_2*r3*(one/3) - (four/3)*aI2*aJ2_3*r2*scale - 4*aI2*aJ2_2*r*scale - 4*aI2*aJ2*scale)*expJ +
                (aI2*aJ2*r2*(one/3) + aI2_2*aJ2*r3*(one/3) + (four/3)*aI2_3*aJ2*r2*scale + 4*aI2_2*aJ2*r*scale + 4*aI2*aJ2*scale)*expI;
        fexp3 = (aI2*aJ2_3*r4*(one/15) + aI2*aJ2_2*r3*(one/5) + aI2*aJ2*r2*(one/5) - (four/15)*aI2*aJ2_4*r3*scale - (eight/5)*aI2*aJ2_3*r2*scale - 4*aI2*aJ2_2*r*scale - 4*scale*aI2*aJ2)*expJ +
                (aI2_3*aJ2*r4*(one/15) + aI2_2*aJ2*r3*(one/5) + aI2*aJ2*r2*(one/5) + (four/15)*aI2_4*aJ2*r3*scale + (eight/5)*aI2_3*aJ2*r2*scale + 4*aI2_2*aJ2*r*scale + 4*scale*aI2*aJ2)*expI;
        fexp4 = (aI2*aJ2_4*r5*(one/105) + (two/35)*aI2*aJ2_3*r4 + aI2*aJ2_2*r3*(one/7) + aI2*aJ2*r2*(one/7) - (four/105)*aI2*aJ2_5*r4*scale - (eight/21)*aI2*aJ2_4*r3*scale - (twelve/7)*aI2*aJ2_3*r2*scale - 4*aI2*aJ2_2*r*scale - 4*aI2*aJ2*scale)*expJ +
                (aI2_4*aJ2*r5*(one/105) + (two/35)*aI2_3*aJ2*r4 + aI2_2*aJ2*r3*(one/7) + aI2*aJ2*r2*(one/7) + (four/105)*aI2_5*aJ2*r4*scale + (eight/21)*aI2_4*aJ2*r3*scale + (twelve/7)*aI2_3*aJ2*r2*scale + 4*aI2_2*aJ2*r*scale + 4*aI2*aJ2*scale)*expI;
        fexp5 = (aI2*aJ2_5*r6*(one/945) + (two/189)*aI2*aJ2_4*r5 + aI2*aJ2_3*r4*(one/21) + aI2*aJ2_2*r3*(one/9) + aI2*aJ2*r2*(one/9) - (four/945)*aI2*aJ2_6*r5*scale - (four/63)*aI2*aJ2_5*r4*scale - (four/9)*aI2*aJ2_4*r3*scale - (sixteen/9)*aI2*aJ2_3*r2*scale - 4*aI2*aJ2_2*r*scale - 4*aI2*aJ2*scale)*expJ +
                (aI2_5*aJ2*r6*(one/945) + (two/189)*aI2_4*aJ2*r5 + aI2_3*aJ2*r4*(one/21) + aI2_2*aJ2*r3*(one/9) + aI2*aJ2*r2*(one/9) + (four/945)*aI2_6*aJ2*r5*scale + (four/63)*aI2_5*aJ2*r4*scale + (four/9)*aI2_4*aJ2*r3*scale + (sixteen/9)*aI2_3*aJ2*r2*scale + 4*aI2_2*aJ2*r*scale + 4*aI2*aJ2*scale)*expI;
    }
    fexp = fexp/r;
    fexp1 = fexp1/r3;
    fexp2 = 3*fexp2/r5;
    fexp3 = 15*fexp3/(r5*r2);
    fexp4 = 105*fexp4/(r5*r4);
    fexp5 = 945*fexp5/(r5*r6);
    fdamp1 = 0.5f*pre*fexp*fexp;
    fdamp3 = pre*fexp*fexp1;
    fdamp5 = pre*(fexp*fexp2 + fexp1*fexp1);
    fdamp7 = pre*(fexp*fexp3 + 3*fexp1*fexp2);
    fdamp9 = pre*(fexp*fexp4 + 4*fexp1*fexp3 + 3*fexp2*fexp2);
    fdamp11 = pre*(fexp*fexp5 + 5*fexp1*fexp4 + 10*fexp2*fexp3);
}

extern "C" __global__ void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const tileflags* __restrict__ exclusions,
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockSize, const unsigned int* __restrict__ interactingAtoms, unsigned int maxSinglePairs,
        const int2* __restrict__ singlePairs
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES
    // used shared memory if the device cannot shuffle
#ifndef ENABLE_SHUFFLE
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];
#endif

    // First loop: process tiles that contain exclusions.

    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real3 force = make_real3(0);
        real3 torque = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
        tileflags excl = exclusions[pos*TILE_SIZE+tgx];
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.
#ifdef ENABLE_SHUFFLE
            real4 shflPosq = posq1;
#else
            localData[threadIdx.x].x = posq1.x;
            localData[threadIdx.x].y = posq1.y;
            localData[threadIdx.x].z = posq1.z;
            localData[threadIdx.x].q = posq1.w;
            LOAD_LOCAL_PARAMETERS_FROM_1
#endif

            // we do not need to fetch parameters from global since this is a symmetric tile
            // instead we can broadcast the values using shuffle
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2;
#ifdef ENABLE_SHUFFLE
                BROADCAST_WARP_DATA
#else   
                posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real rInv = RSQRT(r2);
                real r = r2*rInv;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+j;
                real3 tempForce = make_real3(0);
                real3 tempTorque1 = make_real3(0);
                real3 tempTorque2 = make_real3(0);
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                real tempEnergy = 0.0f;
                const real interactionScale = 0.5f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
                force -= tempForce;
                torque += tempTorque1;
                excl >>= 1;
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            real4 shflPosq = posq[j];
#ifdef ENABLE_SHUFFLE
            real3 shflForce = make_real3(0);
            real3 shflTorque = make_real3(0);
#else
            localData[threadIdx.x].x = shflPosq.x;
            localData[threadIdx.x].y = shflPosq.y;
            localData[threadIdx.x].z = shflPosq.z;
            localData[threadIdx.x].q = shflPosq.w;
            localData[threadIdx.x].fx = 0.0f;
            localData[threadIdx.x].fy = 0.0f;
            localData[threadIdx.x].fz = 0.0f;
            localData[threadIdx.x].tx = 0.0f;
            localData[threadIdx.x].ty = 0.0f;
            localData[threadIdx.x].tz = 0.0f;
#endif
            DECLARE_LOCAL_PARAMETERS
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                real4 posq2 = shflPosq;
#else
                real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real rInv = RSQRT(r2);
                real r = r2*rInv;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+tj;
                real3 tempForce = make_real3(0);
                real3 tempTorque1 = make_real3(0);
                real3 tempTorque2 = make_real3(0);
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                real tempEnergy = 0.0f;
                const real interactionScale = 1.0f;
                COMPUTE_INTERACTION
                energy += tempEnergy;
                force -= tempForce;
                torque + tempTorque1;
#ifdef ENABLE_SHUFFLE
                shflForce += tempForce;
                shflTorque += tempTorque2;
                SHUFFLE_WARP_DATA
                shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                localData[tbx+tj].fx += tempForce.x;
                localData[tbx+tj].fy += tempForce.y;
                localData[tbx+tj].fz += tempForce.z;
                localData[tbx+tj].tx += tempTorque2.x;
                localData[tbx+tj].ty += tempTorque2.y;
                localData[tbx+tj].tz += tempTorque2.z;
#endif
                excl >>= 1;
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
            // write results for off diagonal tiles
#ifdef ENABLE_SHUFFLE
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (shflForce.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (shflTorque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.z*0x100000000)));
#else
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tx*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].ty*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tz*0x100000000)));
#endif
        }
        // Write results for on and off diagonal tiles

        const unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
        atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (torque.x*0x100000000)));
        atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.y*0x100000000)));
        atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.z*0x100000000)));
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        const bool hasExclusions = false;
        real3 force = make_real3(0);
        real3 torque = make_real3(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
#else
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
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;
            // Load atom data for this tile.
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            //const unsigned int localAtomIndex = threadIdx.x;
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
#ifdef ENABLE_SHUFFLE
            DECLARE_LOCAL_PARAMETERS
            real4 shflPosq;
            real3 shflForce = make_real3(0);
            real3 shflTorque = make_real3(0);
#endif
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from from global memory
#ifdef ENABLE_SHUFFLE
                shflPosq = posq[j];
#else
                localData[threadIdx.x].x = posq[j].x;
                localData[threadIdx.x].y = posq[j].y;
                localData[threadIdx.x].z = posq[j].z;
                localData[threadIdx.x].q = posq[j].w;
                localData[threadIdx.x].fx = 0.0f;
                localData[threadIdx.x].fy = 0.0f;
                localData[threadIdx.x].fz = 0.0f;
#endif                
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            else {
#ifdef ENABLE_SHUFFLE
                shflPosq = make_real4(0, 0, 0, 0);
#else
                localData[threadIdx.x].x = 0;
                localData[threadIdx.x].y = 0;
                localData[threadIdx.x].z = 0;
#endif
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.
                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
#ifdef ENABLE_SHUFFLE
                APPLY_PERIODIC_TO_POS_WITH_CENTER(shflPosq, blockCenterX)
#else
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[threadIdx.x], blockCenterX)
#endif
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                    real4 posq2 = shflPosq; 
#else
                    real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real rInv = RSQRT(r2);
                    real r = r2*rInv;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
                    real3 tempForce = make_real3(0);
                    real3 tempTorque1 = make_real3(0);
                    real3 tempTorque2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
                    force -= tempForce;
                    torque += tempTorque1;
#ifdef ENABLE_SHUFFLE
                    shflForce += tempForce;
                    shflTorque += tempTorque2;
                    SHUFFLE_WARP_DATA
                    shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                    shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                    shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                    localData[tbx+tj].fx += tempForce.x;
                    localData[tbx+tj].fy += tempForce.y;
                    localData[tbx+tj].fz += tempForce.z;
                    localData[tbx+tj].tx += tempTorque2.x;
                    localData[tbx+tj].ty += tempTorque2.y;
                    localData[tbx+tj].tz += tempTorque2.z;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
#ifdef ENABLE_SHUFFLE
                    real4 posq2 = shflPosq;
#else
                    real4 posq2 = make_real4(localData[atom2].x, localData[atom2].y, localData[atom2].z, localData[atom2].q);
#endif
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real rInv = RSQRT(r2);
                    real r = r2*rInv;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
                    real3 tempForce = make_real3(0);
                    real3 tempTorque1 = make_real3(0);
                    real3 tempTorque2 = make_real3(0);
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
                    force -= tempForce;
                    torque += tempTorque1;
#ifdef ENABLE_SHUFFLE
                    shflForce += tempForce;
                    shflTorque += tempTorque2;
                    SHUFFLE_WARP_DATA
                    shflTorque.x = real_shfl(shflTorque.x, tgx+1);
                    shflTorque.y = real_shfl(shflTorque.y, tgx+1);
                    shflTorque.z = real_shfl(shflTorque.z, tgx+1);
#else
                    localData[tbx+tj].fx += tempForce.x;
                    localData[tbx+tj].fy += tempForce.y;
                    localData[tbx+tj].fz += tempForce.z;
                    localData[tbx+tj].tx += tempTorque.x;
                    localData[tbx+tj].ty += tempTorque.y;
                    localData[tbx+tj].tz += tempTorque.z;
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results.

            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (torque.z*0x100000000)));
#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[threadIdx.x];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
#ifdef ENABLE_SHUFFLE
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (shflForce.x*0x100000000)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.y*0x100000000)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflForce.z*0x100000000)));
                atomicAdd(&torqueBuffers[atom2], static_cast<unsigned long long>((long long) (shflTorque.x*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.y*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (shflTorque.z*0x100000000)));
#else
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fx*0x100000000)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fy*0x100000000)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fz*0x100000000)));
                atomicAdd(&torqueBuffers[atom2], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tx*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].ty*0x100000000)));
                atomicAdd(&torqueBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].tz*0x100000000)));
#endif
            }
        }
        pos++;
    }
#ifdef INCLUDE_ENERGY
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
#endif
    SAVE_DERIVATIVES
}