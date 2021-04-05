#define TILE_SIZE 32

/**
 * Reduce the Born sums to compute the Born radii.
 */
KERNEL void reduceBornSum(GLOBAL const mm_long* RESTRICT bornSum, GLOBAL const float2* RESTRICT params, GLOBAL real* RESTRICT bornRadii) {
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
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
KERNEL void computeSurfaceAreaForce(GLOBAL mm_long* RESTRICT bornForce, GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const float2* RESTRICT params, GLOBAL const real* RESTRICT bornRadii) {
    mixed energy = 0;
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        real bornRadius = bornRadii[index];
        float radius = params[index].x;
        real r = radius + DIELECTRIC_OFFSET + PROBE_RADIUS;
        real ratio6 = (radius+DIELECTRIC_OFFSET)/bornRadius;
        ratio6 = ratio6*ratio6*ratio6;
        ratio6 = ratio6*ratio6;
        real saTerm = SURFACE_AREA_FACTOR * r * r * ratio6;
        bornForce[index] += (mm_long) (saTerm*0x100000000/bornRadius);
        energy += saTerm;
    }
    energyBuffer[GLOBAL_ID] -= energy/6;
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

DEVICE void computeBornSumOneInteraction(AtomData1& atom1, AtomData1& atom2) {
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
KERNEL void computeBornSum(GLOBAL mm_ulong* RESTRICT bornSum, GLOBAL const real4* RESTRICT posq,
        GLOBAL const float2* RESTRICT params, unsigned int numTiles) {
    unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    unsigned int pos = (unsigned int) (warp*(mm_long)numTiles/totalWarps);
    unsigned int end = (unsigned int) ((warp+1)*(mm_long)numTiles/totalWarps);
    unsigned int lasty = 0xFFFFFFFF;
    LOCAL AtomData1 localData[BORN_SUM_THREAD_BLOCK_SIZE];
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
        const unsigned int tbx = LOCAL_ID - tgx;
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

                localData[LOCAL_ID].pos = data.pos;
                localData[LOCAL_ID].radius = params1.x;
                localData[LOCAL_ID].scaledRadius = params1.y;
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
                    localData[LOCAL_ID].pos = trimTo3(tempPosq);
                    float2 tempParams = params[j];
                    localData[LOCAL_ID].radius = tempParams.x;
                    localData[LOCAL_ID].scaledRadius = tempParams.y;
                }
                localData[LOCAL_ID].bornSum = 0;
                
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
            ATOMIC_ADD(&bornSum[offset], (mm_ulong) ((mm_long) (data.bornSum*0x100000000)));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&bornSum[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].bornSum*0x100000000)));
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

DEVICE void computeOneInteractionF1(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force);
DEVICE void computeOneInteractionF2(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force);
DEVICE void computeOneInteractionT1(AtomData2& atom1, volatile AtomData2& atom2, real3& torque);
DEVICE void computeOneInteractionT2(AtomData2& atom1, volatile AtomData2& atom2, real3& torque);
DEVICE void computeOneInteractionB1B2(AtomData2& atom1, volatile AtomData2& atom2);

inline DEVICE void loadAtomData2(AtomData2& data, int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real3* RESTRICT inducedDipole, GLOBAL const real3* RESTRICT inducedDipolePolar, GLOBAL const real* RESTRICT bornRadius) {
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
    data.inducedDipole = inducedDipole[atom];
    data.inducedDipolePolar = inducedDipolePolar[atom];
    data.bornRadius = bornRadius[atom];
}

inline DEVICE void zeroAtomData(AtomData2& data) {
    data.force = make_real3(0);
    data.bornForce = 0;
}

/**
 * Compute electrostatic interactions.
 */
KERNEL void computeGKForces(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT torqueBuffers, GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL const real4* RESTRICT posq, unsigned int startTileIndex, unsigned int numTileIndices, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real3* RESTRICT inducedDipole, GLOBAL const real3* RESTRICT inducedDipolePolar,
        GLOBAL const real* RESTRICT bornRadii, GLOBAL mm_ulong* RESTRICT bornForce) {
    unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = (unsigned int) (startTileIndex+warp*(mm_long)numTiles/totalWarps);
    unsigned int end = (unsigned int) (startTileIndex+(warp+1)*(mm_long)numTiles/totalWarps);
    mixed energy = 0;
    LOCAL AtomData2 localData[GK_FORCE_THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
        const unsigned int tbx = LOCAL_ID - tgx;
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

                localData[LOCAL_ID].pos = data.pos;
                localData[LOCAL_ID].q = data.q;
                localData[LOCAL_ID].dipole = data.dipole;
                localData[LOCAL_ID].quadrupoleXX = data.quadrupoleXX;
                localData[LOCAL_ID].quadrupoleXY = data.quadrupoleXY;
                localData[LOCAL_ID].quadrupoleXZ = data.quadrupoleXZ;
                localData[LOCAL_ID].quadrupoleYY = data.quadrupoleYY;
                localData[LOCAL_ID].quadrupoleYZ = data.quadrupoleYZ;
                localData[LOCAL_ID].quadrupoleZZ = data.quadrupoleZZ;
                localData[LOCAL_ID].inducedDipole = data.inducedDipole;
                localData[LOCAL_ID].inducedDipolePolar = data.inducedDipolePolar;
                localData[LOCAL_ID].bornRadius = data.bornRadius;
                
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
                ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
                
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
                ATOMIC_ADD(&torqueBuffers[atom1], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
                ATOMIC_ADD(&torqueBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
                ATOMIC_ADD(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
                
                // Compute chain rule terms.
                
                zeroAtomData(data);
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                        computeOneInteractionB1B2(data, localData[tbx+j]);
                }
                ATOMIC_ADD(&bornForce[atom1], (mm_ulong) ((mm_long) (data.bornForce*0x100000000)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData2(localData[LOCAL_ID], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, bornRadii);
                zeroAtomData(localData[LOCAL_ID]);
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
                localData[LOCAL_ID].force *= 0.5f;
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));
                }

                // Compute torques.

                zeroAtomData(data);
                zeroAtomData(localData[LOCAL_ID]);
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
                    ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
                    ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
                    ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
                    ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
                    ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));
                }

                // Compute chain rule terms.

                zeroAtomData(data);
                zeroAtomData(localData[LOCAL_ID]);
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                        computeOneInteractionB1B2(data, localData[tbx+tj]);
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&bornForce[offset], (mm_ulong) ((mm_long) (data.bornForce*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&bornForce[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].bornForce*0x100000000)));
                }
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[GLOBAL_ID] += energy*0.5f;
}


/**
 * Data structure used by computeChainRuleForce().
 */
typedef struct {
    real3 pos, force;
    real radius, scaledRadius, bornRadius, bornForce;
} AtomData3;

inline DEVICE void loadAtomData3(AtomData3& data, int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const float2* RESTRICT params,
        GLOBAL const real* RESTRICT bornRadius, GLOBAL const mm_long* RESTRICT bornForce) {
    data.pos = trimTo3(posq[atom]);
    data.bornRadius = bornRadius[atom];
    float2 params1 = params[atom];
    data.radius = params1.x;
    data.scaledRadius = params1.y;
    data.bornForce = bornForce[atom]/(real) 0x100000000;
}

DEVICE void computeBornChainRuleInteraction(AtomData3& atom1, AtomData3& atom2, real3& force) {
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
KERNEL void computeChainRuleForce(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL const real4* RESTRICT posq, unsigned int startTileIndex, unsigned int numTileIndices,
        GLOBAL const float2* RESTRICT params, GLOBAL const real* RESTRICT bornRadii, GLOBAL const mm_long* RESTRICT bornForce) {
    unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = startTileIndex+warp*numTiles/totalWarps;
    unsigned int end = startTileIndex+(warp+1)*numTiles/totalWarps;
    LOCAL AtomData3 localData[CHAIN_RULE_THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
        const unsigned int tbx = LOCAL_ID - tgx;
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

                localData[LOCAL_ID].pos = data.pos;
                localData[LOCAL_ID].radius = data.radius;
                localData[LOCAL_ID].scaledRadius = data.scaledRadius;
                localData[LOCAL_ID].bornRadius = data.bornRadius;
                localData[LOCAL_ID].bornForce = data.bornForce;
                localData[LOCAL_ID].force = make_real3(0);
                
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
                ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) ((mm_long) ((data.force.x+localData[LOCAL_ID].force.x)*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) ((data.force.y+localData[LOCAL_ID].force.y)*0x100000000)));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) ((data.force.z+localData[LOCAL_ID].force.z)*0x100000000)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                loadAtomData3(localData[LOCAL_ID], j, posq, params, bornRadii, bornForce);
                localData[LOCAL_ID].force = make_real3(0);
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
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));
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

DEVICE void computeOneEDiffInteractionF1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real& outputEnergy, real3& outputForce);
DEVICE void computeOneEDiffInteractionT1(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce);
DEVICE void computeOneEDiffInteractionT3(AtomData4& atom1, volatile AtomData4& atom2, float dScale, float pScale, real3& outputForce);

inline DEVICE void loadAtomData4(AtomData4& data, int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real3* RESTRICT inducedDipole, GLOBAL const real3* RESTRICT inducedDipolePolar,
        GLOBAL const real3* RESTRICT inducedDipoleS, GLOBAL const real3* RESTRICT inducedDipolePolarS, GLOBAL const float2* RESTRICT dampingAndThole) {
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
    data.inducedDipole = inducedDipole[atom];
    data.inducedDipolePolar = inducedDipolePolar[atom];
    data.inducedDipoleS = inducedDipoleS[atom];
    data.inducedDipolePolarS = inducedDipolePolarS[atom];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

DEVICE real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

DEVICE float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

/**
 * Compute electrostatic interactions.
 */
KERNEL void computeEDiffForce(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT torqueBuffers, GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL const real4* RESTRICT posq, GLOBAL const uint2* RESTRICT covalentFlags, GLOBAL const unsigned int* RESTRICT polarizationGroupFlags,
        GLOBAL const int2* RESTRICT exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
        GLOBAL const real* RESTRICT labFrameDipole, GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real3* RESTRICT inducedDipole,
        GLOBAL const real3* RESTRICT inducedDipolePolar, GLOBAL const real3* RESTRICT inducedDipoleS, GLOBAL const real3* RESTRICT inducedDipolePolarS,
        GLOBAL const float2* RESTRICT dampingAndThole) {
    const unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    const unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    mixed energy = 0;
    LOCAL AtomData4 localData[EDIFF_THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
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

            localData[LOCAL_ID].pos = data.pos;
            localData[LOCAL_ID].q = data.q;
            localData[LOCAL_ID].dipole = data.dipole;
            localData[LOCAL_ID].quadrupoleXX = data.quadrupoleXX;
            localData[LOCAL_ID].quadrupoleXY = data.quadrupoleXY;
            localData[LOCAL_ID].quadrupoleXZ = data.quadrupoleXZ;
            localData[LOCAL_ID].quadrupoleYY = data.quadrupoleYY;
            localData[LOCAL_ID].quadrupoleYZ = data.quadrupoleYZ;
            localData[LOCAL_ID].quadrupoleZZ = data.quadrupoleZZ;
            localData[LOCAL_ID].inducedDipole = data.inducedDipole;
            localData[LOCAL_ID].inducedDipolePolar = data.inducedDipolePolar;
            localData[LOCAL_ID].inducedDipoleS = data.inducedDipoleS;
            localData[LOCAL_ID].inducedDipolePolarS = data.inducedDipolePolarS;
            localData[LOCAL_ID].thole = data.thole;
            localData[LOCAL_ID].damp = data.damp;

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
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));

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
            ATOMIC_ADD(&torqueBuffers[atom1], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData4(localData[LOCAL_ID], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[LOCAL_ID].force = make_real3(0);

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
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));

            // Compute torques.

            data.force = make_real3(0);
            localData[LOCAL_ID].force = make_real3(0);
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
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions (by enumerating all of them, since there's no cutoff).

    const unsigned int numTiles = numTileIndices;
    int pos = startTileIndex+warp*numTiles/totalWarps;
    int end = startTileIndex+(warp+1)*numTiles/totalWarps;
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL volatile int skipTiles[EDIFF_THREAD_BLOCK_SIZE];
    skipTiles[LOCAL_ID] = -1;

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
                int2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[LOCAL_ID] = end;
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
            loadAtomData4(localData[LOCAL_ID], atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData4(localData[LOCAL_ID], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[LOCAL_ID].force = make_real3(0);

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
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));

            // Compute torques.

            data.force = make_real3(0);
            localData[LOCAL_ID].force = make_real3(0);
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
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (data.force.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.force.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.x*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.y*0x100000000)));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (localData[LOCAL_ID].force.z*0x100000000)));
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy*ENERGY_SCALE_FACTOR;
}
