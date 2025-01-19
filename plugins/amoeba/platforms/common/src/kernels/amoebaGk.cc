#define TILE_SIZE 32

/**
 * Reduce the Born sums to compute the Born radii.
 */
KERNEL void reduceBornSum(GLOBAL const mm_long* RESTRICT bornSum, GLOBAL const float4* RESTRICT params, GLOBAL real* RESTRICT bornRadii) {

    real RECIP_MAX_RADIUS3 = POW(BIG_RADIUS, (real) -3);
    real PI4_3 = 4 * M_PI / 3;
    real INVERSE_PI4_3 = 1 / PI4_3;
    real ONE_THIRD = 1 / (real) 3;

    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        // Get summed Born data
        real sum = RECIP((real) 0x100000000)*bornSum[index];
        sum = sum * PI4_3;

        real radius = params[index].x;
        real radius3 = radius * radius * radius;
        real ir3 = RECIP(radius3);

        // Scale up the solute integral in account for interstitial spaces.
        if (TANH_RESCALING) {
            // Set up tanh function components
            real rhoi3Psi = radius3 * sum;
            real rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
            real rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
            // If the output of the tanh function is 1.0, then the Born radius will be MaxBornRadius
            real tanh_constant = PI4_3 * (ir3 - RECIP_MAX_RADIUS3);
            sum = tanh_constant * tanh(BETA0 * rhoi3Psi - BETA1 * rhoi6Psi2 + BETA2 * rhoi9Psi3);
        }

        // Now calculate Born radius.
        sum = PI4_3 * ir3 - sum;

        // If the sum is less than zero, set the Born radius to 30.0 Angstroms.
        sum = (sum <= 0 ? BIG_RADIUS : POW(INVERSE_PI4_3 * sum, -ONE_THIRD));
        bornRadii[index] = sum;

        // Born radius should be at least as large as its base radius.
        if (bornRadii[index] < radius) {
            bornRadii[index] = radius;
        }

        // Check if we have exceeded the maximum Born radius.
        if (bornRadii[index] > BIG_RADIUS) {
            bornRadii[index] = BIG_RADIUS;
        }
    }
}

#ifdef SURFACE_AREA_FACTOR
/**
 * Apply the surface area term to the force and energy.
 */
KERNEL void computeSurfaceAreaForce(GLOBAL mm_long* RESTRICT bornForce, GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const float4* RESTRICT params, GLOBAL const real* RESTRICT bornRadii) {
    mixed energy = 0;
    for (unsigned int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        real bornRadius = bornRadii[index];
        float radius = params[index].x;
        real r = radius + DIELECTRIC_OFFSET + PROBE_RADIUS;
        real ratio6 = (radius+DIELECTRIC_OFFSET)/bornRadius;
        ratio6 = ratio6*ratio6*ratio6;
        ratio6 = ratio6*ratio6;
        real saTerm = SURFACE_AREA_FACTOR * r * r * ratio6;
        bornForce[index] += realToFixedPoint(saTerm/bornRadius);
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
    float radius, scaleFactor, descreenRadius, neckFactor, padding;
} AtomData1;

DEVICE real interpolate2D(real x1, real x2, real y1, real y2, real x, real y,
                            real fx1y1, real fx2y1, real fx1y2, real fx2y2) {
    real fxy1 = (x2 - x) / (x2 - x1) * fx1y1 + (x - x1) / (x2 - x1) * fx2y1;
    real fxy2 = (x2 - x) / (x2 - x1) * fx1y2 + (x - x1) / (x2 - x1) * fx2y2;

    return (y2 - y) / (y2 - y1) * fxy1 + (y - y1) / (y2 - y1) * fxy2;
}

DEVICE void getBounds(real rho, int* below, int* above) {
    // Note this function is in Angstroms.
    real MINIMUM_RADIUS = 0.8f;
    real SPACING = 0.05f;
    real calculateIndex = (rho - MINIMUM_RADIUS) / SPACING;
    *below = (int) floor(calculateIndex);
    *above = *below + 1;
    int NUM_POINTS = 45;

    if (*above >= NUM_POINTS) {
        // Extrapolate up from the top table values.
        *below = NUM_POINTS - 1;
        *above = NUM_POINTS - 2;
    }
    else if (*below < 0) {
        // If below is less than 0, extrapolate down from the bottom table values.
        *below = 0;
        *above = 1;
    }
}

DEVICE void getNeckConstants(real rhoDescreened, real rhoDescreening, real* aij, real *bij,
    GLOBAL const float* RESTRICT neckRadii,
    GLOBAL const float* RESTRICT neckA,
    GLOBAL const float* RESTRICT neckB) {

    // Convert the input radii from nm to A.
    rhoDescreened *= 10;
    rhoDescreening *= 10;

    // Determine low and high values for integration
    int lowI, highI, lowJ, highJ;
    getBounds(rhoDescreened, &lowI, &highI);
    getBounds(rhoDescreening, &lowJ, &highJ);
    // Interpolate/Extrapolate Aij and Bij constant values
    int NUM_POINTS = 45;
    int lowOff = lowI * NUM_POINTS;
    int highOff = highI * NUM_POINTS;
    *aij = interpolate2D(neckRadii[lowI], neckRadii[highI], neckRadii[lowJ], neckRadii[highJ],
                               rhoDescreened, rhoDescreening,
                               neckA[lowOff + lowJ], neckA[highOff + lowJ],
                               neckA[lowOff + highJ], neckA[highOff + highJ]);

    *bij = interpolate2D(neckRadii[lowI], neckRadii[highI],neckRadii[lowJ], neckRadii[highJ],
                               rhoDescreened, rhoDescreening,
                               neckB[lowOff + lowJ], neckB[highOff + lowJ],
                               neckB[lowOff + highJ], neckB[highOff + highJ]);

    // Never let Aij be negative.
    if (*aij < 0.0f) {
        *aij = 0.0f;
    }

    // Convert aij from A^(-11) to nm^(-11);
    *aij /= 10;
    // Convert bij from A to nm.
    *bij /= 10;
}

DEVICE real neckDescreen(real r, real radius, real radiusK, real sneck,
        GLOBAL const float* RESTRICT neckRadii,
        GLOBAL const float* RESTRICT neckA,
        GLOBAL const float* RESTRICT neckB) {

    real radiusWater = 0.14f;

    // If atoms are too widely separated there is no neck formed.
    if (r > radius + radiusK + 2 * radiusWater) {
        return 0.0f;
    }

    // Get Aij and Bij based on parameterization by Corrigan et al.
    real aij, bij;
    getNeckConstants(radius, radiusK, &aij, &bij, neckRadii, neckA, neckB);
    real rMinusBij = r - bij;
    real radiiMinusr = radius + radiusK + 2 * radiusWater - r;
    real power1 = rMinusBij * rMinusBij * rMinusBij * rMinusBij;
    real power2 = radiiMinusr * radiiMinusr * radiiMinusr * radiiMinusr;

    // Use Aij and Bij to get neck integral using Equations 13 and 14 from Aguilar/Onufriev 2010 paper
    // Sneck may be based on the number of heavy atoms bound to the atom being descreened.
    real PI4_3 = 4 * M_PI / (real) 3;
    real neckIntegral = PI4_3 * sneck * aij * power1 * power2;

    return neckIntegral;
}


DEVICE real computeBornSumOneInteraction(AtomData1 atom1, AtomData1 atom2,
        GLOBAL const float* RESTRICT neckRadii,
        GLOBAL const float* RESTRICT neckA,
        GLOBAL const float* RESTRICT neckB) {

    if (atom1.radius <= 0)
        return 0; // Ignore this interaction

    float sk = atom2.scaleFactor * atom2.descreenRadius;
    if (sk <= 0.0f)
        return 0; // No descreening.

    real3 delta = atom2.pos - atom1.pos;
    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    real baseRadius = max(atom1.radius, atom1.descreenRadius) + DESCREEN_OFFSET;

    if (baseRadius > r + sk)
        return 0; // No descreening due to atom1 engulfing atom2.

    real sk2 = sk*sk;
    if (baseRadius + r < sk) {
        real lik = baseRadius;
        real uik = sk - r;
        atom1.bornSum -= RECIP(uik*uik*uik) - RECIP(lik*lik*lik);
    }
    real uik = r+sk;
    real lik;
    if (baseRadius+r < sk)
        lik = sk-r;
    else if (r < baseRadius+sk)
        lik = baseRadius;
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
    term = term / (real) 16;

    real neck1 = atom1.neckFactor;
    real neck2 = atom2.neckFactor;
    real mixedNeckScale = 0.5f * (neck1 + neck2);
    if (mixedNeckScale > 0 && atom2.scaleFactor > 0) {
        real ret = neckDescreen(r, baseRadius, atom2.descreenRadius, mixedNeckScale, neckRadii, neckA, neckB);
        real pi43 = 4 * M_PI / (real) 3;
        term += ret / pi43;
    }

    return term;
}

/**
 * Compute the Born sum.
 */
KERNEL void computeBornSum(GLOBAL mm_ulong* RESTRICT bornSum, GLOBAL const real4* RESTRICT posq,
        GLOBAL const float4* RESTRICT params, unsigned int numTiles,
        GLOBAL const float* RESTRICT neckRadii,
        GLOBAL const float* RESTRICT neckA,
        GLOBAL const float* RESTRICT neckB) {
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
            float4 params1 = params[atom1];
            data.radius = params1.x;
            data.scaleFactor = params1.y;
            data.descreenRadius = params1.z;
            data.neckFactor = params1.w;
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.
                localData[LOCAL_ID].pos = data.pos;
                localData[LOCAL_ID].radius = params1.x;
                localData[LOCAL_ID].scaleFactor = params1.y;
                localData[LOCAL_ID].descreenRadius = params1.z;
                localData[LOCAL_ID].neckFactor = params1.w;

                SYNC_WARPS;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
                        real bornSum = computeBornSumOneInteraction(data, localData[tbx+j], neckRadii, neckA, neckB);
                        data.bornSum += bornSum;
                    }
                }
                SYNC_WARPS;
            }
            else {
                // This is an off-diagonal tile.

                if (lasty != y) {
                    unsigned int j = y*TILE_SIZE + tgx;
                    real4 tempPosq = posq[j];
                    localData[LOCAL_ID].pos = trimTo3(tempPosq);
                    float4 paramsJ = params[j];
                    localData[LOCAL_ID].radius = paramsJ.x;
                    localData[LOCAL_ID].scaleFactor = paramsJ.y;
                    localData[LOCAL_ID].descreenRadius = paramsJ.z;
                    localData[LOCAL_ID].neckFactor = paramsJ.w;
                }
                localData[LOCAL_ID].bornSum = 0;
                SYNC_WARPS;
                
                // Compute the full set of interactions in this tile.

                unsigned int tj = tgx;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real bornSum = computeBornSumOneInteraction(data, localData[tbx+tj], neckRadii, neckA, neckB);
                        data.bornSum += bornSum;
                        bornSum = computeBornSumOneInteraction(localData[tbx+tj], data, neckRadii, neckA, neckB);
                        localData[tbx+tj].bornSum += bornSum;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
        }
        
        // Write results.
        
        if (pos < end) {
            const unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&bornSum[offset], (mm_ulong) realToFixedPoint(data.bornSum));
        }
        if (pos < end && x != y) {
            const unsigned int offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&bornSum[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].bornSum));
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

#if defined(USE_HIP)
#define ATOM2_ARG_SPEC
#else
#define ATOM2_ARG_SPEC volatile
#endif

DEVICE void computeOneInteractionF1(AtomData2 atom1, ATOM2_ARG_SPEC AtomData2 atom2, real* outputEnergy, real3* force);
DEVICE void computeOneInteractionF2(AtomData2 atom1, ATOM2_ARG_SPEC AtomData2 atom2, real* outputEnergy, real3* force);
DEVICE void computeOneInteractionT1(AtomData2 atom1, ATOM2_ARG_SPEC AtomData2 atom2, real3* torque);
DEVICE void computeOneInteractionT2(AtomData2 atom1, ATOM2_ARG_SPEC AtomData2 atom2, real3* torque);
DEVICE void computeOneInteractionB1B2(AtomData2 atom1, ATOM2_ARG_SPEC AtomData2 atom2, real* bornForce1, real* bornForce2);

inline DEVICE AtomData2 loadAtomData2(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT inducedDipolePolar, GLOBAL const real* RESTRICT bornRadius) {
    AtomData2 data;
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
    data.inducedDipole = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
    data.inducedDipolePolar = make_real3(inducedDipolePolar[3*atom], inducedDipolePolar[3*atom+1], inducedDipolePolar[3*atom+2]);
    data.bornRadius = bornRadius[atom];
    return data;
}

/**
 * Compute electrostatic interactions.
 */
KERNEL void computeGKForces(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mm_ulong* RESTRICT torqueBuffers, GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL const real4* RESTRICT posq, unsigned int startTileIndex, unsigned int numTileIndices, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT inducedDipolePolar,
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
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            AtomData2 data = loadAtomData2(atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, bornRadii);
            data.force = make_real3(0);
            data.bornForce = 0;
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
                SYNC_WARPS;
                
                // Compute forces.
                
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        real tempEnergy;
                        computeOneInteractionF1(data, localData[tbx+j], &tempEnergy, &tempForce);
                        computeOneInteractionF2(data, localData[tbx+j], &tempEnergy, &tempForce);
                        data.force += tempForce;
                        energy += 0.5f*tempEnergy;
                    }
                }
                SYNC_WARPS;
                data.force *= 0.5f;
                ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) realToFixedPoint(data.force.x));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));

                // Compute torques.
                
                data.force = make_real3(0);
                data.bornForce = 0;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempTorque;
                        computeOneInteractionT1(data, localData[tbx+j], &tempTorque);
                        computeOneInteractionT2(data, localData[tbx+j], &tempTorque);
                        data.force += tempTorque;
                    }
                }
                SYNC_WARPS;
                ATOMIC_ADD(&torqueBuffers[atom1], (mm_ulong) realToFixedPoint(data.force.x));
                ATOMIC_ADD(&torqueBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
                ATOMIC_ADD(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));

                // Compute chain rule terms.
                
                data.force = make_real3(0);
                data.bornForce = 0;
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real bornForce1 = 0, bornForce2 = 0;
                        computeOneInteractionB1B2(data, localData[tbx+j], &bornForce1, &bornForce2);
                        data.bornForce += bornForce1;
                        localData[tbx+j].bornForce += bornForce2;
                        SYNC_WARPS;
                    }
                }
                ATOMIC_ADD(&bornForce[atom1], (mm_ulong) realToFixedPoint(data.bornForce));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                localData[LOCAL_ID] = loadAtomData2(j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, bornRadii);
                localData[LOCAL_ID].force = make_real3(0);
                localData[LOCAL_ID].bornForce = 0;
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        real tempEnergy;
                        computeOneInteractionF1(data, localData[tbx+tj], &tempEnergy, &tempForce);
                        computeOneInteractionF2(data, localData[tbx+tj], &tempEnergy, &tempForce);
                        data.force += tempForce;
                        localData[tbx+tj].force -= tempForce;
                        energy += tempEnergy;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
                data.force *= 0.5f;
                localData[LOCAL_ID].force *= 0.5f;
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
                }

                // Compute torques.

                data.force = make_real3(0);
                data.bornForce = 0;
                localData[LOCAL_ID].force = make_real3(0);
                localData[LOCAL_ID].bornForce = 0;
                SYNC_WARPS;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempTorque;
                        computeOneInteractionT1(data, localData[tbx+tj], &tempTorque);
                        computeOneInteractionT2(data, localData[tbx+tj], &tempTorque);
                        data.force += tempTorque;
                        computeOneInteractionT1(localData[tbx+tj], data, &tempTorque);
                        computeOneInteractionT2(localData[tbx+tj], data, &tempTorque);
                        localData[tbx+tj].force += tempTorque;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
                    ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
                    ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
                    ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
                    ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
                }

                // Compute chain rule terms.

                data.force = make_real3(0);
                data.bornForce = 0;
                localData[LOCAL_ID].force = make_real3(0);
                localData[LOCAL_ID].bornForce = 0;
                SYNC_WARPS;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real bornForce1 = 0, bornForce2 = 0;
                        computeOneInteractionB1B2(data, localData[tbx+tj], &bornForce1, &bornForce2);
                        data.bornForce += bornForce1;
                        localData[tbx+tj].bornForce += bornForce2;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&bornForce[offset], (mm_ulong) realToFixedPoint(data.bornForce));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&bornForce[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].bornForce));
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
    real radius, scaleFactor, descreenRadius, neckFactor, bornSum, bornRadius, bornForce;
} AtomData3;

inline DEVICE AtomData3 loadAtomData3(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const float4* RESTRICT params,
                                      GLOBAL const mm_long* RESTRICT bornSum,
                                      GLOBAL const real* RESTRICT bornRadius, GLOBAL const mm_long* RESTRICT bornForce) {
    AtomData3 data;
    data.pos = trimTo3(posq[atom]);
    data.bornRadius = bornRadius[atom];
    data.bornSum = bornSum[atom]/(real) 0x100000000;
    float4 params1 = params[atom];
    data.radius = params1.x;
    data.scaleFactor = params1.y;
    data.descreenRadius = params1.z;
    data.neckFactor = params1.w;
    data.bornForce = bornForce[atom]/(real) 0x100000000;
    return data;
}

/**
 * Use pairwise descreening to compute derivative of the integral of 1/r^6 with respect to r.
 *
 * @param r separation distance.
 * @param r2 separation distance squared.
 * @param radius base radius of descreened atom.
 * @param scaledRadius scaled radius descreening atom.
 * @return the derivative.
 */
DEVICE real pairIntegralDerivative(real r, real r2, real radius, real scaledRadius) {

        real de = 0.0f;
        // Descreen only if the scaledRadius is greater than zero.
        // and atom I does not engulf atom K.
        if (scaledRadius > 0.0f && (radius < r + scaledRadius)) {
            // Atom i is engulfed by atom k.
            if (radius + r < scaledRadius) {
                real uik = scaledRadius - r;
                real uik2 = uik * uik;
                real uik4 = uik2 * uik2;
                de = -4.0f * M_PI / uik4;
            }

            // Lower integration bound depends on atoms sizes and separation.
            real sk2 = scaledRadius * scaledRadius;
            if (radius + r < scaledRadius) {
                // Atom i is engulfed by atom k.
                real lik = scaledRadius - r;
                real lik2 = lik * lik;
                real lik4 = lik2 * lik2;
                de = de + 0.25f * M_PI * (sk2 - 4.0f * scaledRadius * r + 17.0f * r2) / (r2 * lik4);
            }
            else if (r < radius + scaledRadius) {
                // Atoms are overlapped, begin integration from ri.
                real lik = radius;
                real lik2 = lik * lik;
                real lik4 = lik2 * lik2;
                de = de + 0.25f * M_PI * (2.0f * radius * radius - sk2 - r2) / (r2 * lik4);
            }
            else {
                // No overlap between atoms.
                real lik = r - scaledRadius;
                real lik2 = lik * lik;
                real lik4 = lik2 * lik2;
                de = de + 0.25f * M_PI * (sk2 - 4.0f * scaledRadius * r + r2) / (r2 * lik4);
            }

            // Upper integration bound is always the same.
            real uik = r + scaledRadius;
            real uik2 = uik * uik;
            real uik4 = uik2 * uik2;
            de = de - 0.25f * M_PI * (sk2 + 4.0f * scaledRadius * r + r2) / (r2 * uik4);
        }
        return de;
}

DEVICE real neckDescreenDerivative(real r, real radius, real radiusK, real sneck,
                                   GLOBAL const float* RESTRICT neckRadii,
                                   GLOBAL const float* RESTRICT neckA,
                                   GLOBAL const float* RESTRICT neckB) {

    real radiusWater = 0.14f;

    if (r > radius + radiusK + 2 * radiusWater) {
        return 0;
    }

    // Get Aij and Bij
    real Aij, Bij;
    getNeckConstants(radius, radiusK, &Aij, &Bij, neckRadii, neckA, neckB);

    // Use Aij and Bij to get neck value using derivative of Equation 13 from Aguilar/Onufriev 2010 paper
    real rMinusBij = r - Bij;
    real rMinusBij3 = rMinusBij * rMinusBij * rMinusBij;
    real rMinusBij4 = rMinusBij3 * rMinusBij;
    real radiiMinusr = radius + radiusK + 2 * radiusWater - r;
    real radiiMinusr3 = radiiMinusr * radiiMinusr * radiiMinusr;
    real radiiMinusr4 = radiiMinusr3 * radiiMinusr;

    real PI4_3 = 4 * M_PI / (real) 3;

    return 4 * PI4_3 * (sneck * Aij * rMinusBij3 * radiiMinusr4 - sneck * Aij * rMinusBij4 * radiiMinusr3);
}

DEVICE void computeBornChainRuleInteraction(AtomData3 atom1, AtomData3 atom2, real3* force,
                                            GLOBAL const float* RESTRICT neckRadii,
                                            GLOBAL const float* RESTRICT neckA,
                                            GLOBAL const float* RESTRICT neckB) {
    real third = 1 / (real) 3;
    real pi43 = 4 * third * M_PI;
    real factor = -POW(M_PI, third)*POW((real) 6.0f, (real) 2 * third) / (real) 9;
    real term = pi43/(atom1.bornRadius*atom1.bornRadius*atom1.bornRadius);
    term = factor/POW(term, (real) 4 * third);

    if (TANH_RESCALING) {
        real bornSum = atom1.bornSum * pi43;
        real rhoi = atom1.radius;
        real rhoi3 = rhoi * rhoi * rhoi;
        real rhoi3Psi = rhoi3 * bornSum;
        real rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
        real rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
        real rhoi6Psi = rhoi3 * rhoi3 * bornSum;
        real rhoi9Psi2 = rhoi6Psi2 * rhoi3;
        real tanhTerm = tanh(BETA0 * rhoi3Psi - BETA1 * rhoi6Psi2 + BETA2 * rhoi9Psi3);
        real tanh2 = tanhTerm * tanhTerm;
        real chainRuleTerm = BETA0 * rhoi3 - 2 * BETA1 * rhoi6Psi + 3 * BETA2 * rhoi9Psi2;
        real recipBigRad = 1 / (real) BIG_RADIUS;
        real recipBigRad3 = recipBigRad * recipBigRad * recipBigRad;
        real tanh_constant = pi43 * ((1 / (real) rhoi3) - recipBigRad3);
        term = term * tanh_constant * chainRuleTerm * (1.0f - tanh2);
    }

    real de = 0.0f;
    real3 delta = atom2.pos-atom1.pos;
    real sk = atom2.scaleFactor;
    if (sk <= 0 || atom1.bornRadius >= BIG_RADIUS || atom2.descreenRadius <= 0) {
        *force = delta * de;
        return;
    }

    real baseRadius = max(atom1.radius, atom1.descreenRadius) + DESCREEN_OFFSET;
    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    de = pairIntegralDerivative(r, r2, baseRadius, sk * atom2.descreenRadius);

    real mixedNeckScale = 0.5f * (atom1.neckFactor + atom2.neckFactor);
    if (mixedNeckScale > 0) {
        de += neckDescreenDerivative(r, baseRadius, atom2.descreenRadius, mixedNeckScale, neckRadii, neckA, neckB);
    }

    real dbr = term*de/r;
    de = dbr*atom1.bornForce;
    *force = delta*de;
}

/**
 * Compute chain rule terms.
 */
KERNEL void computeChainRuleForce(
        GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL const real4* RESTRICT posq, unsigned int startTileIndex, unsigned int numTileIndices,
        GLOBAL const float4* RESTRICT params,
        GLOBAL const float* RESTRICT neckRadii,
        GLOBAL const float* RESTRICT neckA,
        GLOBAL const float* RESTRICT neckB,
        GLOBAL const mm_long* RESTRICT bornSum,
        GLOBAL const real* RESTRICT bornRadii, GLOBAL const mm_long* RESTRICT bornForce) {
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
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            AtomData3 data = loadAtomData3(atom1, posq, params, bornSum, bornRadii, bornForce);
            data.force = make_real3(0);
            if (pos >= end)
                ; // This warp is done.
            else if (x == y) {
                // This tile is on the diagonal.
                localData[LOCAL_ID].pos = data.pos;
                localData[LOCAL_ID].radius = data.radius;
                localData[LOCAL_ID].scaleFactor = data.scaleFactor;
                localData[LOCAL_ID].descreenRadius = data.descreenRadius;
                localData[LOCAL_ID].neckFactor = data.neckFactor;
                localData[LOCAL_ID].bornRadius = data.bornRadius;
                localData[LOCAL_ID].bornSum = data.bornSum;
                localData[LOCAL_ID].bornForce = data.bornForce;
                localData[LOCAL_ID].force = make_real3(0);
                SYNC_WARPS;
                
                // Compute forces.
                
                for (unsigned int j = (tgx+1)&(TILE_SIZE-1); j != tgx; j = (j+1)&(TILE_SIZE-1)) {
                    int atom2 = y*TILE_SIZE+j;
                    if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        computeBornChainRuleInteraction(data, localData[tbx+j], &tempForce, neckRadii, neckA, neckB);
                        data.force -= tempForce;
                        localData[tbx+j].force += tempForce;
                    }
                    SYNC_WARPS;
                }
                ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) realToFixedPoint((data.force.x+localData[LOCAL_ID].force.x)));
                ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint((data.force.y+localData[LOCAL_ID].force.y)));
                ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint((data.force.z+localData[LOCAL_ID].force.z)));
            }
            else {
                // This is an off-diagonal tile.

                unsigned int j = y*TILE_SIZE + tgx;
                localData[LOCAL_ID] = loadAtomData3(j, posq, params, bornSum, bornRadii, bornForce);
                localData[LOCAL_ID].force = make_real3(0);
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = y*TILE_SIZE+tj;
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        real3 tempForce;
                        computeBornChainRuleInteraction(data, localData[tbx+tj], &tempForce, neckRadii, neckA, neckB);
                        data.force -= tempForce;
                        localData[tbx+tj].force += tempForce;
                        computeBornChainRuleInteraction(localData[tbx+tj], data, &tempForce, neckRadii, neckA, neckB);
                        data.force += tempForce;
                        localData[tbx+tj].force -= tempForce;
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
                if (pos < end) {
                    unsigned int offset = x*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
                    offset = y*TILE_SIZE + tgx;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
                }
            }
        }
        pos++;
    } while (pos < end);
}

#if defined(USE_HIP)
    #define ALIGN alignas(16)
#else
    #define ALIGN
#endif

typedef struct ALIGN {
    real3 pos;
    real q;
    real3 dipole;
#if defined(USE_HIP)
    real padding0;
#endif
    real3 inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS;
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ, quadrupoleZZ;
    real3 force;
    float thole, damp;
#if defined(USE_HIP) && !defined(USE_DOUBLE_PRECISION)
    real padding1[2]; // Prevent bank conflicts because the aligned size is 128
#endif
} AtomData4;

#if defined(USE_HIP)
#define ATOM2_PTR_ARG_SPEC const
#else
#define ATOM2_PTR_ARG_SPEC volatile
#endif

DEVICE void computeOneEDiffInteractionF1(const AtomData4* atom1, LOCAL_ARG ATOM2_PTR_ARG_SPEC AtomData4* atom2, float dScale, float pScale, real* outputEnergy, real3* outputForce);
DEVICE void computeOneEDiffInteractionT1(const AtomData4* atom1, LOCAL_ARG ATOM2_PTR_ARG_SPEC AtomData4* atom2, float dScale, float pScale, real3* outputForce);
DEVICE void computeOneEDiffInteractionT3(const AtomData4* atom1, LOCAL_ARG ATOM2_PTR_ARG_SPEC AtomData4* atom2, float dScale, float pScale, real3* outputForce);

inline DEVICE AtomData4 loadAtomData4(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT labFrameDipole,
        GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT inducedDipolePolar,
        GLOBAL const real* RESTRICT inducedDipoleS, GLOBAL const real* RESTRICT inducedDipolePolarS, GLOBAL const float2* RESTRICT dampingAndThole) {
    AtomData4 data;
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
    data.inducedDipole = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
    data.inducedDipolePolar = make_real3(inducedDipolePolar[3*atom], inducedDipolePolar[3*atom+1], inducedDipolePolar[3*atom+2]);
    data.inducedDipoleS = make_real3(inducedDipoleS[3*atom], inducedDipoleS[3*atom+1], inducedDipoleS[3*atom+2]);
    data.inducedDipolePolarS = make_real3(inducedDipolePolarS[3*atom], inducedDipolePolarS[3*atom+1], inducedDipolePolarS[3*atom+2]);
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
    return data;
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
        GLOBAL const real* RESTRICT labFrameDipole, GLOBAL const real* RESTRICT labFrameQuadrupole, GLOBAL const real* RESTRICT inducedDipole,
        GLOBAL const real* RESTRICT inducedDipolePolar, GLOBAL const real* RESTRICT inducedDipoleS, GLOBAL const real* RESTRICT inducedDipolePolarS,
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
        unsigned int atom1 = x*TILE_SIZE + tgx;
        AtomData4 data = loadAtomData4(atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
        data.force = make_real3(0);
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

            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneEDiffInteractionF1(&data, &localData[tbx+j], d, p, &tempEnergy, &tempForce);
                    energy += 0.25f*tempEnergy;
                    data.force += tempForce;
                }
            }
            SYNC_WARPS;
            data.force *= ENERGY_SCALE_FACTOR;
            ATOMIC_ADD(&forceBuffers[atom1], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&forceBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));

            // Compute torques.

            data.force = make_real3(0);
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneEDiffInteractionT1(&data, &localData[tbx+j], d, p, &tempTorque);
                    data.force += tempTorque;
                }
            }
            data.force *= ENERGY_SCALE_FACTOR;
            ATOMIC_ADD(&torqueBuffers[atom1], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&torqueBuffers[atom1+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            SYNC_WARPS;
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            localData[LOCAL_ID] = loadAtomData4(j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[LOCAL_ID].force = make_real3(0);
            SYNC_WARPS;

            // Compute forces.

            unsigned int tj = tgx;
            SYNC_WARPS;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneEDiffInteractionF1(&data, &localData[tbx+tj], d, p, &tempEnergy, &tempForce);
                    energy += 0.5f*tempEnergy;
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));

            // Compute torques.

            data.force = make_real3(0);
            localData[LOCAL_ID].force = make_real3(0);
            SYNC_WARPS;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneEDiffInteractionT1(&data, &localData[tbx+tj], d, p, &tempTorque);
                    data.force += tempTorque;
                    computeOneEDiffInteractionT3(&data, &localData[tbx+tj], d, p, &tempTorque);
                    localData[tbx+tj].force += tempTorque;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
            SYNC_WARPS;
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

        SYNC_WARPS;
        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            SYNC_WARPS;
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[LOCAL_ID] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[LOCAL_ID] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
            SYNC_WARPS;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        bool includeTile = (skipTiles[currentSkipIndex] != pos);
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData4 data = loadAtomData4(atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            data.force = make_real3(0);
            localData[LOCAL_ID] = loadAtomData4(atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            unsigned int j = y*TILE_SIZE + tgx;
            localData[LOCAL_ID] = loadAtomData4(j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, inducedDipoleS, inducedDipolePolarS, dampingAndThole);
            localData[LOCAL_ID].force = make_real3(0);
            SYNC_WARPS;

            // Compute forces.

            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    computeOneEDiffInteractionF1(&data, &localData[tbx+tj], 1, 1, &tempEnergy, &tempForce);
                    energy += 0.5f*tempEnergy;
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));

            // Compute torques.

            data.force = make_real3(0);
            localData[LOCAL_ID].force = make_real3(0);
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempTorque;
                    computeOneEDiffInteractionT1(&data, &localData[tbx+tj], 1, 1, &tempTorque);
                    data.force += tempTorque;
                    computeOneEDiffInteractionT3(&data, &localData[tbx+tj], 1, 1, &tempTorque);
                    localData[tbx+tj].force += tempTorque;
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
            data.force *= ENERGY_SCALE_FACTOR;
            localData[LOCAL_ID].force *= ENERGY_SCALE_FACTOR;
            offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            offset = y*TILE_SIZE + tgx;
            ATOMIC_ADD(&torqueBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
            ATOMIC_ADD(&torqueBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
            ATOMIC_ADD(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy*ENERGY_SCALE_FACTOR;
}
