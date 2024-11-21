#define TILE_SIZE 32
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos, force;
    float radius, epsilon, padding;
} AtomData;

inline DEVICE AtomData loadAtomData(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const float2* RESTRICT radiusEpsilon) {
    AtomData data;
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    float2 temp = radiusEpsilon[atom];
    data.radius = temp.x;
    data.epsilon = temp.y;
    return data;
}

DEVICE void initParticleParameters(float radius, float epsilon, real* rmixo, real* rmixh, real* emixo, real* emixh) {
    real sqrtEps = SQRT(epsilon);
    real denominator = SQRT(EPSO) + sqrtEps;
    *emixo = 4*EPSO*epsilon / (denominator*denominator);
    denominator = SQRT(EPSH) + sqrtEps;
    *emixh = 4*EPSH*epsilon / (denominator*denominator);
    real radius2 = radius*radius;
    real rmino2 = RMINO*RMINO; 
    *rmixo = 2*(rmino2*RMINO + radius2*radius) / (rmino2 + radius2);
    real rminh2 = RMINH*RMINH;
    *rmixh = 2*(rminh2*RMINH + radius2*radius) / (rminh2+radius2);
}

DEVICE real integralBeforeRMin(real eps, real r, real r2, real sk2,
        real lik2, real lik3, real lik4, real uik2, real uik3, real uik4) {
        return -eps * (4 * M_PI / (48 * r) * (3 * (lik4 - uik4) - 8 * r * (lik3 - uik3) + 6 * (r2 - sk2) * (lik2 - uik2)));
}

DEVICE real integralBeforeRminDerivative(real ri, real eps, real rmin, real r, real r2,
                                         real r3, real sk, real sk2, real lik, real lik2,
                                         real lik3, real uik, real uik2, real uik3) {
    real dl;
    if (ri > r - sk)
        dl = (-lik2 + 2 * r2 + 2 * sk2) * lik2;
    else
        dl = (-lik3 + 4 * lik2 * r - 6 * lik * r2 + 2 * lik * sk2 + 4 * r3 - 4 * r * sk2) * lik;
    real du;
    if (r + sk > rmin)
        du = -(-uik2 + 2 * r2 + 2 * sk2) * uik2;
    else
        du = -(-uik3 + 4 * uik2 * r - 6 * uik * r2 + 2 * uik * sk2 + 4 * r3 - 4 * r * sk2) * uik;
    return -eps * M_PI * (dl + du) / (4 * r2);
}

DEVICE real integralAfterRmin(real eps, real rmin7, real r, real r2, real sk2,
                               real lik, real lik2, real lik3, real lik4, real lik5, real lik10,
                                real lik11, real lik12, real uik, real uik2, real uik3, real uik4,
                                real uik5, real uik10, real uik11, real uik12) {
    real er7 = eps * rmin7;
    real term = 4 * M_PI / (120 * r * lik5 * uik5)
            * (15 * uik * lik * r * (uik4 - lik4)
               - 10 * uik2 * lik2 * (uik3 - lik3)
               + 6 * (sk2 - r2) * (uik5 - lik5));
    real term2 = 4 * M_PI / (2640 * r * lik12 * uik12)
             * (120 * uik * lik * r * (uik11 - lik11)
                - 66 * uik2 * lik2 * (uik10 - lik10)
                + 55 * (sk2 - r2) * (uik12 - lik12));
    real idisp = -2 * er7 * term;
    real irep = er7 * rmin7 * term2;
    return irep + idisp;
}

DEVICE real integralAfterRminDerivative(real ri, real eps, real rmin, real rmin7, real rmax,
                                        real r, real r2, real r3, real sk, real sk2, real lik,
                                        real lik2, real lik3, real lik5, real lik6, real lik12,
                                        real lik13, real uik, real uik2, real uik3, real uik6,
                                        real uik13) {
    real er7 = eps * rmin7;
    real lowerTerm = lik2 * r + r3 - r * sk2;
    real upperTerm = uik2 * r + r3 - r * sk2;

    real dl;
    if (ri > r - sk || rmax < rmin)
        dl = -(-5 * lik2 + 3 * r2 + 3 * sk2) / lik5;
    else
        dl = (5 * lik3 - 33 * lik * r2 - 3 * lik * sk2 + 15 * lowerTerm) / lik6;
    real du = -(5 * uik3 - 33 * uik * r2 - 3 * uik * sk2 + 15 * upperTerm) / uik6;
    real de = -2 * M_PI * er7 * (dl + du) / (15 * r2);

    if (ri > r - sk || rmax < rmin)
        dl = -(-6 * lik2 + 5 * r2 + 5 * sk2) / lik12;
    else
        dl = (6 * lik3 - 125 * lik * r2 - 5 * lik * sk2 + 60 * lowerTerm) / lik13;
    du = -(6 * uik3 - 125 * uik * r2 - 5 * uik * sk2 + 60 * upperTerm) / uik13;
    de += M_PI * er7 * rmin7 * (dl + du) / (60 * r2);
    return de;
}

DEVICE real interact(real factor, real ri, real sk, real rmix, real emix,
                     real r, real r2, real r3, real3 *force) {
    real sum = 0;
    // Nothing to do if the integral begins beyond r + sk (i.e. atom k does not exclude solvent)
    if (ri < r + sk) {
        // Zero out the derivative contribution of atom k.
        real de = 0;
        real sk2 = sk * sk;
        // Compute the maximum of 1) the beginning of the integral and 2) closest edge of atom K.
        real iStart = ri > r - sk ? ri : r - sk;
        // Use this as the lower limit for integrating the constant eps value below Rmin.
        real lik = iStart;
        // Interaction with water from lik to Rmin; nothing to do if the lower limit is greater than Rmin.
        if (lik < rmix) {
            real lik2 = lik * lik;
            real lik3 = lik2 * lik;
            real lik4 = lik3 * lik;
            // Upper limit is the minimum of Rmin and the farthest edge of atom K.
            real uik = r + sk < rmix ? r + sk : rmix;
            real uik2 = uik * uik;
            real uik3 = uik2 * uik;
            real uik4 = uik3 * uik;
            sum = integralBeforeRMin(emix, r, r2, sk2, lik2, lik3, lik4, uik2, uik3, uik4);
            de = integralBeforeRminDerivative(ri, emix, rmix, r, r2, r3, sk, sk2, lik, lik2, lik3, uik, uik2, uik3);
        }
        // Upper limit the variable part of Uwca always the farthest edge of atom K.
        real uik = r + sk;
        // Interaction with water beyond Rmin, from lik to uik = r + sk.
        if (uik > rmix) {
            // Start the integral at the max of 1) iStart and 2) Rmin.
            lik = iStart > rmix ? iStart : rmix;
            real lik2 = lik * lik;
            real lik3 = lik2 * lik;
            real lik4 = lik3 * lik;
            real lik5 = lik4 * lik;
            real lik6 = lik5 * lik;
            real lik10 = lik5 * lik5;
            real lik11 = lik10 * lik;
            real lik12 = lik11 * lik;
            real uik2 = uik * uik;
            real uik3 = uik2 * uik;
            real uik4 = uik3 * uik;
            real uik5 = uik4 * uik;
            real uik10 = uik5 * uik5;
            real uik11 = uik10 * uik;
            real uik12 = uik11 * uik;
            real rmix3 = rmix * rmix * rmix;
            real rmix7 = rmix3 * rmix3 * rmix;
            sum += integralAfterRmin(emix, rmix7, r, r2, sk2, lik, lik2, lik3, lik4, lik5, lik10, lik11, lik12, uik,
                    uik2, uik3, uik4, uik5, uik10, uik11, uik12);
            real lik13 = lik12 * lik;
            real uik6 = uik5 * uik;
            real uik13 = uik12 * uik;
            de += integralAfterRminDerivative(ri, emix, rmix, rmix7, iStart, r, r2, r3, sk, sk2, lik, lik2, lik3,
                    lik5, lik6, lik12, lik13, uik, uik2, uik3, uik6, uik13);
        }
        // Increment the individual dispersion gradient components.
        de *= factor / r;
        (*force).x += de;
        (*force).y += de;
        (*force).z += de;
    }
    return factor * sum;
}


DEVICE void computeOneInteraction(AtomData atom1, AtomData atom2, real rmixo, real rmixh, real emixo, real emixh, real3* force, real* energy) {
    // get deltaR and r between 2 atoms

    *force = atom1.pos - atom2.pos;
    real r2 = dot(*force, *force);
    if (r2 <= 0) {
        *force = make_real3(0);
        *energy = 0;
        return;
    }
    real rI = RSQRT(r2);
    real r = RECIP(rI);

    real xr = (*force).x;
    real yr = (*force).y;
    real zr = (*force).z;
    real r3 = r2 * r;

    real sk = atom2.radius * SHCTD;

    // Start of integration of dispersion for atom i with water oxygen.
    real riO = rmixo * 0.5f + DISPOFF;
    real nO = 1.0f;

    // Start of integration of dispersion for atom i with water hydrogen.
    real riH = rmixh * 0.5f + DISPOFF;
    real nH = 2.0f;

    *force = make_real3(0);
    *energy = interact(nO, riO, sk, rmixo, emixo, r, r2, r3, force) +
             interact(nH, riH, sk, rmixh, emixh, r, r2, r3, force);

    (*force).x *= AWATER * xr;
    (*force).y *= AWATER * yr;
    (*force).z *= AWATER * zr;
}

/**
 * Compute WCA interaction.
 */
KERNEL void computeWCAForce(GLOBAL mm_ulong* RESTRICT forceBuffers, GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL const real4* RESTRICT posq, unsigned int startTileIndex, unsigned int numTileIndices, GLOBAL const float2* RESTRICT radiusEpsilon) {
    unsigned int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    unsigned int warp = GLOBAL_ID/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = (unsigned int) (startTileIndex+warp*(mm_long)numTiles/totalWarps);
    unsigned int end = (unsigned int) (startTileIndex+(warp+1)*(mm_long)numTiles/totalWarps);
    mixed energy = 0;
    LOCAL AtomData localData[THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        
        const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
        const unsigned int tbx = LOCAL_ID - tgx;
        const unsigned int localGroupIndex = LOCAL_ID/TILE_SIZE;
        int x, y;
        AtomData data;
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            data = loadAtomData(atom1, posq, radiusEpsilon);
            localData[LOCAL_ID] = loadAtomData(y*TILE_SIZE+tgx, posq, radiusEpsilon);
            real emixo, emixh, rmixo, rmixh;
            initParticleParameters(data.radius, data.epsilon, &rmixo, &rmixh, &emixo, &emixh);
            data.force = make_real3(0);
            localData[LOCAL_ID].force = make_real3(0);
            SYNC_WARPS;

            // Compute forces.

            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    computeOneInteraction(data, localData[tbx+tj], rmixo, rmixh, emixo, emixh, &tempForce, &tempEnergy);
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                    energy += (x == y ? 0.5f*tempEnergy : tempEnergy);
                    real emjxo, emjxh, rmjxo, rmjxh;
                    initParticleParameters(localData[tbx+tj].radius, localData[tbx+tj].epsilon, &rmjxo, &rmjxh, &emjxo, &emjxh);
                    computeOneInteraction(localData[tbx+tj], data, rmjxo, rmjxh, emjxo, emjxh, &tempForce, &tempEnergy);
                    data.force -= tempForce;
                    localData[tbx+tj].force += tempForce;
                    energy += (x == y ? 0.5f*tempEnergy : tempEnergy);
                }
                tj = (tj+1) & (TILE_SIZE-1);
                SYNC_WARPS;
            }
            unsigned int offset = x*TILE_SIZE + tgx;
            ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(data.force.x));
            ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.y));
            ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(data.force.z));
            if (x != y) {
                offset = y*TILE_SIZE + tgx;
                ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.x));
                ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.y));
                ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].force.z));
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[GLOBAL_ID] -= AWATER*energy;
}
