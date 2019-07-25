#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real4 posq;
    real3 field, fieldPolar, dipole;
#ifdef INCLUDE_QUADRUPOLES
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ;
    real quadrupoleYY, quadrupoleYZ, quadrupoleZZ;
#endif
    float thole, damp;
#ifdef USE_GK
    real3 gkField;
    real bornRadius;
#endif
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const float2* __restrict__ dampingAndThole) {
    data.posq = posq[atom];
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
#ifdef INCLUDE_QUADRUPOLES
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
    data.quadrupoleZZ = -(data.quadrupoleXX+data.quadrupoleYY);
#endif
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real r2 = dot(deltaR, deltaR);
    if (r2 <= CUTOFF_SQUARED) {
        // calculate the error function damping terms

        real r = SQRT(r2);
        real ralpha = EWALD_ALPHA*r;
        real exp2a = EXP(-(ralpha*ralpha));
#ifdef USE_DOUBLE_PRECISION
        const real erfcAlphaR = erfc(ralpha);
#else
        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 1.5e-7.

        const real t = RECIP(1.0f+0.3275911f*ralpha);
        const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
        real bn0 = erfcAlphaR/r;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn3 = (5*bn2+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        real scale3 = 1;
        real scale5 = 1;
        real scale7 = 1;
        real damp = atom1.damp*atom2.damp;
        if (damp != 0) {
            real ratio = (r/damp);
            ratio = ratio*ratio*ratio;
            real pgamma = (atom1.thole < atom2.thole ? atom1.thole : atom2.thole);
            damp = -pgamma*ratio;
            if (damp > -50) {
                real expdamp = EXP(damp);
                scale3 = 1 - expdamp;
                scale5 = 1 - expdamp*(1-damp);
                scale7 = 1 - expdamp*(1-damp+(0.6f*damp*damp));
            }
        }
        real dsc3 = dScale*scale3;
        real dsc5 = dScale*scale5;
        real dsc7 = dScale*scale7;

        real psc3 = pScale*scale3;
        real psc5 = pScale*scale5;
        real psc7 = pScale*scale7;

        real r3 = r*r2;
        real r5 = r3*r2;
        real r7 = r5*r2;
        real drr3 = (1-dsc3)/r3;
        real drr5 = 3*(1-dsc5)/r5;
        real drr7 = 15*(1-dsc7)/r7;

        real prr3 = (1-psc3)/r3;
        real prr5 = 3*(1-psc5)/r5;
        real prr7 = 15*(1-psc7)/r7;

        real dir = dot(atom1.dipole, deltaR);
        real dkr = dot(atom2.dipole, deltaR);

#ifdef INCLUDE_QUADRUPOLES
        real3 qi;
        qi.x = atom1.quadrupoleXX*deltaR.x + atom1.quadrupoleXY*deltaR.y + atom1.quadrupoleXZ*deltaR.z;
        qi.y = atom1.quadrupoleXY*deltaR.x + atom1.quadrupoleYY*deltaR.y + atom1.quadrupoleYZ*deltaR.z;
        qi.z = atom1.quadrupoleXZ*deltaR.x + atom1.quadrupoleYZ*deltaR.y + atom1.quadrupoleZZ*deltaR.z;
        real qir = dot(qi, deltaR);

        real3 qk;
        qk.x = atom2.quadrupoleXX*deltaR.x + atom2.quadrupoleXY*deltaR.y + atom2.quadrupoleXZ*deltaR.z;
        qk.y = atom2.quadrupoleXY*deltaR.x + atom2.quadrupoleYY*deltaR.y + atom2.quadrupoleYZ*deltaR.z;
        qk.z = atom2.quadrupoleXZ*deltaR.x + atom2.quadrupoleYZ*deltaR.y + atom2.quadrupoleZZ*deltaR.z;
        real qkr = dot(qk, deltaR);

        real3 fim = -deltaR*(bn1*atom2.posq.w-bn2*dkr+bn3*qkr) - bn1*atom2.dipole + 2*bn2*qk;
        real3 fkm = deltaR*(bn1*atom1.posq.w+bn2*dir+bn3*qir) - bn1*atom1.dipole - 2*bn2*qi;
        real3 fid = -deltaR*(drr3*atom2.posq.w-drr5*dkr+drr7*qkr) - drr3*atom2.dipole + 2*drr5*qk;
        real3 fkd = deltaR*(drr3*atom1.posq.w+drr5*dir+drr7*qir) - drr3*atom1.dipole - 2*drr5*qi;
        real3 fip = -deltaR*(prr3*atom2.posq.w-prr5*dkr+prr7*qkr) - prr3*atom2.dipole + 2*prr5*qk;
        real3 fkp = deltaR*(prr3*atom1.posq.w+prr5*dir+prr7*qir) - prr3*atom1.dipole - 2*prr5*qi;
#else
        real3 fim = -deltaR*(bn1*atom2.posq.w-bn2*dkr) - bn1*atom2.dipole;
        real3 fkm = deltaR*(bn1*atom1.posq.w+bn2*dir) - bn1*atom1.dipole;
        real3 fid = -deltaR*(drr3*atom2.posq.w-drr5*dkr) - drr3*atom2.dipole;
        real3 fkd = deltaR*(drr3*atom1.posq.w+drr5*dir) - drr3*atom1.dipole;
        real3 fip = -deltaR*(prr3*atom2.posq.w-prr5*dkr) - prr3*atom2.dipole;
        real3 fkp = deltaR*(prr3*atom1.posq.w+prr5*dir) - prr3*atom1.dipole;
#endif
        // increment the field at each site due to this interaction

        fields[0] = fim-fid;
        fields[1] = fim-fip;
        fields[2] = fkm-fkd;
        fields[3] = fkm-fkp;
    }
    else {
        fields[0] = make_real3(0);
        fields[1] = make_real3(0);
        fields[2] = make_real3(0);
        fields[3] = make_real3(0);
    }
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;

    real rr3 = rI*r2I;
    real rr5 = 3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp = atom1.damp*atom2.damp;
    real dampExp;
    if (damp != 0) {

        // get scaling factors
      
        real ratio = r/damp;
        float pGamma = atom2.thole > atom1.thole ? atom1.thole : atom2.thole; 
        damp = ratio*ratio*ratio*pGamma;
        dampExp = EXP(-damp);
    }
    else
        dampExp = 0;
      
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    rr7 *= 1 - (1+damp+(0.6f*damp*damp))*dampExp;
      
    real rr5_2 = 2*rr5;
 
    real dir = dot(atom1.dipole, deltaR);
    real dkr = dot(atom2.dipole, deltaR);

#ifdef INCLUDE_QUADRUPOLES
    real3 qi;
    qi.x = atom1.quadrupoleXX*deltaR.x + atom1.quadrupoleXY*deltaR.y + atom1.quadrupoleXZ*deltaR.z;
    qi.y = atom1.quadrupoleXY*deltaR.x + atom1.quadrupoleYY*deltaR.y + atom1.quadrupoleYZ*deltaR.z;
    qi.z = atom1.quadrupoleXZ*deltaR.x + atom1.quadrupoleYZ*deltaR.y + atom1.quadrupoleZZ*deltaR.z;
    real qir = dot(qi, deltaR);

    real3 qk;
    qk.x = atom2.quadrupoleXX*deltaR.x + atom2.quadrupoleXY*deltaR.y + atom2.quadrupoleXZ*deltaR.z;
    qk.y = atom2.quadrupoleXY*deltaR.x + atom2.quadrupoleYY*deltaR.y + atom2.quadrupoleYZ*deltaR.z;
    qk.z = atom2.quadrupoleXZ*deltaR.x + atom2.quadrupoleYZ*deltaR.y + atom2.quadrupoleZZ*deltaR.z;
    real qkr = dot(qk, deltaR);

    real factor = -rr3*atom2.posq.w + rr5*dkr - rr7*qkr;
    real3 field1 = deltaR*factor - rr3*atom2.dipole + rr5_2*qk;
    factor = rr3*atom1.posq.w + rr5*dir + rr7*qir;
    real3 field2 = deltaR*factor - rr3*atom1.dipole - rr5_2*qi;
#else
    real factor = -rr3*atom2.posq.w + rr5*dkr;
    real3 field1 = deltaR*factor - rr3*atom2.dipole;
    factor = rr3*atom1.posq.w + rr5*dir;
    real3 field2 = deltaR*factor - rr3*atom1.dipole;
#endif
    fields[0] = dScale*field1;
    fields[1] = pScale*field1;
    fields[2] = dScale*field2;
    fields[3] = pScale*field2;
}
#endif

#ifdef USE_GK
__device__ void computeOneGkInteraction(AtomData& atom1, AtomData& atom2, real3 delta, real3* fields) {
    real a[4][4];
    real gc[5];
    real gux[11],guy[11],guz[11];
    real gqxx[5],gqxy[5];
    real gqxz[5],gqyy[5];
    real gqyz[5],gqzz[5];

    real ci = atom1.posq.w;
    real ck = atom2.posq.w;

    real uxi = atom1.dipole.x;
    real uyi = atom1.dipole.y;
    real uzi = atom1.dipole.z;
    real uxk = atom2.dipole.x;
    real uyk = atom2.dipole.y;
    real uzk = atom2.dipole.z;

#ifdef INCLUDE_QUADRUPOLES
    real qxxi = atom1.quadrupoleXX;
    real qxyi = atom1.quadrupoleXY;
    real qxzi = atom1.quadrupoleXZ;
    real qyyi = atom1.quadrupoleYY;
    real qyzi = atom1.quadrupoleYZ;
    real qzzi = atom1.quadrupoleZZ;
    real qxxk = atom2.quadrupoleXX;
    real qxyk = atom2.quadrupoleXY;
    real qxzk = atom2.quadrupoleXZ;
    real qyyk = atom2.quadrupoleYY;
    real qyzk = atom2.quadrupoleYZ;
    real qzzk = atom2.quadrupoleZZ;
#else
    real qxxi = 0;
    real qxyi = 0;
    real qxzi = 0;
    real qyyi = 0;
    real qyzi = 0;
    real qzzi = 0;
    real qxxk = 0;
    real qxyk = 0;
    real qxzk = 0;
    real qyyk = 0;
    real qyzk = 0;
    real qzzk = 0;
#endif
    real xr2 = delta.x*delta.x;
    real yr2 = delta.y*delta.y;
    real zr2 = delta.z*delta.z;
    real r2 = xr2 + yr2 + zr2;

    real rb2 = atom1.bornRadius*atom2.bornRadius;
    real expterm = EXP(-r2/(GK_C*rb2));
    real expc = expterm / GK_C;
    real dexpc = -2/(GK_C*rb2);
    real gf2 = RECIP(r2+rb2*expterm);
    real gf = SQRT(gf2);
    real gf3 = gf2*gf;
    real gf5 = gf3*gf2;
    real gf7 = gf5*gf2;

    // reaction potential auxiliary terms

    a[0][0] = gf;
    a[1][0] = -gf3;
    a[2][0] = 3*gf5;
    a[3][0] = -15*gf7;

    // reaction potential gradient auxiliary terms

    real expc1 = 1 - expc;
    a[0][1] = expc1*a[1][0];
    a[1][1] = expc1*a[2][0];
    a[2][1] = expc1*a[3][0];

    // dipole second reaction potential gradient auxiliary term

    real expcdexpc = -expc*dexpc;
    a[1][2] = expc1*a[2][1] + expcdexpc*a[2][0];

    // multiply the auxillary terms by dielectric functions;

    a[0][1] = GK_FC*a[0][1];
    a[1][0] = GK_FD*a[1][0];
    a[1][1] = GK_FD*a[1][1];
    a[1][2] = GK_FD*a[1][2];
    a[2][0] = GK_FQ*a[2][0];
    a[2][1] = GK_FQ*a[2][1];

    // unweighted dipole reaction potential tensor

    gux[1] = delta.x*a[1][0];
    guy[1] = delta.y*a[1][0];
    guz[1] = delta.z*a[1][0];

    // unweighted reaction potential gradient tensor

    gc[2] = delta.x*a[0][1];
    gc[3] = delta.y*a[0][1];
    gc[4] = delta.z*a[0][1];
    gux[2] = a[1][0] + xr2*a[1][1];
    gux[3] = delta.x*delta.y*a[1][1];
    gux[4] = delta.x*delta.z*a[1][1];
    guy[2] = gux[3];
    guy[3] = a[1][0] + yr2*a[1][1];
    guy[4] = delta.y*delta.z*a[1][1];
    guz[2] = gux[4];
    guz[3] = guy[4];
    guz[4] = a[1][0] + zr2*a[1][1];
    gqxx[2] = delta.x*(2*a[2][0]+xr2*a[2][1]);
    gqxx[3] = delta.y*xr2*a[2][1];
    gqxx[4] = delta.z*xr2*a[2][1];
    gqyy[2] = delta.x*yr2*a[2][1];
    gqyy[3] = delta.y*(2*a[2][0]+yr2*a[2][1]);
    gqyy[4] = delta.z*yr2*a[2][1];
    gqzz[2] = delta.x*zr2*a[2][1];
    gqzz[3] = delta.y*zr2*a[2][1];
    gqzz[4] = delta.z*(2*a[2][0]+zr2*a[2][1]);
    gqxy[2] = delta.y*(a[2][0]+xr2*a[2][1]);
    gqxy[3] = delta.x*(a[2][0]+yr2*a[2][1]);
    gqxy[4] = delta.z*delta.x*delta.y*a[2][1];
    gqxz[2] = delta.z*(a[2][0]+xr2*a[2][1]);
    gqxz[3] = gqxy[4];
    gqxz[4] = delta.x*(a[2][0]+zr2*a[2][1]);
    gqyz[2] = gqxy[4];
    gqyz[3] = delta.z*(a[2][0]+yr2*a[2][1]);
    gqyz[4] = delta.y*(a[2][0]+zr2*a[2][1]);

    // unweighted dipole second reaction potential gradient tensor

    gux[5] = delta.x*(3*a[1][1]+xr2*a[1][2]);
    gux[6] = delta.y*(a[1][1]+xr2*a[1][2]);
    gux[7] = delta.z*(a[1][1]+xr2*a[1][2]);
    gux[8] = delta.x*(a[1][1]+yr2*a[1][2]);
    gux[9] = delta.z*delta.x*delta.y*a[1][2];
    gux[10] = delta.x*(a[1][1]+zr2*a[1][2]);
    guy[5] = delta.y*(a[1][1]+xr2*a[1][2]);
    guy[6] = delta.x*(a[1][1]+yr2*a[1][2]);
    guy[7] = gux[9];
    guy[8] = delta.y*(3*a[1][1]+yr2*a[1][2]);
    guy[9] = delta.z*(a[1][1]+yr2*a[1][2]);
    guy[10] = delta.y*(a[1][1]+zr2*a[1][2]);
    guz[5] = delta.z*(a[1][1]+xr2*a[1][2]);
    guz[6] = gux[9];
    guz[7] = delta.x*(a[1][1]+zr2*a[1][2]);
    guz[8] = delta.z*(a[1][1]+yr2*a[1][2]);
    guz[9] = delta.y*(a[1][1]+zr2*a[1][2]);
    guz[10] = delta.z*(3*a[1][1]+zr2*a[1][2]);

    // generalized Kirkwood permanent reaction field

    fields[0].x = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                                   + 0.5f*(ck*gux[1] + qxxk*gux[5]
                                   + qyyk*gux[8] + qzzk*gux[10]
                                   + 2*(qxyk*gux[6]+qxzk*gux[7]
                                   + qyzk*gux[9]))
                                   + 0.5f*(ck*gc[2] + qxxk*gqxx[2]
                                   + qyyk*gqyy[2] + qzzk*gqzz[2]
                                   + 2*(qxyk*gqxy[2]+qxzk*gqxz[2]
                                   + qyzk*gqyz[2]));

    fields[0].y = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                                   + 0.5f*(ck*guy[1] + qxxk*guy[5]
                                   + qyyk*guy[8] + qzzk*guy[10]
                                   + 2*(qxyk*guy[6]+qxzk*guy[7]
                                   + qyzk*guy[9]))
                                   + 0.5f*(ck*gc[3] + qxxk*gqxx[3]
                                   + qyyk*gqyy[3] + qzzk*gqzz[3]
                                   + 2*(qxyk*gqxy[3]+qxzk*gqxz[3]
                                   + qyzk*gqyz[3]));

    fields[0].z = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                                   + 0.5f*(ck*guz[1] + qxxk*guz[5]
                                   + qyyk*guz[8] + qzzk*guz[10]
                                   + 2*(qxyk*guz[6]+qxzk*guz[7]
                                   + qyzk*guz[9]))
                                   + 0.5f*(ck*gc[4] + qxxk*gqxx[4]
                                   + qyyk*gqyy[4] + qzzk*gqzz[4]
                                   + 2*(qxyk*gqxy[4]+qxzk*gqxz[4]
                                   + qyzk*gqyz[4]));

    fields[1].x = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                                   - 0.5f*(ci*gux[1] + qxxi*gux[5]
                                   + qyyi*gux[8] + qzzi*gux[10]
                                   + 2*(qxyi*gux[6]+qxzi*gux[7]
                                   + qyzi*gux[9]))
                                   - 0.5f*(ci*gc[2] + qxxi*gqxx[2]
                                   + qyyi*gqyy[2] + qzzi*gqzz[2]
                                   + 2*(qxyi*gqxy[2]+qxzi*gqxz[2]
                                   + qyzi*gqyz[2]));

    fields[1].y = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                                   - 0.5f*(ci*guy[1] + qxxi*guy[5]
                                   + qyyi*guy[8] + qzzi*guy[10]
                                   + 2*(qxyi*guy[6]+qxzi*guy[7]
                                   + qyzi*guy[9]))
                                   - 0.5f*(ci*gc[3]      + qxxi*gqxx[3]
                                   + qyyi*gqyy[3] + qzzi*gqzz[3]
                                   + 2*(qxyi*gqxy[3]+qxzi*gqxz[3]
                                   + qyzi*gqyz[3]));

    fields[1].z = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                                   - 0.5f*(ci*guz[1] + qxxi*guz[5]
                                   + qyyi*guz[8] + qzzi*guz[10]
                                   + 2*(qxyi*guz[6]+qxzi*guz[7]
                                   + qyzi*guz[9]))
                                   - 0.5f*(ci*gc[4] + qxxi*gqxx[4]
                                   + qyyi*gqyy[4] + qzzi*gqzz[4]
                                   + 2*(qxyi*gqxy[4]+qxzi*gqxz[4]
                                   + qyzi*gqyz[4]));
}
#endif

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
 * Compute nonbonded interactions.
 */
extern "C" __global__ void computeFixedField(
        unsigned long long* __restrict__ fieldBuffers, unsigned long long* __restrict__ fieldPolarBuffers, const real4* __restrict__ posq,
        const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags, const ushort2* __restrict__ exclusionTiles,
        unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const unsigned int* __restrict__ interactingAtoms,
#elif defined USE_GK
        const real* __restrict__ bornRadii, unsigned long long* __restrict__ gkFieldBuffers,
#endif
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const float2* __restrict__ dampingAndThole) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        data.field = make_real3(0);
        data.fieldPolar = make_real3(0);
#ifdef USE_GK
        data.gkField = make_real3(0);
#endif
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
#ifdef USE_GK
        data.bornRadius = bornRadii[atom1];
#endif
        uint2 covalent = covalentFlags[pos*TILE_SIZE+tgx];
        unsigned int polarizationGroup = polarizationGroupFlags[pos*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = threadIdx.x;
            localData[localAtomIndex].posq = data.posq;
            localData[localAtomIndex].dipole = data.dipole;
#ifdef INCLUDE_QUADRUPOLES
            localData[localAtomIndex].quadrupoleXX = data.quadrupoleXX;
            localData[localAtomIndex].quadrupoleXY = data.quadrupoleXY;
            localData[localAtomIndex].quadrupoleXZ = data.quadrupoleXZ;
            localData[localAtomIndex].quadrupoleYY = data.quadrupoleYY;
            localData[localAtomIndex].quadrupoleYZ = data.quadrupoleYZ;
            localData[localAtomIndex].quadrupoleZZ = data.quadrupoleZZ;
#endif
            localData[localAtomIndex].thole = data.thole;
            localData[localAtomIndex].damp = data.damp;
#ifdef USE_GK
            localData[localAtomIndex].bornRadius = data.bornRadius;
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+j].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneInteraction(data, localData[tbx+j], delta, d, p, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                }
#ifdef USE_GK
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[2];
                    computeOneGkInteraction(data, localData[tbx+j], delta, fields);
                    data.gkField += fields[0];
                }
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            const unsigned int localAtomIndex = threadIdx.x;
            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData(localData[localAtomIndex], j, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
            localData[localAtomIndex].field = make_real3(0);
            localData[localAtomIndex].fieldPolar = make_real3(0);
#ifdef USE_GK
            localData[localAtomIndex].bornRadius = bornRadii[j];
            localData[localAtomIndex].gkField = make_real3(0);
#endif
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+tj].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneInteraction(data, localData[tbx+tj], delta, d, p, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                    localData[tbx+tj].field += fields[2];
                    localData[tbx+tj].fieldPolar += fields[3];
#ifdef USE_GK
                }
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[2];
                    computeOneGkInteraction(data, localData[tbx+tj], delta, fields);
                    data.gkField += fields[0];
                    localData[tbx+tj].gkField += fields[1];
#endif
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
        }
        
        // Write results.
        
        unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
        atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
        atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
        atomicAdd(&gkFieldBuffers[offset], static_cast<unsigned long long>((long long) (data.gkField.x*0x100000000)));
        atomicAdd(&gkFieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.gkField.y*0x100000000)));
        atomicAdd(&gkFieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.gkField.z*0x100000000)));
#endif
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&gkFieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.x*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.y*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.z*0x100000000)));
#endif
        }
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
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;

    while (pos < end) {
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
#ifdef USE_CUTOFF
        x = tiles[pos];
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

            AtomData data;
            data.field = make_real3(0);
            data.fieldPolar = make_real3(0);
#ifdef USE_GK
            data.gkField = make_real3(0);
#endif
            loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
#ifdef USE_GK
            data.bornRadius = bornRadii[atom1];
#endif
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            const unsigned int localAtomIndex = threadIdx.x;
            loadAtomData(localData[localAtomIndex], j, posq, labFrameDipole, labFrameQuadrupole, dampingAndThole);
            localData[localAtomIndex].field = make_real3(0);
            localData[localAtomIndex].fieldPolar = make_real3(0);
#ifdef USE_GK
            localData[localAtomIndex].bornRadius = bornRadii[j];
            localData[localAtomIndex].gkField = make_real3(0);
#endif

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+tj].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    computeOneInteraction(data, localData[tbx+tj], delta, 1, 1, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                    localData[tbx+tj].field += fields[2];
                    localData[tbx+tj].fieldPolar += fields[3];
#ifdef USE_GK
                }
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[2];
                    computeOneGkInteraction(data, localData[tbx+tj], delta, fields);
                    data.gkField += fields[0];
                    localData[tbx+tj].gkField += fields[1];
#endif
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&gkFieldBuffers[offset], static_cast<unsigned long long>((long long) (data.gkField.x*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.gkField.y*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.gkField.z*0x100000000)));
#endif
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
#ifdef USE_GK
            atomicAdd(&gkFieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.x*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.y*0x100000000)));
            atomicAdd(&gkFieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].gkField.z*0x100000000)));
#endif
        }
        pos++;
    }
}
