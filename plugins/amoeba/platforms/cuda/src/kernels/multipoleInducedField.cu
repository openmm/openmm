#ifndef HIPPO
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos;
    real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
#ifdef EXTRAPOLATED_POLARIZATION
    real fieldGradient[6], fieldGradientPolar[6];
#endif
#ifdef USE_GK
    real3 fieldS, fieldPolarS, inducedDipoleS, inducedDipolePolarS;
    real bornRadius;
    #ifdef EXTRAPOLATED_POLARIZATION
        real fieldGradientS[6], fieldGradientPolarS[6];
    #endif
#endif
    float thole, damp;
} AtomData;

#ifdef USE_GK
inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole, const real* __restrict__ inducedDipoleS,
        const real* __restrict__ inducedDipolePolarS, const real* __restrict__ bornRadii) {
#else
inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
#endif
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
#ifdef USE_GK
    data.inducedDipoleS.x = inducedDipoleS[atom*3];
    data.inducedDipoleS.y = inducedDipoleS[atom*3+1];
    data.inducedDipoleS.z = inducedDipoleS[atom*3+2];
    data.inducedDipolePolarS.x = inducedDipolePolarS[atom*3];
    data.inducedDipolePolarS.y = inducedDipolePolarS[atom*3+1];
    data.inducedDipolePolarS.z = inducedDipolePolarS[atom*3+2];
    data.bornRadius = bornRadii[atom];
#endif
}

inline __device__ void zeroAtomData(AtomData& data) {
    data.field = make_real3(0);
    data.fieldPolar = make_real3(0);
#ifdef USE_GK
    data.fieldS = make_real3(0);
    data.fieldPolarS = make_real3(0);
#endif
#ifdef EXTRAPOLATED_POLARIZATION
    for (int i = 0; i < 6; i++) {
        data.fieldGradient[i] = 0;
        data.fieldGradientPolar[i] = 0;
#ifdef USE_GK
        data.fieldGradientS[i] = 0;
        data.fieldGradientPolarS[i] = 0;
#endif
    }
#endif
}

#ifdef EXTRAPOLATED_POLARIZATION
    #ifdef USE_GK
        #define SAVE_ATOM_DATA(index, data) saveAtomData(index, data, field, fieldPolar, fieldGradient, fieldGradientPolar, fieldS, fieldPolarS, fieldGradientS, fieldGradientPolarS);
    #else
        #define SAVE_ATOM_DATA(index, data) saveAtomData(index, data, field, fieldPolar, fieldGradient, fieldGradientPolar);
    #endif
#else
    #ifdef USE_GK
        #define SAVE_ATOM_DATA(index, data) saveAtomData(index, data, field, fieldPolar, fieldS, fieldPolarS);
    #else
        #define SAVE_ATOM_DATA(index, data) saveAtomData(index, data, field, fieldPolar);
    #endif
#endif

inline __device__ void saveAtomData(int index, AtomData& data, unsigned long long* __restrict__ field, unsigned long long* __restrict__ fieldPolar
#ifdef EXTRAPOLATED_POLARIZATION
        , unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar
#endif
#ifdef USE_GK
        , unsigned long long* __restrict__ fieldS, unsigned long long* __restrict__ fieldPolarS
    #ifdef EXTRAPOLATED_POLARIZATION
        , unsigned long long* __restrict__ fieldGradientS, unsigned long long* __restrict__ fieldGradientPolarS
    #endif
#endif
        ) {
    atomicAdd(&field[index], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
    atomicAdd(&field[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
    atomicAdd(&field[index+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
    atomicAdd(&fieldPolar[index], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
    atomicAdd(&fieldPolar[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
    atomicAdd(&fieldPolar[index+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
    atomicAdd(&fieldS[index], static_cast<unsigned long long>((long long) (data.fieldS.x*0x100000000)));
    atomicAdd(&fieldS[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.y*0x100000000)));
    atomicAdd(&fieldS[index+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldS.z*0x100000000)));
    atomicAdd(&fieldPolarS[index], static_cast<unsigned long long>((long long) (data.fieldPolarS.x*0x100000000)));
    atomicAdd(&fieldPolarS[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.y*0x100000000)));
    atomicAdd(&fieldPolarS[index+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolarS.z*0x100000000)));
#endif
#ifdef EXTRAPOLATED_POLARIZATION
    for (int i = 0; i < 6; i++) {
        atomicAdd(&fieldGradient[6*index+i], static_cast<unsigned long long>((long long) (data.fieldGradient[i]*0x100000000)));
        atomicAdd(&fieldGradientPolar[6*index+i], static_cast<unsigned long long>((long long) (data.fieldGradientPolar[i]*0x100000000)));
#ifdef USE_GK
        atomicAdd(&fieldGradientS[6*index+i], static_cast<unsigned long long>((long long) (data.fieldGradientS[i]*0x100000000)));
        atomicAdd(&fieldGradientPolarS[6*index+i], static_cast<unsigned long long>((long long) (data.fieldGradientPolarS[i]*0x100000000)));
#endif
    }
#endif
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    real scale1, scale2, scale3;
    real r2 = dot(deltaR, deltaR);
    if (r2 < CUTOFF_SQUARED) {
        real rI = RSQRT(r2);
        real r = RECIP(rI);
        real rI2 = rI*rI;

        // calculate the error function damping terms

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
        real bn0 = erfcAlphaR*rI;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)*rI2;

        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)*rI2;

        alsq2n *= alsq2;
        real bn3 = (5*bn2+alsq2n*exp2a)*rI2;

        // compute the error function scaled and unscaled terms

        real damp = atom1.damp*atom2.damp;
        real ratio = (r/damp);
        ratio = ratio*ratio*ratio;
        float pgamma = atom1.thole < atom2.thole ? atom1.thole : atom2.thole;
        damp = damp == 0 ? 0 : -pgamma*ratio;
        real expdamp = EXP(damp);
        real dsc3 = 1 - expdamp;
        real dsc5 = 1 - expdamp*(1-damp);
        real dsc7 = 1 - (1-damp+(0.6f*damp*damp))*expdamp;
        real r3 = (r*r2);
        real r5 = (r3*r2);
        real r7 = (r5*r2);
        real rr3 = (1-dsc3)/r3;
        real rr5 = 3*(1-dsc5)/r5;
        real rr7 = 15*(1-dsc7)/r7;

        scale1 = rr3 - bn1;
        scale2 = bn2 - rr5;
        scale3 = bn3 - rr7;
    }
    else {
        scale1 = 0;
        scale2 = 0;
        scale3 = 0;
    }
    real dDotDelta = scale2*dot(deltaR, atom2.inducedDipole);
    atom1.field += scale1*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += scale1*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipole);
    atom2.field += scale1*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += scale1*atom1.inducedDipolePolar + dDotDelta*deltaR;
#ifdef EXTRAPOLATED_POLARIZATION
    // Compute and store the field gradients for later use.
    
    real3 dipole = atom1.inducedDipole;
    real muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2.fieldGradient[0] -= (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom2.fieldGradient[1] -= (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom2.fieldGradient[2] -= (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom2.fieldGradient[3] -= (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom2.fieldGradient[4] -= (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom2.fieldGradient[5] -= (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom1.inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2.fieldGradientPolar[0] -= (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom2.fieldGradientPolar[1] -= (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom2.fieldGradientPolar[2] -= (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom2.fieldGradientPolar[3] -= (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom2.fieldGradientPolar[4] -= (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom2.fieldGradientPolar[5] -= (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom2.inducedDipole;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1.fieldGradient[0] += (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom1.fieldGradient[1] += (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom1.fieldGradient[2] += (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom1.fieldGradient[3] += (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom1.fieldGradient[4] += (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom1.fieldGradient[5] += (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom2.inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1.fieldGradientPolar[0] += (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom1.fieldGradientPolar[1] += (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom1.fieldGradientPolar[2] += (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom1.fieldGradientPolar[3] += (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom1.fieldGradientPolar[4] += (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom1.fieldGradientPolar[5] += (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;
#endif
}
#elif defined USE_GK
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    real r2 = dot(deltaR, deltaR);
    real r = SQRT(r2);
    if (!isSelfInteraction) {
        real rI = RECIP(r);
        real r2I = rI*rI;
        real rr3 = -rI*r2I;
        real rr5 = -3*rr3*r2I;

        real dampProd = atom1.damp*atom2.damp;
        real ratio = (dampProd != 0 ? r/dampProd : 1);
        float pGamma = (atom1.thole > atom2.thole ? atom2.thole: atom1.thole);
        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 

        rr3 *= 1-dampExp;
        rr5 *= 1-(1+damp)*dampExp;

        real dDotDelta = rr5*dot(deltaR, atom2.inducedDipole);
        atom1.field += rr3*atom2.inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
        atom1.fieldPolar += rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
        atom2.field += rr3*atom1.inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
        atom2.fieldPolar += rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipoleS);
        atom1.fieldS += rr3*atom2.inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolarS);
        atom1.fieldPolarS += rr3*atom2.inducedDipolePolarS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipoleS);
        atom2.fieldS += rr3*atom1.inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolarS);
        atom2.fieldPolarS += rr3*atom1.inducedDipolePolarS + dDotDelta*deltaR;
    }

    real rb2 = atom1.bornRadius*atom2.bornRadius;
    real expterm = EXP(-r2/(GK_C*rb2));
    real expc = expterm/GK_C; 
    real gf2 = RECIP(r2+rb2*expterm);
    real gf = SQRT(gf2);
    real gf3 = gf2*gf;
    real gf5 = gf3*gf2;
    real a10 = -gf3;
    real expc1 = 1 - expc;
    real a11 = expc1 * 3 * gf5;
    real3 gux = GK_FD*make_real3(a10+deltaR.x*deltaR.x*a11, deltaR.x*deltaR.y*a11, deltaR.x*deltaR.z*a11);
    real3 guy = make_real3(gux.y, GK_FD*(a10+deltaR.y*deltaR.y*a11), GK_FD*deltaR.y*deltaR.z*a11);
    real3 guz = make_real3(gux.z, guy.z, GK_FD*(a10+deltaR.z*deltaR.z*a11));
 
    atom1.fieldS += atom2.inducedDipoleS.x*gux+atom2.inducedDipoleS.y*guy+atom2.inducedDipoleS.z*guz;
    atom2.fieldS += atom1.inducedDipoleS.x*gux+atom1.inducedDipoleS.y*guy+atom1.inducedDipoleS.z*guz;
    atom1.fieldPolarS += atom2.inducedDipolePolarS.x*gux+atom2.inducedDipolePolarS.y*guy+atom2.inducedDipolePolarS.z*guz;
    atom2.fieldPolarS += atom1.inducedDipolePolarS.x*gux+atom1.inducedDipolePolarS.y*guy+atom1.inducedDipolePolarS.z*guz;
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;
    real rr3 = -rI*r2I;
    real rr5 = -3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
    real dampProd = atom1.damp*atom2.damp;
    real ratio = (dampProd != 0 ? r/dampProd : 1);
    float pGamma = (atom2.thole > atom1.thole ? atom1.thole: atom2.thole);
    real damp = ratio*ratio*ratio*pGamma;
    real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    rr7 *= 1 - (1+damp+(0.6f*damp*damp))*dampExp;
    real dDotDelta = rr5*dot(deltaR, atom2.inducedDipole);
    atom1.field += rr3*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
    atom2.field += rr3*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
#ifdef EXTRAPOLATED_POLARIZATION
    // Compute and store the field gradients for later use.
    
    real3 dipole = atom1.inducedDipole;
    real muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2.fieldGradient[0] -= (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom2.fieldGradient[1] -= (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom2.fieldGradient[2] -= (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom2.fieldGradient[3] -= (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom2.fieldGradient[4] -= (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom2.fieldGradient[5] -= (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom1.inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2.fieldGradientPolar[0] -= (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom2.fieldGradientPolar[1] -= (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom2.fieldGradientPolar[2] -= (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom2.fieldGradientPolar[3] -= (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom2.fieldGradientPolar[4] -= (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom2.fieldGradientPolar[5] -= (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom2.inducedDipole;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1.fieldGradient[0] += (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom1.fieldGradient[1] += (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom1.fieldGradient[2] += (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom1.fieldGradient[3] += (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom1.fieldGradient[4] += (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom1.fieldGradient[5] += (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom2.inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1.fieldGradientPolar[0] += (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom1.fieldGradientPolar[1] += (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom1.fieldGradientPolar[2] += (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom1.fieldGradientPolar[3] += (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom1.fieldGradientPolar[4] += (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom1.fieldGradientPolar[5] += (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;
#endif
}
#endif

/**
 * Compute the mutual induced field.
 */
extern "C" __global__ void computeInducedField(
        unsigned long long* __restrict__ field, unsigned long long* __restrict__ fieldPolar, const real4* __restrict__ posq, const ushort2* __restrict__ exclusionTiles, 
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef EXTRAPOLATED_POLARIZATION
        unsigned long long* __restrict__ fieldGradient, unsigned long long* __restrict__ fieldGradientPolar,
#endif
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
#elif defined USE_GK
        unsigned long long* __restrict__ fieldS, unsigned long long* __restrict__ fieldPolarS, const real* __restrict__ inducedDipoleS,
        const real* __restrict__ inducedDipolePolarS, const real* __restrict__ bornRadii,
    #ifdef EXTRAPOLATED_POLARIZATION
        unsigned long long* __restrict__ fieldGradientS, unsigned long long* __restrict__ fieldGradientPolarS,
    #endif
#endif
        const float2* __restrict__ dampingAndThole) {
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
        zeroAtomData(data);
        unsigned int atom1 = x*TILE_SIZE + tgx;
#ifdef USE_GK
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].thole = data.thole;
            localData[threadIdx.x].damp = data.damp;
#ifdef USE_GK
            localData[threadIdx.x].inducedDipoleS = data.inducedDipoleS;
            localData[threadIdx.x].inducedDipolePolarS = data.inducedDipolePolarS;
            localData[threadIdx.x].bornRadius = data.bornRadius;
#endif
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+j].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+j], delta, atom1 == atom2);
            }
        }
        else {
            // This is an off-diagonal tile.

#ifdef USE_GK
            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomData(localData[threadIdx.x]);
            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
        }

        // Write results.

        unsigned int offset = x*TILE_SIZE + tgx;
        SAVE_ATOM_DATA(offset, data)
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            SAVE_ATOM_DATA(offset, localData[threadIdx.x])
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
            zeroAtomData(data);
#ifdef USE_GK
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
#ifdef USE_GK
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomData(localData[threadIdx.x]);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            SAVE_ATOM_DATA(offset, data)
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            SAVE_ATOM_DATA(offset, localData[threadIdx.x])
        }
        pos++;
    }
}

extern "C" __global__ void updateInducedFieldBySOR(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ fixedFieldS, const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar,
        real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors) {
    extern __shared__ real2 buffer[];
    const float polarSOR = 0.55f;
#ifdef USE_EWALD
    const real ewaldScale = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI;
#else
    const real ewaldScale = 0;
#endif
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            real previousDipole = inducedDipole[dipoleIndex];
            real previousDipolePolar = inducedDipolePolar[dipoleIndex];
            long long fixedS = (fixedFieldS == NULL ? (long long) 0 : fixedFieldS[fieldIndex]);
            real newDipole = scale*((fixedField[fieldIndex]+fixedS+inducedField[fieldIndex])*fieldScale+ewaldScale*previousDipole);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+fixedS+inducedFieldPolar[fieldIndex])*fieldScale+ewaldScale*previousDipolePolar);
            newDipole = previousDipole + polarSOR*(newDipole-previousDipole);
            newDipolePolar = previousDipolePolar + polarSOR*(newDipolePolar-previousDipolePolar);
            inducedDipole[dipoleIndex] = newDipole;
            inducedDipolePolar[dipoleIndex] = newDipolePolar;
            sumErrors += (newDipole-previousDipole)*(newDipole-previousDipole);
            sumPolarErrors += (newDipolePolar-previousDipolePolar)*(newDipolePolar-previousDipolePolar);
        }
    }
    
    // Sum the errors over threads and store the total for this block.
    
    buffer[threadIdx.x] = make_real2(sumErrors, sumPolarErrors);
    __syncthreads();
    for (int offset = 1; offset < blockDim.x; offset *= 2) {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0) {
            buffer[threadIdx.x].x += buffer[threadIdx.x+offset].x;
            buffer[threadIdx.x].y += buffer[threadIdx.x+offset].y;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        errors[blockIdx.x] = make_float2((float) buffer[0].x, (float) buffer[0].y);
}

extern "C" __global__ void recordInducedDipolesForDIIS(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ fixedFieldS, const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar,
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors,
        real* __restrict__ prevDipoles, real* __restrict__ prevDipolesPolar, real* __restrict__ prevErrors, int iteration, bool recordPrevErrors, real* __restrict__ matrix) {
    extern __shared__ real2 buffer[];
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            if (iteration >= MAX_PREV_DIIS_DIPOLES) {
                // We have filled up the buffer for previous dipoles, so shift them all over by one.
                
                for (int i = 1; i < MAX_PREV_DIIS_DIPOLES; i++) {
                    int index1 = dipoleIndex+(i-1)*NUM_ATOMS*3;
                    int index2 = dipoleIndex+i*NUM_ATOMS*3;
                    prevDipoles[index1] = prevDipoles[index2];
                    prevDipolesPolar[index1] = prevDipolesPolar[index2];
                    if (recordPrevErrors)
                        prevErrors[index1] = prevErrors[index2];
                }
            }
            
            // Compute the new dipole, and record it along with the error.
            
            real oldDipole = inducedDipole[dipoleIndex];
            real oldDipolePolar = inducedDipolePolar[dipoleIndex];
            long long fixedS = (fixedFieldS == NULL ? (long long) 0 : fixedFieldS[fieldIndex]);
            real newDipole = scale*((fixedField[fieldIndex]+fixedS+inducedField[fieldIndex])*fieldScale);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+fixedS+inducedFieldPolar[fieldIndex])*fieldScale);
            int storePrevIndex = dipoleIndex+min(iteration, MAX_PREV_DIIS_DIPOLES-1)*NUM_ATOMS*3;
            prevDipoles[storePrevIndex] = newDipole;
            prevDipolesPolar[storePrevIndex] = newDipolePolar;
            if (recordPrevErrors)
                prevErrors[storePrevIndex] = newDipole-oldDipole;
            sumErrors += (newDipole-oldDipole)*(newDipole-oldDipole);
            sumPolarErrors += (newDipolePolar-oldDipolePolar)*(newDipolePolar-oldDipolePolar);
        }
    }
    
    // Sum the errors over threads and store the total for this block.
    
    buffer[threadIdx.x] = make_real2(sumErrors, sumPolarErrors);
    __syncthreads();
    for (int offset = 1; offset < blockDim.x; offset *= 2) {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0) {
            buffer[threadIdx.x].x += buffer[threadIdx.x+offset].x;
            buffer[threadIdx.x].y += buffer[threadIdx.x+offset].y;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        errors[blockIdx.x] = make_float2((float) buffer[0].x, (float) buffer[0].y);
    
    if (iteration >= MAX_PREV_DIIS_DIPOLES && recordPrevErrors && blockIdx.x == 0) {
        // Shift over the existing matrix elements.
        
        for (int i = 0; i < MAX_PREV_DIIS_DIPOLES-1; i++) {
            if (threadIdx.x < MAX_PREV_DIIS_DIPOLES-1)
                matrix[threadIdx.x+i*MAX_PREV_DIIS_DIPOLES] = matrix[(threadIdx.x+1)+(i+1)*MAX_PREV_DIIS_DIPOLES];
            __syncthreads();
        }
    }
}

extern "C" __global__ void computeDIISMatrix(real* __restrict__ prevErrors, int iteration, real* __restrict__ matrix) {
    extern __shared__ real sumBuffer[];
    int j = min(iteration, MAX_PREV_DIIS_DIPOLES-1);
    for (int i = blockIdx.x; i <= j; i += gridDim.x) {
        // All the threads in this thread block work together to compute a single matrix element.

        real sum = 0;
        for (int index = threadIdx.x; index < NUM_ATOMS*3; index += blockDim.x)
            sum += prevErrors[index+i*NUM_ATOMS*3]*prevErrors[index+j*NUM_ATOMS*3];
        sumBuffer[threadIdx.x] = sum;
        __syncthreads();
        for (int offset = 1; offset < blockDim.x; offset *= 2) { 
            if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
                sumBuffer[threadIdx.x] += sumBuffer[threadIdx.x+offset];
            __syncthreads();
        }
        if (threadIdx.x == 0) {
            matrix[i+MAX_PREV_DIIS_DIPOLES*j] = sumBuffer[0];
            if (i != j)
                matrix[j+MAX_PREV_DIIS_DIPOLES*i] = sumBuffer[0];
        }
    }
}

extern "C" __global__ void solveDIISMatrix(int iteration, const real* __restrict__ matrix, float* __restrict__ coefficients) {
    __shared__ real b[MAX_PREV_DIIS_DIPOLES+1][MAX_PREV_DIIS_DIPOLES+1];
    __shared__ real piv[MAX_PREV_DIIS_DIPOLES+1];
    __shared__ real x[MAX_PREV_DIIS_DIPOLES+1];

    // On the first iteration we don't need to do any calculation.
    
    if (iteration == 0) {
        if (threadIdx.x == 0)
            coefficients[0] = 1;
        return;
    }
    
    // Load the matrix.
    
    int numPrev = min(iteration+1, MAX_PREV_DIIS_DIPOLES);
    int rank = numPrev+1;
    for (int index = threadIdx.x; index < numPrev*numPrev; index += blockDim.x) {
        int i = index/numPrev;
        int j = index-i*numPrev;
        b[i+1][j+1] = matrix[i*MAX_PREV_DIIS_DIPOLES+j];
    }
    for (int i = threadIdx.x; i < rank; i += blockDim.x) {
        b[i][0] = -1;
        piv[i] = i;
    }
    __syncthreads();
    
    // Compute the mean absolute value of the values we just loaded.  We use that for preconditioning it,
    // which is essential for doing the computation in single precision.
    
    if (threadIdx.x == 0) {
        real mean = 0;
        for (int i = 0; i < numPrev; i++)
            for (int j = 0; j < numPrev; j++)
                mean += fabs(b[i+1][j+1]);
        mean /= numPrev*numPrev;
        b[0][0] = 0;
        for (int i = 1; i < rank; i++)
            b[0][i] = -mean;

        // Compute the LU decomposition of the matrix.  This code is adapted from JAMA.
    
        int pivsign = 1;
        for (int j = 0; j < rank; j++) {
            // Apply previous transformations.

            for (int i = 0; i < rank; i++) {
                // Most of the time is spent in the following dot product.

                int kmax = min(i, j);
                real s = 0;
                for (int k = 0; k < kmax; k++)
                    s += b[i][k] * b[k][j];
                b[i][j] -= s;
            }

            // Find pivot and exchange if necessary.

            int p = j;
            for (int i = j+1; i < rank; i++)
                if (abs(b[i][j]) > abs(b[p][j]))
                    p = i;
            if (p != j) {
                int k = 0;
                for (k = 0; k < rank; k++) {
                    real t = b[p][k];
                    b[p][k] = b[j][k];
                    b[j][k] = t;
                }
                k = piv[p];
                piv[p] = piv[j];
                piv[j] = k;
                pivsign = -pivsign;
            }

            // Compute multipliers.

            if ((j < rank) && (b[j][j] != 0))
                for (int i = j+1; i < rank; i++)
                    b[i][j] /= b[j][j];
        }
        for (int i = 0; i < rank; i++)
            if (b[i][i] == 0) {
                // The matrix is singular.
                
                for (int j = 0; j < rank-1; j++)
                    coefficients[j] = 0;
                coefficients[rank-1] = 1;
                return;
            }

        // Solve b*Y = X(piv)
        
        for (int i = 0; i < rank; i++) 
            x[i] = (piv[i] == 0 ? -1 : 0);
        for (int k = 0; k < rank; k++)
            for (int i = k+1; i < rank; i++)
                x[i] -= x[k] * b[i][k];

        // Solve U*X = Y;
        
        for (int k = rank-1; k >= 0; k--) {
            x[k] /= b[k][k];
            for (int i = 0; i < k; i++)
                x[i] -= x[k] * b[i][k];
        }
        
        // Record the coefficients.
        
        real lastCoeff = 1;
        for (int i = 0; i < rank-1; i++) {
            real c = x[i+1]*mean;
            coefficients[i] = c;
            lastCoeff -= c;
        }
        coefficients[rank-1] = lastCoeff;
    }
}

extern "C" __global__ void updateInducedFieldByDIIS(real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, 
        const real* __restrict__ prevDipoles, const real* __restrict__ prevDipolesPolar, const float* __restrict__ coefficients, int numPrev) {
    for (int index = blockIdx.x*blockDim.x + threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real sum = 0;
        real sumPolar = 0;
        for (int i = 0; i < numPrev; i++) {
            sum += coefficients[i]*prevDipoles[i*3*NUM_ATOMS+index];
            sumPolar += coefficients[i]*prevDipolesPolar[i*3*NUM_ATOMS+index];
        }
        inducedDipole[index] = sum;
        inducedDipolePolar[index] = sumPolar;
    }
}
#endif // not HIPPO

extern "C" __global__ void initExtrapolatedDipoles(real* __restrict__ inducedDipole, real* __restrict__ extrapolatedDipole
#ifndef HIPPO
        , real* __restrict__ inducedDipolePolar, real* __restrict__ extrapolatedDipolePolar, long long* __restrict__ inducedDipoleFieldGradient, long long* __restrict__ inducedDipoleFieldGradientPolar
#endif
#ifdef USE_GK
        , real* __restrict__ inducedDipoleGk, real* __restrict__ inducedDipoleGkPolar, real* __restrict__ extrapolatedDipoleGk, real* __restrict__ extrapolatedDipoleGkPolar,
        real* __restrict__ inducedDipoleFieldGradientGk, real* __restrict__ inducedDipoleFieldGradientGkPolar
#endif
        ) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        extrapolatedDipole[index] = inducedDipole[index];
#ifndef HIPPO
        extrapolatedDipolePolar[index] = inducedDipolePolar[index];
#endif
#ifdef USE_GK
        extrapolatedDipoleGk[index] = inducedDipoleGk[index];
        extrapolatedDipoleGkPolar[index] = inducedDipoleGkPolar[index];
#endif
    }
#ifndef HIPPO
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < 6*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        inducedDipoleFieldGradient[index] = 0;
        inducedDipoleFieldGradientPolar[index] = 0;
#ifdef USE_GK
        inducedDipoleFieldGradientGk[index] = 0;
        inducedDipoleFieldGradientGkPolar[index] = 0;
#endif
    }
#endif
}

extern "C" __global__ void iterateExtrapolatedDipoles(int order, real* __restrict__ inducedDipole, real* __restrict__ extrapolatedDipole, long long* __restrict__ inducedDipoleField,
#ifndef HIPPO
        real* __restrict__ inducedDipolePolar, real* __restrict__ extrapolatedDipolePolar, long long* __restrict__ inducedDipoleFieldPolar, long long* __restrict__ inducedDipoleFieldGradient,
        long long* __restrict__ inducedDipoleFieldGradientPolar, real* __restrict__ extrapolatedDipoleFieldGradient, real* __restrict__ extrapolatedDipoleFieldGradientPolar,
#endif
#ifdef USE_GK
        real* __restrict__ inducedDipoleGk, real* __restrict__ inducedDipoleGkPolar, real* __restrict__ extrapolatedDipoleGk, real* __restrict__ extrapolatedDipoleGkPolar,
        real* __restrict__ inducedDipoleFieldGradientGk, real* __restrict__ inducedDipoleFieldGradientGkPolar, long long* __restrict__ inducedDipoleFieldGk,
        long long* __restrict__ inducedDipoleFieldGkPolar, real* __restrict__ extrapolatedDipoleFieldGradientGk, real* __restrict__ extrapolatedDipoleFieldGradientGkPolar,
#endif
#ifdef HIPPO
        const real* __restrict__ polarizability
#else
        const float* __restrict__ polarizability
#endif
    ) {
    const real fieldScale = 1/(real) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        int atom = index/3;
        int component = index-3*atom;
        int fieldIndex = atom+component*PADDED_NUM_ATOMS;
        float polar = polarizability[atom];
        real value = inducedDipoleField[fieldIndex]*fieldScale*polar;
        inducedDipole[index] = value;
        extrapolatedDipole[order*3*NUM_ATOMS+index] = value;
#ifndef HIPPO
        value = inducedDipoleFieldPolar[fieldIndex]*fieldScale*polar;
        inducedDipolePolar[index] = value;
        extrapolatedDipolePolar[order*3*NUM_ATOMS+index] = value;
#endif
#ifdef USE_GK
        value = inducedDipoleFieldGk[fieldIndex]*fieldScale*polar;
        inducedDipoleGk[index] = value;
        extrapolatedDipoleGk[order*3*NUM_ATOMS+index] = value;
        value = inducedDipoleFieldGkPolar[fieldIndex]*fieldScale*polar;
        inducedDipoleGkPolar[index] = value;
        extrapolatedDipoleGkPolar[order*3*NUM_ATOMS+index] = value;
#endif
    }
#ifndef HIPPO
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < 6*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        int index2 = (order-1)*6*NUM_ATOMS+index;
        extrapolatedDipoleFieldGradient[index2] = fieldScale*inducedDipoleFieldGradient[index];
        extrapolatedDipoleFieldGradientPolar[index2] = fieldScale*inducedDipoleFieldGradientPolar[index];
#ifdef USE_GK
        extrapolatedDipoleFieldGradientGk[index2] = fieldScale*inducedDipoleFieldGradientGk[index];
        extrapolatedDipoleFieldGradientGkPolar[index2] = fieldScale*inducedDipoleFieldGradientGkPolar[index];
#endif
    }
#endif
}

extern "C" __global__ void computeExtrapolatedDipoles(real* __restrict__ inducedDipole, real* __restrict__ extrapolatedDipole
#ifndef HIPPO
        , real* __restrict__ inducedDipolePolar, real* __restrict__ extrapolatedDipolePolar
#endif
#ifdef USE_GK
        , real* __restrict__ inducedDipoleGk, real* __restrict__ inducedDipoleGkPolar, real* __restrict__ extrapolatedDipoleGk, real* __restrict__ extrapolatedDipoleGkPolar
#endif
        ) {
    real coeff[] = {EXTRAPOLATION_COEFFICIENTS_SUM};
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real sum = 0, sumPolar = 0, sumGk = 0, sumGkPolar = 0;
        for (int order = 0; order < MAX_EXTRAPOLATION_ORDER; order++) {
            sum += extrapolatedDipole[order*3*NUM_ATOMS+index]*coeff[order];
#ifndef HIPPO
            sumPolar += extrapolatedDipolePolar[order*3*NUM_ATOMS+index]*coeff[order];
#endif
#ifdef USE_GK
            sumGk += extrapolatedDipoleGk[order*3*NUM_ATOMS+index]*coeff[order];
            sumGkPolar += extrapolatedDipoleGkPolar[order*3*NUM_ATOMS+index]*coeff[order];
#endif
        }
        inducedDipole[index] = sum;
#ifndef HIPPO
        inducedDipolePolar[index] = sumPolar;
#endif
#ifdef USE_GK
        inducedDipoleGk[index] = sumGk;
        inducedDipoleGkPolar[index] = sumGkPolar;
#endif
    }
}

extern "C" __global__ void addExtrapolatedFieldGradientToForce(long long* __restrict__ forceBuffers, real* __restrict__ extrapolatedDipole,
        real* __restrict__ extrapolatedDipolePolar, real* __restrict__ extrapolatedDipoleFieldGradient, real* __restrict__ extrapolatedDipoleFieldGradientPolar
#ifdef USE_GK
        , real* __restrict__ extrapolatedDipoleGk, real* __restrict__ extrapolatedDipoleGkPolar,
        real* __restrict__ extrapolatedDipoleFieldGradientGk, real* __restrict__ extrapolatedDipoleFieldGradientGkPolar
#endif
        ) {
    real coeff[] = {EXTRAPOLATION_COEFFICIENTS_SUM};
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real fx = 0, fy = 0, fz = 0;
        for (int l = 0; l < MAX_EXTRAPOLATION_ORDER-1; l++) {
            int index1 = 3*(l*NUM_ATOMS+atom);
            real dipole[] = {extrapolatedDipole[index1], extrapolatedDipole[index1+1], extrapolatedDipole[index1+2]};
            real dipolePolar[] = {extrapolatedDipolePolar[index1], extrapolatedDipolePolar[index1+1], extrapolatedDipolePolar[index1+2]};
#ifdef USE_GK
            real dipoleGk[] = {extrapolatedDipoleGk[index1], extrapolatedDipoleGk[index1+1], extrapolatedDipoleGk[index1+2]};
            real dipoleGkPolar[] = {extrapolatedDipoleGkPolar[index1], extrapolatedDipoleGkPolar[index1+1], extrapolatedDipoleGkPolar[index1+2]};
#endif
            for (int m = 0; m < MAX_EXTRAPOLATION_ORDER-1-l; m++) {
                int index2 = 6*(m*NUM_ATOMS+atom);
                real scale = 0.5f*coeff[l+m+1]*ENERGY_SCALE_FACTOR;
                real gradient[] = {extrapolatedDipoleFieldGradient[index2], extrapolatedDipoleFieldGradient[index2+1], extrapolatedDipoleFieldGradient[index2+2],
                                   extrapolatedDipoleFieldGradient[index2+3], extrapolatedDipoleFieldGradient[index2+4], extrapolatedDipoleFieldGradient[index2+5]};
                real gradientPolar[] = {extrapolatedDipoleFieldGradientPolar[index2], extrapolatedDipoleFieldGradientPolar[index2+1], extrapolatedDipoleFieldGradientPolar[index2+2],
                                        extrapolatedDipoleFieldGradientPolar[index2+3], extrapolatedDipoleFieldGradientPolar[index2+4], extrapolatedDipoleFieldGradientPolar[index2+5]};
                fx += scale*(dipole[0]*gradientPolar[0] + dipole[1]*gradientPolar[3] + dipole[2]*gradientPolar[4]);
                fy += scale*(dipole[0]*gradientPolar[3] + dipole[1]*gradientPolar[1] + dipole[2]*gradientPolar[5]);
                fz += scale*(dipole[0]*gradientPolar[4] + dipole[1]*gradientPolar[5] + dipole[2]*gradientPolar[2]);
                fx += scale*(dipolePolar[0]*gradient[0] + dipolePolar[1]*gradient[3] + dipolePolar[2]*gradient[4]);
                fy += scale*(dipolePolar[0]*gradient[3] + dipolePolar[1]*gradient[1] + dipolePolar[2]*gradient[5]);
                fz += scale*(dipolePolar[0]*gradient[4] + dipolePolar[1]*gradient[5] + dipolePolar[2]*gradient[2]);
#ifdef USE_GK
                real gradientGk[] = {extrapolatedDipoleFieldGradient[index2], extrapolatedDipoleFieldGradient[index2+1], extrapolatedDipoleFieldGradient[index2+2],
                                   extrapolatedDipoleFieldGradient[index2+3], extrapolatedDipoleFieldGradient[index2+4], extrapolatedDipoleFieldGradient[index2+5]};
                real gradientGkPolar[] = {extrapolatedDipoleFieldGradientPolar[index2], extrapolatedDipoleFieldGradientPolar[index2+1], extrapolatedDipoleFieldGradientPolar[index2+2],
                                        extrapolatedDipoleFieldGradientPolar[index2+3], extrapolatedDipoleFieldGradientPolar[index2+4], extrapolatedDipoleFieldGradientPolar[index2+5]};
                fx += scale*(dipoleGk[0]*gradientGkPolar[0] + dipoleGk[1]*gradientGkPolar[3] + dipoleGk[2]*gradientGkPolar[4]);
                fy += scale*(dipoleGk[0]*gradientGkPolar[3] + dipoleGk[1]*gradientGkPolar[1] + dipoleGk[2]*gradientGkPolar[5]);
                fz += scale*(dipoleGk[0]*gradientGkPolar[4] + dipoleGk[1]*gradientGkPolar[5] + dipoleGk[2]*gradientGkPolar[2]);
                fx += scale*(dipoleGkPolar[0]*gradientGk[0] + dipoleGkPolar[1]*gradientGk[3] + dipoleGkPolar[2]*gradientGk[4]);
                fy += scale*(dipoleGkPolar[0]*gradientGk[3] + dipoleGkPolar[1]*gradientGk[1] + dipoleGkPolar[2]*gradientGk[5]);
                fz += scale*(dipoleGkPolar[0]*gradientGk[4] + dipoleGkPolar[1]*gradientGk[5] + dipoleGkPolar[2]*gradientGk[2]);
#endif
            }
        }
        forceBuffers[atom] += (long long) (fx*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (long long) (fy*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS*2] += (long long) (fz*0x100000000);
    }
}

#ifdef HIPPO
extern "C" __global__ void computePolarizationEnergy(mixed* __restrict__ energyBuffer, const real3* __restrict__ inducedDipole,
        const real3* __restrict__ extrapolatedDipole, const real* __restrict__ polarizability) {
    mixed energy = 0;
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x)
        energy -= (ENERGY_SCALE_FACTOR/2)*dot(extrapolatedDipole[atom], inducedDipole[atom])/polarizability[atom];
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}
#endif