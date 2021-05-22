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
inline DEVICE AtomData loadAtomData(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT inducedDipole,
        GLOBAL const real* RESTRICT inducedDipolePolar, GLOBAL const float2* RESTRICT dampingAndThole, GLOBAL const real* RESTRICT inducedDipoleS,
        GLOBAL const real* RESTRICT inducedDipolePolarS, GLOBAL const real* RESTRICT bornRadii) {
#else
inline DEVICE AtomData loadAtomData(int atom, GLOBAL const real4* RESTRICT posq, GLOBAL const real* RESTRICT inducedDipole,
        GLOBAL const real* RESTRICT inducedDipolePolar, GLOBAL const float2* RESTRICT dampingAndThole) {
#endif
    AtomData data;
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.inducedDipole = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
    data.inducedDipolePolar = make_real3(inducedDipolePolar[3*atom], inducedDipolePolar[3*atom+1], inducedDipolePolar[3*atom+2]);
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
#ifdef USE_GK
    data.inducedDipoleS = make_real3(inducedDipoleS[3*atom], inducedDipoleS[3*atom+1], inducedDipoleS[3*atom+2]);
    data.inducedDipolePolarS = make_real3(inducedDipolePolarS[3*atom], inducedDipolePolarS[3*atom+1], inducedDipolePolarS[3*atom+2]);
    data.bornRadius = bornRadii[atom];
#endif
    return data;
}

inline DEVICE void zeroAtomData(AtomData* data) {
    data->field = make_real3(0);
    data->fieldPolar = make_real3(0);
#ifdef USE_GK
    data->fieldS = make_real3(0);
    data->fieldPolarS = make_real3(0);
#endif
#ifdef EXTRAPOLATED_POLARIZATION
    for (int i = 0; i < 6; i++) {
        data->fieldGradient[i] = 0;
        data->fieldGradientPolar[i] = 0;
#ifdef USE_GK
        data->fieldGradientS[i] = 0;
        data->fieldGradientPolarS[i] = 0;
#endif
    }
#endif
}

// OpenCL requires a second version of this function, since the signature depends
// on the address space of the argument.

inline DEVICE void zeroAtomDataLocal(LOCAL_ARG AtomData* data) {
    data->field = make_real3(0);
    data->fieldPolar = make_real3(0);
#ifdef USE_GK
    data->fieldS = make_real3(0);
    data->fieldPolarS = make_real3(0);
#endif
#ifdef EXTRAPOLATED_POLARIZATION
    for (int i = 0; i < 6; i++) {
        data->fieldGradient[i] = 0;
        data->fieldGradientPolar[i] = 0;
#ifdef USE_GK
        data->fieldGradientS[i] = 0;
        data->fieldGradientPolarS[i] = 0;
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

inline DEVICE void saveAtomData(int index, AtomData data, GLOBAL mm_ulong* RESTRICT field, GLOBAL mm_ulong* RESTRICT fieldPolar
#ifdef EXTRAPOLATED_POLARIZATION
        , GLOBAL mm_ulong* RESTRICT fieldGradient, GLOBAL mm_ulong* RESTRICT fieldGradientPolar
#endif
#ifdef USE_GK
        , GLOBAL mm_ulong* RESTRICT fieldS, GLOBAL mm_ulong* RESTRICT fieldPolarS
    #ifdef EXTRAPOLATED_POLARIZATION
        , GLOBAL mm_ulong* RESTRICT fieldGradientS, GLOBAL mm_ulong* RESTRICT fieldGradientPolarS
    #endif
#endif
        ) {
    ATOMIC_ADD(&field[index], (mm_ulong) ((mm_long) (data.field.x*0x100000000)));
    ATOMIC_ADD(&field[index+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.field.y*0x100000000)));
    ATOMIC_ADD(&field[index+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.field.z*0x100000000)));
    ATOMIC_ADD(&fieldPolar[index], (mm_ulong) ((mm_long) (data.fieldPolar.x*0x100000000)));
    ATOMIC_ADD(&fieldPolar[index+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldPolar.y*0x100000000)));
    ATOMIC_ADD(&fieldPolar[index+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_GK
    ATOMIC_ADD(&fieldS[index], (mm_ulong) ((mm_long) (data.fieldS.x*0x100000000)));
    ATOMIC_ADD(&fieldS[index+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldS.y*0x100000000)));
    ATOMIC_ADD(&fieldS[index+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldS.z*0x100000000)));
    ATOMIC_ADD(&fieldPolarS[index], (mm_ulong) ((mm_long) (data.fieldPolarS.x*0x100000000)));
    ATOMIC_ADD(&fieldPolarS[index+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldPolarS.y*0x100000000)));
    ATOMIC_ADD(&fieldPolarS[index+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (data.fieldPolarS.z*0x100000000)));
#endif
#ifdef EXTRAPOLATED_POLARIZATION
    for (int i = 0; i < 6; i++) {
        ATOMIC_ADD(&fieldGradient[6*index+i], (mm_ulong) ((mm_long) (data.fieldGradient[i]*0x100000000)));
        ATOMIC_ADD(&fieldGradientPolar[6*index+i], (mm_ulong) ((mm_long) (data.fieldGradientPolar[i]*0x100000000)));
#ifdef USE_GK
        ATOMIC_ADD(&fieldGradientS[6*index+i], (mm_ulong) ((mm_long) (data.fieldGradientS[i]*0x100000000)));
        ATOMIC_ADD(&fieldGradientPolarS[6*index+i], (mm_ulong) ((mm_long) (data.fieldGradientPolarS[i]*0x100000000)));
#endif
    }
#endif
}

#ifdef USE_EWALD
DEVICE void computeOneInteraction(AtomData* atom1, LOCAL_ARG AtomData* atom2, real3 deltaR, bool isSelfInteraction) {
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

        real damp = atom1->damp*atom2->damp;
        real ratio = (r/damp);
        ratio = ratio*ratio*ratio;
        float pgamma = atom1->thole < atom2->thole ? atom1->thole : atom2->thole;
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
    real dDotDelta = scale2*dot(deltaR, atom2->inducedDipole);
    atom1->field += scale1*atom2->inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom2->inducedDipolePolar);
    atom1->fieldPolar += scale1*atom2->inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1->inducedDipole);
    atom2->field += scale1*atom1->inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1->inducedDipolePolar);
    atom2->fieldPolar += scale1*atom1->inducedDipolePolar + dDotDelta*deltaR;
#ifdef EXTRAPOLATED_POLARIZATION
    // Compute and store the field gradients for later use.
    
    real3 dipole = atom1->inducedDipole;
    real muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2->fieldGradient[0] -= (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom2->fieldGradient[1] -= (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom2->fieldGradient[2] -= (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom2->fieldGradient[3] -= (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom2->fieldGradient[4] -= (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom2->fieldGradient[5] -= (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom1->inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2->fieldGradientPolar[0] -= (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom2->fieldGradientPolar[1] -= (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom2->fieldGradientPolar[2] -= (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom2->fieldGradientPolar[3] -= (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom2->fieldGradientPolar[4] -= (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom2->fieldGradientPolar[5] -= (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom2->inducedDipole;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1->fieldGradient[0] += (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom1->fieldGradient[1] += (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom1->fieldGradient[2] += (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom1->fieldGradient[3] += (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom1->fieldGradient[4] += (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom1->fieldGradient[5] += (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;

    dipole = atom2->inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1->fieldGradientPolar[0] += (muDotR*scale3)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*scale2;
    atom1->fieldGradientPolar[1] += (muDotR*scale3)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*scale2;
    atom1->fieldGradientPolar[2] += (muDotR*scale3)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*scale2;
    atom1->fieldGradientPolar[3] += (muDotR*scale3)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*scale2;
    atom1->fieldGradientPolar[4] += (muDotR*scale3)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*scale2;
    atom1->fieldGradientPolar[5] += (muDotR*scale3)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*scale2;
#endif
}
#elif defined USE_GK
DEVICE void computeOneInteraction(AtomData* atom1, LOCAL_ARG AtomData* atom2, real3 deltaR, bool isSelfInteraction) {
    real r2 = dot(deltaR, deltaR);
    real r = SQRT(r2);
    if (!isSelfInteraction) {
        real rI = RECIP(r);
        real r2I = rI*rI;
        real rr3 = -rI*r2I;
        real rr5 = -3*rr3*r2I;

        real dampProd = atom1->damp*atom2->damp;
        real ratio = (dampProd != 0 ? r/dampProd : 1);
        float pGamma = (atom1->thole > atom2->thole ? atom2->thole: atom1->thole);
        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 

        rr3 *= 1-dampExp;
        rr5 *= 1-(1+damp)*dampExp;

        real dDotDelta = rr5*dot(deltaR, atom2->inducedDipole);
        atom1->field += rr3*atom2->inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2->inducedDipolePolar);
        atom1->fieldPolar += rr3*atom2->inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1->inducedDipole);
        atom2->field += rr3*atom1->inducedDipole + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1->inducedDipolePolar);
        atom2->fieldPolar += rr3*atom1->inducedDipolePolar + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2->inducedDipoleS);
        atom1->fieldS += rr3*atom2->inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom2->inducedDipolePolarS);
        atom1->fieldPolarS += rr3*atom2->inducedDipolePolarS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1->inducedDipoleS);
        atom2->fieldS += rr3*atom1->inducedDipoleS + dDotDelta*deltaR;
        dDotDelta = rr5*dot(deltaR, atom1->inducedDipolePolarS);
        atom2->fieldPolarS += rr3*atom1->inducedDipolePolarS + dDotDelta*deltaR;
    }

    real rb2 = atom1->bornRadius*atom2->bornRadius;
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
 
    atom1->fieldS += atom2->inducedDipoleS.x*gux+atom2->inducedDipoleS.y*guy+atom2->inducedDipoleS.z*guz;
    atom2->fieldS += atom1->inducedDipoleS.x*gux+atom1->inducedDipoleS.y*guy+atom1->inducedDipoleS.z*guz;
    atom1->fieldPolarS += atom2->inducedDipolePolarS.x*gux+atom2->inducedDipolePolarS.y*guy+atom2->inducedDipolePolarS.z*guz;
    atom2->fieldPolarS += atom1->inducedDipolePolarS.x*gux+atom1->inducedDipolePolarS.y*guy+atom1->inducedDipolePolarS.z*guz;
}
#else
DEVICE void computeOneInteraction(AtomData* atom1, LOCAL_ARG AtomData* atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;
    real rr3 = -rI*r2I;
    real rr5 = -3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
    real dampProd = atom1->damp*atom2->damp;
    real ratio = (dampProd != 0 ? r/dampProd : 1);
    float pGamma = (atom2->thole > atom1->thole ? atom1->thole: atom2->thole);
    real damp = ratio*ratio*ratio*pGamma;
    real dampExp = (dampProd != 0 ? EXP(-damp) : 0); 
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    rr7 *= 1 - (1+damp+(0.6f*damp*damp))*dampExp;
    real dDotDelta = rr5*dot(deltaR, atom2->inducedDipole);
    atom1->field += rr3*atom2->inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom2->inducedDipolePolar);
    atom1->fieldPolar += rr3*atom2->inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1->inducedDipole);
    atom2->field += rr3*atom1->inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1->inducedDipolePolar);
    atom2->fieldPolar += rr3*atom1->inducedDipolePolar + dDotDelta*deltaR;
#ifdef EXTRAPOLATED_POLARIZATION
    // Compute and store the field gradients for later use.
    
    real3 dipole = atom1->inducedDipole;
    real muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2->fieldGradient[0] -= (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom2->fieldGradient[1] -= (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom2->fieldGradient[2] -= (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom2->fieldGradient[3] -= (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom2->fieldGradient[4] -= (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom2->fieldGradient[5] -= (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom1->inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom2->fieldGradientPolar[0] -= (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom2->fieldGradientPolar[1] -= (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom2->fieldGradientPolar[2] -= (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom2->fieldGradientPolar[3] -= (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom2->fieldGradientPolar[4] -= (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom2->fieldGradientPolar[5] -= (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom2->inducedDipole;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1->fieldGradient[0] += (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom1->fieldGradient[1] += (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom1->fieldGradient[2] += (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom1->fieldGradient[3] += (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom1->fieldGradient[4] += (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom1->fieldGradient[5] += (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;

    dipole = atom2->inducedDipolePolar;
    muDotR = dipole.x*deltaR.x + dipole.y*deltaR.y + dipole.z*deltaR.z;
    atom1->fieldGradientPolar[0] += (muDotR*rr7)*deltaR.x*deltaR.x - (2*dipole.x*deltaR.x + muDotR)*rr5;
    atom1->fieldGradientPolar[1] += (muDotR*rr7)*deltaR.y*deltaR.y - (2*dipole.y*deltaR.y + muDotR)*rr5;
    atom1->fieldGradientPolar[2] += (muDotR*rr7)*deltaR.z*deltaR.z - (2*dipole.z*deltaR.z + muDotR)*rr5;
    atom1->fieldGradientPolar[3] += (muDotR*rr7)*deltaR.x*deltaR.y - (dipole.x*deltaR.y + dipole.y*deltaR.x)*rr5;
    atom1->fieldGradientPolar[4] += (muDotR*rr7)*deltaR.x*deltaR.z - (dipole.x*deltaR.z + dipole.z*deltaR.x)*rr5;
    atom1->fieldGradientPolar[5] += (muDotR*rr7)*deltaR.y*deltaR.z - (dipole.y*deltaR.z + dipole.z*deltaR.y)*rr5;
#endif
}
#endif

/**
 * Compute the mutual induced field.
 */
KERNEL void computeInducedField(
        GLOBAL mm_ulong* RESTRICT field, GLOBAL mm_ulong* RESTRICT fieldPolar, GLOBAL const real4* RESTRICT posq, GLOBAL const int2* RESTRICT exclusionTiles, 
        GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter, GLOBAL const unsigned int* RESTRICT interactingAtoms,
#elif defined USE_GK
        GLOBAL mm_ulong* RESTRICT fieldS, GLOBAL mm_ulong* RESTRICT fieldPolarS, GLOBAL const real* RESTRICT inducedDipoleS,
        GLOBAL const real* RESTRICT inducedDipolePolarS, GLOBAL const real* RESTRICT bornRadii,
#endif
#ifdef EXTRAPOLATED_POLARIZATION
        GLOBAL mm_ulong* RESTRICT fieldGradient, GLOBAL mm_ulong* RESTRICT fieldGradientPolar,
    #ifdef USE_GK
        GLOBAL mm_ulong* RESTRICT fieldGradientS, GLOBAL mm_ulong* RESTRICT fieldGradientPolarS,
    #endif
#endif
        GLOBAL const float2* RESTRICT dampingAndThole) {
    const unsigned int totalWarps = (GLOBAL_SIZE)/TILE_SIZE;
    const unsigned int warp = (GLOBAL_ID)/TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE-1);
    const unsigned int tbx = LOCAL_ID - tgx;
    LOCAL AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        unsigned int atom1 = x*TILE_SIZE + tgx;
#ifdef USE_GK
        AtomData data = loadAtomData(atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
        AtomData data = loadAtomData(atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
        zeroAtomData(&data);
        if (x == y) {
            // This tile is on the diagonal.

            localData[LOCAL_ID].pos = data.pos;
            localData[LOCAL_ID].inducedDipole = data.inducedDipole;
            localData[LOCAL_ID].inducedDipolePolar = data.inducedDipolePolar;
            localData[LOCAL_ID].thole = data.thole;
            localData[LOCAL_ID].damp = data.damp;
#ifdef USE_GK
            localData[LOCAL_ID].inducedDipoleS = data.inducedDipoleS;
            localData[LOCAL_ID].inducedDipolePolarS = data.inducedDipolePolarS;
            localData[LOCAL_ID].bornRadius = data.bornRadius;
#endif
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+j].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(&data, &localData[tbx+j], delta, atom1 == atom2);
            }
            SYNC_WARPS;
        }
        else {
            // This is an off-diagonal tile.

#ifdef USE_GK
            localData[LOCAL_ID] = loadAtomData(y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            localData[LOCAL_ID] = loadAtomData(y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomDataLocal(&localData[LOCAL_ID]);
            unsigned int tj = tgx;
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(&data, &localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            SYNC_WARPS;
        }

        // Write results.

        unsigned int offset = x*TILE_SIZE + tgx;
        SAVE_ATOM_DATA(offset, data)
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            SAVE_ATOM_DATA(offset, localData[LOCAL_ID])
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(mm_long)numTileIndices/totalWarps : warp*(mm_long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(mm_long)numTileIndices/totalWarps : (warp+1)*(mm_long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(mm_long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(mm_long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    LOCAL int atomIndices[THREAD_BLOCK_SIZE];
    LOCAL volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[LOCAL_ID] = -1;
    
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
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

#ifdef USE_GK
            AtomData data = loadAtomData(atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            AtomData data = loadAtomData(atom1, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomData(&data);
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[LOCAL_ID] = j;
#ifdef USE_GK
            localData[LOCAL_ID] = loadAtomData(j, posq, inducedDipole, inducedDipolePolar, dampingAndThole, inducedDipoleS, inducedDipolePolarS, bornRadii);
#else
            localData[LOCAL_ID] = loadAtomData(j, posq, inducedDipole, inducedDipolePolar, dampingAndThole);
#endif
            zeroAtomDataLocal(&localData[LOCAL_ID]);
            SYNC_WARPS;

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(&data, &localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            SAVE_ATOM_DATA(offset, data)
#ifdef USE_CUTOFF
            offset = atomIndices[LOCAL_ID];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            SAVE_ATOM_DATA(offset, localData[LOCAL_ID])
        }
        pos++;
    }
}

KERNEL void recordInducedDipolesForDIIS(GLOBAL const mm_long* RESTRICT fixedField, GLOBAL const mm_long* RESTRICT fixedFieldPolar,
        GLOBAL const float* RESTRICT polarizability, GLOBAL float2* RESTRICT errors, GLOBAL real* RESTRICT prevErrors, GLOBAL real* RESTRICT matrix,
        GLOBAL const mm_long* RESTRICT fixedFieldS, GLOBAL const mm_long* RESTRICT inducedField, GLOBAL const mm_long* RESTRICT inducedFieldPolar,
        GLOBAL const real* RESTRICT inducedDipole, GLOBAL const real* RESTRICT inducedDipolePolar,
        GLOBAL real* RESTRICT prevDipoles, GLOBAL real* RESTRICT prevDipolesPolar, int iteration, int isGK) {
    LOCAL real2 buffer[64];
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        real scale = polarizability[atom];
        if (iteration >= MAX_PREV_DIIS_DIPOLES) {
            // We have filled up the buffer for previous dipoles, so shift them all over by one.

            for (int i = 1; i < MAX_PREV_DIIS_DIPOLES; i++) {
                int index1 = atom+(i-1)*NUM_ATOMS;
                int index2 = atom+i*NUM_ATOMS;
                for (int j = 0; j < 3; j++) {
                    prevDipoles[3*index1+j] = prevDipoles[3*index2+j];
                    prevDipolesPolar[3*index1+j] = prevDipolesPolar[3*index2+j];
                    if (!isGK)
                        prevErrors[3*index1+j] = prevErrors[3*index2+j];
                }
            }
        }

        // Compute the new dipole, and record it along with the error.

        real3 oldDipole = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
        real3 oldDipolePolar = make_real3(inducedDipolePolar[3*atom], inducedDipolePolar[3*atom+1], inducedDipolePolar[3*atom+2]);
        real3 fixed = make_real3(fixedField[atom], fixedField[atom+PADDED_NUM_ATOMS], fixedField[atom+2*PADDED_NUM_ATOMS])*fieldScale;
        real3 fixedPolar = make_real3(fixedFieldPolar[atom], fixedFieldPolar[atom+PADDED_NUM_ATOMS], fixedFieldPolar[atom+2*PADDED_NUM_ATOMS])*fieldScale;
        real3 induced = make_real3(inducedField[atom], inducedField[atom+PADDED_NUM_ATOMS], inducedField[atom+2*PADDED_NUM_ATOMS])*fieldScale;
        real3 inducedPolar = make_real3(inducedFieldPolar[atom], inducedFieldPolar[atom+PADDED_NUM_ATOMS], inducedFieldPolar[atom+2*PADDED_NUM_ATOMS])*fieldScale;
        real3 fixedS = make_real3(0);
        if (isGK)
            fixedS = make_real3(fixedFieldS[atom], fixedFieldS[atom+PADDED_NUM_ATOMS], fixedFieldS[atom+2*PADDED_NUM_ATOMS])*fieldScale;
        real3 newDipole = scale*(fixed+fixedS+induced);
        real3 newDipolePolar = scale*(fixedPolar+fixedS+inducedPolar);
        int storePrevIndex = atom+min(iteration, MAX_PREV_DIIS_DIPOLES-1)*NUM_ATOMS;
        prevDipoles[3*storePrevIndex] = newDipole.x;
        prevDipoles[3*storePrevIndex+1] = newDipole.y;
        prevDipoles[3*storePrevIndex+2] = newDipole.z;
        prevDipolesPolar[3*storePrevIndex] = newDipolePolar.x;
        prevDipolesPolar[3*storePrevIndex+1] = newDipolePolar.y;
        prevDipolesPolar[3*storePrevIndex+2] = newDipolePolar.z;
        if (!isGK) {
            prevErrors[3*storePrevIndex] = newDipole.x-oldDipole.x;
            prevErrors[3*storePrevIndex+1] = newDipole.y-oldDipole.y;
            prevErrors[3*storePrevIndex+2] = newDipole.z-oldDipole.z;
        }
        real3 errors = (newDipole-oldDipole)*(newDipole-oldDipole);
        real3 errorsPolar = (newDipolePolar-oldDipolePolar)*(newDipolePolar-oldDipolePolar);
        sumErrors += errors.x + errors.y + errors.z;
        sumPolarErrors += errorsPolar.x + errorsPolar.y + errorsPolar.z;
    }
    
    // Sum the errors over threads and store the total for this block.
    
    buffer[LOCAL_ID] = make_real2(sumErrors, sumPolarErrors);
    SYNC_THREADS;
    for (int offset = 1; offset < LOCAL_SIZE; offset *= 2) {
        if (LOCAL_ID+offset < LOCAL_SIZE && (LOCAL_ID&(2*offset-1)) == 0) {
            buffer[LOCAL_ID].x += buffer[LOCAL_ID+offset].x;
            buffer[LOCAL_ID].y += buffer[LOCAL_ID+offset].y;
        }
        SYNC_THREADS;
    }
    if (LOCAL_ID == 0)
        errors[GROUP_ID] = make_float2((float) buffer[0].x, (float) buffer[0].y);
    
    if (iteration >= MAX_PREV_DIIS_DIPOLES && !isGK && GROUP_ID == 0) {
        // Shift over the existing matrix elements.
        
        for (int i = 0; i < MAX_PREV_DIIS_DIPOLES-1; i++) {
            if (LOCAL_ID < MAX_PREV_DIIS_DIPOLES-1)
                matrix[LOCAL_ID+i*MAX_PREV_DIIS_DIPOLES] = matrix[(LOCAL_ID+1)+(i+1)*MAX_PREV_DIIS_DIPOLES];
            SYNC_THREADS;
        }
    }
}

KERNEL void computeDIISMatrix(GLOBAL real* RESTRICT prevErrors, int iteration, GLOBAL real* RESTRICT matrix) {
    LOCAL real sumBuffer[512];
    int j = min(iteration, MAX_PREV_DIIS_DIPOLES-1);
    for (int i = GROUP_ID; i <= j; i += NUM_GROUPS) {
        // All the threads in this thread block work together to compute a single matrix element.

        real sum = 0;
        for (int index = LOCAL_ID; index < 3*NUM_ATOMS; index += LOCAL_SIZE)
            sum += prevErrors[index+i*3*NUM_ATOMS]*prevErrors[index+j*3*NUM_ATOMS];
        sumBuffer[LOCAL_ID] = sum;
        SYNC_THREADS;
        for (int offset = 1; offset < LOCAL_SIZE; offset *= 2) { 
            if (LOCAL_ID+offset < LOCAL_SIZE && (LOCAL_ID&(2*offset-1)) == 0)
                sumBuffer[LOCAL_ID] += sumBuffer[LOCAL_ID+offset];
            SYNC_THREADS;
        }
        if (LOCAL_ID == 0) {
            matrix[i+MAX_PREV_DIIS_DIPOLES*j] = sumBuffer[0];
            if (i != j)
                matrix[j+MAX_PREV_DIIS_DIPOLES*i] = sumBuffer[0];
        }
    }
}

KERNEL void solveDIISMatrix(int iteration, GLOBAL const real* RESTRICT matrix, GLOBAL float* RESTRICT coefficients) {
    LOCAL real b[MAX_PREV_DIIS_DIPOLES+1][MAX_PREV_DIIS_DIPOLES+1];
    LOCAL real piv[MAX_PREV_DIIS_DIPOLES+1];
    LOCAL real x[MAX_PREV_DIIS_DIPOLES+1];

    // On the first iteration we don't need to do any calculation.
    
    if (iteration == 0) {
        if (LOCAL_ID == 0)
            coefficients[0] = 1;
        return;
    }
    
    // Load the matrix.
    
    int numPrev = min(iteration+1, MAX_PREV_DIIS_DIPOLES);
    int rank = numPrev+1;
    for (int index = LOCAL_ID; index < numPrev*numPrev; index += LOCAL_SIZE) {
        int i = index/numPrev;
        int j = index-i*numPrev;
        b[i+1][j+1] = matrix[i*MAX_PREV_DIIS_DIPOLES+j];
    }
    for (int i = LOCAL_ID; i < rank; i += LOCAL_SIZE) {
        b[i][0] = -1;
        piv[i] = i;
    }
    SYNC_THREADS;
    
    // Compute the mean absolute value of the values we just loaded.  We use that for preconditioning it,
    // which is essential for doing the computation in single precision.
    
    if (LOCAL_ID == 0) {
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
                if (fabs(b[i][j]) > fabs(b[p][j]))
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

KERNEL void updateInducedFieldByDIIS(GLOBAL real* RESTRICT inducedDipole, GLOBAL real* RESTRICT inducedDipolePolar, 
        GLOBAL const real* RESTRICT prevDipoles, GLOBAL const real* RESTRICT prevDipolesPolar, GLOBAL const float* RESTRICT coefficients, int numPrev) {
    for (int index = GLOBAL_ID; index < 3*NUM_ATOMS; index += GLOBAL_SIZE) {
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

KERNEL void initExtrapolatedDipoles(GLOBAL real* RESTRICT inducedDipole, GLOBAL real* RESTRICT extrapolatedDipole
#ifndef HIPPO
        , GLOBAL real* RESTRICT inducedDipolePolar, GLOBAL real* RESTRICT extrapolatedDipolePolar, GLOBAL mm_long* RESTRICT inducedDipoleFieldGradient, GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientPolar
#endif
#ifdef USE_GK
        , GLOBAL real* RESTRICT inducedDipoleGk, GLOBAL real* RESTRICT inducedDipoleGkPolar, GLOBAL real* RESTRICT extrapolatedDipoleGk, GLOBAL real* RESTRICT extrapolatedDipoleGkPolar,
        GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientGk, GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientGkPolar
#endif
        ) {
    for (int index = GLOBAL_ID; index < 3*NUM_ATOMS; index += GLOBAL_SIZE) {
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
    for (int index = GLOBAL_ID; index < 6*NUM_ATOMS; index += GLOBAL_SIZE) {
        inducedDipoleFieldGradient[index] = 0;
        inducedDipoleFieldGradientPolar[index] = 0;
#ifdef USE_GK
        inducedDipoleFieldGradientGk[index] = 0;
        inducedDipoleFieldGradientGkPolar[index] = 0;
#endif
    }
#endif
}

KERNEL void iterateExtrapolatedDipoles(int order, GLOBAL real* RESTRICT inducedDipole, GLOBAL real* RESTRICT extrapolatedDipole, GLOBAL mm_long* RESTRICT inducedDipoleField,
#ifndef HIPPO
        GLOBAL real* RESTRICT inducedDipolePolar, GLOBAL real* RESTRICT extrapolatedDipolePolar, GLOBAL mm_long* RESTRICT inducedDipoleFieldPolar, GLOBAL mm_long* RESTRICT inducedDipoleFieldGradient,
        GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientPolar, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradient, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientPolar,
#endif
#ifdef USE_GK
        GLOBAL real* RESTRICT inducedDipoleGk, GLOBAL real* RESTRICT inducedDipoleGkPolar, GLOBAL real* RESTRICT extrapolatedDipoleGk, GLOBAL real* RESTRICT extrapolatedDipoleGkPolar,
        GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientGk, GLOBAL mm_long* RESTRICT inducedDipoleFieldGradientGkPolar, GLOBAL mm_long* RESTRICT inducedDipoleFieldGk,
        GLOBAL mm_long* RESTRICT inducedDipoleFieldGkPolar, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientGk, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientGkPolar,
#endif
#ifdef HIPPO
        GLOBAL const real* RESTRICT polarizability
#else
        GLOBAL const float* RESTRICT polarizability
#endif
    ) {
    const real fieldScale = 1/(real) 0x100000000;
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
        float polar = polarizability[atom];
        real3 value = make_real3(inducedDipoleField[atom], inducedDipoleField[atom+PADDED_NUM_ATOMS], inducedDipoleField[atom+2*PADDED_NUM_ATOMS])*fieldScale*polar;
        inducedDipole[3*atom] = value.x;
        inducedDipole[3*atom+1] = value.y;
        inducedDipole[3*atom+2] = value.z;
        extrapolatedDipole[3*(order*NUM_ATOMS+atom)] = value.x;
        extrapolatedDipole[3*(order*NUM_ATOMS+atom)+1] = value.y;
        extrapolatedDipole[3*(order*NUM_ATOMS+atom)+2] = value.z;
#ifndef HIPPO
        value = make_real3(inducedDipoleFieldPolar[atom], inducedDipoleFieldPolar[atom+PADDED_NUM_ATOMS], inducedDipoleFieldPolar[atom+2*PADDED_NUM_ATOMS])*fieldScale*polar;
        inducedDipolePolar[3*atom] = value.x;
        inducedDipolePolar[3*atom+1] = value.y;
        inducedDipolePolar[3*atom+2] = value.z;
        extrapolatedDipolePolar[3*(order*NUM_ATOMS+atom)] = value.x;
        extrapolatedDipolePolar[3*(order*NUM_ATOMS+atom)+1] = value.y;
        extrapolatedDipolePolar[3*(order*NUM_ATOMS+atom)+2] = value.z;
#endif
#ifdef USE_GK
        value = make_real3(inducedDipoleFieldGk[atom], inducedDipoleFieldGk[atom+PADDED_NUM_ATOMS], inducedDipoleFieldGk[atom+2*PADDED_NUM_ATOMS])*fieldScale*polar;
        inducedDipoleGk[3*atom] = value.x;
        inducedDipoleGk[3*atom+1] = value.y;
        inducedDipoleGk[3*atom+2] = value.z;
        extrapolatedDipoleGk[3*(order*NUM_ATOMS+atom)] = value.x;
        extrapolatedDipoleGk[3*(order*NUM_ATOMS+atom)+1] = value.y;
        extrapolatedDipoleGk[3*(order*NUM_ATOMS+atom)+2] = value.z;
        value = make_real3(inducedDipoleFieldGkPolar[atom], inducedDipoleFieldGkPolar[atom+PADDED_NUM_ATOMS], inducedDipoleFieldGkPolar[atom+2*PADDED_NUM_ATOMS])*fieldScale*polar;
        inducedDipoleGkPolar[3*atom] = value.x;
        inducedDipoleGkPolar[3*atom+1] = value.y;
        inducedDipoleGkPolar[3*atom+2] = value.z;
        extrapolatedDipoleGkPolar[3*(order*NUM_ATOMS+atom)] = value.x;
        extrapolatedDipoleGkPolar[3*(order*NUM_ATOMS+atom)+1] = value.y;
        extrapolatedDipoleGkPolar[3*(order*NUM_ATOMS+atom)+2] = value.z;
#endif
    }
#ifndef HIPPO
    for (int index = GLOBAL_ID; index < 2*NUM_ATOMS; index += GLOBAL_SIZE) {
        int index2 = (order-1)*2*NUM_ATOMS+index;
        for (int i = 0; i < 3; i++) {
            extrapolatedDipoleFieldGradient[3*index2+i] = fieldScale*inducedDipoleFieldGradient[3*index+i];
            extrapolatedDipoleFieldGradientPolar[3*index2+i] = fieldScale*inducedDipoleFieldGradientPolar[3*index+i];
#ifdef USE_GK
            extrapolatedDipoleFieldGradientGk[3*index2+i] = fieldScale*inducedDipoleFieldGradientGk[3*index+i];
            extrapolatedDipoleFieldGradientGkPolar[3*index2+i] = fieldScale*inducedDipoleFieldGradientGkPolar[3*index+i];
#endif
        }
    }
#endif
}

KERNEL void computeExtrapolatedDipoles(GLOBAL real* RESTRICT inducedDipole, GLOBAL real* RESTRICT extrapolatedDipole
#ifndef HIPPO
        , GLOBAL real* RESTRICT inducedDipolePolar, GLOBAL real* RESTRICT extrapolatedDipolePolar
#endif
#ifdef USE_GK
        , GLOBAL real* RESTRICT inducedDipoleGk, GLOBAL real* RESTRICT inducedDipoleGkPolar, GLOBAL real* RESTRICT extrapolatedDipoleGk, GLOBAL real* RESTRICT extrapolatedDipoleGkPolar
#endif
        ) {
    real coeff[] = {EXTRAPOLATION_COEFFICIENTS_SUM};
    for (int index = GLOBAL_ID; index < 3*NUM_ATOMS; index += GLOBAL_SIZE) {
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

KERNEL void addExtrapolatedFieldGradientToForce(GLOBAL mm_long* RESTRICT forceBuffers, GLOBAL real* RESTRICT extrapolatedDipole,
        GLOBAL real* RESTRICT extrapolatedDipolePolar, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradient, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientPolar
#ifdef USE_GK
        , GLOBAL real* RESTRICT extrapolatedDipoleGk, GLOBAL real* RESTRICT extrapolatedDipoleGkPolar,
        GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientGk, GLOBAL real* RESTRICT extrapolatedDipoleFieldGradientGkPolar
#endif
        ) {
    real coeff[] = {EXTRAPOLATION_COEFFICIENTS_SUM};
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE) {
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
        forceBuffers[atom] += (mm_long) (fx*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS] += (mm_long) (fy*0x100000000);
        forceBuffers[atom+PADDED_NUM_ATOMS*2] += (mm_long) (fz*0x100000000);
    }
}

#ifdef HIPPO
KERNEL void computePolarizationEnergy(GLOBAL mixed* RESTRICT energyBuffer, GLOBAL const real* RESTRICT inducedDipole,
        GLOBAL const real* RESTRICT extrapolatedDipole, GLOBAL const real* RESTRICT polarizability) {
    mixed energy = 0;
    for (int atom = GLOBAL_ID; atom < NUM_ATOMS; atom += GLOBAL_SIZE)
        for (int i = 0; i < 3; i++)
            energy -= (ENERGY_SCALE_FACTOR/2)*extrapolatedDipole[3*atom+i]*inducedDipole[3*atom+i]/polarizability[atom];
    energyBuffer[GLOBAL_ID] += energy;
}
#endif