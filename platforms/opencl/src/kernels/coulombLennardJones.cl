#if USE_EWALD
if (r2 < CUTOFF_SQUARED) {
    bool needCorrection = isExcluded && atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS;
    if (!isExcluded || needCorrection) {
        const float prefactor = 138.935485f*posq1.w*posq2.w*invR;
        float alphaR = EWALD_ALPHA*r;
        float erfcAlphaR = erfc(alphaR);
        float tempForce;
        if (needCorrection) {
            // Subtract off the part of this interaction that was included in the reciprocal space contribution.

            tempForce = -prefactor*((1.0f-erfcAlphaR)-alphaR*exp(-alphaR*alphaR)*TWO_OVER_SQRT_PI);
            tempEnergy += -prefactor*(1.0f-erfcAlphaR);
        }
        else {
            float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
            float sig2 = invR*sig;
            sig2 *= sig2;
            float sig6 = sig2*sig2*sig2;
            float eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
            tempForce = eps*(12.0f*sig6 - 6.0f)*sig6 + prefactor*(erfcAlphaR+alphaR*exp(-alphaR*alphaR)*TWO_OVER_SQRT_PI);
            tempEnergy += eps*(sig6 - 1.0f)*sig6 + prefactor*erfcAlphaR;
        }
        dEdR += tempForce*invR*invR;
    }
}
#else
#ifdef USE_CUTOFF
if (!isExcluded && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded) {
#endif
    float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
    float sig2 = invR*sig;
    sig2 *= sig2;
    float sig6 = sig2*sig2*sig2;
    float eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
    float tempForce = eps*(12.0f*sig6 - 6.0f)*sig6;
    tempEnergy += eps*(sig6 - 1.0f)*sig6;
#ifdef USE_CUTOFF
    const float prefactor = 138.935485f*posq1.w*posq2.w;
    tempForce += prefactor*(invR - 2.0f*REACTION_FIELD_K*r2);
    tempEnergy += prefactor*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C);
#else
    const float prefactor = 138.935485f*posq1.w*posq2.w*invR;
    tempForce += prefactor;
    tempEnergy += prefactor;
#endif
    dEdR += tempForce*invR*invR;
}
#endif