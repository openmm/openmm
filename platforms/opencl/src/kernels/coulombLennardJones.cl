#if USE_EWALD
bool needCorrection = isExcluded && atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS;
if (!isExcluded || needCorrection) {
    const float prefactor = 138.935456f*posq1.w*posq2.w*invR;
    float alphaR = EWALD_ALPHA*r;
    float erfcAlphaR = 0.0f;
    if (r2 < CUTOFF_SQUARED) {
        float normalized = ERFC_TABLE_SCALE*alphaR;
        int tableIndex = (int) normalized;
        float fract2 = normalized-tableIndex;
        float fract1 = 1.0f-fract2;
        erfcAlphaR = fract1*erfcTable[tableIndex] + fract2*erfcTable[tableIndex+1];
    }
    else if (needCorrection)
        erfcAlphaR = erfc(alphaR);
    float tempForce = 0.0f;
    if (needCorrection) {
        // Subtract off the part of this interaction that was included in the reciprocal space contribution.

        tempForce = -prefactor*((1.0f-erfcAlphaR)-alphaR*EXP(-alphaR*alphaR)*TWO_OVER_SQRT_PI);
        tempEnergy += -prefactor*(1.0f-erfcAlphaR);
    }
    else if (r2 < CUTOFF_SQUARED) {
#if HAS_LENNARD_JONES
        float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
        float sig2 = invR*sig;
        sig2 *= sig2;
        float sig6 = sig2*sig2*sig2;
        float eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
        tempForce = eps*(12.0f*sig6 - 6.0f)*sig6 + prefactor*(erfcAlphaR+alphaR*EXP(-alphaR*alphaR)*TWO_OVER_SQRT_PI);
        tempEnergy += eps*(sig6 - 1.0f)*sig6 + prefactor*erfcAlphaR;
#else
        tempForce = prefactor*(erfcAlphaR+alphaR*EXP(-alphaR*alphaR)*TWO_OVER_SQRT_PI);
        tempEnergy += prefactor*erfcAlphaR;
#endif
    }
    dEdR += tempForce*invR*invR;
}
#else
{
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (!isExcluded);
#endif
    float tempForce = 0.0f;
  #if HAS_LENNARD_JONES
    float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
    float sig2 = invR*sig;
    sig2 *= sig2;
    float sig6 = sig2*sig2*sig2;
    float eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
    tempForce = eps*(12.0f*sig6 - 6.0f)*sig6;
    tempEnergy += select(0.0f, eps*(sig6 - 1.0f)*sig6, includeInteraction);
  #endif
#if HAS_COULOMB
  #ifdef USE_CUTOFF
    const float prefactor = 138.935456f*posq1.w*posq2.w;
    tempForce += prefactor*(invR - 2.0f*REACTION_FIELD_K*r2);
    tempEnergy += select(0.0f, prefactor*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C), includeInteraction);
  #else
    const float prefactor = 138.935456f*posq1.w*posq2.w*invR;
    tempForce += prefactor;
    tempEnergy += select(0.0f, prefactor, includeInteraction);
  #endif
#endif
    dEdR += select(0.0f, tempForce*invR*invR, includeInteraction);
}
#endif