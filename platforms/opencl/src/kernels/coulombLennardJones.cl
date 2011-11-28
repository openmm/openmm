#if USE_EWALD
bool needCorrection = isExcluded && atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS;
if (!isExcluded || needCorrection) {
    float tempForce = 0.0f;
    if (r2 < CUTOFF_SQUARED || needCorrection) {
        const float alphaR = EWALD_ALPHA*r;
        const float expAlphaRSqr = EXP(-alphaR*alphaR);
        const float prefactor = 138.935456f*posq1.w*posq2.w*invR;

        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 3e-7.

        float t = 1.0f+(0.0705230784f+(0.0422820123f+(0.0092705272f+(0.0001520143f+(0.0002765672f+0.0000430638f*alphaR)*alphaR)*alphaR)*alphaR)*alphaR)*alphaR;
        t *= t;
        t *= t;
        t *= t;
        const float erfcAlphaR = RECIP(t*t);
        if (needCorrection) {
            // Subtract off the part of this interaction that was included in the reciprocal space contribution.

            tempForce = -prefactor*((1.0f-erfcAlphaR)-alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
            tempEnergy += -prefactor*(1.0f-erfcAlphaR);
        }
        else {
#if HAS_LENNARD_JONES
            float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
            float sig2 = invR*sig;
            sig2 *= sig2;
            float sig6 = sig2*sig2*sig2;
            float epssig6 = sig6*(sigmaEpsilon1.y*sigmaEpsilon2.y);
            tempForce = epssig6*(12.0f*sig6 - 6.0f) + prefactor*(erfcAlphaR+alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
            tempEnergy += epssig6*(sig6 - 1.0f) + prefactor*erfcAlphaR;
#else
            tempForce = prefactor*(erfcAlphaR+alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
            tempEnergy += prefactor*erfcAlphaR;
#endif
        }
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
    float epssig6 = sig6*(sigmaEpsilon1.y*sigmaEpsilon2.y);
    tempForce = epssig6*(12.0f*sig6 - 6.0f);
    tempEnergy += select(0.0f, epssig6*(sig6 - 1.0f), includeInteraction);
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