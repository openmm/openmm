
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (!isExcluded);
#endif
    real tempForce = 0.0f;
#if SIGMA_COMBINING_RULE == 1
    real sigma = sigmaEpsilon1.x + sigmaEpsilon2.x;
#elif SIGMA_COMBINING_RULE == 2
    real sigma = 2*SQRT(sigmaEpsilon1.x*sigmaEpsilon2.x);
#else
    real sigma1_2 = sigmaEpsilon1.x*sigmaEpsilon1.x;
    real sigma2_2 = sigmaEpsilon2.x*sigmaEpsilon2.x;
    real sigmasum = sigma1_2+sigma2_2;
    real sigma = (sigmasum == 0.0f ? (real) 0 : 2*(sigmaEpsilon1.x*sigma1_2 + sigmaEpsilon2.x*sigma2_2)/(sigma1_2+sigma2_2));
#endif
#if EPSILON_COMBINING_RULE == 1
    real epsilon = 0.5f*(sigmaEpsilon1.y + sigmaEpsilon2.y);
#elif EPSILON_COMBINING_RULE == 2
    real epsilon = SQRT(sigmaEpsilon1.y*sigmaEpsilon2.y);
#elif EPSILON_COMBINING_RULE == 3
    real epssum = sigmaEpsilon1.y+sigmaEpsilon2.y;
    real epsilon = (epssum == 0.0f ? (real) 0 : 2*(sigmaEpsilon1.y*sigmaEpsilon2.y)/(sigmaEpsilon1.y+sigmaEpsilon2.y));
#else
    real epsilon_s = SQRT(sigmaEpsilon1.y) + SQRT(sigmaEpsilon2.y);
    real epsilon = (epsilon_s == 0.0f ? (real) 0 : 4*sigmaEpsilon1.y*sigmaEpsilon2.y/(epsilon_s*epsilon_s));
#endif

    real softcore = 0.0;
#if USE_SOFTCORE == 1
    // Decouple
    bool both = isAlchemicalI && isAlchemicalJ;
    bool either = isAlchemicalI || isAlchemicalJ;
    bool soft = !both && either;
    if (soft) { 
       epsilon = epsilon * pow(VDW_LAMBDA, VDW_SOFTCORE_POWER);
       softcore = VDW_SOFTCORE_ALPHA * (1.0 - VDW_LAMBDA) * (1.0 - VDW_LAMBDA);
    }
#else if USE_SOFTCORE == 2 
    // Annihilate
    bool soft = isAlchemicalI || isAlchemicalJ;
    if (soft) {
       epsilon = epsilon * pow(VDW_LAMBDA, VDW_SOFTCORE_POWER);
       softcore = VDW_SOFTCORE_ALPHA * (1.0 - VDW_LAMBDA) * (1.0 - VDW_LAMBDA);
    }
@endif

    // The Buffered-14-7 buffering constants.
    real dhal = 0.07;
    real ghal = 0.12;
    // The buffering constant plus one.
    real dhal1 = 1.07;
    real ghal1 = 1.12;

    real rho = r / sigma;
    real rho2 = rho * rho;
    real rho6 = rho2 * rho2 * rho2;
    real rhoplus = rho + dhal;
    real rhodec2 = rhoplus * rhoplus;
    real rhodec = rhodec2 * rhodec2 * rhodec2;
    real s1 = 1.0 / (softcore + rhodec * rhoplus);
    real s2 = 1.0 / (softcore + rho6 * rho + 0.12);
    real point72 = dhal1 * dhal1;
    real t1 = dhal1 * point72 * point72 * point72 * s1;
    real t2 = ghal1 * s2;
    real t2min = t2 - 2;
    real dt1 = -7.0 * rhodec * t1 * s1;
    real dt2 = -7.0 * rho6 * t2 * s2;
    real termEnergy = epsilon * t1 * t2min;
    real deltaE = epsilon * (dt1 * t2min + t1 * dt2) / sigma;

#ifdef USE_CUTOFF
    if (r > TAPER_CUTOFF) {
        real x = r-TAPER_CUTOFF;
        real taper = 1+x*x*x*(TAPER_C3+x*(TAPER_C4+x*TAPER_C5));
        real dtaper = x*x*(3*TAPER_C3+x*(4*TAPER_C4+x*5*TAPER_C5));
        deltaE = termEnergy*dtaper + deltaE*taper;
        termEnergy *= taper;
    }
#endif
    tempEnergy += (includeInteraction ? termEnergy : 0);
    dEdR -= (includeInteraction ? deltaE*invR : 0);
}
