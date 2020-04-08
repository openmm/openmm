// vdwpr
{
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (!isExcluded);
#endif

if (atom2 < NUM_ATOMS) {
    int pairK = condensedTypes1 * NUM_VDW_TYPES + condensedTypes2;
    real sigma = pairSigmaEpsilon[pairK].x;
    real epsilon = pairSigmaEpsilon[pairK].y;

#if USE_LJ
    // LENNARD-JONES
    real pp1 = sigma * invR;
    real pp2 = pp1 * pp1;
    real pp3 = pp2 * pp1;
    real pp6 = pp3 * pp3;
    real pp12 = pp6 * pp6;
    real termEnergy = epsilon * (pp12 - 2.0f * pp6);
    real deltaE = epsilon * (pp12 - pp6) * (-12.0f) * invR;

#else
    // BUFFERED 14-7

    real softcore = 0.0f;
#if VDW_ALCHEMICAL_METHOD == 1
    if (isAlchemical1 != isAlchemical2) {
#elif VDW_ALCHEMICAL_METHOD == 2
    if (isAlchemical1 || isAlchemical2) {
#endif
#if VDW_ALCHEMICAL_METHOD != 0
       real lambda = vdwLambda[0];
       epsilon = epsilon * POW(lambda, VDW_SOFTCORE_POWER);
       softcore = VDW_SOFTCORE_ALPHA * (1.0f - lambda) * (1.0f - lambda);
    }
#endif
    real dhal = 0.07f;
    real ghal = 0.12f;
    real dhal1 = 1.07f;
    real ghal1 = 1.12f;
    real rho = r / sigma;
    real rho2 = rho * rho;
    real rho6 = rho2 * rho2 * rho2;
    real rhoplus = rho + dhal;
    real rhodec2 = rhoplus * rhoplus;
    real rhodec = rhodec2 * rhodec2 * rhodec2;
    real s1 = 1.0f / (softcore + rhodec * rhoplus);
    real s2 = 1.0f / (softcore + rho6 * rho + ghal);
    real point72 = dhal1 * dhal1;
    real t1 = dhal1 * point72 * point72 * point72 * s1;
    real t2 = ghal1 * s2;
    real t2min = t2 - 2.0f;
    real dt1 = -7.0f * rhodec * t1 * s1;
    real dt2 = -7.0f * rho6 * t2 * s2;
    real termEnergy = epsilon * t1 * t2min;
    real deltaE = epsilon * (dt1 * t2min + t1 * dt2) / sigma;
#endif

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
}
