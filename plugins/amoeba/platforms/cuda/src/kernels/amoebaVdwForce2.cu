#ifdef USE_CUTOFF
unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
unsigned int includeInteraction = (!isExcluded);
#endif

if (atom2 < PADDED_NUM_ATOMS) {
int pairK = vdwTypes1 * NUM_VDWPR_TYPES + vdwTypes2;
real sigma = sigmaEpsilon[pairK].x;
real epsilon = sigmaEpsilon[pairK].y;

#if FUNCTIONAL_FORM == 1
#if COUPLE_METHOD == 0
real combinedLambda = (lambdas1 == lambdas2 ? 1.0f : (lambdas1 < lambdas2 ? lambdas1 : lambdas2));
#endif /* decouple */
#if COUPLE_METHOD == 1
real combinedLambda = (lambdas1 < lambdas2 ? lambdas1 : lambdas2);
#endif /* annihilate */
real comblambda2 = combinedLambda * combinedLambda;
epsilon = epsilon * comblambda2 * comblambda2 * combinedLambda;
real r_sigma = RECIP(sigma);
real rho = r * r_sigma;
real rho2 = rho * rho;
real rho6 = rho2 * rho2 * rho2;
real rhoplus = rho + 0.07f;
real rhodec2 = rhoplus * rhoplus;
real rhodec = rhodec2 * rhodec2 * rhodec2;
real scal = 0.7f * (1.0f - combinedLambda) * (1.0f - combinedLambda);
real s1 = RECIP(scal + rhodec * rhoplus);
real s2 = RECIP(scal + rho6 * rho + 0.12f);
real point72 = 1.60578147648f;
real t1 = point72 * s1;
real t2 = 1.12f * s2;
real t2min = t2 - 2.0f;
real dt1 = -7.0f * rhodec * t1 * s1;
real dt2 = -7.0f * rho6 * t2 * s2;
real termEnergy = epsilon * t1 * t2min;
real deltaE = epsilon * (dt1 * t2min + t1 * dt2) * r_sigma;
#endif /* FUNCTIONAL_FORM == 1 BUFFERED-14-7 */

#if FUNCTIONAL_FORM == 2
real pp1 = sigma * invR;
real pp2 = pp1 * pp1;
real pp3 = pp2 * pp1;
real pp6 = pp3 * pp3;
real pp12 = pp6 * pp6;
real termEnergy = epsilon * (pp12 - 2.0f * pp6);
real deltaE = epsilon * (pp12 - pp6) * (-12.0f) * invR;
#endif /* FUNCTIONAL_FORM == 2 LENNARD-JONES */

#ifdef USE_CUTOFF
if (r > TAPER_CUTOFF) {
    real x = r - TAPER_CUTOFF;
    real taper = 1 + x * x * x * (TAPER_C3 + x * (TAPER_C4 + x * TAPER_C5));
    real dtaper = x * x * (3 * TAPER_C3 + x * (4 * TAPER_C4 + x * 5 * TAPER_C5));
    deltaE = termEnergy * dtaper + deltaE * taper;
    termEnergy *= taper;
}
#endif
tempEnergy += (includeInteraction ? termEnergy : 0);
dEdR -= (includeInteraction ? deltaE * invR : 0);
}
