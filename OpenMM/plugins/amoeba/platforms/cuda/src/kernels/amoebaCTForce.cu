#ifdef USE_CUTOFF
unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
unsigned int includeInteraction = (!isExcluded);
#endif

if (atom2 < PADDED_NUM_ATOMS) {
int pairK = CTTypes1 * NUM_CTPR_TYPES + CTTypes2;
real apre = apreBexp[pairK].x;
real bexp = apreBexp[pairK].y;
real termEnergy = -apre*1000.0*EXP(-bexp*r);
real deltaE = -bexp*termEnergy; 
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
