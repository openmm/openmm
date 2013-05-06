#ifdef USE_CUTOFF
if (!isExcluded && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded) {
#endif
    real tempForce = 0;
    COMPUTE_FORCE
#if USE_SWITCH
    if (r > SWITCH_CUTOFF) {
        real x = r-SWITCH_CUTOFF;
        real switchValue = 1+x*x*x*(SWITCH_C3+x*(SWITCH_C4+x*SWITCH_C5));
        real switchDeriv = x*x*(3*SWITCH_C3+x*(4*SWITCH_C4+x*5*SWITCH_C5));
        tempForce = tempForce*switchValue - tempEnergy*switchDeriv;
        tempEnergy *= switchValue;
    }
#endif
    dEdR += tempForce*invR;
}
