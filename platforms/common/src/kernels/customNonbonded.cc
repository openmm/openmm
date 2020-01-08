#ifdef USE_CUTOFF
if (!isExcluded && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded) {
#endif
    real tempForce = 0;
    real switchValue = 1, switchDeriv = 0;
#if USE_SWITCH
    if (r > SWITCH_CUTOFF) {
        real x = r-SWITCH_CUTOFF;
        switchValue = 1+x*x*x*(SWITCH_C3+x*(SWITCH_C4+x*SWITCH_C5));
        switchDeriv = x*x*(3*SWITCH_C3+x*(4*SWITCH_C4+x*5*SWITCH_C5));
    }
#endif
    COMPUTE_FORCE
#if USE_SWITCH
    tempForce = tempForce*switchValue - customEnergy*switchDeriv;
    tempEnergy += customEnergy*switchValue;
#else
    tempEnergy += customEnergy;
#endif
    dEdR += tempForce*invR;
}
