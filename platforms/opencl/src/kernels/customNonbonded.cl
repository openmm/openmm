#ifdef USE_CUTOFF
if (!isExcluded && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded) {
#endif
    real tempForce = 0.0f;
    COMPUTE_FORCE
    dEdR += tempForce*invR;
}
