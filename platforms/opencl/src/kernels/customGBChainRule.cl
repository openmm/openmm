#ifdef USE_CUTOFF
if (!isExcluded && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED) {
#else
if (!isExcluded && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
    float tempForce = 0.0f;
    COMPUTE_FORCE
    dEdR += tempForce*invR;
}
