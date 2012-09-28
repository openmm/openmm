#ifdef USE_CUTOFF
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED) {
#else
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
#ifdef USE_SYMMETRIC
    real tempForce = 0;
#else
    real3 tempForce1 = make_real3(0);
    real3 tempForce2 = make_real3(0);
#endif
    COMPUTE_FORCE
#ifdef USE_SYMMETRIC
    dEdR += tempForce*invR;
#else
    dEdR1 += tempForce1;
    dEdR2 += tempForce2;
#endif
}
