#ifdef USE_CUTOFF
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED) {
#else
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
#ifdef USE_SYMMETRIC
    float tempForce = 0.0f;
#else
    float4 tempForce1 = (float4) 0.0f;
    float4 tempForce2 = (float4) 0.0f;
#endif
    COMPUTE_FORCE
#ifdef USE_SYMMETRIC
    dEdR += tempForce*invR;
#else
    dEdR1 += tempForce1;
    dEdR2 += tempForce2;
#endif
}
