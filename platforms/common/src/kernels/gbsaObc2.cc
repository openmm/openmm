{
#if USE_GBSA_CHAIN_RULE_GUARD
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2);
#endif
    // Only pair arithmetic is conditional; enclosing SIMD communication stays uniform.
    if (includeInteraction) {
#endif
    real invRSquaredOver4 = 0.25f*invR*invR;
    real rScaledRadiusJ = r+OBC_PARAMS2.y;
    real rScaledRadiusI = r+OBC_PARAMS1.y;
#if USE_GBSA_RECIPROCAL_REUSE && !defined(OPENMM_METAL_FLOAT_ACCUMULATORS) && !defined(OPENMM_METAL_REQUIRE_SAFE_MATH)
    real lowerBoundJ = max((real) OBC_PARAMS1.x, fabs(r-OBC_PARAMS2.y));
    real lowerBoundI = max((real) OBC_PARAMS2.x, fabs(r-OBC_PARAMS1.y));
    real l_ijJ = RECIP(lowerBoundJ);
    real l_ijI = RECIP(lowerBoundI);
#else
    real l_ijJ = RECIP(max((real) OBC_PARAMS1.x, fabs(r-OBC_PARAMS2.y)));
    real l_ijI = RECIP(max((real) OBC_PARAMS2.x, fabs(r-OBC_PARAMS1.y)));
#endif
    real u_ijJ = RECIP(rScaledRadiusJ);
    real u_ijI = RECIP(rScaledRadiusI);
    real l_ij2J = l_ijJ*l_ijJ;
    real l_ij2I = l_ijI*l_ijI;
    real u_ij2J = u_ijJ*u_ijJ;
    real u_ij2I = u_ijI*u_ijI;
#if USE_GBSA_RECIPROCAL_REUSE && !defined(OPENMM_METAL_FLOAT_ACCUMULATORS) && !defined(OPENMM_METAL_REQUIRE_SAFE_MATH)
    real t1J = LOG(u_ijJ*lowerBoundJ);
    real t1I = LOG(u_ijI*lowerBoundI);
#else
    real t1J = LOG(u_ijJ*RECIP(l_ijJ));
    real t1I = LOG(u_ijI*RECIP(l_ijI));
#endif
    real t2J = (l_ij2J-u_ij2J);
    real t2I = (l_ij2I-u_ij2I);
    real term1 = (0.5f*(0.25f+OBC_PARAMS2.y*OBC_PARAMS2.y*invRSquaredOver4)*t2J + t1J*invRSquaredOver4)*invR;
    real term2 = (0.5f*(0.25f+OBC_PARAMS1.y*OBC_PARAMS1.y*invRSquaredOver4)*t2I + t1I*invRSquaredOver4)*invR;
    real tempdEdR = (OBC_PARAMS1.x < rScaledRadiusJ ? BORN_FORCE1*term1/0x100000000 : 0);
    tempdEdR += (OBC_PARAMS2.x < rScaledRadiusI ? BORN_FORCE2*term2/0x100000000 : 0);
#if !USE_GBSA_CHAIN_RULE_GUARD
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2);
#endif
#endif
    dEdR += (includeInteraction ? tempdEdR : (real) 0);
#if USE_GBSA_CHAIN_RULE_GUARD
    }
#endif
}
