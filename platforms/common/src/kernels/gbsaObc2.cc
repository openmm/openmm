{
    real invRSquaredOver4 = 0.25f*invR*invR;
    real rScaledRadiusJ = r+OBC_PARAMS2.y;
    real rScaledRadiusI = r+OBC_PARAMS1.y;
    real l_ijJ = RECIP(max((real) OBC_PARAMS1.x, fabs(r-OBC_PARAMS2.y)));
    real l_ijI = RECIP(max((real) OBC_PARAMS2.x, fabs(r-OBC_PARAMS1.y)));
    real u_ijJ = RECIP(rScaledRadiusJ);
    real u_ijI = RECIP(rScaledRadiusI);
    real l_ij2J = l_ijJ*l_ijJ;
    real l_ij2I = l_ijI*l_ijI;
    real u_ij2J = u_ijJ*u_ijJ;
    real u_ij2I = u_ijI*u_ijI;
    real t1J = LOG(u_ijJ*RECIP(l_ijJ));
    real t1I = LOG(u_ijI*RECIP(l_ijI));
    real t2J = (l_ij2J-u_ij2J);
    real t2I = (l_ij2I-u_ij2I);
    real term1 = (0.5f*(0.25f+OBC_PARAMS2.y*OBC_PARAMS2.y*invRSquaredOver4)*t2J + t1J*invRSquaredOver4)*invR;
    real term2 = (0.5f*(0.25f+OBC_PARAMS1.y*OBC_PARAMS1.y*invRSquaredOver4)*t2I + t1I*invRSquaredOver4)*invR;
#ifdef SUPPORTS_64_BIT_ATOMICS
    real tempdEdR = (OBC_PARAMS1.x < rScaledRadiusJ ? BORN_FORCE1*term1/0x100000000 : 0);
    tempdEdR += (OBC_PARAMS2.x < rScaledRadiusI ? BORN_FORCE2*term2/0x100000000 : 0);
#else
    real tempdEdR = (OBC_PARAMS1.x < rScaledRadiusJ ? BORN_FORCE1*term1 : (real) 0);
    tempdEdR += (OBC_PARAMS2.x < rScaledRadiusI ? BORN_FORCE2*term2 : (real) 0);
#endif
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2);
#endif
    dEdR += (includeInteraction ? tempdEdR : (real) 0);
}
