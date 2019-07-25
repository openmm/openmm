{
    real invRSquaredOver4 = 0.25f*invR*invR;
    real rScaledRadiusJ = r+obcParams2.y;
    real rScaledRadiusI = r+obcParams1.y;
    real l_ijJ = RECIP(max((real) obcParams1.x, fabs(r-obcParams2.y)));
    real l_ijI = RECIP(max((real) obcParams2.x, fabs(r-obcParams1.y)));
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
    real term1 = (0.5f*(0.25f+obcParams2.y*obcParams2.y*invRSquaredOver4)*t2J + t1J*invRSquaredOver4)*invR;
    real term2 = (0.5f*(0.25f+obcParams1.y*obcParams1.y*invRSquaredOver4)*t2I + t1I*invRSquaredOver4)*invR;
    real tempdEdR = (obcParams1.x < rScaledRadiusJ ? bornForce1*term1 : (real) 0);
    tempdEdR += (obcParams2.x < rScaledRadiusI ? bornForce2*term2 : (real) 0);
#ifdef USE_CUTOFF
    bool includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED);
#else
    bool includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2);
#endif
    dEdR += (includeInteraction ? tempdEdR : (real) 0);
}
