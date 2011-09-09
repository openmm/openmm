{
    float invRSquaredOver4 = 0.25f*invR*invR;
    float rScaledRadiusJ = r+obcParams2.y;
    float rScaledRadiusI = r+obcParams1.y;
    float l_ijJ = RECIP(max(obcParams1.x, fabs(r-obcParams2.y)));
    float l_ijI = RECIP(max(obcParams2.x, fabs(r-obcParams1.y)));
    float u_ijJ = RECIP(rScaledRadiusJ);
    float u_ijI = RECIP(rScaledRadiusI);
    float l_ij2J = l_ijJ*l_ijJ;
    float l_ij2I = l_ijI*l_ijI;
    float u_ij2J = u_ijJ*u_ijJ;
    float u_ij2I = u_ijI*u_ijI;
    float t1J = LOG(u_ijJ*RECIP(l_ijJ));
    float t1I = LOG(u_ijI*RECIP(l_ijI));
    float t2J = (l_ij2J-u_ij2J);
    float t2I = (l_ij2I-u_ij2I);
    float term1 = (0.5f*(0.25f+obcParams2.y*obcParams2.y*invRSquaredOver4)*t2J + t1J*invRSquaredOver4)*invR;
    float term2 = (0.5f*(0.25f+obcParams1.y*obcParams1.y*invRSquaredOver4)*t2I + t1I*invRSquaredOver4)*invR;
    float tempdEdR = select(0.0f, bornForce1*term1, obcParams1.x < rScaledRadiusJ);
    tempdEdR += select(0.0f, bornForce2*term2, obcParams2.x < rScaledRadiusJ);
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2);
#endif
    dEdR += select(0.0f, tempdEdR, includeInteraction);
}
