#ifdef USE_CUTOFF
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2 && r2 < CUTOFF_SQUARED) {
#else
if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
    float invRSquared = RECIP(r2);
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
    float t1J = LOG(u_ijJ/l_ijJ);
    float t1I = LOG(u_ijI/l_ijI);
    float t2J = (l_ij2J-u_ij2J);
    float t2I = (l_ij2I-u_ij2I);
    float t3J = t2J*invR;
    float t3I = t2I*invR;
    t1J *= invR;
    t1I *= invR;
    float term1 = 0.125f*(1.0f+obcParams2.y*obcParams2.y*invRSquared)*t3J + 0.25f*t1J*invRSquared;
    float term2 = 0.125f*(1.0f+obcParams1.y*obcParams1.y*invRSquared)*t3I + 0.25f*t1I*invRSquared;
    dEdR += (obcParams1.x < rScaledRadiusJ ? bornForce1*term1 : 0.0f);
    dEdR += (obcParams2.x < rScaledRadiusJ ? bornForce2*term2 : 0.0f);
}
