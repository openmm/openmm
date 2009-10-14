#ifdef USE_CUTOFF
if (atom1 < numAtoms && atom2 < numAtoms && atom1 != atom2 && r2 < cutoffSquared) {
#else
if (atom1 < numAtoms && atom2 < numAtoms && atom1 != atom2) {
#endif
    float invRSquared = 1.0f/r2;
    float rScaledRadiusJ = r+obcParams2.y;
    float rScaledRadiusI = r+obcParams1.y;
    float l_ijJ = 1.0f/max(obcParams1.x, fabs(r-obcParams2.y));
    float l_ijI = 1.0f/max(obcParams2.x, fabs(r-obcParams1.y));
    float u_ijJ = 1.0f/rScaledRadiusJ;
    float u_ijI = 1.0f/rScaledRadiusI;
    float l_ij2J = l_ijJ*l_ijJ;
    float l_ij2I = l_ijI*l_ijI;
    float u_ij2J = u_ijJ*u_ijJ;
    float u_ij2I = u_ijI*u_ijI;
    float t1J = log(u_ijJ/l_ijJ);
    float t1I = log(u_ijI/l_ijI);
    float t2J = (l_ij2J-u_ij2J);
    float t2I = (l_ij2I-u_ij2I);
    float t3J = t2J*invR;
    float t3I = t2I*invR;
    t1J *= invR;
    t1I *= invR;
    if (obcParams1.x < rScaledRadiusJ) {
        float term = 0.125f*(1.0f+obcParams2.y*obcParams2.y*invRSquared)*t3J + 0.25f*t1J*invRSquared;
        dEdR += bornForce1*term;
    }
    if (obcParams2.x < rScaledRadiusJ) {
        float term = 0.125f*(1.0f+obcParams1.y*obcParams1.y*invRSquared)*t3I + 0.25f*t1I*invRSquared;
        dEdR += bornForce2*term;
    }
}
