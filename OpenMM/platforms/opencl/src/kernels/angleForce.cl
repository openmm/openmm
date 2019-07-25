real4 v0 = pos2-pos1;
real4 v1 = pos2-pos3;
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(v0)
APPLY_PERIODIC_TO_DELTA(v1)
#endif
real4 cp = cross(v0, v1);
real rp = cp.x*cp.x + cp.y*cp.y + cp.z*cp.z;
rp = max(SQRT(rp), (real) 1.0e-06f);
real r21 = v0.x*v0.x + v0.y*v0.y + v0.z*v0.z;
real r23 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
real dot = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
real cosine = clamp(dot*RSQRT(r21*r23), (real) -1, (real) 1);
real theta = acos(cosine);
COMPUTE_FORCE
real4 force1 = cross(v0, cp)*(dEdAngle/(r21*rp));
real4 force3 = cross(cp, v1)*(dEdAngle/(r23*rp));
real4 force2 = -force1-force3;
