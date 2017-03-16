real4 delta = pos2-pos1;
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(delta)
#endif
real r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
COMPUTE_FORCE
dEdR = (r > 0.0f) ? (dEdR / r) : 0.0f;
delta.xyz *= dEdR;
real4 force1 = delta;
real4 force2 = -delta;