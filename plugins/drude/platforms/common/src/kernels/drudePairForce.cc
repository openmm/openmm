float2 drudeParams = PARAMS[index];
real3 force1 = make_real3(0);
real3 force2 = make_real3(0);
real3 force3 = make_real3(0);
real3 force4 = make_real3(0);

// First pair.

real3 delta = make_real3(pos1.x-pos3.x, pos1.y-pos3.y, pos1.z-pos3.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(delta)
#endif
real rInv = RSQRT(dot(delta, delta));
real r = RECIP(rInv);
real u = drudeParams.x*r;
real screening = 1-(1+0.5f*u)*EXP(-u);
real pairEnergy = drudeParams.y*screening*rInv;
energy += pairEnergy;
real3 f = delta*(drudeParams.y*rInv*rInv*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x));
force1 += f;
force3 -= f;

// Second pair.

delta = make_real3(pos1.x-pos4.x, pos1.y-pos4.y, pos1.z-pos4.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(delta)
#endif
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = -drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(-drudeParams.y*rInv*rInv*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x));
force1 += f;
force4 -= f;

// Third pair.

delta = make_real3(pos2.x-pos3.x, pos2.y-pos3.y, pos2.z-pos3.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(delta)
#endif
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = -drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(-drudeParams.y*rInv*rInv*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x));
force2 += f;
force3 -= f;

// Fourth pair.

delta = make_real3(pos2.x-pos4.x, pos2.y-pos4.y, pos2.z-pos4.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(delta)
#endif
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(drudeParams.y*rInv*rInv*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x));
force2 += f;
force4 -= f;
