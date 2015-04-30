float2 drudeParams = PARAMS[index];
real4 force1 = 0;
real4 force2 = 0;
real4 force3 = 0;
real4 force4 = 0;

// First pair.

real4 delta = (real4) (pos1.xyz-pos3.xyz, 0);
real rInv = RSQRT(dot(delta, delta));
real r = RECIP(rInv);
real u = drudeParams.x*r;
real screening = 1-(1+0.5f*u)*EXP(-u);
real pairEnergy = drudeParams.y*screening*rInv;
energy += pairEnergy;
real4 f = delta*(drudeParams.y*rInv*rInv)*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x);
force1 += f;
force3 -= f;

// Second pair.

delta = (real4) (pos1.xyz-pos4.xyz, 0);
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = -drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(-drudeParams.y*rInv*rInv)*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x);
force1 += f;
force4 -= f;

// Third pair.

delta = (real4) (pos2.xyz-pos3.xyz, 0);
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = -drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(-drudeParams.y*rInv*rInv)*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x);
force2 += f;
force3 -= f;

// Fourth pair.

delta = (real4) (pos2.xyz-pos4.xyz, 0);
rInv = RSQRT(dot(delta, delta));
r = RECIP(rInv);
u = drudeParams.x*r;
screening = 1-(1+0.5f*u)*EXP(-u);
pairEnergy = drudeParams.y*screening*rInv;
energy += pairEnergy;
f = delta*(drudeParams.y*rInv*rInv)*(screening*rInv-0.5f*(1+u)*EXP(-u)*drudeParams.x);
force2 += f;
force4 -= f;
