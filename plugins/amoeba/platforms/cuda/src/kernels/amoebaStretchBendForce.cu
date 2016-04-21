// compute the value of the bond angle

real3 ab = make_real3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z);
real3 cb = make_real3(pos3.x-pos2.x, pos3.y-pos2.y, pos3.z-pos2.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ab)
APPLY_PERIODIC_TO_DELTA(cb)
#endif

real rab = SQRT(ab.x*ab.x + ab.y*ab.y + ab.z*ab.z);
real rcb = SQRT(cb.x*cb.x + cb.y*cb.y + cb.z*cb.z);

real xp = cb.y*ab.z - cb.z*ab.y;
real yp = cb.z*ab.x - cb.x*ab.z;
real zp = cb.x*ab.y - cb.y*ab.x;

real rp = SQRT(xp*xp + yp*yp + zp*zp);

real dotp = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;
real cosine = rab*rcb > 0 ? (dotp / (rab*rcb)) : (real) 1;
cosine = (cosine > 1 ? (real) 1 : cosine);
cosine = (cosine < -1 ? -(real) 1 : cosine);
real angle;
if (cosine > 0.99f || cosine < -0.99f) {
    // Highly unlikely a stretch-bend angle will be near 0 or 180, but just in case...
    real3 cross_prod = cross(make_real3(ab.x, ab.y, ab.z), make_real3(cb.x, cb.y, cb.z));
    angle = ASIN(SQRT(dot(cross_prod, cross_prod))/(rab*rcb))*RAD_TO_DEG;
    if (cosine < 0.0f)
        angle = 180-angle;
}
else
    angle = ACOS(cosine)*RAD_TO_DEG;

// find chain rule terms for the bond angle deviation

float3 parameters = PARAMS[index];
float2 force_constants = FORCE_CONSTANTS[index];

real dt = angle - RAD_TO_DEG*parameters.z;
real terma = rab*rp != 0 ? (-RAD_TO_DEG/(rab*rab*rp)) : (real) 0;
real termc = rcb*rp != 0 ? (RAD_TO_DEG/(rcb*rcb*rp)) : (real) 0;

real ddtdxia = terma * (ab.y*zp-ab.z*yp);
real ddtdyia = terma * (ab.z*xp-ab.x*zp);
real ddtdzia = terma * (ab.x*yp-ab.y*xp);

real ddtdxic = termc * (cb.y*zp-cb.z*yp);
real ddtdyic = termc * (cb.z*xp-cb.x*zp);
real ddtdzic = termc * (cb.x*yp-cb.y*xp);

// find chain rule terms for the bond length deviations

real dr1 = (parameters.x > 0 ? (rab - parameters.x) : (real) 0);
terma = (parameters.x > 0 ? RECIP(rab) : (real) 0);

real dr2 = (parameters.y > 0 ? (rcb - parameters.y) : (real) 0);
termc = (parameters.y > 0 ? RECIP(rcb) : (real) 0);
real frc1 = ((rp != 0) ? force_constants.x : (real) 0);
real frc2 = ((rp != 0) ? force_constants.y : (real) 0);

real drkk = dr1*frc1 + dr2*frc2;

real ddrdxia = terma * ab.x;
real ddrdyia = terma * ab.y;
real ddrdzia = terma * ab.z;

real ddrdxic = termc * cb.x;
real ddrdyic = termc * cb.y;
real ddrdzic = termc * cb.z;

// get the energy and master chain rule terms for derivatives

energy += dt*drkk;

real3 force1 = make_real3(-frc1*dt*ddrdxia-ddtdxia*drkk, -frc1*dt*ddrdyia-ddtdyia*drkk, -frc1*dt*ddrdzia-ddtdzia*drkk);
real3 force3 = make_real3(-frc2*dt*ddrdxic-ddtdxic*drkk, -frc2*dt*ddrdyic-ddtdyic*drkk, -frc2*dt*ddrdzic-ddtdzic*drkk);
real3 force2 = make_real3(-force1.x-force3.x, -force1.y-force3.y, -force1.z-force3.z);
