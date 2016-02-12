// compute the value of the bond angle

real xab = pos1.x - pos2.x;
real yab = pos1.y - pos2.y;
real zab = pos1.z - pos2.z;

real xcb = pos3.x - pos2.x;
real ycb = pos3.y - pos2.y;
real zcb = pos3.z - pos2.z;

real rab = SQRT(xab*xab + yab*yab + zab*zab);
real rcb = SQRT(xcb*xcb + ycb*ycb + zcb*zcb);

real xp = ycb*zab - zcb*yab;
real yp = zcb*xab - xcb*zab;
real zp = xcb*yab - ycb*xab;

real rp = SQRT(xp*xp + yp*yp + zp*zp);

real dotp = xab*xcb + yab*ycb + zab*zcb;
real cosine = rab*rcb > 0 ? (dotp / (rab*rcb)) : (real) 1;
cosine = (cosine > 1 ? (real) 1 : cosine);
cosine = (cosine < -1 ? -(real) 1 : cosine);
real angle;
if (cosine > 0.99f || cosine < -0.99f) {
    // Highly unlikely a stretch-bend angle will be near 0 or 180, but just in case...
    real3 cross_prod = cross(make_real3(xab, yab, zab), make_real3(xcb, ycb, zcb));
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

real ddtdxia = terma * (yab*zp-zab*yp);
real ddtdyia = terma * (zab*xp-xab*zp);
real ddtdzia = terma * (xab*yp-yab*xp);

real ddtdxic = termc * (ycb*zp-zcb*yp);
real ddtdyic = termc * (zcb*xp-xcb*zp);
real ddtdzic = termc * (xcb*yp-ycb*xp);

// find chain rule terms for the bond length deviations

real dr1 = (parameters.x > 0 ? (rab - parameters.x) : (real) 0);
terma = (parameters.x > 0 ? RECIP(rab) : (real) 0);

real dr2 = (parameters.y > 0 ? (rcb - parameters.y) : (real) 0);
termc = (parameters.y > 0 ? RECIP(rcb) : (real) 0);
real frc1 = ((rp != 0) ? force_constants.x : (real) 0);
real frc2 = ((rp != 0) ? force_constants.y : (real) 0);

real drkk = dr1*frc1 + dr2*frc2;

real ddrdxia = terma * xab;
real ddrdyia = terma * yab;
real ddrdzia = terma * zab;

real ddrdxic = termc * xcb;
real ddrdyic = termc * ycb;
real ddrdzic = termc * zcb;

// get the energy and master chain rule terms for derivatives

energy += dt*drkk;

real3 force1 = make_real3(-frc1*dt*ddrdxia-ddtdxia*drkk, -frc1*dt*ddrdyia-ddtdyia*drkk, -frc1*dt*ddrdzia-ddtdzia*drkk);
real3 force3 = make_real3(-frc2*dt*ddrdxic-ddtdxic*drkk, -frc2*dt*ddrdyic-ddtdyic*drkk, -frc2*dt*ddrdzic-ddtdzic*drkk);
real3 force2 = make_real3(-force1.x-force3.x, -force1.y-force3.y, -force1.z-force3.z);
