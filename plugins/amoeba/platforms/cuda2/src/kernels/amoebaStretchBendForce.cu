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

real dot = xab*xcb + yab*ycb + zab*zcb;
real cosine = rab*rcb > 0 ? (dot / (rab*rcb)) : (real) 1;
cosine = (cosine > 1 ? (real) 1 : cosine);
cosine = (cosine < -1 ? -(real) 1 : cosine);
real angle = ACOS(cosine);

// find chain rule terms for the bond angle deviation

float4 parameters = PARAMS[index];

real dt = RAD_TO_DEG*(angle - parameters.z);
real terma = rab*rp != 0 ? (-RAD_TO_DEG/(rab*rab*rp)) : (real) 0;
real termc = rcb*rp != 0 ? (RAD_TO_DEG/(rcb*rcb*rp)) : (real) 0;

real ddtdxia = terma * (yab*zp-zab*yp);
real ddtdyia = terma * (zab*xp-xab*zp);
real ddtdzia = terma * (xab*yp-yab*xp);

real ddtdxic = termc * (ycb*zp-zcb*yp);
real ddtdyic = termc * (zcb*xp-xcb*zp);
real ddtdzic = termc * (xcb*yp-ycb*xp);

// find chain rule terms for the bond length deviations

real dr = (parameters.x > 0 ? (rab - parameters.x) : (real) 0);
terma = (parameters.x > 0 ? RECIP(rab) : (real) 0);

dr += (parameters.y > 0 ? (rcb - parameters.y) : (real) 0);
termc = (parameters.y > 0 ? RECIP(rcb) : (real) 0);

real ddrdxia = terma * xab;
real ddrdyia = terma * yab;
real ddrdzia = terma * zab;

real ddrdxic = termc * xcb;
real ddrdyic = termc * ycb;
real ddrdzic = termc * zcb;

// get the energy and master chain rule terms for derivatives

real term = ((rp != 0) ? parameters.w : (real) 0);
energy += term*dt*dr;

real3 force1 = make_real3(-term*(dt*ddrdxia+ddtdxia*dr), -term*(dt*ddrdyia+ddtdyia*dr), -term*(dt*ddrdzia+ddtdzia*dr));
real3 force3 = make_real3(-term*(dt*ddrdxic+ddtdxic*dr), -term*(dt*ddrdyic+ddtdyic*dr), -term*(dt*ddrdzic+ddtdzic*dr));
real3 force2 = make_real3(-force1.x-force3.x, -force1.y-force3.y, -force1.z-force3.z);
