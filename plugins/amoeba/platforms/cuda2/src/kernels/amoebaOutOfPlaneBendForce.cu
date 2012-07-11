// compute the value of the bond angle

real xab = pos1.x - pos2.x;
real yab = pos1.y - pos2.y;
real zab = pos1.z - pos2.z;

real xcb = pos3.x - pos2.x;
real ycb = pos3.y - pos2.y;
real zcb = pos3.z - pos2.z;

// compute the out-of-plane bending angle

real xdb = pos4.x - pos2.x;
real ydb = pos4.y - pos2.y;
real zdb = pos4.z - pos2.z;

real xad = pos1.x - pos4.x;
real yad = pos1.y - pos4.y;
real zad = pos1.z - pos4.z;

real xcd = pos3.x - pos4.x;
real ycd = pos3.y - pos4.y;
real zcd = pos3.z - pos4.z;

real rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
real rad2 = xad*xad + yad*yad + zad*zad;
real rcd2 = xcd*xcd + ycd*ycd + zcd*zcd;

real ee = xab*(ycb*zdb-zcb*ydb) + yab*(zcb*xdb-xcb*zdb) + zab*(xcb*ydb-ycb*xdb);

real dot = xad*xcd + yad*ycd + zad*zcd;
real cc = rad2*rcd2 - dot*dot;
real bkk2 = (cc != 0 ? (ee*ee)/(cc) : (real) 0);
bkk2 = rdb2 - bkk2;

real adXcd_0 = yad*zcd - zad*ycd;
real adXcd_1 = zad*xcd - xad*zcd;
real adXcd_2 = xad*ycd - yad*xcd;
real adXcd_nrm2 = adXcd_0*adXcd_0 + adXcd_1*adXcd_1 + adXcd_2*adXcd_2;

real adXcd_dot_db = xdb*adXcd_0 + ydb*adXcd_1 + zdb*adXcd_2;
adXcd_dot_db /= SQRT(rdb2*adXcd_nrm2);

real angle = abs(ASIN(adXcd_dot_db));

// find the out-of-plane energy and master chain rule terms

real dt = RAD_TO_DEG*angle;
real dt2 = dt * dt;
real dt3 = dt2 * dt;
real dt4 = dt2 * dt2;
float k = (rdb2 != 0 && cc != 0) ? PARAMS[index] : 0.0f;

energy += k*dt2*(1.0f + CUBIC_K*dt + QUARTIC_K*dt2 + PENTIC_K*dt3 + SEXTIC_K*dt4);

real deddt = k*dt*RAD_TO_DEG*(2.0f + 3.0f*CUBIC_K*dt + 4.0f*QUARTIC_K*dt2 + 5.0f*PENTIC_K*dt3 + 6.0f*SEXTIC_K*dt4);

real eeSign = (ee >= 0 ? 1 : -1);
real dedcos = -deddt*eeSign/SQRT(cc*bkk2);

// chain rule terms for first derivative components

real term = ee / cc;

real dccdxia = (xad*rcd2-xcd*dot) * term;
real dccdyia = (yad*rcd2-ycd*dot) * term;
real dccdzia = (zad*rcd2-zcd*dot) * term;

real dccdxic = (xcd*rad2-xad*dot) * term;
real dccdyic = (ycd*rad2-yad*dot) * term;
real dccdzic = (zcd*rad2-zad*dot) * term;

real dccdxid = -dccdxia - dccdxic;
real dccdyid = -dccdyia - dccdyic;
real dccdzid = -dccdzia - dccdzic;

term = ee / rdb2;

real deedxia = ydb*zcb - zdb*ycb;
real deedyia = zdb*xcb - xdb*zcb;
real deedzia = xdb*ycb - ydb*xcb;

real deedxic = yab*zdb - zab*ydb;
real deedyic = zab*xdb - xab*zdb;
real deedzic = xab*ydb - yab*xdb;

real deedxid = ycb*zab - zcb*yab + xdb*term;
real deedyid = zcb*xab - xcb*zab + ydb*term;
real deedzid = xcb*yab - ycb*xab + zdb*term;

// compute first derivative components for this angle

real3 force1 = make_real3(-dedcos*(dccdxia+deedxia), -dedcos*(dccdyia+deedyia), -dedcos*(dccdzia+deedzia));
real3 force3 = make_real3(-dedcos*(dccdxic+deedxic), -dedcos*(dccdyic+deedyic), -dedcos*(dccdzic+deedzic));
real3 force4 = make_real3(-dedcos*(dccdxid+deedxid), -dedcos*(dccdyid+deedyid), -dedcos*(dccdzid+deedzid));
real3 force2 = make_real3(-force1.x-force3.x-force4.x, -force1.y-force3.y-force4.y, -force1.z-force3.z-force4.z);
