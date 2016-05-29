// compute the value of the bond angle

real3 ab = make_real3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z);
real3 cb = make_real3(pos3.x-pos2.x, pos3.y-pos2.y, pos3.z-pos2.z);
real3 db = make_real3(pos4.x-pos2.x, pos4.y-pos2.y, pos4.z-pos2.z);
real3 ad = make_real3(pos1.x-pos4.x, pos1.y-pos4.y, pos1.z-pos4.z);
real3 cd = make_real3(pos3.x-pos4.x, pos3.y-pos4.y, pos3.z-pos4.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ab)
APPLY_PERIODIC_TO_DELTA(cb)
APPLY_PERIODIC_TO_DELTA(db)
APPLY_PERIODIC_TO_DELTA(ad)
APPLY_PERIODIC_TO_DELTA(cd)
#endif

real rdb2 = db.x*db.x + db.y*db.y + db.z*db.z;
real rad2 = ad.x*ad.x + ad.y*ad.y + ad.z*ad.z;
real rcd2 = cd.x*cd.x + cd.y*cd.y + cd.z*cd.z;

real ee = ab.x*(cb.y*db.z-cb.z*db.y) + ab.y*(cb.z*db.x-cb.x*db.z) + ab.z*(cb.x*db.y-cb.y*db.x);

real dot = ad.x*cd.x + ad.y*cd.y + ad.z*cd.z;
real cc = rad2*rcd2 - dot*dot;
real bkk2 = (cc != 0 ? (ee*ee)/(cc) : (real) 0);
bkk2 = rdb2 - bkk2;

real adXcd_0 = ad.y*cd.z - ad.z*cd.y;
real adXcd_1 = ad.z*cd.x - ad.x*cd.z;
real adXcd_2 = ad.x*cd.y - ad.y*cd.x;
real adXcd_nrm2 = adXcd_0*adXcd_0 + adXcd_1*adXcd_1 + adXcd_2*adXcd_2;

real adXcd_dot_db = db.x*adXcd_0 + db.y*adXcd_1 + db.z*adXcd_2;
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

real dccdxia = (ad.x*rcd2-cd.x*dot) * term;
real dccdyia = (ad.y*rcd2-cd.y*dot) * term;
real dccdzia = (ad.z*rcd2-cd.z*dot) * term;

real dccdxic = (cd.x*rad2-ad.x*dot) * term;
real dccdyic = (cd.y*rad2-ad.y*dot) * term;
real dccdzic = (cd.z*rad2-ad.z*dot) * term;

real dccdxid = -dccdxia - dccdxic;
real dccdyid = -dccdyia - dccdyic;
real dccdzid = -dccdzia - dccdzic;

term = ee / rdb2;

real deedxia = db.y*cb.z - db.z*cb.y;
real deedyia = db.z*cb.x - db.x*cb.z;
real deedzia = db.x*cb.y - db.y*cb.x;

real deedxic = ab.y*db.z - ab.z*db.y;
real deedyic = ab.z*db.x - ab.x*db.z;
real deedzic = ab.x*db.y - ab.y*db.x;

real deedxid = cb.y*ab.z - cb.z*ab.y + db.x*term;
real deedyid = cb.z*ab.x - cb.x*ab.z + db.y*term;
real deedzid = cb.x*ab.y - cb.y*ab.x + db.z*term;

// compute first derivative components for this angle

real3 force1 = make_real3(-dedcos*(dccdxia+deedxia), -dedcos*(dccdyia+deedyia), -dedcos*(dccdzia+deedzia));
real3 force3 = make_real3(-dedcos*(dccdxic+deedxic), -dedcos*(dccdyic+deedyic), -dedcos*(dccdzic+deedzic));
real3 force4 = make_real3(-dedcos*(dccdxid+deedxid), -dedcos*(dccdyid+deedyid), -dedcos*(dccdzid+deedzid));
real3 force2 = make_real3(-force1.x-force3.x-force4.x, -force1.y-force3.y-force4.y, -force1.z-force3.z-force4.z);
