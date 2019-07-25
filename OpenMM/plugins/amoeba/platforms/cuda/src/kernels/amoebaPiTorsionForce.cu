// compute the value of the pi-orbital torsion angle

real3 ad = make_real3(pos1.x-pos4.x, pos1.y-pos4.y, pos1.z-pos4.z);
real3 bd = make_real3(pos2.x-pos4.x, pos2.y-pos4.y, pos2.z-pos4.z);
real3 ec = make_real3(pos5.x-pos3.x, pos5.y-pos3.y, pos5.z-pos3.z);
real3 gc = make_real3(pos6.x-pos3.x, pos6.y-pos3.y, pos6.z-pos3.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ad)
APPLY_PERIODIC_TO_DELTA(bd)
APPLY_PERIODIC_TO_DELTA(ec)
APPLY_PERIODIC_TO_DELTA(gc)
#endif

real xip = ad.y*bd.z - bd.y*ad.z + pos3.x;
real yip = ad.z*bd.x - bd.z*ad.x + pos3.y;
real zip = ad.x*bd.y - bd.x*ad.y + pos3.z;

real xiq = ec.y*gc.z - gc.y*ec.z + pos4.x;
real yiq = ec.z*gc.x - gc.z*ec.x + pos4.y;
real ziq = ec.x*gc.y - gc.x*ec.y + pos4.z;

real xcp = pos3.x - xip;
real ycp = pos3.y - yip;
real zcp = pos3.z - zip;

real xdc = pos4.x - pos3.x;
real ydc = pos4.y - pos3.y;
real zdc = pos4.z - pos3.z;

real xqd = xiq - pos4.x;
real yqd = yiq - pos4.y;
real zqd = ziq - pos4.z;

real xt = ycp*zdc - ydc*zcp;
real yt = zcp*xdc - zdc*xcp;
real zt = xcp*ydc - xdc*ycp;

real xu = ydc*zqd - yqd*zdc;
real yu = zdc*xqd - zqd*xdc;
real zu = xdc*yqd - xqd*ydc;

real xtu = yt*zu - yu*zt;
real ytu = zt*xu - zu*xt;
real ztu = xt*yu - xu*yt;

real rt2 = xt*xt + yt*yt + zt*zt;
real ru2 = xu*xu + yu*yu + zu*zu;

real rtru = sqrtf(rt2 * ru2);
real rdc = sqrtf(xdc*xdc + ydc*ydc + zdc*zdc);

real cosine = rtru > 0.0f ? (xt*xu + yt*yu + zt*zu) / rtru : 0.0f;
real sine = (rtru*rdc) > 0.0f ? (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru) : 0.0f;

// zero energy/force if rtru == 0

float v2 = PARAMS[index];
v2 = (rtru > 0 ? v2 : 0.0f);

// compute the multiple angle trigonometry and the phase terms

real cosine2 = cosine*cosine - sine*sine;
real sine2 = 2.0f * cosine * sine;
real phi2 = 1.0f - cosine2;
real dphi2 = 2.0f * sine2;

// calculate pi-orbital torsion energy and master chain rule term

energy += v2 * phi2;
real dedphi = v2 * dphi2;

// chain rule terms for first derivative components

real xdp = pos4.x - xip;
real ydp = pos4.y - yip;
real zdp = pos4.z - zip;

real xqc = xiq - pos3.x;
real yqc = yiq - pos3.y;
real zqc = ziq - pos3.z;

real dedxt = dedphi * (yt*zdc - ydc*zt) / (rt2*rdc);
real dedyt = dedphi * (zt*xdc - zdc*xt) / (rt2*rdc);
real dedzt = dedphi * (xt*ydc - xdc*yt) / (rt2*rdc);

real dedxu = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc);
real dedyu = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc);
real dedzu = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc);

// compute first derivative components for pi-orbital angle

real dedxip = zdc*dedyt - ydc*dedzt;
real dedyip = xdc*dedzt - zdc*dedxt;
real dedzip = ydc*dedxt - xdc*dedyt;

real dedxic = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu;
real dedyic = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu;
real dedzic = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu;

real dedxid = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu;
real dedyid = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu;
real dedzid = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu;

real dedxiq = zdc*dedyu - ydc*dedzu;
real dedyiq = xdc*dedzu - zdc*dedxu;
real dedziq = ydc*dedxu - xdc*dedyu;

// compute first derivative components for individual atoms

real dedxia = bd.y*dedzip - bd.z*dedyip;
real dedyia = bd.z*dedxip - bd.x*dedzip;
real dedzia = bd.x*dedyip - bd.y*dedxip;

real dedxib = ad.z*dedyip - ad.y*dedzip;
real dedyib = ad.x*dedzip - ad.z*dedxip;
real dedzib = ad.y*dedxip - ad.x*dedyip;

real dedxie = gc.y*dedziq - gc.z*dedyiq;
real dedyie = gc.z*dedxiq - gc.x*dedziq;
real dedzie = gc.x*dedyiq - gc.y*dedxiq;

real dedxig = ec.z*dedyiq - ec.y*dedziq;
real dedyig = ec.x*dedziq - ec.z*dedxiq;
real dedzig = ec.y*dedxiq - ec.x*dedyiq;

dedxic = dedxic + dedxip - dedxie - dedxig;
dedyic = dedyic + dedyip - dedyie - dedyig;
dedzic = dedzic + dedzip - dedzie - dedzig;
dedxid = dedxid + dedxiq - dedxia - dedxib;
dedyid = dedyid + dedyiq - dedyia - dedyib;
dedzid = dedzid + dedziq - dedzia - dedzib;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib, -dedyib, -dedzib);
real3 force3 = make_real3(-dedxic, -dedyic, -dedzic);
real3 force4 = make_real3(-dedxid, -dedyid, -dedzid);
real3 force5 = make_real3(-dedxie, -dedyie, -dedzie);
real3 force6 = make_real3(-dedxig, -dedyig, -dedzig);