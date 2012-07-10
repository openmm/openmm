// compute the value of the pi-orbital torsion angle

real xad = pos1.x - pos4.x;
real yad = pos1.y - pos4.y;
real zad = pos1.z - pos4.z;

real xbd = pos2.x - pos4.x;
real ybd = pos2.y - pos4.y;
real zbd = pos2.z - pos4.z;

real xec = pos5.x - pos3.x;
real yec = pos5.y - pos3.y;
real zec = pos5.z - pos3.z;

real xgc = pos6.x - pos3.x;
real ygc = pos6.y - pos3.y;
real zgc = pos6.z - pos3.z;

real xip = yad*zbd - ybd*zad + pos3.x;
real yip = zad*xbd - zbd*xad + pos3.y;
real zip = xad*ybd - xbd*yad + pos3.z;

real xiq = yec*zgc - ygc*zec + pos4.x;
real yiq = zec*xgc - zgc*xec + pos4.y;
real ziq = xec*ygc - xgc*yec + pos4.z;

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

real dedxia = ybd*dedzip - zbd*dedyip;
real dedyia = zbd*dedxip - xbd*dedzip;
real dedzia = xbd*dedyip - ybd*dedxip;

real dedxib = zad*dedyip - yad*dedzip;
real dedyib = xad*dedzip - zad*dedxip;
real dedzib = yad*dedxip - xad*dedyip;

real dedxie = ygc*dedziq - zgc*dedyiq;
real dedyie = zgc*dedxiq - xgc*dedziq;
real dedzie = xgc*dedyiq - ygc*dedxiq;

real dedxig = zec*dedyiq - yec*dedziq;
real dedyig = xec*dedziq - zec*dedxiq;
real dedzig = yec*dedxiq - xec*dedyiq;

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