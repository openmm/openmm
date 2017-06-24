// compute the value of the angle-torsion energy & force

real xba = pos2.x - pos1.x;
real yba = pos2.y - pos1.y;
real zba = pos2.z - pos1.z;

real xcb = pos3.x - pos2.x;
real ycb = pos3.y - pos2.y;
real zcb = pos3.z - pos2.z;

real xdc = pos4.x - pos3.x;
real ydc = pos4.y - pos3.y;
real zdc = pos4.z - pos3.z;

real rba2 = xba*xba + yba*yba + zba*zba;
real rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
real rdc2 = xdc*xdc + ydc*ydc + zdc*zdc;

real xt = yba*zcb - ycb*zba;
real yt = zba*xcb - zcb*xba;
real zt = xba*ycb - xcb*yba;

real xu = ycb*zdc - ydc*zcb;
real yu = zcb*xdc - zdc*xcb;
real zu = xcb*ydc - xdc*ycb;

real xtu = yt*zu - yu*zt;
real ytu = zt*xu - zu*xt;
real ztu = xt*yu - xu*yt;

real rt2 = xt*xt + yt*yt + zt*zt;
real ru2 = xu*xu + yu*yu + zu*zu;
real rtru = SQRT(rt2 * ru2);

real xca = pos3.x - pos1.x;
real yca = pos3.y - pos1.y;
real zca = pos3.z - pos1.z;

real xdb = pos4.x - pos2.x;
real ydb = pos4.y - pos2.y;
real zdb = pos4.z - pos2.z;

real rcb = SQRT(rcb2);
real cosine = rtru > 0.0f ? (xt*xu + yt*yu + zt*zu) / rtru : 0.0f;
real sine = rtru > 0.0f ? (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru) : 0.0f;

real c1 = 1.0f;
real s1 = 0.0f;
real c2 = -1.0f;
real s2 = 0.0f;
real c3 = 1.0f;
real s3 = 0.0f;

real cosine2 = cosine*cosine - sine*sine;
real sine2 = 2.0f * cosine * sine;
real cosine3 = cosine*cosine2 - sine*sine2;
real sine3 = cosine*sine2 + sine*cosine2;

real phi1 = rtru > 0.0f ? 1.0f + (cosine*c1 + sine*s1) : 0.0f;
real phi2 = rtru > 0.0f ? 1.0f + (cosine2*c2 + sine2*s2) : 0.0f;
real phi3 = rtru > 0.0f ? 1.0f + (cosine3*c3 + sine3*s3) : 0.0f;

real dphi1 = rtru > 0.0f ? cosine*s1 - sine*c1 : 0.0f;
real dphi2 = rtru > 0.0f ? 2.0f * (cosine2*s2 - sine2*c2) : 0.0f;
real dphi3 = rtru > 0.0f ? 3.0f * (cosine3*s3 - sine3*c3) : 0.0f;

float2 parameters = PARAMS[index];
float3 amplitude1 = FORCE_CONSTANTS_FIRST[index];
float3 amplitude2 = FORCE_CONSTANTS_SECOND[index];

real dot = xba*xcb + yba*ycb + zba*zcb;
real cosang = -dot / SQRT(rba2*rcb2);
real angle = RAD_TO_DEG * ACOS(cosang);
real dt = angle - parameters.x;
real e1 = dt * (amplitude1.x*phi1 + amplitude1.y*phi2 + amplitude1.z*phi3);
real dedphi = dt * (amplitude1.x*dphi1 + amplitude1.y*dphi2 + amplitude1.z*dphi3);
real ddt = RAD_TO_DEG * (amplitude1.x*phi1 + amplitude1.y*phi2 + amplitude1.z*phi3);

energy += e1;

real dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb);
real dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb);
real dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb);
real dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb);
real dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb);
real dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb);

real terma = -ddt / (rba2*SQRT(rt2));
real termc = ddt / (rcb2*SQRT(rt2));

real dedxia = terma*(zba*yt-yba*zt) + zcb*dedyt - ycb*dedzt;
real dedyia = terma*(xba*zt-zba*xt) + xcb*dedzt - zcb*dedxt;
real dedzia = terma*(yba*xt-xba*yt) + ycb*dedxt - xcb*dedyt;
real dedxib = terma*(yba*zt-zba*yt) + termc*(zcb*yt-ycb*zt)
	+ yca*dedzt - zca*dedyt
	+ zdc*dedyu - ydc*dedzu;
real dedyib = terma*(zba*xt-xba*zt) + termc*(xcb*zt-zcb*xt)
	+ zca*dedxt - xca*dedzt
	+ xdc*dedzu - zdc*dedxu;
real dedzib = terma*(xba*yt-yba*xt) + termc*(ycb*xt-xcb*yt)
	+ xca*dedyt - yca*dedxt
	+ ydc*dedxu - xdc*dedyu;
real dedxic = termc*(ycb*zt-zcb*yt) + zba*dedyt
	- yba*dedzt + ydb*dedzu - zdb*dedyu;
real dedyic = termc*(zcb*xt-xcb*zt) + xba*dedzt
	- zba*dedxt + zdb*dedxu - xdb*dedzu;
real dedzic = termc*(xcb*yt-ycb*xt) + yba*dedxt
	- xba*dedyt + xdb*dedyu - ydb*dedxu;
real dedxid = zcb*dedyu - ycb*dedzu;
real dedyid = xcb*dedzu - zcb*dedxu;
real dedzid = ycb*dedxu - xcb*dedyu;

dot = dot = xcb*xdc + ycb*ydc + zcb*zdc;
cosang = -dot / SQRT(rcb2*rdc2);
angle = RAD_TO_DEG * ACOS(cosang);
dt = angle - parameters.y;
real e2 = dt * (amplitude2.x*phi1 + amplitude2.y*phi2 + amplitude2.z*phi3);
dedphi = dt * (amplitude2.x*dphi1 + amplitude2.y*dphi2 + amplitude2.z*dphi3);
ddt = RAD_TO_DEG * (amplitude2.x*phi1 + amplitude2.y*phi2 + amplitude2.z*phi3);

energy += e2;

dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb);
dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb);
dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb);
dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb);
dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb);
dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb);

real termb = -ddt / (rcb2*SQRT(ru2));
real termd = ddt / (rdc2*sqrt(ru2));

dedxia += zcb*dedyt - ycb*dedzt;
dedyia += xcb*dedzt - zcb*dedxt;
dedzia += ycb*dedxt - xcb*dedyt;
dedxib += termb*(zcb*yu-ycb*zu) + yca*dedzt
	- zca*dedyt + zdc*dedyu - ydc*dedzu;
dedyib += termb*(xcb*zu-zcb*xu) + zca*dedxt
	- xca*dedzt + xdc*dedzu - zdc*dedxu;
dedzib += termb*(ycb*xu-xcb*yu) + xca*dedyt
	- yca*dedxt + ydc*dedxu - xdc*dedyu;
dedxic += termb*(ycb*zu-zcb*yu)
	+ termd*(zdc*yu-ydc*zu) + zba*dedyt
	- yba*dedzt + ydb*dedzu - zdb*dedyu;
dedyic += termb*(zcb*xu-xcb*zu)
	+ termd*(xdc*zu-zdc*xu) + xba*dedzt
	- zba*dedxt + zdb*dedxu - xdb*dedzu;
dedzic += termb*(xcb*yu-ycb*xu)
	+ termd*(ydc*xu-xdc*yu) + yba*dedxt
	- xba*dedyt + xdb*dedyu - ydb*dedxu;
dedxid += termd*(ydc*zu-zdc*yu)
	+ zcb*dedyu - ycb*dedzu;
dedyid += termd*(zdc*xu-xdc*zu)
	+ xcb*dedzu - zcb*dedxu;
dedzid += termd*(xdc*yu-ydc*xu)
	+ ycb*dedxu - xcb*dedyu;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib, -dedyib, -dedzib);
real3 force3 = make_real3(-dedxic, -dedyic, -dedzic);
real3 force4 = make_real3(-dedxid, -dedyid, -dedzid);
