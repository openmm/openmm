// compute the value of the strech-torsion energy & force

real xba = pos2.x - pos1.x;
real yba = pos2.y - pos1.y;
real zba = pos2.z - pos1.z;

real xcb = pos3.x - pos2.x;
real ycb = pos3.y - pos2.y;
real zcb = pos3.z - pos2.z;

real xdc = pos4.x - pos3.x;
real ydc = pos4.y - pos3.y;
real zdc = pos4.z - pos3.z;

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

real rba = SQRT(xba*xba + yba*yba + zba*zba);
real rcb = SQRT(xcb*xcb + ycb*ycb + zcb*zcb);
real rdc = SQRT(xdc*xdc + ydc*ydc + zdc*zdc);

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

float3 parameters = PARAMS[index];
float3 amplitude1 = FORCE_CONSTANTS_FIRST[index];
float3 amplitude2 = FORCE_CONSTANTS_SECOND[index];
float3 amplitude3 = FORCE_CONSTANTS_THIRD[index];

real dr = rba - parameters.x;
real e1 = dr * (amplitude1.x*phi1 + amplitude1.y*phi2 + amplitude1.z*phi3);
real dedphi = dr * (amplitude1.x*dphi1 + amplitude1.y*dphi2 + amplitude1.z*dphi3);
real ddr = (amplitude1.x*phi1 + amplitude1.y*phi2 + amplitude1.z*phi3) / rba;

energy += e1;

real ddrdx = xba * ddr;
real ddrdy = yba * ddr;
real ddrdz = zba * ddr;

real dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
real dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
real dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);

real dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
real dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
real dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);

real dedxia = zcb*dedyt - ycb*dedzt - ddrdx;
real dedyia = xcb*dedzt - zcb*dedxt - ddrdy;
real dedzia = ycb*dedxt - xcb*dedyt - ddrdz;

real dedxib = yca*dedzt - zca*dedyt + zdc*dedyu
                 - ydc*dedzu + ddrdx;
real dedyib = zca*dedxt - xca*dedzt + xdc*dedzu
                 - zdc*dedxu + ddrdy;
real dedzib = xca*dedyt - yca*dedxt + ydc*dedxu
                 - xdc*dedyu + ddrdz;

real dedxic = zba*dedyt - yba*dedzt + ydb*dedzu
                 - zdb*dedyu;
real dedyic = xba*dedzt - zba*dedxt + zdb*dedxu
                 - xdb*dedzu;
real dedzic = yba*dedxt - xba*dedyt + xdb*dedyu
                 - ydb*dedxu;

real dedxid = zcb*dedyu - ycb*dedzu;
real dedyid = xcb*dedzu - zcb*dedxu;
real dedzid = ycb*dedxu - xcb*dedyu;

dr = rcb - parameters.y;
real e2 = dr * (amplitude2.x*phi1 + amplitude2.y*phi2 + amplitude2.z*phi3);
dedphi = dr * (amplitude2.x*dphi1 + amplitude2.y*dphi2 + amplitude2.z*dphi3);
ddr = (amplitude2.x*phi1 + amplitude2.y*phi2 + amplitude2.z*phi3) / rcb;

energy += e2;

ddrdx = xcb * ddr;
ddrdy = ycb * ddr;
ddrdz = zcb * ddr;

dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);

dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);

dedxia += zcb*dedyt - ycb*dedzt;
dedyia += xcb*dedzt - zcb*dedxt;
dedzia += ycb*dedxt - xcb*dedyt;

dedxib += yca*dedzt - zca*dedyt + zdc*dedyu
             - ydc*dedzu - ddrdx;
dedyib += zca*dedxt - xca*dedzt + xdc*dedzu
             - zdc*dedxu - ddrdy;
dedzib += xca*dedyt - yca*dedxt + ydc*dedxu
             - xdc*dedyu - ddrdz;

dedxic += zba*dedyt - yba*dedzt + ydb*dedzu
             - zdb*dedyu + ddrdx;
dedyic += xba*dedzt - zba*dedxt + zdb*dedxu
             - xdb*dedzu + ddrdy;
dedzic += yba*dedxt - xba*dedyt + xdb*dedyu
             - ydb*dedxu + ddrdz;

dedxid += zcb*dedyu - ycb*dedzu;
dedyid += xcb*dedzu - zcb*dedxu;
dedzid += ycb*dedxu - xcb*dedyu;

dr = rdc - parameters.z;
real e3 = dr * (amplitude3.x*phi1 + amplitude3.y*phi2 + amplitude3.z*phi3);
dedphi = dr * (amplitude3.x*dphi1 + amplitude3.y*dphi2 + amplitude3.z*dphi3);
ddr = (amplitude3.x*phi1 + amplitude3.y*phi2 + amplitude3.z*phi3) / rdc;

energy += e3;

ddrdx = xdc * ddr;
ddrdy = ydc * ddr;
ddrdz = zdc * ddr;

dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);

dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);

dedxia += zcb*dedyt - ycb*dedzt;
dedyia += xcb*dedzt - zcb*dedxt;
dedzia += ycb*dedxt - xcb*dedyt;

dedxib += yca*dedzt - zca*dedyt + zdc*dedyu
             - ydc*dedzu;
dedyib += zca*dedxt - xca*dedzt + xdc*dedzu
             - zdc*dedxu;
dedzib += xca*dedyt - yca*dedxt + ydc*dedxu
             - xdc*dedyu;

dedxic += zba*dedyt - yba*dedzt + ydb*dedzu
             - zdb*dedyu - ddrdx;
dedyic += xba*dedzt - zba*dedxt + zdb*dedxu
             - xdb*dedzu - ddrdy;
dedzic += yba*dedxt - xba*dedyt + xdb*dedyu
             - ydb*dedxu - ddrdz;

dedxid += zcb*dedyu - ycb*dedzu + ddrdx;
dedyid += xcb*dedzu - zcb*dedxu + ddrdy;
dedzid += ycb*dedxu - xcb*dedyu + ddrdz;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib, -dedyib, -dedzib);
real3 force3 = make_real3(-dedxic, -dedyic, -dedzic);
real3 force4 = make_real3(-dedxid, -dedyid, -dedzid);
