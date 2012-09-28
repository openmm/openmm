float2 angleParams = PARAMS[index];
real xad = pos1.x - pos4.x;
real yad = pos1.y - pos4.y;
real zad = pos1.z - pos4.z;

real xbd = pos2.x - pos4.x;
real ybd = pos2.y - pos4.y;
real zbd = pos2.z - pos4.z;

real xcd = pos3.x - pos4.x;
real ycd = pos3.y - pos4.y;
real zcd = pos3.z - pos4.z;

real xt = yad*zcd - zad*ycd;
real yt = zad*xcd - xad*zcd;
real zt = xad*ycd - yad*xcd;

real rt2 = xt*xt + yt*yt + zt*zt;

real delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;

real xip = pos2.x + xt*delta;
real yip = pos2.y + yt*delta;
real zip = pos2.z + zt*delta;

real xap = pos1.x - xip;
real yap = pos1.y - yip;
real zap = pos1.z - zip;

real xcp = pos3.x - xip;
real ycp = pos3.y - yip;
real zcp = pos3.z - zip;

real rap2 = xap*xap + yap*yap + zap*zap;
real rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;

real xm = ycp*zap - zcp*yap;
real ym = zcp*xap - xcp*zap;
real zm = xcp*yap - ycp*xap;

real rm = max(SQRT(xm*xm + ym*ym + zm*zm), (real) 1e-6f);
real dot = xap*xcp + yap*ycp + zap*zcp;
real product = SQRT(rap2*rcp2);
real cosine = (product > 0 ? (dot/product) : 0);
cosine = max(min(cosine, (real) 1), (real) -1);
real angle = ACOS(cosine);

// if product == 0, set force/energy to 0

real deltaIdeal = (product > 0 ? (angle*RAD_TO_DEG - angleParams.x) : 0);
real deltaIdeal2 = deltaIdeal*deltaIdeal;
real deltaIdeal3 = deltaIdeal*deltaIdeal2;
real deltaIdeal4 = deltaIdeal2*deltaIdeal2;

energy += angleParams.y*deltaIdeal2*(1.0f + CUBIC_K*deltaIdeal + QUARTIC_K*deltaIdeal2 + PENTIC_K*deltaIdeal3 + SEXTIC_K*deltaIdeal4);
real dEdAngle = angleParams.y*deltaIdeal*(2.0f + 3.0f*CUBIC_K*deltaIdeal + 4.0f*QUARTIC_K*deltaIdeal2 + 5.0f*PENTIC_K*deltaIdeal3 + 6.0f*SEXTIC_K*deltaIdeal4);
dEdAngle *= RAD_TO_DEG;

real terma = -dEdAngle/(rap2*rm);
real termc = dEdAngle/(rcp2*rm);

real dedxia = terma * (yap*zm-zap*ym);
real dedyia = terma * (zap*xm-xap*zm);
real dedzia = terma * (xap*ym-yap*xm);

real dedxic = termc * (ycp*zm-zcp*ym);
real dedyic = termc * (zcp*xm-xcp*zm);
real dedzic = termc * (xcp*ym-ycp*xm);

real dedxip = -dedxia - dedxic;
real dedyip = -dedyia - dedyic;
real dedzip = -dedzia - dedzic;

real delta2 = 2.0f*delta;
real ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;

real term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
real dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;

term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
real dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;

term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
real dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;

term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
real dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;

term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
real dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;

term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
real dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;

dedxia = dedxia + dpdxia;
dedyia = dedyia + dpdyia;
dedzia = dedzia + dpdzia;

real dedxib = dedxip;
real dedyib = dedyip;
real dedzib = dedzip;

dedxic = dedxic + dpdxic;
dedyic = dedyic + dpdyic;
dedzic = dedzic + dpdzic;

real dedxid = -dedxia - dedxib - dedxic;
real dedyid = -dedyia - dedyib - dedyic;
real dedzid = -dedzia - dedzib - dedzic;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib, -dedyib, -dedzib);
real3 force3 = make_real3(-dedxic, -dedyic, -dedzic);
real3 force4 = make_real3(-dedxid, -dedyid, -dedzid);