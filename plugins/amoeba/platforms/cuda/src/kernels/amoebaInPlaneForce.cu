float2 angleParams = PARAMS[index];
real3 ad = make_real3(pos1.x-pos4.x, pos1.y-pos4.y, pos1.z-pos4.z);
real3 bd = make_real3(pos2.x-pos4.x, pos2.y-pos4.y, pos2.z-pos4.z);
real3 cd = make_real3(pos3.x-pos4.x, pos3.y-pos4.y, pos3.z-pos4.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ad)
APPLY_PERIODIC_TO_DELTA(bd)
APPLY_PERIODIC_TO_DELTA(cd)
#endif

real xt = ad.y*cd.z - ad.z*cd.y;
real yt = ad.z*cd.x - ad.x*cd.z;
real zt = ad.x*cd.y - ad.y*cd.x;

real rt2 = xt*xt + yt*yt + zt*zt;

real delta = -(xt*bd.x + yt*bd.y + zt*bd.z) / rt2;

real xip = pos2.x + xt*delta;
real yip = pos2.y + yt*delta;
real zip = pos2.z + zt*delta;

real3 ap = make_real3(pos1.x-xip, pos1.y-yip, pos1.z-zip);
real3 cp = make_real3(pos3.x-xip, pos3.y-yip, pos3.z-zip);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ap)
APPLY_PERIODIC_TO_DELTA(cp)
#endif

real rap2 = ap.x*ap.x + ap.y*ap.y + ap.z*ap.z;
real rcp2 = cp.x*cp.x + cp.y*cp.y + cp.z*cp.z;

real xm = cp.y*ap.z - cp.z*ap.y;
real ym = cp.z*ap.x - cp.x*ap.z;
real zm = cp.x*ap.y - cp.y*ap.x;

real rm = max(SQRT(xm*xm + ym*ym + zm*zm), (real) 1e-6f);
real dotp = ap.x*cp.x + ap.y*cp.y + ap.z*cp.z;
real product = SQRT(rap2*rcp2);
real cosine = (product > 0 ? (dotp/product) : 0);
cosine = max(min(cosine, (real) 1), (real) -1);
real angle;
if (cosine > 0.99f || cosine < -0.99f) {
    real3 cross_prod = cross(ap, cp);
    angle = ASIN(SQRT(dot(cross_prod, cross_prod)/(rap2*rcp2)))*RAD_TO_DEG;
    if (cosine < 0.0f)
        angle = 180-angle;
}
else
    angle = ACOS(cosine)*RAD_TO_DEG;

// if product == 0, set force/energy to 0

real deltaIdeal = (product > 0 ? (angle - angleParams.x) : 0);
real deltaIdeal2 = deltaIdeal*deltaIdeal;
real deltaIdeal3 = deltaIdeal*deltaIdeal2;
real deltaIdeal4 = deltaIdeal2*deltaIdeal2;

energy += angleParams.y*deltaIdeal2*(1.0f + CUBIC_K*deltaIdeal + QUARTIC_K*deltaIdeal2 + PENTIC_K*deltaIdeal3 + SEXTIC_K*deltaIdeal4);
real dEdAngle = angleParams.y*deltaIdeal*(2.0f + 3.0f*CUBIC_K*deltaIdeal + 4.0f*QUARTIC_K*deltaIdeal2 + 5.0f*PENTIC_K*deltaIdeal3 + 6.0f*SEXTIC_K*deltaIdeal4);
dEdAngle *= RAD_TO_DEG;

real terma = -dEdAngle/(rap2*rm);
real termc = dEdAngle/(rcp2*rm);

real dedxia = terma * (ap.y*zm-ap.z*ym);
real dedyia = terma * (ap.z*xm-ap.x*zm);
real dedzia = terma * (ap.x*ym-ap.y*xm);

real dedxic = termc * (cp.y*zm-cp.z*ym);
real dedyic = termc * (cp.z*xm-cp.x*zm);
real dedzic = termc * (cp.x*ym-cp.y*xm);

real dedxip = -dedxia - dedxic;
real dedyip = -dedyia - dedyic;
real dedzip = -dedzia - dedzic;

real delta2 = 2.0f*delta;
real ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;

real term = (cd.z*bd.y-cd.y*bd.z) + delta2*(yt*cd.z-zt*cd.y);
real dpdxia = delta*(cd.y*dedzip-cd.z*dedyip) + term*ptrt2;

term = (cd.x*bd.z-cd.z*bd.x) + delta2*(zt*cd.x-xt*cd.z);
real dpdyia = delta*(cd.z*dedxip-cd.x*dedzip) + term*ptrt2;

term = (cd.y*bd.x-cd.x*bd.y) + delta2*(xt*cd.y-yt*cd.x);
real dpdzia = delta*(cd.x*dedyip-cd.y*dedxip) + term*ptrt2;

term = (ad.y*bd.z-ad.z*bd.y) + delta2*(zt*ad.y-yt*ad.z);
real dpdxic = delta*(ad.z*dedyip-ad.y*dedzip) + term*ptrt2;

term = (ad.z*bd.x-ad.x*bd.z) + delta2*(xt*ad.z-zt*ad.x);
real dpdyic = delta*(ad.x*dedzip-ad.z*dedxip) + term*ptrt2;

term = (ad.x*bd.y-ad.y*bd.x) + delta2*(yt*ad.x-xt*ad.y);
real dpdzic = delta*(ad.y*dedxip-ad.x*dedyip) + term*ptrt2;

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
