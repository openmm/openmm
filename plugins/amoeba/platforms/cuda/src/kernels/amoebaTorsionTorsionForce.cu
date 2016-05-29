int2 torsionParams = TORSION_PARAMS[index];

real3 ba = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
real3 cb = make_real3(pos3.x-pos2.x, pos3.y-pos2.y, pos3.z-pos2.z);
real3 dc = make_real3(pos4.x-pos3.x, pos4.y-pos3.y, pos4.z-pos3.z);
real3 ed = make_real3(pos5.x-pos4.x, pos5.y-pos4.y, pos5.z-pos4.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ba)
APPLY_PERIODIC_TO_DELTA(cb)
APPLY_PERIODIC_TO_DELTA(dc)
APPLY_PERIODIC_TO_DELTA(ed)
#endif

real xt = ba.y*cb.z - cb.y*ba.z;
real yt = ba.z*cb.x - cb.z*ba.x;
real zt = ba.x*cb.y - cb.x*ba.y;

real xu = cb.y*dc.z - dc.y*cb.z;
real yu = cb.z*dc.x - dc.z*cb.x;
real zu = cb.x*dc.y - dc.x*cb.y;

real rt2 = xt*xt + yt*yt + zt*zt;
real ru2 = xu*xu + yu*yu + zu*zu;

real rtru = SQRT(rt2 * ru2);

real xv = dc.y*ed.z - ed.y*dc.z;
real yv = dc.z*ed.x - ed.z*dc.x;
real zv = dc.x*ed.y - ed.x*dc.y;

real rv2 = xv*xv + yv*yv + zv*zv;
real rurv = SQRT(ru2 * rv2);

real rcb = SQRT(cb.x*cb.x + cb.y*cb.y + cb.z*cb.z);
real cosine1 = (rtru != 0 ? (xt*xu+yt*yu+zt*zu)/rtru : (real) 0);
cosine1 = (cosine1 > 1 ? (real) 1 : cosine1);
cosine1 = (cosine1 < -1 ? (real) -1 : cosine1);
real angle1;
if (cosine1 > 0.99f || cosine1 < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real3 cross_prod = cross(make_real3(xt, yt, zt), make_real3(xu, yu, zu));
    angle1 = RAD_TO_DEG*ASIN(SQRT(dot(cross_prod, cross_prod)/(rt2*ru2)));
    if (cosine1 < 0.0f)
        angle1 = 180-angle1;
}
else
   angle1 = RAD_TO_DEG*ACOS(cosine1);
real sign = ba.x*xu + ba.y*yu + ba.z*zu;
angle1 = (sign < 0 ? -angle1 : angle1);
real value1 = angle1;

real rdc = SQRT(dc.x*dc.x + dc.y*dc.y + dc.z*dc.z);
real cosine2 = (xu*xv + yu*yv + zu*zv) / rurv;
cosine2 = (cosine2 > 1 ? (real) 1 : cosine2);
cosine2 = (cosine2 < -1 ? (real) -1 : cosine2);
real angle2;
if (cosine2 > 0.99f || cosine2 < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real3 cross_prod = cross(make_real3(xu, yu, zu), make_real3(xv, yv, zv));
    angle2 = RAD_TO_DEG*ASIN(SQRT(dot(cross_prod, cross_prod)/(ru2*rv2)));
    if (cosine2 < 0.0f)
        angle2 = 180-angle2;
}
else
   angle2 = RAD_TO_DEG*ACOS(cosine2);
sign = cb.x*xv + cb.y*yv + cb.z*zv;
angle2 = (sign < 0 ? -angle2 : angle2);
real value2 = angle2;

// check for inverted chirality at the central atom

// if atom2.y < 0, then no chiral check required
// sign is set to 1.0 in this case
// use atom5 for the atom index to avoid warp divergence

int chiralAtomIndex = (torsionParams.x > -1 ? torsionParams.x : atom5);
real4 pos6 = posq[chiralAtomIndex];

real3 ac = make_real3(pos6.x-pos3.x, pos6.y-pos3.y, pos6.z-pos3.z);
real3 bc = make_real3(pos2.x-pos3.x, pos2.y-pos3.y, pos2.z-pos3.z);
real3 dc1 = make_real3(pos4.x-pos3.x, pos4.y-pos3.y, pos4.z-pos3.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ac)
APPLY_PERIODIC_TO_DELTA(bc)
APPLY_PERIODIC_TO_DELTA(dc1)
#endif

real c1 = bc.y*dc1.z - bc.z*dc1.y;
real c2 = dc1.y*ac.z - dc1.z*ac.y;
real c3 = ac.y*bc.z - ac.z*bc.y;
real vol = ac.x*c1 + bc.x*c2 + dc1.x*c3;
sign = (vol > 0 ? (real) 1 : (real) -1);
sign = (torsionParams.x < 0 ? (real) 1 : sign);
value1 *= sign;
value2 *= sign;

// use bicubic interpolation to compute spline values
// compute indices into grid based on angles

float4 gridParams = GRID_PARAMS[torsionParams.y];
int index1 = (int) ((value1 - gridParams.y)/gridParams.z + 1.0e-05f);
real fIndex = (real) index1;
real x1l = gridParams.z*fIndex + gridParams.y;
real x1u = x1l + gridParams.z;

int index2 = (int) ((value2 - gridParams.y)/gridParams.z + 1.0e-05f);
fIndex = (real) index2;
real x2l = gridParams.z*fIndex + gridParams.y;
real x2u = x2l + gridParams.z;

int posIndex1 = index2 + index1*(int) gridParams.w;
posIndex1 += (int) gridParams.x;

int posIndex2 = index2 + (index1+1)*(int) gridParams.w;
posIndex2 += (int) gridParams.x;

// load grid points surrounding angle

real4 y;
real4 y1;
real4 y2;
real4 y12;

int localIndex = posIndex1;
y.x = GRID_VALUES[localIndex].x;
y1.x = GRID_VALUES[localIndex].y;
y2.x = GRID_VALUES[localIndex].z;
y12.x = GRID_VALUES[localIndex].w;

localIndex = posIndex2;
y.y = GRID_VALUES[localIndex].x;
y1.y = GRID_VALUES[localIndex].y;
y2.y = GRID_VALUES[localIndex].z;
y12.y = GRID_VALUES[localIndex].w;

localIndex = posIndex2 + 1;
y.z = GRID_VALUES[localIndex].x;
y1.z = GRID_VALUES[localIndex].y;
y2.z = GRID_VALUES[localIndex].z;
y12.z = GRID_VALUES[localIndex].w;

localIndex = posIndex1 + 1;
y.w = GRID_VALUES[localIndex].x;
y1.w = GRID_VALUES[localIndex].y;
y2.w = GRID_VALUES[localIndex].z;
y12.w = GRID_VALUES[localIndex].w;

// perform interpolation

real e;
real dedang1;
real dedang2;

bicubic(y, y1, y2, y12, value1, x1l, x1u, value2, x2l, x2u, &e, &dedang1, &dedang2);
energy += e;
dedang1 *= sign * RAD_TO_DEG;
dedang2 *= sign * RAD_TO_DEG;

// chain rule terms for first angle derivative components

real3 ca = make_real3(pos3.x-pos1.x, pos3.y-pos1.y, pos3.z-pos1.z);
real3 db = make_real3(pos4.x-pos2.x, pos4.y-pos2.y, pos4.z-pos2.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ca)
APPLY_PERIODIC_TO_DELTA(db)
#endif

real dedxt = dedang1 * (yt*cb.z - cb.y*zt) / (rt2*rcb);
real dedyt = dedang1 * (zt*cb.x - cb.z*xt) / (rt2*rcb);
real dedzt = dedang1 * (xt*cb.y - cb.x*yt) / (rt2*rcb);
real dedxu = -dedang1 * (yu*cb.z - cb.y*zu) / (ru2*rcb);
real dedyu = -dedang1 * (zu*cb.x - cb.z*xu) / (ru2*rcb);
real dedzu = -dedang1 * (xu*cb.y - cb.x*yu) / (ru2*rcb);

// compute first derivative components for first angle

real dedxia = cb.z*dedyt - cb.y*dedzt;
real dedyia = cb.x*dedzt - cb.z*dedxt;
real dedzia = cb.y*dedxt - cb.x*dedyt;

real dedxib = ca.y*dedzt - ca.z*dedyt + dc.z*dedyu - dc.y*dedzu;
real dedyib = ca.z*dedxt - ca.x*dedzt + dc.x*dedzu - dc.z*dedxu;
real dedzib = ca.x*dedyt - ca.y*dedxt + dc.y*dedxu - dc.x*dedyu;

real dedxic = ba.z*dedyt - ba.y*dedzt + db.y*dedzu - db.z*dedyu;
real dedyic = ba.x*dedzt - ba.z*dedxt + db.z*dedxu - db.x*dedzu;
real dedzic = ba.y*dedxt - ba.x*dedyt + db.x*dedyu - db.y*dedxu;

real dedxid = cb.z*dedyu - cb.y*dedzu;
real dedyid = cb.x*dedzu - cb.z*dedxu;
real dedzid = cb.y*dedxu - cb.x*dedyu;

// chain rule terms for second angle derivative components

real3 ec = make_real3(pos5.x-pos3.x, pos5.y-pos3.y, pos5.z-pos3.z);

#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(ec)
#endif

real dedxu2 = dedang2 * (yu*dc.z - dc.y*zu) / (ru2*rdc);
real dedyu2 = dedang2 * (zu*dc.x - dc.z*xu) / (ru2*rdc);
real dedzu2 = dedang2 * (xu*dc.y - dc.x*yu) / (ru2*rdc);
real dedxv2 = -dedang2 * (yv*dc.z - dc.y*zv) / (rv2*rdc);
real dedyv2 = -dedang2 * (zv*dc.x - dc.z*xv) / (rv2*rdc);
real dedzv2 = -dedang2 * (xv*dc.y - dc.x*yv) / (rv2*rdc);

// compute first derivative components for second angle

real dedxib2 = dc.z*dedyu2 - dc.y*dedzu2;
real dedyib2 = dc.x*dedzu2 - dc.z*dedxu2;
real dedzib2 = dc.y*dedxu2 - dc.x*dedyu2;
real dedxic2 = db.y*dedzu2 - db.z*dedyu2 + ed.z*dedyv2 - ed.y*dedzv2;
real dedyic2 = db.z*dedxu2 - db.x*dedzu2 + ed.x*dedzv2 - ed.z*dedxv2;
real dedzic2 = db.x*dedyu2 - db.y*dedxu2 + ed.y*dedxv2 - ed.x*dedyv2;
real dedxid2 = cb.z*dedyu2 - cb.y*dedzu2 + ec.y*dedzv2 - ec.z*dedyv2;
real dedyid2 = cb.x*dedzu2 - cb.z*dedxu2 + ec.z*dedxv2 - ec.x*dedzv2;
real dedzid2 = cb.y*dedxu2 - cb.x*dedyu2 + ec.x*dedyv2 - ec.y*dedxv2;
real dedxie2 = dc.z*dedyv2 - dc.y*dedzv2;
real dedyie2 = dc.x*dedzv2 - dc.z*dedxv2;
real dedzie2 = dc.y*dedxv2 - dc.x*dedyv2;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib-dedxib2, -dedyib-dedyib2, -dedzib-dedzib2);
real3 force3 = make_real3(-dedxic-dedxic2, -dedyic-dedyic2, -dedzic-dedzic2);
real3 force4 = make_real3(-dedxid-dedxid2, -dedyid-dedyid2, -dedzid-dedzid2);
real3 force5 = make_real3(-dedxie2, -dedyie2, -dedzie2);