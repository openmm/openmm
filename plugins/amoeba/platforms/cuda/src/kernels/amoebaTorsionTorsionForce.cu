int2 torsionParams = TORSION_PARAMS[index];

real xba = pos2.x - pos1.x;
real yba = pos2.y - pos1.y;
real zba = pos2.z - pos1.z;

real xcb = pos3.x - pos2.x;
real ycb = pos3.y - pos2.y;
real zcb = pos3.z - pos2.z;

real xdc = pos4.x - pos3.x;
real ydc = pos4.y - pos3.y;
real zdc = pos4.z - pos3.z;

real xed = pos5.x - pos4.x;
real yed = pos5.y - pos4.y;
real zed = pos5.z - pos4.z;

real xt = yba*zcb - ycb*zba;
real yt = zba*xcb - zcb*xba;
real zt = xba*ycb - xcb*yba;

real xu = ycb*zdc - ydc*zcb;
real yu = zcb*xdc - zdc*xcb;
real zu = xcb*ydc - xdc*ycb;

real rt2 = xt*xt + yt*yt + zt*zt;
real ru2 = xu*xu + yu*yu + zu*zu;

real rtru = SQRT(rt2 * ru2);

real xv = ydc*zed - yed*zdc;
real yv = zdc*xed - zed*xdc;
real zv = xdc*yed - xed*ydc;

real rv2 = xv*xv + yv*yv + zv*zv;
real rurv = SQRT(ru2 * rv2);

real rcb = SQRT(xcb*xcb + ycb*ycb + zcb*zcb);
real cosine1 = (rtru != 0 ? (xt*xu+yt*yu+zt*zu)/rtru : (real) 0);
cosine1 = (cosine1 > 1 ? (real) 1 : cosine1);
cosine1 = (cosine1 < -1 ? (real) -1 : cosine1);

real angle1 = RAD_TO_DEG * ACOS(cosine1);
real sign = xba*xu + yba*yu + zba*zu;
angle1 = (sign < 0 ? -angle1 : angle1);
real value1 = angle1;

real rdc = SQRT(xdc*xdc + ydc*ydc + zdc*zdc);
real cosine2 = (xu*xv + yu*yv + zu*zv) / rurv;
cosine2 = (cosine2 > 1 ? (real) 1 : cosine2);
cosine2 = (cosine2 < -1 ? (real) -1 : cosine2);
real angle2 = RAD_TO_DEG * ACOS(cosine2);

sign = xcb*xv + ycb*yv + zcb*zv;
angle2 = (sign < 0 ? -angle2 : angle2);
real value2 = angle2;

// check for inverted chirality at the central atom

// if atom2.y < 0, then no chiral check required
// sign is set to 1.0 in this case
// use atom5 for the atom index to avoid warp divergence

int chiralAtomIndex = (torsionParams.x > -1 ? torsionParams.x : atom5);
real4 pos6 = posq[chiralAtomIndex];

real xac = pos6.x - pos3.x;
real yac = pos6.y - pos3.y;
real zac = pos6.z - pos3.z;

real xbc = pos2.x - pos3.x;
real ybc = pos2.y - pos3.y;
real zbc = pos2.z - pos3.z;

// xdc, ydc, zdc appear above

real xdc1 = pos4.x - pos3.x;
real ydc1 = pos4.y - pos3.y;
real zdc1 = pos4.z - pos3.z;

real c1 = ybc*zdc1 - zbc*ydc1;
real c2 = ydc1*zac - zdc1*yac;
real c3 = yac*zbc - zac*ybc;
real vol = xac*c1 + xbc*c2 + xdc1*c3;
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

real xca = pos3.x - pos1.x;
real yca = pos3.y - pos1.y;
real zca = pos3.z - pos1.z;

real xdb = pos4.x - pos2.x;
real ydb = pos4.y - pos2.y;
real zdb = pos4.z - pos2.z;

real dedxt = dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb);
real dedyt = dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb);
real dedzt = dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb);
real dedxu = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb);
real dedyu = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb);
real dedzu = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb);

// compute first derivative components for first angle

real dedxia = zcb*dedyt - ycb*dedzt;
real dedyia = xcb*dedzt - zcb*dedxt;
real dedzia = ycb*dedxt - xcb*dedyt;

real dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
real dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
real dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;

real dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
real dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
real dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;

real dedxid = zcb*dedyu - ycb*dedzu;
real dedyid = xcb*dedzu - zcb*dedxu;
real dedzid = ycb*dedxu - xcb*dedyu;

// chain rule terms for second angle derivative components

real xec = pos5.x - pos3.x;
real yec = pos5.y - pos3.y;
real zec = pos5.z - pos3.z;

real dedxu2 = dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc);
real dedyu2 = dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc);
real dedzu2 = dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc);
real dedxv2 = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc);
real dedyv2 = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc);
real dedzv2 = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc);

// compute first derivative components for second angle

real dedxib2 = zdc*dedyu2 - ydc*dedzu2;
real dedyib2 = xdc*dedzu2 - zdc*dedxu2;
real dedzib2 = ydc*dedxu2 - xdc*dedyu2;
real dedxic2 = ydb*dedzu2 - zdb*dedyu2 + zed*dedyv2 - yed*dedzv2;
real dedyic2 = zdb*dedxu2 - xdb*dedzu2 + xed*dedzv2 - zed*dedxv2;
real dedzic2 = xdb*dedyu2 - ydb*dedxu2 + yed*dedxv2 - xed*dedyv2;
real dedxid2 = zcb*dedyu2 - ycb*dedzu2 + yec*dedzv2 - zec*dedyv2;
real dedyid2 = xcb*dedzu2 - zcb*dedxu2 + zec*dedxv2 - xec*dedzv2;
real dedzid2 = ycb*dedxu2 - xcb*dedyu2 + xec*dedyv2 - yec*dedxv2;
real dedxie2 = zdc*dedyv2 - ydc*dedzv2;
real dedyie2 = xdc*dedzv2 - zdc*dedxv2;
real dedzie2 = ydc*dedxv2 - xdc*dedyv2;

real3 force1 = make_real3(-dedxia, -dedyia, -dedzia);
real3 force2 = make_real3(-dedxib-dedxib2, -dedyib-dedyib2, -dedzib-dedzib2);
real3 force3 = make_real3(-dedxic-dedxic2, -dedyic-dedyic2, -dedzic-dedzic2);
real3 force4 = make_real3(-dedxid-dedxid2, -dedyid-dedyid2, -dedzid-dedzid2);
real3 force5 = make_real3(-dedxie2, -dedyie2, -dedzie2);