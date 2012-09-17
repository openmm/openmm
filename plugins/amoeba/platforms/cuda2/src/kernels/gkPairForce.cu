/**
 * This defines three different closely related functions, depending on which constant (F1, T1, or T3) is defined.
 */

#if defined F1
__device__ void computeOneInteractionF1(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force) {
#elif defined F2
__device__ void computeOneInteractionF2(AtomData2& atom1, volatile AtomData2& atom2, real& outputEnergy, real3& force) {
#elif defined T1
__device__ void computeOneInteractionT1(AtomData2& atom1, volatile AtomData2& atom2, real3& torque) {
#elif defined T2
__device__ void computeOneInteractionT2(AtomData2& atom1, volatile AtomData2& atom2, real3& torque) {
#elif defined B1 && defined B2
__device__ void computeOneInteractionB1B2(AtomData2& atom1, volatile AtomData2& atom2) {
#endif

    const real fc = EPSILON_FACTOR*GK_FC;
    const real fd = EPSILON_FACTOR*GK_FD;
    const real fq = EPSILON_FACTOR*GK_FQ;

#if defined F2 || defined B2
    real sxi = atom1.inducedDipole.x + atom1.inducedDipolePolar.x;
    real syi = atom1.inducedDipole.y + atom1.inducedDipolePolar.y;
    real szi = atom1.inducedDipole.z + atom1.inducedDipolePolar.z;
#endif

#if defined F2 || defined T2 || defined B2
    real sxk = atom2.inducedDipole.x + atom2.inducedDipolePolar.x;
    real syk = atom2.inducedDipole.y + atom2.inducedDipolePolar.y;
    real szk = atom2.inducedDipole.z + atom2.inducedDipolePolar.z;
#endif

    // decide whether to compute the current interaction;

    real xr = atom2.pos.x - atom1.pos.x;
    real yr = atom2.pos.y - atom1.pos.y;
    real zr = atom2.pos.z - atom1.pos.z;

    real xr2 = xr*xr;
    real yr2 = yr*yr;
    real zr2 = zr*zr;
    real r2 = xr2 + yr2 + zr2;

    real rb2 = atom1.bornRadius*atom2.bornRadius;

    real expterm = EXP(-r2/(GK_C*rb2));
    real expc = expterm/GK_C;
    real expcr = r2*expterm/(GK_C*GK_C*rb2*rb2);
    real dexpc = -2 / (GK_C*rb2);
    real dexpcr = 2 / (GK_C*rb2*rb2);
    real dgfdr = 0.5f*expterm*(1 + r2/(rb2*GK_C));
    real gf2 = 1 / (r2 + rb2*expterm);

    real gf = SQRT(gf2);
    real gf3 = gf2*gf;
    real gf5 = gf3*gf2;
    real gf7 = gf5*gf2;
    real gf9 = gf7*gf2;
    real gf11 = gf9*gf2;

    // reaction potential auxiliary terms;

    real a00 = gf;
    real a10 = -gf3;
    real a20 = 3*gf5;
    real a30 = -15*gf7;
    real a40 = 105*gf9;
    real a50 = -945*gf11;

    // Born radii derivatives of reaction potential auxiliary terms;

    real b00 = dgfdr*a10;
    real b10 = dgfdr*a20;
    real b20 = dgfdr*a30;
    real b30 = dgfdr*a40;
    real b40 = dgfdr*a50;

    // reaction potential gradient auxiliary terms;

    real expc1 = 1 - expc;
    real a01 = expc1*a10;
    real a11 = expc1*a20;
    real a21 = expc1*a30;
    real a31 = expc1*a40;
    real a41 = expc1*a50;

    // Born radii derivs of reaction potential gradient auxiliary terms;

    real b01 = b10 - expcr*a10 - expc*b10;
    real b11 = b20 - expcr*a20 - expc*b20;
    real b21 = b30 - expcr*a30 - expc*b30;
    real b31 = b40 - expcr*a40 - expc*b40;

    // 2nd reaction potential gradient auxiliary terms;

    real expcdexpc = -expc*dexpc;
    real a02 = expc1*a11 + expcdexpc*a10;
    real a12 = expc1*a21 + expcdexpc*a20;
    real a22 = expc1*a31 + expcdexpc*a30;
    real a32 = expc1*a41 + expcdexpc*a40;

    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     real b02 = b11 - (expcr*(a11 + dexpc*a10) + expc*(b11 + dexpcr*a10 + dexpc*b10));
     real b12 = b21 - (expcr*(a21 + dexpc*a20) + expc*(b21 + dexpcr*a20 + dexpc*b20));
     real b22 = b31 - (expcr*(a31 + dexpc*a30) + expc*(b31 + dexpcr*a30 + dexpc*b30));

    // 3rd reaction potential gradient auxiliary terms

    expcdexpc = 2*expcdexpc;
    real a03 = expc1*a12 + expcdexpc*a11;
    real a13 = expc1*a22 + expcdexpc*a21;
    real a23 = expc1*a32 + expcdexpc*a31;

    expcdexpc = -expc*dexpc*dexpc;
    a03 = a03 + expcdexpc*a10;
    a13 = a13 + expcdexpc*a20;
    a23 = a23 + expcdexpc*a30;

    // multiply the auxillary terms by their dieletric functions;

    a00 *= fc;
    a01 *= fc;
    a02 *= fc;
    a03 *= fc;

    b00 *= fc;
    b01 *= fc;
    b02 *= fc;

    a10 *= fd;
    a11 *= fd;
    a12 *= fd;
    a13 *= fd;

    b10 *= fd;
    b11 *= fd;
    b12 *= fd;

    a20 *= fq;
    a21 *= fq;
    a22 *= fq;
    a23 *= fq;

    b20 *= fq;
    b21 *= fq;
    b22 *= fq;

    // unweighted reaction potential tensor

#if defined F2
    real energy = -a10*atom2.q*(atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr);
    energy += a10*atom1.q*(atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr);
#endif

#if defined F1
    real energy = 2*atom1.q*atom2.q*a00;
    energy += -a10*atom2.q*(atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr);
    energy += a10*atom1.q*(atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr);
    energy += a20*atom2.q*(atom1.quadrupoleXX*xr2 + atom1.quadrupoleYY*yr2 + atom1.quadrupoleZZ*zr2 + 2*(atom1.quadrupoleXY*xr*yr + atom1.quadrupoleXZ*xr*zr + atom1.quadrupoleYZ*yr*zr));
    energy += a20*atom1.q*(atom2.quadrupoleXX*xr2 + atom2.quadrupoleYY*yr2 + atom2.quadrupoleZZ*zr2 + 2*(atom2.quadrupoleXY*xr*yr + atom2.quadrupoleXZ*xr*zr + atom2.quadrupoleYZ*yr*zr));
#endif

    // Born radii derivs of unweighted reaction potential tensor

#if defined B1
    real dsumdrB1 = 2*(atom1.q*atom2.q*b00);
    dsumdrB1 -= b10*atom2.q*(atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr);
    dsumdrB1 += b10*atom1.q*(atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr);
#endif
#if defined B2
    real dsumdrB2 = -b10*atom2.q*(sxi*xr + syi*yr + szi*zr);
    dsumdrB2 += b10*atom1.q*(sxk*xr + syk*yr + szk*zr);
#endif

#if defined B1
    real gqxx21 = xr*xr;
    real gqyy21 = yr*yr;
    real gqzz21 = zr*zr;

    real gqxy21 = xr*yr;
    real gqxz21 = xr*zr;
    real gqyz21 = yr*zr;
    dsumdrB1 += b20*atom2.q*(atom1.quadrupoleXX*gqxx21 + atom1.quadrupoleYY*gqyy21 + atom1.quadrupoleZZ*gqzz21 + 2*(atom1.quadrupoleXY*gqxy21 + atom1.quadrupoleXZ*gqxz21 + atom1.quadrupoleYZ*gqyz21));
    dsumdrB1 += b20*atom1.q*(atom2.quadrupoleXX*gqxx21 + atom2.quadrupoleYY*gqyy21 + atom2.quadrupoleZZ*gqzz21 + 2*(atom2.quadrupoleXY*gqxy21 + atom2.quadrupoleXZ*gqxz21 + atom2.quadrupoleYZ*gqyz21));
#endif

#if defined F1
    energy += a01*atom1.q*(atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr);
    energy -= a01*atom2.q*(atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr);
    real factor = a01*2*atom1.q*atom2.q;
    real dedx = factor*xr;
    real dedy = factor*yr;
    real dedz = factor*zr;
#endif
#if defined F2 
    energy += a01*atom1.q*(atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr);
    energy -= a01*atom2.q*(atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr);
#endif

#if defined F1 || defined F2 || defined T1 || defined T2
    real gux2 = a10 + xr*xr*a11;
    real gux3 = xr*yr*a11;
    real gux4 = xr*zr*a11;

    real guy3 = a10 + yr*yr*a11;
    real guy4 = yr*zr*a11;
    real guz4 = a10 + zr*zr*a11;
#if defined T1
    real guy2 = gux3;
    real guz2 = gux4;
    real guz3 = guy4;
#endif
#if defined T2
    real fid1 = sxk*gux2 + syk*gux3 + szk*gux4;
    real fid2 = sxk*gux3 + syk*guy3 + szk*guy4;
    real fid3 = sxk*gux4 + syk*guy4 + szk*guz4;

    real trqi1 = atom1.dipole.y*fid3 - atom1.dipole.z*fid2;
    real trqi2 = atom1.dipole.z*fid1 - atom1.dipole.x*fid3;
    real trqi3 = atom1.dipole.x*fid2 - atom1.dipole.y*fid1;
#endif

#if defined F1
    energy -= 2*(atom1.dipole.x*(atom2.dipole.x*gux2 + atom2.dipole.y*gux3 + atom2.dipole.z*gux4) +
                             atom1.dipole.y*(atom2.dipole.x*gux3 + atom2.dipole.y*guy3 + atom2.dipole.z*guy4) +
                             atom1.dipole.z*(atom2.dipole.x*gux4 + atom2.dipole.y*guy4 + atom2.dipole.z*guz4));

    dedx -= atom2.q*(atom1.dipole.x*gux2 + atom1.dipole.y*gux3 + atom1.dipole.z*gux4);
    dedx += atom1.q*(atom2.dipole.x*gux2 + atom2.dipole.y*gux3 + atom2.dipole.z*gux4);

    dedy -= atom2.q*(atom1.dipole.x*gux3 + atom1.dipole.y*guy3 + atom1.dipole.z*guy4);
    dedy += atom1.q*(atom2.dipole.x*gux3 + atom2.dipole.y*guy3 + atom2.dipole.z*guy4);

    dedz -= atom2.q*(atom1.dipole.x*gux4 + atom1.dipole.y*guy4 + atom1.dipole.z*guz4);
    dedz += atom1.q*(atom2.dipole.x*gux4 + atom2.dipole.y*guy4 + atom2.dipole.z*guz4);
#endif
#if defined F2
    energy -= 2*(
                       atom1.dipole.x*(atom2.inducedDipole.x*gux2 + atom2.inducedDipole.y*gux3 + atom2.inducedDipole.z*gux4) + 
                       atom1.dipole.y*(atom2.inducedDipole.x*gux3 + atom2.inducedDipole.y*guy3 + atom2.inducedDipole.z*guy4) + 
                       atom1.dipole.z*(atom2.inducedDipole.x*gux4 + atom2.inducedDipole.y*guy4 + atom2.inducedDipole.z*guz4) + 
                       atom2.dipole.x*(atom1.inducedDipole.x*gux2 + atom1.inducedDipole.y*gux3 + atom1.inducedDipole.z*gux4) + 
                       atom2.dipole.y*(atom1.inducedDipole.x*gux3 + atom1.inducedDipole.y*guy3 + atom1.inducedDipole.z*guy4) + 
                       atom2.dipole.z*(atom1.inducedDipole.x*gux4 + atom1.inducedDipole.y*guy4 + atom1.inducedDipole.z*guz4));

    real dpdx = atom1.q*(sxk*gux2 + syk*gux3 + szk*gux4);
    dpdx -= atom2.q*(sxi*gux2 + syi*gux3 + szi*gux4);

    real dpdy = atom1.q*(sxk*gux3 + syk*guy3 + szk*guy4);
    dpdy -= atom2.q*(sxi*gux3 + syi*guy3 + szi*guy4);

    real dpdz = atom1.q*(sxk*gux4 + syk*guy4 + szk*guz4);
    dpdz -= atom2.q*(sxi*gux4 + syi*guy4 + szi*guz4);

#endif
    real gqxx2 = xr*(2*a20 + xr*xr*a21);
    real gqxx3 = yr*xr*xr*a21;
    real gqxx4 = zr*xr*xr*a21;
    real gqyy2 = xr*yr*yr*a21;
    real gqyy3 = yr*(2*a20 + yr*yr*a21);
    real gqyy4 = zr*yr*yr*a21;
    real gqzz2 = xr*zr*zr*a21;
    real gqzz3 = yr*zr*zr*a21;
    real gqzz4 = zr*(2*a20 + zr*zr*a21);
    real gqxy2 = yr*(a20 + xr*xr*a21);
    real gqxy3 = xr*(a20 + yr*yr*a21);
    real gqxy4 = zr*xr*yr*a21;
    real gqxz2 = zr*(a20 + xr*xr*a21);
    real gqxz4 = xr*(a20 + zr*zr*a21);
    real gqyz3 = zr*(a20 + yr*yr*a21);
    real gqyz4 = yr*(a20 + zr*zr*a21);
#if defined T1 || defined T2
    real gqxz3 = gqxy4;
    real gqyz2 = gqxy4;
#endif

#if defined F1
    energy += atom2.dipole.x*(atom1.quadrupoleXX*gqxx2 + atom1.quadrupoleYY*gqyy2 + atom1.quadrupoleZZ*gqzz2 + 2*(atom1.quadrupoleXY*gqxy2 + atom1.quadrupoleXZ*gqxz2 + atom1.quadrupoleYZ*gqxy4)) +
                       atom2.dipole.y*(atom1.quadrupoleXX*gqxx3 + atom1.quadrupoleYY*gqyy3 + atom1.quadrupoleZZ*gqzz3 + 2*(atom1.quadrupoleXY*gqxy3 + atom1.quadrupoleXZ*gqxy4 + atom1.quadrupoleYZ*gqyz3)) +
                       atom2.dipole.z*(atom1.quadrupoleXX*gqxx4 + atom1.quadrupoleYY*gqyy4 + atom1.quadrupoleZZ*gqzz4 + 2*(atom1.quadrupoleXY*gqxy4 + atom1.quadrupoleXZ*gqxz4 + atom1.quadrupoleYZ*gqyz4));
    energy -= atom1.dipole.x*(atom2.quadrupoleXX*gqxx2 + atom2.quadrupoleYY*gqyy2 + atom2.quadrupoleZZ*gqzz2 + 2*(atom2.quadrupoleXY*gqxy2 + atom2.quadrupoleXZ*gqxz2 + atom2.quadrupoleYZ*gqxy4)) +
                       atom1.dipole.y*(atom2.quadrupoleXX*gqxx3 + atom2.quadrupoleYY*gqyy3 + atom2.quadrupoleZZ*gqzz3 + 2*(atom2.quadrupoleXY*gqxy3 + atom2.quadrupoleXZ*gqxy4 + atom2.quadrupoleYZ*gqyz3)) +
                       atom1.dipole.z*(atom2.quadrupoleXX*gqxx4 + atom2.quadrupoleYY*gqyy4 + atom2.quadrupoleZZ*gqzz4 + 2*(atom2.quadrupoleXY*gqxy4 + atom2.quadrupoleXZ*gqxz4 + atom2.quadrupoleYZ*gqyz4));

    dedx += atom2.q*(atom1.quadrupoleXX*gqxx2 + atom1.quadrupoleYY*gqyy2 + atom1.quadrupoleZZ*gqzz2 + 2*(atom1.quadrupoleXY*gqxy2 + atom1.quadrupoleXZ*gqxz2 + atom1.quadrupoleYZ*gqxy4));
    dedx += atom1.q*(atom2.quadrupoleXX*gqxx2 + atom2.quadrupoleYY*gqyy2 + atom2.quadrupoleZZ*gqzz2 + 2*(atom2.quadrupoleXY*gqxy2 + atom2.quadrupoleXZ*gqxz2 + atom2.quadrupoleYZ*gqxy4));

    dedy += atom2.q*(atom1.quadrupoleXX*gqxx3 + atom1.quadrupoleYY*gqyy3 + atom1.quadrupoleZZ*gqzz3 + 2*(atom1.quadrupoleXY*gqxy3 + atom1.quadrupoleXZ*gqxy4 + atom1.quadrupoleYZ*gqyz3));
    dedy += atom1.q*(atom2.quadrupoleXX*gqxx3 + atom2.quadrupoleYY*gqyy3 + atom2.quadrupoleZZ*gqzz3 + 2*(atom2.quadrupoleXY*gqxy3 + atom2.quadrupoleXZ*gqxy4 + atom2.quadrupoleYZ*gqyz3));

    dedz += atom2.q*(atom1.quadrupoleXX*gqxx4 + atom1.quadrupoleYY*gqyy4 + atom1.quadrupoleZZ*gqzz4 + 2*(atom1.quadrupoleXY*gqxy4 + atom1.quadrupoleXZ*gqxz4 + atom1.quadrupoleYZ*gqyz4));
    dedz += atom1.q*(atom2.quadrupoleXX*gqxx4 + atom2.quadrupoleYY*gqyy4 + atom2.quadrupoleZZ*gqzz4 + 2*(atom2.quadrupoleXY*gqxy4 + atom2.quadrupoleXZ*gqxz4 + atom2.quadrupoleYZ*gqyz4));
#endif

#if defined F2
    energy += atom2.inducedDipole.x*(atom1.quadrupoleXX*gqxx2 + atom1.quadrupoleYY*gqyy2 + atom1.quadrupoleZZ*gqzz2 + 2*(atom1.quadrupoleXY*gqxy2 + atom1.quadrupoleXZ*gqxz2 + atom1.quadrupoleYZ*gqxy4)) +
              atom2.inducedDipole.y*(atom1.quadrupoleXX*gqxx3 + atom1.quadrupoleYY*gqyy3 + atom1.quadrupoleZZ*gqzz3 + 2*(atom1.quadrupoleXY*gqxy3 + atom1.quadrupoleXZ*gqxy4 + atom1.quadrupoleYZ*gqyz3)) +
              atom2.inducedDipole.z*(atom1.quadrupoleXX*gqxx4 + atom1.quadrupoleYY*gqyy4 + atom1.quadrupoleZZ*gqzz4 + 2*(atom1.quadrupoleXY*gqxy4 + atom1.quadrupoleXZ*gqxz4 + atom1.quadrupoleYZ*gqyz4));

    energy -= atom1.inducedDipole.x*(atom2.quadrupoleXX*gqxx2 + atom2.quadrupoleYY*gqyy2 + atom2.quadrupoleZZ*gqzz2 + 2*(atom2.quadrupoleXY*gqxy2 + atom2.quadrupoleXZ*gqxz2 + atom2.quadrupoleYZ*gqxy4)) +
              atom1.inducedDipole.y*(atom2.quadrupoleXX*gqxx3 + atom2.quadrupoleYY*gqyy3 + atom2.quadrupoleZZ*gqzz3 + 2*(atom2.quadrupoleXY*gqxy3 + atom2.quadrupoleXZ*gqxy4 + atom2.quadrupoleYZ*gqyz3)) +
              atom1.inducedDipole.z*(atom2.quadrupoleXX*gqxx4 + atom2.quadrupoleYY*gqyy4 + atom2.quadrupoleZZ*gqzz4 + 2*(atom2.quadrupoleXY*gqxy4 + atom2.quadrupoleXZ*gqxz4 + atom2.quadrupoleYZ*gqyz4));

#endif
#endif

    // Born derivs of the unweighted reaction potential gradient tensor

#if defined B1
    dsumdrB1 += b01*atom1.q*(atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr);
    dsumdrB1 -= b01*atom2.q*(atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr);
#endif
#if defined B2
    dsumdrB2 += b01*atom1.q*(sxk*xr+ syk*yr + szk*zr);
    dsumdrB2 -= b01*atom2.q*(sxi*xr+ syi*yr + szi*zr);
#endif

#if defined B1 || defined B2
    real gux22 = b10 + xr2*b11;
    real gux23 = xr*yr*b11;
    real gux24 = xr*zr*b11;
    real guy22 = gux23;
    real guy23 = b10 + yr2*b11;
    real guy24 = yr*zr*b11;
    real guz22 = gux24;
    real guz23 = guy24;
    real guz24 = b10 + zr2*b11;
#if defined B1
    dsumdrB1 -= 2*(atom1.dipole.x*(atom2.dipole.x*gux22 + atom2.dipole.y*guy22 + atom2.dipole.z*guz22) +
                              atom1.dipole.y*(atom2.dipole.x*gux23 + atom2.dipole.y*guy23 + atom2.dipole.z*guz23) +
                              atom1.dipole.z*(atom2.dipole.x*gux24 + atom2.dipole.y*guy24 + atom2.dipole.z*guz24));
#endif
#if defined B2
    dsumdrB2 -= 2*(atom1.dipole.x*(sxk*gux22 + syk*guy22 + szk*guz22) +
                            atom1.dipole.y*(sxk*gux23 + syk*guy23 + szk*guz23) +
                            atom1.dipole.z*(sxk*gux24 + syk*guy24 + szk*guz24) +
                            atom2.dipole.x*(sxi*gux22 + syi*guy22 + szi*guz22) +
                            atom2.dipole.y*(sxi*gux23 + syi*guy23 + szi*guz23) +
                            atom2.dipole.z*(sxi*gux24 + syi*guy24 + szi*guz24));

#ifndef DIRECT_POLARIZATION
    dsumdrB2 -= 2*(atom1.inducedDipole.x*(atom2.inducedDipolePolar.x*gux22 + atom2.inducedDipolePolar.y*gux23 + atom2.inducedDipolePolar.z*gux24)
                            + atom1.inducedDipole.y*(atom2.inducedDipolePolar.x*guy22 + atom2.inducedDipolePolar.y*guy23 + atom2.inducedDipolePolar.z*guy24)
                            + atom1.inducedDipole.z*(atom2.inducedDipolePolar.x*guz22 + atom2.inducedDipolePolar.y*guz23 + atom2.inducedDipolePolar.z*guz24)
                            + atom2.inducedDipole.x*(atom1.inducedDipolePolar.x*gux22 + atom1.inducedDipolePolar.y*gux23 + atom1.inducedDipolePolar.z*gux24)
                            + atom2.inducedDipole.y*(atom1.inducedDipolePolar.x*guy22 + atom1.inducedDipolePolar.y*guy23 + atom1.inducedDipolePolar.z*guy24)
                            + atom2.inducedDipole.z*(atom1.inducedDipolePolar.x*guz22 + atom1.inducedDipolePolar.y*guz23 + atom1.inducedDipolePolar.z*guz24));
#endif
#endif
    real gqxx22 = xr*(2*b20 + xr2*b21);
    real gqxx23 = yr*xr2*b21;
    real gqxx24 = zr*xr2*b21;
    real gqyy22 = xr*yr2*b21;
    real gqyy23 = yr*(2*b20 + yr2*b21);
    real gqyy24 = zr*yr2*b21;
    real gqzz22 = xr*zr2*b21;
    real gqzz23 = yr*zr2*b21;
    real gqzz24 = zr*(2*b20 + zr2*b21);
    real gqxy22 = yr*(b20 + xr2*b21);
    real gqxy23 = xr*(b20 + yr2*b21);
    real gqxy24 = zr*xr*yr*b21;
    real gqxz22 = zr*(b20 + xr2*b21);
    real gqxz23 = gqxy24;
    real gqxz24 = xr*(b20 + zr2*b21);
    real gqyz22 = gqxy24;
    real gqyz23 = zr*(b20 + yr2*b21);
    real gqyz24 = yr*(b20 + zr2*b21);
#if defined B1
    dsumdrB1 += atom2.dipole.x*(atom1.quadrupoleXX*gqxx22 + atom1.quadrupoleYY*gqyy22 + atom1.quadrupoleZZ*gqzz22 + 2*(atom1.quadrupoleXY*gqxy22 + atom1.quadrupoleXZ*gqxz22 + atom1.quadrupoleYZ*gqyz22)) +
                atom2.dipole.y*(atom1.quadrupoleXX*gqxx23 + atom1.quadrupoleYY*gqyy23 + atom1.quadrupoleZZ*gqzz23 + 2*(atom1.quadrupoleXY*gqxy23 + atom1.quadrupoleXZ*gqxz23 + atom1.quadrupoleYZ*gqyz23)) +
                atom2.dipole.z*(atom1.quadrupoleXX*gqxx24 + atom1.quadrupoleYY*gqyy24 + atom1.quadrupoleZZ*gqzz24 + 2*(atom1.quadrupoleXY*gqxy24 + atom1.quadrupoleXZ*gqxz24 + atom1.quadrupoleYZ*gqyz24));
    dsumdrB1 -= atom1.dipole.x*(atom2.quadrupoleXX*gqxx22 + atom2.quadrupoleYY*gqyy22 + atom2.quadrupoleZZ*gqzz22 + 2*(atom2.quadrupoleXY*gqxy22 + atom2.quadrupoleXZ*gqxz22 + atom2.quadrupoleYZ*gqyz22)) +
                atom1.dipole.y*(atom2.quadrupoleXX*gqxx23 + atom2.quadrupoleYY*gqyy23 + atom2.quadrupoleZZ*gqzz23 + 2*(atom2.quadrupoleXY*gqxy23 + atom2.quadrupoleXZ*gqxz23 + atom2.quadrupoleYZ*gqyz23)) +
                atom1.dipole.z*(atom2.quadrupoleXX*gqxx24 + atom2.quadrupoleYY*gqyy24 + atom2.quadrupoleZZ*gqzz24 + 2*(atom2.quadrupoleXY*gqxy24 + atom2.quadrupoleXZ*gqxz24 + atom2.quadrupoleYZ*gqyz24));
#endif
#if defined B2

    dsumdrB2 += sxk*(atom1.quadrupoleXX*gqxx22 + atom1.quadrupoleYY*gqyy22 + atom1.quadrupoleZZ*gqzz22 + 2*(atom1.quadrupoleXY*gqxy22 + atom1.quadrupoleXZ*gqxz22 + atom1.quadrupoleYZ*gqyz22)) +
                syk*(atom1.quadrupoleXX*gqxx23 + atom1.quadrupoleYY*gqyy23 + atom1.quadrupoleZZ*gqzz23 + 2*(atom1.quadrupoleXY*gqxy23 + atom1.quadrupoleXZ*gqxz23 + atom1.quadrupoleYZ*gqyz23)) +
                szk*(atom1.quadrupoleXX*gqxx24 + atom1.quadrupoleYY*gqyy24 + atom1.quadrupoleZZ*gqzz24 + 2*(atom1.quadrupoleXY*gqxy24 + atom1.quadrupoleXZ*gqxz24 + atom1.quadrupoleYZ*gqyz24));
    dsumdrB2 -= sxi*(atom2.quadrupoleXX*gqxx22 + atom2.quadrupoleYY*gqyy22 + atom2.quadrupoleZZ*gqzz22 + 2*(atom2.quadrupoleXY*gqxy22 + atom2.quadrupoleXZ*gqxz22 + atom2.quadrupoleYZ*gqyz22)) +
                syi*(atom2.quadrupoleXX*gqxx23 + atom2.quadrupoleYY*gqyy23 + atom2.quadrupoleZZ*gqzz23 + 2*(atom2.quadrupoleXY*gqxy23 + atom2.quadrupoleXZ*gqxz23 + atom2.quadrupoleYZ*gqyz23)) +
                szi*(atom2.quadrupoleXX*gqxx24 + atom2.quadrupoleYY*gqyy24 + atom2.quadrupoleZZ*gqzz24 + 2*(atom2.quadrupoleXY*gqxy24 + atom2.quadrupoleXZ*gqxz24 + atom2.quadrupoleYZ*gqyz24));

#endif
#endif

    // unweighted 2nd reaction potential gradient tensor;

#if defined F1 || defined F2 || defined T1
    real gc5 = a01 + xr2*a02;
    real gc6 = xr*yr*a02;
    real gc7 = xr*zr*a02;
    real gc8 = a01 + yr2*a02;
    real gc9 = yr*zr*a02;
    real gc10 = a01 + zr2*a02;
#if defined F1
    energy += atom1.q*(atom2.quadrupoleXX*gc5 + atom2.quadrupoleYY*gc8 + atom2.quadrupoleZZ*gc10 + 2*(atom2.quadrupoleXY*gc6 + atom2.quadrupoleXZ*gc7 + atom2.quadrupoleYZ*gc9));
    energy += atom2.q*(atom1.quadrupoleXX*gc5 + atom1.quadrupoleYY*gc8 + atom1.quadrupoleZZ*gc10 + 2*(atom1.quadrupoleXY*gc6 + atom1.quadrupoleXZ*gc7 + atom1.quadrupoleYZ*gc9));

    dedx += atom1.q*(atom2.dipole.x*gc5 + atom2.dipole.y*gc6 + atom2.dipole.z*gc7);
    dedx -= atom2.q*(atom1.dipole.x*gc5 + atom1.dipole.y*gc6 + atom1.dipole.z*gc7);

    dedy += atom1.q*(atom2.dipole.x*gc6 + atom2.dipole.y*gc8 + atom2.dipole.z*gc9);
    dedy -= atom2.q*(atom1.dipole.x*gc6 + atom1.dipole.y*gc8 + atom1.dipole.z*gc9);

    dedz += atom1.q*(atom2.dipole.x*gc7 + atom2.dipole.y*gc9 + atom2.dipole.z*gc10);
    dedz -= atom2.q*(atom1.dipole.x*gc7 + atom1.dipole.y*gc9 + atom1.dipole.z*gc10);
#endif

#if defined F2
    dpdx += atom1.q*(sxk*gc5 + syk*gc6 + szk*gc7);
    dpdx -= atom2.q*(sxi*gc5 + syi*gc6 + szi*gc7);
    dpdy += atom1.q*(sxk*gc6 + syk*gc8 + szk*gc9);
    dpdy -= atom2.q*(sxi*gc6 + syi*gc8 + szi*gc9);
    dpdz += atom1.q*(sxk*gc7 + syk*gc9 + szk*gc10);
    dpdz -= atom2.q*(sxi*gc7 + syi*gc9 + szi*gc10);
#endif

#endif

#if defined F1 || defined F2 || defined T1 || defined T2
    real gux5 = xr*(3*a11 + xr2*a12);
    real gux6 = yr*(a11 + xr2*a12);
    real gux7 = zr*(a11 + xr2*a12);
    real gux8 = xr*(a11 + yr2*a12);
    real gux9 = zr*xr*yr*a12;
    real gux10 = xr*(a11 + zr2*a12);
    real guy5 = yr*(a11 + xr2*a12);
    real guy6 = xr*(a11 + yr2*a12);
    real guy7 = gux9;
    real guy8 = yr*(3*a11 + yr2*a12);
    real guy9 = zr*(a11 + yr2*a12);
    real guy10 = yr*(a11 + zr2*a12);
    real guz5 = zr*(a11 + xr2*a12);
    real guz6 = gux9;
    real guz7 = xr*(a11 + zr2*a12);
    real guz8 = zr*(a11 + yr2*a12);
    real guz9 = yr*(a11 + zr2*a12);
    real guz10 = zr*(3*a11 + zr2*a12);
#if defined F1
    energy -= atom1.dipole.x*(atom2.quadrupoleXX*gux5 + atom2.quadrupoleYY*gux8 + atom2.quadrupoleZZ*gux10 + 2*(atom2.quadrupoleXY*gux6 + atom2.quadrupoleXZ*gux7 + atom2.quadrupoleYZ*gux9)) +
                       atom1.dipole.y*(atom2.quadrupoleXX*guy5 + atom2.quadrupoleYY*guy8 + atom2.quadrupoleZZ*guy10 + 2*(atom2.quadrupoleXY*guy6 + atom2.quadrupoleXZ*guy7 + atom2.quadrupoleYZ*guy9)) +
                       atom1.dipole.z*(atom2.quadrupoleXX*guz5 + atom2.quadrupoleYY*guz8 + atom2.quadrupoleZZ*guz10 + 2*(atom2.quadrupoleXY*guz6 + atom2.quadrupoleXZ*guz7 + atom2.quadrupoleYZ*guz9));

    energy += atom2.dipole.x*(atom1.quadrupoleXX*gux5 + atom1.quadrupoleYY*gux8 + atom1.quadrupoleZZ*gux10 + 2*(atom1.quadrupoleXY*gux6 + atom1.quadrupoleXZ*gux7 + atom1.quadrupoleYZ*gux9)) +
                       atom2.dipole.y*(atom1.quadrupoleXX*guy5 + atom1.quadrupoleYY*guy8 + atom1.quadrupoleZZ*guy10 + 2*(atom1.quadrupoleXY*guy6 + atom1.quadrupoleXZ*guy7 + atom1.quadrupoleYZ*guy9)) +
                       atom2.dipole.z*(atom1.quadrupoleXX*guz5 + atom1.quadrupoleYY*guz8 + atom1.quadrupoleZZ*guz10 + 2*(atom1.quadrupoleXY*guz6 + atom1.quadrupoleXZ*guz7 + atom1.quadrupoleYZ*guz9));

    dedx -= 2*(atom1.dipole.x*(atom2.dipole.x*gux5 + atom2.dipole.y*guy5 + atom2.dipole.z*guz5) +
                              atom1.dipole.y*(atom2.dipole.x*gux6 + atom2.dipole.y*guy6 + atom2.dipole.z*guz6) +
                              atom1.dipole.z*(atom2.dipole.x*gux7 + atom2.dipole.y*guy7 + atom2.dipole.z*guz7));

    dedy -= 2*(atom1.dipole.x*(atom2.dipole.x*gux6 + atom2.dipole.y*guy6 + atom2.dipole.z*guz6) +
                              atom1.dipole.y*(atom2.dipole.x*gux8 + atom2.dipole.y*guy8 + atom2.dipole.z*guz8) +
                              atom1.dipole.z*(atom2.dipole.x*gux9 + atom2.dipole.y*guy9 + atom2.dipole.z*guz9));

    dedz -= 2*(atom1.dipole.x*(atom2.dipole.x*gux7 + atom2.dipole.y*guy7 + atom2.dipole.z*guz7) +
                              atom1.dipole.y*(atom2.dipole.x*gux9 + atom2.dipole.y*guy9 + atom2.dipole.z*guz9) +
                              atom1.dipole.z*(atom2.dipole.x*gux10 + atom2.dipole.y*guy10 + atom2.dipole.z*guz10));

#endif

#if defined F2
    energy -= atom1.inducedDipole.x*(atom2.quadrupoleXX*gux5 + atom2.quadrupoleYY*gux8 + atom2.quadrupoleZZ*gux10 + 2*(atom2.quadrupoleXY*gux6 + atom2.quadrupoleXZ*gux7 + atom2.quadrupoleYZ*gux9)) +
              atom1.inducedDipole.y*(atom2.quadrupoleXX*guy5 + atom2.quadrupoleYY*guy8 + atom2.quadrupoleZZ*guy10 + 2*(atom2.quadrupoleXY*guy6 + atom2.quadrupoleXZ*guy7 + atom2.quadrupoleYZ*guy9)) +
              atom1.inducedDipole.z*(atom2.quadrupoleXX*guz5 + atom2.quadrupoleYY*guz8 + atom2.quadrupoleZZ*guz10 + 2*(atom2.quadrupoleXY*guz6 + atom2.quadrupoleXZ*guz7 + atom2.quadrupoleYZ*guz9));

    energy += atom2.inducedDipole.x*(atom1.quadrupoleXX*gux5 + atom1.quadrupoleYY*gux8 + atom1.quadrupoleZZ*gux10 + 2*(atom1.quadrupoleXY*gux6 + atom1.quadrupoleXZ*gux7 + atom1.quadrupoleYZ*gux9)) +
              atom2.inducedDipole.y*(atom1.quadrupoleXX*guy5 + atom1.quadrupoleYY*guy8 + atom1.quadrupoleZZ*guy10 + 2*(atom1.quadrupoleXY*guy6 + atom1.quadrupoleXZ*guy7 + atom1.quadrupoleYZ*guy9)) +
              atom2.inducedDipole.z*(atom1.quadrupoleXX*guz5 + atom1.quadrupoleYY*guz8 + atom1.quadrupoleZZ*guz10 + 2*(atom1.quadrupoleXY*guz6 + atom1.quadrupoleXZ*guz7 + atom1.quadrupoleYZ*guz9));

    dpdx -= 2*(atom1.dipole.x*(sxk*gux5 + syk*guy5 + szk*guz5) + atom1.dipole.y*(sxk*gux6 + syk*guy6 + szk*guz6) + atom1.dipole.z*(sxk*gux7 + syk*guy7 + szk*guz7) +
                    atom2.dipole.x*(sxi*gux5 + syi*guy5 + szi*guz5) + atom2.dipole.y*(sxi*gux6 + syi*guy6 + szi*guz6) + atom2.dipole.z*(sxi*gux7 + syi*guy7 + szi*guz7));

    dpdy -= 2*(atom1.dipole.x*(sxk*gux6 + syk*guy6 + szk*guz6) + atom1.dipole.y*(sxk*gux8 + syk*guy8 + szk*guz8) + atom1.dipole.z*(sxk*gux9 + syk*guy9 + szk*guz9) +
                    atom2.dipole.x*(sxi*gux6 + syi*guy6 + szi*guz6) + atom2.dipole.y*(sxi*gux8 + syi*guy8 + szi*guz8) + atom2.dipole.z*(sxi*gux9 + syi*guy9 + szi*guz9));

    dpdz -= 2*(atom1.dipole.x*(sxk*gux7 + syk*guy7 + szk*guz7) + atom1.dipole.y*(sxk*gux9 + syk*guy9 + szk*guz9) + atom1.dipole.z*(sxk*gux10 + syk*guy10 + szk*guz10) +
                    atom2.dipole.x*(sxi*gux7 + syi*guy7 + szi*guz7) + atom2.dipole.y*(sxi*gux9 + syi*guy9 + szi*guz9) + atom2.dipole.z*(sxi*gux10 + syi*guy10 + szi*guz10));

#ifndef DIRECT_POLARIZATION
        dpdx -= 2*(atom1.inducedDipole.x*(atom2.inducedDipolePolar.x*gux5 + atom2.inducedDipolePolar.y*gux6 + atom2.inducedDipolePolar.z*gux7)
                            + atom1.inducedDipole.y*(atom2.inducedDipolePolar.x*guy5 + atom2.inducedDipolePolar.y*guy6 + atom2.inducedDipolePolar.z*guy7)
                            + atom1.inducedDipole.z*(atom2.inducedDipolePolar.x*guz5 + atom2.inducedDipolePolar.y*guz6 + atom2.inducedDipolePolar.z*guz7)
                            + atom2.inducedDipole.x*(atom1.inducedDipolePolar.x*gux5 + atom1.inducedDipolePolar.y*gux6 + atom1.inducedDipolePolar.z*gux7)
                            + atom2.inducedDipole.y*(atom1.inducedDipolePolar.x*guy5 + atom1.inducedDipolePolar.y*guy6 + atom1.inducedDipolePolar.z*guy7)
                            + atom2.inducedDipole.z*(atom1.inducedDipolePolar.x*guz5 + atom1.inducedDipolePolar.y*guz6 + atom1.inducedDipolePolar.z*guz7));

        dpdy -= 2*(atom1.inducedDipole.x*(atom2.inducedDipolePolar.x*gux6 + atom2.inducedDipolePolar.y*gux8 + atom2.inducedDipolePolar.z*gux9)
                            + atom1.inducedDipole.y*(atom2.inducedDipolePolar.x*guy6 + atom2.inducedDipolePolar.y*guy8 + atom2.inducedDipolePolar.z*guy9)
                            + atom1.inducedDipole.z*(atom2.inducedDipolePolar.x*guz6 + atom2.inducedDipolePolar.y*guz8 + atom2.inducedDipolePolar.z*guz9)
                            + atom2.inducedDipole.x*(atom1.inducedDipolePolar.x*gux6 + atom1.inducedDipolePolar.y*gux8 + atom1.inducedDipolePolar.z*gux9)
                            + atom2.inducedDipole.y*(atom1.inducedDipolePolar.x*guy6 + atom1.inducedDipolePolar.y*guy8 + atom1.inducedDipolePolar.z*guy9)
                            + atom2.inducedDipole.z*(atom1.inducedDipolePolar.x*guz6 + atom1.inducedDipolePolar.y*guz8 + atom1.inducedDipolePolar.z*guz9));

        dpdz -= 2*(atom1.inducedDipole.x*(atom2.inducedDipolePolar.x*gux7 + atom2.inducedDipolePolar.y*gux9 + atom2.inducedDipolePolar.z*gux10)
                            + atom1.inducedDipole.y*(atom2.inducedDipolePolar.x*guy7 + atom2.inducedDipolePolar.y*guy9 + atom2.inducedDipolePolar.z*guy10)
                            + atom1.inducedDipole.z*(atom2.inducedDipolePolar.x*guz7 + atom2.inducedDipolePolar.y*guz9 + atom2.inducedDipolePolar.z*guz10)
                            + atom2.inducedDipole.x*(atom1.inducedDipolePolar.x*gux7 + atom1.inducedDipolePolar.y*gux9 + atom1.inducedDipolePolar.z*gux10)
                            + atom2.inducedDipole.y*(atom1.inducedDipolePolar.x*guy7 + atom1.inducedDipolePolar.y*guy9 + atom1.inducedDipolePolar.z*guy10)
                            + atom2.inducedDipole.z*(atom1.inducedDipolePolar.x*guz7 + atom1.inducedDipolePolar.y*guz9 + atom1.inducedDipolePolar.z*guz10));
#endif
#endif
#endif

#if defined F1 || defined F2 || defined T1
    real gqxx5 = 2*a20 + xr2*(5*a21 + xr2*a22);
    real gqxx6 = yr*xr*(2*a21 + xr2*a22);
    real gqxx7 = zr*xr*(2*a21 + xr2*a22);
    real gqxx8 = xr2*(a21 + yr2*a22);
    real gqxx9 = zr*yr*xr2*a22;
    real gqxx10 = xr2*(a21 + zr2*a22);
    real gqyy5 = yr2*(a21 + xr2*a22);
    real gqyy6 = xr*yr*(2*a21 + yr2*a22);
    real gqyy7 = xr*zr*yr2*a22;
    real gqyy8 = 2*a20 + yr2*(5*a21 + yr2*a22);
    real gqyy9 = yr*zr*(2*a21 + yr2*a22);
    real gqyy10 = yr2*(a21 + zr2*a22);
    real gqzz5 = zr2*(a21 + xr2*a22);
    real gqzz6 = xr*yr*zr2*a22;
    real gqzz7 = xr*zr*(2*a21 + zr2*a22);
    real gqzz8 = zr2*(a21 + yr2*a22);
    real gqzz9 = yr*zr*(2*a21 + zr2*a22);
    real gqzz10 = 2*a20 + zr2*(5*a21 + zr2*a22);
    real gqxy5 = xr*yr*(3*a21 + xr2*a22);
    real gqxy6 = a20 + (xr2 + yr2)*a21 + xr2*yr2*a22;
    real gqxy7 = zr*yr*(a21 + xr2*a22);
    real gqxy8 = xr*yr*(3*a21 + yr2*a22);
    real gqxy9 = zr*xr*(a21 + yr2*a22);
    real gqxy10 = xr*yr*(a21 + zr2*a22);
    real gqxz5 = xr*zr*(3*a21 + xr2*a22);
    real gqxz6 = yr*zr*(a21 + xr2*a22);
    real gqxz7 = a20 + (xr2 + zr2)*a21 + xr2*zr2*a22;
    real gqxz8 = xr*zr*(a21 + yr2*a22);
    real gqxz9 = xr*yr*(a21 + zr2*a22);
    real gqxz10 = xr*zr*(3*a21 + zr2*a22);
    real gqyz5 = zr*yr*(a21 + xr2*a22);
    real gqyz6 = xr*zr*(a21 + yr2*a22);
    real gqyz7 = xr*yr*(a21 + zr2*a22);
    real gqyz8 = yr*zr*(3*a21 + yr2*a22);
    real gqyz9 = a20 + (yr2 + zr2)*a21 + yr2*zr2*a22;
    real gqyz10 = yr*zr*(3*a21 + zr2*a22);
#if defined F1
    energy += atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqxx8 + atom2.quadrupoleZZ*gqxx10 + 2*(atom2.quadrupoleXY*gqxx6 + atom2.quadrupoleXZ*gqxx7 + atom2.quadrupoleYZ*gqxx9))
              + atom1.quadrupoleYY*(atom2.quadrupoleXX*gqyy5 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqyy10 + 2*(atom2.quadrupoleXY*gqyy6 + atom2.quadrupoleXZ*gqyy7 + atom2.quadrupoleYZ*gqyy9))
              + atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqzz5 + atom2.quadrupoleYY*gqzz8 + atom2.quadrupoleZZ*gqzz10 + 2*(atom2.quadrupoleXY*gqzz6 + atom2.quadrupoleXZ*gqzz7 + atom2.quadrupoleYZ*gqzz9))
              + 2*(atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxy5 + atom2.quadrupoleYY*gqxy8 + atom2.quadrupoleZZ*gqxy10
              + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxy7 + atom2.quadrupoleYZ*gqxy9))
              + atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxz5 + atom2.quadrupoleYY*gqxz8 + atom2.quadrupoleZZ*gqxz10
              + 2*(atom2.quadrupoleXY*gqxz6 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqxz9))
              + atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqyz5 + atom2.quadrupoleYY*gqyz8 + atom2.quadrupoleZZ*gqyz10
              + 2*(atom2.quadrupoleXY*gqyz6 + atom2.quadrupoleXZ*gqyz7 + atom2.quadrupoleYZ*gqyz9)));

    energy += atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqyy5 + atom2.quadrupoleZZ*gqzz5 + 2*(atom2.quadrupoleXY*gqxy5 + atom2.quadrupoleXZ*gqxz5 + atom2.quadrupoleYZ*gqyz5))
             + atom1.quadrupoleYY*(atom2.quadrupoleXX*gqxx8 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqzz8
             + 2*(atom2.quadrupoleXY*gqxy8 + atom2.quadrupoleXZ*gqxz8 + atom2.quadrupoleYZ*gqyz8))
             + atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqxx10 + atom2.quadrupoleYY*gqyy10 + atom2.quadrupoleZZ*gqzz10
             + 2*(atom2.quadrupoleXY*gqxy10 + atom2.quadrupoleXZ*gqxz10 + atom2.quadrupoleYZ*gqyz10))
             + 2*(atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6
             + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6))
             + atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7
             + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7))
             + atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9
             + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9)));

    dedx += atom2.dipole.x*(atom1.quadrupoleXX*gqxx5 + atom1.quadrupoleYY*gqyy5 + atom1.quadrupoleZZ*gqzz5 + 2*(atom1.quadrupoleXY*gqxy5 + atom1.quadrupoleXZ*gqxz5 + atom1.quadrupoleYZ*gqyz5)) +
              atom2.dipole.y*(atom1.quadrupoleXX*gqxx6 + atom1.quadrupoleYY*gqyy6 + atom1.quadrupoleZZ*gqzz6 + 2*(atom1.quadrupoleXY*gqxy6 + atom1.quadrupoleXZ*gqxz6 + atom1.quadrupoleYZ*gqyz6)) +
              atom2.dipole.z*(atom1.quadrupoleXX*gqxx7 + atom1.quadrupoleYY*gqyy7 + atom1.quadrupoleZZ*gqzz7 + 2*(atom1.quadrupoleXY*gqxy7 + atom1.quadrupoleXZ*gqxz7 + atom1.quadrupoleYZ*gqyz7));

    dedx -= atom1.dipole.x*(atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqyy5 + atom2.quadrupoleZZ*gqzz5 + 2*(atom2.quadrupoleXY*gqxy5 + atom2.quadrupoleXZ*gqxz5 + atom2.quadrupoleYZ*gqyz5)) +
              atom1.dipole.y*(atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6 + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6)) +
              atom1.dipole.z*(atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7 + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7));

    dedy += atom2.dipole.x*(atom1.quadrupoleXX*gqxx6 + atom1.quadrupoleYY*gqyy6 + atom1.quadrupoleZZ*gqzz6 + 2*(atom1.quadrupoleXY*gqxy6 + atom1.quadrupoleXZ*gqxz6 + atom1.quadrupoleYZ*gqyz6)) +
              atom2.dipole.y*(atom1.quadrupoleXX*gqxx8 + atom1.quadrupoleYY*gqyy8 + atom1.quadrupoleZZ*gqzz8 + 2*(atom1.quadrupoleXY*gqxy8 + atom1.quadrupoleXZ*gqxz8 + atom1.quadrupoleYZ*gqyz8)) +
              atom2.dipole.z*(atom1.quadrupoleXX*gqxx9 + atom1.quadrupoleYY*gqyy9 + atom1.quadrupoleZZ*gqzz9 + 2*(atom1.quadrupoleXY*gqxy9 + atom1.quadrupoleXZ*gqxz9 + atom1.quadrupoleYZ*gqyz9));

    dedy -= atom1.dipole.x*(atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6 + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6)) +
              atom1.dipole.y*(atom2.quadrupoleXX*gqxx8 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqzz8 + 2*(atom2.quadrupoleXY*gqxy8 + atom2.quadrupoleXZ*gqxz8 + atom2.quadrupoleYZ*gqyz8)) +
              atom1.dipole.z*(atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9 + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9));

    dedz += atom2.dipole.x*(atom1.quadrupoleXX*gqxx7 + atom1.quadrupoleYY*gqyy7 + atom1.quadrupoleZZ*gqzz7 + 2*(atom1.quadrupoleXY*gqxy7 + atom1.quadrupoleXZ*gqxz7 + atom1.quadrupoleYZ*gqyz7)) +
              atom2.dipole.y*(atom1.quadrupoleXX*gqxx9 + atom1.quadrupoleYY*gqyy9 + atom1.quadrupoleZZ*gqzz9 + 2*(atom1.quadrupoleXY*gqxy9 + atom1.quadrupoleXZ*gqxz9 + atom1.quadrupoleYZ*gqyz9)) +
              atom2.dipole.z*(atom1.quadrupoleXX*gqxx10 + atom1.quadrupoleYY*gqyy10 + atom1.quadrupoleZZ*gqzz10 + 2*(atom1.quadrupoleXY*gqxy10 + atom1.quadrupoleXZ*gqxz10 + atom1.quadrupoleYZ*gqyz10));

    dedz -= atom1.dipole.x*(atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7 + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7)) +
              atom1.dipole.y*(atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9 + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9)) +
              atom1.dipole.z*(atom2.quadrupoleXX*gqxx10 + atom2.quadrupoleYY*gqyy10 + atom2.quadrupoleZZ*gqzz10 + 2*(atom2.quadrupoleXY*gqxy10 + atom2.quadrupoleXZ*gqxz10 + atom2.quadrupoleYZ*gqyz10));
#endif
#if defined F2
    dpdx += sxk*(atom1.quadrupoleXX*gqxx5 + atom1.quadrupoleYY*gqyy5 + atom1.quadrupoleZZ*gqzz5 + 2*(atom1.quadrupoleXY*gqxy5 + atom1.quadrupoleXZ*gqxz5 + atom1.quadrupoleYZ*gqyz5)) +
            syk*(atom1.quadrupoleXX*gqxx6 + atom1.quadrupoleYY*gqyy6 + atom1.quadrupoleZZ*gqzz6 + 2*(atom1.quadrupoleXY*gqxy6 + atom1.quadrupoleXZ*gqxz6 + atom1.quadrupoleYZ*gqyz6)) +
            szk*(atom1.quadrupoleXX*gqxx7 + atom1.quadrupoleYY*gqyy7 + atom1.quadrupoleZZ*gqzz7 + 2*(atom1.quadrupoleXY*gqxy7 + atom1.quadrupoleXZ*gqxz7 + atom1.quadrupoleYZ*gqyz7));

    dpdx -= sxi*(atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqyy5 + atom2.quadrupoleZZ*gqzz5 + 2*(atom2.quadrupoleXY*gqxy5 + atom2.quadrupoleXZ*gqxz5 + atom2.quadrupoleYZ*gqyz5)) +
            syi*(atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6 + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6)) +
            szi*(atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7 + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7));

    dpdy += sxk*(atom1.quadrupoleXX*gqxx6 + atom1.quadrupoleYY*gqyy6 + atom1.quadrupoleZZ*gqzz6 + 2*(atom1.quadrupoleXY*gqxy6 + atom1.quadrupoleXZ*gqxz6 + atom1.quadrupoleYZ*gqyz6)) +
            syk*(atom1.quadrupoleXX*gqxx8 + atom1.quadrupoleYY*gqyy8 + atom1.quadrupoleZZ*gqzz8 + 2*(atom1.quadrupoleXY*gqxy8 + atom1.quadrupoleXZ*gqxz8 + atom1.quadrupoleYZ*gqyz8)) +
            szk*(atom1.quadrupoleXX*gqxx9 + atom1.quadrupoleYY*gqyy9 + atom1.quadrupoleZZ*gqzz9 + 2*(atom1.quadrupoleXY*gqxy9 + atom1.quadrupoleXZ*gqxz9 + atom1.quadrupoleYZ*gqyz9));

    dpdy -= sxi*(atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6 + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6)) +
            syi*(atom2.quadrupoleXX*gqxx8 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqzz8 + 2*(atom2.quadrupoleXY*gqxy8 + atom2.quadrupoleXZ*gqxz8 + atom2.quadrupoleYZ*gqyz8)) +
            szi*(atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9 + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9));

    dpdz -= sxi*(atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7 + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7)) +
            syi*(atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9 + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9)) +
            szi*(atom2.quadrupoleXX*gqxx10 + atom2.quadrupoleYY*gqyy10 + atom2.quadrupoleZZ*gqzz10 + 2*(atom2.quadrupoleXY*gqxy10 + atom2.quadrupoleXZ*gqxz10 + atom2.quadrupoleYZ*gqyz10));

    dpdz += sxk*(atom1.quadrupoleXX*gqxx7 + atom1.quadrupoleYY*gqyy7 + atom1.quadrupoleZZ*gqzz7 + 2*(atom1.quadrupoleXY*gqxy7 + atom1.quadrupoleXZ*gqxz7 + atom1.quadrupoleYZ*gqyz7)) +
            syk*(atom1.quadrupoleXX*gqxx9 + atom1.quadrupoleYY*gqyy9 + atom1.quadrupoleZZ*gqzz9 + 2*(atom1.quadrupoleXY*gqxy9 + atom1.quadrupoleXZ*gqxz9 + atom1.quadrupoleYZ*gqyz9)) +
            szk*(atom1.quadrupoleXX*gqxx10 + atom1.quadrupoleYY*gqyy10 + atom1.quadrupoleZZ*gqzz10 + 2*(atom1.quadrupoleXY*gqxy10 + atom1.quadrupoleXZ*gqxz10 + atom1.quadrupoleYZ*gqyz10));
#endif
#endif

    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

#if defined B1
    real gc25 = b01 + xr2*b02;
    real gc26 = xr*yr*b02;
    real gc27 = xr*zr*b02;
    real gc28 = b01 + yr2*b02;
    real gc29 = yr*zr*b02;
    real gc30 = b01 + zr2*b02;
    dsumdrB1 += atom1.q*(atom2.quadrupoleXX*gc25 + atom2.quadrupoleYY*gc28 + atom2.quadrupoleZZ*gc30 + 2*(atom2.quadrupoleXY*gc26 + atom2.quadrupoleXZ*gc27 + atom2.quadrupoleYZ*gc29));
    dsumdrB1 += atom2.q*(atom1.quadrupoleXX*gc25 + atom1.quadrupoleYY*gc28 + atom1.quadrupoleZZ*gc30 + 2*(atom1.quadrupoleXY*gc26 + atom1.quadrupoleXZ*gc27 + atom1.quadrupoleYZ*gc29));
#endif
#if defined B1 || defined B2
    real gux25 = xr*(3*b11 + xr2*b12);
    real gux26 = yr*(b11 + xr2*b12);
    real gux27 = zr*(b11 + xr2*b12);
    real gux28 = xr*(b11 + yr2*b12);
    real gux29 = zr*xr*yr*b12;
    real gux30 = xr*(b11 + zr2*b12);
    real guy25 = yr*(b11 + xr2*b12);
    real guy26 = xr*(b11 + yr2*b12);
    real guy27 = gux29;
    real guy28 = yr*(3*b11 + yr2*b12);
    real guy29 = zr*(b11 + yr2*b12);
    real guy30 = yr*(b11 + zr2*b12);
    real guz25 = zr*(b11 + xr2*b12);
    real guz26 = gux29;
    real guz27 = xr*(b11 + zr2*b12);
    real guz28 = zr*(b11 + yr2*b12);
    real guz29 = yr*(b11 + zr2*b12);
    real guz30 = zr*(3*b11 + zr2*b12);
#endif
#if defined B2
    dsumdrB2 -= sxi*(atom2.quadrupoleXX*gux25 + atom2.quadrupoleYY*gux28 + atom2.quadrupoleZZ*gux30 + 2*(atom2.quadrupoleXY*gux26 + atom2.quadrupoleXZ*gux27 + atom2.quadrupoleYZ*gux29)) +
                syi*(atom2.quadrupoleXX*guy25 + atom2.quadrupoleYY*guy28 + atom2.quadrupoleZZ*guy30 + 2*(atom2.quadrupoleXY*guy26 + atom2.quadrupoleXZ*guy27 + atom2.quadrupoleYZ*guy29)) +
                szi*(atom2.quadrupoleXX*guz25 + atom2.quadrupoleYY*guz28 + atom2.quadrupoleZZ*guz30 + 2*(atom2.quadrupoleXY*guz26 + atom2.quadrupoleXZ*guz27 + atom2.quadrupoleYZ*guz29));
    dsumdrB2 += sxk*(atom1.quadrupoleXX*gux25 + atom1.quadrupoleYY*gux28 + atom1.quadrupoleZZ*gux30 + 2*(atom1.quadrupoleXY*gux26 + atom1.quadrupoleXZ*gux27 + atom1.quadrupoleYZ*gux29)) +
                syk*(atom1.quadrupoleXX*guy25 + atom1.quadrupoleYY*guy28 + atom1.quadrupoleZZ*guy30 + 2*(atom1.quadrupoleXY*guy26 + atom1.quadrupoleXZ*guy27 + atom1.quadrupoleYZ*guy29)) +
                szk*(atom1.quadrupoleXX*guz25 + atom1.quadrupoleYY*guz28 + atom1.quadrupoleZZ*guz30 + 2*(atom1.quadrupoleXY*guz26 + atom1.quadrupoleXZ*guz27 + atom1.quadrupoleYZ*guz29));
#endif
#if defined B1
    dsumdrB1 -= atom1.dipole.x*(atom2.quadrupoleXX*gux25 + atom2.quadrupoleYY*gux28 + atom2.quadrupoleZZ*gux30 + 2*(atom2.quadrupoleXY*gux26 + atom2.quadrupoleXZ*gux27 + atom2.quadrupoleYZ*gux29)) +
                atom1.dipole.y*(atom2.quadrupoleXX*guy25 + atom2.quadrupoleYY*guy28 + atom2.quadrupoleZZ*guy30 + 2*(atom2.quadrupoleXY*guy26 + atom2.quadrupoleXZ*guy27 + atom2.quadrupoleYZ*guy29)) +
                atom1.dipole.z*(atom2.quadrupoleXX*guz25 + atom2.quadrupoleYY*guz28 + atom2.quadrupoleZZ*guz30 + 2*(atom2.quadrupoleXY*guz26 + atom2.quadrupoleXZ*guz27 + atom2.quadrupoleYZ*guz29));
    dsumdrB1 += atom2.dipole.x*(atom1.quadrupoleXX*gux25 + atom1.quadrupoleYY*gux28 + atom1.quadrupoleZZ*gux30 + 2*(atom1.quadrupoleXY*gux26 + atom1.quadrupoleXZ*gux27 + atom1.quadrupoleYZ*gux29)) +
                atom2.dipole.y*(atom1.quadrupoleXX*guy25 + atom1.quadrupoleYY*guy28 + atom1.quadrupoleZZ*guy30 + 2*(atom1.quadrupoleXY*guy26 + atom1.quadrupoleXZ*guy27 + atom1.quadrupoleYZ*guy29)) +
                atom2.dipole.z*(atom1.quadrupoleXX*guz25 + atom1.quadrupoleYY*guz28 + atom1.quadrupoleZZ*guz30 + 2*(atom1.quadrupoleXY*guz26 + atom1.quadrupoleXZ*guz27 + atom1.quadrupoleYZ*guz29));

    real gqxx25 = 2*b20 + xr2*(5*b21 + xr2*b22);
    real gqxx26 = yr*xr*(2*b21 + xr2*b22);
    real gqxx27 = zr*xr*(2*b21 + xr2*b22);
    real gqxx28 = xr2*(b21 + yr2*b22);
    real gqxx29 = zr*yr*xr2*b22;
    real gqxx30 = xr2*(b21 + zr2*b22);
    real gqyy25 = yr2*(b21 + xr2*b22);
    real gqyy26 = xr*yr*(2*b21 + yr2*b22);
    real gqyy27 = xr*zr*yr2*b22;
    real gqyy28 = 2*b20 + yr2*(5*b21 + yr2*b22);
    real gqyy29 = yr*zr*(2*b21 + yr2*b22);
    real gqyy30 = yr2*(b21 + zr2*b22);
    real gqzz25 = zr2*(b21 + xr2*b22);
    real gqzz26 = xr*yr*zr2*b22;
    real gqzz27 = xr*zr*(2*b21 + zr2*b22);
    real gqzz28 = zr2*(b21 + yr2*b22);
    real gqzz29 = yr*zr*(2*b21 + zr2*b22);
    real gqzz30 = 2*b20 + zr2*(5*b21 + zr2*b22);
    real gqxy25 = xr*yr*(3*b21 + xr2*b22);
    real gqxy26 = b20 + (xr2 + yr2)*b21 + xr2*yr2*b22;
    real gqxy27 = zr*yr*(b21 + xr2*b22);
    real gqxy28 = xr*yr*(3*b21 + yr2*b22);
    real gqxy29 = zr*xr*(b21 + yr2*b22);
    real gqxy30 = xr*yr*(b21 + zr2*b22);
    real gqxz25 = xr*zr*(3*b21 + xr2*b22);
    real gqxz26 = yr*zr*(b21 + xr2*b22);
    real gqxz27 = b20 + (xr2 + zr2)*b21 + xr2*zr2*b22;
    real gqxz28 = xr*zr*(b21 + yr2*b22);
    real gqxz29 = xr*yr*(b21 + zr2*b22);
    real gqxz30 = xr*zr*(3*b21 + zr2*b22);
    real gqyz25 = zr*yr*(b21 + xr2*b22);
    real gqyz26 = xr*zr*(b21 + yr2*b22);
    real gqyz27 = xr*yr*(b21 + zr2*b22);
    real gqyz28 = yr*zr*(3*b21 + yr2*b22);
    real gqyz29 = b20 + (yr2 + zr2)*b21 + yr2*zr2*b22;
    real gqyz30 = yr*zr*(3*b21 + zr2*b22);

    dsumdrB1 +=
        atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx25 + atom2.quadrupoleYY*gqxx28 + atom2.quadrupoleZZ*gqxx30 + 2*(atom2.quadrupoleXY*gqxx26 + atom2.quadrupoleXZ*gqxx27 + atom2.quadrupoleYZ*gqxx29)) +
        atom1.quadrupoleYY*(atom2.quadrupoleXX*gqyy25 + atom2.quadrupoleYY*gqyy28 + atom2.quadrupoleZZ*gqyy30 + 2*(atom2.quadrupoleXY*gqyy26 + atom2.quadrupoleXZ*gqyy27 + atom2.quadrupoleYZ*gqyy29)) +
        atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqzz25 + atom2.quadrupoleYY*gqzz28 + atom2.quadrupoleZZ*gqzz30 + 2*(atom2.quadrupoleXY*gqzz26 + atom2.quadrupoleXZ*gqzz27 + atom2.quadrupoleYZ*gqzz29));

    dsumdrB1 += 2*(
        atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxy25 + atom2.quadrupoleYY*gqxy28 + atom2.quadrupoleZZ*gqxy30 + 2*(atom2.quadrupoleXY*gqxy26 + atom2.quadrupoleXZ*gqxy27 + atom2.quadrupoleYZ*gqxy29)) +
        atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxz25 + atom2.quadrupoleYY*gqxz28 + atom2.quadrupoleZZ*gqxz30 + 2*(atom2.quadrupoleXY*gqxz26 + atom2.quadrupoleXZ*gqxz27 + atom2.quadrupoleYZ*gqxz29)) +
        atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqyz25 + atom2.quadrupoleYY*gqyz28 + atom2.quadrupoleZZ*gqyz30 + 2*(atom2.quadrupoleXY*gqyz26 + atom2.quadrupoleXZ*gqyz27 + atom2.quadrupoleYZ*gqyz29)));

    dsumdrB1 +=
        atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx25 + atom2.quadrupoleYY*gqyy25 + atom2.quadrupoleZZ*gqzz25 + 2*(atom2.quadrupoleXY*gqxy25 + atom2.quadrupoleXZ*gqxz25 + atom2.quadrupoleYZ*gqyz25)) +
        atom1.quadrupoleYY*(atom2.quadrupoleXX*gqxx28 + atom2.quadrupoleYY*gqyy28 + atom2.quadrupoleZZ*gqzz28 + 2*(atom2.quadrupoleXY*gqxy28 + atom2.quadrupoleXZ*gqxz28 + atom2.quadrupoleYZ*gqyz28)) +
        atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqxx30 + atom2.quadrupoleYY*gqyy30 + atom2.quadrupoleZZ*gqzz30 + 2*(atom2.quadrupoleXY*gqxy30 + atom2.quadrupoleXZ*gqxz30 + atom2.quadrupoleYZ*gqyz30));

    dsumdrB1 += 2*(
        atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxx26 + atom2.quadrupoleYY*gqyy26 + atom2.quadrupoleZZ*gqzz26 + 2*(atom2.quadrupoleXY*gqxy26 + atom2.quadrupoleXZ*gqxz26 + atom2.quadrupoleYZ*gqyz26)) +
        atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxx27 + atom2.quadrupoleYY*gqyy27 + atom2.quadrupoleZZ*gqzz27 + 2*(atom2.quadrupoleXY*gqxy27 + atom2.quadrupoleXZ*gqxz27 + atom2.quadrupoleYZ*gqyz27)) +
        atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqxx29 + atom2.quadrupoleYY*gqyy29 + atom2.quadrupoleZZ*gqzz29 + 2*(atom2.quadrupoleXY*gqxy29 + atom2.quadrupoleXZ*gqxz29 + atom2.quadrupoleYZ*gqyz29)));

    dsumdrB1 *= 0.5f;
    atom1.bornForce += atom2.bornRadius*dsumdrB1;
    atom2.bornForce += atom1.bornRadius*dsumdrB1;
#endif

    // unweighted 3rd reaction potential gradient tensor;

#if defined F1
    real gc11 = xr*(3*a02 + xr2*a03);
    real gc12 = yr*(a02 + xr2*a03);
    real gc13 = zr*(a02 + xr2*a03);
    real gc14 = xr*(a02 + yr2*a03);
    real gc15 = xr*yr*zr*a03;
    real gc16 = xr*(a02 + zr2*a03);
    real gc17 = yr*(3*a02 + yr2*a03);
    real gc18 = zr*(a02 + yr2*a03);
    real gc19 = yr*(a02 + zr2*a03);
    real gc20 = zr*(3*a02 + zr2*a03);
    dedx += atom1.q*(atom2.quadrupoleXX*gc11 + atom2.quadrupoleYY*gc14 + atom2.quadrupoleZZ*gc16 + 2*(atom2.quadrupoleXY*gc12 + atom2.quadrupoleXZ*gc13 + atom2.quadrupoleYZ*gc15));
    dedx += atom2.q*(atom1.quadrupoleXX*gc11 + atom1.quadrupoleYY*gc14 + atom1.quadrupoleZZ*gc16 + 2*(atom1.quadrupoleXY*gc12 + atom1.quadrupoleXZ*gc13 + atom1.quadrupoleYZ*gc15));
    dedy += atom1.q*(atom2.quadrupoleXX*gc12 + atom2.quadrupoleYY*gc17 + atom2.quadrupoleZZ*gc19 + 2*(atom2.quadrupoleXY*gc14 + atom2.quadrupoleXZ*gc15 + atom2.quadrupoleYZ*gc18));
    dedy += atom2.q*(atom1.quadrupoleXX*gc12 + atom1.quadrupoleYY*gc17 + atom1.quadrupoleZZ*gc19 + 2*(atom1.quadrupoleXY*gc14 + atom1.quadrupoleXZ*gc15 + atom1.quadrupoleYZ*gc18));
    dedz += atom1.q*(atom2.quadrupoleXX*gc13 + atom2.quadrupoleYY*gc18 + atom2.quadrupoleZZ*gc20 + 2*(atom2.quadrupoleXY*gc15 + atom2.quadrupoleXZ*gc16 + atom2.quadrupoleYZ*gc19));
    dedz += atom2.q*(atom1.quadrupoleXX*gc13 + atom1.quadrupoleYY*gc18 + atom1.quadrupoleZZ*gc20 + 2*(atom1.quadrupoleXY*gc15 + atom1.quadrupoleXZ*gc16 + atom1.quadrupoleYZ*gc19));
#endif
#if defined F1 || defined F2
    real gux11 = 3*a11 + xr2*(6*a12 + xr2*a13);
    real gux12 = xr*yr*(3*a12 + xr2*a13);
    real gux13 = xr*zr*(3*a12 + xr2*a13);
    real gux14 = a11 + (xr2 + yr2)*a12 + xr2*yr2*a13;
    real gux15 = yr*zr*(a12 + xr2*a13);
    real gux16 = a11 + (xr2 + zr2)*a12 + xr2*zr2*a13;
    real gux17 = xr*yr*(3*a12 + yr2*a13);
    real gux18 = xr*zr*(a12 + yr2*a13);
    real gux19 = xr*yr*(a12 + zr2*a13);
    real gux20 = xr*zr*(3*a12 + zr2*a13);
    real guy11 = gux12;
    real guy12 = gux14;
    real guy13 = gux15;
    real guy14 = gux17;
    real guy15 = gux18;
    real guy16 = gux19;
    real guy17 = 3*a11 + yr2*(6*a12 + yr2*a13);
    real guy18 = yr*zr*(3*a12 + yr2*a13);
    real guy19 = a11 + (yr2 + zr2)*a12 + yr2*zr2*a13;
    real guy20 = yr*zr*(3*a12 + zr2*a13);
    real guz11 = gux13;
    real guz12 = gux15;
    real guz13 = gux16;
    real guz14 = gux18;
    real guz15 = gux19;
    real guz16 = gux20;
    real guz17 = guy18;
    real guz18 = guy19;
    real guz19 = guy20;
    real guz20 = 3*a11 + zr2*(6*a12 + zr2*a13);
#if defined F1
    dedx -= atom1.dipole.x*(atom2.quadrupoleXX*gux11 + atom2.quadrupoleYY*gux14 + atom2.quadrupoleZZ*gux16 + 2*(atom2.quadrupoleXY*gux12 + atom2.quadrupoleXZ*gux13 + atom2.quadrupoleYZ*gux15)) +
                       atom1.dipole.y*(atom2.quadrupoleXX*guy11 + atom2.quadrupoleYY*guy14 + atom2.quadrupoleZZ*guy16 + 2*(atom2.quadrupoleXY*guy12 + atom2.quadrupoleXZ*guy13 + atom2.quadrupoleYZ*guy15)) +
                       atom1.dipole.z*(atom2.quadrupoleXX*guz11 + atom2.quadrupoleYY*guz14 + atom2.quadrupoleZZ*guz16 + 2*(atom2.quadrupoleXY*guz12 + atom2.quadrupoleXZ*guz13 + atom2.quadrupoleYZ*guz15));

    dedx += atom2.dipole.x*(atom1.quadrupoleXX*gux11 + atom1.quadrupoleYY*gux14 + atom1.quadrupoleZZ*gux16 + 2*(atom1.quadrupoleXY*gux12 + atom1.quadrupoleXZ*gux13 + atom1.quadrupoleYZ*gux15)) +
                       atom2.dipole.y*(atom1.quadrupoleXX*guy11 + atom1.quadrupoleYY*guy14 + atom1.quadrupoleZZ*guy16 + 2*(atom1.quadrupoleXY*guy12 + atom1.quadrupoleXZ*guy13 + atom1.quadrupoleYZ*guy15)) +
                       atom2.dipole.z*(atom1.quadrupoleXX*guz11 + atom1.quadrupoleYY*guz14 + atom1.quadrupoleZZ*guz16 + 2*(atom1.quadrupoleXY*guz12 + atom1.quadrupoleXZ*guz13 + atom1.quadrupoleYZ*guz15));

    dedy -= atom1.dipole.x*(atom2.quadrupoleXX*gux12 + atom2.quadrupoleYY*gux17 + atom2.quadrupoleZZ*gux19 + 2*(atom2.quadrupoleXY*gux14 + atom2.quadrupoleXZ*gux15 + atom2.quadrupoleYZ*gux18)) +
                       atom1.dipole.y*(atom2.quadrupoleXX*guy12 + atom2.quadrupoleYY*guy17 + atom2.quadrupoleZZ*guy19 + 2*(atom2.quadrupoleXY*guy14 + atom2.quadrupoleXZ*guy15 + atom2.quadrupoleYZ*guy18)) +
                       atom1.dipole.z*(atom2.quadrupoleXX*guz12 + atom2.quadrupoleYY*guz17 + atom2.quadrupoleZZ*guz19 + 2*(atom2.quadrupoleXY*guz14 + atom2.quadrupoleXZ*guz15 + atom2.quadrupoleYZ*guz18));

    dedy += atom2.dipole.x*(atom1.quadrupoleXX*gux12 + atom1.quadrupoleYY*gux17 + atom1.quadrupoleZZ*gux19 + 2*(atom1.quadrupoleXY*gux14 + atom1.quadrupoleXZ*gux15 + atom1.quadrupoleYZ*gux18)) +
                       atom2.dipole.y*(atom1.quadrupoleXX*guy12 + atom1.quadrupoleYY*guy17 + atom1.quadrupoleZZ*guy19 + 2*(atom1.quadrupoleXY*guy14 + atom1.quadrupoleXZ*guy15 + atom1.quadrupoleYZ*guy18)) +
                       atom2.dipole.z*(atom1.quadrupoleXX*guz12 + atom1.quadrupoleYY*guz17 + atom1.quadrupoleZZ*guz19 + 2*(atom1.quadrupoleXY*guz14 + atom1.quadrupoleXZ*guz15 + atom1.quadrupoleYZ*guz18));

    dedz -= atom1.dipole.x*(atom2.quadrupoleXX*gux13 + atom2.quadrupoleYY*gux18 + atom2.quadrupoleZZ*gux20 + 2*(atom2.quadrupoleXY*gux15 + atom2.quadrupoleXZ*gux16 + atom2.quadrupoleYZ*gux19)) +
                       atom1.dipole.y*(atom2.quadrupoleXX*guy13 + atom2.quadrupoleYY*guy18 + atom2.quadrupoleZZ*guy20 + 2*(atom2.quadrupoleXY*guy15 + atom2.quadrupoleXZ*guy16 + atom2.quadrupoleYZ*guy19)) +
                       atom1.dipole.z*(atom2.quadrupoleXX*guz13 + atom2.quadrupoleYY*guz18 + atom2.quadrupoleZZ*guz20 + 2*(atom2.quadrupoleXY*guz15 + atom2.quadrupoleXZ*guz16 + atom2.quadrupoleYZ*guz19));

    dedz += atom2.dipole.x*(atom1.quadrupoleXX*gux13 + atom1.quadrupoleYY*gux18 + atom1.quadrupoleZZ*gux20 + 2*(atom1.quadrupoleXY*gux15 + atom1.quadrupoleXZ*gux16 + atom1.quadrupoleYZ*gux19)) +
                       atom2.dipole.y*(atom1.quadrupoleXX*guy13 + atom1.quadrupoleYY*guy18 + atom1.quadrupoleZZ*guy20 + 2*(atom1.quadrupoleXY*guy15 + atom1.quadrupoleXZ*guy16 + atom1.quadrupoleYZ*guy19)) +
                       atom2.dipole.z*(atom1.quadrupoleXX*guz13 + atom1.quadrupoleYY*guz18 + atom1.quadrupoleZZ*guz20 + 2*(atom1.quadrupoleXY*guz15 + atom1.quadrupoleXZ*guz16 + atom1.quadrupoleYZ*guz19));
#endif
#if defined F2
    dpdx -= sxi*(atom2.quadrupoleXX*gux11 + atom2.quadrupoleYY*gux14 + atom2.quadrupoleZZ*gux16 + 2*(atom2.quadrupoleXY*gux12 + atom2.quadrupoleXZ*gux13 + atom2.quadrupoleYZ*gux15)) +
            syi*(atom2.quadrupoleXX*guy11 + atom2.quadrupoleYY*guy14 + atom2.quadrupoleZZ*guy16 + 2*(atom2.quadrupoleXY*guy12 + atom2.quadrupoleXZ*guy13 + atom2.quadrupoleYZ*guy15)) +
            szi*(atom2.quadrupoleXX*guz11 + atom2.quadrupoleYY*guz14 + atom2.quadrupoleZZ*guz16 + 2*(atom2.quadrupoleXY*guz12 + atom2.quadrupoleXZ*guz13 + atom2.quadrupoleYZ*guz15));

    dpdx += sxk*(atom1.quadrupoleXX*gux11 + atom1.quadrupoleYY*gux14 + atom1.quadrupoleZZ*gux16 + 2*(atom1.quadrupoleXY*gux12 + atom1.quadrupoleXZ*gux13 + atom1.quadrupoleYZ*gux15)) +
            syk*(atom1.quadrupoleXX*guy11 + atom1.quadrupoleYY*guy14 + atom1.quadrupoleZZ*guy16 + 2*(atom1.quadrupoleXY*guy12 + atom1.quadrupoleXZ*guy13 + atom1.quadrupoleYZ*guy15)) +
            szk*(atom1.quadrupoleXX*guz11 + atom1.quadrupoleYY*guz14 + atom1.quadrupoleZZ*guz16 + 2*(atom1.quadrupoleXY*guz12 + atom1.quadrupoleXZ*guz13 + atom1.quadrupoleYZ*guz15));

    dpdy -= sxi*(atom2.quadrupoleXX*gux12 + atom2.quadrupoleYY*gux17 + atom2.quadrupoleZZ*gux19 + 2*(atom2.quadrupoleXY*gux14 + atom2.quadrupoleXZ*gux15 + atom2.quadrupoleYZ*gux18)) +
            syi*(atom2.quadrupoleXX*guy12 + atom2.quadrupoleYY*guy17 + atom2.quadrupoleZZ*guy19 + 2*(atom2.quadrupoleXY*guy14 + atom2.quadrupoleXZ*guy15 + atom2.quadrupoleYZ*guy18)) +
            szi*(atom2.quadrupoleXX*guz12 + atom2.quadrupoleYY*guz17 + atom2.quadrupoleZZ*guz19 + 2*(atom2.quadrupoleXY*guz14 + atom2.quadrupoleXZ*guz15 + atom2.quadrupoleYZ*guz18));

    dpdy += sxk*(atom1.quadrupoleXX*gux12 + atom1.quadrupoleYY*gux17 + atom1.quadrupoleZZ*gux19 + 2*(atom1.quadrupoleXY*gux14 + atom1.quadrupoleXZ*gux15 + atom1.quadrupoleYZ*gux18)) +
            syk*(atom1.quadrupoleXX*guy12 + atom1.quadrupoleYY*guy17 + atom1.quadrupoleZZ*guy19 + 2*(atom1.quadrupoleXY*guy14 + atom1.quadrupoleXZ*guy15 + atom1.quadrupoleYZ*guy18)) +
            szk*(atom1.quadrupoleXX*guz12 + atom1.quadrupoleYY*guz17 + atom1.quadrupoleZZ*guz19 + 2*(atom1.quadrupoleXY*guz14 + atom1.quadrupoleXZ*guz15 + atom1.quadrupoleYZ*guz18));

    dpdz -= sxi*(atom2.quadrupoleXX*gux13 + atom2.quadrupoleYY*gux18 + atom2.quadrupoleZZ*gux20 + 2*(atom2.quadrupoleXY*gux15 + atom2.quadrupoleXZ*gux16 + atom2.quadrupoleYZ*gux19)) +
            syi*(atom2.quadrupoleXX*guy13 + atom2.quadrupoleYY*guy18 + atom2.quadrupoleZZ*guy20 + 2*(atom2.quadrupoleXY*guy15 + atom2.quadrupoleXZ*guy16 + atom2.quadrupoleYZ*guy19)) +
            szi*(atom2.quadrupoleXX*guz13 + atom2.quadrupoleYY*guz18 + atom2.quadrupoleZZ*guz20 + 2*(atom2.quadrupoleXY*guz15 + atom2.quadrupoleXZ*guz16 + atom2.quadrupoleYZ*guz19));

    dpdz += sxk*(atom1.quadrupoleXX*gux13 + atom1.quadrupoleYY*gux18 + atom1.quadrupoleZZ*gux20 + 2*(atom1.quadrupoleXY*gux15 + atom1.quadrupoleXZ*gux16 + atom1.quadrupoleYZ*gux19)) +
            syk*(atom1.quadrupoleXX*guy13 + atom1.quadrupoleYY*guy18 + atom1.quadrupoleZZ*guy20 + 2*(atom1.quadrupoleXY*guy15 + atom1.quadrupoleXZ*guy16 + atom1.quadrupoleYZ*guy19)) +
            szk*(atom1.quadrupoleXX*guz13 + atom1.quadrupoleYY*guz18 + atom1.quadrupoleZZ*guz20 + 2*(atom1.quadrupoleXY*guz15 + atom1.quadrupoleXZ*guz16 + atom1.quadrupoleYZ*guz19));
#endif

#endif

#if defined F1
    real gqxx11 = xr*(12*a21 + xr2*(9*a22 + xr2*a23));
    real gqxx12 = yr*(2*a21 + xr2*(5*a22 + xr2*a23));
    real gqxx13 = zr*(2*a21 + xr2*(5*a22 + xr2*a23));
    real gqxx14 = xr*(2*a21 + yr2*2*a22 + xr2*(a22 + yr2*a23));
    real gqxx15 = xr*yr*zr*(2*a22 + xr2*a23);
    real gqxx16 = xr*(2*a21 + zr2*2*a22 + xr2*(a22 + zr2*a23));
    real gqxx17 = yr*xr2*(3*a22 + yr2*a23);
    real gqxx18 = zr*xr2*(a22 + yr2*a23);
    real gqxx19 = yr*xr2*(a22 + zr2*a23);
    real gqxx20 = zr*xr2*(3*a22 + zr2*a23);

    real gqxy11 = yr*(3*a21 + xr2*(6*a22 + xr2*a23));
    real gqxy12 = xr*(3*(a21 + yr2*a22) + xr2*(a22 + yr2*a23));

    real gqxy13 = xr*yr*zr*(3*a22 + xr2*a23);

    real gqxy14 = yr*(3*(a21 + xr2*a22) + yr2*(a22 + xr2*a23));
    real gqxy15 = zr*(a21 + (yr2 + xr2)*a22 + yr2*xr2*a23);
    real gqxy16 = yr*(a21 + (xr2 + zr2)*a22 + xr2*zr2*a23);
    real gqxy17 = xr*(3*(a21 + yr2*a22) + yr2*(3*a22 + yr2*a23));
    real gqxy18 = xr*yr*zr*(3*a22 + yr2*a23);
    real gqxy19 = xr*(a21 + (yr2 + zr2)*a22 + yr2*zr2*a23);
    real gqxy20 = xr*yr*zr*(3*a22 + zr2*a23);
    real gqxz11 = zr*(3*a21 + xr2*(6*a22 + xr2*a23));

    real gqxz12 = xr*yr*zr*(3*a22 + xr2*a23);

    real gqxz13 = xr*(3*(a21 + zr2*a22) + xr2*(a22 + zr2*a23));
    real gqxz14 = zr*(a21 + (xr2 + yr2)*a22 + xr2*yr2*a23);
    real gqxz15 = yr*(a21 + (xr2 + zr2)*a22 + zr2*xr2*a23);
    real gqxz16 = zr*(3*(a21 + xr2*a22) + zr2*(a22 + xr2*a23));
    real gqxz17 = xr*yr*zr*(3*a22 + yr2*a23);
    real gqxz18 = xr*(a21 + (zr2 + yr2)*a22 + zr2*yr2*a23);
    real gqxz19 = xr*yr*zr*(3*a22 + zr2*a23);
    real gqxz20 = xr*(3*a21 + zr2*(6*a22 + zr2*a23));
    real gqyy11 = xr*yr2*(3*a22 + xr2*a23);
    real gqyy12 = yr*(2*a21 + xr2*2*a22 + yr2*(a22 + xr2*a23));
    real gqyy13 = zr*yr2*(a22 + xr2*a23);
    real gqyy14 = xr*(2*a21 + yr2*(5*a22 + yr2*a23));
    real gqyy15 = xr*yr*zr*(2*a22 + yr2*a23);
    real gqyy16 = xr*yr2*(a22 + zr2*a23);
    real gqyy17 = yr*(12*a21 + yr2*(9*a22 + yr2*a23));
    real gqyy18 = zr*(2*a21 + yr2*(5*a22 + yr2*a23));
    real gqyy19 = yr*(2*a21 + zr2*2*a22 + yr2*(a22 + zr2*a23));
    real gqyy20 = zr*yr2*(3*a22 + zr2*a23);
    real gqyz11 = xr*yr*zr*(3*a22 + xr2*a23);
    real gqyz12 = zr*(a21 + (xr2 + yr2)*a22 + xr2*yr2*a23);
    real gqyz13 = yr*(a21 + (xr2 + zr2)*a22 + xr2*zr2*a23);
    real gqyz14 = xr*yr*zr*(3*a22 + yr2*a23);
    real gqyz15 = xr*(a21 + (yr2 + zr2)*a22 + yr2*zr2*a23);
    real gqyz16 = xr*yr*zr*(3*a22 + zr2*a23);
    real gqyz17 = zr*(3*a21 + yr2*(6*a22 + yr2*a23));
    real gqyz18 = yr*(3*(a21 + zr2*a22) + yr2*(a22 + zr2*a23));
    real gqyz19 = zr*(3*(a21 + yr2*a22) + zr2*(a22 + yr2*a23));
    real gqyz20 = yr*(3*a21 + zr2*(6*a22 + zr2*a23));
    real gqzz11 = xr*zr2*(3*a22 + xr2*a23);
    real gqzz12 = yr*(zr2*a22 + xr2*(zr2*a23));
    real gqzz13 = zr*(2*a21 + xr2*2*a22 + zr2*(a22 + xr2*a23));
    real gqzz14 = xr*zr2*(a22 + yr2*a23);
    real gqzz15 = xr*yr*zr*(2*a22 + zr2*a23);
    real gqzz16 = xr*(2*a21 + zr2*(5*a22 + zr2*a23));
    real gqzz17 = yr*zr2*(3*a22 + yr2*a23);
    real gqzz18 = zr*(2*a21 + yr2*2*a22 + zr2*(a22 + yr2*a23));
    real gqzz19 = yr*(2*a21 + zr2*(5*a22 + zr2*a23));
    real gqzz20 = zr*(12*a21 + zr2*(9*a22 + zr2*a23));

    dedx += atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx11 + atom2.quadrupoleYY*gqxx14 + atom2.quadrupoleZZ*gqxx16 + 2*(atom2.quadrupoleXY*gqxx12 + atom2.quadrupoleXZ*gqxx13 + atom2.quadrupoleYZ*gqxx15)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqyy11 + atom2.quadrupoleYY*gqyy14 + atom2.quadrupoleZZ*gqyy16 + 2*(atom2.quadrupoleXY*gqyy12 + atom2.quadrupoleXZ*gqyy13 + atom2.quadrupoleYZ*gqyy15)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqzz11 + atom2.quadrupoleYY*gqzz14 + atom2.quadrupoleZZ*gqzz16 + 2*(atom2.quadrupoleXY*gqzz12 + atom2.quadrupoleXZ*gqzz13 + atom2.quadrupoleYZ*gqzz15)) +
          2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxy11 + atom2.quadrupoleYY*gqxy14 + atom2.quadrupoleZZ*gqxy16 + 2*(atom2.quadrupoleXY*gqxy12 + atom2.quadrupoleXZ*gqxy13 + atom2.quadrupoleYZ*gqxy15)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxz11 + atom2.quadrupoleYY*gqxz14 + atom2.quadrupoleZZ*gqxz16 + 2*(atom2.quadrupoleXY*gqxz12 + atom2.quadrupoleXZ*gqxz13 + atom2.quadrupoleYZ*gqxz15)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqyz11 + atom2.quadrupoleYY*gqyz14 + atom2.quadrupoleZZ*gqyz16 + 2*(atom2.quadrupoleXY*gqyz12 + atom2.quadrupoleXZ*gqyz13 + atom2.quadrupoleYZ*gqyz15))) +

            atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx11 + atom2.quadrupoleYY*gqyy11 + atom2.quadrupoleZZ*gqzz11 + 2*(atom2.quadrupoleXY*gqxy11 + atom2.quadrupoleXZ*gqxz11 + atom2.quadrupoleYZ*gqyz11)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqxx14 + atom2.quadrupoleYY*gqyy14 + atom2.quadrupoleZZ*gqzz14 + 2*(atom2.quadrupoleXY*gqxy14 + atom2.quadrupoleXZ*gqxz14 + atom2.quadrupoleYZ*gqyz14)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqxx16 + atom2.quadrupoleYY*gqyy16 + atom2.quadrupoleZZ*gqzz16 + 2*(atom2.quadrupoleXY*gqxy16 + atom2.quadrupoleXZ*gqxz16 + atom2.quadrupoleYZ*gqyz16)) +

          2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxx12 + atom2.quadrupoleYY*gqyy12 + atom2.quadrupoleZZ*gqzz12 + 2*(atom2.quadrupoleXY*gqxy12 + atom2.quadrupoleXZ*gqxz12 + atom2.quadrupoleYZ*gqyz12)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxx13 + atom2.quadrupoleYY*gqyy13 + atom2.quadrupoleZZ*gqzz13 + 2*(atom2.quadrupoleXY*gqxy13 + atom2.quadrupoleXZ*gqxz13 + atom2.quadrupoleYZ*gqyz13)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqxx15 + atom2.quadrupoleYY*gqyy15 + atom2.quadrupoleZZ*gqzz15 + 2*(atom2.quadrupoleXY*gqxy15 + atom2.quadrupoleXZ*gqxz15 + atom2.quadrupoleYZ*gqyz15)));

    dedy += atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx12 + atom2.quadrupoleYY*gqxx17 + atom2.quadrupoleZZ*gqxx19 + 2*(atom2.quadrupoleXY*gqxx14 + atom2.quadrupoleXZ*gqxx15 + atom2.quadrupoleYZ*gqxx18)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqyy12 + atom2.quadrupoleYY*gqyy17 + atom2.quadrupoleZZ*gqyy19 + 2*(atom2.quadrupoleXY*gqyy14 + atom2.quadrupoleXZ*gqyy15 + atom2.quadrupoleYZ*gqyy18)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqzz12 + atom2.quadrupoleYY*gqzz17 + atom2.quadrupoleZZ*gqzz19 + 2*(atom2.quadrupoleXY*gqzz14 + atom2.quadrupoleXZ*gqzz15 + atom2.quadrupoleYZ*gqzz18)) +

          2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxy12 + atom2.quadrupoleYY*gqxy17 + atom2.quadrupoleZZ*gqxy19 + 2*(atom2.quadrupoleXY*gqxy14 + atom2.quadrupoleXZ*gqxy15 + atom2.quadrupoleYZ*gqxy18)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxz12 + atom2.quadrupoleYY*gqxz17 + atom2.quadrupoleZZ*gqxz19 + 2*(atom2.quadrupoleXY*gqxz14 + atom2.quadrupoleXZ*gqxz15 + atom2.quadrupoleYZ*gqxz18)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqyz12 + atom2.quadrupoleYY*gqyz17 + atom2.quadrupoleZZ*gqyz19 + 2*(atom2.quadrupoleXY*gqyz14 + atom2.quadrupoleXZ*gqyz15 + atom2.quadrupoleYZ*gqyz18))) +

            atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx12 + atom2.quadrupoleYY*gqyy12 + atom2.quadrupoleZZ*gqzz12 + 2*(atom2.quadrupoleXY*gqxy12 + atom2.quadrupoleXZ*gqxz12 + atom2.quadrupoleYZ*gqyz12)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqxx17 + atom2.quadrupoleYY*gqyy17 + atom2.quadrupoleZZ*gqzz17 + 2*(atom2.quadrupoleXY*gqxy17 + atom2.quadrupoleXZ*gqxz17 + atom2.quadrupoleYZ*gqyz17)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqxx19 + atom2.quadrupoleYY*gqyy19 + atom2.quadrupoleZZ*gqzz19 + 2*(atom2.quadrupoleXY*gqxy19 + atom2.quadrupoleXZ*gqxz19 + atom2.quadrupoleYZ*gqyz19)) +

          2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxx14 + atom2.quadrupoleYY*gqyy14 + atom2.quadrupoleZZ*gqzz14 + 2*(atom2.quadrupoleXY*gqxy14 + atom2.quadrupoleXZ*gqxz14 + atom2.quadrupoleYZ*gqyz14)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxx15 + atom2.quadrupoleYY*gqyy15 + atom2.quadrupoleZZ*gqzz15 + 2*(atom2.quadrupoleXY*gqxy15 + atom2.quadrupoleXZ*gqxz15 + atom2.quadrupoleYZ*gqyz15)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqxx18 + atom2.quadrupoleYY*gqyy18 + atom2.quadrupoleZZ*gqzz18 + 2*(atom2.quadrupoleXY*gqxy18 + atom2.quadrupoleXZ*gqxz18 + atom2.quadrupoleYZ*gqyz18)));

    dedz += atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx13 + atom2.quadrupoleYY*gqxx18 + atom2.quadrupoleZZ*gqxx20 + 2*(atom2.quadrupoleXY*gqxx15 + atom2.quadrupoleXZ*gqxx16 + atom2.quadrupoleYZ*gqxx19)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqyy13 + atom2.quadrupoleYY*gqyy18 + atom2.quadrupoleZZ*gqyy20 + 2*(atom2.quadrupoleXY*gqyy15 + atom2.quadrupoleXZ*gqyy16 + atom2.quadrupoleYZ*gqyy19)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqzz13 + atom2.quadrupoleYY*gqzz18 + atom2.quadrupoleZZ*gqzz20 + 2*(atom2.quadrupoleXY*gqzz15 + atom2.quadrupoleXZ*gqzz16 + atom2.quadrupoleYZ*gqzz19)) +

           2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxy13 + atom2.quadrupoleYY*gqxy18 + atom2.quadrupoleZZ*gqxy20 + 2*(atom2.quadrupoleXY*gqxy15 + atom2.quadrupoleXZ*gqxy16 + atom2.quadrupoleYZ*gqxy19)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxz13 + atom2.quadrupoleYY*gqxz18 + atom2.quadrupoleZZ*gqxz20 + 2*(atom2.quadrupoleXY*gqxz15 + atom2.quadrupoleXZ*gqxz16 + atom2.quadrupoleYZ*gqxz19)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqyz13 + atom2.quadrupoleYY*gqyz18 + atom2.quadrupoleZZ*gqyz20 + 2*(atom2.quadrupoleXY*gqyz15 + atom2.quadrupoleXZ*gqyz16 + atom2.quadrupoleYZ*gqyz19))) +

            atom1.quadrupoleXX*(atom2.quadrupoleXX*gqxx13 + atom2.quadrupoleYY*gqyy13 + atom2.quadrupoleZZ*gqzz13 + 2*(atom2.quadrupoleXY*gqxy13 + atom2.quadrupoleXZ*gqxz13 + atom2.quadrupoleYZ*gqyz13)) +
            atom1.quadrupoleYY*(atom2.quadrupoleXX*gqxx18 + atom2.quadrupoleYY*gqyy18 + atom2.quadrupoleZZ*gqzz18 + 2*(atom2.quadrupoleXY*gqxy18 + atom2.quadrupoleXZ*gqxz18 + atom2.quadrupoleYZ*gqyz18)) +
            atom1.quadrupoleZZ*(atom2.quadrupoleXX*gqxx20 + atom2.quadrupoleYY*gqyy20 + atom2.quadrupoleZZ*gqzz20 + 2*(atom2.quadrupoleXY*gqxy20 + atom2.quadrupoleXZ*gqxz20 + atom2.quadrupoleYZ*gqyz20)) +

           2*(
            atom1.quadrupoleXY*(atom2.quadrupoleXX*gqxx15 + atom2.quadrupoleYY*gqyy15 + atom2.quadrupoleZZ*gqzz15 + 2*(atom2.quadrupoleXY*gqxy15 + atom2.quadrupoleXZ*gqxz15 + atom2.quadrupoleYZ*gqyz15)) +
            atom1.quadrupoleXZ*(atom2.quadrupoleXX*gqxx16 + atom2.quadrupoleYY*gqyy16 + atom2.quadrupoleZZ*gqzz16 + 2*(atom2.quadrupoleXY*gqxy16 + atom2.quadrupoleXZ*gqxz16 + atom2.quadrupoleYZ*gqyz16)) +
            atom1.quadrupoleYZ*(atom2.quadrupoleXX*gqxx19 + atom2.quadrupoleYY*gqyy19 + atom2.quadrupoleZZ*gqzz19 + 2*(atom2.quadrupoleXY*gqxy19 + atom2.quadrupoleXZ*gqxz19 + atom2.quadrupoleYZ*gqyz19)));
#endif


#if defined T1
    if (xr != 0 || yr != 0 || zr != 0) {

        real gux1 = xr*a10;
        real guy1 = yr*a10;
        real guz1 = zr*a10;

        real gc2 = xr*a01;
        real gc3 = yr*a01;
        real gc4 = zr*a01;
        real fid1 = atom2.dipole.x*gux2 + atom2.dipole.y*gux3 + atom2.dipole.z*gux4 + 0.5f*(atom2.q*gux1 + atom2.quadrupoleXX*gux5 + atom2.quadrupoleYY*gux8 + atom2.quadrupoleZZ*gux10 +
                           2*(atom2.quadrupoleXY*gux6 + atom2.quadrupoleXZ*gux7 + atom2.quadrupoleYZ*gux9) +
                          atom2.q*gc2 + atom2.quadrupoleXX*gqxx2 + atom2.quadrupoleYY*gqyy2 + atom2.quadrupoleZZ*gqzz2 +
                           2*(atom2.quadrupoleXY*gqxy2 + atom2.quadrupoleXZ*gqxz2 + atom2.quadrupoleYZ*gqyz2));

        real fid2 = atom2.dipole.x*guy2 + atom2.dipole.y*guy3 + atom2.dipole.z*guy4 + 0.5f*(atom2.q*guy1 + atom2.quadrupoleXX*guy5 + atom2.quadrupoleYY*guy8 + atom2.quadrupoleZZ*guy10 +
                           2*(atom2.quadrupoleXY*guy6 + atom2.quadrupoleXZ*guy7 + atom2.quadrupoleYZ*guy9) +
                          atom2.q*gc3 + atom2.quadrupoleXX*gqxx3 + atom2.quadrupoleYY*gqyy3 + atom2.quadrupoleZZ*gqzz3 + 
                           2*(atom2.quadrupoleXY*gqxy3 + atom2.quadrupoleXZ*gqxz3 + atom2.quadrupoleYZ*gqyz3));

        real fid3 = atom2.dipole.x*guz2 + atom2.dipole.y*guz3 + atom2.dipole.z*guz4 + 0.5f*(atom2.q*guz1 + atom2.quadrupoleXX*guz5 + atom2.quadrupoleYY*guz8 + atom2.quadrupoleZZ*guz10 +
                           2*(atom2.quadrupoleXY*guz6 + atom2.quadrupoleXZ*guz7 + atom2.quadrupoleYZ*guz9) +
                           atom2.q*gc4 + atom2.quadrupoleXX*gqxx4 + atom2.quadrupoleYY*gqyy4 + atom2.quadrupoleZZ*gqzz4 +
                           2*(atom2.quadrupoleXY*gqxy4 + atom2.quadrupoleXZ*gqxz4 + atom2.quadrupoleYZ*gqyz4));

        real trq1 = atom1.dipole.y*fid3 - atom1.dipole.z*fid2;
        real trq2 = atom1.dipole.z*fid1 - atom1.dipole.x*fid3;
        real trq3 = atom1.dipole.x*fid2 - atom1.dipole.y*fid1;

        // torque on quadrupoles due to permanent reaction field gradient

        real fidg11 =
                (atom2.q*xr2*a20 + atom2.dipole.x*gqxx2 + atom2.dipole.y*gqxx3 + atom2.dipole.z*gqxx4
                       + atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqxx8 + atom2.quadrupoleZZ*gqxx10
                       + 2*(atom2.quadrupoleXY*gqxx6 + atom2.quadrupoleXZ*gqxx7 + atom2.quadrupoleYZ*gqxx9)
                       + atom2.q*gc5 + atom2.dipole.x*gux5 + atom2.dipole.y*guy5 + atom2.dipole.z*guz5
                       + atom2.quadrupoleXX*gqxx5 + atom2.quadrupoleYY*gqyy5 + atom2.quadrupoleZZ*gqzz5
                       + 2*(atom2.quadrupoleXY*gqxy5 + atom2.quadrupoleXZ*gqxz5 + atom2.quadrupoleYZ*gqyz5));

        real fidg12 =
                (atom2.q*xr*yr*a20 + atom2.dipole.x*gqxy2 + atom2.dipole.y*gqxy3 + atom2.dipole.z*gqxy4
                       + atom2.quadrupoleXX*gqxy5 + atom2.quadrupoleYY*gqxy8 + atom2.quadrupoleZZ*gqxy10
                       + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxy7 + atom2.quadrupoleYZ*gqxy9)
                       + atom2.q*gc6 + atom2.dipole.x*gux6 + atom2.dipole.y*guy6 + atom2.dipole.z*guz6
                       + atom2.quadrupoleXX*gqxx6 + atom2.quadrupoleYY*gqyy6 + atom2.quadrupoleZZ*gqzz6
                       + 2*(atom2.quadrupoleXY*gqxy6 + atom2.quadrupoleXZ*gqxz6 + atom2.quadrupoleYZ*gqyz6));

        real fidg13 =
                (atom2.q*xr*zr*a20 + atom2.dipole.x*gqxz2 + atom2.dipole.y*gqxz3 + atom2.dipole.z*gqxz4
                       + atom2.quadrupoleXX*gqxz5 + atom2.quadrupoleYY*gqxz8 + atom2.quadrupoleZZ*gqxz10
                       + 2*(atom2.quadrupoleXY*gqxz6 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqxz9)
                       + atom2.q*gc7 + atom2.dipole.x*gux7 + atom2.dipole.y*guy7 + atom2.dipole.z*guz7
                       + atom2.quadrupoleXX*gqxx7 + atom2.quadrupoleYY*gqyy7 + atom2.quadrupoleZZ*gqzz7
                       + 2*(atom2.quadrupoleXY*gqxy7 + atom2.quadrupoleXZ*gqxz7 + atom2.quadrupoleYZ*gqyz7));

        real fidg22 =
                (atom2.q*yr2*a20 + atom2.dipole.x*gqyy2 + atom2.dipole.y*gqyy3 + atom2.dipole.z*gqyy4
                       + atom2.quadrupoleXX*gqyy5 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqyy10
                       + 2*(atom2.quadrupoleXY*gqyy6 + atom2.quadrupoleXZ*gqyy7 + atom2.quadrupoleYZ*gqyy9)
                       + atom2.q*gc8 + atom2.dipole.x*gux8 + atom2.dipole.y*guy8 + atom2.dipole.z*guz8
                       + atom2.quadrupoleXX*gqxx8 + atom2.quadrupoleYY*gqyy8 + atom2.quadrupoleZZ*gqzz8
                       + 2*(atom2.quadrupoleXY*gqxy8 + atom2.quadrupoleXZ*gqxz8 + atom2.quadrupoleYZ*gqyz8));

        real fidg23 =
                (atom2.q*yr*zr*a20 + atom2.dipole.x*gqyz2 + atom2.dipole.y*gqyz3 + atom2.dipole.z*gqyz4
                       + atom2.quadrupoleXX*gqyz5 + atom2.quadrupoleYY*gqyz8 + atom2.quadrupoleZZ*gqyz10
                       + 2*(atom2.quadrupoleXY*gqyz6 + atom2.quadrupoleXZ*gqyz7 + atom2.quadrupoleYZ*gqyz9)
                       + atom2.q*gc9 + atom2.dipole.x*gux9 + atom2.dipole.y*guy9 + atom2.dipole.z*guz9
                       + atom2.quadrupoleXX*gqxx9 + atom2.quadrupoleYY*gqyy9 + atom2.quadrupoleZZ*gqzz9
                       + 2*(atom2.quadrupoleXY*gqxy9 + atom2.quadrupoleXZ*gqxz9 + atom2.quadrupoleYZ*gqyz9));

        real fidg33 =
                (atom2.q*zr2*a20 + atom2.dipole.x*gqzz2 + atom2.dipole.y*gqzz3 + atom2.dipole.z*gqzz4
                       + atom2.quadrupoleXX*gqzz5 + atom2.quadrupoleYY*gqzz8 + atom2.quadrupoleZZ*gqzz10
                       + 2*(atom2.quadrupoleXY*gqzz6 + atom2.quadrupoleXZ*gqzz7 + atom2.quadrupoleYZ*gqzz9)
                       + atom2.q*gc10 + atom2.dipole.x*gux10 + atom2.dipole.y*guy10 + atom2.dipole.z*guz10
                       + atom2.quadrupoleXX*gqxx10 + atom2.quadrupoleYY*gqyy10 + atom2.quadrupoleZZ*gqzz10
                    + 2*(atom2.quadrupoleXY*gqxy10 + atom2.quadrupoleXZ*gqxz10 + atom2.quadrupoleYZ*gqyz10));

        trq1 -= (atom1.quadrupoleXY*fidg13 + atom1.quadrupoleYY*fidg23 + atom1.quadrupoleYZ*fidg33 -atom1.quadrupoleXZ*fidg12-atom1.quadrupoleYZ*fidg22-atom1.quadrupoleZZ*fidg23);
        trq2 -= (atom1.quadrupoleXZ*fidg11 + atom1.quadrupoleYZ*fidg12 + atom1.quadrupoleZZ*fidg13 -atom1.quadrupoleXX*fidg13-atom1.quadrupoleXY*fidg23-atom1.quadrupoleXZ*fidg33);
        trq3 -= (atom1.quadrupoleXX*fidg12 + atom1.quadrupoleXY*fidg22 + atom1.quadrupoleXZ*fidg23 -atom1.quadrupoleXY*fidg11-atom1.quadrupoleYY*fidg12-atom1.quadrupoleYZ*fidg13);

        torque.x = trq1;
        torque.y = trq2;
        torque.z = trq3;

    } else {
        torque.x = 0;
        torque.y = 0;
        torque.z = 0;
    }
#endif

#if defined B2 
    dsumdrB2 *= 0.5f;
    atom1.bornForce += 0.5f*atom2.bornRadius*dsumdrB2;
    atom2.bornForce += 0.5f*atom1.bornRadius*dsumdrB2;
#endif

#if defined T2
    // torque due to induced reaction field gradient on quadrupoles;

    real fidg11 = sxk*gqxx2 + syk*gqxx3 + szk*gqxx4 + sxk*gux5 + syk*guy5 + szk*guz5;
    real fidg12 = sxk*gqxy2 + syk*gqxy3 + szk*gqxy4 + sxk*gux6 + syk*guy6 + szk*guz6;
    real fidg13 = sxk*gqxz2 + syk*gqxz3 + szk*gqxz4 + sxk*gux7 + syk*guy7 + szk*guz7;
    real fidg22 = sxk*gqyy2 + syk*gqyy3 + szk*gqyy4 + sxk*gux8 + syk*guy8 + szk*guz8;
    real fidg23 = sxk*gqyz2 + syk*gqyz3 + szk*gqyz4 + sxk*gux9 + syk*guy9 + szk*guz9;
    real fidg33 = sxk*gqzz2 + syk*gqzz3 + szk*gqzz4 + sxk*gux10 + syk*guy10 + szk*guz10;

    trqi1 -= atom1.quadrupoleXY*fidg13 + atom1.quadrupoleYY*fidg23 + atom1.quadrupoleYZ*fidg33
                                -atom1.quadrupoleXZ*fidg12 - atom1.quadrupoleYZ*fidg22 - atom1.quadrupoleZZ*fidg23;

    trqi2 -= atom1.quadrupoleXZ*fidg11 + atom1.quadrupoleYZ*fidg12 + atom1.quadrupoleZZ*fidg13
                                -atom1.quadrupoleXX*fidg13 - atom1.quadrupoleXY*fidg23 - atom1.quadrupoleXZ*fidg33;

    trqi3 -= atom1.quadrupoleXX*fidg12 + atom1.quadrupoleXY*fidg22 + atom1.quadrupoleXZ*fidg23
                                -atom1.quadrupoleXY*fidg11 - atom1.quadrupoleYY*fidg12 - atom1.quadrupoleYZ*fidg13;

    torque.x += 0.5f*trqi1;
    torque.y += 0.5f*trqi2;
    torque.z += 0.5f*trqi3;
#endif

#if defined F1
    
    outputEnergy += energy;

    if ((xr != 0 || yr != 0 || zr != 0)) {
        force.x = dedx;
        force.y = dedy;
        force.z = dedz;
    } else {
        force.x = force.y = force.z = 0;
    }

#endif

#if defined F2
    outputEnergy += 0.5f*energy;

    dpdx *= 0.5f;
    dpdy *= 0.5f;
    dpdz *= 0.5f;

    if ((xr != 0 || yr != 0 || zr != 0)) {
        force.x += dpdx;
        force.y += dpdy;
        force.z += dpdz;
    }
#endif
}
