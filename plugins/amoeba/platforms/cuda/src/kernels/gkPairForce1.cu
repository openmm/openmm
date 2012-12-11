/**
 * This defines three different closely related functions, depending on which constant (F1, T1, or T3) is defined.
 * To work around limitations in Visual Studio, this kernel is split into two files.
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
