KERNEL void applyPositionDeltas(int numAtoms, GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT posDelta
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        ) {
    for (unsigned int index = GLOBAL_ID; index < numAtoms; index += GLOBAL_SIZE) {
#ifdef USE_MIXED_PRECISION
        real4 pos1 = posq[index];
        real4 pos2 = posqCorrection[index];
        mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
        mixed4 pos = posq[index];
#endif
        pos.x += posDelta[index].x;
        pos.y += posDelta[index].y;
        pos.z += posDelta[index].z;
#ifdef USE_MIXED_PRECISION
        posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
        posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
        posq[index] = pos;
#endif
    }
}

/**
 * Load the position of a particle.
 */
inline DEVICE mixed4 loadPos(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
inline DEVICE void storePos(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

/**
 * Enforce position constraints on a single SHAKE cluster.
 */
inline DEVICE void applyShakeToPositionsCluster(int4 atoms, float4 params, mixed tol,
        mixed4 pos, mixed4 pos1, mixed4 pos2, mixed4 pos3,
        PRIVATE mixed4* xp0Ptr, PRIVATE mixed4* xp1Ptr, PRIVATE mixed4* xp2Ptr, PRIVATE mixed4* xp3Ptr) {
    mixed4 xpi = *xp0Ptr;
    mixed4 xpj1 = *xp1Ptr;
    mixed4 xpj2 = *xp2Ptr;
    mixed4 xpj3 = *xp3Ptr;
    float invMassCentral = params.x;
    float avgMass = params.y;
    float d2 = params.z;
    float invMassPeripheral = params.w;

    // Precompute quantities.

    mixed3 rij1 = make_mixed3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
    mixed3 rij2 = make_mixed3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
    mixed3 rij3 = make_mixed3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
    mixed rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
    mixed rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
    mixed rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
    mixed ld1 = d2-rij1sq;
    mixed ld2 = d2-rij2sq;
    mixed ld3 = d2-rij3sq;

    // Iterate until convergence.

    bool converged = false;
    int iteration = 0;
    while (iteration < 15 && !converged) {
        converged = true;
        mixed3 rpij = make_mixed3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
        mixed rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
        mixed rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
        mixed diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
        if (diff >= 1.0f) {
            mixed acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
            mixed3 dr = rij1*acor;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj1.x -= dr.x*invMassPeripheral;
            xpj1.y -= dr.y*invMassPeripheral;
            xpj1.z -= dr.z*invMassPeripheral;
            converged = false;
        }
        if (atoms.z != -1) {
            rpij = make_mixed3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
            rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            rrpr = rij2.x*rpij.x + rij2.y*rpij.y + rij2.z*rpij.z;
            diff = fabs(ld2-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                mixed acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                mixed3 dr = rij2*acor;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj2.x -= dr.x*invMassPeripheral;
                xpj2.y -= dr.y*invMassPeripheral;
                xpj2.z -= dr.z*invMassPeripheral;
                converged = false;
            }
        }
        if (atoms.w != -1) {
            rpij = make_mixed3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
            rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            rrpr = rij3.x*rpij.x + rij3.y*rpij.y + rij3.z*rpij.z;
            diff = fabs(ld3 - 2.0f*rrpr - rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                mixed acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                mixed3 dr = rij3*acor;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj3.x -= dr.x*invMassPeripheral;
                xpj3.y -= dr.y*invMassPeripheral;
                xpj3.z -= dr.z*invMassPeripheral;
                converged = false;
            }
        }
        iteration++;
    }
    *xp0Ptr = xpi;
    *xp1Ptr = xpj1;
    *xp2Ptr = xpj2;
    *xp3Ptr = xpj3;
}

/**
 * Enforce velocity constraints on a single SHAKE cluster.
 */
inline DEVICE void applyShakeToVelocitiesCluster(int4 atoms, float4 params, mixed tol,
        mixed4 pos, mixed4 pos1, mixed4 pos2, mixed4 pos3,
        PRIVATE mixed4* v0Ptr, PRIVATE mixed4* v1Ptr, PRIVATE mixed4* v2Ptr, PRIVATE mixed4* v3Ptr) {
    mixed4 xpi = *v0Ptr;
    mixed4 xpj1 = *v1Ptr;
    mixed4 xpj2 = *v2Ptr;
    mixed4 xpj3 = *v3Ptr;
    float invMassCentral = params.x;
    float avgMass = params.y;
    float invMassPeripheral = params.w;

    // Precompute quantities.

    mixed3 rij1 = make_mixed3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
    mixed3 rij2 = make_mixed3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
    mixed3 rij3 = make_mixed3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
    mixed rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
    mixed rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
    mixed rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;

    // Iterate until convergence.

    bool converged = false;
    int iteration = 0;
    while (iteration < 15 && !converged) {
        converged = true;
        mixed3 rpij = make_mixed3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
        mixed rrpr = rpij.x*rij1.x + rpij.y*rij1.y + rpij.z*rij1.z;
        mixed delta = -2.0f*avgMass*rrpr/rij1sq;
        mixed3 dr = rij1*delta;
        xpi.x += dr.x*invMassCentral;
        xpi.y += dr.y*invMassCentral;
        xpi.z += dr.z*invMassCentral;
        xpj1.x -= dr.x*invMassPeripheral;
        xpj1.y -= dr.y*invMassPeripheral;
        xpj1.z -= dr.z*invMassPeripheral;
        if (fabs(delta) > tol)
            converged = false;
        if (atoms.z != -1) {
            rpij = make_mixed3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
            rrpr = rpij.x*rij2.x + rpij.y*rij2.y + rpij.z*rij2.z;
            delta = -2.0f*avgMass*rrpr/rij2sq;
            dr = rij2*delta;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj2.x -= dr.x*invMassPeripheral;
            xpj2.y -= dr.y*invMassPeripheral;
            xpj2.z -= dr.z*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
        }
        if (atoms.w != -1) {
            rpij = make_mixed3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
            rrpr = rpij.x*rij3.x + rpij.y*rij3.y + rpij.z*rij3.z;
            delta = -2.0f*avgMass*rrpr/rij3sq;
            dr = rij3*delta;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj3.x -= dr.x*invMassPeripheral;
            xpj3.y -= dr.y*invMassPeripheral;
            xpj3.z -= dr.z*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
        }
        iteration++;
    }
    *v0Ptr = xpi;
    *v1Ptr = xpj1;
    *v2Ptr = xpj2;
    *v3Ptr = xpj3;
}

/**
 * Enforce position constraints on a single SETTLE cluster.
 */
inline DEVICE void applySettleToPositionsCluster(float2 params, mixed4 apos0, mixed4 apos1, mixed4 apos2,
        mixed m0, mixed m1, mixed m2, PRIVATE mixed4* xp0Ptr, PRIVATE mixed4* xp1Ptr, PRIVATE mixed4* xp2Ptr) {
    mixed4 xp0 = *xp0Ptr;
    mixed4 xp1 = *xp1Ptr;
    mixed4 xp2 = *xp2Ptr;

    // Apply the SETTLE algorithm.

    mixed xb0 = apos1.x-apos0.x;
    mixed yb0 = apos1.y-apos0.y;
    mixed zb0 = apos1.z-apos0.z;
    mixed xc0 = apos2.x-apos0.x;
    mixed yc0 = apos2.y-apos0.y;
    mixed zc0 = apos2.z-apos0.z;

    mixed invTotalMass = 1/(m0+m1+m2);
    mixed xcom = (xp0.x*m0 + (xb0+xp1.x)*m1 + (xc0+xp2.x)*m2) * invTotalMass;
    mixed ycom = (xp0.y*m0 + (yb0+xp1.y)*m1 + (yc0+xp2.y)*m2) * invTotalMass;
    mixed zcom = (xp0.z*m0 + (zb0+xp1.z)*m1 + (zc0+xp2.z)*m2) * invTotalMass;

    mixed xa1 = xp0.x - xcom;
    mixed ya1 = xp0.y - ycom;
    mixed za1 = xp0.z - zcom;
    mixed xb1 = xb0 + xp1.x - xcom;
    mixed yb1 = yb0 + xp1.y - ycom;
    mixed zb1 = zb0 + xp1.z - zcom;
    mixed xc1 = xc0 + xp2.x - xcom;
    mixed yc1 = yc0 + xp2.y - ycom;
    mixed zc1 = zc0 + xp2.z - zcom;

    mixed xaksZd = yb0*zc0 - zb0*yc0;
    mixed yaksZd = zb0*xc0 - xb0*zc0;
    mixed zaksZd = xb0*yc0 - yb0*xc0;
    mixed xaksXd = ya1*zaksZd - za1*yaksZd;
    mixed yaksXd = za1*xaksZd - xa1*zaksZd;
    mixed zaksXd = xa1*yaksZd - ya1*xaksZd;
    mixed xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
    mixed yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
    mixed zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

    mixed axlng = sqrt(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
    mixed aylng = sqrt(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
    mixed azlng = sqrt(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
    mixed trns11 = xaksXd / axlng;
    mixed trns21 = yaksXd / axlng;
    mixed trns31 = zaksXd / axlng;
    mixed trns12 = xaksYd / aylng;
    mixed trns22 = yaksYd / aylng;
    mixed trns32 = zaksYd / aylng;
    mixed trns13 = xaksZd / azlng;
    mixed trns23 = yaksZd / azlng;
    mixed trns33 = zaksZd / azlng;

    mixed xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
    mixed yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
    mixed xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
    mixed yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
    mixed za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
    mixed xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
    mixed yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
    mixed zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
    mixed xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
    mixed yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
    mixed zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

    //                                        --- Step2  A2' ---

    float rc = 0.5f*params.y;
    mixed rb = sqrt(params.x*params.x-rc*rc);
    mixed ra = rb*(m1+m2)*invTotalMass;
    rb -= ra;
    mixed sinphi = za1d/ra;
    mixed cosphi = sqrt(1-sinphi*sinphi);
    mixed sinpsi = (zb1d-zc1d) / (2*rc*cosphi);
    mixed cospsi = sqrt(1-sinpsi*sinpsi);

    mixed ya2d =   ra*cosphi;
    mixed xb2d = - rc*cospsi;
    mixed yb2d = - rb*cosphi - rc*sinpsi*sinphi;
    mixed yc2d = - rb*cosphi + rc*sinpsi*sinphi;
    mixed xb2d2 = xb2d*xb2d;
    mixed hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
    mixed deltx = 2.0f*xb2d + sqrt(4.0f*xb2d2 - hh2 + params.y*params.y);
    xb2d -= deltx*0.5f;

    //                                        --- Step3  al,be,ga ---

    mixed alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
    mixed beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
    mixed gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

    mixed al2be2 = alpha*alpha + beta*beta;
    mixed sintheta = (alpha*gamma - beta*sqrt(al2be2 - gamma*gamma)) / al2be2;

    //                                        --- Step4  A3' ---

    mixed costheta = sqrt(1-sintheta*sintheta);
    mixed xa3d = - ya2d*sintheta;
    mixed ya3d =   ya2d*costheta;
    mixed za3d = za1d;
    mixed xb3d =   xb2d*costheta - yb2d*sintheta;
    mixed yb3d =   xb2d*sintheta + yb2d*costheta;
    mixed zb3d = zb1d;
    mixed xc3d = - xb2d*costheta - yc2d*sintheta;
    mixed yc3d = - xb2d*sintheta + yc2d*costheta;
    mixed zc3d = zc1d;

    //                                        --- Step5  A3 ---

    mixed xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
    mixed ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
    mixed za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
    mixed xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
    mixed yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
    mixed zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
    mixed xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
    mixed yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
    mixed zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

    xp0.x = xcom + xa3;
    xp0.y = ycom + ya3;
    xp0.z = zcom + za3;
    xp1.x = xcom + xb3 - xb0;
    xp1.y = ycom + yb3 - yb0;
    xp1.z = zcom + zb3 - zb0;
    xp2.x = xcom + xc3 - xc0;
    xp2.y = ycom + yc3 - yc0;
    xp2.z = zcom + zc3 - zc0;

    *xp0Ptr = xp0;
    *xp1Ptr = xp1;
    *xp2Ptr = xp2;
}

/**
 * Enforce velocity constraints on a single SETTLE cluster.
 */
inline DEVICE void applySettleToVelocitiesCluster(mixed4 apos0, mixed4 apos1, mixed4 apos2,
        PRIVATE mixed4* v0Ptr, PRIVATE mixed4* v1Ptr, PRIVATE mixed4* v2Ptr) {
    mixed4 v0 = *v0Ptr;
    mixed4 v1 = *v1Ptr;
    mixed4 v2 = *v2Ptr;

    // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
    // and the angle cosines and sines.

    mixed mA = 1/v0.w;
    mixed mB = 1/v1.w;
    mixed mC = 1/v2.w;
    mixed3 eAB = make_mixed3(apos1.x-apos0.x, apos1.y-apos0.y, apos1.z-apos0.z);
    mixed3 eBC = make_mixed3(apos2.x-apos1.x, apos2.y-apos1.y, apos2.z-apos1.z);
    mixed3 eCA = make_mixed3(apos0.x-apos2.x, apos0.y-apos2.y, apos0.z-apos2.z);
    eAB *= RSQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z);
    eBC *= RSQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z);
    eCA *= RSQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z);
    mixed vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
    mixed vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
    mixed vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
    mixed cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
    mixed cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
    mixed cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
    mixed s2A = 1-cA*cA;
    mixed s2B = 1-cB*cB;
    mixed s2C = 1-cC*cC;

    // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
    // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
    // making that assumption).  We allow all three atoms to have different masses.

    mixed mABCinv = 1/(mA*mB*mC);
    mixed denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
    mixed tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
    mixed tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
    mixed tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
    v0.x += (tab*eAB.x - tca*eCA.x)*v0.w;
    v0.y += (tab*eAB.y - tca*eCA.y)*v0.w;
    v0.z += (tab*eAB.z - tca*eCA.z)*v0.w;
    v1.x += (tbc*eBC.x - tab*eAB.x)*v1.w;
    v1.y += (tbc*eBC.y - tab*eAB.y)*v1.w;
    v1.z += (tbc*eBC.z - tab*eAB.z)*v1.w;
    v2.x += (tca*eCA.x - tbc*eBC.x)*v2.w;
    v2.y += (tca*eCA.y - tbc*eBC.y)*v2.w;
    v2.z += (tca*eCA.z - tbc*eBC.z)*v2.w;
    *v0Ptr = v0;
    *v1Ptr = v1;
    *v2Ptr = v2;
}
