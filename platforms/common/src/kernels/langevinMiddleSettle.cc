// Experimental CUDA/mixed-only SETTLE position + exact LangevinMiddle Part3.
// SETTLE body copied verbatim from integrationUtilities.cc except signature,
// original-delta capture and the final stores. Keep them synchronized.
#ifndef USE_MIXED_PRECISION
#error LangevinMiddle SETTLE fusion requires mixed precision
#endif

inline DEVICE mixed4 loadPos(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, int index) {
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
}

inline DEVICE mixed settleMiddleStoredDouble(mixed value) {
    // Volatile inline assembly is opaque to front-end contraction across this
    // boundary without requiring a global scratch store. Not a bit-equivalence
    // claim: final generated instructions still require numerical validation.
    mixed rounded;
    asm volatile("mov.b64 %0, %1;" : "=d"(rounded) : "d"(value));
    return rounded;
}

inline DEVICE void finishSettleMiddleAtom(int index, mixed4 constrained,
#ifdef RELOAD_LANGEVIN_MIDDLE_ORIGINAL_DELTA
        GLOBAL const volatile mixed4* RESTRICT originalDelta,
#else
        mixed3 original,
#endif
        mixed invDt,
        GLOBAL real4* RESTRICT posq, GLOBAL mixed4* RESTRICT velm, GLOBAL real4* RESTRICT posqCorrection) {
    mixed4 delta = make_mixed4(settleMiddleStoredDouble(constrained.x),
            settleMiddleStoredDouble(constrained.y), settleMiddleStoredDouble(constrained.z), 0);
#ifdef RELOAD_LANGEVIN_MIDDLE_ORIGINAL_DELTA
    // SETTLE fusion never overwrites posDelta for its disjoint fast atoms.
    // Re-read exact stored xyz instead of retaining nine original doubles
    // across SETTLE. Volatile requests actual loads; it is not synchronization
    // or a guarantee about generated register allocation or load scheduling.
    mixed3 original;
    original.x = originalDelta[index].x;
    original.y = originalDelta[index].y;
    original.z = originalDelta[index].z;
#endif
    mixed4 velocity = velm[index];
    if (velocity.w != 0.0) {
        velocity.x += (delta.x-original.x)*invDt;
        velocity.y += (delta.y-original.y)*invDt;
        velocity.z += (delta.z-original.z)*invDt;
        velm[index] = velocity;
        // Keep the original Part3 position reconstruction and w handling.
        // Rereading positions limits live ranges; this is the 120 B/atom scope.
        real4 pos1 = posq[index];
        real4 pos2 = posqCorrection[index];
        mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
        pos.x += delta.x;
        pos.y += delta.y;
        pos.z += delta.z;
        posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
        posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
    }
}

KERNEL void applySettleAndLangevinMiddlePart3(int numClusters, mixed tol, GLOBAL real4* RESTRICT oldPos,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT velm, GLOBAL const int4* RESTRICT clusterAtoms,
        GLOBAL const float2* RESTRICT clusterParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL real4* RESTRICT posqCorrection
#endif
        , GLOBAL const mixed2* RESTRICT dt
#ifdef FUSE_LANGEVIN_MIDDLE_RESIDUAL_TAIL
        , int numResidualAtoms, GLOBAL const int* RESTRICT residualAtoms, GLOBAL const mixed4* RESTRICT oldDelta
#endif
        ) {
    mixed invDt = 1/dt[0].y;
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    int index = GLOBAL_ID;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        mixed4 apos0 = loadPos(oldPos, posqCorrection, atoms.x);
        mixed4 xp0 = posDelta[atoms.x];
        mixed4 apos1 = loadPos(oldPos, posqCorrection, atoms.y);
        mixed4 xp1 = posDelta[atoms.y];
        mixed4 apos2 = loadPos(oldPos, posqCorrection, atoms.z);
        mixed4 xp2 = posDelta[atoms.z];
#ifndef RELOAD_LANGEVIN_MIDDLE_ORIGINAL_DELTA
        mixed3 original0 = make_mixed3(xp0.x, xp0.y, xp0.z);
        mixed3 original1 = make_mixed3(xp1.x, xp1.y, xp1.z);
        mixed3 original2 = make_mixed3(xp2.x, xp2.y, xp2.z);
#endif
        mixed m0 = 1/velm[atoms.x].w;
        mixed m1 = 1/velm[atoms.y].w;
        mixed m2 = 1/velm[atoms.z].w;

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

        // Record the new positions.

        // Replace the original global mixed4 output boundary with an opaque
        // CUDA f64 register boundary before the unchanged Part3 arithmetic.
#ifdef RELOAD_LANGEVIN_MIDDLE_ORIGINAL_DELTA
        finishSettleMiddleAtom(atoms.x, xp0, posDelta, invDt, oldPos, velm, posqCorrection);
        finishSettleMiddleAtom(atoms.y, xp1, posDelta, invDt, oldPos, velm, posqCorrection);
        finishSettleMiddleAtom(atoms.z, xp2, posDelta, invDt, oldPos, velm, posqCorrection);
#else
        finishSettleMiddleAtom(atoms.x, xp0, original0, invDt, oldPos, velm, posqCorrection);
        finishSettleMiddleAtom(atoms.y, xp1, original1, invDt, oldPos, velm, posqCorrection);
        finishSettleMiddleAtom(atoms.z, xp2, original2, invDt, oldPos, velm, posqCorrection);
#endif
        index += GLOBAL_SIZE;
    }
#ifdef FUSE_LANGEVIN_MIDDLE_RESIDUAL_TAIL
    // SHAKE/CCMA completed on the same stream before this launch. The list is
    // disjoint from SETTLE, so no cross-workgroup barrier is needed here.
    for (int residual = GLOBAL_ID; residual < numResidualAtoms; residual += GLOBAL_SIZE) {
        int atom = residualAtoms[residual];
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[atom];
            velocity.x += (delta.x-oldDelta[atom].x)*invDt;
            velocity.y += (delta.y-oldDelta[atom].y)*invDt;
            velocity.z += (delta.z-oldDelta[atom].z)*invDt;
            velm[atom] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = oldPos[atom];
            real4 pos2 = posqCorrection[atom];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = oldPos[atom];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            oldPos[atom] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[atom] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            oldPos[atom] = pos;
#endif
        }
    }
#endif
}
