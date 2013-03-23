mixed4 loadPos(__global const real4* restrict posq, __global const real4* restrict posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Enforce constraints on SETTLE clusters
 */

__kernel void applySettle(int numClusters, mixed tol, __global const real4* restrict oldPos, __global const real4* restrict posCorrection, __global mixed4* restrict posDelta, __global const mixed4* restrict velm, __global const int4* restrict clusterAtoms, __global const float2* restrict clusterParams) {
    int index = get_global_id(0);
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        mixed4 apos0 = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 xp0 = posDelta[atoms.x];
        mixed4 apos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 xp1 = posDelta[atoms.y];
        mixed4 apos2 = loadPos(oldPos, posCorrection, atoms.z);
        mixed4 xp2 = posDelta[atoms.z];
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

        mixed invTotalMass = 1.0f/(m0+m1+m2);
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

        float rc = 0.5*params.y;
        mixed rb = sqrt(params.x*params.x-rc*rc);
        mixed ra = rb*(m1+m2)*invTotalMass;
        rb -= ra;
        mixed sinphi = za1d / ra;
        mixed cosphi = sqrt(1.0f - sinphi*sinphi);
        mixed sinpsi = (zb1d - zc1d) / (2*rc*cosphi);
        mixed cospsi = sqrt(1.0f - sinpsi*sinpsi);

        mixed ya2d =   ra*cosphi;
        mixed xb2d = - rc*cospsi;
        mixed yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        mixed yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        mixed xb2d2 = xb2d*xb2d;
        mixed hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        mixed deltx = 2.0f*xb2d + sqrt(4.0f*xb2d2 - hh2 + params.y*params.y);
        xb2d -= deltx*0.5;

        //                                        --- Step3  al,be,ga ---

        mixed alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        mixed beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        mixed gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        mixed al2be2 = alpha*alpha + beta*beta;
        mixed sintheta = (alpha*gamma - beta*sqrt(al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        mixed costheta = sqrt(1.0f - sintheta*sintheta);
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

        posDelta[atoms.x] = xp0;
        posDelta[atoms.y] = xp1;
        posDelta[atoms.z] = xp2;
        index += get_global_size(0);
    }
}

/**
 * Enforce velocity constraints on SETTLE clusters
 */

__kernel void constrainVelocities(int numClusters, mixed tol, __global const real4* restrict oldPos, __global const real4* restrict posCorrection, __global mixed4* restrict posDelta, __global mixed4* restrict velm, __global const int4* restrict clusterAtoms, __global const float2* restrict clusterParams) {
    for (int index = get_global_id(0); index < numClusters; index += get_global_size(0)) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        mixed4 apos0 = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 apos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 apos2 = loadPos(oldPos, posCorrection, atoms.z);
        mixed4 v0 = velm[atoms.x];
        mixed4 v1 = velm[atoms.y];
        mixed4 v2 = velm[atoms.z];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        mixed mA = 1/v0.w;
        mixed mB = 1/v1.w;
        mixed mC = 1/v2.w;
        mixed4 eAB = apos1-apos0;
        mixed4 eBC = apos2-apos1;
        mixed4 eCA = apos0-apos2;
        eAB.xyz /= sqrt(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z);
        eBC.xyz /= sqrt(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z);
        eCA.xyz /= sqrt(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z);
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
        v0.xyz += (tab*eAB.xyz - tca*eCA.xyz)*v0.w;
        v1.xyz += (tbc*eBC.xyz - tab*eAB.xyz)*v1.w;
        v2.xyz += (tca*eCA.xyz - tbc*eBC.xyz)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}