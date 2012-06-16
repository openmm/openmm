/**
 * Enforce constraints on SETTLE clusters
 */

extern "C" __global__ void applySettle(int numClusters, float tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, const real4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        real4 apos0 = oldPos[atoms.x];
        real4 xp0 = posDelta[atoms.x];
        real4 apos1 = oldPos[atoms.y];
        real4 xp1 = posDelta[atoms.y];
        real4 apos2 = oldPos[atoms.z];
        real4 xp2 = posDelta[atoms.z];
        real m0 = RECIP(velm[atoms.x].w);
        real m1 = RECIP(velm[atoms.y].w);
        real m2 = RECIP(velm[atoms.z].w);

        // Apply the SETTLE algorithm.

        real xb0 = apos1.x-apos0.x;
        real yb0 = apos1.y-apos0.y;
        real zb0 = apos1.z-apos0.z;
        real xc0 = apos2.x-apos0.x;
        real yc0 = apos2.y-apos0.y;
        real zc0 = apos2.z-apos0.z;

        real invTotalMass = RECIP(m0+m1+m2);
        real xcom = (xp0.x*m0 + (xb0+xp1.x)*m1 + (xc0+xp2.x)*m2) * invTotalMass;
        real ycom = (xp0.y*m0 + (yb0+xp1.y)*m1 + (yc0+xp2.y)*m2) * invTotalMass;
        real zcom = (xp0.z*m0 + (zb0+xp1.z)*m1 + (zc0+xp2.z)*m2) * invTotalMass;

        real xa1 = xp0.x - xcom;
        real ya1 = xp0.y - ycom;
        real za1 = xp0.z - zcom;
        real xb1 = xb0 + xp1.x - xcom;
        real yb1 = yb0 + xp1.y - ycom;
        real zb1 = zb0 + xp1.z - zcom;
        real xc1 = xc0 + xp2.x - xcom;
        real yc1 = yc0 + xp2.y - ycom;
        real zc1 = zc0 + xp2.z - zcom;

        real xaksZd = yb0*zc0 - zb0*yc0;
        real yaksZd = zb0*xc0 - xb0*zc0;
        real zaksZd = xb0*yc0 - yb0*xc0;
        real xaksXd = ya1*zaksZd - za1*yaksZd;
        real yaksXd = za1*xaksZd - xa1*zaksZd;
        real zaksXd = xa1*yaksZd - ya1*xaksZd;
        real xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
        real yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
        real zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

        real axlng = SQRT(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
        real aylng = SQRT(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
        real azlng = SQRT(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
        real trns11 = xaksXd / axlng;
        real trns21 = yaksXd / axlng;
        real trns31 = zaksXd / axlng;
        real trns12 = xaksYd / aylng;
        real trns22 = yaksYd / aylng;
        real trns32 = zaksYd / aylng;
        real trns13 = xaksZd / azlng;
        real trns23 = yaksZd / azlng;
        real trns33 = zaksZd / azlng;

        real xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
        real yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
        real xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
        real yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
        real za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
        real xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
        real yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
        real zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
        real xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
        real yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
        real zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

        //                                        --- Step2  A2' ---

        float rc = 0.5f*params.y;
        float rb = SQRT(params.x*params.x-rc*rc);
        real ra = rb*(m1+m2)*invTotalMass;
        rb -= ra;
        real sinphi = za1d/ra;
        real cosphi = SQRT(1-sinphi*sinphi);
        real sinpsi = (zb1d-zc1d) / (2*rc*cosphi);
        real cospsi = SQRT(1-sinpsi*sinpsi);

        real ya2d =   ra*cosphi;
        real xb2d = - rc*cospsi;
        real yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        real yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        real xb2d2 = xb2d*xb2d;
        real hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        real deltx = 2.0f*xb2d + SQRT(4.0f*xb2d2 - hh2 + params.y*params.y);
        xb2d -= deltx*0.5f;

        //                                        --- Step3  al,be,ga ---

        real alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        real beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        real gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        real al2be2 = alpha*alpha + beta*beta;
        real sintheta = (alpha*gamma - beta*SQRT(al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        real costheta = SQRT(1-sintheta*sintheta);
        real xa3d = - ya2d*sintheta;
        real ya3d =   ya2d*costheta;
        real za3d = za1d;
        real xb3d =   xb2d*costheta - yb2d*sintheta;
        real yb3d =   xb2d*sintheta + yb2d*costheta;
        real zb3d = zb1d;
        real xc3d = - xb2d*costheta - yc2d*sintheta;
        real yc3d = - xb2d*sintheta + yc2d*costheta;
        real zc3d = zc1d;

        //                                        --- Step5  A3 ---

        real xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
        real ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
        real za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
        real xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
        real yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
        real zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
        real xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
        real yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
        real zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

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
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Enforce velocity constraints on SETTLE clusters
 */

extern "C" __global__ void constrainVelocities(int numClusters, float tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, real4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numClusters; index += blockDim.x*gridDim.x) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        real4 apos0 = oldPos[atoms.x];
        real4 apos1 = oldPos[atoms.y];
        real4 apos2 = oldPos[atoms.z];
        real4 v0 = velm[atoms.x];
        real4 v1 = velm[atoms.y];
        real4 v2 = velm[atoms.z];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        real mA = RECIP(v0.w);
        real mB = RECIP(v1.w);
        real mC = RECIP(v2.w);
        real3 eAB = make_real3(apos1.x-apos0.x, apos1.y-apos0.y, apos1.z-apos0.z);
        real3 eBC = make_real3(apos2.x-apos1.x, apos2.y-apos1.y, apos2.z-apos1.z);
        real3 eCA = make_real3(apos0.x-apos2.x, apos0.y-apos2.y, apos0.z-apos2.z);
        eAB *= RECIP(SQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z));
        eBC *= RECIP(SQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z));
        eCA *= RECIP(SQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z));
        real vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
        real vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
        real vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
        real cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
        real cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
        real cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
        real s2A = 1-cA*cA;
        real s2B = 1-cB*cB;
        real s2C = 1-cC*cC;
        
        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.
        
        real mABCinv = RECIP(mA*mB*mC);
        real denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
        real tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
        real tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        real tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0.x += (tab*eAB.x - tca*eCA.x)*v0.w;
        v0.y += (tab*eAB.y - tca*eCA.y)*v0.w;
        v0.z += (tab*eAB.z - tca*eCA.z)*v0.w;
        v1.x += (tbc*eBC.x - tab*eAB.x)*v1.w;
        v1.y += (tbc*eBC.y - tab*eAB.y)*v1.w;
        v1.z += (tbc*eBC.z - tab*eAB.z)*v1.w;
        v2.x += (tca*eCA.x - tbc*eBC.x)*v2.w;
        v2.y += (tca*eCA.y - tbc*eBC.y)*v2.w;
        v2.z += (tca*eCA.z - tbc*eBC.z)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}