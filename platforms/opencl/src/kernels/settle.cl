/**
 * Enforce constraints on SETTLE clusters
 */

__kernel void applySettle(int numClusters, float tol, __global const float4* restrict oldPos, __global float4* restrict posDelta, __global const float4* restrict velm, __global const int4* restrict clusterAtoms, __global const float2* restrict clusterParams) {
    int index = get_global_id(0);
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        float4 apos0 = oldPos[atoms.x];
        float4 xp0 = posDelta[atoms.x];
        float4 apos1 = oldPos[atoms.y];
        float4 xp1 = posDelta[atoms.y];
        float4 apos2 = oldPos[atoms.z];
        float4 xp2 = posDelta[atoms.z];
        float m0 = RECIP(velm[atoms.x].w);
        float m1 = RECIP(velm[atoms.y].w);
        float m2 = RECIP(velm[atoms.z].w);

        // Translate the molecule to the origin to improve numerical precision.

        float4 center = apos0;
        apos0.xyz -= center.xyz;
        apos1.xyz -= center.xyz;
        apos2.xyz -= center.xyz;

        // Apply the SETTLE algorithm.

        float xb0 = apos1.x-apos0.x;
        float yb0 = apos1.y-apos0.y;
        float zb0 = apos1.z-apos0.z;
        float xc0 = apos2.x-apos0.x;
        float yc0 = apos2.y-apos0.y;
        float zc0 = apos2.z-apos0.z;

        float totalMass = m0+m1+m2;
        float xcom = ((apos0.x+xp0.x)*m0 + (apos1.x+xp1.x)*m1 + (apos2.x+xp2.x)*m2) / totalMass;
        float ycom = ((apos0.y+xp0.y)*m0 + (apos1.y+xp1.y)*m1 + (apos2.y+xp2.y)*m2) / totalMass;
        float zcom = ((apos0.z+xp0.z)*m0 + (apos1.z+xp1.z)*m1 + (apos2.z+xp2.z)*m2) / totalMass;

        float xa1 = apos0.x + xp0.x - xcom;
        float ya1 = apos0.y + xp0.y - ycom;
        float za1 = apos0.z + xp0.z - zcom;
        float xb1 = apos1.x + xp1.x - xcom;
        float yb1 = apos1.y + xp1.y - ycom;
        float zb1 = apos1.z + xp1.z - zcom;
        float xc1 = apos2.x + xp2.x - xcom;
        float yc1 = apos2.y + xp2.y - ycom;
        float zc1 = apos2.z + xp2.z - zcom;

        float xaksZd = yb0*zc0 - zb0*yc0;
        float yaksZd = zb0*xc0 - xb0*zc0;
        float zaksZd = xb0*yc0 - yb0*xc0;
        float xaksXd = ya1*zaksZd - za1*yaksZd;
        float yaksXd = za1*xaksZd - xa1*zaksZd;
        float zaksXd = xa1*yaksZd - ya1*xaksZd;
        float xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
        float yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
        float zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

        float axlng = sqrt(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
        float aylng = sqrt(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
        float azlng = sqrt(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
        float trns11 = xaksXd / axlng;
        float trns21 = yaksXd / axlng;
        float trns31 = zaksXd / axlng;
        float trns12 = xaksYd / aylng;
        float trns22 = yaksYd / aylng;
        float trns32 = zaksYd / aylng;
        float trns13 = xaksZd / azlng;
        float trns23 = yaksZd / azlng;
        float trns33 = zaksZd / azlng;

        float xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
        float yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
        float xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
        float yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
        float za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
        float xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
        float yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
        float zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
        float xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
        float yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
        float zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

        //                                        --- Step2  A2' ---

        float rc = 0.5*params.y;
        float rb = sqrt(params.x*params.x-rc*rc);
        float ra = rb*(m1+m2)/totalMass;
        rb -= ra;
        float sinphi = za1d / ra;
        float cosphi = sqrt(1.0f - sinphi*sinphi);
        float sinpsi = (zb1d - zc1d) / (2*rc*cosphi);
        float cospsi = sqrt(1.0f - sinpsi*sinpsi);

        float ya2d =   ra*cosphi;
        float xb2d = - rc*cospsi;
        float yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        float yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        float xb2d2 = xb2d*xb2d;
        float hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        float deltx = 2.0f*xb2d + sqrt(4.0f*xb2d2 - hh2 + params.y*params.y);
        xb2d -= deltx*0.5;

        //                                        --- Step3  al,be,ga ---

        float alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        float beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        float gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        float al2be2 = alpha*alpha + beta*beta;
        float sintheta = (alpha*gamma - beta*sqrt (al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        float costheta = sqrt (1.0f - sintheta*sintheta);
        float xa3d = - ya2d*sintheta;
        float ya3d =   ya2d*costheta;
        float za3d = za1d;
        float xb3d =   xb2d*costheta - yb2d*sintheta;
        float yb3d =   xb2d*sintheta + yb2d*costheta;
        float zb3d = zb1d;
        float xc3d = - xb2d*costheta - yc2d*sintheta;
        float yc3d = - xb2d*sintheta + yc2d*costheta;
        float zc3d = zc1d;

        //                                        --- Step5  A3 ---

        float xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
        float ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
        float za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
        float xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
        float yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
        float zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
        float xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
        float yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
        float zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

        xp0.x = xcom + xa3 - apos0.x;
        xp0.y = ycom + ya3 - apos0.y;
        xp0.z = zcom + za3 - apos0.z;
        xp1.x = xcom + xb3 - apos1.x;
        xp1.y = ycom + yb3 - apos1.y;
        xp1.z = zcom + zb3 - apos1.z;
        xp2.x = xcom + xc3 - apos2.x;
        xp2.y = ycom + yc3 - apos2.y;
        xp2.z = zcom + zc3 - apos2.z;

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

__kernel void constrainVelocities(int numClusters, float tol, __global const float4* restrict oldPos, __global float4* restrict posDelta, __global float4* restrict velm, __global const int4* restrict clusterAtoms, __global const float2* restrict clusterParams) {
    for (int index = get_global_id(0); index < numClusters; index++) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 apos0 = oldPos[atoms.x];
        float4 apos1 = oldPos[atoms.y];
        float4 apos2 = oldPos[atoms.z];
        float4 v0 = velm[atoms.x];
        float4 v1 = velm[atoms.y];
        float4 v2 = velm[atoms.z];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        float mA = RECIP(v0.w);
        float mB = RECIP(v1.w);
        float mC = RECIP(v2.w);
        float4 eAB = apos1-apos0;
        float4 eBC = apos2-apos1;
        float4 eCA = apos0-apos2;
        eAB.xyz /= SQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z);
        eBC.xyz /= SQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z);
        eCA.xyz /= SQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z);
        float vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
        float vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
        float vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
        float cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
        float cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
        float cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
        float s2A = 1-cA*cA;
        float s2B = 1-cB*cB;
        float s2C = 1-cC*cC;
        
        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.
        
        float mABC = mA*mB*mC;
        float denom = ((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB);
        float tab = (mABC*(cB*cC*mA-cA*mB-cA*mC)*vCA + mABC*(cA*cC*mB-cB*mC-cB*mA)*vBC + (mABC*(mA+mB+mC)+s2C*mA*mA*mB*mB)*vAB)/denom;
        float tbc = (mABC*(cA*cB*mC-cC*mB-cC*mA)*vCA + ((s2A*mB*mB)*mC*mC+mABC*(mA+mB+mC))*vBC + mABC*(cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        float tca = (((s2B*mA*mA)*mC*mC+mABC*(mA+mB+mC))*vCA + mABC*(cA*cB*mC-cC*mB-cC*mA)*vBC + mABC*(cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0.xyz += (tab*eAB.xyz - tca*eCA.xyz)*v0.w;
        v1.xyz += (tbc*eBC.xyz - tab*eAB.xyz)*v1.w;
        v2.xyz += (tca*eCA.xyz - tbc*eBC.xyz)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}