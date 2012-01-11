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
        float m0 = RECIP(v0.w);
        float m1 = RECIP(v1.w);
        float m2 = RECIP(v2.w);

        // Compute the center of mass position and velocity, and atom positions and
        // velocities relative to it.

        float4 com = (apos0*m0 + apos1*m1 + apos2*m2)/(m0+m1+m2);
        float4 deltax0 = apos0-com;
        float4 deltax1 = apos1-com;
        float4 deltax2 = apos2-com;
        float4 p = v0*m0 + v1*m1 + v2*m2;
        float4 v = p/(m0+m1+m2);
        float4 deltav0 = v0-v;
        float4 deltav1 = v1-v;
        float4 deltav2 = v2-v;

        // Compute linear and angular momentum, and the inertia tensor.

        float4 L = cross(deltax0, deltav0)*m0 + cross(deltax1, deltav1)*m1 + cross(deltax2, deltav2)*m2;
        float Ixx = m0*(deltax0.y*deltax0.y+deltax0.z*deltax0.z) + m1*(deltax1.y*deltax1.y+deltax1.z*deltax1.z) + m2*(deltax2.y*deltax2.y+deltax2.z*deltax2.z);
        float Iyy = m0*(deltax0.x*deltax0.x+deltax0.z*deltax0.z) + m1*(deltax1.x*deltax1.x+deltax1.z*deltax1.z) + m2*(deltax2.x*deltax2.x+deltax2.z*deltax2.z);
        float Izz = m0*(deltax0.x*deltax0.x+deltax0.y*deltax0.y) + m1*(deltax1.x*deltax1.x+deltax1.y*deltax1.y) + m2*(deltax2.x*deltax2.x+deltax2.y*deltax2.y);
        float Ixy = -m0*deltax0.x*deltax0.y - m1*deltax1.x*deltax1.y - m2*deltax2.x*deltax2.y;
        float Ixz = -m0*deltax0.x*deltax0.z - m1*deltax1.x*deltax1.z - m2*deltax2.x*deltax2.z;
        float Iyz = -m0*deltax0.y*deltax0.z - m1*deltax1.y*deltax1.z - m2*deltax2.y*deltax2.z;
        float Iyx = Ixy;
        float Izx = Ixz;
        float Izy = Iyz;
        
        // Invert the inertia tensor and use it to compute the angular velocity.
        
        float Axx = Iyy*Izz - Iyz*Izy;
        float Axy = Iyz*Izx - Iyx*Izz;
        float Axz = Iyx*Izy - Iyy*Izx;
        float Ayx = Izy*Ixz - Izz*Ixy;
        float Ayy = Izz*Ixx - Izx*Ixz;
        float Ayz = Izx*Ixy - Izy*Ixx;
        float Azx = Ixy*Iyz - Ixz*Iyy;
        float Azy = Ixz*Iyx - Ixx*Iyz;
        float Azz = Ixx*Iyy - Ixy*Iyx;
        float invDet = 1.0f/(Ixx*Axx + Ixy*Axy + Ixz*Axz);
        float4 w = (float4) (invDet*(L.x*Axx + L.y*Axy + L.z*Axz),
                             invDet*(L.x*Ayx + L.y*Ayy + L.z*Ayz),
                             invDet*(L.x*Azx + L.y*Azy + L.z*Azz), 0.0f);
        
        // Compute the particle velocities from the molecule's linear and angular velocities.
        
        v0.xyz = v.xyz+cross(w, deltax0).xyz;
        v1.xyz = v.xyz+cross(w, deltax1).xyz;
        v2.xyz = v.xyz+cross(w, deltax2).xyz;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}