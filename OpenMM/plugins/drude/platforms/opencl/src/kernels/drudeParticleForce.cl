real4 delta = (real4) (pos1.xyz-pos2.xyz, 0);
real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
float4 drudeParams = PARAMS[index];
float k1 = drudeParams.x;
float k2 = drudeParams.y;
float k3 = drudeParams.z;

// Compute the isotropic force.

energy += 0.5f*k3*r2;
real4 force1 = -delta*k3;
real4 force2 = delta*k3;
real4 force3 = 0;
real4 force4 = 0;
real4 force5 = 0;

// Compute the first anisotropic force.

if (k1 != 0) {
    real4 dir = (real4) (pos2.xyz-pos3.xyz, 0);
    real invDist = RSQRT(dot(dir, dir));
    dir *= invDist;
    real rprime = dot(dir, delta);
    energy += 0.5f*k1*rprime*rprime;
    real4 f1 = dir*(k1*rprime); 
    real4 f2 = (delta-dir*rprime)*(k1*rprime*invDist);
    force1 -= f1;
    force2 += f1-f2;
    force3 += f2;
}

// Compute the second anisotropic force.

if (k2 != 0) {
    real4 dir = (real4) (pos4.xyz-pos5.xyz, 0);
    real invDist = RSQRT(dot(dir, dir));
    dir *= invDist;
    real rprime = dot(dir, delta);
    energy += 0.5f*k2*rprime*rprime;
    real4 f1 = dir*(k2*rprime);
    real4 f2 = (delta-dir*rprime)*(k2*rprime*invDist);
    force1 -= f1;
    force2 += f1;
    force4 -= f2;
    force5 += f2;
}
