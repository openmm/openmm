real3 delta = make_real3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z);
real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
float4 drudeParams = PARAMS[index];
float k1 = drudeParams.x;
float k2 = drudeParams.y;
float k3 = drudeParams.z;

// Compute the isotropic force.

energy += 0.5f*k3*r2;
real3 force1 = -delta*k3;
real3 force2 = delta*k3;
real3 force3 = make_real3(0);
real3 force4 = make_real3(0);
real3 force5 = make_real3(0);

// Compute the first anisotropic force.

if (k1 != 0) {
    real3 dir = make_real3(pos2.x-pos3.x, pos2.y-pos3.y, pos2.z-pos3.z);
    real invDist = RSQRT(dot(dir, dir));
    dir *= invDist;
    real rprime = dot(dir, delta);
    energy += 0.5f*k1*rprime*rprime;
    real3 f1 = dir*(k1*rprime); 
    real3 f2 = (delta-dir*rprime)*(k1*rprime*invDist);
    force1 -= f1;
    force2 += f1-f2;
    force3 += f2;
}

// Compute the second anisotropic force.

if (k2 != 0) {
    real3 dir = make_real3(pos4.x-pos5.x, pos4.y-pos5.y, pos4.z-pos5.z);
    real invDist = RSQRT(dot(dir, dir));
    dir *= invDist;
    real rprime = dot(dir, delta);
    energy += 0.5f*k2*rprime*rprime;
    real3 f1 = dir*(k2*rprime);
    real3 f2 = (delta-dir*rprime)*(k2*rprime*invDist);
    force1 -= f1;
    force2 += f1;
    force4 -= f2;
    force5 += f2;
}
