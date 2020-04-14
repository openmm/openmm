const real PI = (real) 3.14159265358979323846;

// Compute the first angle.

real3 v0a = make_real3(pos1.x-pos2.x, pos1.y-pos2.y, pos1.z-pos2.z);
real3 v1a = make_real3(pos3.x-pos2.x, pos3.y-pos2.y, pos3.z-pos2.z);
real3 v2a = make_real3(pos3.x-pos4.x, pos3.y-pos4.y, pos3.z-pos4.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(v0a)
APPLY_PERIODIC_TO_DELTA(v1a)
APPLY_PERIODIC_TO_DELTA(v2a)
#endif
real3 cp0a = cross(v0a, v1a);
real3 cp1a = cross(v1a, v2a);
real cosangle = dot(normalize(cp0a), normalize(cp1a));
real angleA;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real3 cross_prod = cross(cp0a, cp1a);
    real scale = dot(cp0a, cp0a)*dot(cp1a, cp1a);
    angleA = ASIN(SQRT(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0.0f)
        angleA = PI-angleA;
}
else
   angleA = ACOS(cosangle);
angleA = (dot(v0a, cp1a) >= 0 ? angleA : -angleA);
angleA = fmod(angleA+2.0f*PI, 2.0f*PI);

// Compute the second angle.

real3 v0b = make_real3(pos5.x-pos6.x, pos5.y-pos6.y, pos5.z-pos6.z);
real3 v1b = make_real3(pos7.x-pos6.x, pos7.y-pos6.y, pos7.z-pos6.z);
real3 v2b = make_real3(pos7.x-pos8.x, pos7.y-pos8.y, pos7.z-pos8.z);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(v0b)
APPLY_PERIODIC_TO_DELTA(v1b)
APPLY_PERIODIC_TO_DELTA(v2b)
#endif
real3 cp0b = cross(v0b, v1b);
real3 cp1b = cross(v1b, v2b);
cosangle = dot(normalize(cp0b), normalize(cp1b));
real angleB;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real3 cross_prod = cross(cp0b, cp1b);
    real scale = dot(cp0b, cp0b)*dot(cp1b, cp1b);
    angleB = ASIN(SQRT(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0.0f)
        angleB = PI-angleB;
}
else
   angleB = ACOS(cosangle);
angleB = (dot(v0b, cp1b) >= 0 ? angleB : -angleB);
angleB = fmod(angleB+2.0f*PI, 2.0f*PI);

// Identify which patch this is in.

int2 pos = MAP_POS[MAPS[index]];
int size = pos.y;
real delta = 2*PI/size;
int s = (int) fmin(angleA/delta, size-1);
int t = (int) fmin(angleB/delta, size-1);
float4 c[4];
int coeffIndex = pos.x+4*(s+size*t);
c[0] = COEFF[coeffIndex];
c[1] = COEFF[coeffIndex+1];
c[2] = COEFF[coeffIndex+2];
c[3] = COEFF[coeffIndex+3];
real da = angleA/delta-s;
real db = angleB/delta-t;

// Evaluate the spline to determine the energy and gradients.

real torsionEnergy = 0.0f;
real dEdA = 0.0f;
real dEdB = 0.0f;
torsionEnergy = da*torsionEnergy + ((c[3].w*db + c[3].z)*db + c[3].y)*db + c[3].x;
dEdA = db*dEdA + (3.0f*c[3].w*da + 2.0f*c[2].w)*da + c[1].w;
dEdB = da*dEdB + (3.0f*c[3].w*db + 2.0f*c[3].z)*db + c[3].y;
torsionEnergy = da*torsionEnergy + ((c[2].w*db + c[2].z)*db + c[2].y)*db + c[2].x;
dEdA = db*dEdA + (3.0f*c[3].z*da + 2.0f*c[2].z)*da + c[1].z;
dEdB = da*dEdB + (3.0f*c[2].w*db + 2.0f*c[2].z)*db + c[2].y;
torsionEnergy = da*torsionEnergy + ((c[1].w*db + c[1].z)*db + c[1].y)*db + c[1].x;
dEdA = db*dEdA + (3.0f*c[3].y*da + 2.0f*c[2].y)*da + c[1].y;
dEdB = da*dEdB + (3.0f*c[1].w*db + 2.0f*c[1].z)*db + c[1].y;
torsionEnergy = da*torsionEnergy + ((c[0].w*db + c[0].z)*db + c[0].y)*db + c[0].x;
dEdA = db*dEdA + (3.0f*c[3].x*da + 2.0f*c[2].x)*da + c[1].x;
dEdB = da*dEdB + (3.0f*c[0].w*db + 2.0f*c[0].z)*db + c[0].y;
dEdA /= delta;
dEdB /= delta;
energy += torsionEnergy;

// Apply the force to the first torsion.

real normCross1 = dot(cp0a, cp0a);
real normSqrBC = dot(v1a, v1a);
real normBC = SQRT(normSqrBC);
real normCross2 = dot(cp1a, cp1a);
real dp = RECIP(normSqrBC);
real4 ff = make_real4((-dEdA*normBC)/normCross1, dot(v0a, v1a)*dp, dot(v2a, v1a)*dp, (dEdA*normBC)/normCross2);
real3 force1 = ff.x*cp0a;
real3 force4 = ff.w*cp1a;
real3 d = ff.y*force1 - ff.z*force4;
real3 force2 = d-force1;
real3 force3 = -d-force4;

// Apply the force to the second torsion.

normCross1 = dot(cp0b, cp0b);
normSqrBC = dot(v1b, v1b);
normBC = SQRT(normSqrBC);
normCross2 = dot(cp1b, cp1b);
dp = RECIP(normSqrBC);
ff = make_real4((-dEdB*normBC)/normCross1, dot(v0b, v1b)*dp, dot(v2b, v1b)*dp, (dEdB*normBC)/normCross2);
real3 force5 = ff.x*cp0b;
real3 force8 = ff.w*cp1b;
d = ff.y*force5 - ff.z*force8;
real3 force6 = d-force5;
real3 force7 = -d-force8;
