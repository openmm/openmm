const real PI = 3.14159265358979323846f;
real4 v0 = (real4) (pos1.xyz-pos2.xyz, 0.0f);
real4 v1 = (real4) (pos3.xyz-pos2.xyz, 0.0f);
real4 v2 = (real4) (pos3.xyz-pos4.xyz, 0.0f);
#if APPLY_PERIODIC
APPLY_PERIODIC_TO_DELTA(v0)
APPLY_PERIODIC_TO_DELTA(v1)
APPLY_PERIODIC_TO_DELTA(v2)
#endif
real4 cp0 = cross(v0, v1);
real4 cp1 = cross(v1, v2);
real cosangle = dot(normalize(cp0), normalize(cp1));
real theta;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real4 cross_prod = cross(cp0, cp1);
    real scale = dot(cp0, cp0)*dot(cp1, cp1);
    theta = asin(SQRT(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0)
        theta = PI-theta;
}
else
   theta = acos(cosangle);
theta = (dot(v0, cp1) >= 0 ? theta : -theta);
COMPUTE_FORCE
real normCross1 = dot(cp0, cp0);
real normSqrBC = dot(v1, v1);
real normBC = SQRT(normSqrBC);
real normCross2 = dot(cp1, cp1);
real dp = 1.0f/normSqrBC;
real4 ff = (real4) ((-dEdAngle*normBC)/normCross1, dot(v0, v1)*dp, dot(v2, v1)*dp, (dEdAngle*normBC)/normCross2);
real4 force1 = ff.x*cp0;
real4 force4 = ff.w*cp1;
real4 s = ff.y*force1 - ff.z*force4;
real4 force2 = s-force1;
real4 force3 = -s-force4;
