const float PI = 3.14159265358979323846f;
float4 v0 = (float4) (pos1.xyz-pos2.xyz, 0.0f);
float4 v1 = (float4) (pos3.xyz-pos2.xyz, 0.0f);
float4 v2 = (float4) (pos3.xyz-pos4.xyz, 0.0f);
float4 cp0 = cross(v0, v1);
float4 cp1 = cross(v1, v2);
float cosangle = dot(normalize(cp0), normalize(cp1));
float theta;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    float4 cross_prod = cross(cp0, cp1);
    float scale = dot(cp0, cp0)*dot(cp1, cp1);
    theta = asin(sqrt(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0.0f)
        theta = PI-theta;
}
else
   theta = acos(cosangle);
theta = (dot(v0, cp1) >= 0 ? theta : -theta);
COMPUTE_FORCE
float normCross1 = dot(cp0, cp0);
float normSqrBC = dot(v1, v1);
float normBC = sqrt(normSqrBC);
float normCross2 = dot(cp1, cp1);
float dp = 1.0f/normSqrBC;
float4 ff = (float4) ((-dEdAngle*normBC)/normCross1, dot(v0, v1)*dp, dot(v2, v1)*dp, (dEdAngle*normBC)/normCross2);
float4 force1 = ff.x*cp0;
float4 force4 = ff.w*cp1;
float4 s = ff.y*force1 - ff.z*force4;
float4 force2 = s-force1;
float4 force3 = -s-force4;
