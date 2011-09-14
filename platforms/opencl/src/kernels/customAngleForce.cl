float4 v0 = pos2-pos1;
float4 v1 = pos2-pos3;
float4 cp = cross(v0, v1);
float rp = cp.x*cp.x + cp.y*cp.y + cp.z*cp.z;
rp = max(sqrt(rp), 1.0e-06f);
float r21 = v0.x*v0.x + v0.y*v0.y + v0.z*v0.z;
float r23 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
float dot = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
float cosine = clamp(dot/sqrt(r21*r23), -1.0f, 1.0f);
float theta = acos(cosine);
COMPUTE_FORCE
float4 force1 = cross(v0, cp)*(dEdAngle/(r21*rp));
float4 force3 = cross(cp, v1)*(dEdAngle/(r23*rp));
float4 force2 = -force1-force3;
