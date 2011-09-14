float4 delta = pos2-pos1;
float r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
COMPUTE_FORCE
delta.xyz *= -dEdR/r;
float4 force1 = -delta;
float4 force2 = delta;
