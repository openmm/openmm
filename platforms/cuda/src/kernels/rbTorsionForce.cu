float4 torsionParams1 = PARAMS1[index];
float2 torsionParams2 = PARAMS2[index];
if (theta < 0)
    theta += PI;
else
    theta -= PI;
cosangle = -cosangle;
real cosFactor = cosangle;
real dEdAngle = -torsionParams1.y;
real rbEnergy = torsionParams1.x;
rbEnergy += torsionParams1.y*cosFactor;
dEdAngle -= 2.0f*torsionParams1.z*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 3.0f*torsionParams1.w*cosFactor;
rbEnergy += torsionParams1.z*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 4.0f*torsionParams2.x*cosFactor;
rbEnergy += torsionParams1.w*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 5.0f*torsionParams2.y*cosFactor;
rbEnergy += torsionParams2.x*cosFactor;
rbEnergy += torsionParams2.y*cosFactor*cosangle;
energy += rbEnergy;
dEdAngle *= SIN(theta);
