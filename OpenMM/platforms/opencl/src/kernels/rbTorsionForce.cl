float8 torsionParams = PARAMS[index];
if (theta < 0.0f)
    theta += PI;
else
    theta -= PI;
cosangle = -cosangle;
real cosFactor = cosangle;
real dEdAngle = -torsionParams.s1;
real rbEnergy = torsionParams.s0;
rbEnergy += torsionParams.s1*cosFactor;
dEdAngle -= 2.0f*torsionParams.s2*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 3.0f*torsionParams.s3*cosFactor;
rbEnergy += torsionParams.s2*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 4.0f*torsionParams.s4*cosFactor;
rbEnergy += torsionParams.s3*cosFactor;
cosFactor *= cosangle;
dEdAngle -= 5.0f*torsionParams.s5*cosFactor;
rbEnergy += torsionParams.s4*cosFactor;
rbEnergy += torsionParams.s5*cosFactor*cosangle;
energy += rbEnergy;
dEdAngle *= sin(theta);
