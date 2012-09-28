float4 torsionParams = PARAMS[index];
real deltaAngle = torsionParams.z*theta-torsionParams.y;
energy += torsionParams.x*(1.0f+COS(deltaAngle));
real sinDeltaAngle = SIN(deltaAngle);
real dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle;
