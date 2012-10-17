float4 torsionParams = PARAMS[index];
real deltaAngle = torsionParams.z*theta-torsionParams.y;
energy += torsionParams.x*(1.0f+cos(deltaAngle));
real sinDeltaAngle = sin(deltaAngle);
real dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle;
