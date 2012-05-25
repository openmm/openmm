float4 torsionParams = PARAMS[index];
float deltaAngle = torsionParams.z*theta-torsionParams.y;
energy += torsionParams.x*(1.0f+cos(deltaAngle));
float sinDeltaAngle = sin(deltaAngle);
float dEdAngle = -torsionParams.x*torsionParams.z*sinDeltaAngle;
