float2 angleParams = PARAMS[index];
real deltaIdeal = theta-angleParams.x;
energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal;
real dEdAngle = angleParams.y*deltaIdeal;
