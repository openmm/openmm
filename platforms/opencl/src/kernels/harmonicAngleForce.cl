float2 angleParams = PARAMS[index];
float deltaIdeal = theta-angleParams.x;
energy += 0.5f*angleParams.y*deltaIdeal*deltaIdeal;
float dEdAngle = angleParams.y*deltaIdeal;
