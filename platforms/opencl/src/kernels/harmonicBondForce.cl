float2 bondParams = PARAMS[index];
float deltaIdeal = r-bondParams.x;
energy += 0.5f * bondParams.y*deltaIdeal*deltaIdeal;
float dEdR = bondParams.y * deltaIdeal;
