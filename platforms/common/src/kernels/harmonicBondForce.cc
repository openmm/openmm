float2 bondParams = PARAMS[index];
real deltaIdeal = r-bondParams.x;
energy += 0.5f * bondParams.y*deltaIdeal*deltaIdeal;
real dEdR = bondParams.y * deltaIdeal;
