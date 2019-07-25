float2 angleParams = PARAMS[index];
real deltaIdeal = theta*RAD_TO_DEG-angleParams.x;
real deltaIdeal2 = deltaIdeal*deltaIdeal;
real deltaIdeal3 = deltaIdeal*deltaIdeal2;
real deltaIdeal4 = deltaIdeal2*deltaIdeal2;
energy += angleParams.y*deltaIdeal2*(1.0f + CUBIC_K*deltaIdeal + QUARTIC_K*deltaIdeal2 + PENTIC_K*deltaIdeal3 + SEXTIC_K*deltaIdeal4);
real dEdAngle = angleParams.y*deltaIdeal*(2.0f + 3.0f*CUBIC_K*deltaIdeal + 4.0f*QUARTIC_K*deltaIdeal2 + 5.0f*PENTIC_K*deltaIdeal3 + 6.0f*SEXTIC_K*deltaIdeal4);
dEdAngle *= RAD_TO_DEG;