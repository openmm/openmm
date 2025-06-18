const real exceptionScale = PARAMS[index];
real3 delta = make_real3(pos2.x - pos1.x, pos2.y - pos1.y, pos2.z - pos1.z);
#if APPLY_PERIODIC
    APPLY_PERIODIC_TO_DELTA(delta)
#endif

const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
const real invR = RSQRT(r2);

const real tempEnergy = exceptionScale * invR;
const real tempForce = tempEnergy * invR * invR;

energy += tempEnergy;
delta *= tempForce;
real3 force1 = -delta;
real3 force2 = delta;
