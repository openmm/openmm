const real exclusionScale = PARAMS[index];
real3 delta = make_real3(pos2.x - pos1.x, pos2.y - pos1.y, pos2.z - pos1.z);
#if APPLY_PERIODIC
    APPLY_PERIODIC_TO_DELTA(delta)
#endif

const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
const real r = SQRT(r2);
const real invR = RECIP(r);

const real alphaR = EWALD_ALPHA * r;
real tempForce = 0.0f;
if (alphaR > 1e-6f) {
    const real erfAlphaR = ERF(alphaR);
    const real prefactor = exclusionScale * invR;
    tempForce = prefactor * (erfAlphaR - TWO_OVER_SQRT_PI * alphaR * EXP(-alphaR * alphaR)) * invR * invR;

    energy -= prefactor * erfAlphaR;
}
else {
    energy -= TWO_OVER_SQRT_PI * EWALD_ALPHA * exclusionScale;
}
delta *= tempForce;
real3 force1 = delta;
real3 force2 = -delta;
