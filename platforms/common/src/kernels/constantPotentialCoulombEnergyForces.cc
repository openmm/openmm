// The approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
// the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
// error of 1.5e-7.

if (!isExcluded && r2 < CUTOFF_SQUARED) {
    const real prefactor = ONE_4PI_EPS0 * CHARGE1 * CHARGE2 * invR;

    const real alphaR = EWALD_ALPHA * r;
    const real expAlphaRSqr = EXP(-alphaR * alphaR);
#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(alphaR);
#else
    const real tAlpha = RECIP(1.0f+0.3275911f*alphaR);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*tAlpha)*tAlpha)*tAlpha)*tAlpha)*tAlpha*expAlphaRSqr;
#endif

    real tempForceScale = erfcAlphaR + TWO_OVER_SQRT_PI * alphaR * expAlphaRSqr;
    real tempEnergyScale = erfcAlphaR;

    if (SYSELEC1 != -1 || SYSELEC2 != -1) {
        const real4 params1 = PARAMS[SYSELEC1 + 1];
        const real4 params2 = PARAMS[SYSELEC2 + 1];

        const real etaR = r / SQRT(params1.y * params1.y + params2.y * params2.y);
        const real expEtaRSqr = EXP(-etaR * etaR);
#ifdef USE_DOUBLE_PRECISION
        const real erfcEtaR = erfc(etaR);
#else
        const real tEta = RECIP(1.0f+0.3275911f*etaR);
        const real erfcEtaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*tEta)*tEta)*tEta)*tEta)*tEta*expEtaRSqr;
#endif

        tempForceScale -= erfcEtaR + TWO_OVER_SQRT_PI * etaR * expEtaRSqr;
        tempEnergyScale -= erfcEtaR;
    }

    tempEnergy += prefactor * tempEnergyScale;
    dEdR += prefactor * tempForceScale * invR * invR;
}
