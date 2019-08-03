real fdamp3, fdamp5;
computeMutualFieldDampingFactors(alpha1, alpha2, r, fdamp3, fdamp5);
#ifdef COMPUTING_EXCEPTIONS
fdamp3 *= scale;
fdamp5 *= scale;
#endif
real invR2 = invR*invR;
real invR3 = invR*invR2;
#if USE_EWALD
real ralpha = PME_ALPHA*r;
real exp2a = EXP(-(ralpha*ralpha));
#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(ralpha);
#else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    const real t = RECIP(1.0f+0.3275911f*ralpha);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
real bn0 = erfcAlphaR*invR;
real alsq2 = 2*PME_ALPHA*PME_ALPHA;
real alsq2n = 1/(SQRT_PI*PME_ALPHA);
alsq2n *= alsq2;
real bn1 = (bn0+alsq2n*exp2a)*invR2;
alsq2n *= alsq2;
real bn2 = (3*bn1+alsq2n*exp2a)*invR2;
real scale3 = -bn1 + (1-fdamp3)*invR3;
real scale5 = bn2 - 3*(1-fdamp5)*invR3*invR2;
#else
real scale3 = -fdamp3*invR3;
real scale5 = 3*fdamp5*invR3*invR2;
#endif
tempField1 = inducedDipole2*scale3 + delta*scale5*dot(inducedDipole2, delta);
tempField2 = inducedDipole1*scale3 + delta*scale5*dot(inducedDipole1, delta);
