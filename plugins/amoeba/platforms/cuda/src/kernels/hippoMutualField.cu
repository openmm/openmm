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
real bn0 = erfc(ralpha)*invR;
real alsq2 = 2*PME_ALPHA*PME_ALPHA;
real alsq2n = 1/(SQRT_PI*PME_ALPHA);
real exp2a = EXP(-(ralpha*ralpha));
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
