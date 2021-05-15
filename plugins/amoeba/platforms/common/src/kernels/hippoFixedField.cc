real invR2 = invR*invR;
real invR3 = invR*invR2;
real invR5 = invR3*invR2;
real invR7 = invR5*invR2;

#if USE_EWALD
// Calculate the error function damping terms.

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
alsq2n *= alsq2;
real bn3 = (5*bn2+alsq2n*exp2a)*invR2;
#endif

// Calculate the field at particle 1 due to multipoles at particle 2

real fdamp3, fdamp5, fdamp7;
computeDirectFieldDampingFactors(alpha2, r, fdamp3, fdamp5, fdamp7);
#ifndef COMPUTING_EXCEPTIONS
real scale = 1;
#endif
#ifdef USE_EWALD
real rr3 = bn1 - (1-scale)*invR3;
real rr3j = bn1 - (1-scale*fdamp3)*invR3;
real rr5j = bn2 - (1-scale*fdamp5)*3*invR5;
real rr7j = bn3 - (1-scale*fdamp7)*15*invR7;
#else
real rr3 = scale*invR3;
real rr3j = scale*fdamp3*invR3;
real rr5j = scale*3*fdamp5*invR5;
real rr7j = scale*15*fdamp7*invR7;
#endif
real qZZ2 = -qXX2-qYY2;
real3 qDotDelta2 = make_real3(delta.x*qXX2 + delta.y*qXY2 + delta.z*qXZ2,
                              delta.x*qXY2 + delta.y*qYY2 + delta.z*qYZ2,
                              delta.x*qXZ2 + delta.y*qYZ2 + delta.z*qZZ2);
real dipoleDelta2 = dot(dipole2, delta);
real qdpoleDelta2 = dot(qDotDelta2, delta);
real factor2 = rr3*coreCharge2 + rr3j*valenceCharge2 - rr5j*dipoleDelta2 + rr7j*qdpoleDelta2;
tempField1 = -delta*factor2 - dipole2*rr3j + qDotDelta2*2*rr5j;

// Calculate the field at particle 2 due to multipoles at particle 1

computeDirectFieldDampingFactors(alpha1, r, fdamp3, fdamp5, fdamp7);
#ifdef USE_EWALD
real rr3i = bn1 - (1-scale*fdamp3)*invR3;
real rr5i = bn2 - (1-scale*fdamp5)*3*invR5;
real rr7i = bn3 - (1-scale*fdamp7)*15*invR7;
#else
real rr3i = scale*fdamp3*invR3;
real rr5i = scale*3*fdamp5*invR5;
real rr7i = scale*15*fdamp7*invR7;
#endif
real qZZ1 = -qXX1-qYY1;
real3 qDotDelta1 = make_real3(delta.x*qXX1 + delta.y*qXY1 + delta.z*qXZ1,
                              delta.x*qXY1 + delta.y*qYY1 + delta.z*qYZ1,
                              delta.x*qXZ1 + delta.y*qYZ1 + delta.z*qZZ1);
real dipoleDelta1 = dot(dipole1, delta);
real qdpoleDelta1 = dot(qDotDelta1, delta);
real factor1 = rr3*coreCharge1 + rr3i*valenceCharge1 + rr5i*dipoleDelta1 + rr7i*qdpoleDelta1;
tempField2 = delta*factor1 - dipole1*rr3i - qDotDelta1*2*rr5i;
