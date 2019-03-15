real invR2 = invR*invR;
real invR3 = invR*invR2;
real invR5 = invR3*invR2;
real invR7 = invR5*invR2;

// field at particle 1 due to multipoles at particle 2

real qZZ2 = -qXX2-qYY2;
real3 qDotDelta2 = make_real3(delta.x*qXX2 + delta.y*qXY2 + delta.z*qXZ2,
                              delta.x*qXY2 + delta.y*qYY2 + delta.z*qYZ2,
                              delta.x*qXZ2 + delta.y*qYZ2 + delta.z*qZZ2);
real dipoleDelta2 = dot(dipole2, delta);
real qdpoleDelta2 = dot(qDotDelta2, delta);
real fdamp3, fdamp5, fdamp7;
computeDirectFieldDampingFactors(alpha2, r, fdamp3, fdamp5, fdamp7);
real factor2 = invR3*coreCharge2 + fdamp3*invR3*valenceCharge2 - 3*fdamp5*invR5*dipoleDelta2 + 15*fdamp7*invR7*qdpoleDelta2;
tempField1 = -delta*factor2 - dipole2*fdamp3*invR3 + qDotDelta2*6*fdamp5*invR5;

// field at particle 2 due to multipoles at particle 1

real qZZ1 = -qXX1-qYY1;
real3 qDotDelta1 = make_real3(delta.x*qXX1 + delta.y*qXY1 + delta.z*qXZ1,
                              delta.x*qXY1 + delta.y*qYY1 + delta.z*qYZ1,
                              delta.x*qXZ1 + delta.y*qYZ1 + delta.z*qZZ1);
real dipoleDelta1 = dot(dipole1, delta);
real qdpoleDelta1 = dot(qDotDelta1, delta);
computeDirectFieldDampingFactors(alpha1, r, fdamp3, fdamp5, fdamp7);
real factor1 = invR3*coreCharge1 + fdamp3*invR3*valenceCharge1 + 3*fdamp5*invR5*dipoleDelta1 + 15*fdamp7*invR7*qdpoleDelta1;
tempField2 = delta*factor1 - dipole1*fdamp3*invR3 - qDotDelta1*6*fdamp5*invR5;
#ifdef COMPUTING_EXCEPTIONS
tempField1 *= scale;
tempField2 *= scale;
#endif
