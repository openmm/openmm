real fdamp3, fdamp5;
computeMutualFieldDampingFactors(alpha1, alpha2, r, fdamp3, fdamp5);
#ifdef COMPUTING_EXCEPTIONS
fdamp3 *= dipoleDipoleScale;
fdamp5 *= dipoleDipoleScale;
#endif
real invR2 = invR*invR;
real invR3 = invR*invR2;
real scale3 = -fdamp3*invR3;
real scale5 = 3*fdamp5*invR3*invR2;
tempField1 = inducedDipole2*scale3 + delta*scale5*dot(inducedDipole2, delta);
tempField2 = inducedDipole1*scale3 + delta*scale5*dot(inducedDipole1, delta);
