#ifdef USE_CUTOFF
unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
unsigned int includeInteraction = (!isExcluded);
#endif

// Convert to quasi-internal coordinates.

real mat[3][3];
formQIRotationMatrix(delta, invR, mat);

// Compute the dispersion force and energy.

{
real dispEnergy = -c61*c62/(r2*r2*r2);
real dispForce = -6*dispEnergy*invR;
#ifdef COMPUTING_EXCEPTIONS
dispEnergy *= dispersionScale;
dispForce *= dispersionScale;
#endif
real fdamp, ddamp;
computeDispersionDampingFactors(alpha1, alpha2, r, fdamp, ddamp);
dispForce = dispForce*fdamp*fdamp + 2*dispEnergy*fdamp*ddamp;
dispEnergy *= fdamp*fdamp;
tempEnergy += includeInteraction ? dispEnergy : 0;
tempForce.z += includeInteraction ? dispForce : 0;
}

// Compute the charge transfer force and energy.

{
real term1 = epsilon1*EXP(-damping2*r);
real term2 = epsilon2*EXP(-damping1*r);
real ctEnergy = -(term1+term2);
real ctForce = (term1*damping2 + term2*damping1)*invR;
#ifdef USE_CUTOFF
if (r > SWITCH_CUTOFF) {
    real x = r-SWITCH_CUTOFF;
    real switchValue = 1+x*x*x*(SWITCH_C3+x*(SWITCH_C4+x*SWITCH_C5));
    real switchDeriv = x*x*(3*SWITCH_C3+x*(4*SWITCH_C4+x*5*SWITCH_C5));
    ctForce = ctForce*switchValue - ctEnergy*switchDeriv*r;
    ctEnergy *= switchValue;
}
#endif
tempEnergy += includeInteraction ? ctEnergy : 0;
tempForce.z += includeInteraction ? ctForce*r : 0;
}

// Rotate back to the lab frame.

if (includeInteraction)
    tempForce = rotateVectorFromQI(tempForce, mat);
