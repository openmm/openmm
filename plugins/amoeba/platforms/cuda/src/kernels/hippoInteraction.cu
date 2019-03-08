#ifdef USE_CUTOFF
unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
unsigned int includeInteraction = (!isExcluded);
#endif

// Convert to quasi-internal coordinates.

real mat[3][3];
formQIRotationMatrix(delta, rInv, mat);
real3 qiDipole1 = rotateVectorToQI(dipole1, mat);
real3 qiDipole2 = rotateVectorToQI(dipole2, mat);
real qiQXX1, qiQXY1, qiQXZ1, qiQYY1, qiQYZ1, qiQZZ1;
real qiQXX2, qiQXY2, qiQXZ2, qiQYY2, qiQYZ2, qiQZZ2;
rotateQuadrupoleToQI(qXX1, qXY1, qXZ1, qYY1, qYZ1, qiQXX1, qiQXY1, qiQXZ1, qiQYY1, qiQYZ1, qiQZZ1, mat);
rotateQuadrupoleToQI(qXX2, qXY2, qXZ2, qYY2, qYZ2, qiQXX2, qiQXY2, qiQXZ2, qiQYY2, qiQYZ2, qiQZZ2, mat);

// Compute the repulsion interaction.

{
real dir = qiDipole1.z*r;
real3 qxI = make_real3(qiQXX1, qiQXY1, qiQXZ1);
real3 qyI = make_real3(qiQXY1, qiQYY1, qiQYZ1);
real3 qzI = make_real3(qiQXZ1, qiQYZ1, qiQZZ1);
real3 qi = r*make_real3(qiQXZ1, qiQYZ1, qiQZZ1);
real qir = qi.z*r;
real dkr = qiDipole2.z*r;
real3 qxK = make_real3(qiQXX2, qiQXY2, qiQXZ2);
real3 qyK = make_real3(qiQXY2, qiQYY2, qiQYZ2);
real3 qzK = make_real3(qiQXZ2, qiQYZ2, qiQZZ2);
real3 qk = r*make_real3(qiQXZ2, qiQYZ2, qiQZZ2);
real qkr = qk.z*r;
real dik = dot(qiDipole1, qiDipole2);
real qik = dot(qi, qk);
real diqk = dot(qiDipole1, qk);
real dkqi = dot(qiDipole2, qi);
real qiqk = 2*(qxI.y*qxK.y+qxI.z*qxK.z+qyI.z*qyK.z) + qxI.x*qxK.x + qyI.y*qyK.y + qzI.z*qzK.z;

// Additional intermediates involving moments and distance.

real3 dirCross = make_real3(qiDipole1.y*r, -qiDipole1.x*r, 0);
real3 dkrCross = make_real3(qiDipole2.y*r, -qiDipole2.x*r, 0);
real3 dikCross = cross(qiDipole1, qiDipole2);
real3 qirCross = make_real3(qi.y*r, -qi.x*r, 0);
real3 qkrCross = make_real3(qk.y*r, -qk.x*r, 0);
real3 qikCross = cross(qk, qi);
real3 qikTemp = make_real3(dot(qxI, qk), dot(qyI, qk), dot(qzI, qk));
real3 qkiTemp = make_real3(dot(qxK, qi), dot(qyK, qi), dot(qzK, qi));
real3 qikrCross = make_real3(-r*qikTemp.y, r*qikTemp.x, 0);
real3 qkirCross = make_real3(-r*qkiTemp.y, r*qkiTemp.x, 0);
real3 diqkTemp = make_real3(dot(qiDipole1, qxK), dot(qiDipole1, qyK), dot(qiDipole1, qzK));
real3 dkqiTemp = make_real3(dot(qiDipole2, qxI), dot(qiDipole2, qyI), dot(qiDipole2, qzI));
real3 diqkrCross = make_real3(-r*diqkTemp.y, r*diqkTemp.x, 0);
real3 dkqirCross = make_real3(-r*dkqiTemp.y, r*dkqiTemp.x, 0);
real3 dqik = cross(qiDipole1, qk) + cross(qiDipole2, qi) - 2*(cross(qxI, qxK) + cross(qyI, qyK) + cross(qzI, qzK));

// Get reciprocal distance terms for this interaction.

real rInv2 = rInv*rInv;
real rr1 = rInv;
real rr3 = rr1*rInv2;

// Compute damping coefficients.

real fdamp1, fdamp3, fdamp5, fdamp7, fdamp9, fdamp11;
computeRepulsionDampingFactors(pauliAlpha1, pauliAlpha2, r, fdamp1, fdamp3, fdamp5, fdamp7, fdamp9, fdamp11);

// Calculate intermediate terms needed for the energy

real eterm1 = pauliQ1*pauliQ2;
real eterm2 = pauliQ2*dir - pauliQ1*dkr + dik;
real eterm3 = pauliQ1*qkr + pauliQ2*qir - dir*dkr + 2*(dkqi-diqk+qiqk);
real eterm4 = dir*qkr - dkr*qir - 4*qik;
real eterm5 = qir*qkr;
real eterm = eterm1*fdamp1 + eterm2*fdamp3 + eterm3*fdamp5 + eterm4*fdamp7 + eterm5*fdamp9;

// Compute the energy.

real sizik = pauliK1*pauliK2;
#ifdef COMPUTING_EXCEPTIONS
sizik *= repulsionScale;
#endif
real repEnergy = sizik*eterm*rr1;

// Calculate intermediate terms for force and torque

real de = eterm1*fdamp3 + eterm2*fdamp5 + eterm3*fdamp7 + eterm4*fdamp9 + eterm5*fdamp11;
real term1 = -pauliQ2*fdamp3 + dkr*fdamp5 - qkr*fdamp7;
real term2 = pauliQ1*fdamp3 + dir*fdamp5 + qir*fdamp7;
real term3 = 2*fdamp5;
real term4 = 2*(-pauliQ2*fdamp5 + dkr*fdamp7 - qkr*fdamp9);
real term5 = 2*(-pauliQ1*fdamp5 - dir*fdamp7 - qir*fdamp9);
real term6 = 4*fdamp7;

// Compute the force and torque.

real3 repForce = make_real3(0, 0, de*r) + term1*qiDipole1 + term2*qiDipole2 + term3*(diqkTemp-dkqiTemp)
        + term4*qi + term5*qk + term6*(qikTemp+qkiTemp);
repForce = -sizik*(repForce*rr1 + make_real3(0, 0, eterm*rr3*r));
real3 tI = -fdamp3*dikCross + term1*dirCross + term3*(dqik+dkqirCross) + term4*qirCross - term6*(qikrCross+qikCross);
real3 tK = fdamp3*dikCross + term2*dkrCross - term3*(dqik+diqkrCross) + term5*qkrCross - term6*(qkirCross-qikCross);
tI *= -sizik*rr1;
tK *= -sizik*rr1;
#ifdef USE_CUTOFF
if (r > SWITCH_CUTOFF) {
    real x = r-SWITCH_CUTOFF;
    real switchValue = 1+x*x*x*(SWITCH_C3+x*(SWITCH_C4+x*SWITCH_C5));
    real switchDeriv = x*x*(3*SWITCH_C3+x*(4*SWITCH_C4+x*5*SWITCH_C5));
    repForce *= switchValue;
    repForce.z += repEnergy*switchDeriv;
    repEnergy *= switchValue;
    tI *= switchValue;
    tK *= switchValue;
}
#endif
tempEnergy += includeInteraction ? repEnergy : 0;
tempForce += includeInteraction ? repForce : make_real3(0);
tempTorque1 += includeInteraction ? tI : make_real3(0);
tempTorque2 += includeInteraction ? tK : make_real3(0);
}

// Compute the dispersion force and energy.

{
real dispEnergy = -c61*c62/(r2*r2*r2);
real dispForce = -6*dispEnergy*rInv;
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
real ctForce = (term1*damping2 + term2*damping1)*rInv;
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

if (includeInteraction) {
    tempForce = rotateVectorFromQI(tempForce, mat);
    tempTorque1 = rotateVectorFromQI(tempTorque1, mat);
    tempTorque2 = rotateVectorFromQI(tempTorque2, mat);
}
