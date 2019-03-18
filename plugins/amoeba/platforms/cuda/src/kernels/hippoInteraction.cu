#ifdef USE_CUTOFF
unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
unsigned int includeInteraction = (!isExcluded);
#endif

// Convert to quasi-internal coordinates.

real mat[3][3];
formQIRotationMatrix(delta, rInv, mat);
real3 labForce = make_real3(0);
real3 qiDipole1 = rotateVectorToQI(dipole1, mat);
real3 qiDipole2 = rotateVectorToQI(dipole2, mat);
real3 qiInducedDipole1 = rotateVectorToQI(inducedDipole1, mat);
real3 qiInducedDipole2 = rotateVectorToQI(inducedDipole2, mat);
real qiQXX1, qiQXY1, qiQXZ1, qiQYY1, qiQYZ1, qiQZZ1;
real qiQXX2, qiQXY2, qiQXZ2, qiQYY2, qiQYZ2, qiQZZ2;
rotateQuadrupoleToQI(qXX1, qXY1, qXZ1, qYY1, qYZ1, qiQXX1, qiQXY1, qiQXZ1, qiQYY1, qiQYZ1, qiQZZ1, mat);
rotateQuadrupoleToQI(qXX2, qXY2, qXZ2, qYY2, qYZ2, qiQXX2, qiQXY2, qiQXZ2, qiQYY2, qiQYZ2, qiQZZ2, mat);

// Compute intermediates used for dipole and quadrupole calculations.

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
real rr5 = 3*rr3*rInv2;
real rr7 = 5*rr5*rInv2;
real rr9 = 7*rr7*rInv2;
real rr11 = 9*rr9*rInv2;

// Compute the fixed multipole interaction.

{
    // Find damped multipole intermediates and energy value.

    real term1 = coreCharge1*coreCharge2;
    real term1i = coreCharge2*valenceCharge1;
    real term2i = coreCharge2*dir;
    real term3i = coreCharge2*qir;
    real term1k = coreCharge1*valenceCharge2;
    real term2k = -coreCharge1*dkr;
    real term3k = coreCharge1*qkr;
    real term1ik = valenceCharge1*valenceCharge2;
    real term2ik = valenceCharge2*dir - valenceCharge1*dkr + dik;
    real term3ik = valenceCharge1*qkr + valenceCharge2*qir - dir*dkr + 2*(dkqi-diqk+qiqk);
    real term4ik = dir*qkr - dkr*qir - 4*qik;
    real term5ik = qir*qkr;
    real fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    real fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    real fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(alpha1, alpha2, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    real rr1i = fdampI1*rr1;
    real rr3i = fdampI3*rr3;
    real rr5i = fdampI5*rr5;
    real rr7i = fdampI7*rr7;
    real rr1k = fdampK1*rr1;
    real rr3k = fdampK3*rr3;
    real rr5k = fdampK5*rr5;
    real rr7k = fdampK7*rr7;
    real rr1ik = fdampIK1*rr1;
    real rr3ik = fdampIK3*rr3;
    real rr5ik = fdampIK5*rr5;
    real rr7ik = fdampIK7*rr7;
    real rr9ik = fdampIK9*rr9;
    real rr11ik = fdampIK11*rr11;
    real scale = ENERGY_SCALE_FACTOR;
#ifdef COMPUTING_EXCEPTIONS
    scale *= multipoleMultipoleScale;
#endif
    real elecEnergy = scale*(term1*rr1 + term4ik*rr7ik + term5ik*rr9ik +
                             term1i*rr1i + term1k*rr1k + term1ik*rr1ik +
                             term2i*rr3i + term2k*rr3k + term2ik*rr3ik +
                             term3i*rr5i + term3k*rr5k + term3ik*rr5ik);

    // Find damped multipole intermediates for force and torque.

    real de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik +
              term1i*rr3i + term1k*rr3k + term1ik*rr3ik +
              term2i*rr5i + term2k*rr5k + term2ik*rr5ik +
              term3i*rr7i + term3k*rr7k + term3ik*rr7ik;
    term1 = -coreCharge2*rr3i - valenceCharge2*rr3ik + dkr*rr5ik - qkr*rr7ik;
    real term2 = coreCharge1*rr3k + valenceCharge1*rr3ik + dir*rr5ik + qir*rr7ik;
    real term3 = 2*rr5ik;
    real term4 = -2*(coreCharge2*rr5i+valenceCharge2*rr5ik-dkr*rr7ik+qkr*rr9ik);
    real term5 = -2*(coreCharge1*rr5k+valenceCharge1*rr5ik+dir*rr7ik+qir*rr9ik);
    real term6 = 4*rr7ik;

    // Compute the force and torque.

    real3 elecForce = -scale*(de*make_real3(0, 0, r) + term1*qiDipole1 + term2*qiDipole2 +
            term3*(diqkTemp-dkqiTemp) + term4*qi + term5*qk + term6*(qikTemp+qkiTemp));
    real3 tI = -scale*(-rr3ik*dikCross + term1*dirCross + term3*(dqik+dkqirCross) + term4*qirCross - term6*(qikrCross+qikCross));
    real3 tK = -scale*(rr3ik*dikCross + term2*dkrCross - term3*(dqik+diqkrCross) + term5*qkrCross - term6*(qkirCross-qikCross));
    tempEnergy += includeInteraction ? elecEnergy : 0;
    tempForce += includeInteraction ? elecForce : make_real3(0);
    tempTorque1 += includeInteraction ? tI : make_real3(0);
    tempTorque2 += includeInteraction ? tK : make_real3(0);
}

// Compute the induced dipole interactions.

{
    real uir = qiInducedDipole1.z*r;
    real ukr = qiInducedDipole2.z*r;

    // Apply charge penetration damping to scale factors.

    real fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    real fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    real fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(alpha1, alpha2, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    real scale = ENERGY_SCALE_FACTOR;
#ifdef COMPUTING_EXCEPTIONS
    scale *= dipoleMultipoleScale;
#endif
    real dsr3i = rr3*fdampI3*scale;
    real dsr5i = rr5*fdampI5*scale;
    real dsr7i = rr7*fdampI7*scale;
    real dsr3k = rr3*fdampK3*scale;
    real dsr5k = rr5*fdampK5*scale;
    real dsr7k = rr7*fdampK7*scale;

    // Get the induced dipole field used for dipole torques.

    real3 torqueField1 = dsr3i*qiInducedDipole2;
    torqueField1.z -= dsr5i*ukr*r;
    real3 torqueField2 = dsr3k*qiInducedDipole1;
    torqueField2.z -= dsr5k*uir*r;

    // Get induced dipole field gradient used for quadrupole torques.

    real3 dtorqueField1 = 2*r*dsr5i*qiInducedDipole2;
    dtorqueField1.z -= r2*dsr7i*ukr;
    real3 dtorqueField2 = -2*r*dsr5k*qiInducedDipole1;
    dtorqueField2.z += r2*dsr7k*uir;

    // Get the field gradient for direct polarization force

    real t1XX = valenceCharge1*rr3*fdampI3 + coreCharge1*rr3 + dir*rr5*fdampI5 - qxI.x*2*rr5*fdampI5 + qi.z*r*rr7*fdampI7;
    real t2XX = valenceCharge2*rr3*fdampK3 + coreCharge2*rr3 - dkr*rr5*fdampK5 - qxK.x*2*rr5*fdampK5 + qk.z*r*rr7*fdampK7;
    real t1YY = valenceCharge1*rr3*fdampI3 + coreCharge1*rr3 + dir*rr5*fdampI5 - qyI.y*2*rr5*fdampI5 + qi.z*r*rr7*fdampI7;
    real t2YY = valenceCharge2*rr3*fdampK3 + coreCharge2*rr3 - dkr*rr5*fdampK5 - qyK.y*2*rr5*fdampK5 + qk.z*r*rr7*fdampK7;
    real t1ZZ = valenceCharge1*(rr3*fdampI3-rr5*fdampI5*r2) + coreCharge1*(rr3-rr5*r2) + qiDipole1.z*2*rr5*fdampI5*r -
              dir*(rr7*fdampI7*r2-rr5*fdampI5) - qzI.z*2*rr5*fdampI5 + qi.z*5*rr7*fdampI7*r - qir*rr9*fdampI9*r2;
    real t2ZZ = valenceCharge2*(rr3*fdampK3-rr5*fdampK5*r2) + coreCharge2*(rr3-rr5*r2) - qiDipole2.z*2*rr5*fdampK5*r +
              dkr*(rr7*fdampK7*r2-rr5*fdampK5) - qzK.z*2*rr5*fdampK5 + qk.z*5*rr7*fdampK7*r - qkr*rr9*fdampK9*r2;
    real t1XY = -qxI.y*2*rr5*fdampI5;
    real t2XY = -qxK.y*2*rr5*fdampK5;
    real t1XZ = qiDipole1.x*rr5*fdampI5*r - qxI.z*2*rr5*fdampI5 + qi.x*2*rr7*fdampI7*r;
    real t2XZ = -qiDipole2.x*rr5*fdampK5*r - qxK.z*2*rr5*fdampK5 + qk.x*2*rr7*fdampK7*r;
    real t1YZ = qiDipole1.y*rr5*fdampI5*r - qyI.z*2*rr5*fdampI5 + qi.y*2*rr7*fdampI7*r;
    real t2YZ = -qiDipole2.y*rr5*fdampK5*r - qyK.z*2*rr5*fdampK5 + qk.y*2*rr7*fdampK7*r;

    // Get the dEp/dR terms for chgpen direct polarization force.

    real depx = t1XX*qiInducedDipole2.x + t1XY*qiInducedDipole2.y + t1XZ*qiInducedDipole2.z - t2XX*qiInducedDipole1.x - t2XY*qiInducedDipole1.y - t2XZ*qiInducedDipole1.z;
    real depy = t1XY*qiInducedDipole2.x + t1YY*qiInducedDipole2.y + t1YZ*qiInducedDipole2.z - t2XY*qiInducedDipole1.x - t2YY*qiInducedDipole1.y - t2YZ*qiInducedDipole1.z;
    real depz = t1XZ*qiInducedDipole2.x + t1YZ*qiInducedDipole2.y + t1ZZ*qiInducedDipole2.z - t2XZ*qiInducedDipole1.x - t2YZ*qiInducedDipole1.y - t2ZZ*qiInducedDipole1.z;
    real3 indForce = scale*make_real3(depx, depy, depz);

    // Torque is induced field and gradient cross permanent moments.

    real3 tI = cross(torqueField1, qiDipole1);
    tI.x += -qxI.y*dtorqueField1.x - 2*qyI.z*dtorqueField1.z + (qzI.z-qyI.y)*dtorqueField1.y;
    tI.y += qxI.y*dtorqueField1.y + 2*qxI.z*dtorqueField1.z + (qxI.x-qzI.z)*dtorqueField1.x;
    tI.z += qyI.z*dtorqueField1.x - qxI.z*dtorqueField1.y;
    real3 tK = cross(torqueField2, qiDipole2);
    tK.x += -qxK.y*dtorqueField2.x - 2*qyK.z*dtorqueField2.z + (qzK.z-qyK.y)*dtorqueField2.y;
    tK.y += qxK.y*dtorqueField2.y + 2*qxK.z*dtorqueField2.z + (qxK.x-qzK.z)*dtorqueField2.x;
    tK.z += qyK.z*dtorqueField2.x - qxK.z*dtorqueField2.y;
    tempForce -= includeInteraction ? indForce : make_real3(0);
    tempTorque1 -= includeInteraction ? tI : make_real3(0);
    tempTorque2 -= includeInteraction ? tK : make_real3(0);

    // Get the dtau/dr terms used for OPT polarization force.

#ifdef COMPUTING_EXCEPTIONS
    real ddscale = dipoleDipoleScale*ENERGY_SCALE_FACTOR/2;
#else
    real ddscale = ENERGY_SCALE_FACTOR/2;
#endif
    real coeff[] = {EXTRAPOLATION_COEFFICIENTS_SUM};
    for (int j = 0; j < MAX_EXTRAPOLATION_ORDER-1; j++) {
        real3 extDipole1 = (atom1 < NUM_ATOMS ? extrapolatedDipole[j*NUM_ATOMS+atom1] : make_real3(0));
        real uirm = dot(extDipole1, delta);
        for (int m = 0; m < MAX_EXTRAPOLATION_ORDER-1-j; m++) {
            real3 extDipole2 = (atom2 < NUM_ATOMS ? extrapolatedDipole[m*NUM_ATOMS+atom2] : make_real3(0));
            real ukrm = dot(extDipole2, delta);
            real term1 = 2*fdampIK5*rr5;
            real term2 = term1*delta.x;
            real term3 = rr5*fdampIK5 - rr7*fdampIK7*delta.x*delta.x;
            real tixx = extDipole1.x*term2 + uirm*term3;
            real tkxx = extDipole2.x*term2 + ukrm*term3;
            term2 = term1*delta.y;
            term3 = rr5*fdampIK5 - rr7*fdampIK7*delta.y*delta.y;
            real tiyy = extDipole1.y*term2 + uirm*term3;
            real tkyy = extDipole2.y*term2 + ukrm*term3;
            term2 = term1*delta.z;
            term3 = rr5*fdampIK5 - rr7*fdampIK7*delta.z*delta.z;
            real tizz = extDipole1.z*term2 + uirm*term3;
            real tkzz = extDipole2.z*term2 + ukrm*term3;
            term1 = rr5*fdampIK5*delta.y;
            term2 = rr5*fdampIK5*delta.x;
            term3 = delta.y * (rr7*fdampIK7*delta.x);
            real tixy = extDipole1.x*term1 + extDipole1.y*term2 - uirm*term3;
            real tkxy = extDipole2.x*term1 + extDipole2.y*term2 - ukrm*term3;
            term1 = rr5 *fdampIK5 * delta.z;
            term3 = delta.z * (rr7*fdampIK7*delta.x);
            real tixz = extDipole1.x*term1 + extDipole1.z*term2 - uirm*term3;
            real tkxz = extDipole2.x*term1 + extDipole2.z*term2 - ukrm*term3;
            term2 = rr5*fdampIK5*delta.y;
            term3 = delta.z * (rr7*fdampIK7*delta.y);
            real tiyz = extDipole1.y*term1 + extDipole1.z*term2 - uirm*term3;
            real tkyz = extDipole2.y*term1 + extDipole2.z*term2 - ukrm*term3;
            real depx = tixx*extDipole2.x + tkxx*extDipole1.x + tixy*extDipole2.y + tkxy*extDipole1.y + tixz*extDipole2.z + tkxz*extDipole1.z;
            real depy = tixy*extDipole2.x + tkxy*extDipole1.x + tiyy*extDipole2.y + tkyy*extDipole1.y + tiyz*extDipole2.z + tkyz*extDipole1.z;
            real depz = tixz*extDipole2.x + tkxz*extDipole1.x + tiyz*extDipole2.y + tkyz*extDipole1.y + tizz*extDipole2.z + tkzz*extDipole1.z;
            labForce += ddscale*coeff[j+m+1]*make_real3(depx, depy, depz);
        }
    }
}

// Compute the repulsion interaction.

{
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
#ifdef COMPUTING_EXCEPTIONS
    ctForce *= multipoleMultipoleScale;
    ctEnergy *= multipoleMultipoleScale;
#endif
    tempEnergy += includeInteraction ? ctEnergy : 0;
    tempForce.z += includeInteraction ? ctForce*r : 0;
}

// Rotate back to the lab frame.

if (includeInteraction) {
    tempForce = rotateVectorFromQI(tempForce, mat) - labForce;
    tempTorque1 = rotateVectorFromQI(tempTorque1, mat);
    tempTorque2 = rotateVectorFromQI(tempTorque2, mat);
}
