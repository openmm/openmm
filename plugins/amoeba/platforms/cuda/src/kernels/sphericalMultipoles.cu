__device__ void buildQIRotationMatrix(real3 deltaR, real rInv, real (&rotationMatrix)[3][3]) {
    real3 vectorZ = deltaR*rInv;
    real3 vectorX = vectorZ;
    if (deltaR.y != 0 || deltaR.z != 0)
        vectorX.x += 1;
    else
        vectorX.y += 1;

    vectorX -= vectorZ*dot(vectorX, vectorZ);
    vectorX = normalize(vectorX);
    real3 vectorY = cross(vectorZ, vectorX);

    // Reorder the Cartesian {x,y,z} dipole rotation matrix, to account
    // for spherical harmonic ordering {z,x,y}.
    rotationMatrix[0][0] = vectorZ.z;
    rotationMatrix[0][1] = vectorZ.x;
    rotationMatrix[0][2] = vectorZ.y;
    rotationMatrix[1][0] = vectorX.z;
    rotationMatrix[1][1] = vectorX.x;
    rotationMatrix[1][2] = vectorX.y;
    rotationMatrix[2][0] = vectorY.z;
    rotationMatrix[2][1] = vectorY.x;
    rotationMatrix[2][2] = vectorY.y;
}

__device__ real3 rotateDipole(real3& dipole, const real (&rotationMatrix)[3][3]) {
    return make_real3(rotationMatrix[0][0]*dipole.x + rotationMatrix[0][1]*dipole.y + rotationMatrix[0][2]*dipole.z,
                      rotationMatrix[1][0]*dipole.x + rotationMatrix[1][1]*dipole.y + rotationMatrix[1][2]*dipole.z,
                      rotationMatrix[2][0]*dipole.x + rotationMatrix[2][1]*dipole.y + rotationMatrix[2][2]*dipole.z);
}


__device__ void rotateQuadupoles(const real (&rotationMatrix)[3][3], const real* quad1, const real* quad2, real* rotated1, real* rotated2) {
    real sqrtThree = SQRT((real) 3);
    real element;
    element = 0.5f*(3.0f*rotationMatrix[0][0]*rotationMatrix[0][0] - 1.0f);
    rotated1[0] += quad1[0]*element;
    rotated2[0] += quad2[0]*element;
    element = sqrtThree*rotationMatrix[0][0]*rotationMatrix[0][1];
    rotated1[0] += quad1[1]*element;
    rotated2[0] += quad2[1]*element;
    element = sqrtThree*rotationMatrix[0][0]*rotationMatrix[0][2];
    rotated1[0] += quad1[2]*element;
    rotated2[0] += quad2[2]*element;
    element = 0.5f*sqrtThree*(rotationMatrix[0][1]*rotationMatrix[0][1] - rotationMatrix[0][2]*rotationMatrix[0][2]);
    rotated1[0] += quad1[3]*element;
    rotated2[0] += quad2[3]*element;
    element = sqrtThree*rotationMatrix[0][1]*rotationMatrix[0][2];
    rotated1[0] += quad1[4]*element;
    rotated2[0] += quad2[4]*element;
    element = sqrtThree*rotationMatrix[0][0]*rotationMatrix[1][0];
    rotated1[1] += quad1[0]*element;
    rotated2[1] += quad2[0]*element;
    element = rotationMatrix[1][0]*rotationMatrix[0][1] + rotationMatrix[0][0]*rotationMatrix[1][1];
    rotated1[1] += quad1[1]*element;
    rotated2[1] += quad2[1]*element;
    element = rotationMatrix[1][0]*rotationMatrix[0][2] + rotationMatrix[0][0]*rotationMatrix[1][2];
    rotated1[1] += quad1[2]*element;
    rotated2[1] += quad2[2]*element;
    element = rotationMatrix[0][1]*rotationMatrix[1][1] - rotationMatrix[0][2]*rotationMatrix[1][2];
    rotated1[1] += quad1[3]*element;
    rotated2[1] += quad2[3]*element;
    element = rotationMatrix[1][1]*rotationMatrix[0][2] + rotationMatrix[0][1]*rotationMatrix[1][2];
    rotated1[1] += quad1[4]*element;
    rotated2[1] += quad2[4]*element;
    element = sqrtThree*rotationMatrix[0][0]*rotationMatrix[2][0];
    rotated1[2] += quad1[0]*element;
    rotated2[2] += quad2[0]*element;
    element = rotationMatrix[2][0]*rotationMatrix[0][1] + rotationMatrix[0][0]*rotationMatrix[2][1];
    rotated1[2] += quad1[1]*element;
    rotated2[2] += quad2[1]*element;
    element = rotationMatrix[2][0]*rotationMatrix[0][2] + rotationMatrix[0][0]*rotationMatrix[2][2];
    rotated1[2] += quad1[2]*element;
    rotated2[2] += quad2[2]*element;
    element = rotationMatrix[0][1]*rotationMatrix[2][1] - rotationMatrix[0][2]*rotationMatrix[2][2];
    rotated1[2] += quad1[3]*element;
    rotated2[2] += quad2[3]*element;
    element = rotationMatrix[2][1]*rotationMatrix[0][2] + rotationMatrix[0][1]*rotationMatrix[2][2];
    rotated1[2] += quad1[4]*element;
    rotated2[2] += quad2[4]*element;
    element = 0.5f*sqrtThree*(rotationMatrix[1][0]*rotationMatrix[1][0] - rotationMatrix[2][0]*rotationMatrix[2][0]);
    rotated1[3] += quad1[0]*element;
    rotated2[3] += quad2[0]*element;
    element = rotationMatrix[1][0]*rotationMatrix[1][1] - rotationMatrix[2][0]*rotationMatrix[2][1];
    rotated1[3] += quad1[1]*element;
    rotated2[3] += quad2[1]*element;
    element = rotationMatrix[1][0]*rotationMatrix[1][2] - rotationMatrix[2][0]*rotationMatrix[2][2];
    rotated1[3] += quad1[2]*element;
    rotated2[3] += quad2[2]*element;
    element = 0.5f*(rotationMatrix[1][1]*rotationMatrix[1][1] - rotationMatrix[2][1]*rotationMatrix[2][1] - rotationMatrix[1][2]*rotationMatrix[1][2] + rotationMatrix[2][2]*rotationMatrix[2][2]);
    rotated1[3] += quad1[3]*element;
    rotated2[3] += quad2[3]*element;
    element = rotationMatrix[1][1]*rotationMatrix[1][2] - rotationMatrix[2][1]*rotationMatrix[2][2];
    rotated1[3] += quad1[4]*element;
    rotated2[3] += quad2[4]*element;
    element = sqrtThree*rotationMatrix[1][0]*rotationMatrix[2][0];
    rotated1[4] += quad1[0]*element;
    rotated2[4] += quad2[0]*element;
    element = rotationMatrix[2][0]*rotationMatrix[1][1] + rotationMatrix[1][0]*rotationMatrix[2][1];
    rotated1[4] += quad1[1]*element;
    rotated2[4] += quad2[1]*element;
    element = rotationMatrix[2][0]*rotationMatrix[1][2] + rotationMatrix[1][0]*rotationMatrix[2][2];
    rotated1[4] += quad1[2]*element;
    rotated2[4] += quad2[2]*element;
    element = rotationMatrix[1][1]*rotationMatrix[2][1] - rotationMatrix[1][2]*rotationMatrix[2][2];
    rotated1[4] += quad1[3]*element;
    rotated2[4] += quad2[3]*element;
    element = rotationMatrix[2][1]*rotationMatrix[1][2] + rotationMatrix[1][1]*rotationMatrix[2][2];
    rotated1[4] += quad1[4]*element;
    rotated2[4] += quad2[4]*element;
}
