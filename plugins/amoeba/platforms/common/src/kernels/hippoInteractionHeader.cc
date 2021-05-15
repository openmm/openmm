// Functions that are called from hippoInteraction.cu.

__device__ void formQIRotationMatrix(real3 deltaR, real rInv, real (&rotationMatrix)[3][3]) {
    real3 vectorZ = deltaR*rInv;
    real3 vectorX = make_real3(0);
    if (fabs(vectorZ.y) > fabs(vectorZ.x))
        vectorX.x = 1;
    else
        vectorX.y = 1;

    vectorX -= vectorZ*dot(vectorZ, vectorX);
    vectorX = normalize(vectorX);
    real3 vectorY = cross(vectorZ, vectorX);

    rotationMatrix[0][0] = vectorX.x;
    rotationMatrix[0][1] = vectorX.y;
    rotationMatrix[0][2] = vectorX.z;
    rotationMatrix[1][0] = vectorY.x;
    rotationMatrix[1][1] = vectorY.y;
    rotationMatrix[1][2] = vectorY.z;
    rotationMatrix[2][0] = vectorZ.x;
    rotationMatrix[2][1] = vectorZ.y;
    rotationMatrix[2][2] = vectorZ.z;
}

__device__ real3 rotateVectorToQI(real3 v, const real (&mat)[3][3]) {
    return make_real3(mat[0][0]*v.x + mat[0][1]*v.y + mat[0][2]*v.z,
                      mat[1][0]*v.x + mat[1][1]*v.y + mat[1][2]*v.z,
                      mat[2][0]*v.x + mat[2][1]*v.y + mat[2][2]*v.z);
}

__device__ real3 rotateVectorFromQI(real3 v, const real (&mat)[3][3]) {
    return make_real3(mat[0][0]*v.x + mat[1][0]*v.y + mat[2][0]*v.z,
                      mat[0][1]*v.x + mat[1][1]*v.y + mat[2][1]*v.z,
                      mat[0][2]*v.x + mat[1][2]*v.y + mat[2][2]*v.z);
}

__device__ void rotateQuadrupoleToQI(real qXX, real qXY, real qXZ, real qYY, real qYZ, 
            real &qiQXX, real &qiQXY, real &qiQXZ, real &qiQYY, real &qiQYZ, real &qiQZZ, const real (&mat)[3][3]) {
    real qZZ = -qXX-qYY;
    qiQXX = mat[0][0]*(mat[0][0]*qXX + 2*(mat[0][1]*qXY + mat[0][2]*qXZ)) + mat[0][1]*(mat[0][1]*qYY + 2*mat[0][2]*qYZ) + mat[0][2]*mat[0][2]*qZZ;
    qiQYY = mat[1][0]*(mat[1][0]*qXX + 2*(mat[1][1]*qXY + mat[1][2]*qXZ)) + mat[1][1]*(mat[1][1]*qYY + 2*mat[1][2]*qYZ) + mat[1][2]*mat[1][2]*qZZ;
    qiQXY = mat[0][0]*mat[1][0]*qXX + mat[0][1]*mat[1][1]*qYY + mat[0][2]*mat[1][2]*qZZ + (mat[0][0]*mat[1][1] + mat[0][1]*mat[1][0])*qXY + (mat[0][0]*mat[1][2] + mat[0][2]*mat[1][0])*qXZ + (mat[0][1]*mat[1][2] + mat[0][2]*mat[1][1])*qYZ;
    qiQXZ = mat[0][0]*mat[2][0]*qXX + mat[0][1]*mat[2][1]*qYY + mat[0][2]*mat[2][2]*qZZ + (mat[0][0]*mat[2][1] + mat[0][1]*mat[2][0])*qXY + (mat[0][0]*mat[2][2] + mat[0][2]*mat[2][0])*qXZ + (mat[0][1]*mat[2][2] + mat[0][2]*mat[2][1])*qYZ;
    qiQYZ = mat[1][0]*mat[2][0]*qXX + mat[1][1]*mat[2][1]*qYY + mat[1][2]*mat[2][2]*qZZ + (mat[1][0]*mat[2][1] + mat[1][1]*mat[2][0])*qXY + (mat[1][0]*mat[2][2] + mat[1][2]*mat[2][0])*qXZ + (mat[1][1]*mat[2][2] + mat[1][2]*mat[2][1])*qYZ;
    qiQZZ = -qiQXX-qiQYY;
}

__device__ void computeOverlapDampingFactors(real alphaI, real alphaJ, real r,
                                  real& fdampI1, real& fdampI3, real& fdampI5, real& fdampI7, real& fdampI9,
                                  real& fdampJ1, real& fdampJ3, real& fdampJ5, real& fdampJ7, real& fdampJ9,
                                  real& fdampIJ1, real& fdampIJ3, real& fdampIJ5, real& fdampIJ7, real& fdampIJ9, real& fdampIJ11) {
    real arI = alphaI*r;
    real arI2 = arI*arI;
    real arI3 = arI2*arI;
    real arI4 = arI2*arI2;
    real arI5 = arI3*arI2;
    real arI6 = arI3*arI3;
    real expARI = EXP(-arI);
    real one = 1;
    real two = 2;
    real three = 3;
    real four = 4;
    real five = 5;
    real seven = 7;
    real eleven = 11;
    fdampI1 = 1 - (1 + arI*(one/2))*expARI;
    fdampI3 = 1 - (1 + arI + arI2*(one/2))*expARI;
    fdampI5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6))*expARI;
    fdampI7 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/30))*expARI;
    fdampI9 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(four/105) + arI5*(one/210))*expARI;
    if (alphaI == alphaJ) {
        fdampJ1 = fdampI1;
        fdampJ3 = fdampI3;
        fdampJ5 = fdampI5;
        fdampJ7 = fdampI7;
        fdampJ9 = fdampI9;
        real arI7 = arI4*arI3;
        real arI8 = arI4*arI4;
        fdampIJ1 = 1 - (1 + arI*(eleven/16) + arI2*(three/16) + arI3*(one/48))*expARI;
        fdampIJ3 = 1 - (1 + arI + arI2*(one/2) + arI3*(seven/48) + arI4*(one/48))*expARI;
        fdampIJ5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/144))*expARI;
        fdampIJ7 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/120) + arI6*(one/720))*expARI;
        fdampIJ9 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/120) + arI6*(one/720) + arI7*(one/5040))*expARI;
        fdampIJ11 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/120) + arI6*(one/720) + arI7*(one/5040) + arI8*(one/45360))*expARI;
    }
    else {
        real arJ = alphaJ*r;
        real arJ2 = arJ*arJ;
        real arJ3 = arJ2*arJ;
        real arJ4 = arJ2*arJ2;
        real arJ5 = arJ3*arJ2;
        real arJ6 = arJ3*arJ3;
        real expARJ = EXP(-arJ);
        real aI2 = alphaI*alphaI;
        real aJ2 = alphaJ*alphaJ;
        real A = aJ2/(aJ2-aI2);
        real B = aI2/(aI2-aJ2);
        real A2expARI = A*A*expARI;
        real B2expARJ = B*B*expARJ;
        fdampJ1 = 1 - (1 + arJ*(one/2))*expARJ;
        fdampJ3 = 1 - (1 + arJ + arJ2*(one/2))*expARJ;
        fdampJ5 = 1 - (1 + arJ + arJ2*(one/2) + arJ3*(one/6))*expARJ;
        fdampJ7 = 1 - (1 + arJ + arJ2*(one/2) + arJ3*(one/6) + arJ4*(one/30))*expARJ;
        fdampJ9 = 1 - (1 + arJ + arJ2*(one/2) + arJ3*(one/6) + 4*arJ4*(one/105) + arJ5*(one/210))*expARJ;
        fdampIJ1 = 1 - (1 + 2*B + arI*(one/2))*A2expARI -
                       (1 + 2*A + arJ*(one/2))*B2expARJ;
        fdampIJ3 = 1 - (1 + arI + arI2*(one/2))*A2expARI -
                       (1 + arJ + arJ2*(one/2))*B2expARJ -
                       (1 + arI)*2*B*A2expARI -
                       (1 + arJ)*2*A*B2expARJ;
        fdampIJ5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6))*A2expARI -
                       (1 + arJ + arJ2*(one/2) + arJ3*(one/6))*B2expARJ -
                       (1 + arI + arI2*(one/3))*2*B*A2expARI -
                       (1 + arJ + arJ2*(one/3))*2*A*B2expARJ;
        fdampIJ7 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/30))*A2expARI -
                       (1 + arJ + arJ2*(one/2) + arJ3*(one/6) + arJ4*(one/30))*B2expARJ -
                       (1 + arI + arI2*(two/5) + arI3*(one/15))*2*B*A2expARI -
                       (1 + arJ + arJ2*(two/5) + arJ3*(one/15))*2*A*B2expARJ;
        fdampIJ9 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*4*(one/105) + arI5*(one/210))*A2expARI -
                       (1 + arJ + arJ2*(one/2) + arJ3*(one/6) + arJ4*4*(one/105) + arJ5*(one/210))*B2expARJ -
                       (1 + arI + arI2*(three/7) + arI3*(two/21) + arI4*(one/105))*2*B*A2expARI -
                       (1 + arJ + arJ2*(three/7) + arJ3*(two/21) + arJ4*(one/105))*2*A*B2expARJ;
        fdampIJ11 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(five/126) + arI5*(two/315) + arI6*(one/1890))*A2expARI -
                        (1 + arJ + arJ2*(one/2) + arJ3*(one/6) + arJ4*(five/126) + arJ5*(two/315) + arJ6*(one/1890))*B2expARJ -
                        (1 + arI + arI2*(four/9) + arI3*(one/9) + arI4*(one/63) + arI5*(one/945))*2*B*A2expARI -
                        (1 + arJ + arJ2*(four/9) + arJ3*(one/9) + arJ4*(one/63) + arJ5*(one/945))*2*A*B2expARJ;
    }
}

__device__ void computeDispersionDampingFactors(real alphaI, real alphaJ, real r, real& fdamp, real& ddamp) {
    real arI = alphaI*r;
    real arI2 = arI*arI;
    real arI3 = arI2*arI;
    real expARI = EXP(-arI);
    real fdamp3, fdamp5;
    real one = 1;
    real seven = 7;
    if (alphaI == alphaJ) {
        real arI4 = arI3*arI;
        real arI5 = arI4*arI;
        fdamp3 = 1 - (1 + arI + arI2*(one/2) + arI3*(seven/48) + arI4*(one/48))*expARI;
        fdamp5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6) + arI4*(one/24) + arI5*(one/144))*expARI;
        ddamp = alphaI*(arI5 - 3*arI3 - 3*arI2)*expARI*(one/96);
    }
    else {
        real arJ = alphaJ*r;
        real arJ2 = arJ*arJ;
        real arJ3 = arJ2*arJ;
        real expARJ = EXP(-arJ);
        real aI2 = alphaI*alphaI;
        real aJ2 = alphaJ*alphaJ;
        real A = aJ2/(aJ2-aI2);
        real B = aI2/(aI2-aJ2);
        real A2expARI = A*A*expARI;
        real B2expARJ = B*B*expARJ;
        fdamp3 = 1 - (1 + arI + arI2*(one/2))*A2expARI -
                     (1 + arJ + arJ2*(one/2))*B2expARJ -
                     (1 + arI)*2*B*A2expARI -
                     (1 + arJ)*2*A*B2expARJ;
        fdamp5 = 1 - (1 + arI + arI2*(one/2) + arI3*(one/6))*A2expARI -
                     (1 + arJ + arJ2*(one/2) + arJ3*(one/6))*B2expARJ -
                     (1 + arI + arI2*(one/3))*2*B*A2expARI -
                     (1 + arJ + arJ2*(one/3))*2*A*B2expARJ;
        ddamp = (arI2*alphaI*A2expARI*(r*alphaI + 4*B - 1) +
                (arJ2*alphaJ*B2expARJ*(r*alphaJ + 4*A - 1)))*(one/4);
    }
    fdamp = 1.5f*fdamp5 - 0.5f*fdamp3;
}

__device__ void computeRepulsionDampingFactors(real pauliAlphaI, real pauliAlphaJ, real r,
            real& fdamp1, real& fdamp3, real& fdamp5, real& fdamp7, real& fdamp9, real& fdamp11) {
    real r2 = r*r;
    real r3 = r2*r;
    real r4 = r2*r2;
    real r5 = r3*r2;
    real r6 = r3*r3;
    real aI2 = 0.5f*pauliAlphaI;
    real arI2 = aI2*r;
    real expI = EXP(-arI2);
    real aI2_2 = aI2*aI2;
    real aI2_3 = aI2_2*aI2;
    real aI2_4 = aI2_2*aI2_2;
    real aI2_5 = aI2_3*aI2_2;
    real aI2_6 = aI2_3*aI2_3;
    real fexp, fexp1, fexp2, fexp3, fexp4, fexp5, pre;
    real one = 1;
    real two = 2;
    real four = 4;
    real eight = 8;
    real twelve = 12;
    real sixteen = 16;
    if (pauliAlphaI == pauliAlphaJ) {
        real r7 = r4*r3;
        real r8 = r4*r4;
        real aI2_7 = aI2_4*aI2_3;
        pre = 128;
        fexp = (r + aI2*r2 + aI2_2*r3*(one/3))*expI;
        fexp1 = (aI2_2*r3 + aI2_3*r4)*expI*(one/3);
        fexp2 = aI2_4*expI*r5*(one/9);
        fexp3 = aI2_5*expI*r6*(one/45);
        fexp4 = (aI2_5*r6 + aI2_6*r7)*expI*(one/315);
        fexp5 = (aI2_5*r6 + aI2_6*r7 + aI2_7*r8*(one/3))*expI*(one/945);
    }
    else {
        real aJ2 = 0.5f*pauliAlphaJ;
        real arJ2 = aJ2*r;
        real expJ = EXP(-arJ2);
        real aJ2_2 = aJ2*aJ2;
        real aJ2_3 = aJ2_2*aJ2;
        real aJ2_4 = aJ2_2*aJ2_2;
        real aJ2_5 = aJ2_3*aJ2_2;
        real scale = 1/(aI2_2-aJ2_2);
        real aI2aJ2expI = aI2*aJ2*expI;
        real aI2aJ2expJ = aI2*aJ2*expJ;
        pre = 8192*aI2_3*aJ2_3*(scale*scale*scale*scale);
        real tmp = 4*aI2*aJ2*scale;
        fexp = (arI2-tmp)*expJ + (arJ2+tmp)*expI;
        fexp1 = (r2 - (4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (r2 + (4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp2 = (r2*(one/3) + aJ2*r3*(one/3) - ((four/3)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (r2*(one/3) + aI2*r3*(one/3) + ((four/3)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp3 = (aJ2_2*r4*(one/15) + aJ2*r3*(one/5) + r2*(one/5) - ((four/15)*aJ2_3*r3 + (eight/5)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_2*r4*(one/15) + aI2*r3*(one/5) + r2*(one/5) + ((four/15)*aI2_3*r3 + (eight/5)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp4 = (aJ2_3*r5*(one/105) + (two/35)*aJ2_2*r4 + aJ2*r3*(one/7) + r2*(one/7) - ((four/105)*aJ2_4*r4 + (eight/21)*aJ2_3*r3 + (twelve/7)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_3*r5*(one/105) + (two/35)*aI2_2*r4 + aI2*r3*(one/7) + r2*(one/7) + ((four/105)*aI2_4*r4 + (eight/21)*aI2_3*r3 + (twelve/7)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp5 = (aJ2_4*r6*(one/945) + (two/189)*aJ2_3*r5 + aJ2_2*r4*(one/21) + aJ2*r3*(one/9) + r2*(one/9) - ((four/945)*aJ2_5*r5 + (four/63)*aJ2_4*r4 + (four/9)*aJ2_3*r3 + (sixteen/9)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_4*r6*(one/945) + (two/189)*aI2_3*r5 + aI2_2*r4*(one/21) + aI2*r3*(one/9) + r2*(one/9) + ((four/945)*aI2_5*r5 + (four/63)*aI2_4*r4 + (four/9)*aI2_3*r3 + (sixteen/9)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
    }
    fexp = fexp/r;
    fexp1 = fexp1/r3;
    fexp2 = 3*fexp2/r5;
    fexp3 = 15*fexp3/(r5*r2);
    fexp4 = 105*fexp4/(r5*r4);
    fexp5 = 945*fexp5/(r5*r6);
    fdamp1 = 0.5f*pre*fexp*fexp;
    fdamp3 = pre*fexp*fexp1;
    fdamp5 = pre*(fexp*fexp2 + fexp1*fexp1);
    fdamp7 = pre*(fexp*fexp3 + 3*fexp1*fexp2);
    fdamp9 = pre*(fexp*fexp4 + 4*fexp1*fexp3 + 3*fexp2*fexp2);
    fdamp11 = pre*(fexp*fexp5 + 5*fexp1*fexp4 + 10*fexp2*fexp3);
}

