
/* Portions copyright (c) 2006-2019 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "AmoebaReferenceHippoNonbondedForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "ReferencePME.h"
#include "jama_svd.h"
#include <algorithm>

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using std::vector;
using namespace OpenMM;

AmoebaReferenceHippoNonbondedForce::AmoebaReferenceHippoNonbondedForce(const HippoNonbondedForce& force) : _electric(138.9354558456) {
    _numParticles = force.getNumParticles();
    particleData.resize(_numParticles);
    std::vector<double> dipoles, quadrupoles;
    for (int i = 0; i < _numParticles; i++) {
        particleData[i].index = i;

        MultipoleParticleData& p = particleData[i];
        double charge;
        force.getParticleParameters(i, charge, dipoles, quadrupoles, p.coreCharge, p.alpha, p.epsilon, p.damping, p.c6,
                p.pauliK, p.pauliQ, p.pauliAlpha, p.polarizability, p.axisType, p.multipoleAtomZ, p.multipoleAtomX, p.multipoleAtomY);
        p.valenceCharge = charge-p.coreCharge;

        particleData[i].localDipole[0]            = dipoles[0];
        particleData[i].localDipole[1]            = dipoles[1];
        particleData[i].localDipole[2]            = dipoles[2];

        particleData[i].localQuadrupole[QXX]      = quadrupoles[0];
        particleData[i].localQuadrupole[QXY]      = quadrupoles[1];
        particleData[i].localQuadrupole[QXZ]      = quadrupoles[2];
        particleData[i].localQuadrupole[QYY]      = quadrupoles[4];
        particleData[i].localQuadrupole[QYZ]      = quadrupoles[5];
        particleData[i].localQuadrupole[QZZ]      = quadrupoles[8];
    }
    for (int i = 0; i < force.getNumExceptions(); i++) {
        Exception e;
        force.getExceptionParameters(i, e.particle1, e.particle2, e.multipoleMultipoleScale, e.dipoleMultipoleScale,
                                     e.dipoleDipoleScale, e.dispersionScale, e.repulsionScale, e.chargeTransferScale);
        exceptions[make_pair(e.particle1, e.particle2)] = e;
        exceptions[make_pair(e.particle2, e.particle1)] = e;
    }

    setExtrapolationCoefficients(force.getExtrapolationCoefficients());
    _nonbondedMethod = force.getNonbondedMethod();
    _cutoffDistance = force.getCutoffDistance();
    _cutoffDistanceSquared = _cutoffDistance*_cutoffDistance;
    _switchingDistance = force.getSwitchingDistance();
    checkChiral();
}

HippoNonbondedForce::NonbondedMethod AmoebaReferenceHippoNonbondedForce::getNonbondedMethod() const {
    return _nonbondedMethod;
}

void AmoebaReferenceHippoNonbondedForce::setExtrapolationCoefficients(const std::vector<double> &coefficients) {
    _maxPTOrder = coefficients.size(); // This accounts for the zero-based counting; actual highest order is 1 less
    _extrapolationCoefficients = coefficients;
    _extPartCoefficients.resize(_maxPTOrder);
    for (int i = 0; i < _maxPTOrder; ++i) {
        _extPartCoefficients[i] = 0.0;
        for (int j = i; j < _maxPTOrder; ++j)
            _extPartCoefficients[i] += _extrapolationCoefficients[j];
    }
}

double AmoebaReferenceHippoNonbondedForce::normalizeVec3(Vec3& vectorToNormalize) const {
    double norm = sqrt(vectorToNormalize.dot(vectorToNormalize));
    if (norm > 0.0)
        vectorToNormalize *= (1.0/norm);
    return norm;
}

void AmoebaReferenceHippoNonbondedForce::initializeVec3Vector(vector<Vec3>& vectorToInitialize) const {
    vectorToInitialize.resize(_numParticles);
    Vec3 zeroVec(0.0, 0.0, 0.0);
    std::fill(vectorToInitialize.begin(), vectorToInitialize.end(), zeroVec);
}

void AmoebaReferenceHippoNonbondedForce::loadParticleData(const vector<Vec3>& particlePositions) {
    for (int i = 0; i < _numParticles; i++)
        particleData[i].position = particlePositions[i];
}

void AmoebaReferenceHippoNonbondedForce::checkChiralCenterAtParticle(MultipoleParticleData& particleI, int axisType,
                                                                MultipoleParticleData& particleZ, MultipoleParticleData& particleX,
                                                                MultipoleParticleData& particleY) {

    if (axisType != HippoNonbondedForce::ZThenX || particleY.index == -1)
        return;

    Vec3 deltaAD   = particleI.position - particleY.position;
    Vec3 deltaBD   = particleZ.position - particleY.position;
    Vec3 deltaCD   = particleX.position - particleY.position;

    Vec3 deltaC    = deltaBD.cross(deltaCD);
    double volume = deltaC.dot(deltaAD);

    if (volume < 0.0) {
        particleI.dipole[1] *= -1.0; // pole(3,i)
        particleI.quadrupole[QXY] *= -1.0; // pole(6,i) && pole(8,i)
        particleI.quadrupole[QYZ] *= -1.0; // pole(10,i) && pole(12,i)
    }
}

void AmoebaReferenceHippoNonbondedForce::checkChiral() {
    for (int ii = 0; ii < _numParticles; ii++) {
        MultipoleParticleData& p = particleData[ii];
        if (p.multipoleAtomY > -1) {
            checkChiralCenterAtParticle(p, p.axisType,
                                        particleData[p.multipoleAtomZ],
                                        particleData[p.multipoleAtomX],
                                        particleData[p.multipoleAtomY]);
        }
    }
}

void AmoebaReferenceHippoNonbondedForce::applyRotationMatrixToParticle( MultipoleParticleData& particleI,
                                                                  const MultipoleParticleData* particleZ,
                                                                  const MultipoleParticleData* particleX,
                                                                        MultipoleParticleData* particleY,
                                                                        int axisType) const {

    // handle case where rotation matrix is identity (e.g. single ion)

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom


    Vec3 vectorX, vectorY;
    Vec3 vectorZ = particleZ->position - particleI.position;
    normalizeVec3(vectorZ);

    // branch based on axis type

    if (axisType == HippoNonbondedForce::ZOnly) {

        // z-only

        if (fabs(vectorZ[0]) < 0.866)
            vectorX = Vec3(1.0, 0.0, 0.0);
        else
            vectorX = Vec3(0.0, 1.0, 0.0);
    }
    else {
        vectorX = particleX->position - particleI.position;
        if (axisType == HippoNonbondedForce::Bisector) {

            // bisector

            // dx = dx1 + dx2 (in TINKER code)

            normalizeVec3(vectorX);
            vectorZ += vectorX;
            normalizeVec3(vectorZ);
        }
        else if (axisType == HippoNonbondedForce::ZBisect) {

            // z-bisect

            // dx = dx1 + dx2 (in TINKER code)

            normalizeVec3(vectorX);

            vectorY  = particleY->position - particleI.position;
            normalizeVec3(vectorY);

            vectorX += vectorY;
            normalizeVec3(vectorX);
        }
        else if (axisType == HippoNonbondedForce::ThreeFold) {

            // 3-fold

            // dx = dx1 + dx2 + dx3 (in TINKER code)

            normalizeVec3(vectorX);

            vectorY   = particleY->position - particleI.position;
            normalizeVec3(vectorY);

            vectorZ  += vectorX +  vectorY;
            normalizeVec3(vectorZ);
        }
    }

    double dot = vectorZ.dot(vectorX);
    vectorX -= vectorZ*dot;

    normalizeVec3(vectorX);
    vectorY = vectorZ.cross(vectorX);

    Vec3 rotationMatrix[3];
    rotationMatrix[0] = vectorX;
    rotationMatrix[1] = vectorY;
    rotationMatrix[2] = vectorZ;

    Vec3 labDipole;
    for (int ii = 0; ii < 3; ii++) {
        labDipole[ii] = particleI.localDipole[0]*rotationMatrix[0][ii];
        for (int jj = 1; jj < 3; jj++) {
            labDipole[ii] += particleI.localDipole[jj]*rotationMatrix[jj][ii];
        }
    }
    particleI.dipole = labDipole;

    double mPole[3][3];
    double rPole[3][3] = { { 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0 } };

    mPole[0][0] = particleI.localQuadrupole[QXX];
    mPole[0][1] = particleI.localQuadrupole[QXY];
    mPole[0][2] = particleI.localQuadrupole[QXZ];

    mPole[1][0] = particleI.localQuadrupole[QXY];
    mPole[1][1] = particleI.localQuadrupole[QYY];
    mPole[1][2] = particleI.localQuadrupole[QYZ];

    mPole[2][0] = particleI.localQuadrupole[QXZ];
    mPole[2][1] = particleI.localQuadrupole[QYZ];
    mPole[2][2] = particleI.localQuadrupole[QZZ];

    for (int ii = 0; ii < 3; ii++) {
       for (int jj = ii; jj < 3; jj++) {
          for (int kk = 0; kk < 3; kk++) {
             for (int mm = 0; mm < 3; mm++) {
                 rPole[ii][jj] += rotationMatrix[kk][ii]*rotationMatrix[mm][jj]*mPole[kk][mm];
             }
          }
       }
    }

    particleI.quadrupole[QXX] = rPole[0][0];
    particleI.quadrupole[QXY] = rPole[0][1];
    particleI.quadrupole[QXZ] = rPole[0][2];

    particleI.quadrupole[QYY] = rPole[1][1];
    particleI.quadrupole[QYZ] = rPole[1][2];
    particleI.quadrupole[QZZ] = rPole[2][2];
}

void AmoebaReferenceHippoNonbondedForce::formQIRotationMatrix(const Vec3 &deltaR, double r, double (&rotationMatrix)[3][3]) const {
    Vec3 vectorZ = deltaR/r;
    Vec3 vectorX;
    if (fabs(vectorZ[1]) > fabs(vectorZ[0]))
        vectorX[0] = 1;
    else
        vectorX[1] = 1;

    double dot = vectorZ.dot(vectorX);
    vectorX -= vectorZ*dot;
    normalizeVec3(vectorX);
    Vec3 vectorY = vectorZ.cross(vectorX);

    rotationMatrix[0][0] = vectorX[0];
    rotationMatrix[0][1] = vectorX[1];
    rotationMatrix[0][2] = vectorX[2];
    rotationMatrix[1][0] = vectorY[0];
    rotationMatrix[1][1] = vectorY[1];
    rotationMatrix[1][2] = vectorY[2];
    rotationMatrix[2][0] = vectorZ[0];
    rotationMatrix[2][1] = vectorZ[1];
    rotationMatrix[2][2] = vectorZ[2];
}

Vec3 AmoebaReferenceHippoNonbondedForce::rotateVectorToQI(const Vec3 v, const double (&mat)[3][3]) const {
    Vec3 rotated;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            rotated[i] += mat[i][j]*v[j];
    return rotated;
}

Vec3 AmoebaReferenceHippoNonbondedForce::rotateVectorFromQI(const Vec3 v, const double (&mat)[3][3]) const {
    Vec3 rotated;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            rotated[i] += mat[j][i]*v[j];
    return rotated;
}

void AmoebaReferenceHippoNonbondedForce::rotateQuadrupoleToQI(const double (&q)[6], double (&r)[6], const double (&mat)[3][3]) const {
    r[QXX] = mat[0][0]*(mat[0][0]*q[QXX] + 2*(mat[0][1]*q[QXY] + mat[0][2]*q[QXZ])) + mat[0][1]*(mat[0][1]*q[QYY] + 2*mat[0][2]*q[QYZ]) + mat[0][2]*mat[0][2]*q[QZZ];
    r[QYY] = mat[1][0]*(mat[1][0]*q[QXX] + 2*(mat[1][1]*q[QXY] + mat[1][2]*q[QXZ])) + mat[1][1]*(mat[1][1]*q[QYY] + 2*mat[1][2]*q[QYZ]) + mat[1][2]*mat[1][2]*q[QZZ];
    r[QXY] = mat[0][0]*mat[1][0]*q[QXX] + mat[0][1]*mat[1][1]*q[QYY] + mat[0][2]*mat[1][2]*q[QZZ] + (mat[0][0]*mat[1][1] + mat[0][1]*mat[1][0])*q[QXY] + (mat[0][0]*mat[1][2] + mat[0][2]*mat[1][0])*q[QXZ] + (mat[0][1]*mat[1][2] + mat[0][2]*mat[1][1])*q[QYZ];
    r[QXZ] = mat[0][0]*mat[2][0]*q[QXX] + mat[0][1]*mat[2][1]*q[QYY] + mat[0][2]*mat[2][2]*q[QZZ] + (mat[0][0]*mat[2][1] + mat[0][1]*mat[2][0])*q[QXY] + (mat[0][0]*mat[2][2] + mat[0][2]*mat[2][0])*q[QXZ] + (mat[0][1]*mat[2][2] + mat[0][2]*mat[2][1])*q[QYZ];
    r[QYZ] = mat[1][0]*mat[2][0]*q[QXX] + mat[1][1]*mat[2][1]*q[QYY] + mat[1][2]*mat[2][2]*q[QZZ] + (mat[1][0]*mat[2][1] + mat[1][1]*mat[2][0])*q[QXY] + (mat[1][0]*mat[2][2] + mat[1][2]*mat[2][0])*q[QXZ] + (mat[1][1]*mat[2][2] + mat[1][2]*mat[2][1])*q[QYZ];
    r[QZZ] = -r[QXX]-r[QYY];
}

void AmoebaReferenceHippoNonbondedForce::applyRotationMatrix() {
    for (int ii = 0; ii < _numParticles; ii++) {
        MultipoleParticleData& p = particleData[ii];
        if (p.multipoleAtomZ >= 0) {
            applyRotationMatrixToParticle(p, &particleData[p.multipoleAtomZ],
                                          p.multipoleAtomX > -1 ? &particleData[p.multipoleAtomX] : NULL,
                                          p.multipoleAtomY > -1 ? &particleData[p.multipoleAtomY] : NULL, p.axisType);
        }
    }
}

void AmoebaReferenceHippoNonbondedForce::computeDirectFieldDampingFactors(const MultipoleParticleData& particle, double r, double& fdamp3, double& fdamp5, double& fdamp7) const {
    double ar = particle.alpha*r;
    double ar2 = ar*ar;
    double ar3 = ar2*ar;
    double ar4 = ar2*ar2;
    double expAR = exp(-ar);
    fdamp3 = 1 - (1 + ar + ar2*(1.0/2))*expAR;
    fdamp5 = 1 - (1 + ar + ar2*(1.0/2) + ar3*(1.0/6))*expAR;
    fdamp7 = 1 - (1 + ar + ar2*(1.0/2) + ar3*(1.0/6) + ar4*(1.0/30))*expAR;
}

void AmoebaReferenceHippoNonbondedForce::computeMutualFieldDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp3, double& fdamp5) const {
    double arI = particleI.alpha*r;
    double arI2 = arI*arI;
    double arI3 = arI2*arI;
    double expARI = exp(-arI);
    if (particleI.alpha == particleJ.alpha) {
        double arI4 = arI3*arI;
        double arI5 = arI4*arI;
        fdamp3 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(7.0/48) + arI4*(1.0/48))*expARI;
        fdamp5 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/144))*expARI;
    }
    else {
        double arJ = particleJ.alpha*r;
        double arJ2 = arJ*arJ;
        double arJ3 = arJ2*arJ;
        double expARJ = exp(-arJ);
        double aI2 = particleI.alpha*particleI.alpha;
        double aJ2 = particleJ.alpha*particleJ.alpha;
        double A = aJ2/(aJ2-aI2);
        double B = aI2/(aI2-aJ2);
        double A2 = A*A;
        double B2 = B*B;
        fdamp3 = 1 - A2*(1 + arI + arI2*(1.0/2))*expARI -
                     B2*(1 + arJ + arJ2*(1.0/2))*expARJ -
                     2*A2*B*(1 + arI)*expARI -
                     2*B2*A*(1 + arJ)*expARJ;
        fdamp5 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6))*expARI -
                     B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6))*expARJ -
                     2*A2*B*(1 + arI + arI2*(1.0/3))*expARI -
                     2*B2*A*(1 + arJ + arJ2*(1.0/3))*expARJ;
    }
}

void AmoebaReferenceHippoNonbondedForce::computeOverlapDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r,
                                double& fdampI1, double& fdampI3, double& fdampI5, double& fdampI7, double& fdampI9,
                                double& fdampJ1, double& fdampJ3, double& fdampJ5, double& fdampJ7, double& fdampJ9,
                                double& fdampIJ1, double& fdampIJ3, double& fdampIJ5, double& fdampIJ7, double& fdampIJ9, double& fdampIJ11) const {
    double arI = particleI.alpha*r;
    double arI2 = arI*arI;
    double arI3 = arI2*arI;
    double arI4 = arI2*arI2;
    double arI5 = arI3*arI2;
    double arI6 = arI3*arI3;
    double expARI = exp(-arI);
    fdampI1 = 1 - (1 + arI*(1.0/2))*expARI;
    fdampI3 = 1 - (1 + arI + arI2*(1.0/2))*expARI;
    fdampI5 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6))*expARI;
    fdampI7 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/30))*expARI;
    fdampI9 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(4.0/105) + arI5*(1.0/210))*expARI;
    if (particleI.alpha == particleJ.alpha) {
        fdampJ1 = fdampI1;
        fdampJ3 = fdampI3;
        fdampJ5 = fdampI5;
        fdampJ7 = fdampI7;
        fdampJ9 = fdampI9;
        double arI7 = arI4*arI3;
        double arI8 = arI4*arI4;
        fdampIJ1 = 1 - (1 + arI*(11.0/16) + arI2*(3.0/16) + arI3*(1.0/48))*expARI;
        fdampIJ3 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(7.0/48) + arI4*(1.0/48))*expARI;
        fdampIJ5 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/144))*expARI;
        fdampIJ7 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/120) + arI6*(1.0/720))*expARI;
        fdampIJ9 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/120) + arI6*(1.0/720) + arI7*(1.0/5040))*expARI;
        fdampIJ11 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/120) + arI6*(1.0/720) + arI7*(1.0/5040) + arI8*(1.0/45360))*expARI;
    }
    else {
        double arJ = particleJ.alpha*r;
        double arJ2 = arJ*arJ;
        double arJ3 = arJ2*arJ;
        double arJ4 = arJ2*arJ2;
        double arJ5 = arJ3*arJ2;
        double arJ6 = arJ3*arJ3;
        double expARJ = exp(-arJ);
        double aI2 = particleI.alpha*particleI.alpha;
        double aJ2 = particleJ.alpha*particleJ.alpha;
        double A = aJ2/(aJ2-aI2);
        double B = aI2/(aI2-aJ2);
        double A2 = A*A;
        double B2 = B*B;
        fdampJ1 = 1 - (1 + arJ*(1.0/2))*expARJ;
        fdampJ3 = 1 - (1 + arJ + arJ2*(1.0/2))*expARJ;
        fdampJ5 = 1 - (1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6))*expARJ;
        fdampJ7 = 1 - (1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6) + arJ4*(1.0/30))*expARJ;
        fdampJ9 = 1 - (1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6) + 4*arJ4*(1.0/105) + arJ5*(1.0/210))*expARJ;
        fdampIJ1 = 1 - A2*(1 + 2*B + arI*(1.0/2))*expARI -
                       B2*(1 + 2*A + arJ*(1.0/2))*expARJ;
        fdampIJ3 = 1 - A2*(1 + arI + arI2*(1.0/2))*expARI -
                       B2*(1 + arJ + arJ2*(1.0/2))*expARJ -
                       2*A2*B*(1 + arI)*expARI -
                       2*B2*A*(1 + arJ)*expARJ;
        fdampIJ5 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6))*expARI -
                       B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6))*expARJ -
                       2*A2*B*(1 + arI + arI2*(1.0/3))*expARI -
                       2*B2*A*(1 + arJ + arJ2*(1.0/3))*expARJ;
        fdampIJ7 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/30))*expARI -
                       B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6) + arJ4*(1.0/30))*expARJ -
                       2*A2*B*(1 + arI + arI2*(2.0/5) + arI3*(1.0/15))*expARI -
                       2*B2*A*(1 + arJ + arJ2*(2.0/5) + arJ3*(1.0/15))*expARJ;
        fdampIJ9 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*4*(1.0/105) + arI5*(1.0/210))*expARI -
                       B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6) + arJ4*4*(1.0/105) + arJ5*(1.0/210))*expARJ -
                       2*A2*B*(1 + arI + arI2*(3.0/7) + arI3*(2.0/21) + arI4*(1.0/105))*expARI -
                       2*B2*A*(1 + arJ + arJ2*(3.0/7) + arJ3*(2.0/21) + arJ4*(1.0/105))*expARJ;
        fdampIJ11 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(5.0/126) + arI5*(2.0/315) + arI6*(1.0/1890))*expARI -
                        B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6) + arJ4*(5.0/126) + arJ5*(2.0/315) + arJ6*(1.0/1890))*expARJ -
                        2*A2*B*(1 + arI + arI2*(4.0/9) + arI3*(1.0/9) + arI4*(1.0/63) + arI5*(1.0/945))*expARI -
                        2*B2*A*(1 + arJ + arJ2*(4.0/9) + arJ3*(1.0/9) + arJ4*(1.0/63) + arJ5*(1.0/945))*expARJ;
    }
}

void AmoebaReferenceHippoNonbondedForce::computeDispersionDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp, double& ddamp) const {
    double arI = particleI.alpha*r;
    double arI2 = arI*arI;
    double arI3 = arI2*arI;
    double expARI = exp(-arI);
    double fdamp3, fdamp5;
    if (particleI.alpha == particleJ.alpha) {
        double arI4 = arI3*arI;
        double arI5 = arI4*arI;
        fdamp3 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(7.0/48) + arI4*(1.0/48))*expARI;
        fdamp5 = 1 - (1 + arI + arI2*(1.0/2) + arI3*(1.0/6) + arI4*(1.0/24) + arI5*(1.0/144))*expARI;
        ddamp = particleI.alpha*(arI5 - 3*arI3 - 3*arI2)*expARI*(1.0/96);
    }
    else {
        double arJ = particleJ.alpha*r;
        double arJ2 = arJ*arJ;
        double arJ3 = arJ2*arJ;
        double expARJ = exp(-arJ);
        double aI2 = particleI.alpha*particleI.alpha;
        double aJ2 = particleJ.alpha*particleJ.alpha;
        double A = aJ2/(aJ2-aI2);
        double B = aI2/(aI2-aJ2);
        double A2 = A*A;
        double B2 = B*B;
        fdamp3 = 1 - A2*(1 + arI + arI2*(1.0/2))*expARI -
                     B2*(1 + arJ + arJ2*(1.0/2))*expARJ -
                     2*A2*B*(1 + arI)*expARI -
                     2*B2*A*(1 + arJ)*expARJ;
        fdamp5 = 1 - A2*(1 + arI + arI2*(1.0/2) + arI3*(1.0/6))*expARI -
                     B2*(1 + arJ + arJ2*(1.0/2) + arJ3*(1.0/6))*expARJ -
                     2*A2*B*(1 + arI + arI2*(1.0/3))*expARI -
                     2*B2*A*(1 + arJ + arJ2*(1.0/3))*expARJ;
        ddamp = (A2*arI2*particleI.alpha*expARI*(r*particleI.alpha + 4*B - 1) +
                (B2*arJ2*particleJ.alpha*expARJ*(r*particleJ.alpha + 4*A - 1)))*(1.0/4);
    }
    fdamp = 1.5*fdamp5 - 0.5*fdamp3;
}

void AmoebaReferenceHippoNonbondedForce::computeRepulsionDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r,
            double& fdamp1, double& fdamp3, double& fdamp5, double& fdamp7, double& fdamp9, double& fdamp11) const {
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;
    double r5 = r3*r2;
    double r6 = r3*r3;
    double aI2 = 0.5*particleI.pauliAlpha;
    double arI2 = aI2*r;
    double expI = exp(-arI2);
    double aI2_2 = aI2*aI2;
    double aI2_3 = aI2_2*aI2;
    double aI2_4 = aI2_2*aI2_2;
    double aI2_5 = aI2_3*aI2_2;
    double aI2_6 = aI2_3*aI2_3;
    double fexp, fexp1, fexp2, fexp3, fexp4, fexp5, pre;
    if (particleI.pauliAlpha == particleJ.pauliAlpha) {
        double r7 = r4*r3;
        double r8 = r4*r4;
        double aI2_7 = aI2_4*aI2_3;
        pre = 128;
        fexp = (r + aI2*r2 + aI2_2*r3*(1.0/3))*expI;
        fexp1 = (aI2_2*r3 + aI2_3*r4)*expI*(1.0/3);
        fexp2 = aI2_4*expI*r5*(1.0/9);
        fexp3 = aI2_5*expI*r6*(1.0/45);
        fexp4 = (aI2_5*r6 + aI2_6*r7)*expI*(1.0/315);
        fexp5 = (aI2_5*r6 + aI2_6*r7 + aI2_7*r8*(1.0/3))*expI*(1.0/945);
    }
    else {
        double aJ2 = 0.5*particleJ.pauliAlpha;
        double arJ2 = aJ2*r;
        double expJ = exp(-arJ2);
        double aJ2_2 = aJ2*aJ2;
        double aJ2_3 = aJ2_2*aJ2;
        double aJ2_4 = aJ2_2*aJ2_2;
        double aJ2_5 = aJ2_3*aJ2_2;
        double scale = 1/(aI2_2-aJ2_2);
        double aI2aJ2expI = aI2*aJ2*expI;
        double  aI2aJ2expJ = aI2*aJ2*expJ;
        pre = 8192*aI2_3*aJ2_3*(scale*scale*scale*scale);
        double tmp = 4*aI2*aJ2*scale;
        fexp = (arI2-tmp)*expJ + (arJ2+tmp)*expI;
        fexp1 = (r2 - (4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (r2 + (4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp2 = (r2*(1.0/3) + aJ2*r3*(1.0/3) - ((4.0/3)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (r2*(1.0/3) + aI2*r3*(1.0/3) + ((4.0/3)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp3 = (aJ2_2*r4*(1.0/15) + aJ2*r3*(1.0/5) + r2*(1.0/5) - ((4.0/15)*aJ2_3*r3 + (8.0/5)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_2*r4*(1.0/15) + aI2*r3*(1.0/5) + r2*(1.0/5) + ((4.0/15)*aI2_3*r3 + (8.0/5)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp4 = (aJ2_3*r5*(1.0/105) + (2.0/35)*aJ2_2*r4 + aJ2*r3*(1.0/7) + r2*(1.0/7) - ((4.0/105)*aJ2_4*r4 + (8.0/21)*aJ2_3*r3 + (12.0/7)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_3*r5*(1.0/105) + (2.0/35)*aI2_2*r4 + aI2*r3*(1.0/7) + r2*(1.0/7) + ((4.0/105)*aI2_4*r4 + (8.0/21)*aI2_3*r3 + (12.0/7)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
        fexp5 = (aJ2_4*r6*(1.0/945) + (2.0/189)*aJ2_3*r5 + aJ2_2*r4*(1.0/21) + aJ2*r3*(1.0/9) + r2*(1.0/9) - ((4.0/945)*aJ2_5*r5 + (4.0/63)*aJ2_4*r4 + (4.0/9)*aJ2_3*r3 + (16.0/9)*aJ2_2*r2 + 4*aJ2*r + 4)*scale)*aI2aJ2expJ +
                (aI2_4*r6*(1.0/945) + (2.0/189)*aI2_3*r5 + aI2_2*r4*(1.0/21) + aI2*r3*(1.0/9) + r2*(1.0/9) + ((4.0/945)*aI2_5*r5 + (4.0/63)*aI2_4*r4 + (4.0/9)*aI2_3*r3 + (16.0/9)*aI2_2*r2 + 4*aI2*r + 4)*scale)*aI2aJ2expI;
    }
    fexp = fexp/r;
    fexp1 = fexp1/r3;
    fexp2 = 3*fexp2/r5;
    fexp3 = 15*fexp3/(r5*r2);
    fexp4 = 105*fexp4/(r5*r4);
    fexp5 = 945*fexp5/(r5*r6);
    fdamp1 = 0.5*pre*fexp*fexp;
    fdamp3 = pre*fexp*fexp1;
    fdamp5 = pre*(fexp*fexp2 + fexp1*fexp1);
    fdamp7 = pre*(fexp*fexp3 + 3*fexp1*fexp2);
    fdamp9 = pre*(fexp*fexp4 + 4*fexp1*fexp3 + 3*fexp2*fexp2);
    fdamp11 = pre*(fexp*fexp5 + 5*fexp1*fexp4 + 10*fexp2*fexp3);
}

void AmoebaReferenceHippoNonbondedForce::calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI,
                                                                             const MultipoleParticleData& particleJ) {
    Vec3 deltaR = particleJ.position - particleI.position;
    double r = sqrt(deltaR.dot(deltaR));
    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rInv3 = rInv*rInv2;
    double rInv5 = rInv3*rInv2;
    double rInv7 = rInv5*rInv2;

    // Calculate the field at particle I due to multipoles at particle J.

    double fdamp3, fdamp5, fdamp7;
    computeDirectFieldDampingFactors(particleJ, r, fdamp3, fdamp5, fdamp7);
    double scale = 1;
    auto exception = exceptions.find(make_pair(particleI.index, particleJ.index));
    if (exception != exceptions.end())
        scale = exception->second.dipoleMultipoleScale;
    double rr3 = scale*rInv3;
    double rr3j = scale*fdamp3*rInv3;
    double rr5j = scale*3*fdamp5*rInv5;
    double rr7j = scale*15*fdamp7*rInv7;
    Vec3 qDotDelta(deltaR[0]*particleJ.quadrupole[QXX] + deltaR[1]*particleJ.quadrupole[QXY] + deltaR[2]*particleJ.quadrupole[QXZ],
                   deltaR[0]*particleJ.quadrupole[QXY] + deltaR[1]*particleJ.quadrupole[QYY] + deltaR[2]*particleJ.quadrupole[QYZ],
                   deltaR[0]*particleJ.quadrupole[QXZ] + deltaR[1]*particleJ.quadrupole[QYZ] + deltaR[2]*particleJ.quadrupole[QZZ]);
    double dipoleDelta = particleJ.dipole.dot(deltaR);
    double qdpoleDelta = qDotDelta.dot(deltaR);
    double factor = rr3*particleJ.coreCharge + rr3j*particleJ.valenceCharge - rr5j*dipoleDelta + rr7j*qdpoleDelta;
    Vec3 field = deltaR*factor + particleJ.dipole*rr3j - qDotDelta*2*rr5j;
    _fixedMultipoleField[particleI.index] -= field;
}

void AmoebaReferenceHippoNonbondedForce::calculateFixedMultipoleField() {
    for (int i = 0; i < _numParticles; i++)
        for (int j = 0; j < _numParticles; j++)
            if (i != j)
                calculateFixedMultipoleFieldPairIxn(particleData[i], particleData[j]);
}

void AmoebaReferenceHippoNonbondedForce::initializeInducedDipoles() {
    _inducedDipole.resize(_numParticles);
    for (int i = 0; i < _numParticles; i++)
        _inducedDipole[i] = _fixedMultipoleField[i];
}

void AmoebaReferenceHippoNonbondedForce::calculateInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                                                        const MultipoleParticleData& particleJ) {
    int i = particleI.index;
    int j = particleJ.index;
    if (i == j)
        return;

    Vec3 deltaR = particleJ.position - particleI.position;
    double r = sqrt(deltaR.dot(deltaR));
    double fdamp3, fdamp5;
    computeMutualFieldDampingFactors(particleI, particleJ, r, fdamp3, fdamp5);
    auto exception = exceptions.find(make_pair(particleI.index, particleJ.index));
    if (exception != exceptions.end()) {
        fdamp3 *= exception->second.dipoleDipoleScale;
        fdamp5 *= exception->second.dipoleDipoleScale;
    }
    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rInv3 = rInv*rInv2;
    double scale3 = -fdamp3*rInv3;
    double scale5 = 3*fdamp5*rInv3*rInv2;
    _inducedDipoleField[i] += _inducedDipole[j]*scale3 + deltaR*scale5*(_inducedDipole[j].dot(deltaR));
    _inducedDipoleField[j] += _inducedDipole[i]*scale3 + deltaR*scale5*(_inducedDipole[i].dot(deltaR));
}

void AmoebaReferenceHippoNonbondedForce::calculateInducedDipoleFields(const vector<MultipoleParticleData>& particleData, int optOrder) {
    // Initialize the fields to zero.

    Vec3 zeroVec(0.0, 0.0, 0.0);
    std::fill(_inducedDipoleField.begin(), _inducedDipoleField.end(), zeroVec);

    // Add fields from all induced dipoles.

    for (int ii = 0; ii < _numParticles; ii++)
        for (int jj = ii; jj < _numParticles; jj++)
            calculateInducedDipolePairIxns(particleData[ii], particleData[jj]);
}

void AmoebaReferenceHippoNonbondedForce::convergeInduceDipolesByExtrapolation(const vector<MultipoleParticleData>& particleData) {
    // Start by storing the direct dipoles as PT0

    _ptDipoleD.resize(_maxPTOrder);
    _ptDipoleD[0].resize(_numParticles);
    for (int atom = 0; atom < _numParticles; ++atom)
        _ptDipoleD[0][atom] = _inducedDipole[atom];

    // Recursively apply alpha.Tau to the µ_(n) components to generate µ_(n+1), and store the result

    vector<double> zeros(6, 0.0);
    for (int order = 1; order < _maxPTOrder; ++order) {
        calculateInducedDipoleFields(particleData, order-1);
        _ptDipoleD[order].resize(_numParticles);
        for (int atom = 0; atom < _numParticles; ++atom) {
            _inducedDipole[atom] = _inducedDipoleField[atom] * particleData[atom].polarizability;
            _ptDipoleD[order][atom] = _inducedDipole[atom];
        }
    }

    // Take a linear combination of the µ_(n) components to form the total dipole
    
    _inducedDipole = vector<Vec3>(_numParticles, Vec3());
    for (int order = 0; order < _maxPTOrder; ++order)
        for (int atom = 0; atom < _numParticles; ++atom)
            _inducedDipole[atom] += _ptDipoleD[order][atom] * _extPartCoefficients[order];
    calculateInducedDipoleFields(particleData, _maxPTOrder-1);
}

void AmoebaReferenceHippoNonbondedForce::calculateInducedDipoles() {
    // calculate fixed electric fields

    initializeVec3Vector(_fixedMultipoleField);
    calculateFixedMultipoleField();

    // initialize inducedDipoles

    for (int i = 0; i < _numParticles; i++)
        _fixedMultipoleField[i] *= particleData[i].polarizability;
    _inducedDipole.resize(_numParticles);
    _inducedDipoleField.resize(_fixedMultipoleField.size());
    initializeInducedDipoles();
    convergeInduceDipolesByExtrapolation(particleData);
}

double AmoebaReferenceHippoNonbondedForce::calculateElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                         double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const {
    double dir = particleI.qiDipole[2]*r;
    Vec3 qxI = Vec3(particleI.qiQuadrupole[QXX], particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QXZ]);
    Vec3 qyI = Vec3(particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QYY], particleI.qiQuadrupole[QYZ]);
    Vec3 qzI = Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    Vec3 qi = r*Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    double qir = qi[2]*r;
    double dkr = particleK.qiDipole[2]*r;
    Vec3 qxK = Vec3(particleK.qiQuadrupole[QXX], particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QXZ]);
    Vec3 qyK = Vec3(particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QYY], particleK.qiQuadrupole[QYZ]);
    Vec3 qzK = Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    Vec3 qk = r*Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    double qkr = qk[2]*r;
    double dik = particleI.qiDipole.dot(particleK.qiDipole);
    double qik = qi.dot(qk);
    double diqk = particleI.qiDipole.dot(qk);
    double dkqi = particleK.qiDipole.dot(qi);
    double qiqk = 2*(qxI[1]*qxK[1]+qxI[2]*qxK[2]+qyI[2]*qyK[2]) + qxI[0]*qxK[0] + qyI[1]*qyK[1] + qzI[2]*qzK[2];

    // Additional intermediates involving moments and distance.

    Vec3 dirCross(particleI.qiDipole[1]*r, -particleI.qiDipole[0]*r, 0);
    Vec3 dkrCross(particleK.qiDipole[1]*r, -particleK.qiDipole[0]*r, 0);
    Vec3 dikCross = particleI.qiDipole.cross(particleK.qiDipole);
    Vec3 qirCross(qi[1]*r, -qi[0]*r, 0);
    Vec3 qkrCross(qk[1]*r, -qk[0]*r, 0);
    Vec3 qikCross = qk.cross(qi);
    Vec3 qikTemp(qxI.dot(qk), qyI.dot(qk), qzI.dot(qk));
    Vec3 qkiTemp(qxK.dot(qi), qyK.dot(qi), qzK.dot(qi));
    Vec3 qikrCross(-r*qikTemp[1], r*qikTemp[0], 0);
    Vec3 qkirCross(-r*qkiTemp[1], r*qkiTemp[0], 0);
    Vec3 diqkTemp(particleI.qiDipole.dot(qxK), particleI.qiDipole.dot(qyK), particleI.qiDipole.dot(qzK));
    Vec3 dkqiTemp(particleK.qiDipole.dot(qxI), particleK.qiDipole.dot(qyI), particleK.qiDipole.dot(qzI));
    Vec3 diqkrCross(-r*diqkTemp[1], r*diqkTemp[0], 0);
    Vec3 dkqirCross(-r*dkqiTemp[1], r*dkqiTemp[0], 0);
    Vec3 dqik = particleI.qiDipole.cross(qk) + particleK.qiDipole.cross(qi) - 2*(qxI.cross(qxK) + qyI.cross(qyK) + qzI.cross(qzK));

    // Get reciprocal distance terms for this interaction.

    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rr1 = rInv;
    double rr3 = rr1*rInv2;
    double rr5 = 3*rr3*rInv2;
    double rr7 = 5*rr5*rInv2;
    double rr9 = 7*rr7*rInv2;
    double rr11 = 9*rr9*rInv2;

    // Find damped multipole intermediates and energy value.

    double term1 = particleI.coreCharge*particleK.coreCharge;
    double term1i = particleK.coreCharge*particleI.valenceCharge;
    double term2i = particleK.coreCharge*dir;
    double term3i = particleK.coreCharge*qir;
    double term1k = particleI.coreCharge*particleK.valenceCharge;
    double term2k = -particleI.coreCharge*dkr;
    double term3k = particleI.coreCharge*qkr;
    double term1ik = particleI.valenceCharge*particleK.valenceCharge;
    double term2ik = particleK.valenceCharge*dir - particleI.valenceCharge*dkr + dik;
    double term3ik = particleI.valenceCharge*qkr + particleK.valenceCharge*qir - dir*dkr + 2*(dkqi-diqk+qiqk);
    double term4ik = dir*qkr - dkr*qir - 4*qik;
    double term5ik = qir*qkr;
    double fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    double fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    double fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(particleI, particleK, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    double rr1i = fdampI1*rr1;
    double rr3i = fdampI3*rr3;
    double rr5i = fdampI5*rr5;
    double rr7i = fdampI7*rr7;
    double rr1k = fdampK1*rr1;
    double rr3k = fdampK3*rr3;
    double rr5k = fdampK5*rr5;
    double rr7k = fdampK7*rr7;
    double rr1ik = fdampIK1*rr1;
    double rr3ik = fdampIK3*rr3;
    double rr5ik = fdampIK5*rr5;
    double rr7ik = fdampIK7*rr7;
    double rr9ik = fdampIK9*rr9;
    double rr11ik = fdampIK11*rr11;
    double scale = _electric;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end())
        scale *= exception->second.multipoleMultipoleScale;
    double energy = scale*(term1*rr1 + term4ik*rr7ik + term5ik*rr9ik +
                           term1i*rr1i + term1k*rr1k + term1ik*rr1ik +
                           term2i*rr3i + term2k*rr3k + term2ik*rr3ik +
                           term3i*rr5i + term3k*rr5k + term3ik*rr5ik);

    // Find damped multipole intermediates for force and torque.

    double de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik +
                term1i*rr3i + term1k*rr3k + term1ik*rr3ik +
                term2i*rr5i + term2k*rr5k + term2ik*rr5ik +
                term3i*rr7i + term3k*rr7k + term3ik*rr7ik;
    term1 = -particleK.coreCharge*rr3i - particleK.valenceCharge*rr3ik + dkr*rr5ik - qkr*rr7ik;
    double term2 = particleI.coreCharge*rr3k + particleI.valenceCharge*rr3ik + dir*rr5ik + qir*rr7ik;
    double term3 = 2 * rr5ik;
    double term4 = -2 * (particleK.coreCharge*rr5i+particleK.valenceCharge*rr5ik-dkr*rr7ik+qkr*rr9ik);
    double term5 = -2 * (particleI.coreCharge*rr5k+particleI.valenceCharge*rr5ik+dir*rr7ik+qir*rr9ik);
    double term6 = 4 * rr7ik;

    // Compute the force and torque.

    force += scale*(Vec3(0, 0, de*r) + term1*particleI.qiDipole + term2*particleK.qiDipole +
            term3*(diqkTemp-dkqiTemp) + term4*qi + term5*qk + term6*(qikTemp+qkiTemp));
    torqueI += scale*(-rr3ik*dikCross + term1*dirCross + term3*(dqik+dkqirCross) + term4*qirCross - term6*(qikrCross+qikCross));
    torqueK += scale*(rr3ik*dikCross + term2*dkrCross - term3*(dqik+diqkrCross) + term5*qkrCross - term6*(qkirCross-qikCross));
    return energy;
}

void AmoebaReferenceHippoNonbondedForce::calculateInducedDipolePairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                       Vec3 deltaR, double r, Vec3& force, Vec3& torqueI, Vec3& torqueK, Vec3& labForce) const {
    double r2 = r*r;

    // Intermediates involving moments and separation distance

    double dir = particleI.qiDipole[2]*r;
    Vec3 qxI = Vec3(particleI.qiQuadrupole[QXX], particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QXZ]);
    Vec3 qyI = Vec3(particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QYY], particleI.qiQuadrupole[QYZ]);
    Vec3 qzI = Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    Vec3 qi = r*Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    double qir = qi[2]*r;
    double dkr = particleK.qiDipole[2]*r;
    Vec3 qxK = Vec3(particleK.qiQuadrupole[QXX], particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QXZ]);
    Vec3 qyK = Vec3(particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QYY], particleK.qiQuadrupole[QYZ]);
    Vec3 qzK = Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    Vec3 qk = r*Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    double qkr = qk[2]*r;
    const Vec3& ui = particleI.qiInducedDipole;
    const Vec3& uk = particleK.qiInducedDipole;
    double uir = ui[2]*r;
    double ukr = uk[2]*r;

    // Get reciprocal distance terms for this interaction.

    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rr1 = rInv;
    double rr3 = rr1*rInv2;
    double rr5 = 3*rr3*rInv2;
    double rr7 = 5*rr5*rInv2;
    double rr9 = 7*rr7*rInv2;

    // Apply charge penetration damping to scale factors.

    double fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    double fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    double fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(particleI, particleK, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    double dmScale = _electric;
    double ddScale = _electric;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end()) {
        dmScale *= exception->second.dipoleMultipoleScale;
        ddScale *= exception->second.dipoleDipoleScale;
    }
    double rr3core = rr3*dmScale;
    double rr5core = rr5*dmScale;
    double rr3i = rr3*fdampI3*dmScale;
    double rr5i = rr5*fdampI5*dmScale;
    double rr7i = rr7*fdampI7*dmScale;
    double rr9i = rr9*fdampI9*dmScale;
    double rr3k = rr3*fdampK3*dmScale;
    double rr5k = rr5*fdampK5*dmScale;
    double rr7k = rr7*fdampK7*dmScale;
    double rr9k = rr9*fdampK9*dmScale;
    double rr5ik = rr5*fdampIK5*ddScale;
    double rr7ik = rr7*fdampIK7*ddScale;

    // Get the induced dipole field used for dipole torques.

    Vec3 torqueFieldI = rr3i*uk;
    torqueFieldI[2] -= rr5i*ukr*r;
    Vec3 torqueFieldK = rr3k*ui;
    torqueFieldK[2] -= rr5k*uir*r;

    // Get induced dipole field gradient used for quadrupole torques.

    Vec3 dtorqueFieldI = 2*r*rr5i*uk;
    dtorqueFieldI[2] -= r2*rr7i*ukr;
    Vec3 dtorqueFieldK = -2*r*rr5k*ui;
    dtorqueFieldK[2] += r2*rr7k*uir;

    // Get the field gradient for direct polarization force

    double ti[6], tk[6];
    ti[QXX] = particleI.valenceCharge*rr3i + particleI.coreCharge*rr3core + dir*rr5i - qxI[0]*2*rr5i + qi[2]*r*rr7i;
    tk[QXX] = particleK.valenceCharge*rr3k + particleK.coreCharge*rr3core - dkr*rr5k - qxK[0]*2*rr5k + qk[2]*r*rr7k;
    ti[QYY] = particleI.valenceCharge*rr3i + particleI.coreCharge*rr3core + dir*rr5i - qyI[1]*2*rr5i + qi[2]*r*rr7i;
    tk[QYY] = particleK.valenceCharge*rr3k + particleK.coreCharge*rr3core - dkr*rr5k - qyK[1]*2*rr5k + qk[2]*r*rr7k;
    ti[QZZ] = particleI.valenceCharge*(rr3i-rr5i*r2) + particleI.coreCharge*(rr3core-rr5core*r2) + particleI.qiDipole[2]*2*rr5i*r -
              dir*(rr7i*r2-rr5i) - qzI[2]*2*rr5i + qi[2]*5*rr7i*r - qir*rr9i*r2;
    tk[QZZ] = particleK.valenceCharge*(rr3k-rr5k*r2) + particleK.coreCharge*(rr3core-rr5core*r2) - particleK.qiDipole[2]*2*rr5k*r +
              dkr*(rr7k*r2-rr5k) - qzK[2]*2*rr5k + qk[2]*5*rr7k*r - qkr*rr9k*r2;
    ti[QXY] = -qxI[1]*2*rr5i;
    tk[QXY] = -qxK[1]*2*rr5k;
    ti[QXZ] = particleI.qiDipole[0]*rr5i*r - qxI[2]*2*rr5i + qi[0]*2*rr7i*r;
    tk[QXZ] = -particleK.qiDipole[0]*rr5k*r - qxK[2]*2*rr5k + qk[0]*2*rr7k*r;
    ti[QYZ] = particleI.qiDipole[1]*rr5i*r - qyI[2]*2*rr5i + qi[1]*2*rr7i*r;
    tk[QYZ] = -particleK.qiDipole[1]*rr5k*r - qyK[2]*2*rr5k + qk[1]*2*rr7k*r;

    // Get the dEp/dR terms for chgpen direct polarization force.

    double depx = ti[QXX]*uk[0] + ti[QXY]*uk[1] + ti[QXZ]*uk[2] - tk[QXX]*ui[0] - tk[QXY]*ui[1] - tk[QXZ]*ui[2];
    double depy = ti[QXY]*uk[0] + ti[QYY]*uk[1] + ti[QYZ]*uk[2] - tk[QXY]*ui[0] - tk[QYY]*ui[1] - tk[QYZ]*ui[2];
    double depz = ti[QXZ]*uk[0] + ti[QYZ]*uk[1] + ti[QZZ]*uk[2] - tk[QXZ]*ui[0] - tk[QYZ]*ui[1] - tk[QZZ]*ui[2];
    force += Vec3(depx, depy, depz);

    // Get the dtau/dr terms used for OPT polarization force.

    int i = particleI.index;
    int k = particleK.index;
    for (int j = 0; j < _maxPTOrder-1; j++) {
        Vec3 extDipole1 = _ptDipoleD[j][i];
        double uirm = extDipole1.dot(deltaR);
        for (int m = 0; m < _maxPTOrder-1-j; m++) {
            Vec3 extDipole2 = _ptDipoleD[m][k];
            double ukrm = extDipole2.dot(deltaR);
            double term1 = 2*rr5ik;
            double term2 = term1*deltaR[0];
            double term3 = rr5ik - rr7ik*deltaR[0]*deltaR[0];
            double tixx = extDipole1[0]*term2 + uirm*term3;
            double tkxx = extDipole2[0]*term2 + ukrm*term3;
            term2 = term1*deltaR[1];
            term3 = rr5ik - rr7ik*deltaR[1]*deltaR[1];
            double tiyy = extDipole1[1]*term2 + uirm*term3;
            double tkyy = extDipole2[1]*term2 + ukrm*term3;
            term2 = term1*deltaR[2];
            term3 = rr5ik - rr7ik*deltaR[2]*deltaR[2];
            double tizz = extDipole1[2]*term2 + uirm*term3;
            double tkzz = extDipole2[2]*term2 + ukrm*term3;
            term1 = rr5ik*deltaR[1];
            term2 = rr5ik*deltaR[0];
            term3 = deltaR[1] * (rr7ik*deltaR[0]);
            double tixy = extDipole1[0]*term1 + extDipole1[1]*term2 - uirm*term3;
            double tkxy = extDipole2[0]*term1 + extDipole2[1]*term2 - ukrm*term3;
            term1 = rr5ik * deltaR[2];
            term3 = deltaR[2] * (rr7ik*deltaR[0]);
            double tixz = extDipole1[0]*term1 + extDipole1[2]*term2 - uirm*term3;
            double tkxz = extDipole2[0]*term1 + extDipole2[2]*term2 - ukrm*term3;
            term2 = rr5ik*deltaR[1];
            term3 = deltaR[2] * (rr7ik*deltaR[1]);
            double tiyz = extDipole1[1]*term1 + extDipole1[2]*term2 - uirm*term3;
            double tkyz = extDipole2[1]*term1 + extDipole2[2]*term2 - ukrm*term3;
            double depx = tixx*extDipole2[0] + tkxx*extDipole1[0]
                 + tixy*extDipole2[1] + tkxy*extDipole1[1]
                 + tixz*extDipole2[2] + tkxz*extDipole1[2];
            double depy = tixy*extDipole2[0] + tkxy*extDipole1[0]
                 + tiyy*extDipole2[1] + tkyy*extDipole1[1]
                 + tiyz*extDipole2[2] + tkyz*extDipole1[2];
            double depz = tixz*extDipole2[0] + tkxz*extDipole1[0]
                 + tiyz*extDipole2[1] + tkyz*extDipole1[1]
                 + tizz*extDipole2[2] + tkzz*extDipole1[2];
            labForce += 0.5*_extPartCoefficients[j+m+1]*Vec3(depx, depy, depz);
        }
    }

    // Torque is induced field and gradient cross permanent moments.

    torqueI += torqueFieldI.cross(particleI.qiDipole);
    torqueI[0] += -qxI[1]*dtorqueFieldI[0] - 2*qyI[2]*dtorqueFieldI[2] + (qzI[2]-qyI[1])*dtorqueFieldI[1];
    torqueI[1] += qxI[1]*dtorqueFieldI[1] + 2*qxI[2]*dtorqueFieldI[2] + (qxI[0]-qzI[2])*dtorqueFieldI[0];
    torqueI[2] += qyI[2]*dtorqueFieldI[0] - qxI[2]*dtorqueFieldI[1];
    torqueK += torqueFieldK.cross(particleK.qiDipole);
    torqueK[0] += -qxK[1]*dtorqueFieldK[0] - 2*qyK[2]*dtorqueFieldK[2] + (qzK[2]-qyK[1])*dtorqueFieldK[1];
    torqueK[1] += qxK[1]*dtorqueFieldK[1] + 2*qxK[2]*dtorqueFieldK[2] + (qxK[0]-qzK[2])*dtorqueFieldK[0];
    torqueK[2] += qyK[2]*dtorqueFieldK[0] - qxK[2]*dtorqueFieldK[1];
}

double AmoebaReferenceHippoNonbondedForce::calculateDispersionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                      double r, Vec3& force) const {
    // Compute the undamped force and energy.
    
    double r2 = r*r;
    double energy = -particleI.c6*particleK.c6/(r2*r2*r2);
    double dEnergydR = -6*energy/r;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end()) {
        energy *= exception->second.dispersionScale;
        dEnergydR *= exception->second.dispersionScale;
    }
    
    // Apply the damping function.
    
    double fdamp, ddamp;
    computeDispersionDampingFactors(particleI, particleK, r, fdamp, ddamp);
    dEnergydR = dEnergydR*fdamp*fdamp + 2*energy*fdamp*ddamp;
    energy *= fdamp*fdamp;

    // Accumulate the forces.
    
    force[2] -= dEnergydR;
    return energy;
}

double AmoebaReferenceHippoNonbondedForce::calculateRepulsionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                     double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const {
    double dir = particleI.qiDipole[2]*r;
    Vec3 qxI = Vec3(particleI.qiQuadrupole[QXX], particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QXZ]);
    Vec3 qyI = Vec3(particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QYY], particleI.qiQuadrupole[QYZ]);
    Vec3 qzI = Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    Vec3 qi = r*Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    double qir = qi[2]*r;
    double dkr = particleK.qiDipole[2]*r;
    Vec3 qxK = Vec3(particleK.qiQuadrupole[QXX], particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QXZ]);
    Vec3 qyK = Vec3(particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QYY], particleK.qiQuadrupole[QYZ]);
    Vec3 qzK = Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    Vec3 qk = r*Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    double qkr = qk[2]*r;
    double dik = particleI.qiDipole.dot(particleK.qiDipole);
    double qik = qi.dot(qk);
    double diqk = particleI.qiDipole.dot(qk);
    double dkqi = particleK.qiDipole.dot(qi);
    double qiqk = 2*(qxI[1]*qxK[1]+qxI[2]*qxK[2]+qyI[2]*qyK[2]) + qxI[0]*qxK[0] + qyI[1]*qyK[1] + qzI[2]*qzK[2];

    // Additional intermediates involving moments and distance.

    Vec3 dirCross(particleI.qiDipole[1]*r, -particleI.qiDipole[0]*r, 0);
    Vec3 dkrCross(particleK.qiDipole[1]*r, -particleK.qiDipole[0]*r, 0);
    Vec3 dikCross = particleI.qiDipole.cross(particleK.qiDipole);
    Vec3 qirCross(qi[1]*r, -qi[0]*r, 0);
    Vec3 qkrCross(qk[1]*r, -qk[0]*r, 0);
    Vec3 qikCross = qk.cross(qi);
    Vec3 qikTemp(qxI.dot(qk), qyI.dot(qk), qzI.dot(qk));
    Vec3 qkiTemp(qxK.dot(qi), qyK.dot(qi), qzK.dot(qi));
    Vec3 qikrCross(-r*qikTemp[1], r*qikTemp[0], 0);
    Vec3 qkirCross(-r*qkiTemp[1], r*qkiTemp[0], 0);
    Vec3 diqkTemp(particleI.qiDipole.dot(qxK), particleI.qiDipole.dot(qyK), particleI.qiDipole.dot(qzK));
    Vec3 dkqiTemp(particleK.qiDipole.dot(qxI), particleK.qiDipole.dot(qyI), particleK.qiDipole.dot(qzI));
    Vec3 diqkrCross(-r*diqkTemp[1], r*diqkTemp[0], 0);
    Vec3 dkqirCross(-r*dkqiTemp[1], r*dkqiTemp[0], 0);
    Vec3 dqik = particleI.qiDipole.cross(qk) + particleK.qiDipole.cross(qi) - 2*(qxI.cross(qxK) + qyI.cross(qyK) + qzI.cross(qzK));

    // Get reciprocal distance terms for this interaction.

    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rr1 = rInv;
    double rr3 = rr1*rInv2;

    // Compute damping coefficients.

    double fdamp1, fdamp3, fdamp5, fdamp7, fdamp9, fdamp11;
    computeRepulsionDampingFactors(particleI, particleK, r, fdamp1, fdamp3, fdamp5, fdamp7, fdamp9, fdamp11);

    // Calculate intermediate terms needed for the energy

    double eterm1 = particleI.pauliQ*particleK.pauliQ;
    double eterm2 = particleK.pauliQ*dir - particleI.pauliQ*dkr + dik;
    double eterm3 = particleI.pauliQ*qkr + particleK.pauliQ*qir - dir*dkr + 2*(dkqi-diqk+qiqk);
    double eterm4 = dir*qkr - dkr*qir - 4*qik;
    double eterm5 = qir*qkr;
    double eterm = eterm1*fdamp1 + eterm2*fdamp3 + eterm3*fdamp5 + eterm4*fdamp7 + eterm5*fdamp9;

    // Compute the energy.

    double sizik = particleI.pauliK*particleK.pauliK;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end())
        sizik *= exception->second.repulsionScale;
    double energy = sizik*eterm*rr1;

    // Calculate intermediate terms for force and torque

    double de = eterm1*fdamp3 + eterm2*fdamp5 + eterm3*fdamp7 + eterm4*fdamp9 + eterm5*fdamp11;
    double term1 = -particleK.pauliQ*fdamp3 + dkr*fdamp5 - qkr*fdamp7;
    double term2 = particleI.pauliQ*fdamp3 + dir*fdamp5 + qir*fdamp7;
    double term3 = 2*fdamp5;
    double term4 = 2*(-particleK.pauliQ*fdamp5 + dkr*fdamp7 - qkr*fdamp9);
    double term5 = 2*(-particleI.pauliQ*fdamp5 - dir*fdamp7 - qir*fdamp9);
    double term6 = 4*fdamp7;

    // Compute the force and torque.

    Vec3 f = Vec3(0, 0, de*r) + term1*particleI.qiDipole + term2*particleK.qiDipole + term3*(diqkTemp-dkqiTemp)
            + term4*qi + term5*qk + term6*(qikTemp+qkiTemp);
    f = sizik*(f*rr1 + Vec3(0, 0, eterm*rr3*r));
    Vec3 tI = -fdamp3*dikCross + term1*dirCross + term3*(dqik+dkqirCross) + term4*qirCross - term6*(qikrCross+qikCross);
    Vec3 tK = fdamp3*dikCross + term2*dkrCross - term3*(dqik+diqkrCross) + term5*qkrCross - term6*(qkirCross-qikCross);
    tI *= sizik*rr1;
    tK *= sizik*rr1;
    if (_nonbondedMethod == HippoNonbondedForce::PME && r > _switchingDistance) {
        double t = (r-_switchingDistance)/(_cutoffDistance-_switchingDistance);
        double switchValue = 1+t*t*t*(-10+t*(15-t*6));
        double switchDeriv = t*t*(-30+t*(60-t*30))/(_cutoffDistance-_switchingDistance);
        f *= switchValue;
        f[2] += energy*switchDeriv;
        energy *= switchValue;
        tI *= switchValue;
        tK *= switchValue;
    }
    force += f;
    torqueI += tI;
    torqueK += tK;
    return energy;
}

double AmoebaReferenceHippoNonbondedForce::calculateChargeTransferPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                          double r, Vec3& force) const {
    // Compute the force and energy.

    double term1 = particleI.epsilon*exp(-particleK.damping*r);
    double term2 = particleK.epsilon*exp(-particleI.damping*r);
    double energy = -(term1+term2);
    double dEnergydR = -(term1*particleK.damping + term2*particleI.damping);
    if (_nonbondedMethod == HippoNonbondedForce::PME && r > _switchingDistance) {
        double t = (r-_switchingDistance)/(_cutoffDistance-_switchingDistance);
        double switchValue = 1+t*t*t*(-10+t*(15-t*6));
        double switchDeriv = t*t*(-30+t*(60-t*30))/(_cutoffDistance-_switchingDistance);
        dEnergydR = dEnergydR*switchValue + energy*switchDeriv;
        energy *= switchValue;
    }
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end()) {
        energy *= exception->second.chargeTransferScale;
        dEnergydR *= exception->second.chargeTransferScale;
    }
    force[2] += dEnergydR;
    return energy;
}

void AmoebaReferenceHippoNonbondedForce::mapTorqueToForceForParticle(const MultipoleParticleData& particleI,
                                                                const MultipoleParticleData& particleU,
                                                                const MultipoleParticleData& particleV,
                                                                      MultipoleParticleData* particleW,
                                                                      int axisType, const Vec3& torque,
                                                                      vector<Vec3>& forces) {
    static const int U                  = 0;
    static const int V                  = 1;
    static const int W                  = 2;
    static const int R                  = 3;
    static const int S                  = 4;
    static const int UV                 = 5;
    static const int UW                 = 6;
    static const int VW                 = 7;
    static const int UR                 = 8;
    static const int US                 = 9;
    static const int VS                 = 10;
    static const int WS                 = 11;
    static const int LastVectorIndex    = 12;

    double norms[LastVectorIndex];
    double angles[LastVectorIndex][2];

    // ---------------------------------------------------------------------------------------

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    if (axisType == HippoNonbondedForce::NoAxisType)
        return;

    Vec3 vectorU = particleI.position - particleU.position;
    norms[U] = normalizeVec3(vectorU);
    Vec3 vectorV = particleI.position - particleV.position;
    norms[V] = normalizeVec3(vectorV);
    Vec3 vectorW;
    if (particleW && (axisType == HippoNonbondedForce::ZBisect || axisType == HippoNonbondedForce::ThreeFold))
        vectorW = particleW->position - particleI.position;
    else
        vectorW = vectorU.cross(vectorV);
    norms[W]  = normalizeVec3(vectorW);

    Vec3 vectorUV, vectorUW, vectorVW;
    vectorUV = vectorV.cross(vectorU);
    vectorUW = vectorW.cross(vectorU);
    vectorVW = vectorW.cross(vectorV);

    norms[UV]                     = normalizeVec3(vectorUV);
    norms[UW]                     = normalizeVec3(vectorUW);
    norms[VW]                     = normalizeVec3(vectorVW);

    // angles[][0] is cosine of angle
    // angles[][1] is sine   of angle

    angles[UV][0]                 = vectorU.dot(vectorV);
    angles[UV][1]                 = sqrt(1.0 - angles[UV][0]*angles[UV][0]);

    angles[UW][0]                 = vectorU.dot(vectorW);
    angles[UW][1]                 = sqrt(1.0 - angles[UW][0]*angles[UW][0]);

    angles[VW][0]                 = vectorV.dot(vectorW);
    angles[VW][1]                 = sqrt(1.0 - angles[VW][0]*angles[VW][0]);

    Vec3 dphi;
    dphi[U]                       = vectorU.dot(torque);
    dphi[V]                       = vectorV.dot(torque);
    dphi[W]                       = vectorW.dot(torque);
    dphi                         *= -1.0;

    // branch based on axis type

    if (axisType == HippoNonbondedForce::ZThenX || axisType == HippoNonbondedForce::Bisector) {
        double factor1 =  dphi[V]/(norms[U]*angles[UV][1]);
        double factor2 =  dphi[W]/(norms[U]);
        double factor3 = -dphi[U]/(norms[V]*angles[UV][1]);
        double factor4;
        if (axisType == HippoNonbondedForce::Bisector) {
            factor2 *= 0.5;
            factor4 = 0.5*dphi[W]/(norms[V]);
        }
        else
            factor4 = 0.0;
        Vec3 forceU              =  vectorUV*factor1 + factor2*vectorUW;
        forces[particleU.index] +=  forceU;
        Vec3 forceV              =  vectorUV*factor3 + factor4*vectorVW;
        forces[particleV.index] +=  forceV;
        forces[particleI.index] -=  (forceU + forceV);

    }
    else if (axisType == HippoNonbondedForce::ZBisect) {
        Vec3 vectorR          = vectorV + vectorW;
        Vec3 vectorS          = vectorU.cross(vectorR);

        norms[R]              = normalizeVec3(vectorR);
        norms[S]              = normalizeVec3(vectorS);

        Vec3 vectorUR         =  vectorR.cross(vectorU);
        Vec3 vectorUS         =  vectorS.cross(vectorU);
        Vec3 vectorVS         =  vectorS.cross(vectorV);
        Vec3 vectorWS         =  vectorS.cross(vectorW);

        norms[UR]             = normalizeVec3(vectorUR);
        norms[US]             = normalizeVec3(vectorUS);
        norms[VS]             = normalizeVec3(vectorVS);
        norms[WS]             = normalizeVec3(vectorWS);

        angles[UR][0]         = vectorU.dot(vectorR);
        angles[UR][1]         = sqrt(1.0 - angles[UR][0]*angles[UR][0]);

        angles[US][0]         = vectorU.dot(vectorS);
        angles[US][1]         = sqrt(1.0 - angles[US][0]*angles[US][0]);

        angles[VS][0]         = vectorV.dot(vectorS);
        angles[VS][1]         = sqrt(1.0 - angles[VS][0]*angles[VS][0]);

        angles[WS][0]         = vectorW.dot(vectorS);
        angles[WS][1]         = sqrt(1.0 - angles[WS][0]*angles[WS][0]);

        Vec3 t1               = vectorV - vectorS*angles[VS][0];
        Vec3 t2               = vectorW - vectorS*angles[WS][0];

        double notUsed        = normalizeVec3(t1);
              notUsed         = normalizeVec3(t2);

        double ut1cos         = vectorU.dot(t1);
        double ut1sin         = sqrt(1.0 - ut1cos*ut1cos);

        double ut2cos         = vectorU.dot(t2);
        double ut2sin         = sqrt(1.0 - ut2cos*ut2cos);

        double dphiR          = vectorR.dot(torque)*(-1.0);
        double dphiS          = vectorS.dot(torque)*(-1.0);

        double factor1        = dphiR/(norms[U]*angles[UR][1]);
        double factor2        = dphiS/(norms[U]);
        double factor3        = dphi[U]/(norms[V]*(ut1sin+ut2sin));
        double factor4        = dphi[U]/(norms[W]*(ut1sin+ut2sin));

        Vec3 forceU =  vectorUR*factor1 + vectorUS*factor2;
        forces[particleU.index] += forceU;
        Vec3 forceV = (vectorS*angles[VS][1] - t1*angles[VS][0])*factor3;
        forces[particleV.index] += forceV;
        Vec3 forceW = (vectorS*angles[WS][1] - t2*angles[WS][0])*factor4;
        forces[particleW->index] += forceW;
        forces[particleI.index] -= (forceU + forceV + forceW);

    }
    else if (axisType == HippoNonbondedForce::ThreeFold) {
        Vec3 du =  vectorUW*dphi[W]/(norms[U]*angles[UW][1]) +
                   vectorUV*dphi[V]/(norms[U]*angles[UV][1]) -
                   vectorUW*dphi[U]/(norms[U]*angles[UW][1]) -
                   vectorUV*dphi[U]/(norms[U]*angles[UV][1]);

        Vec3 dv =  vectorVW*dphi[W]/(norms[V]*angles[VW][1]) -
                   vectorUV*dphi[U]/(norms[V]*angles[UV][1]) -
                   vectorVW*dphi[V]/(norms[V]*angles[VW][1]) +
                   vectorUV*dphi[V]/(norms[V]*angles[UV][1]);

        Vec3 dw = -vectorUW*dphi[U]/(norms[W]*angles[UW][1]) -
                   vectorVW*dphi[V]/(norms[W]*angles[VW][1]) +
                   vectorUW*dphi[W]/(norms[W]*angles[UW][1]) +
                   vectorVW*dphi[W]/(norms[W]*angles[VW][1]);

        du /= 3.0;
        dv /= 3.0;
        dw /= 3.0;

        forces[particleU.index] += du;
        forces[particleV.index] += dv;
        if (particleW)
            forces[particleW->index] += dw;
        forces[particleI.index] -= (du + dv + dw);
    }
    else if (axisType == HippoNonbondedForce::ZOnly) {
        Vec3 du = vectorUV*dphi[V]/(norms[U]*angles[UV][1]) + vectorUW*dphi[W]/norms[U];
        forces[particleU.index] += du;
        forces[particleI.index] -= du;
    }
}

void AmoebaReferenceHippoNonbondedForce::mapTorqueToForce(vector<Vec3>& torques,
                                                          vector<Vec3>& forces) {
    // map torques to forces

    for (int ii = 0; ii < _numParticles; ii++) {
        MultipoleParticleData& p = particleData[ii];
        if (p.axisType != HippoNonbondedForce::NoAxisType) {
             mapTorqueToForceForParticle(p,
                                         particleData[p.multipoleAtomZ], particleData[p.multipoleAtomX],
                                         p.multipoleAtomY > -1 ? &particleData[p.multipoleAtomY] : NULL,
                                         p.axisType, torques[ii], forces);
        }
    }
}

double AmoebaReferenceHippoNonbondedForce::calculateInteractions(vector<Vec3>& torques, vector<Vec3>& forces) {

    // main loop over particle pairs

    double energy = 0.0;
    for (int i = 0; i < _numParticles; i++) {
        for (int j = i+1; j < _numParticles; j++) {
            Vec3 deltaR = particleData[j].position - particleData[i].position;
            if (_nonbondedMethod == HippoNonbondedForce::PME)
                getPeriodicDelta(deltaR);
            double r2 = deltaR.dot(deltaR);
            if (_nonbondedMethod == HippoNonbondedForce::PME && r2 > _cutoffDistanceSquared)
                continue;
            double r = sqrt(r2);
            double mat[3][3];
            formQIRotationMatrix(deltaR, r, mat);
            particleData[i].qiDipole = rotateVectorToQI(particleData[i].dipole, mat);
            particleData[j].qiDipole = rotateVectorToQI(particleData[j].dipole, mat);
            particleData[i].qiInducedDipole = rotateVectorToQI(_inducedDipole[i], mat);
            particleData[j].qiInducedDipole = rotateVectorToQI(_inducedDipole[j], mat);
            rotateQuadrupoleToQI(particleData[i].quadrupole, particleData[i].qiQuadrupole, mat);
            rotateQuadrupoleToQI(particleData[j].quadrupole, particleData[j].qiQuadrupole, mat);
            Vec3 force, labForce, torqueI, torqueJ;
            energy += calculateElectrostaticPairIxn(particleData[i], particleData[j], r, force, torqueI, torqueJ);
            calculateInducedDipolePairIxn(particleData[i], particleData[j], deltaR, r, force, torqueI, torqueJ, labForce);
            energy += calculateDispersionPairIxn(particleData[i], particleData[j], r, force);
            energy += calculateRepulsionPairIxn(particleData[i], particleData[j], r, force, torqueI, torqueJ);
            energy += calculateChargeTransferPairIxn(particleData[i], particleData[j], r, force);
            force = rotateVectorFromQI(force, mat);
            torqueI = rotateVectorFromQI(torqueI, mat);
            torqueJ = rotateVectorFromQI(torqueJ, mat);
            forces[i] -= force+labForce;
            forces[j] += force+labForce;
            torques[i] += torqueI;
            torques[j] += torqueJ;
        }
    }
    for (int i = 0; i < _numParticles; i++)
        energy -= (0.5*_electric/particleData[i].polarizability)*_ptDipoleD[0][i].dot(_inducedDipole[i]);
    
    return energy;
}

void AmoebaReferenceHippoNonbondedForce::setup(const vector<Vec3>& particlePositions) {
    loadParticleData(particlePositions);
    applyRotationMatrix();
    calculateInducedDipoles();
}

double AmoebaReferenceHippoNonbondedForce::calculateForceAndEnergy(const vector<Vec3>& particlePositions,
                                                                   vector<Vec3>& forces) {
    // setup, including calculating induced dipoles
    // calculate electrostatic ixns including torques
    // map torques to forces

    setup(particlePositions);

    vector<Vec3> torques;
    initializeVec3Vector(torques);
    double energy = calculateInteractions(torques, forces);

    mapTorqueToForce(torques, forces);

    return energy;
}

void AmoebaReferenceHippoNonbondedForce::calculateInducedDipoles(const vector<Vec3>& particlePositions,
                                                                 vector<Vec3>& outputInducedDipoles) {
    // setup, including calculating induced dipoles

    setup(particlePositions);
    outputInducedDipoles = _inducedDipole;
}




void AmoebaReferenceHippoNonbondedForce::calculateLabFramePermanentDipoles(const vector<Vec3>& particlePositions,
                                                                      vector<Vec3>& outputRotatedPermanentDipoles) {
    // setup, including calculating permanent dipoles

    setup(particlePositions);
    outputRotatedPermanentDipoles.resize(_numParticles);
    for (int i = 0; i < _numParticles; i++)
        outputRotatedPermanentDipoles[i] = particleData[i].dipole;
}

const int AmoebaReferencePmeHippoNonbondedForce::AMOEBA_PME_ORDER = 5;

const double AmoebaReferencePmeHippoNonbondedForce::SQRT_PI = sqrt(M_PI);

AmoebaReferencePmeHippoNonbondedForce::AmoebaReferencePmeHippoNonbondedForce(const HippoNonbondedForce& force, const System& system) :
               AmoebaReferenceHippoNonbondedForce(force) {
    _fftplan = NULL;
    force.getPMEParameters(_alphaEwald, _pmeGridDimensions[0], _pmeGridDimensions[1], _pmeGridDimensions[2]);
    force.getDPMEParameters(_dalphaEwald, _dpmeGridDimensions[0], _dpmeGridDimensions[1], _dpmeGridDimensions[2]);
    if (_alphaEwald == 0.0 || _dalphaEwald == 0.0) {
        NonbondedForce nb;
        nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
        nb.setCutoffDistance(force.getCutoffDistance());
        if (_alphaEwald == 0.0)
            NonbondedForceImpl::calcPMEParameters(system, nb, _alphaEwald, _pmeGridDimensions[0], _pmeGridDimensions[1], _pmeGridDimensions[2], false);
        if (_dalphaEwald == 0.0)
            NonbondedForceImpl::calcPMEParameters(system, nb, _dalphaEwald, _dpmeGridDimensions[0], _dpmeGridDimensions[1], _dpmeGridDimensions[2], true);
    }
    fftpack_init_3d(&_fftplan, _pmeGridDimensions[0], _pmeGridDimensions[1], _pmeGridDimensions[2]);
    initializeBSplineModuli();
}

AmoebaReferencePmeHippoNonbondedForce::~AmoebaReferencePmeHippoNonbondedForce() {
    if (_fftplan != NULL)
        fftpack_destroy(_fftplan);
};

double AmoebaReferencePmeHippoNonbondedForce::getCutoffDistance() const {
     return _cutoffDistance;
};

double AmoebaReferencePmeHippoNonbondedForce::getAlphaEwald() const {
     return _alphaEwald;
};

double AmoebaReferencePmeHippoNonbondedForce::getDispersionAlphaEwald() const {
     return _dalphaEwald;
};

void AmoebaReferencePmeHippoNonbondedForce::getPmeGridDimensions(vector<int>& pmeGridDimensions) const {
    pmeGridDimensions.resize(3);

    pmeGridDimensions[0] = _pmeGridDimensions[0];
    pmeGridDimensions[1] = _pmeGridDimensions[1];
    pmeGridDimensions[2] = _pmeGridDimensions[2];
};

void AmoebaReferencePmeHippoNonbondedForce::getDispersionPmeGridDimensions(vector<int>& pmeGridDimensions) const {
    pmeGridDimensions.resize(3);

    pmeGridDimensions[0] = _dpmeGridDimensions[0];
    pmeGridDimensions[1] = _dpmeGridDimensions[1];
    pmeGridDimensions[2] = _dpmeGridDimensions[2];
};

void AmoebaReferencePmeHippoNonbondedForce::setPmeGridDimensions(vector<int>& pmeGridDimensions) {
    if ((pmeGridDimensions[0] == _pmeGridDimensions[0]) &&
        (pmeGridDimensions[1] == _pmeGridDimensions[1]) &&
        (pmeGridDimensions[2] == _pmeGridDimensions[2]))
        return;

    if (_fftplan) {
        fftpack_destroy(_fftplan);
    }
    fftpack_init_3d(&_fftplan,pmeGridDimensions[0], pmeGridDimensions[1], pmeGridDimensions[2]);

    _pmeGridDimensions[0] = pmeGridDimensions[0];
    _pmeGridDimensions[1] = pmeGridDimensions[1];
    _pmeGridDimensions[2] = pmeGridDimensions[2];

    initializeBSplineModuli();
};

void AmoebaReferencePmeHippoNonbondedForce::setDispersionPmeGridDimensions(vector<int>& pmeGridDimensions) {
    if ((pmeGridDimensions[0] == _dpmeGridDimensions[0]) &&
        (pmeGridDimensions[1] == _dpmeGridDimensions[1]) &&
        (pmeGridDimensions[2] == _dpmeGridDimensions[2]))
        return;

    if (_fftplan) {
        fftpack_destroy(_fftplan);
    }
    fftpack_init_3d(&_fftplan,pmeGridDimensions[0], pmeGridDimensions[1], pmeGridDimensions[2]);

    _dpmeGridDimensions[0] = pmeGridDimensions[0];
    _dpmeGridDimensions[1] = pmeGridDimensions[1];
    _dpmeGridDimensions[2] = pmeGridDimensions[2];

    initializeBSplineModuli();
};

void AmoebaReferencePmeHippoNonbondedForce::setPeriodicBoxSize(OpenMM::Vec3* vectors) {
    if (vectors[0][0] == 0.0 || vectors[1][1] == 0.0 || vectors[2][2] == 0.0) {
        std::stringstream message;
        message << "Box size of zero is invalid.";
        throw OpenMMException(message.str());
    }

    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
    double determinant = vectors[0][0]*vectors[1][1]*vectors[2][2];
    assert(determinant > 0);
    double scale = 1.0/determinant;
    _recipBoxVectors[0] = Vec3(vectors[1][1]*vectors[2][2], 0, 0)*scale;
    _recipBoxVectors[1] = Vec3(-vectors[1][0]*vectors[2][2], vectors[0][0]*vectors[2][2], 0)*scale;
    _recipBoxVectors[2] = Vec3(vectors[1][0]*vectors[2][1]-vectors[1][1]*vectors[2][0], -vectors[0][0]*vectors[2][1], vectors[0][0]*vectors[1][1])*scale;
};

void AmoebaReferencePmeHippoNonbondedForce::resizePmeArrays() {
    int gridSize = _pmeGridDimensions[0]*_pmeGridDimensions[1]*_pmeGridDimensions[2];
    _pmeGrid.resize(gridSize);

    for (int ii = 0; ii < 3; ii++) {
       _pmeBsplineModuli[ii].resize(_pmeGridDimensions[ii]);
       _thetai[ii].resize(AMOEBA_PME_ORDER*_numParticles);
    }

    _iGrid.resize(_numParticles);
    _phi.resize(20*_numParticles);
    _phidp.resize(20*_numParticles);
    optPhi.resize(_maxPTOrder, vector<double>(10*_numParticles));
}

void AmoebaReferencePmeHippoNonbondedForce::initializePmeGrid() {
    for (int jj = 0; jj < _pmeGrid.size(); jj++)
        _pmeGrid[jj].re = _pmeGrid[jj].im = 0.0;
}

void AmoebaReferencePmeHippoNonbondedForce::getPeriodicDelta(Vec3& deltaR) const {
    deltaR -= _periodicBoxVectors[2]*floor(deltaR[2]*_recipBoxVectors[2][2]+0.5);
    deltaR -= _periodicBoxVectors[1]*floor(deltaR[1]*_recipBoxVectors[1][1]+0.5);
    deltaR -= _periodicBoxVectors[0]*floor(deltaR[0]*_recipBoxVectors[0][0]+0.5);
}

void AmoebaReferencePmeHippoNonbondedForce::initializeBSplineModuli() {
    // Initialize the b-spline moduli.

    int maxSize = -1;
    for (int ii = 0; ii < 3; ii++) {
       _pmeBsplineModuli[ii].resize(_pmeGridDimensions[ii]);
        maxSize = maxSize  > _pmeGridDimensions[ii] ? maxSize : _pmeGridDimensions[ii];
    }

    double array[AMOEBA_PME_ORDER];
    double x = 0.0;
    array[0]     = 1.0 - x;
    array[1]     = x;
    for (int k = 2; k < AMOEBA_PME_ORDER; k++) {
        double denom = 1.0/k;
        array[k] = x*array[k-1]*denom;
        for (int i = 1; i < k; i++) {
            array[k-i] = ((x+i)*array[k-i-1] + ((k-i+1)-x)*array[k-i])*denom;
        }
        array[0] = (1.0-x)*array[0]*denom;
    }

    vector<double> bsarray(maxSize+1, 0.0);
    for (int i = 2; i <= AMOEBA_PME_ORDER+1; i++) {
        bsarray[i] = array[i-2];
    }
    for (int dim = 0; dim < 3; dim++) {

        int size = _pmeGridDimensions[dim];

        // get the modulus of the discrete Fourier transform

        double factor = 2.0*M_PI/size;
        for (int i = 0; i < size; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int j = 1; j <= size; j++) {
                double arg = factor*i*(j-1);
                sum1 += bsarray[j]*cos(arg);
                sum2 += bsarray[j]*sin(arg);
            }
            _pmeBsplineModuli[dim][i] = (sum1*sum1 + sum2*sum2);
        }

        // fix for exponential Euler spline interpolation failure

        double eps = 1.0e-7;
        if (_pmeBsplineModuli[dim][0] < eps) {
            _pmeBsplineModuli[dim][0] = 0.5*_pmeBsplineModuli[dim][1];
        }
        for (int i = 1; i < size-1; i++) {
            if (_pmeBsplineModuli[dim][i] < eps) {
                _pmeBsplineModuli[dim][i] = 0.5*(_pmeBsplineModuli[dim][i-1]+_pmeBsplineModuli[dim][i+1]);
            }
        }
        if (_pmeBsplineModuli[dim][size-1] < eps) {
            _pmeBsplineModuli[dim][size-1] = 0.5*_pmeBsplineModuli[dim][size-2];
        }

        // compute and apply the optimal zeta coefficient

        int jcut = 50;
        for (int i = 1; i <= size; i++) {
            int k = i - 1;
            if (i > size/2)
                k = k - size;
            double zeta;
            if (k == 0)
                zeta = 1.0;
            else {
                double sum1 = 1.0;
                double sum2 = 1.0;
                factor = M_PI*k/size;
                for (int j = 1; j <= jcut; j++) {
                    double arg = factor/(factor+M_PI*j);
                    sum1 = sum1 + pow(arg,   AMOEBA_PME_ORDER);
                    sum2 = sum2 + pow(arg, 2*AMOEBA_PME_ORDER);
                }
                for (int j = 1; j <= jcut; j++) {
                    double arg  = factor/(factor-M_PI*j);
                    sum1 += pow(arg,   AMOEBA_PME_ORDER);
                    sum2 += pow(arg, 2*AMOEBA_PME_ORDER);
                }
                zeta = sum2/sum1;
            }
            _pmeBsplineModuli[dim][i-1] = _pmeBsplineModuli[dim][i-1]*(zeta*zeta);
        }
    }
}

void AmoebaReferencePmeHippoNonbondedForce::calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI,
                                                                                const MultipoleParticleData& particleJ) {
    // compute the real space portion of the Ewald summation

    Vec3 deltaR = particleJ.position - particleI.position;
    getPeriodicDelta(deltaR);
    double r2 = deltaR.dot(deltaR);
    if (r2 > _cutoffDistanceSquared)
        return;
    double r = sqrt(r2);
    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rInv3 = rInv*rInv2;
    double rInv5 = rInv3*rInv2;
    double rInv7 = rInv5*rInv2;

    // calculate the error function damping terms

    double ralpha = _alphaEwald*r;
    double bn0 = erfc(ralpha)*rInv;
    double alsq2 = 2*_alphaEwald*_alphaEwald;
    double alsq2n = 1/(SQRT_PI*_alphaEwald);
    double exp2a = exp(-(ralpha*ralpha));
    alsq2n *= alsq2;
    double bn1 = (bn0+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn2 = (3*bn1+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn3 = (5*bn2+alsq2n*exp2a)*rInv2;

    // Calculate the field at particle I due to multipoles at particle J.

    double fdamp3, fdamp5, fdamp7;
    computeDirectFieldDampingFactors(particleJ, r, fdamp3, fdamp5, fdamp7);
    double scale = 1;
    auto exception = exceptions.find(make_pair(particleI.index, particleJ.index));
    if (exception != exceptions.end())
        scale = exception->second.dipoleMultipoleScale;
    double rr3 = bn1 - (1-scale)*rInv3;
    double rr3j = bn1 - (1-scale*fdamp3)*rInv3;
    double rr5j = bn2 - (1-scale*fdamp5)*3*rInv5;
    double rr7j = bn3 - (1-scale*fdamp7)*15*rInv7;
    Vec3 qDotDelta(deltaR[0]*particleJ.quadrupole[QXX] + deltaR[1]*particleJ.quadrupole[QXY] + deltaR[2]*particleJ.quadrupole[QXZ],
                   deltaR[0]*particleJ.quadrupole[QXY] + deltaR[1]*particleJ.quadrupole[QYY] + deltaR[2]*particleJ.quadrupole[QYZ],
                   deltaR[0]*particleJ.quadrupole[QXZ] + deltaR[1]*particleJ.quadrupole[QYZ] + deltaR[2]*particleJ.quadrupole[QZZ]);
    double dipoleDelta = particleJ.dipole.dot(deltaR);
    double qdpoleDelta = qDotDelta.dot(deltaR);
    double factor = rr3*particleJ.coreCharge + rr3j*particleJ.valenceCharge - rr5j*dipoleDelta + rr7j*qdpoleDelta;
    Vec3 field = deltaR*factor + particleJ.dipole*rr3j - qDotDelta*2*rr5j;
    _fixedMultipoleField[particleI.index] -= field;
}

void AmoebaReferencePmeHippoNonbondedForce::calculateFixedMultipoleField() {
    // first calculate reciprocal space fixed multipole fields

    resizePmeArrays();
    computeAmoebaBsplines(particleData);
    initializePmeGrid();
    spreadFixedMultipolesOntoGrid(particleData);
    fftpack_exec_3d(_fftplan, FFTPACK_FORWARD, _pmeGrid.data(), _pmeGrid.data());
    performAmoebaReciprocalConvolution();
    fftpack_exec_3d(_fftplan, FFTPACK_BACKWARD, _pmeGrid.data(), _pmeGrid.data());
    computeFixedPotentialFromGrid();
    recordFixedMultipoleField();

    // include self-energy portion of the multipole field

    double term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (int j = 0; j < _numParticles; j++)
        _fixedMultipoleField[j] += particleData[j].dipole*term;

    // include direct space fixed multipole fields

    AmoebaReferenceHippoNonbondedForce::calculateFixedMultipoleField();
}

#define ARRAY(x,y) array[(x)-1+((y)-1)*AMOEBA_PME_ORDER]

/**
 * This is called from computeBsplines().  It calculates the spline coefficients for a single atom along a single axis.
 */
void AmoebaReferencePmeHippoNonbondedForce::computeBSplinePoint(vector<HippoDouble4>& thetai, double w) {
    double array[AMOEBA_PME_ORDER*AMOEBA_PME_ORDER];

    // initialization to get to 2nd order recursion

    ARRAY(2,2) = w;
    ARRAY(2,1) = 1.0 - w;

    // perform one pass to get to 3rd order recursion

    ARRAY(3,3) = 0.5 * w * ARRAY(2,2);
    ARRAY(3,2) = 0.5 * ((1.0+w)*ARRAY(2,1)+(2.0-w)*ARRAY(2,2));
    ARRAY(3,1) = 0.5 * (1.0-w) * ARRAY(2,1);

    // compute standard B-spline recursion to desired order

    for (int i = 4; i <= AMOEBA_PME_ORDER; i++) {
        int k = i - 1;
        double denom = 1.0 / k;
        ARRAY(i,i) = denom * w * ARRAY(k,k);
        for (int j = 1; j <= i-2; j++)
            ARRAY(i,i-j) = denom * ((w+j)*ARRAY(k,i-j-1)+(i-j-w)*ARRAY(k,i-j));
        ARRAY(i,1) = denom * (1.0-w) * ARRAY(k,1);
    }

    // get coefficients for the B-spline first derivative

    int k = AMOEBA_PME_ORDER - 1;
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline second derivative

    k = AMOEBA_PME_ORDER - 2;
    ARRAY(k,AMOEBA_PME_ORDER-1) = ARRAY(k,AMOEBA_PME_ORDER-2);
    for (int i = AMOEBA_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline third derivative

    k = AMOEBA_PME_ORDER - 3;
    ARRAY(k,AMOEBA_PME_ORDER-2) = ARRAY(k,AMOEBA_PME_ORDER-3);
    for (int i = AMOEBA_PME_ORDER-3; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER-1) = ARRAY(k,AMOEBA_PME_ORDER-2);
    for (int i = AMOEBA_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,AMOEBA_PME_ORDER) = ARRAY(k,AMOEBA_PME_ORDER-1);
    for (int i = AMOEBA_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // copy coefficients from temporary to permanent storage

    for (int i = 1; i <= AMOEBA_PME_ORDER; i++)
        thetai[i-1] = HippoDouble4(ARRAY(AMOEBA_PME_ORDER,i), ARRAY(AMOEBA_PME_ORDER-1,i), ARRAY(AMOEBA_PME_ORDER-2,i), ARRAY(AMOEBA_PME_ORDER-3,i));
}

/**
 * Compute b-spline coefficients.
 */
void AmoebaReferencePmeHippoNonbondedForce::computeAmoebaBsplines(const vector<MultipoleParticleData>& particleData) {
    //  get the B-spline coefficients for each multipole site

    for (int ii = 0; ii < _numParticles; ii++) {
        Vec3 position  = particleData[ii].position;
        getPeriodicDelta(position);
        int igrid[3];
        for (int jj = 0; jj < 3; jj++) {

            double w  = position[0]*_recipBoxVectors[0][jj]+position[1]*_recipBoxVectors[1][jj]+position[2]*_recipBoxVectors[2][jj];
            double fr = _pmeGridDimensions[jj]*(w-(int)(w+0.5)+0.5);
            int ifr   = static_cast<int>(floor(fr));
            w         = fr - ifr;
            igrid[jj] = ifr - AMOEBA_PME_ORDER + 1;
            igrid[jj] += igrid[jj] < 0 ? _pmeGridDimensions[jj] : 0;
            vector<HippoDouble4> thetaiTemp(AMOEBA_PME_ORDER);
            computeBSplinePoint(thetaiTemp, w);
            for (int kk = 0; kk < AMOEBA_PME_ORDER; kk++)
                _thetai[jj][ii*AMOEBA_PME_ORDER+kk] = thetaiTemp[kk];
        }

        // Record the grid point.

        _iGrid[ii] = {igrid[0], igrid[1], igrid[2]};
    }
}

void AmoebaReferencePmeHippoNonbondedForce::transformMultipolesToFractionalCoordinates(const vector<MultipoleParticleData>& particleData) {
    // Build matrices for transforming the dipoles and quadrupoles.

    Vec3 a[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a[j][i] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    int index1[] = {0, 0, 0, 1, 1, 2};
    int index2[] = {0, 1, 2, 1, 2, 2};
    double b[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            b[i][j] = a[index1[i]][index1[j]]*a[index2[i]][index2[j]];
            if (index1[i] != index2[i])
                b[i][j] += a[index1[i]][index2[j]]*a[index2[i]][index1[j]];
        }
    }

    // Transform the multipoles.

    _transformed.resize(_numParticles);
    double quadScale[] = {1, 2, 2, 1, 2, 1};
    for (int i = 0; i < _numParticles; i++) {
        _transformed[i].charge = particleData[i].coreCharge + particleData[i].valenceCharge;
        _transformed[i].dipole = Vec3();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                _transformed[i].dipole[j] += a[j][k]*particleData[i].dipole[k];
        for (int j = 0; j < 6; j++) {
            _transformed[i].quadrupole[j] = 0;
            for (int k = 0; k < 6; k++)
                _transformed[i].quadrupole[j] += quadScale[k]*b[j][k]*particleData[i].quadrupole[k];
        }
    }
}

void AmoebaReferencePmeHippoNonbondedForce::transformPotentialToCartesianCoordinates(const vector<double>& fphi, vector<double>& cphi) const {
    // Build a matrix for transforming the potential.

    Vec3 a[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    int index1[] = {0, 1, 2, 0, 0, 1};
    int index2[] = {0, 1, 2, 1, 2, 2};
    double b[6][6];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            b[i][j] = a[index1[i]][index1[j]]*a[index2[i]][index2[j]];
            if (index1[j] != index2[j])
                b[i][j] *= 2;
        }
    }
    for (int i = 3; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            b[i][j] = a[index1[i]][index1[j]]*a[index2[i]][index2[j]];
            if (index1[j] != index2[j])
                b[i][j] += a[index1[i]][index2[j]]*a[index2[i]][index1[j]];
        }
    }

    // Transform the potential.

    for (int i = 0; i < _numParticles; i++) {
        cphi[10*i] = fphi[20*i];
        cphi[10*i+1] = a[0][0]*fphi[20*i+1] + a[0][1]*fphi[20*i+2] + a[0][2]*fphi[20*i+3];
        cphi[10*i+2] = a[1][0]*fphi[20*i+1] + a[1][1]*fphi[20*i+2] + a[1][2]*fphi[20*i+3];
        cphi[10*i+3] = a[2][0]*fphi[20*i+1] + a[2][1]*fphi[20*i+2] + a[2][2]*fphi[20*i+3];
        for (int j = 0; j < 6; j++) {
            cphi[10*i+4+j] = 0;
            for (int k = 0; k < 6; k++)
                cphi[10*i+4+j] += b[j][k]*fphi[20*i+4+k];
        }
    }
}

void AmoebaReferencePmeHippoNonbondedForce::spreadFixedMultipolesOntoGrid(const vector<MultipoleParticleData>& particleData) {
    transformMultipolesToFractionalCoordinates(particleData);

    // Clear the grid.

    for (int gridIndex = 0; gridIndex < _pmeGrid.size(); gridIndex++)
        _pmeGrid[gridIndex] = t_complex(0, 0);

    // Loop over atoms and spread them on the grid.

    for (int atomIndex = 0; atomIndex < _numParticles; atomIndex++) {
        double atomCharge = _transformed[atomIndex].charge;
        Vec3 atomDipole = Vec3(_transformed[atomIndex].dipole[0],
                               _transformed[atomIndex].dipole[1],
                               _transformed[atomIndex].dipole[2]);

        double atomQuadrupoleXX = _transformed[atomIndex].quadrupole[QXX];
        double atomQuadrupoleXY = _transformed[atomIndex].quadrupole[QXY];
        double atomQuadrupoleXZ = _transformed[atomIndex].quadrupole[QXZ];
        double atomQuadrupoleYY = _transformed[atomIndex].quadrupole[QYY];
        double atomQuadrupoleYZ = _transformed[atomIndex].quadrupole[QYZ];
        double atomQuadrupoleZZ = _transformed[atomIndex].quadrupole[QZZ];
        array<int,3>& gridPoint = _iGrid[atomIndex];
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
            int x = (gridPoint[0]+ix) % _pmeGridDimensions[0];
            HippoDouble4 t = _thetai[0][atomIndex*AMOEBA_PME_ORDER+ix];
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int y = (gridPoint[1]+iy) % _pmeGridDimensions[1];
                HippoDouble4 u = _thetai[1][atomIndex*AMOEBA_PME_ORDER+iy];
                double term0 = atomCharge*t[0]*u[0] + atomDipole[1]*t[0]*u[1] + atomQuadrupoleYY*t[0]*u[2] + atomDipole[0]*t[1]*u[0] + atomQuadrupoleXY*t[1]*u[1] + atomQuadrupoleXX*t[2]*u[0];
                double term1 = atomDipole[2]*t[0]*u[0] + atomQuadrupoleYZ*t[0]*u[1] + atomQuadrupoleXZ*t[1]*u[0];
                double term2 = atomQuadrupoleZZ*t[0]*u[0];
                for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
                    int z = (gridPoint[2]+iz) % _pmeGridDimensions[2];
                    HippoDouble4 v = _thetai[2][atomIndex*AMOEBA_PME_ORDER+iz];
                    t_complex& gridValue = _pmeGrid[x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z];
                    gridValue.re += term0*v[0] + term1*v[1] + term2*v[2];
                }
            }
        }
    }
}

void AmoebaReferencePmeHippoNonbondedForce::performAmoebaReciprocalConvolution() {
    double expFactor = (M_PI*M_PI)/(_alphaEwald*_alphaEwald);
    double scaleFactor = 1.0/(M_PI*_periodicBoxVectors[0][0]*_periodicBoxVectors[1][1]*_periodicBoxVectors[2][2]);

    for (int index = 0; index < _pmeGrid.size(); index++) {
        int kx = index/(_pmeGridDimensions[1]*_pmeGridDimensions[2]);
        int remainder = index-kx*_pmeGridDimensions[1]*_pmeGridDimensions[2];
        int ky = remainder/_pmeGridDimensions[2];
        int kz = remainder-ky*_pmeGridDimensions[2];

        if (kx == 0 && ky == 0 && kz == 0) {
            _pmeGrid[index].re = _pmeGrid[index].im = 0.0;
            continue;
        }

        int mx = (kx < (_pmeGridDimensions[0]+1)/2) ? kx : (kx-_pmeGridDimensions[0]);
        int my = (ky < (_pmeGridDimensions[1]+1)/2) ? ky : (ky-_pmeGridDimensions[1]);
        int mz = (kz < (_pmeGridDimensions[2]+1)/2) ? kz : (kz-_pmeGridDimensions[2]);

        double mhx = mx*_recipBoxVectors[0][0];
        double mhy = mx*_recipBoxVectors[1][0]+my*_recipBoxVectors[1][1];
        double mhz = mx*_recipBoxVectors[2][0]+my*_recipBoxVectors[2][1]+mz*_recipBoxVectors[2][2];

        double bx = _pmeBsplineModuli[0][kx];
        double by = _pmeBsplineModuli[1][ky];
        double bz = _pmeBsplineModuli[2][kz];

        double m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        double denom = m2*bx*by*bz;
        double eterm = scaleFactor*exp(-expFactor*m2)/denom;

        _pmeGrid[index].re *= eterm;
        _pmeGrid[index].im *= eterm;
    }
}

void AmoebaReferencePmeHippoNonbondedForce::computeFixedPotentialFromGrid() {
    // extract the permanent multipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        array<int,3>& gridPoint = _iGrid[m];
        double tuv000 = 0.0;
        double tuv001 = 0.0;
        double tuv010 = 0.0;
        double tuv100 = 0.0;
        double tuv200 = 0.0;
        double tuv020 = 0.0;
        double tuv002 = 0.0;
        double tuv110 = 0.0;
        double tuv101 = 0.0;
        double tuv011 = 0.0;
        double tuv300 = 0.0;
        double tuv030 = 0.0;
        double tuv003 = 0.0;
        double tuv210 = 0.0;
        double tuv201 = 0.0;
        double tuv120 = 0.0;
        double tuv021 = 0.0;
        double tuv102 = 0.0;
        double tuv012 = 0.0;
        double tuv111 = 0.0;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            HippoDouble4 v = _thetai[2][m*AMOEBA_PME_ORDER+iz];
            double tu00 = 0.0;
            double tu10 = 0.0;
            double tu01 = 0.0;
            double tu20 = 0.0;
            double tu11 = 0.0;
            double tu02 = 0.0;
            double tu30 = 0.0;
            double tu21 = 0.0;
            double tu12 = 0.0;
            double tu03 = 0.0;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                HippoDouble4 u = _thetai[1][m*AMOEBA_PME_ORDER+iy];
                HippoDouble4 t = HippoDouble4(0.0, 0.0, 0.0, 0.0);
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    double tq = _pmeGrid[gridIndex].re;
                    HippoDouble4 tadd = _thetai[0][m*AMOEBA_PME_ORDER+ix];
                    t[0] += tq*tadd[0];
                    t[1] += tq*tadd[1];
                    t[2] += tq*tadd[2];
                    t[3] += tq*tadd[3];
                }
                tu00 += t[0]*u[0];
                tu10 += t[1]*u[0];
                tu01 += t[0]*u[1];
                tu20 += t[2]*u[0];
                tu11 += t[1]*u[1];
                tu02 += t[0]*u[2];
                tu30 += t[3]*u[0];
                tu21 += t[2]*u[1];
                tu12 += t[1]*u[2];
                tu03 += t[0]*u[3];
            }
            tuv000 += tu00*v[0];
            tuv100 += tu10*v[0];
            tuv010 += tu01*v[0];
            tuv001 += tu00*v[1];
            tuv200 += tu20*v[0];
            tuv020 += tu02*v[0];
            tuv002 += tu00*v[2];
            tuv110 += tu11*v[0];
            tuv101 += tu10*v[1];
            tuv011 += tu01*v[1];
            tuv300 += tu30*v[0];
            tuv030 += tu03*v[0];
            tuv003 += tu00*v[3];
            tuv210 += tu21*v[0];
            tuv201 += tu20*v[1];
            tuv120 += tu12*v[0];
            tuv021 += tu02*v[1];
            tuv102 += tu10*v[2];
            tuv012 += tu01*v[2];
            tuv111 += tu11*v[1];
        }
        _phi[20*m] = tuv000;
        _phi[20*m+1] = tuv100;
        _phi[20*m+2] = tuv010;
        _phi[20*m+3] = tuv001;
        _phi[20*m+4] = tuv200;
        _phi[20*m+5] = tuv020;
        _phi[20*m+6] = tuv002;
        _phi[20*m+7] = tuv110;
        _phi[20*m+8] = tuv101;
        _phi[20*m+9] = tuv011;
        _phi[20*m+10] = tuv300;
        _phi[20*m+11] = tuv030;
        _phi[20*m+12] = tuv003;
        _phi[20*m+13] = tuv210;
        _phi[20*m+14] = tuv201;
        _phi[20*m+15] = tuv120;
        _phi[20*m+16] = tuv021;
        _phi[20*m+17] = tuv102;
        _phi[20*m+18] = tuv012;
        _phi[20*m+19] = tuv111;
    }
}

void AmoebaReferencePmeHippoNonbondedForce::spreadInducedDipolesOnGrid(const vector<Vec3>& inputInducedDipole) {
    // Create the matrix to convert from Cartesian to fractional coordinates.

    Vec3 cartToFrac[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cartToFrac[j][i] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];

    // Clear the grid.

    for (int gridIndex = 0; gridIndex < _pmeGrid.size(); gridIndex++)
        _pmeGrid[gridIndex] = t_complex(0, 0);

    // Loop over atoms and spread them on the grid.

    for (int atomIndex = 0; atomIndex < _numParticles; atomIndex++) {
        Vec3 inducedDipole = Vec3(inputInducedDipole[atomIndex][0]*cartToFrac[0][0] + inputInducedDipole[atomIndex][1]*cartToFrac[0][1] + inputInducedDipole[atomIndex][2]*cartToFrac[0][2],
                                  inputInducedDipole[atomIndex][0]*cartToFrac[1][0] + inputInducedDipole[atomIndex][1]*cartToFrac[1][1] + inputInducedDipole[atomIndex][2]*cartToFrac[1][2],
                                  inputInducedDipole[atomIndex][0]*cartToFrac[2][0] + inputInducedDipole[atomIndex][1]*cartToFrac[2][1] + inputInducedDipole[atomIndex][2]*cartToFrac[2][2]);
        array<int,3>& gridPoint = _iGrid[atomIndex];
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
            int x = (gridPoint[0]+ix) % _pmeGridDimensions[0];
            HippoDouble4 t = _thetai[0][atomIndex*AMOEBA_PME_ORDER+ix];
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int y = (gridPoint[1]+iy) % _pmeGridDimensions[1];
                HippoDouble4 u = _thetai[1][atomIndex*AMOEBA_PME_ORDER+iy];
                double term01 = inducedDipole[1]*t[0]*u[1] + inducedDipole[0]*t[1]*u[0];
                double term11 = inducedDipole[2]*t[0]*u[0];
                for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
                    int z = (gridPoint[2]+iz) % _pmeGridDimensions[2];
                    HippoDouble4 v = _thetai[2][atomIndex*AMOEBA_PME_ORDER+iz];
                    t_complex& gridValue = _pmeGrid[x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z];
                    gridValue.re += term01*v[0] + term11*v[1];
                }
            }
        }
    }
}

void AmoebaReferencePmeHippoNonbondedForce::computeInducedPotentialFromGrid() {
    // extract the induced dipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        array<int,3>& gridPoint = _iGrid[m];
        double tuv000 = 0.0;
        double tuv001 = 0.0;
        double tuv010 = 0.0;
        double tuv100 = 0.0;
        double tuv200 = 0.0;
        double tuv020 = 0.0;
        double tuv002 = 0.0;
        double tuv110 = 0.0;
        double tuv101 = 0.0;
        double tuv011 = 0.0;
        double tuv300 = 0.0;
        double tuv030 = 0.0;
        double tuv003 = 0.0;
        double tuv210 = 0.0;
        double tuv201 = 0.0;
        double tuv120 = 0.0;
        double tuv021 = 0.0;
        double tuv102 = 0.0;
        double tuv012 = 0.0;
        double tuv111 = 0.0;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            HippoDouble4 v = _thetai[2][m*AMOEBA_PME_ORDER+iz];
            double tu00 = 0.0;
            double tu10 = 0.0;
            double tu01 = 0.0;
            double tu20 = 0.0;
            double tu11 = 0.0;
            double tu02 = 0.0;
            double tu30 = 0.0;
            double tu21 = 0.0;
            double tu12 = 0.0;
            double tu03 = 0.0;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                HippoDouble4 u = _thetai[1][m*AMOEBA_PME_ORDER+iy];
                double t0 = 0.0;
                double t1 = 0.0;
                double t2 = 0.0;
                double t3 = 0.0;
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    t_complex tq = _pmeGrid[gridIndex];
                    HippoDouble4 tadd = _thetai[0][m*AMOEBA_PME_ORDER+ix];
                    t0 += tq.re*tadd[0];
                    t1 += tq.re*tadd[1];
                    t2 += tq.re*tadd[2];
                    t3 += tq.re*tadd[3];
                }
                tu00 += t0*u[0];
                tu10 += t1*u[0];
                tu01 += t0*u[1];
                tu20 += t2*u[0];
                tu11 += t1*u[1];
                tu02 += t0*u[2];
                tu30 += t3*u[0];
                tu21 += t2*u[1];
                tu12 += t1*u[2];
                tu03 += t0*u[3];
            }
            tuv000 += tu00*v[0];
            tuv100 += tu10*v[0];
            tuv010 += tu01*v[0];
            tuv001 += tu00*v[1];
            tuv200 += tu20*v[0];
            tuv020 += tu02*v[0];
            tuv002 += tu00*v[2];
            tuv110 += tu11*v[0];
            tuv101 += tu10*v[1];
            tuv011 += tu01*v[1];
            tuv300 += tu30*v[0];
            tuv030 += tu03*v[0];
            tuv003 += tu00*v[3];
            tuv210 += tu21*v[0];
            tuv201 += tu20*v[1];
            tuv120 += tu12*v[0];
            tuv021 += tu02*v[1];
            tuv102 += tu10*v[2];
            tuv012 += tu01*v[2];
            tuv111 += tu11*v[1];
        }
        _phidp[20*m] = tuv000;
        _phidp[20*m+1] = tuv100;
        _phidp[20*m+2] = tuv010;
        _phidp[20*m+3] = tuv001;
        _phidp[20*m+4] = tuv200;
        _phidp[20*m+5] = tuv020;
        _phidp[20*m+6] = tuv002;
        _phidp[20*m+7] = tuv110;
        _phidp[20*m+8] = tuv101;
        _phidp[20*m+9] = tuv011;
        _phidp[20*m+10] = tuv300;
        _phidp[20*m+11] = tuv030;
        _phidp[20*m+12] = tuv003;
        _phidp[20*m+13] = tuv210;
        _phidp[20*m+14] = tuv201;
        _phidp[20*m+15] = tuv120;
        _phidp[20*m+16] = tuv021;
        _phidp[20*m+17] = tuv102;
        _phidp[20*m+18] = tuv012;
        _phidp[20*m+19] = tuv111;
    }
}

double AmoebaReferencePmeHippoNonbondedForce::computeReciprocalSpaceFixedMultipoleForceAndEnergy(const vector<MultipoleParticleData>& particleData,
                                                                                                 vector<Vec3>& forces, vector<Vec3>& torques) const {
    double multipole[10];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    vector<double> cphi(10*_numParticles);
    transformPotentialToCartesianCoordinates(_phi, cphi);
    Vec3 fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    double energy = 0.0;
    for (int i = 0; i < _numParticles; i++) {

        // Compute the torque.

        multipole[0] = particleData[i].coreCharge + particleData[i].valenceCharge;

        multipole[1] =  particleData[i].dipole[0];
        multipole[2] =  particleData[i].dipole[1];
        multipole[3] =  particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];

        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        const double* phi = &cphi[10*i];
        torques[i][0] += _electric*(multipole[3]*phi[2] - multipole[2]*phi[3]
                      + 2.0*(multipole[6]-multipole[5])*phi[9]
                      + multipole[8]*phi[7] + multipole[9]*phi[5]
                      - multipole[7]*phi[8] - multipole[9]*phi[6]);

        torques[i][1] += _electric*(multipole[1]*phi[3] - multipole[3]*phi[1]
                      + 2.0*(multipole[4]-multipole[6])*phi[8]
                      + multipole[7]*phi[9] + multipole[8]*phi[6]
                      - multipole[8]*phi[4] - multipole[9]*phi[7]);

        torques[i][2] += _electric*(multipole[2]*phi[1] - multipole[1]*phi[2]
                      + 2.0*(multipole[5]-multipole[4])*phi[7]
                      + multipole[7]*phi[4] + multipole[9]*phi[8]
                      - multipole[7]*phi[5] - multipole[8]*phi[9]);

        // Compute the force and energy.

        multipole[1] = _transformed[i].dipole[0];
        multipole[2] = _transformed[i].dipole[1];
        multipole[3] = _transformed[i].dipole[2];
        multipole[4] = _transformed[i].quadrupole[QXX];
        multipole[5] = _transformed[i].quadrupole[QYY];
        multipole[6] = _transformed[i].quadrupole[QZZ];
        multipole[7] = _transformed[i].quadrupole[QXY];
        multipole[8] = _transformed[i].quadrupole[QXZ];
        multipole[9] = _transformed[i].quadrupole[QYZ];

        Vec3 f = Vec3(0.0, 0.0, 0.0);
        for (int k = 0; k < 10; k++) {
            energy += multipole[k]*_phi[20*i+k];
            f[0]   -= multipole[k]*_phi[20*i+deriv1[k]];
            f[1]   -= multipole[k]*_phi[20*i+deriv2[k]];
            f[2]   -= multipole[k]*_phi[20*i+deriv3[k]];
        }
        f *= _electric;
        forces[i] += Vec3(f[0]*fracToCart[0][0] + f[1]*fracToCart[0][1] + f[2]*fracToCart[0][2],
                          f[0]*fracToCart[1][0] + f[1]*fracToCart[1][1] + f[2]*fracToCart[1][2],
                          f[0]*fracToCart[2][0] + f[1]*fracToCart[2][1] + f[2]*fracToCart[2][2]);
    }
    return (0.5*_electric*energy);
}

/**
 * Compute the forces due to the reciprocal space PME calculation for induced dipoles.
 */
void AmoebaReferencePmeHippoNonbondedForce::computeReciprocalSpaceInducedDipoleForce(const vector<MultipoleParticleData>& particleData,
                                                                                     vector<Vec3>& forces, vector<Vec3>& torques) const {
    double multipole[10];
    Vec3 inducedDipole;
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    vector<double> cphi(10*_numParticles);
    transformPotentialToCartesianCoordinates(_phidp, cphi);
    Vec3 cartToFrac[3], fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cartToFrac[j][i] = fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    for (int i = 0; i < _numParticles; i++) {
        // Compute the torque.

        multipole[0] = particleData[i].coreCharge + particleData[i].valenceCharge;

        multipole[1] = particleData[i].dipole[0];
        multipole[2] = particleData[i].dipole[1];
        multipole[3] = particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];
        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        const double* phi = &cphi[10*i];
        torques[i][0] += _electric*(multipole[3]*phi[2] - multipole[2]*phi[3]
                      + 2.0*(multipole[6]-multipole[5])*phi[9]
                      + multipole[8]*phi[7] + multipole[9]*phi[5]
                      - multipole[7]*phi[8] - multipole[9]*phi[6]);

        torques[i][1] += _electric*(multipole[1]*phi[3] - multipole[3]*phi[1]
                      + 2.0*(multipole[4]-multipole[6])*phi[8]
                      + multipole[7]*phi[9] + multipole[8]*phi[6]
                      - multipole[8]*phi[4] - multipole[9]*phi[7]);

        torques[i][2] += _electric*(multipole[2]*phi[1] - multipole[1]*phi[2]
                      + 2.0*(multipole[5]-multipole[4])*phi[7]
                      + multipole[7]*phi[4] + multipole[9]*phi[8]
                      - multipole[7]*phi[5] - multipole[8]*phi[9]);

        // Compute the force and energy.

        multipole[1] = _transformed[i].dipole[0];
        multipole[2] = _transformed[i].dipole[1];
        multipole[3] = _transformed[i].dipole[2];
        multipole[4] = _transformed[i].quadrupole[QXX];
        multipole[5] = _transformed[i].quadrupole[QYY];
        multipole[6] = _transformed[i].quadrupole[QZZ];
        multipole[7] = _transformed[i].quadrupole[QXY];
        multipole[8] = _transformed[i].quadrupole[QXZ];
        multipole[9] = _transformed[i].quadrupole[QYZ];

        inducedDipole[0] = _inducedDipole[i][0]*cartToFrac[0][0] + _inducedDipole[i][1]*cartToFrac[0][1] + _inducedDipole[i][2]*cartToFrac[0][2];
        inducedDipole[1] = _inducedDipole[i][0]*cartToFrac[1][0] + _inducedDipole[i][1]*cartToFrac[1][1] + _inducedDipole[i][2]*cartToFrac[1][2];
        inducedDipole[2] = _inducedDipole[i][0]*cartToFrac[2][0] + _inducedDipole[i][1]*cartToFrac[2][1] + _inducedDipole[i][2]*cartToFrac[2][2];

        Vec3 f = Vec3(0.0, 0.0, 0.0);

        for (int k = 0; k < 3; k++) {
            f[0] -= inducedDipole[k]*_phi[20*i+deriv1[k+1]];
            f[1] -= inducedDipole[k]*_phi[20*i+deriv2[k+1]];
            f[2] -= inducedDipole[k]*_phi[20*i+deriv3[k+1]];
        }

        for (int k = 0; k < 10; k++) {
            f[0] -= multipole[k]*_phidp[20*i+deriv1[k]];
            f[1] -= multipole[k]*_phidp[20*i+deriv2[k]];
            f[2] -= multipole[k]*_phidp[20*i+deriv3[k]];
        }

        // Account for dipole response terms in the OPT method

        for (int j = 0; j < _maxPTOrder-1; j++) {
            for (int m = 0; m < _maxPTOrder-1-j; m++) {
                Vec3 optDipole(_ptDipoleD[m][i][0]*cartToFrac[0][0] + _ptDipoleD[m][i][1]*cartToFrac[0][1] + _ptDipoleD[m][i][2]*cartToFrac[0][2],
                               _ptDipoleD[m][i][0]*cartToFrac[1][0] + _ptDipoleD[m][i][1]*cartToFrac[1][1] + _ptDipoleD[m][i][2]*cartToFrac[1][2],
                               _ptDipoleD[m][i][0]*cartToFrac[2][0] + _ptDipoleD[m][i][1]*cartToFrac[2][1] + _ptDipoleD[m][i][2]*cartToFrac[2][2]);
                Vec3 h;
                for (int k = 0; k < 3; k++) {
                    h[0] += optDipole[k]*optPhi[j][10*i+deriv1[k+1]];
                    h[1] += optDipole[k]*optPhi[j][10*i+deriv2[k+1]];
                    h[2] += optDipole[k]*optPhi[j][10*i+deriv3[k+1]];
                }
                f -= _extPartCoefficients[j+m+1]*h;
            }
        }
        forces[i] += _electric*Vec3(f[0]*fracToCart[0][0] + f[1]*fracToCart[0][1] + f[2]*fracToCart[0][2],
                                    f[0]*fracToCart[1][0] + f[1]*fracToCart[1][1] + f[2]*fracToCart[1][2],
                                    f[0]*fracToCart[2][0] + f[1]*fracToCart[2][1] + f[2]*fracToCart[2][2]);
    }
}

void AmoebaReferencePmeHippoNonbondedForce::recordFixedMultipoleField() {
    Vec3 fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    for (int i = 0; i < _numParticles; i++) {
        _fixedMultipoleField[i][0] = -(_phi[20*i+1]*fracToCart[0][0] + _phi[20*i+2]*fracToCart[0][1] + _phi[20*i+3]*fracToCart[0][2]);
        _fixedMultipoleField[i][1] = -(_phi[20*i+1]*fracToCart[1][0] + _phi[20*i+2]*fracToCart[1][1] + _phi[20*i+3]*fracToCart[1][2]);
        _fixedMultipoleField[i][2] = -(_phi[20*i+1]*fracToCart[2][0] + _phi[20*i+2]*fracToCart[2][1] + _phi[20*i+3]*fracToCart[2][2]);
    }
}

void AmoebaReferencePmeHippoNonbondedForce::initializeInducedDipoles() {
    AmoebaReferenceHippoNonbondedForce::initializeInducedDipoles();
    calculateReciprocalSpaceInducedDipoleField();
}

void AmoebaReferencePmeHippoNonbondedForce::recordInducedDipoleField(vector<Vec3>& field) {
    Vec3 fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    for (int i = 0; i < _numParticles; i++) {
        field[i][0] -= _phidp[20*i+1]*fracToCart[0][0] + _phidp[20*i+2]*fracToCart[0][1] + _phidp[20*i+3]*fracToCart[0][2];
        field[i][1] -= _phidp[20*i+1]*fracToCart[1][0] + _phidp[20*i+2]*fracToCart[1][1] + _phidp[20*i+3]*fracToCart[1][2];
        field[i][2] -= _phidp[20*i+1]*fracToCart[2][0] + _phidp[20*i+2]*fracToCart[2][1] + _phidp[20*i+3]*fracToCart[2][2];
    }
}

void AmoebaReferencePmeHippoNonbondedForce::calculateReciprocalSpaceInducedDipoleField() {
    // Perform PME for the induced dipoles.

    initializePmeGrid();
    spreadInducedDipolesOnGrid(_inducedDipole);
    fftpack_exec_3d(_fftplan, FFTPACK_FORWARD, _pmeGrid.data(), _pmeGrid.data());
    performAmoebaReciprocalConvolution();
    fftpack_exec_3d(_fftplan, FFTPACK_BACKWARD, _pmeGrid.data(), _pmeGrid.data());
    computeInducedPotentialFromGrid();
    recordInducedDipoleField(_inducedDipoleField);
}

void AmoebaReferencePmeHippoNonbondedForce::calculateInducedDipoleFields(const vector<MultipoleParticleData>& particleData, int optOrder) {
    // Initialize the fields to zero.

    Vec3 zeroVec(0.0, 0.0, 0.0);
    std::fill(_inducedDipoleField.begin(), _inducedDipoleField.end(), zeroVec);

    // Add fields from direct space interactions.

    for (int i = 0; i < _numParticles; i++)
        for (int j = i+1; j < _numParticles; j++)
            calculateDirectInducedDipolePairIxns(particleData[i], particleData[j]);

    // reciprocal space ixns

    calculateReciprocalSpaceInducedDipoleField();
    for (int i = 0; i < _numParticles; i++)
        for (int j = 0; j < 10; j++)
            optPhi[optOrder][10*i+j] = _phidp[20*i+j];

    // self ixn

    double term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (int j = 0; j < _numParticles; j++)
        _inducedDipoleField[j] += _inducedDipole[j]*term;
}

void AmoebaReferencePmeHippoNonbondedForce::calculateDirectInducedDipolePairIxn(int iIndex, int jIndex,
                                                                                double preFactor1, double preFactor2,
                                                                                const Vec3& delta,
                                                                                const vector<Vec3>& inducedDipole,
                                                                                vector<Vec3>& field) const {
    // field at i due induced dipole at j

    double dur  = inducedDipole[jIndex].dot(delta);
    field[iIndex]  += delta*(dur*preFactor2) + inducedDipole[jIndex]*preFactor1;

    // field at j due induced dipole at i

               dur  = inducedDipole[iIndex].dot(delta);
    field[jIndex]  += delta*(dur*preFactor2) + inducedDipole[iIndex]*preFactor1;
}

void AmoebaReferencePmeHippoNonbondedForce::calculateDirectInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                                                                 const MultipoleParticleData& particleJ) {
    int i = particleI.index;
    int j = particleJ.index;
    if (i == j)
        return;

    Vec3 deltaR = particleJ.position - particleI.position;
    getPeriodicDelta(deltaR);
    double r2 = deltaR.dot(deltaR);
    if (r2 > _cutoffDistanceSquared)
        return;
    double r = sqrt(r2);
    double fdamp3, fdamp5;
    computeMutualFieldDampingFactors(particleI, particleJ, r, fdamp3, fdamp5);
    auto exception = exceptions.find(make_pair(i, j));
    if (exception != exceptions.end()) {
        fdamp3 *= exception->second.dipoleDipoleScale;
        fdamp5 *= exception->second.dipoleDipoleScale;
    }
    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rInv3 = rInv*rInv2;

    // calculate the error function damping terms

    double ralpha = _alphaEwald*r;
    double bn0 = erfc(ralpha)*rInv;
    double alsq2 = 2*_alphaEwald*_alphaEwald;
    double alsq2n = 1/(SQRT_PI*_alphaEwald);
    double exp2a = exp(-(ralpha*ralpha));
    alsq2n *= alsq2;
    double bn1 = (bn0+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn2 = (3*bn1+alsq2n*exp2a)*rInv2;
    double scale3 = -bn1 + (1-fdamp3)*rInv3;
    double scale5 = bn2 - 3*(1-fdamp5)*rInv3*rInv2;
    _inducedDipoleField[i] += _inducedDipole[j]*scale3 + deltaR*scale5*(_inducedDipole[j].dot(deltaR));
    _inducedDipoleField[j] += _inducedDipole[i]*scale3 + deltaR*scale5*(_inducedDipole[i].dot(deltaR));
}

double AmoebaReferencePmeHippoNonbondedForce::calculatePmeSelfEnergy(const vector<MultipoleParticleData>& particleData) const {
    double cii = 0.0;
    double dii = 0.0;
    double qii = 0.0;
    double c6ii = 0.0;
    for (int i = 0; i < _numParticles; i++) {
        const MultipoleParticleData& particleI = particleData[i];
        double charge = particleI.coreCharge + particleI.valenceCharge;
        cii += charge*charge;
        dii += particleI.dipole.dot(particleI.dipole);
        qii += particleI.quadrupole[QXX]*particleI.quadrupole[QXX] +
               particleI.quadrupole[QYY]*particleI.quadrupole[QYY] +
               particleI.quadrupole[QZZ]*particleI.quadrupole[QZZ] +
               2*(particleI.quadrupole[QXY]*particleI.quadrupole[QXY] +
                  particleI.quadrupole[QXZ]*particleI.quadrupole[QXZ] +
                  particleI.quadrupole[QYZ]*particleI.quadrupole[QYZ]);
        c6ii += particleI.c6*particleI.c6;
    }
    double term = 2*_alphaEwald*_alphaEwald;
    double fterm = -_electric*_alphaEwald/SQRT_PI;
    double alpha3 = _dalphaEwald*_dalphaEwald*_dalphaEwald;
    return fterm*(cii + term*(dii/3+2*term*qii/5)) + alpha3*alpha3*c6ii/12;
}

void AmoebaReferencePmeHippoNonbondedForce::calculatePmeSelfTorque(const vector<MultipoleParticleData>& particleData,
                                                                   vector<Vec3>& torques) const {
    double term = (4.0/3.0)*_electric*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (int i = 0; i < _numParticles; i++) {
        const MultipoleParticleData& particleI = particleData[i];
        Vec3 torque = particleI.dipole.cross(_inducedDipole[i])*term;
        torques[i] += torque;
    }
}

double AmoebaReferencePmeHippoNonbondedForce::calculateElectrostaticPairIxn(const MultipoleParticleData& particleI,
                                                                            const MultipoleParticleData& particleK,
                                                                            double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const {
    double dir = particleI.qiDipole[2]*r;
    Vec3 qxI = Vec3(particleI.qiQuadrupole[QXX], particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QXZ]);
    Vec3 qyI = Vec3(particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QYY], particleI.qiQuadrupole[QYZ]);
    Vec3 qzI = Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    Vec3 qi = r*Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    double qir = qi[2]*r;
    double dkr = particleK.qiDipole[2]*r;
    Vec3 qxK = Vec3(particleK.qiQuadrupole[QXX], particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QXZ]);
    Vec3 qyK = Vec3(particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QYY], particleK.qiQuadrupole[QYZ]);
    Vec3 qzK = Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    Vec3 qk = r*Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    double qkr = qk[2]*r;
    double dik = particleI.qiDipole.dot(particleK.qiDipole);
    double qik = qi.dot(qk);
    double diqk = particleI.qiDipole.dot(qk);
    double dkqi = particleK.qiDipole.dot(qi);
    double qiqk = 2*(qxI[1]*qxK[1]+qxI[2]*qxK[2]+qyI[2]*qyK[2]) + qxI[0]*qxK[0] + qyI[1]*qyK[1] + qzI[2]*qzK[2];

    // Additional intermediates involving moments and distance.

    Vec3 dirCross(particleI.qiDipole[1]*r, -particleI.qiDipole[0]*r, 0);
    Vec3 dkrCross(particleK.qiDipole[1]*r, -particleK.qiDipole[0]*r, 0);
    Vec3 dikCross = particleI.qiDipole.cross(particleK.qiDipole);
    Vec3 qirCross(qi[1]*r, -qi[0]*r, 0);
    Vec3 qkrCross(qk[1]*r, -qk[0]*r, 0);
    Vec3 qikCross = qk.cross(qi);
    Vec3 qikTemp(qxI.dot(qk), qyI.dot(qk), qzI.dot(qk));
    Vec3 qkiTemp(qxK.dot(qi), qyK.dot(qi), qzK.dot(qi));
    Vec3 qikrCross(-r*qikTemp[1], r*qikTemp[0], 0);
    Vec3 qkirCross(-r*qkiTemp[1], r*qkiTemp[0], 0);
    Vec3 diqkTemp(particleI.qiDipole.dot(qxK), particleI.qiDipole.dot(qyK), particleI.qiDipole.dot(qzK));
    Vec3 dkqiTemp(particleK.qiDipole.dot(qxI), particleK.qiDipole.dot(qyI), particleK.qiDipole.dot(qzI));
    Vec3 diqkrCross(-r*diqkTemp[1], r*diqkTemp[0], 0);
    Vec3 dkqirCross(-r*dkqiTemp[1], r*dkqiTemp[0], 0);
    Vec3 dqik = particleI.qiDipole.cross(qk) + particleK.qiDipole.cross(qi) - 2*(qxI.cross(qxK) + qyI.cross(qyK) + qzI.cross(qzK));

    // Get reciprocal distance terms for this interaction.

    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rr1 = rInv;
    double rr3 = rr1*rInv2;
    double rr5 = 3*rr3*rInv2;
    double rr7 = 5*rr5*rInv2;
    double rr9 = 7*rr7*rInv2;
    double rr11 = 9*rr9*rInv2;

    // Calculate the error function damping terms.

    double ralpha = _alphaEwald*r;
    double bn0 = erfc(ralpha)*rInv;
    double alsq2 = 2*_alphaEwald*_alphaEwald;
    double alsq2n = 1/(SQRT_PI*_alphaEwald);
    double exp2a = exp(-(ralpha*ralpha));
    alsq2n *= alsq2;
    double bn1 = (bn0+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn2 = (3*bn1+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn3 = (5*bn2+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn4 = (7*bn3+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn5 = (9*bn4+alsq2n*exp2a)*rInv2;

    // Find damped multipole intermediates and energy value.

    double term1 = particleI.coreCharge*particleK.coreCharge;
    double term1i = particleK.coreCharge*particleI.valenceCharge;
    double term2i = particleK.coreCharge*dir;
    double term3i = particleK.coreCharge*qir;
    double term1k = particleI.coreCharge*particleK.valenceCharge;
    double term2k = -particleI.coreCharge*dkr;
    double term3k = particleI.coreCharge*qkr;
    double term1ik = particleI.valenceCharge*particleK.valenceCharge;
    double term2ik = particleK.valenceCharge*dir - particleI.valenceCharge*dkr + dik;
    double term3ik = particleI.valenceCharge*qkr + particleK.valenceCharge*qir - dir*dkr + 2*(dkqi-diqk+qiqk);
    double term4ik = dir*qkr - dkr*qir - 4*qik;
    double term5ik = qir*qkr;
    double fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    double fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    double fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(particleI, particleK, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    double scale = 1;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end())
        scale = exception->second.multipoleMultipoleScale;
    double rr1i = bn0 - (1-scale*fdampI1)*rr1;
    double rr3i = bn1 - (1-scale*fdampI3)*rr3;
    double rr5i = bn2 - (1-scale*fdampI5)*rr5;
    double rr7i = bn3 - (1-scale*fdampI7)*rr7;
    double rr1k = bn0 - (1-scale*fdampK1)*rr1;
    double rr3k = bn1 - (1-scale*fdampK3)*rr3;
    double rr5k = bn2 - (1-scale*fdampK5)*rr5;
    double rr7k = bn3 - (1-scale*fdampK7)*rr7;
    double rr1ik = bn0 - (1-scale*fdampIK1)*rr1;
    double rr3ik = bn1 - (1-scale*fdampIK3)*rr3;
    double rr5ik = bn2 - (1-scale*fdampIK5)*rr5;
    double rr7ik = bn3 - (1-scale*fdampIK7)*rr7;
    double rr9ik = bn4 - (1-scale*fdampIK9)*rr9;
    double rr11ik = bn5 - (1-scale*fdampIK11)*rr11;
    double rr1s = bn0 - (1-scale)*rr1;
    double rr3s = bn1 - (1-scale)*rr3;
    double energy = _electric*(term1*rr1s + term4ik*rr7ik + term5ik*rr9ik +
                               term1i*rr1i + term1k*rr1k + term1ik*rr1ik +
                               term2i*rr3i + term2k*rr3k + term2ik*rr3ik +
                               term3i*rr5i + term3k*rr5k + term3ik*rr5ik);

    // Find damped multipole intermediates for force and torque.

    double de = term1*rr3s + term4ik*rr9ik + term5ik*rr11ik +
                term1i*rr3i + term1k*rr3k + term1ik*rr3ik +
                term2i*rr5i + term2k*rr5k + term2ik*rr5ik +
                term3i*rr7i + term3k*rr7k + term3ik*rr7ik;
    term1 = -particleK.coreCharge*rr3i - particleK.valenceCharge*rr3ik + dkr*rr5ik - qkr*rr7ik;
    double term2 = particleI.coreCharge*rr3k + particleI.valenceCharge*rr3ik + dir*rr5ik + qir*rr7ik;
    double term3 = 2 * rr5ik;
    double term4 = -2 * (particleK.coreCharge*rr5i+particleK.valenceCharge*rr5ik-dkr*rr7ik+qkr*rr9ik);
    double term5 = -2 * (particleI.coreCharge*rr5k+particleI.valenceCharge*rr5ik+dir*rr7ik+qir*rr9ik);
    double term6 = 4 * rr7ik;

    // Compute the force and torque.

    force += _electric*(Vec3(0, 0, de*r) + term1*particleI.qiDipole + term2*particleK.qiDipole +
            term3*(diqkTemp-dkqiTemp) + term4*qi + term5*qk + term6*(qikTemp+qkiTemp));
    torqueI += _electric*(-rr3ik*dikCross + term1*dirCross + term3*(dqik+dkqirCross) + term4*qirCross - term6*(qikrCross+qikCross));
    torqueK += _electric*(rr3ik*dikCross + term2*dkrCross - term3*(dqik+diqkrCross) + term5*qkrCross - term6*(qkirCross-qikCross));
    return energy;
}

void AmoebaReferencePmeHippoNonbondedForce::calculateInducedDipolePairIxn(const MultipoleParticleData& particleI,
                                                                          const MultipoleParticleData& particleK,
                                                                          Vec3 deltaR, double r, Vec3& force, Vec3& torqueI, Vec3& torqueK, Vec3& labForce) const {
    double r2 = r*r;

    // Intermediates involving moments and separation distance

    double dir = particleI.qiDipole[2]*r;
    Vec3 qxI = Vec3(particleI.qiQuadrupole[QXX], particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QXZ]);
    Vec3 qyI = Vec3(particleI.qiQuadrupole[QXY], particleI.qiQuadrupole[QYY], particleI.qiQuadrupole[QYZ]);
    Vec3 qzI = Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    Vec3 qi = r*Vec3(particleI.qiQuadrupole[QXZ], particleI.qiQuadrupole[QYZ], particleI.qiQuadrupole[QZZ]);
    double qir = qi[2]*r;
    double dkr = particleK.qiDipole[2]*r;
    Vec3 qxK = Vec3(particleK.qiQuadrupole[QXX], particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QXZ]);
    Vec3 qyK = Vec3(particleK.qiQuadrupole[QXY], particleK.qiQuadrupole[QYY], particleK.qiQuadrupole[QYZ]);
    Vec3 qzK = Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    Vec3 qk = r*Vec3(particleK.qiQuadrupole[QXZ], particleK.qiQuadrupole[QYZ], particleK.qiQuadrupole[QZZ]);
    double qkr = qk[2]*r;
    const Vec3& ui = particleI.qiInducedDipole;
    const Vec3& uk = particleK.qiInducedDipole;
    double uir = ui[2]*r;
    double ukr = uk[2]*r;

    // Get reciprocal distance terms for this interaction.

    double rInv = 1/r;
    double rInv2 = rInv*rInv;
    double rr1 = rInv;
    double rr3 = rr1*rInv2;
    double rr5 = 3*rr3*rInv2;
    double rr7 = 5*rr5*rInv2;
    double rr9 = 7*rr7*rInv2;

    // calculate the error function damping terms

    double ralpha = _alphaEwald*r;
    double bn0 = erfc(ralpha)*rInv;
    double alsq2 = 2*_alphaEwald*_alphaEwald;
    double alsq2n = 1/(SQRT_PI*_alphaEwald);
    double exp2a = exp(-(ralpha*ralpha));
    alsq2n *= alsq2;
    double bn1 = (bn0+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn2 = (3*bn1+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn3 = (5*bn2+alsq2n*exp2a)*rInv2;
    alsq2n *= alsq2;
    double bn4 = (7*bn3+alsq2n*exp2a)*rInv2;

    // Apply charge penetration damping to scale factors.

    double fdampI1, fdampI3, fdampI5, fdampI7, fdampI9;
    double fdampK1, fdampK3, fdampK5, fdampK7, fdampK9;
    double fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11;
    computeOverlapDampingFactors(particleI, particleK, r, fdampI1, fdampI3, fdampI5, fdampI7, fdampI9, fdampK1, fdampK3, fdampK5, fdampK7, fdampK9,
                                 fdampIK1, fdampIK3, fdampIK5, fdampIK7, fdampIK9, fdampIK11);
    double dipoleMultipoleScale = 1, dipoleDipoleScale = 1;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end()) {
        dipoleMultipoleScale = exception->second.dipoleMultipoleScale;
        dipoleDipoleScale = exception->second.dipoleDipoleScale;
    }
    double rr3core = _electric*(bn1 - (1-dipoleMultipoleScale)*rr3);
    double rr5core = _electric*(bn2 - (1-dipoleMultipoleScale)*rr5);
    double rr3i = _electric*(bn1 - (1-dipoleMultipoleScale*fdampI3)*rr3);
    double rr5i = _electric*(bn2 - (1-dipoleMultipoleScale*fdampI5)*rr5);
    double rr7i = _electric*(bn3 - (1-dipoleMultipoleScale*fdampI7)*rr7);
    double rr9i = _electric*(bn4 - (1-dipoleMultipoleScale*fdampI9)*rr9);
    double rr3k = _electric*(bn1 - (1-dipoleMultipoleScale*fdampK3)*rr3);
    double rr5k = _electric*(bn2 - (1-dipoleMultipoleScale*fdampK5)*rr5);
    double rr7k = _electric*(bn3 - (1-dipoleMultipoleScale*fdampK7)*rr7);
    double rr9k = _electric*(bn4 - (1-dipoleMultipoleScale*fdampK9)*rr9);
    double rr5ik = _electric*(bn2 - (1-dipoleDipoleScale*fdampIK5)*rr5);
    double rr7ik = _electric*(bn3 - (1-dipoleDipoleScale*fdampIK7)*rr7);

    // Get the induced dipole field used for dipole torques.

    Vec3 torqueFieldI = rr3i*uk;
    torqueFieldI[2] -= rr5i*ukr*r;
    Vec3 torqueFieldK = rr3k*ui;
    torqueFieldK[2] -= rr5k*uir*r;

    // Get induced dipole field gradient used for quadrupole torques.

    Vec3 dtorqueFieldI = 2*r*rr5i*uk;
    dtorqueFieldI[2] -= r2*rr7i*ukr;
    Vec3 dtorqueFieldK = -2*r*rr5k*ui;
    dtorqueFieldK[2] += r2*rr7k*uir;

    // Get the field gradient for direct polarization force

    double ti[6], tk[6];
    ti[QXX] = particleI.valenceCharge*rr3i + particleI.coreCharge*rr3core + dir*rr5i - qxI[0]*2*rr5i + qi[2]*r*rr7i;
    tk[QXX] = particleK.valenceCharge*rr3k + particleK.coreCharge*rr3core - dkr*rr5k - qxK[0]*2*rr5k + qk[2]*r*rr7k;
    ti[QYY] = particleI.valenceCharge*rr3i + particleI.coreCharge*rr3core + dir*rr5i - qyI[1]*2*rr5i + qi[2]*r*rr7i;
    tk[QYY] = particleK.valenceCharge*rr3k + particleK.coreCharge*rr3core - dkr*rr5k - qyK[1]*2*rr5k + qk[2]*r*rr7k;
    ti[QZZ] = particleI.valenceCharge*(rr3i-rr5i*r2) + particleI.coreCharge*(rr3core-rr5core*r2) + particleI.qiDipole[2]*2*rr5i*r -
              dir*(rr7i*r2-rr5i) - qzI[2]*2*rr5i + qi[2]*5*rr7i*r - qir*rr9i*r2;
    tk[QZZ] = particleK.valenceCharge*(rr3k-rr5k*r2) + particleK.coreCharge*(rr3core-rr5core*r2) - particleK.qiDipole[2]*2*rr5k*r +
              dkr*(rr7k*r2-rr5k) - qzK[2]*2*rr5k + qk[2]*5*rr7k*r - qkr*rr9k*r2;
    ti[QXY] = -qxI[1]*2*rr5i;
    tk[QXY] = -qxK[1]*2*rr5k;
    ti[QXZ] = particleI.qiDipole[0]*rr5i*r - qxI[2]*2*rr5i + qi[0]*2*rr7i*r;
    tk[QXZ] = -particleK.qiDipole[0]*rr5k*r - qxK[2]*2*rr5k + qk[0]*2*rr7k*r;
    ti[QYZ] = particleI.qiDipole[1]*rr5i*r - qyI[2]*2*rr5i + qi[1]*2*rr7i*r;
    tk[QYZ] = -particleK.qiDipole[1]*rr5k*r - qyK[2]*2*rr5k + qk[1]*2*rr7k*r;

    // Get the dEp/dR terms for chgpen direct polarization force.

    double depx = ti[QXX]*uk[0] + ti[QXY]*uk[1] + ti[QXZ]*uk[2] - tk[QXX]*ui[0] - tk[QXY]*ui[1] - tk[QXZ]*ui[2];
    double depy = ti[QXY]*uk[0] + ti[QYY]*uk[1] + ti[QYZ]*uk[2] - tk[QXY]*ui[0] - tk[QYY]*ui[1] - tk[QYZ]*ui[2];
    double depz = ti[QXZ]*uk[0] + ti[QYZ]*uk[1] + ti[QZZ]*uk[2] - tk[QXZ]*ui[0] - tk[QYZ]*ui[1] - tk[QZZ]*ui[2];
    force += Vec3(depx, depy, depz);

    // Get the dtau/dr terms used for OPT polarization force.

    int i = particleI.index;
    int k = particleK.index;
    for (int j = 0; j < _maxPTOrder-1; j++) {
        Vec3 extDipole1 = _ptDipoleD[j][i];
        double uirm = extDipole1.dot(deltaR);
        for (int m = 0; m < _maxPTOrder-1-j; m++) {
            Vec3 extDipole2 = _ptDipoleD[m][k];
            double ukrm = extDipole2.dot(deltaR);
            double term1 = 2*rr5ik;
            double term2 = term1*deltaR[0];
            double term3 = rr5ik - rr7ik*deltaR[0]*deltaR[0];
            double tixx = extDipole1[0]*term2 + uirm*term3;
            double tkxx = extDipole2[0]*term2 + ukrm*term3;
            term2 = term1*deltaR[1];
            term3 = rr5ik - rr7ik*deltaR[1]*deltaR[1];
            double tiyy = extDipole1[1]*term2 + uirm*term3;
            double tkyy = extDipole2[1]*term2 + ukrm*term3;
            term2 = term1*deltaR[2];
            term3 = rr5ik - rr7ik*deltaR[2]*deltaR[2];
            double tizz = extDipole1[2]*term2 + uirm*term3;
            double tkzz = extDipole2[2]*term2 + ukrm*term3;
            term1 = rr5ik*deltaR[1];
            term2 = rr5ik*deltaR[0];
            term3 = deltaR[1] * (rr7ik*deltaR[0]);
            double tixy = extDipole1[0]*term1 + extDipole1[1]*term2 - uirm*term3;
            double tkxy = extDipole2[0]*term1 + extDipole2[1]*term2 - ukrm*term3;
            term1 = rr5ik * deltaR[2];
            term3 = deltaR[2] * (rr7ik*deltaR[0]);
            double tixz = extDipole1[0]*term1 + extDipole1[2]*term2 - uirm*term3;
            double tkxz = extDipole2[0]*term1 + extDipole2[2]*term2 - ukrm*term3;
            term2 = rr5ik*deltaR[1];
            term3 = deltaR[2] * (rr7ik*deltaR[1]);
            double tiyz = extDipole1[1]*term1 + extDipole1[2]*term2 - uirm*term3;
            double tkyz = extDipole2[1]*term1 + extDipole2[2]*term2 - ukrm*term3;
            double depx = tixx*extDipole2[0] + tkxx*extDipole1[0]
                 + tixy*extDipole2[1] + tkxy*extDipole1[1]
                 + tixz*extDipole2[2] + tkxz*extDipole1[2];
            double depy = tixy*extDipole2[0] + tkxy*extDipole1[0]
                 + tiyy*extDipole2[1] + tkyy*extDipole1[1]
                 + tiyz*extDipole2[2] + tkyz*extDipole1[2];
            double depz = tixz*extDipole2[0] + tkxz*extDipole1[0]
                 + tiyz*extDipole2[1] + tkyz*extDipole1[1]
                 + tizz*extDipole2[2] + tkzz*extDipole1[2];
            labForce += 0.5*_extPartCoefficients[j+m+1]*Vec3(depx, depy, depz);
        }
    }

    // Torque is induced field and gradient cross permanent moments.

    torqueI += torqueFieldI.cross(particleI.qiDipole);
    torqueI[0] += -qxI[1]*dtorqueFieldI[0] - 2*qyI[2]*dtorqueFieldI[2] + (qzI[2]-qyI[1])*dtorqueFieldI[1];
    torqueI[1] += qxI[1]*dtorqueFieldI[1] + 2*qxI[2]*dtorqueFieldI[2] + (qxI[0]-qzI[2])*dtorqueFieldI[0];
    torqueI[2] += qyI[2]*dtorqueFieldI[0] - qxI[2]*dtorqueFieldI[1];
    torqueK += torqueFieldK.cross(particleK.qiDipole);
    torqueK[0] += -qxK[1]*dtorqueFieldK[0] - 2*qyK[2]*dtorqueFieldK[2] + (qzK[2]-qyK[1])*dtorqueFieldK[1];
    torqueK[1] += qxK[1]*dtorqueFieldK[1] + 2*qxK[2]*dtorqueFieldK[2] + (qxK[0]-qzK[2])*dtorqueFieldK[0];
    torqueK[2] += qyK[2]*dtorqueFieldK[0] - qxK[2]*dtorqueFieldK[1];
}

double AmoebaReferencePmeHippoNonbondedForce::calculateDispersionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                                         double r, Vec3& force) const {
    double r2 = r*r;
    double rInv6 = 1/(r2*r2*r2);
    double ralpha2 = r2*_dalphaEwald*_dalphaEwald;
    double expterm = exp(-ralpha2);
    double expa = expterm * (1 + ralpha2 + 0.5*ralpha2*ralpha2);
    
    // Compute the damping function.
    
    double fdamp, ddamp;
    computeDispersionDampingFactors(particleI, particleK, r, fdamp, ddamp);
    
    // Compute the force and energy.

    double dispersionScale = 1;
    auto exception = exceptions.find(make_pair(particleI.index, particleK.index));
    if (exception != exceptions.end())
        dispersionScale = exception->second.dispersionScale;
    double scale = dispersionScale*fdamp*fdamp - 1;
    double cick = particleI.c6*particleK.c6;
    double energy = -cick*(expa+scale)*rInv6;
    double rterm = -ralpha2*ralpha2*ralpha2*expterm/r;
    double dEnergydR = -6*energy/r - cick*rInv6*(rterm + 2*dispersionScale*fdamp*ddamp);
    force[2] -= dEnergydR;
    return energy;
}

double AmoebaReferencePmeHippoNonbondedForce::computeReciprocalSpaceDispersionForceAndEnergy(const vector<MultipoleParticleData>& particleData, vector<Vec3>& forces) const {
    pme_t pmedata;
    pme_init(&pmedata, _dalphaEwald, _numParticles, _dpmeGridDimensions, 5, 1);
    vector<double> charges(_numParticles);
    vector<Vec3> dpmeforces(_numParticles, Vec3()), coords;
    for (int i = 0; i < _numParticles; i++) {
        charges[i] = particleData[i].c6;
        coords.push_back(particleData[i].position);
    }
    double recipDispersionEnergy;
    pme_exec_dpme(pmedata, coords, dpmeforces, charges, _periodicBoxVectors, &recipDispersionEnergy);
    pme_destroy(pmedata);
    for (int i = 0; i < _numParticles; i++)
        forces[i] += dpmeforces[i];
    return recipDispersionEnergy;
}

double AmoebaReferencePmeHippoNonbondedForce::calculateInteractions(vector<Vec3>& torques, vector<Vec3>& forces) {
    double energy = AmoebaReferenceHippoNonbondedForce::calculateInteractions(torques, forces);
    calculatePmeSelfTorque(particleData, torques);
    computeReciprocalSpaceInducedDipoleForce(particleData, forces, torques);
    energy += computeReciprocalSpaceFixedMultipoleForceAndEnergy(particleData, forces, torques);
    energy += computeReciprocalSpaceDispersionForceAndEnergy(particleData, forces);
    energy += calculatePmeSelfEnergy(particleData);
    return energy;
}
