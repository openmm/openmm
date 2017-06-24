
/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
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

#include "AmoebaReferenceMultipoleForce.h"
#include "jama_svd.h"
#include <algorithm>

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

using std::vector;
using namespace OpenMM;

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce() :
                                                   _nonbondedMethod(NoCutoff),
                                                   _numParticles(0),
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324)
{
    initialize();
}

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce(NonbondedMethod nonbondedMethod) :
                                                   _nonbondedMethod(nonbondedMethod),
                                                   _numParticles(0),
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324)
{
    initialize();
}

void AmoebaReferenceMultipoleForce::initialize()
{

    unsigned int index    = 0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.4;
    _mScale[index++]      = 0.8;

    index                 = 0;
    _dScale[index++]      = 0.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;

    index                 = 0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 1.0;
    _pScale[index++]      = 1.0;

    index                 = 0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
}

AmoebaReferenceMultipoleForce::NonbondedMethod AmoebaReferenceMultipoleForce::getNonbondedMethod() const
{
    return _nonbondedMethod;
}

void AmoebaReferenceMultipoleForce::setNonbondedMethod(AmoebaReferenceMultipoleForce::NonbondedMethod nonbondedMethod)
{
    _nonbondedMethod = nonbondedMethod;
}

AmoebaReferenceMultipoleForce::PolarizationType AmoebaReferenceMultipoleForce::getPolarizationType() const
{
    return _polarizationType;
}

void AmoebaReferenceMultipoleForce::setPolarizationType(AmoebaReferenceMultipoleForce::PolarizationType polarizationType)
{
    _polarizationType = polarizationType;
}

int AmoebaReferenceMultipoleForce::getMutualInducedDipoleConverged() const
{
    return _mutualInducedDipoleConverged;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleConverged(int mutualInducedDipoleConverged)
{
    _mutualInducedDipoleConverged = mutualInducedDipoleConverged;
}

int AmoebaReferenceMultipoleForce::getMutualInducedDipoleIterations() const
{
    return _mutualInducedDipoleIterations;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleIterations(int mutualInducedDipoleIterations)
{
    _mutualInducedDipoleIterations = mutualInducedDipoleIterations;
}

RealOpenMM AmoebaReferenceMultipoleForce::getMutualInducedDipoleEpsilon() const
{
    return _mutualInducedDipoleEpsilon;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleEpsilon(RealOpenMM mutualInducedDipoleEpsilon)
{
    _mutualInducedDipoleEpsilon = mutualInducedDipoleEpsilon;
}

int AmoebaReferenceMultipoleForce::getMaximumMutualInducedDipoleIterations() const
{
    return _maximumMutualInducedDipoleIterations;
}

void AmoebaReferenceMultipoleForce::setMaximumMutualInducedDipoleIterations(int maximumMutualInducedDipoleIterations)
{
    _maximumMutualInducedDipoleIterations = maximumMutualInducedDipoleIterations;
}

RealOpenMM AmoebaReferenceMultipoleForce::getMutualInducedDipoleTargetEpsilon() const
{
    return _mutualInducedDipoleTargetEpsilon;
}

void AmoebaReferenceMultipoleForce::setExtrapolationCoefficients(const std::vector<RealOpenMM> &coefficients)
{
    _maxPTOrder = coefficients.size(); // This accounts for the zero-based counting; actual highest order is 1 less
    _extrapolationCoefficients = coefficients;
    _extPartCoefficients.resize(_maxPTOrder);
    for (int i = 0; i < _maxPTOrder; ++i) {
        _extPartCoefficients[i] = 0.0;
        for (int j = i; j < _maxPTOrder; ++j)
            _extPartCoefficients[i] += _extrapolationCoefficients[j];
    }
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleTargetEpsilon(RealOpenMM mutualInducedDipoleTargetEpsilon)
{
    _mutualInducedDipoleTargetEpsilon = mutualInducedDipoleTargetEpsilon;
}

void AmoebaReferenceMultipoleForce::setupScaleMaps(const vector< vector< vector<int> > >& multipoleParticleCovalentInfo)
{

    /* Setup for scaling maps:
     *
     *     _scaleMaps[particleIndex][ScaleType] = map, where map[covalentIndex] = scaleFactor
     *     _maxScaleIndex[particleIndex]        = max covalent index for particleIndex
     *
     *     multipoleParticleCovalentInfo[ii][jj], jj =0,1,2,3 contains covalent indices (c12, c13, c14, c15)
     *     multipoleParticleCovalentInfo[ii][jj], jj =4,5,6,7 contains covalent indices (p11, p12, p13, p14)
     *
     *     only including covalent particles w/ index >= ii
     */

    _scaleMaps.resize(multipoleParticleCovalentInfo.size());
    _maxScaleIndex.resize(multipoleParticleCovalentInfo.size());

    for (unsigned int ii = 0; ii < multipoleParticleCovalentInfo.size(); ii++) {

        _scaleMaps[ii].resize(LAST_SCALE_TYPE_INDEX);
        _maxScaleIndex[ii] = 0;
        const vector< vector<int> >& covalentInfo = multipoleParticleCovalentInfo[ii];
        const vector<int> covalentListP11              = covalentInfo[AmoebaMultipoleForce::PolarizationCovalent11];

        // pScale & mScale

        for (unsigned jj = 0; jj < AmoebaMultipoleForce::PolarizationCovalent11; jj++) {
            const vector<int> covalentList    = covalentInfo[jj];
            for (unsigned int kk = 0; kk < covalentList.size(); kk++) {
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if (covalentIndex < ii)continue;

                // handle 0.5 factor for p14

                int hit = 0;
                if (jj == AmoebaMultipoleForce::Covalent14) {
                    for (unsigned int mm = 0; mm < covalentListP11.size() && hit == 0; mm++) {
                        if (covalentListP11[mm]  == covalentIndex) {
                            hit = 1;
                        }
                    }
                }

                _scaleMaps[ii][P_SCALE][covalentIndex] = hit ? 0.5*_pScale[jj+1] : _pScale[jj+1];
                _scaleMaps[ii][M_SCALE][covalentIndex] = _mScale[jj+1];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }

        // dScale & uScale

        for (unsigned jj = AmoebaMultipoleForce::PolarizationCovalent11; jj < covalentInfo.size(); jj++) {
            const vector<int> covalentList = covalentInfo[jj];
            for (unsigned int kk = 0; kk < covalentList.size(); kk++) {
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if (covalentIndex < ii)continue;
                _scaleMaps[ii][D_SCALE][covalentIndex] = _dScale[jj-4];
                _scaleMaps[ii][U_SCALE][covalentIndex] = _uScale[jj-4];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }
    }
}

RealOpenMM AmoebaReferenceMultipoleForce::getMultipoleScaleFactor(unsigned int particleI, unsigned int particleJ, ScaleType scaleType) const
{

    MapIntRealOpenMM  scaleMap   = _scaleMaps[particleI][scaleType];
    MapIntRealOpenMMCI isPresent = scaleMap.find(particleJ);
    if (isPresent != scaleMap.end()) {
        return isPresent->second;
    } else {
        return 1.0;
    }
}

void AmoebaReferenceMultipoleForce::getDScaleAndPScale(unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale) const
{
    dScale = getMultipoleScaleFactor(particleI, particleJ, D_SCALE);
    pScale = getMultipoleScaleFactor(particleI, particleJ, P_SCALE);
}

void AmoebaReferenceMultipoleForce::getMultipoleScaleFactors(unsigned int particleI, unsigned int particleJ, vector<RealOpenMM>& scaleFactors) const
{
    scaleFactors[D_SCALE] = getMultipoleScaleFactor(particleI, particleJ, D_SCALE);
    scaleFactors[P_SCALE] = getMultipoleScaleFactor(particleI, particleJ, P_SCALE);
    scaleFactors[M_SCALE] = getMultipoleScaleFactor(particleI, particleJ, M_SCALE);
    scaleFactors[U_SCALE] = getMultipoleScaleFactor(particleI, particleJ, U_SCALE);
}

RealOpenMM AmoebaReferenceMultipoleForce::normalizeRealVec(RealVec& vectorToNormalize) const
{
    RealOpenMM norm = SQRT(vectorToNormalize.dot(vectorToNormalize));
    if (norm > 0.0) {
        vectorToNormalize *= (1.0/norm);
    }
    return norm;
}

void AmoebaReferenceMultipoleForce::initializeRealOpenMMVector(vector<RealOpenMM>& vectorToInitialize) const
{
    RealOpenMM zero = 0.0;
    vectorToInitialize.resize(_numParticles);
    std::fill(vectorToInitialize.begin(), vectorToInitialize.end(), zero);
}

void AmoebaReferenceMultipoleForce::initializeRealVecVector(vector<RealVec>& vectorToInitialize) const
{
    vectorToInitialize.resize(_numParticles);
    RealVec zeroVec(0.0, 0.0, 0.0);
    std::fill(vectorToInitialize.begin(), vectorToInitialize.end(), zeroVec);
}

void AmoebaReferenceMultipoleForce::copyRealVecVector(const vector<OpenMM::RealVec>& inputVector, vector<OpenMM::RealVec>& outputVector) const
{
    outputVector.resize(inputVector.size());
    for (unsigned int ii = 0; ii < inputVector.size(); ii++) {
        outputVector[ii] = inputVector[ii];
    }
}

void AmoebaReferenceMultipoleForce::loadParticleData(const vector<RealVec>& particlePositions,
                                                     const vector<RealOpenMM>& charges,
                                                     const vector<RealOpenMM>& dipoles,
                                                     const vector<RealOpenMM>& quadrupoles,
                                                     const vector<RealOpenMM>& tholes,
                                                     const vector<RealOpenMM>& dampingFactors,
                                                     const vector<RealOpenMM>& polarity,
                                                     vector<MultipoleParticleData>& particleData) const
{

    particleData.resize(_numParticles);
    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        particleData[ii].particleIndex        = ii;

        particleData[ii].position             = particlePositions[ii];
        particleData[ii].charge               = charges[ii];

        particleData[ii].dipole[0]            = dipoles[3*ii+0];
        particleData[ii].dipole[1]            = dipoles[3*ii+1];
        particleData[ii].dipole[2]            = dipoles[3*ii+2];

        particleData[ii].quadrupole[QXX]      = quadrupoles[9*ii+0];
        particleData[ii].quadrupole[QXY]      = quadrupoles[9*ii+1];
        particleData[ii].quadrupole[QXZ]      = quadrupoles[9*ii+2];
        particleData[ii].quadrupole[QYY]      = quadrupoles[9*ii+4];
        particleData[ii].quadrupole[QYZ]      = quadrupoles[9*ii+5];
        particleData[ii].quadrupole[QZZ]      = quadrupoles[9*ii+8];

        // Form spherical harmonic dipoles from Cartesian moments.
        particleData[ii].sphericalDipole[0]  = dipoles[3*ii+2]; // z -> Q_10
        particleData[ii].sphericalDipole[1]  = dipoles[3*ii+0]; // x -> Q_11c
        particleData[ii].sphericalDipole[2]  = dipoles[3*ii+1]; // y -> Q_11s

        // Form spherical harmonic quadrupoles from Cartesian moments.
        particleData[ii].sphericalQuadrupole[0] = quadrupoles[9*ii+8]*3.0; // zz -> Q_20
        particleData[ii].sphericalQuadrupole[1] = sqrtFourThirds * quadrupoles[9*ii+2]*3.0; // xz -> Q_21c
        particleData[ii].sphericalQuadrupole[2] = sqrtFourThirds * quadrupoles[9*ii+5]*3.0; // yz -> Q_21s
        particleData[ii].sphericalQuadrupole[3] = sqrtOneThird * (quadrupoles[9*ii+0] - quadrupoles[9*ii+4])*3.0; // xx-yy -> Q_22c
        particleData[ii].sphericalQuadrupole[4] = sqrtFourThirds * quadrupoles[9*ii+1]*3.0; // xy -> Q_22s

        particleData[ii].thole                = tholes[ii];
        particleData[ii].dampingFactor        = dampingFactors[ii];
        particleData[ii].polarity             = polarity[ii];

    }
}

void AmoebaReferenceMultipoleForce::zeroFixedMultipoleFields()
{
    initializeRealVecVector(_fixedMultipoleField);
    initializeRealVecVector(_fixedMultipoleFieldPolar);
}

void AmoebaReferenceMultipoleForce::checkChiralCenterAtParticle(MultipoleParticleData& particleI, int axisType,
                                                                MultipoleParticleData& particleZ, MultipoleParticleData& particleX,
                                                                MultipoleParticleData& particleY) const
{

    if (axisType == AmoebaMultipoleForce::ZThenX || particleY.particleIndex == -1) {
        return;
    }

    RealVec deltaAD   = particleI.position - particleY.position;
    RealVec deltaBD   = particleZ.position - particleY.position;
    RealVec deltaCD   = particleX.position - particleY.position;

    RealVec deltaC    = deltaBD.cross(deltaCD);
    RealOpenMM volume = deltaC.dot(deltaAD);

    if (volume < 0.0) {
        particleI.dipole[1] *= -1.0; // pole(3,i)
        particleI.quadrupole[QXY] *= -1.0; // pole(6,i) && pole(8,i)
        particleI.quadrupole[QYZ] *= -1.0; // pole(10,i) && pole(12,i)
        particleI.sphericalDipole[2]     *= -1.0;   // y
        particleI.sphericalQuadrupole[2] *= -1.0;   // yz
        particleI.sphericalQuadrupole[4] *= -1.0;   // xy
    }
}

void AmoebaReferenceMultipoleForce::checkChiral(vector<MultipoleParticleData>& particleData,
                                                const vector<int>& multipoleAtomXs,
                                                const vector<int>& multipoleAtomYs,
                                                const vector<int>& multipoleAtomZs,
                                                const vector<int>& axisTypes) const
{
    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        if (multipoleAtomYs[ii] > -1) {
            checkChiralCenterAtParticle(particleData[ii], axisTypes[ii],
                                        particleData[multipoleAtomZs[ii]],
                                        particleData[multipoleAtomXs[ii]],
                                        particleData[multipoleAtomYs[ii]]);
        }
    }
}

void AmoebaReferenceMultipoleForce::applyRotationMatrixToParticle(      MultipoleParticleData& particleI,
                                                                  const MultipoleParticleData& particleZ,
                                                                  const MultipoleParticleData& particleX,
                                                                        MultipoleParticleData* particleY,
                                                                        int axisType) const
{

    // handle case where rotation matrix is identity (e.g. single ion)

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom


    RealVec vectorY;
    RealVec vectorZ = particleZ.position - particleI.position;
    RealVec vectorX = particleX.position - particleI.position;

    normalizeRealVec(vectorZ);

    // branch based on axis type

    if (axisType == AmoebaMultipoleForce::Bisector) {

        // bisector

        // dx = dx1 + dx2 (in TINKER code)

        normalizeRealVec(vectorX);
        vectorZ      += vectorX;
        normalizeRealVec(vectorZ);

    } else if (axisType == AmoebaMultipoleForce::ZBisect) {

        // z-bisect

        // dx = dx1 + dx2 (in TINKER code)

        normalizeRealVec(vectorX);

        vectorY  = particleY->position - particleI.position;
        normalizeRealVec(vectorY);

        vectorX += vectorY;
        normalizeRealVec(vectorX);

    } else if (axisType == AmoebaMultipoleForce::ThreeFold) {

        // 3-fold

        // dx = dx1 + dx2 + dx3 (in TINKER code)

        normalizeRealVec(vectorX);

        vectorY   = particleY->position - particleI.position;
        normalizeRealVec(vectorY);

        vectorZ  += vectorX +  vectorY;
        normalizeRealVec(vectorZ);

    } else if (axisType == AmoebaMultipoleForce::ZOnly) {

        // z-only

        vectorX = RealVec(0.1, 0.1, 0.1);

    }

    RealOpenMM dot      = vectorZ.dot(vectorX);
    vectorX            -= vectorZ*dot;

    normalizeRealVec(vectorX);
    vectorY = vectorZ.cross(vectorX);

    RealVec rotationMatrix[3];
    rotationMatrix[0] = vectorX;
    rotationMatrix[1] = vectorY;
    rotationMatrix[2] = vectorZ;

    RealVec labDipole;
    for (int ii = 0; ii < 3; ii++) {
        labDipole[ii] = particleI.dipole[0]*rotationMatrix[0][ii];
        for (int jj = 1; jj < 3; jj++) {
            labDipole[ii] += particleI.dipole[jj]*rotationMatrix[jj][ii];
        }
    }
    particleI.dipole = labDipole;

    RealOpenMM mPole[3][3];
    RealOpenMM rPole[3][3] = { { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 } };

    mPole[0][0] = particleI.quadrupole[QXX];
    mPole[0][1] = particleI.quadrupole[QXY];
    mPole[0][2] = particleI.quadrupole[QXZ];

    mPole[1][0] = particleI.quadrupole[QXY];
    mPole[1][1] = particleI.quadrupole[QYY];
    mPole[1][2] = particleI.quadrupole[QYZ];

    mPole[2][0] = particleI.quadrupole[QXZ];
    mPole[2][1] = particleI.quadrupole[QYZ];
    mPole[2][2] = particleI.quadrupole[QZZ];

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

    RealOpenMM dipoleRotationMatrix[3][3];

    // Reorder the Cartesian {x,y,z} dipole rotation matrix, to account
    // for spherical harmonic ordering {z,x,y}.
    dipoleRotationMatrix[0][0] = vectorZ[2];
    dipoleRotationMatrix[0][1] = vectorX[2];
    dipoleRotationMatrix[0][2] = vectorY[2];
    dipoleRotationMatrix[1][0] = vectorZ[0];
    dipoleRotationMatrix[1][1] = vectorX[0];
    dipoleRotationMatrix[1][2] = vectorY[0];
    dipoleRotationMatrix[2][0] = vectorZ[1];
    dipoleRotationMatrix[2][1] = vectorX[1];
    dipoleRotationMatrix[2][2] = vectorY[1];

    RealOpenMM quadrupoleRotationMatrix[5][5];
    buildSphericalQuadrupoleRotationMatrix(dipoleRotationMatrix, quadrupoleRotationMatrix);

    // Rotate the dipoles
    RealOpenMM rotatedDipole[3];
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM val = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            val += dipoleRotationMatrix[ii][jj] * particleI.sphericalDipole[jj];
        }
        rotatedDipole[ii] = val;
    }
    for (int ii = 0; ii < 3; ii++)
        particleI.sphericalDipole[ii] = rotatedDipole[ii];
    // Rotate the quadrupoles
    RealOpenMM rotatedQuadrupole[5];
    for (int ii = 0; ii < 5; ii++) {
        RealOpenMM val = 0.0;
        for (int jj = 0; jj < 5; jj++) {
            val += quadrupoleRotationMatrix[ii][jj] * particleI.sphericalQuadrupole[jj];
        }
        rotatedQuadrupole[ii] = val;
    }
    for (int ii = 0; ii < 5; ii++)
        particleI.sphericalQuadrupole[ii] = rotatedQuadrupole[ii];
}

void AmoebaReferenceMultipoleForce::formQIRotationMatrix(const RealVec& iPosition,
                                                         const RealVec& jPosition,
                                                         const RealVec &deltaR,
                                                         RealOpenMM r,
                                                         RealOpenMM (&rotationMatrix)[3][3]) const
{
    RealVec vectorZ = (deltaR)/r;
    RealVec vectorX(vectorZ);
    if ((iPosition[1] != jPosition[1]) || (iPosition[2] != jPosition[2])) {
        vectorX[0] += 1.0;
    }else{
        vectorX[1] += 1.0;
    }
    RealVec vectorY;

    RealOpenMM dot      = vectorZ.dot(vectorX);
    vectorX            -= vectorZ*dot;
    normalizeRealVec(vectorX);
    vectorY = vectorZ.cross(vectorX);

    // Reorder the Cartesian {x,y,z} dipole rotation matrix, to account
    // for spherical harmonic ordering {z,x,y}.
    rotationMatrix[0][0] = vectorZ[2];
    rotationMatrix[0][1] = vectorZ[0];
    rotationMatrix[0][2] = vectorZ[1];
    rotationMatrix[1][0] = vectorX[2];
    rotationMatrix[1][1] = vectorX[0];
    rotationMatrix[1][2] = vectorX[1];
    rotationMatrix[2][0] = vectorY[2];
    rotationMatrix[2][1] = vectorY[0];
    rotationMatrix[2][2] = vectorY[1];
}




void AmoebaReferenceMultipoleForce::buildSphericalQuadrupoleRotationMatrix(const RealOpenMM (&D1)[3][3], RealOpenMM (&D2)[5][5]) const
{
    D2[0][0] = 0.5*(3.0*D1[0][0]*D1[0][0] - 1.0);
    D2[1][0] = sqrtThree*D1[0][0]*D1[1][0];
    D2[2][0] = sqrtThree*D1[0][0]*D1[2][0];
    D2[3][0] = 0.5*sqrtThree*(D1[1][0]*D1[1][0] - D1[2][0]*D1[2][0]);
    D2[4][0] = sqrtThree*D1[1][0]*D1[2][0];
    D2[0][1] = sqrtThree*D1[0][0]*D1[0][1];
    D2[1][1] = D1[1][0]*D1[0][1] + D1[0][0]*D1[1][1];
    D2[2][1] = D1[2][0]*D1[0][1] + D1[0][0]*D1[2][1];
    D2[3][1] = D1[1][0]*D1[1][1] - D1[2][0]*D1[2][1];
    D2[4][1] = D1[2][0]*D1[1][1] + D1[1][0]*D1[2][1];
    D2[0][2] = sqrtThree*D1[0][0]*D1[0][2];
    D2[1][2] = D1[1][0]*D1[0][2] + D1[0][0]*D1[1][2];
    D2[2][2] = D1[2][0]*D1[0][2] + D1[0][0]*D1[2][2];
    D2[3][2] = D1[1][0]*D1[1][2] - D1[2][0]*D1[2][2];
    D2[4][2] = D1[2][0]*D1[1][2] + D1[1][0]*D1[2][2];
    D2[0][3] = 0.5*sqrtThree*(D1[0][1]*D1[0][1] - D1[0][2]*D1[0][2]);
    D2[1][3] = D1[0][1]*D1[1][1] - D1[0][2]*D1[1][2];
    D2[2][3] = D1[0][1]*D1[2][1] - D1[0][2]*D1[2][2];
    D2[3][3] = 0.5*(D1[1][1]*D1[1][1] - D1[2][1]*D1[2][1] - D1[1][2]*D1[1][2] + D1[2][2]*D1[2][2]);
    D2[4][3] = D1[1][1]*D1[2][1] - D1[1][2]*D1[2][2];
    D2[0][4] = sqrtThree*D1[0][1]*D1[0][2];
    D2[1][4] = D1[1][1]*D1[0][2] + D1[0][1]*D1[1][2];
    D2[2][4] = D1[2][1]*D1[0][2] + D1[0][1]*D1[2][2];
    D2[3][4] = D1[1][1]*D1[1][2] - D1[2][1]*D1[2][2];
    D2[4][4] = D1[2][1]*D1[1][2] + D1[1][1]*D1[2][2];
}

void AmoebaReferenceMultipoleForce::buildPartialSphericalQuadrupoleRotationMatrix(const RealOpenMM (&D1)[3][3], RealOpenMM (&D2)[3][5]) const
{
    D2[0][0] = 0.5*(3.0*D1[0][0]*D1[0][0] - 1.0);
    D2[0][1] = sqrtThree*D1[0][0]*D1[0][1];
    D2[0][2] = sqrtThree*D1[0][0]*D1[0][2];
    D2[0][3] = 0.5*sqrtThree*(D1[0][1]*D1[0][1] - D1[0][2]*D1[0][2]);
    D2[0][4] = sqrtThree*D1[0][1]*D1[0][2];
    D2[1][0] = sqrtThree*D1[0][0]*D1[1][0];
    D2[1][1] = D1[1][0]*D1[0][1] + D1[0][0]*D1[1][1];
    D2[1][2] = D1[1][0]*D1[0][2] + D1[0][0]*D1[1][2];
    D2[1][3] = D1[0][1]*D1[1][1] - D1[0][2]*D1[1][2];
    D2[1][4] = D1[1][1]*D1[0][2] + D1[0][1]*D1[1][2];
    D2[2][0] = sqrtThree*D1[0][0]*D1[2][0];
    D2[2][1] = D1[2][0]*D1[0][1] + D1[0][0]*D1[2][1];
    D2[2][2] = D1[2][0]*D1[0][2] + D1[0][0]*D1[2][2];
    D2[2][3] = D1[0][1]*D1[2][1] - D1[0][2]*D1[2][2];
    D2[2][4] = D1[2][1]*D1[0][2] + D1[0][1]*D1[2][2];
}

void AmoebaReferenceMultipoleForce::applyRotationMatrix(vector<MultipoleParticleData>& particleData,
                                                        const vector<int>& multipoleAtomXs,
                                                        const vector<int>& multipoleAtomYs,
                                                        const vector<int>& multipoleAtomZs,
                                                        const vector<int>& axisTypes) const
{

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        if (multipoleAtomZs[ii] >= 0 && multipoleAtomXs[ii] >= 0) {
            applyRotationMatrixToParticle(particleData[ii], particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                           multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL, axisTypes[ii]);
        }
    }
}

void AmoebaReferenceMultipoleForce::getAndScaleInverseRs(RealOpenMM dampI, RealOpenMM dampJ,
                                                         RealOpenMM tholeI, RealOpenMM tholeJ,
                                                         RealOpenMM r, vector<RealOpenMM>& rrI) const
{

    RealOpenMM rI             =  1.0/r;
    RealOpenMM r2I            =  rI*rI;

    rrI[0]                    = rI*r2I;
    RealOpenMM constantFactor = 3.0;
    for (unsigned int ii  = 1; ii < rrI.size(); ii++) {
       rrI[ii]         = constantFactor*rrI[ii-1]*r2I;
       constantFactor += 2.0;
    }

    RealOpenMM damp      = dampI*dampJ;
    if (damp != 0.0) {
        RealOpenMM pgamma    = tholeI < tholeJ ? tholeI : tholeJ;
        RealOpenMM ratio     = (r/damp);
                   ratio     = ratio*ratio*ratio;
                   damp      = -pgamma*ratio;

        if (damp > -50.0) {
            RealOpenMM dampExp   = EXP(damp);

            rrI[0]              *= 1.0 - dampExp;
            rrI[1]              *= 1.0 - (1.0 - damp)*dampExp;
            if (rrI.size() > 2) {
                rrI[2]          *= 1.0 - (1.0 - damp + (0.6*damp*damp))*dampExp;
            }
       }
    }
}

void AmoebaReferenceMultipoleForce::calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI,
                                                                        const MultipoleParticleData& particleJ,
                                                                        RealOpenMM dScale, RealOpenMM pScale)
{

    if (particleI.particleIndex == particleJ.particleIndex)
        return;

    RealVec deltaR    = particleJ.position - particleI.position;
    RealOpenMM r      = SQRT(deltaR.dot(deltaR));

    vector<RealOpenMM> rrI(3);

    // get scaling factors, if needed

    getAndScaleInverseRs(particleI.dampingFactor, particleJ.dampingFactor, particleI.thole, particleJ.thole, r, rrI);

    RealOpenMM rr3    = rrI[0];
    RealOpenMM rr5    = rrI[1];
    RealOpenMM rr7    = rrI[2];
    RealOpenMM rr5_2  = 2.0*rr5;

    // field at particle I due multipoles at particle J

    RealVec qDotDelta;
    qDotDelta[0]                                = deltaR[0]*particleJ.quadrupole[QXX] + deltaR[1]*particleJ.quadrupole[QXY] + deltaR[2]*particleJ.quadrupole[QXZ];
    qDotDelta[1]                                = deltaR[0]*particleJ.quadrupole[QXY] + deltaR[1]*particleJ.quadrupole[QYY] + deltaR[2]*particleJ.quadrupole[QYZ];
    qDotDelta[2]                                = deltaR[0]*particleJ.quadrupole[QXZ] + deltaR[1]*particleJ.quadrupole[QYZ] + deltaR[2]*particleJ.quadrupole[QZZ];

    RealOpenMM dipoleDelta                      = particleJ.dipole.dot(deltaR);
    RealOpenMM qdpoleDelta                      = qDotDelta.dot(deltaR);
    RealOpenMM factor                           = rr3*particleJ.charge - rr5*dipoleDelta + rr7*qdpoleDelta;

    RealVec field                               = deltaR*factor + particleJ.dipole*rr3 - qDotDelta*rr5_2;

    unsigned int particleIndex                  = particleI.particleIndex;
    _fixedMultipoleField[particleIndex]        -= field*dScale;
    _fixedMultipoleFieldPolar[particleIndex]   -= field*pScale;

    // field at particle J due multipoles at particle I

    qDotDelta[0]                                = deltaR[0]*particleI.quadrupole[QXX] + deltaR[1]*particleI.quadrupole[QXY] + deltaR[2]*particleI.quadrupole[QXZ];
    qDotDelta[1]                                = deltaR[0]*particleI.quadrupole[QXY] + deltaR[1]*particleI.quadrupole[QYY] + deltaR[2]*particleI.quadrupole[QYZ];
    qDotDelta[2]                                = deltaR[0]*particleI.quadrupole[QXZ] + deltaR[1]*particleI.quadrupole[QYZ] + deltaR[2]*particleI.quadrupole[QZZ];

    dipoleDelta                                 = particleI.dipole.dot(deltaR);
    qdpoleDelta                                 = qDotDelta.dot(deltaR);
    factor                                      = rr3*particleI.charge + rr5*dipoleDelta + rr7*qdpoleDelta;

    field                                       = deltaR*factor - particleI.dipole*rr3 - qDotDelta*rr5_2;
    particleIndex                               = particleJ.particleIndex;
    _fixedMultipoleField[particleIndex]        += field*dScale;
    _fixedMultipoleFieldPolar[particleIndex]   += field*pScale;
}

void AmoebaReferenceMultipoleForce::calculateFixedMultipoleField(const vector<MultipoleParticleData>& particleData)
{

    // calculate fixed multipole fields

    // loop includes diagonal term ii == jj for GK ixn; other calculateFixedMultipoleFieldPairIxn() methods
    // skip calculations for this case

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        for (unsigned int jj = ii; jj < _numParticles; jj++) {

            // if site jj is less than max covalent scaling index then get/apply scaling constants
            // otherwise add unmodified field and fieldPolar to particle fields

            RealOpenMM dScale, pScale;
            if (jj <= _maxScaleIndex[ii]) {
                getDScaleAndPScale(ii, jj, dScale, pScale);
            } else {
                dScale = pScale = 1.0;
            }
            calculateFixedMultipoleFieldPairIxn(particleData[ii], particleData[jj], dScale, pScale);
        }
    }
}

void AmoebaReferenceMultipoleForce::initializeInducedDipoles(vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    // initialize inducedDipoles

    _inducedDipole.resize(_numParticles);
    _inducedDipolePolar.resize(_numParticles);

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        _inducedDipole[ii]       = _fixedMultipoleField[ii];
        _inducedDipolePolar[ii]  = _fixedMultipoleFieldPolar[ii];
    }
}

void AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxn(unsigned int particleI,
                                                                  unsigned int particleJ,
                                                                  RealOpenMM rr3,
                                                                  RealOpenMM rr5,
                                                                  const RealVec& deltaR,
                                                                  const vector<RealVec>& inducedDipole,
                                                                  vector<RealVec>& field) const
{

    RealOpenMM dDotDelta            = rr5*(inducedDipole[particleJ].dot(deltaR));
    field[particleI]               += inducedDipole[particleJ]*rr3 + deltaR*dDotDelta;
    dDotDelta                       = rr5*(inducedDipole[particleI].dot(deltaR));
    field[particleJ]               += inducedDipole[particleI]*rr3 + deltaR*dDotDelta;
}

void AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                                                   const MultipoleParticleData& particleJ,
                                                                   vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

   if (particleI.particleIndex == particleJ.particleIndex)
       return;

    RealVec deltaR       = particleJ.position - particleI.position;
    RealOpenMM r         =  SQRT(deltaR.dot(deltaR));
    vector<RealOpenMM> rrI(2);
    // If we're using the extrapolation algorithm, we need to compute the field gradient, so ask for one more rrI value.
    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated)
        rrI.push_back(0.0);
  
    getAndScaleInverseRs(particleI.dampingFactor, particleJ.dampingFactor,
                          particleI.thole, particleJ.thole, r, rrI);

    RealOpenMM rr3       = -rrI[0];
    RealOpenMM rr5       =  rrI[1];

    for (unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++) {
        calculateInducedDipolePairIxn(particleI.particleIndex, particleJ.particleIndex, rr3, rr5, deltaR,
                                       *updateInducedDipoleFields[ii].inducedDipoles, updateInducedDipoleFields[ii].inducedDipoleField);
        if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
            // Compute and store the field gradient for later use.
            RealOpenMM dx = deltaR[0];
            RealOpenMM dy = deltaR[1];
            RealOpenMM dz = deltaR[2];

            OpenMM::RealVec &dipolesI = (*updateInducedDipoleFields[ii].inducedDipoles)[particleI.particleIndex];
            RealOpenMM xDipole = dipolesI[0];
            RealOpenMM yDipole = dipolesI[1];
            RealOpenMM zDipole = dipolesI[2];
            RealOpenMM muDotR = xDipole*dx + yDipole*dy + zDipole*dz;
            RealOpenMM Exx = muDotR*dx*dx*rrI[2] - (2.0*xDipole*dx + muDotR)*rrI[1];
            RealOpenMM Eyy = muDotR*dy*dy*rrI[2] - (2.0*yDipole*dy + muDotR)*rrI[1];
            RealOpenMM Ezz = muDotR*dz*dz*rrI[2] - (2.0*zDipole*dz + muDotR)*rrI[1];
            RealOpenMM Exy = muDotR*dx*dy*rrI[2] - (xDipole*dy + yDipole*dx)*rrI[1];
            RealOpenMM Exz = muDotR*dx*dz*rrI[2] - (xDipole*dz + zDipole*dx)*rrI[1];
            RealOpenMM Eyz = muDotR*dy*dz*rrI[2] - (yDipole*dz + zDipole*dy)*rrI[1];

            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][0] -= Exx;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][1] -= Eyy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][2] -= Ezz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][3] -= Exy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][4] -= Exz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][5] -= Eyz;

            OpenMM::RealVec &dipolesJ = (*updateInducedDipoleFields[ii].inducedDipoles)[particleJ.particleIndex];
            xDipole = dipolesJ[0];
            yDipole = dipolesJ[1];
            zDipole = dipolesJ[2];
            muDotR = xDipole*dx + yDipole*dy + zDipole*dz;
            Exx = muDotR*dx*dx*rrI[2] - (2.0*xDipole*dx + muDotR)*rrI[1];
            Eyy = muDotR*dy*dy*rrI[2] - (2.0*yDipole*dy + muDotR)*rrI[1];
            Ezz = muDotR*dz*dz*rrI[2] - (2.0*zDipole*dz + muDotR)*rrI[1];
            Exy = muDotR*dx*dy*rrI[2] - (xDipole*dy + yDipole*dx)*rrI[1];
            Exz = muDotR*dx*dz*rrI[2] - (xDipole*dz + zDipole*dx)*rrI[1];
            Eyz = muDotR*dy*dz*rrI[2] - (yDipole*dz + zDipole*dy)*rrI[1];

            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][0] += Exx;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][1] += Eyy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][2] += Ezz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][3] += Exy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][4] += Exz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][5] += Eyz;
        }
    }
}

void AmoebaReferenceMultipoleForce::calculateInducedDipoleFields(const vector<MultipoleParticleData>& particleData, vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields) {
    // Initialize the fields to zero.

    RealVec zeroVec(0.0, 0.0, 0.0);
    for (unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++)
        std::fill(updateInducedDipoleFields[ii].inducedDipoleField.begin(), updateInducedDipoleFields[ii].inducedDipoleField.end(), zeroVec);

    // Add fields from all induced dipoles.

    for (unsigned int ii = 0; ii < particleData.size(); ii++)
        for (unsigned int jj = ii; jj < particleData.size(); jj++)
            calculateInducedDipolePairIxns(particleData[ii], particleData[jj], updateInducedDipoleFields);
}

RealOpenMM AmoebaReferenceMultipoleForce::updateInducedDipoleFields(const vector<MultipoleParticleData>& particleData,
                                                                    vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{
    // Calculate the fields coming from induced dipoles.

    calculateInducedDipoleFields(particleData, updateInducedDipoleFields);

    // Update the induced dipoles and calculate the convergence factor, maxEpsilon

    RealOpenMM maxEpsilon = 0.0;
    for (unsigned int kk = 0; kk < updateInducedDipoleFields.size(); kk++) {
        RealOpenMM epsilon = updateInducedDipole(particleData,
                                                 *updateInducedDipoleFields[kk].fixedMultipoleField,
                                                 updateInducedDipoleFields[kk].inducedDipoleField,
                                                 *updateInducedDipoleFields[kk].inducedDipoles);

        maxEpsilon = epsilon > maxEpsilon ? epsilon : maxEpsilon;
    }

    return maxEpsilon;
}

RealOpenMM AmoebaReferenceMultipoleForce::updateInducedDipole(const vector<MultipoleParticleData>& particleData,
                                                              const vector<RealVec>& fixedMultipoleField,
                                                              const vector<RealVec>& inducedDipoleField,
                                                              vector<RealVec>& inducedDipole)
{

    RealOpenMM epsilon                    = 0.0;
    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        RealVec    oldValue               = inducedDipole[ii];
        RealVec    newValue               = fixedMultipoleField[ii] + inducedDipoleField[ii]*particleData[ii].polarity;
        RealVec    delta                  = newValue - oldValue;
        inducedDipole[ii]                 = oldValue + delta*_polarSOR;
        epsilon                          += delta.dot(delta);
    }
    return epsilon;
}

void AmoebaReferenceMultipoleForce::convergeInduceDipolesBySOR(const vector<MultipoleParticleData>& particleData,
                                                               vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleField)
{

    bool done                 = false;
    setMutualInducedDipoleConverged(false);
    int iteration             = 0;
    RealOpenMM currentEpsilon = 1.0e+50;

    // loop until (1) induced dipoles are converged or
    //            (2) iterations == max iterations or
    //            (3) convergence factor (spsilon) increases

    while (!done) {

        RealOpenMM epsilon = updateInducedDipoleFields(particleData, updateInducedDipoleField);
                   epsilon = _polarSOR*_debye*SQRT(epsilon/(static_cast<RealOpenMM>(_numParticles)));

        if (epsilon < getMutualInducedDipoleTargetEpsilon()) {
            setMutualInducedDipoleConverged(true);
            done = true;
        } else if (currentEpsilon < epsilon || iteration >= getMaximumMutualInducedDipoleIterations()) {
            done = true;
        }

        currentEpsilon = epsilon;
        iteration++;
    }
    setMutualInducedDipoleEpsilon(currentEpsilon);
    setMutualInducedDipoleIterations(iteration);
}


void AmoebaReferenceMultipoleForce::convergeInduceDipolesByExtrapolation(const vector<MultipoleParticleData>& particleData, vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleField) {
    // Start by storing the direct dipoles as PT0

    int numFields = updateInducedDipoleField.size();
    for (int i = 0; i < numFields; i++) {
        UpdateInducedDipoleFieldStruct& field = updateInducedDipoleField[i];
        field.extrapolatedDipoles->resize(_maxPTOrder);
        (*field.extrapolatedDipoles)[0].resize(_numParticles);
        for (int atom = 0; atom < _numParticles; ++atom)
            (*field.extrapolatedDipoles)[0][atom] = (*field.inducedDipoles)[atom];
        field.inducedDipoleFieldGradient.resize(_numParticles);
    }

    // Recursively apply alpha.Tau to the µ_(n) components to generate µ_(n+1), and store the result

    vector<RealOpenMM> zeros(6, 0.0);
    for (int order = 1; order < _maxPTOrder; ++order) {
        for (int i = 0; i < numFields; i++)
            std::fill(updateInducedDipoleField[i].inducedDipoleFieldGradient.begin(), updateInducedDipoleField[i].inducedDipoleFieldGradient.end(), zeros);
        calculateInducedDipoleFields(particleData, updateInducedDipoleField);
        for (int i = 0; i < numFields; i++) {
            UpdateInducedDipoleFieldStruct& field = updateInducedDipoleField[i];
            (*field.extrapolatedDipoles)[order].resize(_numParticles);
            for (int atom = 0; atom < _numParticles; ++atom) {
                (*field.inducedDipoles)[atom] = field.inducedDipoleField[atom] * particleData[atom].polarity;
                (*field.extrapolatedDipoles)[order][atom] = (*field.inducedDipoles)[atom];
            }
            vector<RealOpenMM> fieldGrad(6*_numParticles, 0.0);
            for (int atom = 0; atom < _numParticles; ++atom)
                for (int component = 0; component < 6; ++component)
                    fieldGrad[6*atom + component] = field.inducedDipoleFieldGradient[atom][component];
            field.extrapolatedDipoleFieldGradient->push_back(fieldGrad);
        }
    }

    // Take a linear combination of the µ_(n) components to form the total dipole
    
    for (int i = 0; i < numFields; i++) {
        UpdateInducedDipoleFieldStruct& field = updateInducedDipoleField[i];
        *field.inducedDipoles = vector<RealVec>(_numParticles, RealVec());
        for (int order = 0; order < _maxPTOrder; ++order)
            for (int atom = 0; atom < _numParticles; ++atom)
                (*field.inducedDipoles)[atom] += (*field.extrapolatedDipoles)[order][atom] * _extPartCoefficients[order];
    }
    calculateInducedDipoleFields(particleData, updateInducedDipoleField);
    setMutualInducedDipoleConverged(true);
}

void AmoebaReferenceMultipoleForce::convergeInduceDipolesByDIIS(const vector<MultipoleParticleData>& particleData, vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleField) {
    int numFields = updateInducedDipoleField.size();
    vector<vector<vector<RealVec> > > prevDipoles(numFields);
    vector<vector<RealVec> > prevErrors;
    setMutualInducedDipoleConverged(false);
    int maxPrevious = 20;
    for (int iteration = 0; ; iteration++) {
        // Compute the field from the induced dipoles.

        calculateInducedDipoleFields(particleData, updateInducedDipoleField);

        // Record the current dipoles and the errors in them.

        RealOpenMM maxEpsilon = 0;
        prevErrors.push_back(vector<RealVec>());
        prevErrors.back().resize(_numParticles);
        for (int k = 0; k < numFields; k++) {
            UpdateInducedDipoleFieldStruct& field = updateInducedDipoleField[k];
            prevDipoles[k].push_back(vector<RealVec>());
            prevDipoles[k].back().resize(_numParticles);
            RealOpenMM epsilon = 0;
            for (int i = 0; i < _numParticles; i++) {
                prevDipoles[k].back()[i] = (*field.inducedDipoles)[i];
                RealVec newDipole = (*field.fixedMultipoleField)[i] + field.inducedDipoleField[i]*particleData[i].polarity;
                RealVec error = newDipole-(*field.inducedDipoles)[i];
                prevDipoles[k].back()[i] = newDipole;
                if (k == 0)
                    prevErrors.back()[i] = error;
                epsilon += error.dot(error);
            }
            if (epsilon > maxEpsilon)
                maxEpsilon = epsilon;
        }
        maxEpsilon = _debye*SQRT(maxEpsilon/_numParticles);

        // Decide whether to stop or continue iterating.

        if (maxEpsilon < getMutualInducedDipoleTargetEpsilon())
            setMutualInducedDipoleConverged(true);
        if (maxEpsilon < getMutualInducedDipoleTargetEpsilon() || iteration == getMaximumMutualInducedDipoleIterations()) {
            setMutualInducedDipoleEpsilon(maxEpsilon);
            setMutualInducedDipoleIterations(iteration);
            return;
        }

        // Select the new dipoles.

        if (prevErrors.size() > maxPrevious) {
            prevErrors.erase(prevErrors.begin());
            for (int k = 0; k < numFields; k++)
                prevDipoles[k].erase(prevDipoles[k].begin());
        }
        int numPrevious = prevErrors.size();
        vector<RealOpenMM> coefficients(numPrevious);
        computeDIISCoefficients(prevErrors, coefficients);
        for (int k = 0; k < numFields; k++) {
            UpdateInducedDipoleFieldStruct& field = updateInducedDipoleField[k];
            for (int i = 0; i < _numParticles; i++) {
                RealVec dipole(0.0, 0.0, 0.0);
                for (int j = 0; j < numPrevious; j++)
                    dipole += prevDipoles[k][j][i]*coefficients[j];
                (*field.inducedDipoles)[i] = dipole;
            }
        }
    }

}

void AmoebaReferenceMultipoleForce::computeDIISCoefficients(const vector<vector<RealVec> >& prevErrors, vector<RealOpenMM>& coefficients) const {
    int steps = coefficients.size();
    if (steps == 1) {
        coefficients[0] = 1;
        return;
    }

    // Create the DIIS matrix.

    int rank = steps+1;
    Array2D<double> b(rank, rank);
    b[0][0] = 0;
    for (int i = 0; i < steps; i++)
        b[i+1][0] = b[0][i+1] = -1;
    for (int i = 0; i < steps; i++)
        for (int j = i; j < steps; j++) {
            double sum = 0;
            for (int k = 0; k < _numParticles; k++)
                sum += prevErrors[i][k].dot(prevErrors[j][k]);
            b[i+1][j+1] = b[j+1][i+1] = sum;
        }

    // Solve using SVD.  Since the right hand side is (-1, 0, 0, 0, ...), this is simpler than the general case.

    JAMA::SVD<double> svd(b);
    Array2D<double> u, v;
    svd.getU(u);
    svd.getV(v);
    Array1D<double> s;
    svd.getSingularValues(s);
    int effectiveRank = svd.rank();
    for (int i = 1; i < rank; i++) {
        double d = 0;
        for (int j = 0; j < effectiveRank; j++)
            d -= u[0][j]*v[i][j]/s[j];
        coefficients[i-1] = d;
    }
}

void AmoebaReferenceMultipoleForce::calculateInducedDipoles(const vector<MultipoleParticleData>& particleData)
{

    // calculate fixed electric fields

    zeroFixedMultipoleFields();
    calculateFixedMultipoleField(particleData);

    // initialize inducedDipoles
    // if polarization type is 'Direct', then return after initializing; otherwise
    // converge induced dipoles.

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        _fixedMultipoleField[ii]      *= particleData[ii].polarity;
        _fixedMultipoleFieldPolar[ii] *= particleData[ii].polarity;
    }

    _inducedDipole.resize(_numParticles);
    _inducedDipolePolar.resize(_numParticles);
    vector<UpdateInducedDipoleFieldStruct> updateInducedDipoleField;
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(_fixedMultipoleField, _inducedDipole, _ptDipoleD, _ptDipoleFieldGradientD));
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(_fixedMultipoleFieldPolar, _inducedDipolePolar, _ptDipoleP, _ptDipoleFieldGradientP));

    initializeInducedDipoles(updateInducedDipoleField);

    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Direct) {
        setMutualInducedDipoleConverged(true);
        return;
    }

    // UpdateInducedDipoleFieldStruct contains induced dipole, fixed multipole fields and fields
    // due to other induced dipoles at each site
    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Mutual)
        convergeInduceDipolesByDIIS(particleData, updateInducedDipoleField);
    else if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated)
        convergeInduceDipolesByExtrapolation(particleData, updateInducedDipoleField);
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateElectrostaticPairIxn(const MultipoleParticleData& particleI,
                                                                        const MultipoleParticleData& particleK,
                                                                        const vector<RealOpenMM>& scalingFactors,
                                                                        vector<RealVec>& forces,
                                                                        vector<RealVec>& torque) const
{
    unsigned int iIndex = particleI.particleIndex;
    unsigned int kIndex = particleK.particleIndex;

    RealVec deltaR = particleK.position - particleI.position;
    RealOpenMM r2  = deltaR.dot(deltaR);
    RealOpenMM r   = SQRT(r2);

    // Start by constructing rotation matrices to put dipoles and
    // quadrupoles into the QI frame, from the lab frame.
    RealOpenMM qiRotationMatrix1[3][3];
    formQIRotationMatrix(particleI.position, particleK.position, deltaR, r, qiRotationMatrix1);
    RealOpenMM qiRotationMatrix2[5][5];
    buildSphericalQuadrupoleRotationMatrix(qiRotationMatrix1, qiRotationMatrix2);
    // The force rotation matrix rotates the QI forces into the lab
    // frame, and makes sure the result is in {x,y,z} ordering. Its
    // transpose is used to rotate the induced dipoles to the QI frame.
    RealOpenMM forceRotationMatrix[3][3];
    forceRotationMatrix[0][0] = qiRotationMatrix1[1][1];
    forceRotationMatrix[0][1] = qiRotationMatrix1[2][1];
    forceRotationMatrix[0][2] = qiRotationMatrix1[0][1];
    forceRotationMatrix[1][0] = qiRotationMatrix1[1][2];
    forceRotationMatrix[1][1] = qiRotationMatrix1[2][2];
    forceRotationMatrix[1][2] = qiRotationMatrix1[0][2];
    forceRotationMatrix[2][0] = qiRotationMatrix1[1][0];
    forceRotationMatrix[2][1] = qiRotationMatrix1[2][0];
    forceRotationMatrix[2][2] = qiRotationMatrix1[0][0];
    // For efficiency, we go ahead and cache that transposed version
    // now, because we need to do 4 rotations in total (I,J, and p,d).
    // We also fold in the factor of 0.5 needed to average the p and d
    // components.
    RealOpenMM inducedDipoleRotationMatrix[3][3];
    inducedDipoleRotationMatrix[0][0] = 0.5*qiRotationMatrix1[0][1];
    inducedDipoleRotationMatrix[0][1] = 0.5*qiRotationMatrix1[0][2];
    inducedDipoleRotationMatrix[0][2] = 0.5*qiRotationMatrix1[0][0];
    inducedDipoleRotationMatrix[1][0] = 0.5*qiRotationMatrix1[1][1];
    inducedDipoleRotationMatrix[1][1] = 0.5*qiRotationMatrix1[1][2];
    inducedDipoleRotationMatrix[1][2] = 0.5*qiRotationMatrix1[1][0];
    inducedDipoleRotationMatrix[2][0] = 0.5*qiRotationMatrix1[2][1];
    inducedDipoleRotationMatrix[2][1] = 0.5*qiRotationMatrix1[2][2];
    inducedDipoleRotationMatrix[2][2] = 0.5*qiRotationMatrix1[2][0];

    // Rotate the induced dipoles to the QI frame.
    RealOpenMM qiUindI[3], qiUindJ[3], qiUinpI[3], qiUinpJ[3];
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM valIP = 0.0;
        RealOpenMM valID = 0.0;
        RealOpenMM valJP = 0.0;
        RealOpenMM valJD = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            valIP += inducedDipoleRotationMatrix[ii][jj] * _inducedDipolePolar[iIndex][jj];
            valID += inducedDipoleRotationMatrix[ii][jj] * _inducedDipole[iIndex][jj];
            valJP += inducedDipoleRotationMatrix[ii][jj] * _inducedDipolePolar[kIndex][jj];
            valJD += inducedDipoleRotationMatrix[ii][jj] * _inducedDipole[kIndex][jj];
        }
        qiUindI[ii] = valID;
        qiUinpI[ii] = valIP;
        qiUindJ[ii] = valJD;
        qiUinpJ[ii] = valJP;
    }

    // The Qtilde intermediates (QI frame multipoles) for atoms I and J
    RealOpenMM qiQI[9], qiQJ[9];
    // Rotate the permanent multipoles to the QI frame.
    qiQI[0] = particleI.charge;
    qiQJ[0] = particleK.charge;
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM valI = 0.0;
        RealOpenMM valJ = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            valI += qiRotationMatrix1[ii][jj] * particleI.sphericalDipole[jj];
            valJ += qiRotationMatrix1[ii][jj] * particleK.sphericalDipole[jj];
        }
        qiQI[ii+1] = valI;
        qiQJ[ii+1] = valJ;
    }
    for (int ii = 0; ii < 5; ii++) {
        RealOpenMM valI = 0.0;
        RealOpenMM valJ = 0.0;
        for (int jj = 0; jj < 5; jj++) {
            valI += qiRotationMatrix2[ii][jj] * particleI.sphericalQuadrupole[jj];
            valJ += qiRotationMatrix2[ii][jj] * particleK.sphericalQuadrupole[jj];
        }
        qiQI[ii+4] = valI;
        qiQJ[ii+4] = valJ;
    }

    // The Qtilde{x,y,z} torque intermediates for atoms I and J, which are used to obtain the torques on the permanent moments.
    RealOpenMM qiQIX[9] = {0.0, qiQI[3], 0.0, -qiQI[1], sqrtThree*qiQI[6], qiQI[8], -sqrtThree*qiQI[4] - qiQI[7], qiQI[6], -qiQI[5]};
    RealOpenMM qiQIY[9] = {0.0, -qiQI[2], qiQI[1], 0.0, -sqrtThree*qiQI[5], sqrtThree*qiQI[4] - qiQI[7], -qiQI[8], qiQI[5], qiQI[6]};
    RealOpenMM qiQIZ[9] = {0.0, 0.0, -qiQI[3], qiQI[2], 0.0, -qiQI[6], qiQI[5], -2.0*qiQI[8], 2.0*qiQI[7]};
    RealOpenMM qiQJX[9] = {0.0, qiQJ[3], 0.0, -qiQJ[1], sqrtThree*qiQJ[6], qiQJ[8], -sqrtThree*qiQJ[4] - qiQJ[7], qiQJ[6], -qiQJ[5]};
    RealOpenMM qiQJY[9] = {0.0, -qiQJ[2], qiQJ[1], 0.0, -sqrtThree*qiQJ[5], sqrtThree*qiQJ[4] - qiQJ[7], -qiQJ[8], qiQJ[5], qiQJ[6]};
    RealOpenMM qiQJZ[9] = {0.0, 0.0, -qiQJ[3], qiQJ[2], 0.0, -qiQJ[6], qiQJ[5], -2.0*qiQJ[8], 2.0*qiQJ[7]};

    // The field derivatives at I due to permanent and induced moments on J, and vice-versa.
    // Also, their derivatives w.r.t. R, which are needed for force calculations
    RealOpenMM Vij[9], Vji[9], VjiR[9], VijR[9];
    // The field derivatives at I due to only permanent moments on J, and vice-versa.
    RealOpenMM Vijp[3], Vijd[3], Vjip[3], Vjid[3];
    RealOpenMM rInvVec[7];

    RealOpenMM prefac = (_electric/_dielectric);
    RealOpenMM rInv = 1.0 / r;

    // The rInvVec array is defined such that the ith element is R^-i, with the
    // dieleectric constant folded in, to avoid conversions later.
    rInvVec[1] = prefac * rInv;
    for (int i = 2; i < 7; ++i)
        rInvVec[i] = rInvVec[i-1] * rInv;

    RealOpenMM mScale = scalingFactors[M_SCALE];
    RealOpenMM dScale = scalingFactors[D_SCALE];
    RealOpenMM pScale = scalingFactors[P_SCALE];
    RealOpenMM uScale = scalingFactors[U_SCALE];

    RealOpenMM dmp = particleI.dampingFactor*particleK.dampingFactor;
    RealOpenMM a = particleI.thole < particleK.thole ? particleI.thole : particleK.thole;
    RealOpenMM u = std::abs(dmp) > 1.0E-5 ? r/dmp : 1E10;
    RealOpenMM au3 = a*u*u*u;
    RealOpenMM expau3 = au3 < 50.0 ? EXP(-au3) : 0.0;
    RealOpenMM a2u6 = au3*au3;
    RealOpenMM a3u9 = a2u6*au3;
    // Thole damping factors for energies
    RealOpenMM thole_c  = 1.0 - expau3;
    RealOpenMM thole_d0 = 1.0 - expau3*(1.0 + 1.5*au3);
    RealOpenMM thole_d1 = 1.0 - expau3;
    RealOpenMM thole_q0 = 1.0 - expau3*(1.0 + au3 + a2u6);
    RealOpenMM thole_q1 = 1.0 - expau3*(1.0 + au3);
    // Thole damping factors for derivatives
    RealOpenMM dthole_c  = 1.0 - expau3*(1.0 + 1.5*au3);
    RealOpenMM dthole_d0 = 1.0 - expau3*(1.0 + au3 + 1.5*a2u6);
    RealOpenMM dthole_d1 = 1.0 - expau3*(1.0 + au3);
    RealOpenMM dthole_q0 = 1.0 - expau3*(1.0 + au3 + 0.25*a2u6 + 0.75*a3u9);
    RealOpenMM dthole_q1 = 1.0 - expau3*(1.0 + au3 + 0.75*a2u6);

    // Now we compute the (attenuated) Coulomb operator and its derivatives, contracted with
    // permanent moments and induced dipoles.  Note that the coefficient of the permanent force
    // terms is half of the expected value; this is because we compute the interaction of I with
    // the sum of induced and permanent moments on J, as well as the interaction of J with I's
    // permanent and induced moments; doing so double counts the permanent-permanent interaction.
    RealOpenMM ePermCoef, dPermCoef, eUIndCoef, dUIndCoef, eUInpCoef, dUInpCoef;

    // C-C terms (m=0)
    ePermCoef = rInvVec[1]*mScale;
    dPermCoef = -0.5*mScale*rInvVec[2];
    Vij[0]  = ePermCoef*qiQJ[0];
    Vji[0]  = ePermCoef*qiQI[0];
    VijR[0] = dPermCoef*qiQJ[0];
    VjiR[0] = dPermCoef*qiQI[0];

    // C-D and C-Uind terms (m=0)
    ePermCoef = rInvVec[2]*mScale;
    eUIndCoef = rInvVec[2]*pScale*thole_c;
    eUInpCoef = rInvVec[2]*dScale*thole_c;
    dPermCoef = -rInvVec[3]*mScale;
    dUIndCoef = -2.0*rInvVec[3]*pScale*dthole_c;
    dUInpCoef = -2.0*rInvVec[3]*dScale*dthole_c;
    Vij[0]  += -(ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0]);
    Vji[1]   = -(ePermCoef*qiQI[0]);
    VijR[0] += -(dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0]);
    VjiR[1]  = -(dPermCoef*qiQI[0]);
    Vjip[0]  = -(eUInpCoef*qiQI[0]);
    Vjid[0]  = -(eUIndCoef*qiQI[0]);
    // D-C and Uind-C terms (m=0)
    Vij[1]   = ePermCoef*qiQJ[0];
    Vji[0]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1]  = dPermCoef*qiQJ[0];
    VjiR[0] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0]  = eUInpCoef*qiQJ[0];
    Vijd[0]  = eUIndCoef*qiQJ[0];

    // D-D and D-Uind terms (m=0)
    ePermCoef = -2.0*rInvVec[3]*mScale;
    eUIndCoef = -2.0*rInvVec[3]*pScale*thole_d0;
    eUInpCoef = -2.0*rInvVec[3]*dScale*thole_d0;
    dPermCoef = 3.0*rInvVec[4]*mScale;
    dUIndCoef = 6.0*rInvVec[4]*pScale*dthole_d0;
    dUInpCoef = 6.0*rInvVec[4]*dScale*dthole_d0;
    Vij[1]  += ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0];
    Vji[1]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1] += dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0];
    VjiR[1] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0] += eUInpCoef*qiQJ[1];
    Vijd[0] += eUIndCoef*qiQJ[1];
    Vjip[0] += eUInpCoef*qiQI[1];
    Vjid[0] += eUIndCoef*qiQI[1];
    // D-D and D-Uind terms (m=1)
    ePermCoef = rInvVec[3]*mScale;
    eUIndCoef = rInvVec[3]*pScale*thole_d1;
    eUInpCoef = rInvVec[3]*dScale*thole_d1;
    dPermCoef = -1.5*rInvVec[4]*mScale;
    dUIndCoef = -3.0*rInvVec[4]*pScale*dthole_d1;
    dUInpCoef = -3.0*rInvVec[4]*dScale*dthole_d1;
    Vij[2]  = ePermCoef*qiQJ[2] + eUIndCoef*qiUindJ[1] + eUInpCoef*qiUinpJ[1];
    Vji[2]  = ePermCoef*qiQI[2] + eUIndCoef*qiUindI[1] + eUInpCoef*qiUinpI[1];
    VijR[2] = dPermCoef*qiQJ[2] + dUIndCoef*qiUindJ[1] + dUInpCoef*qiUinpJ[1];
    VjiR[2] = dPermCoef*qiQI[2] + dUIndCoef*qiUindI[1] + dUInpCoef*qiUinpI[1];
    Vij[3]  = ePermCoef*qiQJ[3] + eUIndCoef*qiUindJ[2] + eUInpCoef*qiUinpJ[2];
    Vji[3]  = ePermCoef*qiQI[3] + eUIndCoef*qiUindI[2] + eUInpCoef*qiUinpI[2];
    VijR[3] = dPermCoef*qiQJ[3] + dUIndCoef*qiUindJ[2] + dUInpCoef*qiUinpJ[2];
    VjiR[3] = dPermCoef*qiQI[3] + dUIndCoef*qiUindI[2] + dUInpCoef*qiUinpI[2];
    Vijp[1] = eUInpCoef*qiQJ[2];
    Vijd[1] = eUIndCoef*qiQJ[2];
    Vjip[1] = eUInpCoef*qiQI[2];
    Vjid[1] = eUIndCoef*qiQI[2];
    Vijp[2] = eUInpCoef*qiQJ[3];
    Vijd[2] = eUIndCoef*qiQJ[3];
    Vjip[2] = eUInpCoef*qiQI[3];
    Vjid[2] = eUIndCoef*qiQI[3];

    // C-Q terms (m=0)
    ePermCoef = mScale*rInvVec[3];
    dPermCoef = -1.5*rInvVec[4]*mScale;
    Vij[0]  += ePermCoef*qiQJ[4];
    Vji[4]   = ePermCoef*qiQI[0];
    VijR[0] += dPermCoef*qiQJ[4];
    VjiR[4]  = dPermCoef*qiQI[0];
    // Q-C terms (m=0)
    Vij[4]   = ePermCoef*qiQJ[0];
    Vji[0]  += ePermCoef*qiQI[4];
    VijR[4]  = dPermCoef*qiQJ[0];
    VjiR[0] += dPermCoef*qiQI[4];

    // D-Q and Uind-Q terms (m=0)
    ePermCoef = rInvVec[4]*3.0*mScale;
    eUIndCoef = rInvVec[4]*3.0*pScale*thole_q0;
    eUInpCoef = rInvVec[4]*3.0*dScale*thole_q0;
    dPermCoef = -6.0*rInvVec[5]*mScale;
    dUIndCoef = -12.0*rInvVec[5]*pScale*dthole_q0;
    dUInpCoef = -12.0*rInvVec[5]*dScale*dthole_q0;
    Vij[1]  += ePermCoef*qiQJ[4];
    Vji[4]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1] += dPermCoef*qiQJ[4];
    VjiR[4] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0] += eUInpCoef*qiQJ[4];
    Vijd[0] += eUIndCoef*qiQJ[4];
    // Q-D and Q-Uind terms (m=0)
    Vij[4]  += -(ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0]);
    Vji[1]  += -(ePermCoef*qiQI[4]);
    VijR[4] += -(dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0]);
    VjiR[1] += -(dPermCoef*qiQI[4]);
    Vjip[0] += -(eUInpCoef*qiQI[4]);
    Vjid[0] += -(eUIndCoef*qiQI[4]);

    // D-Q and Uind-Q terms (m=1)
    ePermCoef = -sqrtThree*rInvVec[4]*mScale;
    eUIndCoef = -sqrtThree*rInvVec[4]*pScale*thole_q1;
    eUInpCoef = -sqrtThree*rInvVec[4]*dScale*thole_q1;
    dPermCoef = 2.0*sqrtThree*rInvVec[5]*mScale;
    dUIndCoef = 4.0*sqrtThree*rInvVec[5]*pScale*dthole_q1;
    dUInpCoef = 4.0*sqrtThree*rInvVec[5]*dScale*dthole_q1;
    Vij[2]  += ePermCoef*qiQJ[5];
    Vji[5]   = ePermCoef*qiQI[2] + eUIndCoef*qiUindI[1] + eUInpCoef*qiUinpI[1];
    VijR[2] += dPermCoef*qiQJ[5];
    VjiR[5]  = dPermCoef*qiQI[2] + dUIndCoef*qiUindI[1] + dUInpCoef*qiUinpI[1];
    Vij[3]  += ePermCoef*qiQJ[6];
    Vji[6]   = ePermCoef*qiQI[3] + eUIndCoef*qiUindI[2] + eUInpCoef*qiUinpI[2];
    VijR[3] += dPermCoef*qiQJ[6];
    VjiR[6]  = dPermCoef*qiQI[3] + dUIndCoef*qiUindI[2] + dUInpCoef*qiUinpI[2];
    Vijp[1] += eUInpCoef*qiQJ[5];
    Vijd[1] += eUIndCoef*qiQJ[5];
    Vijp[2] += eUInpCoef*qiQJ[6];
    Vijd[2] += eUIndCoef*qiQJ[6];
    // D-Q and Uind-Q terms (m=1)
    Vij[5]   = -(ePermCoef*qiQJ[2] + eUIndCoef*qiUindJ[1] + eUInpCoef*qiUinpJ[1]);
    Vji[2]  += -(ePermCoef*qiQI[5]);
    VijR[5]  = -(dPermCoef*qiQJ[2] + dUIndCoef*qiUindJ[1] + dUInpCoef*qiUinpJ[1]);
    VjiR[2] += -(dPermCoef*qiQI[5]);
    Vij[6]   = -(ePermCoef*qiQJ[3] + eUIndCoef*qiUindJ[2] + eUInpCoef*qiUinpJ[2]);
    Vji[3]  += -(ePermCoef*qiQI[6]);
    VijR[6]  = -(dPermCoef*qiQJ[3] + dUIndCoef*qiUindJ[2] + dUInpCoef*qiUinpJ[2]);
    VjiR[3] += -(dPermCoef*qiQI[6]);
    Vjip[1] += -(eUInpCoef*qiQI[5]);
    Vjid[1] += -(eUIndCoef*qiQI[5]);
    Vjip[2] += -(eUInpCoef*qiQI[6]);
    Vjid[2] += -(eUIndCoef*qiQI[6]);

    // Q-Q terms (m=0)
    ePermCoef = 6.0*rInvVec[5]*mScale;
    dPermCoef = -15.0*rInvVec[6]*mScale;
    Vij[4]  += ePermCoef*qiQJ[4];
    Vji[4]  += ePermCoef*qiQI[4];
    VijR[4] += dPermCoef*qiQJ[4];
    VjiR[4] += dPermCoef*qiQI[4];
    // Q-Q terms (m=1)
    ePermCoef = -4.0*rInvVec[5]*mScale;
    dPermCoef = 10.0*rInvVec[6]*mScale;
    Vij[5]  += ePermCoef*qiQJ[5];
    Vji[5]  += ePermCoef*qiQI[5];
    VijR[5] += dPermCoef*qiQJ[5];
    VjiR[5] += dPermCoef*qiQI[5];
    Vij[6]  += ePermCoef*qiQJ[6];
    Vji[6]  += ePermCoef*qiQI[6];
    VijR[6] += dPermCoef*qiQJ[6];
    VjiR[6] += dPermCoef*qiQI[6];
    // Q-Q terms (m=2)
    ePermCoef = rInvVec[5]*mScale;
    dPermCoef = -2.5*rInvVec[6]*mScale;
    Vij[7]  = ePermCoef*qiQJ[7];
    Vji[7]  = ePermCoef*qiQI[7];
    VijR[7] = dPermCoef*qiQJ[7];
    VjiR[7] = dPermCoef*qiQI[7];
    Vij[8]  = ePermCoef*qiQJ[8];
    Vji[8]  = ePermCoef*qiQI[8];
    VijR[8] = dPermCoef*qiQJ[8];
    VjiR[8] = dPermCoef*qiQI[8];

    // Evaluate the energies, forces and torques due to permanent+induced moments
    // interacting with just the permanent moments.
    RealOpenMM energy = 0.5*(qiQI[0]*Vij[0] + qiQJ[0]*Vji[0]);
    RealOpenMM fIZ = qiQI[0]*VijR[0];
    RealOpenMM fJZ = qiQJ[0]*VjiR[0];
    RealOpenMM EIX = 0.0, EIY = 0.0, EIZ = 0.0, EJX = 0.0, EJY = 0.0, EJZ = 0.0;
    for (int i = 1; i < 9; ++i) {
        energy += 0.5*(qiQI[i]*Vij[i] + qiQJ[i]*Vji[i]);
        fIZ += qiQI[i]*VijR[i];
        fJZ += qiQJ[i]*VjiR[i];
        EIX += qiQIX[i]*Vij[i];
        EIY += qiQIY[i]*Vij[i];
        EIZ += qiQIZ[i]*Vij[i];
        EJX += qiQJX[i]*Vji[i];
        EJY += qiQJY[i]*Vji[i];
        EJZ += qiQJZ[i]*Vji[i];
    }

    // Define the torque intermediates for the induced dipoles. These are simply the induced dipole torque
    // intermediates dotted with the field due to permanent moments only, at each center. We inline the
    // induced dipole torque intermediates here, for simplicity. N.B. There are no torques on the dipoles
    // themselves, so we accumulate the torque intermediates into separate variables to allow them to be
    // used only in the force calculation.
    //
    // The torque about the x axis (needed to obtain the y force on the induced dipoles, below)
    //    qiUindIx[0] = qiQUindI[2];    qiUindIx[1] = 0;    qiUindIx[2] = -qiQUindI[0]
    RealOpenMM iEIX = qiUinpI[2]*Vijp[0] + qiUindI[2]*Vijd[0] - qiUinpI[0]*Vijp[2] - qiUindI[0]*Vijd[2];
    RealOpenMM iEJX = qiUinpJ[2]*Vjip[0] + qiUindJ[2]*Vjid[0] - qiUinpJ[0]*Vjip[2] - qiUindJ[0]*Vjid[2];
    // The torque about the y axis (needed to obtain the x force on the induced dipoles, below)
    //    qiUindIy[0] = -qiQUindI[1];   qiUindIy[1] = qiQUindI[0];    qiUindIy[2] = 0
    RealOpenMM iEIY = qiUinpI[0]*Vijp[1] + qiUindI[0]*Vijd[1] - qiUinpI[1]*Vijp[0] - qiUindI[1]*Vijd[0];
    RealOpenMM iEJY = qiUinpJ[0]*Vjip[1] + qiUindJ[0]*Vjid[1] - qiUinpJ[1]*Vjip[0] - qiUindJ[1]*Vjid[0];

    // Add in the induced-induced terms, if needed.
    if(getPolarizationType() == AmoebaReferenceMultipoleForce::Mutual) {
        // Uind-Uind terms (m=0)
        RealOpenMM eCoef = -4.0*rInvVec[3]*uScale*thole_d0;
        RealOpenMM dCoef = 6.0*rInvVec[4]*uScale*dthole_d0;
        iEIX += eCoef*(qiUinpI[2]*qiUindJ[0] + qiUindI[2]*qiUinpJ[0]);
        iEJX += eCoef*(qiUinpJ[2]*qiUindI[0] + qiUindJ[2]*qiUinpI[0]);
        iEIY -= eCoef*(qiUinpI[1]*qiUindJ[0] + qiUindI[1]*qiUinpJ[0]);
        iEJY -= eCoef*(qiUinpJ[1]*qiUindI[0] + qiUindJ[1]*qiUinpI[0]);
        fIZ += dCoef*(qiUinpI[0]*qiUindJ[0] + qiUindI[0]*qiUinpJ[0]);
        fIZ += dCoef*(qiUinpJ[0]*qiUindI[0] + qiUindJ[0]*qiUinpI[0]);
        // Uind-Uind terms (m=1)
        eCoef = 2.0*rInvVec[3]*uScale*thole_d1;
        dCoef = -3.0*rInvVec[4]*uScale*dthole_d1;
        iEIX -= eCoef*(qiUinpI[0]*qiUindJ[2] + qiUindI[0]*qiUinpJ[2]);
        iEJX -= eCoef*(qiUinpJ[0]*qiUindI[2] + qiUindJ[0]*qiUinpI[2]);
        iEIY += eCoef*(qiUinpI[0]*qiUindJ[1] + qiUindI[0]*qiUinpJ[1]);
        iEJY += eCoef*(qiUinpJ[0]*qiUindI[1] + qiUindJ[0]*qiUinpI[1]);
        fIZ += dCoef*(qiUinpI[1]*qiUindJ[1] + qiUindI[1]*qiUinpJ[1] + qiUinpI[2]*qiUindJ[2] + qiUindI[2]*qiUinpJ[2]);
        fIZ += dCoef*(qiUinpJ[1]*qiUindI[1] + qiUindJ[1]*qiUinpI[1] + qiUinpJ[2]*qiUindI[2] + qiUindJ[2]*qiUinpI[2]);
    }

    // The quasi-internal frame forces and torques.  Note that the induced torque intermediates are
    // used in the force expression, but not in the torques; the induced dipoles are isotropic.
    RealOpenMM qiForce[3] = {rInv*(EIY+EJY+iEIY+iEJY), -rInv*(EIX+EJX+iEIX+iEJX), -(fJZ+fIZ)};
    RealOpenMM qiTorqueI[3] = {-EIX, -EIY, -EIZ};
    RealOpenMM qiTorqueJ[3] = {-EJX, -EJY, -EJZ};

    // Rotate the forces and torques back to the lab frame
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM forceVal = 0.0;
        RealOpenMM torqueIVal = 0.0;
        RealOpenMM torqueJVal = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            forceVal   += forceRotationMatrix[ii][jj] * qiForce[jj];
            torqueIVal += forceRotationMatrix[ii][jj] * qiTorqueI[jj];
            torqueJVal += forceRotationMatrix[ii][jj] * qiTorqueJ[jj];
        }
        torque[iIndex][ii] += torqueIVal;
        torque[kIndex][ii] += torqueJVal;
        forces[iIndex][ii] -= forceVal;
        forces[kIndex][ii] += forceVal;
    }
    return energy;
}

void AmoebaReferenceMultipoleForce::mapTorqueToForceForParticle(const MultipoleParticleData& particleI,
                                                                const MultipoleParticleData& particleU,
                                                                const MultipoleParticleData& particleV,
                                                                      MultipoleParticleData* particleW,
                                                                      int axisType, const Vec3& torque,
                                                                      vector<RealVec>& forces) const
{

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

    static const int X                  = 0;
    static const int Y                  = 1;
    static const int Z                  = 2;
    static const int I                  = 3;

    RealOpenMM norms[LastVectorIndex];
    RealOpenMM angles[LastVectorIndex][2];

    // ---------------------------------------------------------------------------------------

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    if (axisType == AmoebaMultipoleForce::NoAxisType) {
        return;
    }

    RealVec vectorU = particleU.position - particleI.position;
    norms[U]        = normalizeRealVec(vectorU);

    RealVec vectorV = particleV.position - particleI.position;
    norms[V]        = normalizeRealVec(vectorV);

    RealVec vectorW;
    if (particleW && (axisType == AmoebaMultipoleForce::ZBisect || axisType == AmoebaMultipoleForce::ThreeFold)) {
         vectorW = particleW->position - particleI.position;
    } else {
         vectorW = vectorU.cross(vectorV);
    }
    norms[W]  = normalizeRealVec(vectorW);

    RealVec vectorUV, vectorUW, vectorVW;
    vectorUV = vectorV.cross(vectorU);
    vectorUW = vectorW.cross(vectorU);
    vectorVW = vectorW.cross(vectorV);

    norms[UV]                     = normalizeRealVec(vectorUV);
    norms[UW]                     = normalizeRealVec(vectorUW);
    norms[VW]                     = normalizeRealVec(vectorVW);

    // angles[][0] is cosine of angle
    // angles[][1] is sine   of angle

    angles[UV][0]                 = vectorU.dot(vectorV);
    angles[UV][1]                 = SQRT(1.0 - angles[UV][0]*angles[UV][0]);

    angles[UW][0]                 = vectorU.dot(vectorW);
    angles[UW][1]                 = SQRT(1.0 - angles[UW][0]*angles[UW][0]);

    angles[VW][0]                 = vectorV.dot(vectorW);
    angles[VW][1]                 = SQRT(1.0 - angles[VW][0]*angles[VW][0]);

    RealVec dphi;
    dphi[U]                       = vectorU.dot(torque);
    dphi[V]                       = vectorV.dot(torque);
    dphi[W]                       = vectorW.dot(torque);
    dphi                         *= -1.0;

    // branch based on axis type

    if (axisType == AmoebaMultipoleForce::ZThenX || axisType == AmoebaMultipoleForce::Bisector) {

        RealOpenMM factor1;
        RealOpenMM factor2;
        RealOpenMM factor3;
        RealOpenMM factor4;
        RealOpenMM half = 0.5;

        factor1                 =  dphi[V]/(norms[U]*angles[UV][1]);
        factor2                 =  dphi[W]/(norms[U]);
        factor3                 = -dphi[U]/(norms[V]*angles[UV][1]);

        if (axisType == AmoebaMultipoleForce::Bisector) {
            factor2    *= half;
            factor4     = half*dphi[W]/(norms[V]);
        } else {
            factor4     = 0.0;
        }

        for (int ii = 0; ii < 3; ii++) {
            double forceU                                        =  vectorUV[ii]*factor1 + factor2*vectorUW[ii];
            forces[particleU.particleIndex][ii]                 -=  forceU;

            double forceV                                        =  vectorUV[ii]*factor3 + factor4*vectorVW[ii];
            forces[particleV.particleIndex][ii]                 -=  forceV;

            forces[particleI.particleIndex][ii]                 +=  (forceU + forceV);
        }

    } else if (axisType == AmoebaMultipoleForce::ZBisect) {

        RealVec vectorR           = vectorV + vectorW;
        RealVec vectorS           = vectorU.cross(vectorR);

        norms[R]                  = normalizeRealVec(vectorR);
        norms[S]                  = normalizeRealVec(vectorS);

        RealVec vectorUR          =  vectorR.cross(vectorU);
        RealVec vectorUS          =  vectorS.cross(vectorU);
        RealVec vectorVS          =  vectorS.cross(vectorV);
        RealVec vectorWS          =  vectorS.cross(vectorW);

        norms[UR]                 = normalizeRealVec(vectorUR);
        norms[US]                 = normalizeRealVec(vectorUS);
        norms[VS]                 = normalizeRealVec(vectorVS);
        norms[WS]                 = normalizeRealVec(vectorWS);

        angles[UR][0]             = vectorU.dot(vectorR);
        angles[UR][1]             = SQRT(1.0 - angles[UR][0]*angles[UR][0]);

        angles[US][0]             = vectorU.dot(vectorS);
        angles[US][1]             = SQRT(1.0 - angles[US][0]*angles[US][0]);

        angles[VS][0]             = vectorV.dot(vectorS);
        angles[VS][1]             = SQRT(1.0 - angles[VS][0]*angles[VS][0]);

        angles[WS][0]             = vectorW.dot(vectorS);
        angles[WS][1]             = SQRT(1.0 - angles[WS][0]*angles[WS][0]);

        RealVec t1                = vectorV - vectorS*angles[VS][0];
        RealVec t2                = vectorW - vectorS*angles[WS][0];

        RealOpenMM notUsed        = normalizeRealVec(t1);
              notUsed             = normalizeRealVec(t2);

        RealOpenMM ut1cos         = vectorU.dot(t1);
        RealOpenMM ut1sin         = SQRT(1.0 - ut1cos*ut1cos);

        RealOpenMM ut2cos         = vectorU.dot(t2);
        RealOpenMM ut2sin         = SQRT(1.0 - ut2cos*ut2cos);

        RealOpenMM dphiR          = vectorR.dot(torque)*(-1.0);
        RealOpenMM dphiS          = vectorS.dot(torque)*(-1.0);

        RealOpenMM factor1        = dphiR/(norms[U]*angles[UR][1]);
        RealOpenMM factor2        = dphiS/(norms[U]);
        RealOpenMM factor3        = dphi[U]/(norms[V]*(ut1sin+ut2sin));
        RealOpenMM factor4        = dphi[U]/(norms[W]*(ut1sin+ut2sin));

        RealVec forceU            =  vectorUR*factor1 + vectorUS*factor2;
        forces[particleU.particleIndex]        -= forceU;

        RealVec forceV            = (vectorS*angles[VS][1] - t1*angles[VS][0])*factor3;
        forces[particleV.particleIndex]        -= forceV;

        RealVec forceW            = (vectorS*angles[WS][1] - t2*angles[WS][0])*factor4;
        forces[particleW->particleIndex]       -= forceW;

        forces[particleI.particleIndex]        += (forceU + forceV + forceW);

    } else if (axisType == AmoebaMultipoleForce::ThreeFold) {

        // 3-fold

        for (int ii = 0; ii < 3; ii++) {

            RealOpenMM du =  vectorUW[ii]*dphi[W]/(norms[U]*angles[UW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]) -
                             vectorUW[ii]*dphi[U]/(norms[U]*angles[UW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[U]*angles[UV][1]);

            RealOpenMM dv =  vectorVW[ii]*dphi[W]/(norms[V]*angles[VW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[V]*angles[UV][1]) -
                             vectorVW[ii]*dphi[V]/(norms[V]*angles[VW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[V]*angles[UV][1]);

            RealOpenMM dw = -vectorUW[ii]*dphi[U]/(norms[W]*angles[UW][1]) -
                             vectorVW[ii]*dphi[V]/(norms[W]*angles[VW][1]) +
                             vectorUW[ii]*dphi[W]/(norms[W]*angles[UW][1]) +
                             vectorVW[ii]*dphi[W]/(norms[W]*angles[VW][1]);

            du /= 3.0;
            dv /= 3.0;
            dw /= 3.0;

            forces[particleU.particleIndex][ii] -= du;
            forces[particleV.particleIndex][ii] -= dv;
            if (particleW)
                forces[particleW->particleIndex][ii] -= dw;
            forces[particleI.particleIndex][ii] += (du + dv + dw);
        }

    } else if (axisType == AmoebaMultipoleForce::ZOnly) {

        // z-only

        for (int ii = 0; ii < 3; ii++) {
            RealOpenMM du                               = vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]);
            forces[particleU.particleIndex][ii]        -= du;
            forces[particleI.particleIndex][ii]        += du;
        }
    }
}

void AmoebaReferenceMultipoleForce::mapTorqueToForce(vector<MultipoleParticleData>& particleData,
                                                     const vector<int>& multipoleAtomXs,
                                                     const vector<int>& multipoleAtomYs,
                                                     const vector<int>& multipoleAtomZs,
                                                     const vector<int>& axisTypes,
                                                     vector<RealVec>& torques,
                                                     vector<RealVec>& forces) const
{

    // map torques to forces

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        if (axisTypes[ii] != AmoebaMultipoleForce::NoAxisType) {
             mapTorqueToForceForParticle(particleData[ii],
                                         particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                         multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL,
                                         axisTypes[ii], torques[ii], forces);
        }
    }
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateElectrostatic(const vector<MultipoleParticleData>& particleData,
                                                                 vector<RealVec>& torques,
                                                                 vector<RealVec>& forces)
{
    RealOpenMM energy = 0.0;
    vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for (unsigned int kk = 0; kk < scaleFactors.size(); kk++) {
        scaleFactors[kk] = 1.0;
    }

    // main loop over particle pairs

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii+1; jj < particleData.size(); jj++) {

            if (jj <= _maxScaleIndex[ii]) {
                getMultipoleScaleFactors(ii, jj, scaleFactors);
            }

            energy += calculateElectrostaticPairIxn(particleData[ii], particleData[jj], scaleFactors, forces, torques);

            if (jj <= _maxScaleIndex[ii]) {
                for (unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++) {
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }
    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
        RealOpenMM prefac = (_electric/_dielectric);
        for (int i = 0; i < _numParticles; i++) {
            // Compute the µ(m) T µ(n) force contributions here
            for (int l = 0; l < _maxPTOrder-1; ++l) {
                for (int m = 0; m < _maxPTOrder-1-l; ++m) {
                    RealOpenMM p = _extPartCoefficients[l+m+1];
                    if(std::fabs(p) < 1e-6) continue;
                    forces[i][0] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+0]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+1]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+4]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+5]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+2]);
                    forces[i][0] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+0]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+1]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+4]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+5]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+2]);
                }
            }
        }
    }

    return energy;
}

void AmoebaReferenceMultipoleForce::setup(const vector<RealVec>& particlePositions,
                                          const vector<RealOpenMM>& charges,
                                          const vector<RealOpenMM>& dipoles,
                                          const vector<RealOpenMM>& quadrupoles,
                                          const vector<RealOpenMM>& tholes,
                                          const vector<RealOpenMM>& dampingFactors,
                                          const vector<RealOpenMM>& polarity,
                                          const vector<int>& axisTypes,
                                          const vector<int>& multipoleAtomZs,
                                          const vector<int>& multipoleAtomXs,
                                          const vector<int>& multipoleAtomYs,
                                          const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                          vector<MultipoleParticleData>& particleData)
{


    // load particle parameters into vector of MultipoleParticleData
    // check for inverted chiral centers
    // apply rotation matrix to get lab frame dipole and quadrupoles
    // setup scaling factors
    // get induced dipoles
    // check if induced dipoles converged

    _numParticles = particlePositions.size();
    loadParticleData(particlePositions, charges, dipoles, quadrupoles,
                      tholes, dampingFactors, polarity, particleData);

    checkChiral(particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes);

    applyRotationMatrix(particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes);

    setupScaleMaps(multipoleAtomCovalentInfo);

    calculateInducedDipoles(particleData);

    if (!getMutualInducedDipoleConverged()) {
        std::stringstream message;
        message << "Induced dipoles did not converge: ";
        message << " iterations="      << getMutualInducedDipoleIterations();
        message << " eps="             << getMutualInducedDipoleEpsilon();
        throw OpenMMException(message.str());
    }
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateForceAndEnergy(const vector<RealVec>& particlePositions,
                                                                  const vector<RealOpenMM>& charges,
                                                                  const vector<RealOpenMM>& dipoles,
                                                                  const vector<RealOpenMM>& quadrupoles,
                                                                  const vector<RealOpenMM>& tholes,
                                                                  const vector<RealOpenMM>& dampingFactors,
                                                                  const vector<RealOpenMM>& polarity,
                                                                  const vector<int>& axisTypes,
                                                                  const vector<int>& multipoleAtomZs,
                                                                  const vector<int>& multipoleAtomXs,
                                                                  const vector<int>& multipoleAtomYs,
                                                                  const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                                  vector<RealVec>& forces)
{

    // setup, including calculating induced dipoles
    // calculate electrostatic ixns including torques
    // map torques to forces

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);

    vector<RealVec> torques;
    initializeRealVecVector(torques);
    RealOpenMM energy = calculateElectrostatic(particleData, torques, forces);

    mapTorqueToForce(particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes, torques, forces);

    return energy;
}

void AmoebaReferenceMultipoleForce::calculateInducedDipoles(const vector<RealVec>& particlePositions,
                                                            const vector<RealOpenMM>& charges,
                                                            const vector<RealOpenMM>& dipoles,
                                                            const vector<RealOpenMM>& quadrupoles,
                                                            const vector<RealOpenMM>& tholes,
                                                            const vector<RealOpenMM>& dampingFactors,
                                                            const vector<RealOpenMM>& polarity,
                                                            const vector<int>& axisTypes,
                                                            const vector<int>& multipoleAtomZs,
                                                            const vector<int>& multipoleAtomXs,
                                                            const vector<int>& multipoleAtomYs,
                                                            const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                            vector<RealVec>& outputInducedDipoles) {
    // setup, including calculating induced dipoles

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);
    outputInducedDipoles = _inducedDipole;
}




void AmoebaReferenceMultipoleForce::calculateLabFramePermanentDipoles(const vector<RealVec>& particlePositions,
                                                                      const vector<RealOpenMM>& charges,
                                                                      const vector<RealOpenMM>& dipoles,
                                                                      const vector<RealOpenMM>& quadrupoles,
                                                                      const vector<RealOpenMM>& tholes,
                                                                      const vector<RealOpenMM>& dampingFactors,
                                                                      const vector<RealOpenMM>& polarity,
                                                                      const vector<int>& axisTypes,
                                                                      const vector<int>& multipoleAtomZs,
                                                                      const vector<int>& multipoleAtomXs,
                                                                      const vector<int>& multipoleAtomYs,
                                                                      const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                                      vector<RealVec>& outputRotatedPermanentDipoles) {
    // setup, including calculating permanent dipoles

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);
    outputRotatedPermanentDipoles.resize(_numParticles);
    for (int i = 0; i < _numParticles; i++)
      {
      outputRotatedPermanentDipoles[i] = particleData[i].dipole;
      }
}

void AmoebaReferenceMultipoleForce::calculateTotalDipoles(const vector<RealVec>& particlePositions,
                                                          const vector<RealOpenMM>& charges,
                                                          const vector<RealOpenMM>& dipoles,
                                                          const vector<RealOpenMM>& quadrupoles,
                                                          const vector<RealOpenMM>& tholes,
                                                          const vector<RealOpenMM>& dampingFactors,
                                                          const vector<RealOpenMM>& polarity,
                                                          const vector<int>& axisTypes,
                                                          const vector<int>& multipoleAtomZs,
                                                          const vector<int>& multipoleAtomXs,
                                                          const vector<int>& multipoleAtomYs,
                                                          const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                          vector<RealVec>& outputTotalDipoles) {
    // setup, including calculating permanent dipoles

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);
    outputTotalDipoles.resize(_numParticles);
    for (int i = 0; i < _numParticles; i++)
      {
      for (int j = 0; j < 3; j++)
        {
        outputTotalDipoles[i][j] = particleData[i].dipole[j] + _inducedDipole[i][j];
        }
      }
}

void AmoebaReferenceMultipoleForce::calculateAmoebaSystemMultipoleMoments(const vector<RealOpenMM>& masses,
                                                                          const vector<RealVec>& particlePositions,
                                                                          const vector<RealOpenMM>& charges,
                                                                          const vector<RealOpenMM>& dipoles,
                                                                          const vector<RealOpenMM>& quadrupoles,
                                                                          const vector<RealOpenMM>& tholes,
                                                                          const vector<RealOpenMM>& dampingFactors,
                                                                          const vector<RealOpenMM>& polarity,
                                                                          const vector<int>& axisTypes,
                                                                          const vector<int>& multipoleAtomZs,
                                                                          const vector<int>& multipoleAtomXs,
                                                                          const vector<int>& multipoleAtomYs,
                                                                          const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                                          vector<RealOpenMM>& outputMultipoleMoments)
{

    // setup, including calculating induced dipoles
    // remove center of mass
    // calculate system moments

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);

    RealOpenMM totalMass  = 0.0;
    RealVec centerOfMass  = RealVec(0.0, 0.0, 0.0);
    for (unsigned int ii  = 0; ii < _numParticles; ii++) {
        RealOpenMM mass   = masses[ii];
        totalMass        += mass;
        centerOfMass     += particleData[ii].position*mass;
    }
    vector<RealVec> localPositions(_numParticles);
    if (totalMass > 0.0) {
        centerOfMass  *= 1.0/totalMass;
    }
    for (unsigned int ii  = 0; ii < _numParticles; ii++) {
        localPositions[ii] = particleData[ii].position - centerOfMass;
    }

    RealOpenMM netchg  = 0.0;

    RealVec dpl        = RealVec(0.0, 0.0, 0.0);

    RealOpenMM xxqdp   = 0.0;
    RealOpenMM xyqdp   = 0.0;
    RealOpenMM xzqdp   = 0.0;

    RealOpenMM yyqdp   = 0.0;
    RealOpenMM yzqdp   = 0.0;

    RealOpenMM zzqdp   = 0.0;

    for (unsigned int ii  = 0; ii < _numParticles; ii++) {

        RealOpenMM charge         = particleData[ii].charge;
        RealVec position          = localPositions[ii];
        netchg                   += charge;

        RealVec netDipole         = (particleData[ii].dipole  + _inducedDipole[ii]);

        dpl                      += position*charge + netDipole;

        xxqdp                    += position[0]*position[0]*charge + 2.0*position[0]*netDipole[0];
        xyqdp                    += position[0]*position[1]*charge + position[0]*netDipole[1] + position[1]*netDipole[0];
        xzqdp                    += position[0]*position[2]*charge + position[0]*netDipole[2] + position[2]*netDipole[0];

        yyqdp                    += position[1]*position[1]*charge + 2.0*position[1]*netDipole[1];
        yzqdp                    += position[1]*position[2]*charge + position[1]*netDipole[2] + position[2]*netDipole[1];

        zzqdp                    += position[2]*position[2]*charge + 2.0*position[2]*netDipole[2];

    }

    // convert the quadrupole from traced to traceless form

    outputMultipoleMoments.resize(13);
    RealOpenMM qave                  = (xxqdp + yyqdp + zzqdp)/3.0;
    outputMultipoleMoments[4]        = 0.5*(xxqdp-qave);
    outputMultipoleMoments[5]        = 0.5*xyqdp;
    outputMultipoleMoments[6]        = 0.5*xzqdp;
    outputMultipoleMoments[8]        = 0.5*(yyqdp-qave);
    outputMultipoleMoments[9]        = 0.5*yzqdp;
    outputMultipoleMoments[12]       = 0.5*(zzqdp-qave);

    // add the traceless atomic quadrupoles to total quadrupole

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        outputMultipoleMoments[4]  += particleData[ii].quadrupole[QXX];
        outputMultipoleMoments[5]  += particleData[ii].quadrupole[QXY];
        outputMultipoleMoments[6]  += particleData[ii].quadrupole[QXZ];
        outputMultipoleMoments[8]  += particleData[ii].quadrupole[QYY];
        outputMultipoleMoments[9]  += particleData[ii].quadrupole[QYZ];
        outputMultipoleMoments[12] += particleData[ii].quadrupole[QZZ];
    }
    outputMultipoleMoments[7]  = outputMultipoleMoments[5];
    outputMultipoleMoments[10] = outputMultipoleMoments[6];
    outputMultipoleMoments[11] = outputMultipoleMoments[9];

    RealOpenMM debye           = 4.80321;

    outputMultipoleMoments[0]  = netchg;

    dpl                       *= 10.0*debye;
    outputMultipoleMoments[1]  = dpl[0];
    outputMultipoleMoments[2]  = dpl[1];
    outputMultipoleMoments[3]  = dpl[2];

    debye *= 3.0;
    for (unsigned int ii = 4; ii < 13; ii++) {
        outputMultipoleMoments[ii] *= 100.0*debye;
    }
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateElectrostaticPotentialForParticleGridPoint(const MultipoleParticleData& particleI, const RealVec& gridPoint) const
{

    RealVec deltaR           = particleI.position - gridPoint;

    getPeriodicDelta(deltaR);

    RealOpenMM r2            = deltaR.dot(deltaR);
    RealOpenMM r             = SQRT(r2);

    RealOpenMM rr1           = 1.0/r;
    RealOpenMM rr2           = rr1*rr1;
    RealOpenMM rr3           = rr1*rr2;
    RealOpenMM potential     = particleI.charge*rr1;

    RealOpenMM scd           = particleI.dipole.dot(deltaR);
    RealOpenMM scu           = _inducedDipole[particleI.particleIndex].dot(deltaR);
    potential               -= (scd + scu)*rr3;

    RealOpenMM rr5           = 3.0*rr3*rr2;
    RealOpenMM scq           = deltaR[0]*(particleI.quadrupole[QXX]*deltaR[0] + particleI.quadrupole[QXY]*deltaR[1] + particleI.quadrupole[QXZ]*deltaR[2]);
          scq               += deltaR[1]*(particleI.quadrupole[QXY]*deltaR[0] + particleI.quadrupole[QYY]*deltaR[1] + particleI.quadrupole[QYZ]*deltaR[2]);
          scq               += deltaR[2]*(particleI.quadrupole[QXZ]*deltaR[0] + particleI.quadrupole[QYZ]*deltaR[1] + particleI.quadrupole[QZZ]*deltaR[2]);
    potential               += scq*rr5;
    return potential;
}

void AmoebaReferenceMultipoleForce::calculateElectrostaticPotential(const vector<RealVec>& particlePositions,
                                                                    const vector<RealOpenMM>& charges,
                                                                    const vector<RealOpenMM>& dipoles,
                                                                    const vector<RealOpenMM>& quadrupoles,
                                                                    const vector<RealOpenMM>& tholes,
                                                                    const vector<RealOpenMM>& dampingFactors,
                                                                    const vector<RealOpenMM>& polarity,
                                                                    const vector<int>& axisTypes,
                                                                    const vector<int>& multipoleAtomZs,
                                                                    const vector<int>& multipoleAtomXs,
                                                                    const vector<int>& multipoleAtomYs,
                                                                    const vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                                                    const vector<RealVec>& grid,
                                                                    vector<RealOpenMM>& potential)
{

    // setup, including calculating induced dipoles
    // initialize potential
    // calculate contribution of each particle to potential at grid point
    // apply prefactor

    vector<MultipoleParticleData> particleData;
    setup(particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData);

    potential.resize(grid.size());
    for (unsigned int ii = 0; ii < grid.size(); ii++) {
        potential[ii] = 0.0;
    }

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        for (unsigned int jj = 0; jj < grid.size(); jj++) {
            potential[jj] += calculateElectrostaticPotentialForParticleGridPoint(particleData[ii], grid[jj]);
        }
    }

    RealOpenMM term = _electric/_dielectric;
    for (unsigned int ii = 0; ii < grid.size(); ii++) {
        potential[ii] *= term;
    }
}

AmoebaReferenceMultipoleForce::UpdateInducedDipoleFieldStruct::UpdateInducedDipoleFieldStruct(vector<OpenMM::RealVec>& inputFixed_E_Field, vector<OpenMM::RealVec>& inputInducedDipoles, vector<vector<RealVec> >& extrapolatedDipoles, vector<vector<RealOpenMM> >& extrapolatedDipoleFieldGradient) :
        fixedMultipoleField(&inputFixed_E_Field), inducedDipoles(&inputInducedDipoles), extrapolatedDipoles(&extrapolatedDipoles), extrapolatedDipoleFieldGradient(&extrapolatedDipoleFieldGradient) { 
    inducedDipoleField.resize(fixedMultipoleField->size());
}

AmoebaReferenceGeneralizedKirkwoodMultipoleForce::AmoebaReferenceGeneralizedKirkwoodMultipoleForce(AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce) :
               AmoebaReferenceMultipoleForce(NoCutoff),
               _amoebaReferenceGeneralizedKirkwoodForce(amoebaReferenceGeneralizedKirkwoodForce),
               _gkc(2.455)
{

    RealOpenMM solventDielectric  = _amoebaReferenceGeneralizedKirkwoodForce->getSolventDielectric();

    _fc                           =      (1.0 - solventDielectric)/(solventDielectric);
    _fd                           =  2.0*(1.0 - solventDielectric)/(1.0 + 2.0*solventDielectric);
    _fq                           =  3.0*(1.0 - solventDielectric)/(2.0 + 3.0*solventDielectric);

    _amoebaReferenceGeneralizedKirkwoodForce->getGrycukBornRadii(_bornRadii);
    _amoebaReferenceGeneralizedKirkwoodForce->getAtomicRadii(_atomicRadii);
    _amoebaReferenceGeneralizedKirkwoodForce->getScaleFactors(_scaledRadii);

    _includeCavityTerm            = _amoebaReferenceGeneralizedKirkwoodForce->getIncludeCavityTerm();
    _probeRadius                  = _amoebaReferenceGeneralizedKirkwoodForce->getProbeRadius();
    _surfaceAreaFactor            = _amoebaReferenceGeneralizedKirkwoodForce->getSurfaceAreaFactor();
    _dielectricOffset             = _amoebaReferenceGeneralizedKirkwoodForce->getDielectricOffset();

    for (unsigned int ii = 0; ii < _scaledRadii.size(); ii++) {
        _scaledRadii[ii] *= _atomicRadii[ii];
    }
}

AmoebaReferenceGeneralizedKirkwoodMultipoleForce::~AmoebaReferenceGeneralizedKirkwoodMultipoleForce()
{
     delete _amoebaReferenceGeneralizedKirkwoodForce;
};

int AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getIncludeCavityTerm() const
{
     return _includeCavityTerm;
};

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getProbeRadius() const
{
     return _probeRadius;
};

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getSurfaceAreaFactor() const
{
     return _surfaceAreaFactor;
};

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getDielectricOffset() const
{
     return _dielectricOffset;
};

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::zeroFixedMultipoleFields()
{
    this->AmoebaReferenceMultipoleForce::zeroFixedMultipoleFields();
    initializeRealVecVector(_gkField);
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI,
                                                                                           const MultipoleParticleData& particleJ,
                                                                                           RealOpenMM dScale, RealOpenMM pScale)
{

    this->AmoebaReferenceMultipoleForce::calculateFixedMultipoleFieldPairIxn(particleI, particleJ, dScale, pScale);

    // get deltaR, R2, and R between 2 atoms

    RealVec deltaR               = particleJ.position - particleI.position;
    RealOpenMM r                 = SQRT(deltaR.dot(deltaR));
    RealOpenMM rb2               = _bornRadii[particleI.particleIndex]*_bornRadii[particleJ.particleIndex];

    RealOpenMM ci                = particleI.charge;

    RealOpenMM uxi               = particleI.dipole[0];
    RealOpenMM uyi               = particleI.dipole[1];
    RealOpenMM uzi               = particleI.dipole[2];

    RealOpenMM qxxi              = particleI.quadrupole[QXX];
    RealOpenMM qxyi              = particleI.quadrupole[QXY];
    RealOpenMM qxzi              = particleI.quadrupole[QXZ];
    RealOpenMM qyyi              = particleI.quadrupole[QYY];
    RealOpenMM qyzi              = particleI.quadrupole[QYZ];
    RealOpenMM qzzi              = particleI.quadrupole[QZZ];

    RealOpenMM xr                = deltaR[0];
    RealOpenMM yr                = deltaR[1];
    RealOpenMM zr                = deltaR[2];

    RealOpenMM ck                = particleJ.charge;

    RealOpenMM xr2               = xr*xr;
    RealOpenMM yr2               = yr*yr;
    RealOpenMM zr2               = zr*zr;
    RealOpenMM r2                = xr2 + yr2 + zr2;

    RealOpenMM uxk               = particleJ.dipole[0];
    RealOpenMM uyk               = particleJ.dipole[1];
    RealOpenMM uzk               = particleJ.dipole[2];

    RealOpenMM qxxk              = particleJ.quadrupole[QXX];
    RealOpenMM qxyk              = particleJ.quadrupole[QXY];
    RealOpenMM qxzk              = particleJ.quadrupole[QXZ];
    RealOpenMM qyyk              = particleJ.quadrupole[QYY];
    RealOpenMM qyzk              = particleJ.quadrupole[QYZ];
    RealOpenMM qzzk              = particleJ.quadrupole[QZZ];

    RealOpenMM expterm           = EXP(-r2/(_gkc*rb2));
    RealOpenMM expc              = expterm/_gkc;
    RealOpenMM dexpc             = -2.0/(_gkc*rb2);
    RealOpenMM gf2               = 1.0/(r2+rb2*expterm);
    RealOpenMM gf                = SQRT(gf2);
    RealOpenMM gf3               = gf2*gf;
    RealOpenMM gf5               = gf3*gf2;
    RealOpenMM gf7               = gf5*gf2;

    // reaction potential auxiliary terms

    RealOpenMM a[4][4];
    a[0][0]                      = gf;
    a[1][0]                      = -gf3;
    a[2][0]                      = 3.0*gf5;
    a[3][0]                      = -15.0*gf7;

    // reaction potential gradient auxiliary terms

    RealOpenMM expc1             = 1.0 - expc;
    a[0][1]                      = expc1*a[1][0];
    a[1][1]                      = expc1*a[2][0];
    a[2][1]                      = expc1*a[3][0];

    // dipole second reaction potential gradient auxiliary term

    RealOpenMM expcdexpc         = -expc*dexpc;
    a[1][2]                      = expc1*a[2][1] + expcdexpc*a[2][0];

    // multiply the auxillary terms by dielectric functions;

    a[0][1]                      = _fc*a[0][1];
    a[1][0]                      = _fd*a[1][0];
    a[1][1]                      = _fd*a[1][1];
    a[1][2]                      = _fd*a[1][2];
    a[2][0]                      = _fq*a[2][0];
    a[2][1]                      = _fq*a[2][1];

    // unweighted dipole reaction potential tensor

    RealOpenMM gux[11],guy[11],guz[11];
    gux[1]                       = xr*a[1][0];
    guy[1]                       = yr*a[1][0];
    guz[1]                       = zr*a[1][0];

    // unweighted reaction potential gradient tensor

    RealOpenMM gc[5];
    gc[2]                        = xr*a[0][1];
    gc[3]                        = yr*a[0][1];
    gc[4]                        = zr*a[0][1];

    gux[2]                       = a[1][0] + xr2*a[1][1];
    gux[3]                       = xr*yr*a[1][1];
    gux[4]                       = xr*zr*a[1][1];
    guy[2]                       = gux[3];
    guy[3]                       = a[1][0] + yr2*a[1][1];
    guy[4]                       = yr*zr*a[1][1];
    guz[2]                       = gux[4];
    guz[3]                       = guy[4];
    guz[4]                       = a[1][0] + zr2*a[1][1];

    RealOpenMM gqxx[5],gqxy[5];
    RealOpenMM gqxz[5],gqyy[5];
    RealOpenMM gqyz[5],gqzz[5];

    gqxx[2]                      = xr*(2.0*a[2][0]+xr2*a[2][1]);
    gqxx[3]                      = yr*xr2*a[2][1];
    gqxx[4]                      = zr*xr2*a[2][1];

    gqyy[2]                      = xr*yr2*a[2][1];
    gqyy[3]                      = yr*(2.0*a[2][0]+yr2*a[2][1]);
    gqyy[4]                      = zr*yr2*a[2][1];

    gqzz[2]                      = xr*zr2*a[2][1];
    gqzz[3]                      = yr*zr2*a[2][1];
    gqzz[4]                      = zr*(2.0*a[2][0]+zr2*a[2][1]);

    gqxy[2]                      = yr*(a[2][0]+xr2*a[2][1]);
    gqxy[3]                      = xr*(a[2][0]+yr2*a[2][1]);
    gqxy[4]                      = zr*xr*yr*a[2][1];

    gqxz[2]                      = zr*(a[2][0]+xr2*a[2][1]);
    gqxz[3]                      = gqxy[4];
    gqxz[4]                      = xr*(a[2][0]+zr2*a[2][1]);

    gqyz[2]                      = gqxy[4];
    gqyz[3]                      = zr*(a[2][0]+yr2*a[2][1]);
    gqyz[4]                      = yr*(a[2][0]+zr2*a[2][1]);

    // unweighted dipole second reaction potential gradient tensor

    gux[5]                       = xr*(3.0*a[1][1]+xr2*a[1][2]);
    gux[6]                       = yr*(a[1][1]+xr2*a[1][2]);
    gux[7]                       = zr*(a[1][1]+xr2*a[1][2]);

    gux[8]                       = xr*(a[1][1]+yr2*a[1][2]);
    gux[9]                       = zr*xr*yr*a[1][2];
    gux[10]                      = xr*(a[1][1]+zr2*a[1][2]);

    guy[5]                       = yr*(a[1][1]+xr2*a[1][2]);
    guy[6]                       = xr*(a[1][1]+yr2*a[1][2]);
    guy[7]                       = gux[9];

    guy[8]                       = yr*(3.0*a[1][1]+yr2*a[1][2]);
    guy[9]                       = zr*(a[1][1]+yr2*a[1][2]);
    guy[10]                      = yr*(a[1][1]+zr2*a[1][2]);
    guz[5]                       = zr*(a[1][1]+xr2*a[1][2]);
    guz[6]                       = gux[9];
    guz[7]                       = xr*(a[1][1]+zr2*a[1][2]);
    guz[8]                       = zr*(a[1][1]+yr2*a[1][2]);
    guz[9]                       = yr*(a[1][1]+zr2*a[1][2]);
    guz[10]                      = zr*(3.0*a[1][1]+zr2*a[1][2]);

    // generalized Kirkwood permanent reaction field

    RealVec fid;
    RealVec fjd;
    fid[0] = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                                   + 0.5*(ck*gux[1] + qxxk*gux[5]
                                   + qyyk*gux[8] + qzzk*gux[10]
                                   + 2.0*(qxyk*gux[6]+qxzk*gux[7]
                                   + qyzk*gux[9]))
                                   + 0.5*(ck*gc[2] + qxxk*gqxx[2]
                                   + qyyk*gqyy[2] + qzzk*gqzz[2]
                                   + 2.0*(qxyk*gqxy[2]+qxzk*gqxz[2]
                                   + qyzk*gqyz[2]));

    fid[1] = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                                   + 0.5*(ck*guy[1] + qxxk*guy[5]
                                   + qyyk*guy[8] + qzzk*guy[10]
                                   + 2.0*(qxyk*guy[6]+qxzk*guy[7]
                                   + qyzk*guy[9]))
                                   + 0.5*(ck*gc[3] + qxxk*gqxx[3]
                                   + qyyk*gqyy[3] + qzzk*gqzz[3]
                                   + 2.0*(qxyk*gqxy[3]+qxzk*gqxz[3]
                                   + qyzk*gqyz[3]));

    fid[2] = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                                   + 0.5*(ck*guz[1] + qxxk*guz[5]
                                   + qyyk*guz[8] + qzzk*guz[10]
                                   + 2.0*(qxyk*guz[6]+qxzk*guz[7]
                                   + qyzk*guz[9]))
                                   + 0.5*(ck*gc[4] + qxxk*gqxx[4]
                                   + qyyk*gqyy[4] + qzzk*gqzz[4]
                                   + 2.0*(qxyk*gqxy[4]+qxzk*gqxz[4]
                                   + qyzk*gqyz[4]));

    fjd[0] = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                                   - 0.5*(ci*gux[1] + qxxi*gux[5]
                                   + qyyi*gux[8] + qzzi*gux[10]
                                   + 2.0*(qxyi*gux[6]+qxzi*gux[7]
                                   + qyzi*gux[9]))
                                   - 0.5*(ci*gc[2] + qxxi*gqxx[2]
                                   + qyyi*gqyy[2] + qzzi*gqzz[2]
                                   + 2.0*(qxyi*gqxy[2]+qxzi*gqxz[2]
                                   + qyzi*gqyz[2]));

    fjd[1] = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                                   - 0.5*(ci*guy[1] + qxxi*guy[5]
                                   + qyyi*guy[8] + qzzi*guy[10]
                                   + 2.0*(qxyi*guy[6]+qxzi*guy[7]
                                   + qyzi*guy[9]))
                                   - 0.5*(ci*gc[3]      + qxxi*gqxx[3]
                                   + qyyi*gqyy[3] + qzzi*gqzz[3]
                                   + 2.0*(qxyi*gqxy[3]+qxzi*gqxz[3]
                                   + qyzi*gqyz[3]));

    fjd[2] = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                                   - 0.5*(ci*guz[1] + qxxi*guz[5]
                                   + qyyi*guz[8] + qzzi*guz[10]
                                   + 2.0*(qxyi*guz[6]+qxzi*guz[7]
                                   + qyzi*guz[9]))
                                   - 0.5*(ci*gc[4] + qxxi*gqxx[4]
                                   + qyyi*gqyy[4] + qzzi*gqzz[4]
                                   + 2.0*(qxyi*gqxy[4]+qxzi*gqxz[4]
                                   + qyzi*gqyz[4]));

    _gkField[particleI.particleIndex] += fid;
    if (particleI.particleIndex != particleJ.particleIndex) {
        _gkField[particleJ.particleIndex] += fjd;
    }
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipolePairGkIxn(const MultipoleParticleData& particleI,
                                                                                       const MultipoleParticleData& particleJ,
                                                                                       const vector<OpenMM::RealVec>& inputFields,
                                                                                       vector<OpenMM::RealVec>& outputFields) const
{

    RealOpenMM a[3][3];

    RealVec deltaR               = particleJ.position - particleI.position;
    RealOpenMM r                 = SQRT(deltaR.dot(deltaR));

    RealOpenMM xr                = deltaR[0];
    RealOpenMM yr                = deltaR[1];
    RealOpenMM zr                = deltaR[2];

    RealOpenMM xr2               = xr*xr;
    RealOpenMM yr2               = yr*yr;
    RealOpenMM zr2               = zr*zr;

    unsigned int iIndex          = particleI.particleIndex;
    unsigned int jIndex          = particleJ.particleIndex;
    RealOpenMM rb2               = _bornRadii[iIndex]*_bornRadii[jIndex];

    RealOpenMM r2                = xr2 + yr2 + zr2;
    RealOpenMM expterm           = EXP(-r2/(_gkc*rb2));
    RealOpenMM expc              = expterm /_gkc;

    RealOpenMM gf2               = 1.0/(r2+rb2*expterm);
    RealOpenMM gf                = SQRT(gf2);
    RealOpenMM gf3               = gf2*gf;
    RealOpenMM gf5               = gf3*gf2;

    RealVec duis                 = inputFields[iIndex];
    RealVec duks                 = inputFields[jIndex];

    // reaction potential auxiliary terms

    a[1][0]                      = -gf3;
    a[2][0]                      = 3.0*gf5;

    // reaction potential gradient auxiliary terms

    RealOpenMM expc1             = 1.0 - expc;
    a[1][1]                      = expc1*a[2][0];

    // unweighted dipole reaction potential gradient tensor

    RealVec gux, guy, guz;
    gux[0]                       = (a[1][0] + xr2*a[1][1]);
    gux[1]                       = xr*yr*a[1][1];
    gux[2]                       = xr*zr*a[1][1];

    guy[0]                       = gux[1];
    guy[1]                       = (a[1][0] + yr2*a[1][1]);
    guy[2]                       = yr*zr*a[1][1];

    guz[0]                       = gux[2];
    guz[1]                       = guy[2];
    guz[2]                       = (a[1][0] + zr2*a[1][1]);

    outputFields[iIndex][0]     += _fd*duks.dot(gux);
    outputFields[iIndex][1]     += _fd*duks.dot(guy);
    outputFields[iIndex][2]     += _fd*duks.dot(guz);

    // skip i == j, i.e., do not include contribution twice

    if (iIndex !=jIndex) {
        outputFields[jIndex][0] += _fd*duis.dot(gux);
        outputFields[jIndex][1] += _fd*duis.dot(guy);
        outputFields[jIndex][2] += _fd*duis.dot(guz);
    }
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                                                                      const MultipoleParticleData& particleJ,
                                                                                      vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxns(particleI, particleJ, updateInducedDipoleFields);

    // include GK contribution

    for (unsigned int ii = 2; ii < updateInducedDipoleFields.size(); ii++) {
        calculateInducedDipolePairGkIxn(particleI, particleJ, *updateInducedDipoleFields[ii].inducedDipoles, updateInducedDipoleFields[ii].inducedDipoleField);
    }
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipoles(const vector<MultipoleParticleData>& particleData)
{

    // calculate fixed electric fields

    zeroFixedMultipoleFields();
    calculateFixedMultipoleField(particleData);

    // initialize inducedDipoles

    _inducedDipole.resize(_numParticles);
    _inducedDipolePolar.resize(_numParticles);
    _inducedDipoleS.resize(_numParticles);
    _inducedDipolePolarS.resize(_numParticles);
    vector<RealVec> gkFieldPolar(_numParticles);

    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        _fixedMultipoleField[ii]      *= particleData[ii].polarity;
        _fixedMultipoleFieldPolar[ii] *= particleData[ii].polarity;
        _gkField[ii]                  *= particleData[ii].polarity;

        _inducedDipole[ii]             =  _fixedMultipoleField[ii];
        _inducedDipolePolar[ii]        =  _fixedMultipoleFieldPolar[ii];

        _inducedDipoleS[ii]            = (_fixedMultipoleField[ii]       + _gkField[ii]);
        _inducedDipolePolarS[ii]       = (_fixedMultipoleFieldPolar[ii]  + _gkField[ii]);

        _gkField[ii]                   = _inducedDipoleS[ii];
         gkFieldPolar[ii]              = _inducedDipolePolarS[ii];
    }

    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Direct) {
        setMutualInducedDipoleConverged(true);
        return;
    }

    vector<UpdateInducedDipoleFieldStruct> updateInducedDipoleField;
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(_fixedMultipoleField, _inducedDipole, _ptDipoleD, _ptDipoleFieldGradientD));
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(_fixedMultipoleFieldPolar, _inducedDipolePolar, _ptDipoleP, _ptDipoleFieldGradientP));
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(_gkField, _inducedDipoleS, _ptDipoleDS, _ptDipoleFieldGradientDS));
    updateInducedDipoleField.push_back(UpdateInducedDipoleFieldStruct(gkFieldPolar, _inducedDipolePolarS, _ptDipolePS, _ptDipoleFieldGradientPS));

    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Mutual)
        convergeInduceDipolesByDIIS(particleData, updateInducedDipoleField);
    else if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated)
        convergeInduceDipolesByExtrapolation(particleData, updateInducedDipoleField);
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateKirkwoodPairIxn(const MultipoleParticleData& particleI,
                                                                                      const MultipoleParticleData& particleJ,
                                                                                      vector<RealVec>& forces,
                                                                                      vector<RealVec>& torques,
                                                                                      vector<RealOpenMM>& dBorn) const
{

    RealOpenMM e,ei;
    RealOpenMM xr,yr,zr;
    RealOpenMM xr2,yr2,zr2;
    RealOpenMM sxi,syi,szi;
    RealOpenMM sxk,syk,szk;
    RealOpenMM r2,rb2;
    RealOpenMM dedx,dedy,dedz;
    RealOpenMM drbi;
    RealOpenMM drbk;
    RealOpenMM dpdx,dpdy,dpdz;
    RealOpenMM dpbi;
    RealOpenMM dpbk;
    RealOpenMM fc,fd,fq;
    RealOpenMM expterm;
    RealOpenMM gf,gf2,gf3,gf5;
    RealOpenMM gf7,gf9,gf11;
    RealOpenMM expc,dexpc;
    RealOpenMM expc1,expcdexpc;
    RealOpenMM expcr,dexpcr;
    RealOpenMM dgfdr;
    RealOpenMM esym,ewi,ewk;
    RealOpenMM desymdx,dewidx,dewkdx;
    RealOpenMM desymdy,dewidy,dewkdy;
    RealOpenMM desymdz,dewidz,dewkdz;
    RealOpenMM dsumdr,desymdr;
    RealOpenMM dewidr,dewkdr;
    RealOpenMM dsymdr;
    RealOpenMM esymi,ewii,ewki;
    RealOpenMM dpwidx,dpwkdx;
    RealOpenMM dpsymdy,dpwidy,dpwkdy;
    RealOpenMM dpsymdz,dpwidz,dpwkdz;
    RealOpenMM dwipdr,dwkpdr;
    RealOpenMM duvdr;

    unsigned int iIndex   = particleI.particleIndex;
    unsigned int jIndex   = particleJ.particleIndex;

    // set the bulk dielectric constant to the water value

    fc           = _electric*_fc;
    fd           = _electric*_fd;
    fq           = _electric*_fq;

    sxi          = _inducedDipoleS[iIndex][0] + _inducedDipolePolarS[iIndex][0];
    syi          = _inducedDipoleS[iIndex][1] + _inducedDipolePolarS[iIndex][1];
    szi          = _inducedDipoleS[iIndex][2] + _inducedDipolePolarS[iIndex][2];

    // decide whether to compute the current interaction

    RealVec deltaR = particleJ.position - particleI.position;
    RealOpenMM r = SQRT(deltaR.dot(deltaR));

    xr           = deltaR[0];
    yr           = deltaR[1];
    zr           = deltaR[2];

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    sxk          = _inducedDipoleS[jIndex][0] + _inducedDipolePolarS[jIndex][0];
    syk          = _inducedDipoleS[jIndex][1] + _inducedDipolePolarS[jIndex][1];
    szk          = _inducedDipoleS[jIndex][2] + _inducedDipolePolarS[jIndex][2];

    rb2          = _bornRadii[iIndex]*_bornRadii[jIndex];

    expterm      = EXP(-r2/(_gkc*rb2));
    expc         = expterm/_gkc;
    expcr        = r2*expterm/(_gkc*_gkc*rb2*rb2);
    dexpc        = -2.0 / (_gkc*rb2);
    dexpcr       = 2.0 / (_gkc*rb2*rb2);
    dgfdr        = 0.5 * expterm * (1.0+r2/(rb2*_gkc));
    gf2          = 1.0 / (r2+rb2*expterm);
    gf           = SQRT(gf2);
    gf3          = gf2*gf;
    gf5          = gf3*gf2;
    gf7          = gf5*gf2;
    gf9          = gf7*gf2;
    gf11         = gf9*gf2;

    // reaction potential auxiliary terms;

    RealOpenMM a00      =        gf;
    RealOpenMM a10      =       -gf3;
    RealOpenMM a20      =  3.0*gf5;
    RealOpenMM a30      =  -15.0*gf7;
    RealOpenMM a40      =  105.0*gf9;
    RealOpenMM a50      = -945.0*gf11;

    // Born radii derivatives of reaction potential auxiliary terms;

    RealOpenMM b00      = dgfdr*a10;
    RealOpenMM b10      = dgfdr*a20;
    RealOpenMM b20      = dgfdr*a30;
    RealOpenMM b30      = dgfdr*a40;
    RealOpenMM b40      = dgfdr*a50;

    // reaction potential gradient auxiliary terms;

    expc1               = 1.0 - expc;
    RealOpenMM a01      = expc1*a10;
    RealOpenMM a11      = expc1*a20;
    RealOpenMM a21      = expc1*a30;
    RealOpenMM a31      = expc1*a40;
    RealOpenMM a41      = expc1*a50;

    // Born radii derivs of reaction potential gradient auxiliary terms;

    RealOpenMM b01      = b10 - expcr*a10 - expc*b10;
    RealOpenMM b11      = b20 - expcr*a20 - expc*b20;
    RealOpenMM b21      = b30 - expcr*a30 - expc*b30;
    RealOpenMM b31      = b40 - expcr*a40 - expc*b40;

    // 2nd reaction potential gradient auxiliary terms;

    expcdexpc           = -expc*dexpc;
    RealOpenMM a02      = expc1*a11 + expcdexpc*a10;
    RealOpenMM a12      = expc1*a21 + expcdexpc*a20;
    RealOpenMM a22      = expc1*a31 + expcdexpc*a30;
    RealOpenMM a32      = expc1*a41 + expcdexpc*a40;

    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     RealOpenMM b02     = b11 - (expcr*(a11 + dexpc*a10)
                           + expc*(b11 + dexpcr*a10+dexpc*b10));

     RealOpenMM b12     = b21 - (expcr*(a21 + dexpc*a20)
                           +   expc*(b21 + dexpcr*a20+dexpc*b20));

     RealOpenMM b22     = b31 - (expcr*(a31 + dexpc*a30)
                           +   expc*(b31 + dexpcr*a30+dexpc*b30));

    // 3rd reaction potential gradient auxiliary terms

    expcdexpc           = 2.0*expcdexpc;
    RealOpenMM a03      = expc1*a12 + expcdexpc*a11;
    RealOpenMM a13      = expc1*a22 + expcdexpc*a21;
    RealOpenMM a23      = expc1*a32 + expcdexpc*a31;

    expcdexpc      = -expc*dexpc*dexpc;
    a03            = a03 + expcdexpc*a10;
    a13            = a13 + expcdexpc*a20;
    a23            = a23 + expcdexpc*a30;

    // multiply the auxillary terms by their dieletric functions;

    a00            *= fc;
    a01            *= fc;
    a02            *= fc;
    a03            *= fc;

    b00            *= fc;
    b01            *= fc;
    b02            *= fc;

    a10            *= fd;
    a11            *= fd;
    a12            *= fd;
    a13            *= fd;

    b10            *= fd;
    b11            *= fd;
    b12            *= fd;

    a20            *= fq;
    a21            *= fq;
    a22            *= fq;
    a23            *= fq;

    b20            *= fq;
    b21            *= fq;
    b22            *= fq;

    // unweighted reaction potential tensor

    RealOpenMM gc1        = a00;

    RealOpenMM gux1       = xr*a10;
    RealOpenMM guy1       = yr*a10;
    RealOpenMM guz1       = zr*a10;

    RealOpenMM gqxx1      = xr2*a20;
    RealOpenMM gqyy1      = yr2*a20;
    RealOpenMM gqzz1      = zr2*a20;

    RealOpenMM gqxy1      = xr*yr*a20;
    RealOpenMM gqxz1      = xr*zr*a20;
    RealOpenMM gqyz1      = yr*zr*a20;

    // Born radii derivs of unweighted reaction potential tensor

    RealOpenMM gc21       = b00;

    RealOpenMM gux21      = xr*b10;
    RealOpenMM guy21      = yr*b10;
    RealOpenMM guz21      = zr*b10;

    RealOpenMM gqxx21     = xr2*b20;
    RealOpenMM gqyy21     = yr2*b20;
    RealOpenMM gqzz21     = zr2*b20;

    RealOpenMM gqxy21     = xr*yr*b20;
    RealOpenMM gqxz21     = xr*zr*b20;
    RealOpenMM gqyz21     = yr*zr*b20;

    // unweighted reaction potential gradient tensor;

    RealOpenMM gc2        = xr*a01;
    RealOpenMM gc3        = yr*a01;
    RealOpenMM gc4        = zr*a01;

    RealOpenMM gux2       = a10 + xr2*a11;
    RealOpenMM gux3       = xr*yr*a11;
    RealOpenMM gux4       = xr*zr*a11;

    RealOpenMM guy2       = gux3;
    RealOpenMM guy3       = a10 + yr2*a11;
    RealOpenMM guy4       = yr*zr*a11;
    RealOpenMM guz2       = gux4;
    RealOpenMM guz3       = guy4;
    RealOpenMM guz4       = a10 + zr2*a11;
    RealOpenMM gqxx2      = xr*(2.0*a20+xr2*a21);
    RealOpenMM gqxx3      = yr*xr2*a21;
    RealOpenMM gqxx4      = zr*xr2*a21;
    RealOpenMM gqyy2      = xr*yr2*a21;
    RealOpenMM gqyy3      = yr*(2.0*a20+yr2*a21);
    RealOpenMM gqyy4      = zr*yr2*a21;
    RealOpenMM gqzz2      = xr*zr2*a21;
    RealOpenMM gqzz3      = yr*zr2*a21;
    RealOpenMM gqzz4      = zr*(2.0*a20+zr2*a21);
    RealOpenMM gqxy2      = yr*(a20+xr2*a21);
    RealOpenMM gqxy3      = xr*(a20+yr2*a21);
    RealOpenMM gqxy4      = zr*xr*yr*a21;
    RealOpenMM gqxz2      = zr*(a20+xr2*a21);
    RealOpenMM gqxz3      = gqxy4;
    RealOpenMM gqxz4      = xr*(a20+zr2*a21);
    RealOpenMM gqyz2      = gqxy4;
    RealOpenMM gqyz3      = zr*(a20+yr2*a21);
    RealOpenMM gqyz4      = yr*(a20+zr2*a21);

    // Born derivs of the unweighted reaction potential gradient tensor

    RealOpenMM gc22       = xr*b01;
    RealOpenMM gc23       = yr*b01;
    RealOpenMM gc24       = zr*b01;
    RealOpenMM gux22      = b10 + xr2*b11;
    RealOpenMM gux23      = xr * yr * b11;
    RealOpenMM gux24      = xr*zr*b11;
    RealOpenMM guy22      = gux23;
    RealOpenMM guy23      = b10 + yr2*b11;
    RealOpenMM guy24      = yr*zr*b11;
    RealOpenMM guz22      = gux24;
    RealOpenMM guz23      = guy24;
    RealOpenMM guz24      = b10 + zr2*b11;
    RealOpenMM gqxx22     = xr*(2.0*b20+xr2*b21);
    RealOpenMM gqxx23     = yr*xr2*b21;
    RealOpenMM gqxx24     = zr*xr2*b21;
    RealOpenMM gqyy22     = xr*yr2*b21;
    RealOpenMM gqyy23     = yr*(2.0*b20+yr2*b21);
    RealOpenMM gqyy24     = zr*yr2*b21;
    RealOpenMM gqzz22     = xr*zr2*b21;
    RealOpenMM gqzz23     = yr*zr2*b21;
    RealOpenMM gqzz24     = zr*(2.0*b20 + zr2*b21);
    RealOpenMM gqxy22     = yr*(b20+xr2*b21);
    RealOpenMM gqxy23     = xr*(b20+yr2*b21);
    RealOpenMM gqxy24     = zr*xr*yr*b21;
    RealOpenMM gqxz22     = zr*(b20+xr2*b21);
    RealOpenMM gqxz23     = gqxy24;
    RealOpenMM gqxz24     = xr*(b20+zr2*b21);
    RealOpenMM gqyz22     = gqxy24;
    RealOpenMM gqyz23     = zr*(b20+yr2*b21);
    RealOpenMM gqyz24     = yr*(b20+zr2*b21);

    // unweighted 2nd reaction potential gradient tensor;

    RealOpenMM gc5        = a01 + xr2*a02;
    RealOpenMM gc6        = xr*yr*a02;
    RealOpenMM gc7        = xr*zr*a02;
    RealOpenMM gc8        = a01 + yr2*a02;
    RealOpenMM gc9        = yr*zr*a02;
    RealOpenMM gc10       = a01 + zr2*a02;
    RealOpenMM gux5       = xr*(3.0*a11+xr2*a12);
    RealOpenMM gux6       = yr*(a11+xr2*a12);
    RealOpenMM gux7       = zr*(a11+xr2*a12);
    RealOpenMM gux8       = xr*(a11+yr2*a12);
    RealOpenMM gux9       = zr*xr*yr*a12;
    RealOpenMM gux10      = xr*(a11+zr2*a12);
    RealOpenMM guy5       = yr*(a11+xr2*a12);
    RealOpenMM guy6       = xr*(a11+yr2*a12);
    RealOpenMM guy7       = gux9;
    RealOpenMM guy8       = yr*(3.0*a11+yr2*a12);
    RealOpenMM guy9       = zr*(a11+yr2*a12);
    RealOpenMM guy10      = yr*(a11+zr2*a12);
    RealOpenMM guz5       = zr*(a11+xr2*a12);
    RealOpenMM guz6       = gux9;
    RealOpenMM guz7       = xr*(a11+zr2*a12);
    RealOpenMM guz8       = zr*(a11+yr2*a12);
    RealOpenMM guz9       = yr*(a11+zr2*a12);
    RealOpenMM guz10      = zr*(3.0*a11+zr2*a12);
    RealOpenMM gqxx5      = 2.0*a20 + xr2*(5.0*a21+xr2*a22);
    RealOpenMM gqxx6      = yr*xr*(2.0*a21+xr2*a22);
    RealOpenMM gqxx7      = zr*xr*(2.0*a21+xr2*a22);
    RealOpenMM gqxx8      = xr2*(a21+yr2*a22);
    RealOpenMM gqxx9      = zr*yr*xr2*a22;
    RealOpenMM gqxx10     = xr2*(a21+zr2*a22);
    RealOpenMM gqyy5      = yr2*(a21+xr2*a22);
    RealOpenMM gqyy6      = xr*yr*(2.0*a21+yr2*a22);
    RealOpenMM gqyy7      = xr*zr*yr2*a22;
    RealOpenMM gqyy8      = 2.0*a20 + yr2*(5.0*a21+yr2*a22);
    RealOpenMM gqyy9      = yr*zr*(2.0*a21+yr2*a22);
    RealOpenMM gqyy10     = yr2*(a21+zr2*a22);
    RealOpenMM gqzz5      = zr2*(a21+xr2*a22);
    RealOpenMM gqzz6      = xr*yr*zr2*a22;
    RealOpenMM gqzz7      = xr*zr*(2.0*a21+zr2*a22);
    RealOpenMM gqzz8      = zr2*(a21+yr2*a22);
    RealOpenMM gqzz9      = yr*zr*(2.0*a21+zr2*a22);
    RealOpenMM gqzz10     = 2.0*a20 + zr2*(5.0*a21+zr2*a22);
    RealOpenMM gqxy5      = xr*yr*(3.0*a21+xr2*a22);
    RealOpenMM gqxy6      = a20 + (xr2+yr2)*a21 + xr2*yr2*a22;
    RealOpenMM gqxy7      = zr*yr*(a21+xr2*a22);
    RealOpenMM gqxy8      = xr*yr*(3.0*a21+yr2*a22);
    RealOpenMM gqxy9      = zr*xr*(a21+yr2*a22);
    RealOpenMM gqxy10     = xr*yr*(a21+zr2*a22);
    RealOpenMM gqxz5      = xr*zr*(3.0*a21+xr2*a22);
    RealOpenMM gqxz6      = yr*zr*(a21+xr2*a22);
    RealOpenMM gqxz7      = a20 + (xr2+zr2)*a21 + xr2*zr2*a22;
    RealOpenMM gqxz8      = xr*zr*(a21+yr2*a22);
    RealOpenMM gqxz9      = xr*yr*(a21+zr2*a22);
    RealOpenMM gqxz10     = xr*zr*(3.0*a21+zr2*a22);
    RealOpenMM gqyz5      = zr*yr*(a21+xr2*a22);
    RealOpenMM gqyz6      = xr*zr*(a21+yr2*a22);
    RealOpenMM gqyz7      = xr*yr*(a21+zr2*a22);
    RealOpenMM gqyz8      = yr*zr*(3.0*a21+yr2*a22);
    RealOpenMM gqyz9      = a20 + (yr2+zr2)*a21 + yr2*zr2*a22;
    RealOpenMM gqyz10     = yr*zr*(3.0*a21+zr2*a22);

    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

    RealOpenMM gc25       = b01 + xr2*b02;
    RealOpenMM gc26       = xr*yr*b02;
    RealOpenMM gc27       = xr*zr*b02;
    RealOpenMM gc28       = b01 + yr2*b02;
    RealOpenMM gc29       = yr*zr*b02;
    RealOpenMM gc30       = b01 + zr2*b02;
    RealOpenMM gux25      = xr*(3.0*b11+xr2*b12);
    RealOpenMM gux26      = yr*(b11+xr2*b12);
    RealOpenMM gux27      = zr*(b11+xr2*b12);
    RealOpenMM gux28      = xr*(b11+yr2*b12);
    RealOpenMM gux29      = zr*xr*yr*b12;
    RealOpenMM gux30      = xr*(b11+zr2*b12);
    RealOpenMM guy25      = yr*(b11+xr2*b12);
    RealOpenMM guy26      = xr*(b11+yr2*b12);
    RealOpenMM guy27      = gux29;
    RealOpenMM guy28      = yr*(3.0*b11+yr2*b12);
    RealOpenMM guy29      = zr*(b11+yr2*b12);
    RealOpenMM guy30      = yr*(b11+zr2*b12);
    RealOpenMM guz25      = zr*(b11+xr2*b12);
    RealOpenMM guz26      = gux29;
    RealOpenMM guz27      = xr*(b11+zr2*b12);
    RealOpenMM guz28      = zr*(b11+yr2*b12);
    RealOpenMM guz29      = yr*(b11+zr2*b12);
    RealOpenMM guz30      = zr*(3.0*b11+zr2*b12);
    RealOpenMM gqxx25     = 2.0*b20 + xr2*(5.0*b21+xr2*b22);
    RealOpenMM gqxx26     = yr*xr*(2.0*b21+xr2*b22);
    RealOpenMM gqxx27     = zr*xr*(2.0*b21+xr2*b22);
    RealOpenMM gqxx28     = xr2*(b21+yr2*b22);
    RealOpenMM gqxx29     = zr*yr*xr2*b22;
    RealOpenMM gqxx30     = xr2*(b21+zr2*b22);
    RealOpenMM gqyy25     = yr2*(b21+xr2*b22);
    RealOpenMM gqyy26     = xr*yr*(2.0*b21+yr2*b22);
    RealOpenMM gqyy27     = xr*zr*yr2*b22;
    RealOpenMM gqyy28     = 2.0*b20 + yr2*(5.0*b21+yr2*b22);
    RealOpenMM gqyy29     = yr*zr*(2.0*b21+yr2*b22);
    RealOpenMM gqyy30     = yr2*(b21+zr2*b22);
    RealOpenMM gqzz25     = zr2*(b21+xr2*b22);
    RealOpenMM gqzz26     = xr*yr*zr2*b22;
    RealOpenMM gqzz27     = xr*zr*(2.0*b21+zr2*b22);
    RealOpenMM gqzz28     = zr2*(b21+yr2*b22);
    RealOpenMM gqzz29     = yr*zr*(2.0*b21+zr2*b22);
    RealOpenMM gqzz30     = 2.0*b20 + zr2*(5.0*b21+zr2*b22);
    RealOpenMM gqxy25     = xr*yr*(3.0*b21 + xr2*b22);
    RealOpenMM gqxy26     = b20 + (xr2+yr2)*b21 + xr2*yr2*b22;
    RealOpenMM gqxy27     = zr*yr*(b21+xr2*b22);
    RealOpenMM gqxy28     = xr*yr*(3.0*b21+yr2*b22);
    RealOpenMM gqxy29     = zr*xr*(b21+yr2*b22);
    RealOpenMM gqxy30     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqxz25     = xr*zr*(3.0*b21+xr2*b22);
    RealOpenMM gqxz26     = yr*zr*(b21+xr2*b22);
    RealOpenMM gqxz27     = b20 + (xr2+zr2)*b21 + xr2*zr2*b22;
    RealOpenMM gqxz28     = xr*zr*(b21+yr2*b22);
    RealOpenMM gqxz29     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqxz30     = xr*zr*(3.0*b21+zr2*b22);
    RealOpenMM gqyz25     = zr*yr*(b21+xr2*b22);
    RealOpenMM gqyz26     = xr*zr*(b21+yr2*b22);
    RealOpenMM gqyz27     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqyz28     = yr*zr*(3.0*b21+yr2*b22);
    RealOpenMM gqyz29     = b20 + (yr2+zr2)*b21 + yr2*zr2*b22;
    RealOpenMM gqyz30     = yr*zr*(3.0*b21+zr2*b22);

    // unweighted 3rd reaction potential gradient tensor;

    RealOpenMM gc11       = xr*(3.0*a02+xr2*a03);
    RealOpenMM gc12       = yr*(a02+xr2*a03);
    RealOpenMM gc13       = zr*(a02+xr2*a03);
    RealOpenMM gc14       = xr*(a02+yr2*a03);
    RealOpenMM gc15       = xr*yr*zr*a03;
    RealOpenMM gc16       = xr*(a02+zr2*a03);
    RealOpenMM gc17       = yr*(3.0*a02+yr2*a03);
    RealOpenMM gc18       = zr*(a02+yr2*a03);
    RealOpenMM gc19       = yr*(a02+zr2*a03);
    RealOpenMM gc20       = zr*(3.0*a02+zr2*a03);
    RealOpenMM gux11      = 3.0*a11 + xr2*(6.0*a12+xr2*a13);
    RealOpenMM gux12      = xr*yr*(3.0*a12+xr2*a13);
    RealOpenMM gux13      = xr*zr*(3.0*a12+xr2*a13);
    RealOpenMM gux14      = a11 + (xr2+yr2)*a12 + xr2*yr2*a13;
    RealOpenMM gux15      = yr*zr*(a12+xr2*a13);
    RealOpenMM gux16      = a11 + (xr2+zr2)*a12 + xr2*zr2*a13;
    RealOpenMM gux17      = xr*yr*(3.0*a12+yr2*a13);
    RealOpenMM gux18      = xr*zr*(a12+yr2*a13);
    RealOpenMM gux19      = xr*yr*(a12+zr2*a13);
    RealOpenMM gux20      = xr*zr*(3.0*a12+zr2*a13);
    RealOpenMM guy11      = gux12;
    RealOpenMM guy12      = gux14;
    RealOpenMM guy13      = gux15;
    RealOpenMM guy14      = gux17;
    RealOpenMM guy15      = gux18;
    RealOpenMM guy16      = gux19;
    RealOpenMM guy17      = 3.0*a11 + yr2*(6.0*a12+yr2*a13);
    RealOpenMM guy18      = yr*zr*(3.0*a12+yr2*a13);
    RealOpenMM guy19      = a11 + (yr2+zr2)*a12 + yr2*zr2*a13;
    RealOpenMM guy20      = yr*zr*(3.0*a12+zr2*a13);
    RealOpenMM guz11      = gux13;
    RealOpenMM guz12      = gux15;
    RealOpenMM guz13      = gux16;
    RealOpenMM guz14      = gux18;
    RealOpenMM guz15      = gux19;
    RealOpenMM guz16      = gux20;
    RealOpenMM guz17      = guy18;
    RealOpenMM guz18      = guy19;
    RealOpenMM guz19      = guy20;
    RealOpenMM guz20      = 3.0*a11 + zr2*(6.0*a12+zr2*a13);

    RealOpenMM gqxx11     = xr*(12.0*a21+xr2*(9.0*a22 + xr2*a23));
    RealOpenMM gqxx12     = yr*(2.0*a21+xr2*(5.0*a22  + xr2*a23));
    RealOpenMM gqxx13     = zr*(2.0*a21+xr2*(5.0*a22  + xr2*a23));
    RealOpenMM gqxx14     = xr*(2.0*a21+yr2*2.0*a22   +xr2*(a22+yr2*a23));
    RealOpenMM gqxx15     = xr*yr*zr*(2.0*a22+xr2*a23);
    RealOpenMM gqxx16     = xr*(2.0*a21+zr2*2.0*a22 +xr2*(a22+zr2*a23));
    RealOpenMM gqxx17     = yr*xr2*(3.0*a22+yr2*a23);
    RealOpenMM gqxx18     = zr*xr2*(a22+yr2*a23);
    RealOpenMM gqxx19     = yr*xr2*(a22+zr2*a23);
    RealOpenMM gqxx20     = zr*xr2*(3.0*a22+zr2*a23);

    RealOpenMM gqxy11     = yr*(3.0*a21+xr2*(6.0*a22 +xr2*a23));
    RealOpenMM gqxy12     = xr*(3.0*(a21+yr2*a22) +xr2*(a22+yr2*a23));
    RealOpenMM gqxy13     = xr*yr*zr*(3.0*a22+xr2*a23);
    RealOpenMM gqxy14     = yr*(3.0*(a21+xr2*a22) +yr2*(a22+xr2*a23));
    RealOpenMM gqxy15     = zr*(a21+(yr2+xr2)*a22 +yr2*xr2*a23);
    RealOpenMM gqxy16     = yr*(a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    RealOpenMM gqxy17     = xr*(3.0*(a21+yr2*a22) +yr2*(3.0*a22+yr2*a23));
    RealOpenMM gqxy18     = xr*yr*zr*(3.0*a22+yr2*a23);
    RealOpenMM gqxy19     = xr*(a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    RealOpenMM gqxy20     = xr*yr*zr*(3.0*a22+zr2*a23);
    RealOpenMM gqxz11     = zr*(3.0*a21+xr2*(6.0*a22 +xr2*a23));
    RealOpenMM gqxz12     = xr*yr*zr*(3.0*a22+xr2*a23);
    RealOpenMM gqxz13     = xr*(3.0*(a21+zr2*a22) +xr2*(a22+zr2*a23));
    RealOpenMM gqxz14     = zr*(a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    RealOpenMM gqxz15     = yr*(a21+(xr2+zr2)*a22 +zr2*xr2*a23);
    RealOpenMM gqxz16     = zr*(3.0*(a21+xr2*a22) +zr2*(a22+xr2*a23));
    RealOpenMM gqxz17     = xr*yr*zr*(3.0*a22+yr2*a23);
    RealOpenMM gqxz18     = xr*(a21+(zr2+yr2)*a22 +zr2*yr2*a23);
    RealOpenMM gqxz19     = xr*yr*zr*(3.0*a22+zr2*a23);
    RealOpenMM gqxz20     = xr*(3.0*a21+zr2*(6.0*a22 +zr2*a23));
    RealOpenMM gqyy11     = xr*yr2*(3.0*a22+xr2*a23);
    RealOpenMM gqyy12     = yr*(2.0*a21+xr2*2.0*a22 +yr2*(a22+xr2*a23));
    RealOpenMM gqyy13     = zr*yr2*(a22+xr2*a23);
    RealOpenMM gqyy14     = xr*(2.0*a21+yr2*(5.0*a22 +yr2*a23));
    RealOpenMM gqyy15     = xr*yr*zr*(2.0*a22+yr2*a23);
    RealOpenMM gqyy16     = xr*yr2*(a22+zr2*a23);
    RealOpenMM gqyy17     = yr*(12.0*a21+yr2*(9.0*a22 +yr2*a23));
    RealOpenMM gqyy18     = zr*(2.0*a21+yr2*(5.0*a22 +yr2*a23));
    RealOpenMM gqyy19     = yr*(2.0*a21+zr2*2.0*a22 +yr2*(a22+zr2*a23));
    RealOpenMM gqyy20     = zr*yr2*(3.0*a22+zr2*a23);
    RealOpenMM gqyz11     = xr*yr*zr*(3.0*a22+xr2*a23);
    RealOpenMM gqyz12     = zr*(a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    RealOpenMM gqyz13     = yr*(a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    RealOpenMM gqyz14     = xr*yr*zr*(3.0*a22+yr2*a23);
    RealOpenMM gqyz15     = xr*(a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    RealOpenMM gqyz16     = xr*yr*zr*(3.0*a22+zr2*a23);
    RealOpenMM gqyz17     = zr*(3.0*a21+yr2*(6.0*a22 +yr2*a23));
    RealOpenMM gqyz18     = yr*(3.0*(a21+zr2*a22) +yr2*(a22+zr2*a23));
    RealOpenMM gqyz19     = zr*(3.0*(a21+yr2*a22) +zr2*(a22+yr2*a23));
    RealOpenMM gqyz20     = yr*(3.0*a21+zr2*(6.0*a22 +zr2*a23));
    RealOpenMM gqzz11     = xr*zr2*(3.0*a22+xr2*a23);
    RealOpenMM gqzz12     = yr*(zr2*a22+xr2*(zr2*a23));
    RealOpenMM gqzz13     = zr*(2.0*a21+xr2*2.0*a22 +zr2*(a22+xr2*a23));
    RealOpenMM gqzz14     = xr*zr2*(a22+yr2*a23);
    RealOpenMM gqzz15     = xr*yr*zr*(2.0*a22+zr2*a23);
    RealOpenMM gqzz16     = xr*(2.0*a21+zr2*(5.0*a22 +zr2*a23));
    RealOpenMM gqzz17     = yr*zr2*(3.0*a22+yr2*a23);
    RealOpenMM gqzz18     = zr*(2.0*a21+yr2*2.0*a22 +zr2*(a22+yr2*a23));
    RealOpenMM gqzz19     = yr*(2.0*a21+zr2*(5.0*a22 +zr2*a23));
    RealOpenMM gqzz20     = zr*(12.0*a21+zr2*(9.0*a22 +zr2*a23));

    // electrostatic solvation energy of the permanent multipoles
    // in their own GK reaction potential

    esym = particleI.charge*particleJ.charge*gc1 - (particleI.dipole[0]*(particleJ.dipole[0]*gux2+particleJ.dipole[1]*guy2+particleJ.dipole[2]*guz2)
                                   +  particleI.dipole[1]*(particleJ.dipole[0]*gux3+particleJ.dipole[1]*guy3+particleJ.dipole[2]*guz3)
                                   +  particleI.dipole[2]*(particleJ.dipole[0]*gux4+particleJ.dipole[1]*guy4+particleJ.dipole[2]*guz4));

    ewi =  particleI.charge*(particleJ.dipole[0]*gc2+particleJ.dipole[1]*gc3+particleJ.dipole[2]*gc4)
          -particleJ.charge*(particleI.dipole[0]*gux1+particleI.dipole[1]*guy1+particleI.dipole[2]*guz1)
           +particleI.charge*(particleJ.quadrupole[QXX]*gc5+particleJ.quadrupole[QYY]*gc8+particleJ.quadrupole[QZZ]*gc10
              +2.0*(particleJ.quadrupole[QXY]*gc6+particleJ.quadrupole[QXZ]*gc7+particleJ.quadrupole[QYZ]*gc9))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx1+particleI.quadrupole[QYY]*gqyy1+particleI.quadrupole[QZZ]*gqzz1
              +2.0*(particleI.quadrupole[QXY]*gqxy1+particleI.quadrupole[QXZ]*gqxz1+particleI.quadrupole[QYZ]*gqyz1))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
              +2.0*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
              +2.0*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
              +2.0*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
              +2.0*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
              +2.0*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
              +2.0*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqxx8+particleJ.quadrupole[QZZ]*gqxx10
              +2.0*(particleJ.quadrupole[QXY]*gqxx6+particleJ.quadrupole[QXZ]*gqxx7+particleJ.quadrupole[QYZ]*gqxx9))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy5+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqyy10
              +2.0*(particleJ.quadrupole[QXY]*gqyy6+particleJ.quadrupole[QXZ]*gqyy7+particleJ.quadrupole[QYZ]*gqyy9))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz5+particleJ.quadrupole[QYY]*gqzz8+particleJ.quadrupole[QZZ]*gqzz10
              +2.0*(particleJ.quadrupole[QXY]*gqzz6+particleJ.quadrupole[QXZ]*gqzz7+particleJ.quadrupole[QYZ]*gqzz9))
              + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy5+particleJ.quadrupole[QYY]*gqxy8+particleJ.quadrupole[QZZ]*gqxy10
              +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxy7+particleJ.quadrupole[QYZ]*gqxy9))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz5+particleJ.quadrupole[QYY]*gqxz8+particleJ.quadrupole[QZZ]*gqxz10
              +2.0*(particleJ.quadrupole[QXY]*gqxz6+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqxz9))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz5+particleJ.quadrupole[QYY]*gqyz8+particleJ.quadrupole[QZZ]*gqyz10
              +2.0*(particleJ.quadrupole[QXY]*gqyz6+particleJ.quadrupole[QXZ]*gqyz7+particleJ.quadrupole[QYZ]*gqyz9)));

    ewk = particleI.charge*(particleJ.dipole[0]*gux1+particleJ.dipole[1]*guy1+particleJ.dipole[2]*guz1)
                      -particleJ.charge*(particleI.dipole[0]*gc2+particleI.dipole[1]*gc3+particleI.dipole[2]*gc4)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx1+particleJ.quadrupole[QYY]*gqyy1+particleJ.quadrupole[QZZ]*gqzz1
              +2.0*(particleJ.quadrupole[QXY]*gqxy1+particleJ.quadrupole[QXZ]*gqxz1+particleJ.quadrupole[QYZ]*gqyz1))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc5+particleI.quadrupole[QYY]*gc8+particleI.quadrupole[QZZ]*gc10
              +2.0*(particleI.quadrupole[QXY]*gc6+particleI.quadrupole[QXZ]*gc7+particleI.quadrupole[QYZ]*gc9))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
              +2.0*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
              +2.0*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
              +2.0*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
              +2.0*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
              +2.0*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
              +2.0*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
              +2.0*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
              +2.0*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
              +2.0*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9)));

    desymdx = particleI.charge*particleJ.charge*gc2 - (particleI.dipole[0]*(particleJ.dipole[0]*gux5+particleJ.dipole[1]*guy5+particleJ.dipole[2]*guz5)
                              +  particleI.dipole[1]*(particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6)
                              +  particleI.dipole[2]*(particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7));

    dewidx = particleI.charge*(particleJ.dipole[0]*gc5+particleJ.dipole[1]*gc6+particleJ.dipole[2]*gc7)
                      -particleJ.charge*(particleI.dipole[0]*gux2+particleI.dipole[1]*guy2+particleI.dipole[2]*guz2)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc11+particleJ.quadrupole[QYY]*gc14+particleJ.quadrupole[QZZ]*gc16
              +2.0*(particleJ.quadrupole[QXY]*gc12+particleJ.quadrupole[QXZ]*gc13+particleJ.quadrupole[QYZ]*gc15))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
              +2.0*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux11+particleJ.quadrupole[QYY]*gux14+particleJ.quadrupole[QZZ]*gux16
              +2.0*(particleJ.quadrupole[QXY]*gux12+particleJ.quadrupole[QXZ]*gux13+particleJ.quadrupole[QYZ]*gux15))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy11+particleJ.quadrupole[QYY]*guy14+particleJ.quadrupole[QZZ]*guy16
              +2.0*(particleJ.quadrupole[QXY]*guy12+particleJ.quadrupole[QXZ]*guy13+particleJ.quadrupole[QYZ]*guy15))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz11+particleJ.quadrupole[QYY]*guz14+particleJ.quadrupole[QZZ]*guz16
              +2.0*(particleJ.quadrupole[QXY]*guz12+particleJ.quadrupole[QXZ]*guz13+particleJ.quadrupole[QYZ]*guz15))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
              +2.0*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
              +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
              +2.0*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx11+particleJ.quadrupole[QYY]*gqxx14+particleJ.quadrupole[QZZ]*gqxx16
              +2.0*(particleJ.quadrupole[QXY]*gqxx12+particleJ.quadrupole[QXZ]*gqxx13+particleJ.quadrupole[QYZ]*gqxx15))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy11+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqyy16
              +2.0*(particleJ.quadrupole[QXY]*gqyy12+particleJ.quadrupole[QXZ]*gqyy13+particleJ.quadrupole[QYZ]*gqyy15))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz11+particleJ.quadrupole[QYY]*gqzz14+particleJ.quadrupole[QZZ]*gqzz16
              +2.0*(particleJ.quadrupole[QXY]*gqzz12+particleJ.quadrupole[QXZ]*gqzz13+particleJ.quadrupole[QYZ]*gqzz15))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy11+particleJ.quadrupole[QYY]*gqxy14+particleJ.quadrupole[QZZ]*gqxy16
              +2.0*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxy13+particleJ.quadrupole[QYZ]*gqxy15))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz11+particleJ.quadrupole[QYY]*gqxz14+particleJ.quadrupole[QZZ]*gqxz16
              +2.0*(particleJ.quadrupole[QXY]*gqxz12+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqxz15))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz11+particleJ.quadrupole[QYY]*gqyz14+particleJ.quadrupole[QZZ]*gqyz16
              +2.0*(particleJ.quadrupole[QXY]*gqyz12+particleJ.quadrupole[QXZ]*gqyz13+particleJ.quadrupole[QYZ]*gqyz15)));

    dewkdx = particleI.charge*(particleJ.dipole[0]*gux2+particleJ.dipole[1]*guy2+particleJ.dipole[2]*guz2)
                      -particleJ.charge*(particleI.dipole[0]*gc5+particleI.dipole[1]*gc6+particleI.dipole[2]*gc7)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
              +2.0*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc11+particleI.quadrupole[QYY]*gc14+particleI.quadrupole[QZZ]*gc16
              +2.0*(particleI.quadrupole[QXY]*gc12+particleI.quadrupole[QXZ]*gc13+particleI.quadrupole[QYZ]*gc15))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
              +2.0*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux11+particleI.quadrupole[QYY]*gux14+particleI.quadrupole[QZZ]*gux16
              +2.0*(particleI.quadrupole[QXY]*gux12+particleI.quadrupole[QXZ]*gux13+particleI.quadrupole[QYZ]*gux15))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy11+particleI.quadrupole[QYY]*guy14+particleI.quadrupole[QZZ]*guy16
              +2.0*(particleI.quadrupole[QXY]*guy12+particleI.quadrupole[QXZ]*guy13+particleI.quadrupole[QYZ]*guy15))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz11+particleI.quadrupole[QYY]*guz14+particleI.quadrupole[QZZ]*guz16
              +2.0*(particleI.quadrupole[QXY]*guz12+particleI.quadrupole[QXZ]*guz13+particleI.quadrupole[QYZ]*guz15))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx11+particleJ.quadrupole[QYY]*gqyy11+particleJ.quadrupole[QZZ]*gqzz11
              +2.0*(particleJ.quadrupole[QXY]*gqxy11+particleJ.quadrupole[QXZ]*gqxz11+particleJ.quadrupole[QYZ]*gqyz11))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx14+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqzz14
              +2.0*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxz14+particleJ.quadrupole[QYZ]*gqyz14))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx16+particleJ.quadrupole[QYY]*gqyy16+particleJ.quadrupole[QZZ]*gqzz16
              +2.0*(particleJ.quadrupole[QXY]*gqxy16+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqyz16))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqyy12+particleJ.quadrupole[QZZ]*gqzz12
              +2.0*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxz12+particleJ.quadrupole[QYZ]*gqyz12))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqyy13+particleJ.quadrupole[QZZ]*gqzz13
              +2.0*(particleJ.quadrupole[QXY]*gqxy13+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqyz13))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +2.0*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15)));

    dedx = desymdx + 0.5*(dewidx + dewkdx);

    desymdy = particleI.charge*particleJ.charge*gc3
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux8+particleJ.dipole[1]*guy8+particleJ.dipole[2]*guz8)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9));

    dewidy = particleI.charge*(particleJ.dipole[0]*gc6+particleJ.dipole[1]*gc8+particleJ.dipole[2]*gc9)
                      -particleJ.charge*(particleI.dipole[0]*gux3+particleI.dipole[1]*guy3+particleI.dipole[2]*guz3)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc12+particleJ.quadrupole[QYY]*gc17+particleJ.quadrupole[QZZ]*gc19
              +2.0*(particleJ.quadrupole[QXY]*gc14+particleJ.quadrupole[QXZ]*gc15+particleJ.quadrupole[QYZ]*gc18))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
              +2.0*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux12+particleJ.quadrupole[QYY]*gux17+particleJ.quadrupole[QZZ]*gux19
              +2.0*(particleJ.quadrupole[QXY]*gux14+particleJ.quadrupole[QXZ]*gux15+particleJ.quadrupole[QYZ]*gux18))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy12+particleJ.quadrupole[QYY]*guy17+particleJ.quadrupole[QZZ]*guy19
              +2.0*(particleJ.quadrupole[QXY]*guy14+particleJ.quadrupole[QXZ]*guy15+particleJ.quadrupole[QYZ]*guy18))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz12+particleJ.quadrupole[QYY]*guz17+particleJ.quadrupole[QZZ]*guz19
              +2.0*(particleJ.quadrupole[QXY]*guz14+particleJ.quadrupole[QXZ]*guz15+particleJ.quadrupole[QYZ]*guz18))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
              +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
              +2.0*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
              +2.0*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqxx17+particleJ.quadrupole[QZZ]*gqxx19
              +2.0*(particleJ.quadrupole[QXY]*gqxx14+particleJ.quadrupole[QXZ]*gqxx15+particleJ.quadrupole[QYZ]*gqxx18))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy12+particleJ.quadrupole[QYY]*gqyy17+particleJ.quadrupole[QZZ]*gqyy19
              +2.0*(particleJ.quadrupole[QXY]*gqyy14+particleJ.quadrupole[QXZ]*gqyy15+particleJ.quadrupole[QYZ]*gqyy18))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz12+particleJ.quadrupole[QYY]*gqzz17+particleJ.quadrupole[QZZ]*gqzz19
              +2.0*(particleJ.quadrupole[QXY]*gqzz14+particleJ.quadrupole[QXZ]*gqzz15+particleJ.quadrupole[QYZ]*gqzz18))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy12+particleJ.quadrupole[QYY]*gqxy17+particleJ.quadrupole[QZZ]*gqxy19
              +2.0*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxy15+particleJ.quadrupole[QYZ]*gqxy18))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz12+particleJ.quadrupole[QYY]*gqxz17+particleJ.quadrupole[QZZ]*gqxz19
              +2.0*(particleJ.quadrupole[QXY]*gqxz14+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqxz18))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz12+particleJ.quadrupole[QYY]*gqyz17+particleJ.quadrupole[QZZ]*gqyz19
              +2.0*(particleJ.quadrupole[QXY]*gqyz14+particleJ.quadrupole[QXZ]*gqyz15+particleJ.quadrupole[QYZ]*gqyz18)));

    dewkdy = particleI.charge*(particleJ.dipole[0]*gux3+particleJ.dipole[1]*guy3+particleJ.dipole[2]*guz3)
                      -particleJ.charge*(particleI.dipole[0]*gc6+particleI.dipole[1]*gc8+particleI.dipole[2]*gc9)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
              +2.0*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc12+particleI.quadrupole[QYY]*gc17+particleI.quadrupole[QZZ]*gc19
              +2.0*(particleI.quadrupole[QXY]*gc14+particleI.quadrupole[QXZ]*gc15+particleI.quadrupole[QYZ]*gc18))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
              +2.0*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux12+particleI.quadrupole[QYY]*gux17+particleI.quadrupole[QZZ]*gux19
              +2.0*(particleI.quadrupole[QXY]*gux14+particleI.quadrupole[QXZ]*gux15+particleI.quadrupole[QYZ]*gux18))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy12+particleI.quadrupole[QYY]*guy17+particleI.quadrupole[QZZ]*guy19
              +2.0*(particleI.quadrupole[QXY]*guy14+particleI.quadrupole[QXZ]*guy15+particleI.quadrupole[QYZ]*guy18))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz12+particleI.quadrupole[QYY]*guz17+particleI.quadrupole[QZZ]*guz19
              +2.0*(particleI.quadrupole[QXY]*guz14+particleI.quadrupole[QXZ]*guz15+particleI.quadrupole[QYZ]*guz18))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqyy12+particleJ.quadrupole[QZZ]*gqzz12
              +2.0*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxz12+particleJ.quadrupole[QYZ]*gqyz12))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx17+particleJ.quadrupole[QYY]*gqyy17+particleJ.quadrupole[QZZ]*gqzz17
              +2.0*(particleJ.quadrupole[QXY]*gqxy17+particleJ.quadrupole[QXZ]*gqxz17+particleJ.quadrupole[QYZ]*gqyz17))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx19+particleJ.quadrupole[QYY]*gqyy19+particleJ.quadrupole[QZZ]*gqzz19
              +2.0*(particleJ.quadrupole[QXY]*gqxy19+particleJ.quadrupole[QXZ]*gqxz19+particleJ.quadrupole[QYZ]*gqyz19))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx14+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqzz14
              +2.0*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxz14+particleJ.quadrupole[QYZ]*gqyz14))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +2.0*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx18+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqzz18
              +2.0*(particleJ.quadrupole[QXY]*gqxy18+particleJ.quadrupole[QXZ]*gqxz18+particleJ.quadrupole[QYZ]*gqyz18)));

    dedy = desymdy + 0.5*(dewidy + dewkdy);

    desymdz = particleI.charge*particleJ.charge*gc4
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux10+particleJ.dipole[1]*guy10+particleJ.dipole[2]*guz10));

    dewidz = particleI.charge*(particleJ.dipole[0]*gc7+particleJ.dipole[1]*gc9+particleJ.dipole[2]*gc10)
                      -particleJ.charge*(particleI.dipole[0]*gux4+particleI.dipole[1]*guy4+particleI.dipole[2]*guz4)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc13+particleJ.quadrupole[QYY]*gc18+particleJ.quadrupole[QZZ]*gc20
              +2.0*(particleJ.quadrupole[QXY]*gc15+particleJ.quadrupole[QXZ]*gc16+particleJ.quadrupole[QYZ]*gc19))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
              +2.0*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux13+particleJ.quadrupole[QYY]*gux18+particleJ.quadrupole[QZZ]*gux20
              +2.0*(particleJ.quadrupole[QXY]*gux15+particleJ.quadrupole[QXZ]*gux16+particleJ.quadrupole[QYZ]*gux19))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy13+particleJ.quadrupole[QYY]*guy18+particleJ.quadrupole[QZZ]*guy20
              +2.0*(particleJ.quadrupole[QXY]*guy15+particleJ.quadrupole[QXZ]*guy16+particleJ.quadrupole[QYZ]*guy19))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz13+particleJ.quadrupole[QYY]*guz18+particleJ.quadrupole[QZZ]*guz20
              +2.0*(particleJ.quadrupole[QXY]*guz15+particleJ.quadrupole[QXZ]*guz16+particleJ.quadrupole[QYZ]*guz19))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
              +2.0*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
              +2.0*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
              +2.0*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqxx18+particleJ.quadrupole[QZZ]*gqxx20
              +2.0*(particleJ.quadrupole[QXY]*gqxx15+particleJ.quadrupole[QXZ]*gqxx16+particleJ.quadrupole[QYZ]*gqxx19))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy13+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqyy20
              +2.0*(particleJ.quadrupole[QXY]*gqyy15+particleJ.quadrupole[QXZ]*gqyy16+particleJ.quadrupole[QYZ]*gqyy19))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz13+particleJ.quadrupole[QYY]*gqzz18+particleJ.quadrupole[QZZ]*gqzz20
              +2.0*(particleJ.quadrupole[QXY]*gqzz15+particleJ.quadrupole[QXZ]*gqzz16+particleJ.quadrupole[QYZ]*gqzz19))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy13+particleJ.quadrupole[QYY]*gqxy18+particleJ.quadrupole[QZZ]*gqxy20
              +2.0*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxy16+particleJ.quadrupole[QYZ]*gqxy19))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz13+particleJ.quadrupole[QYY]*gqxz18+particleJ.quadrupole[QZZ]*gqxz20
              +2.0*(particleJ.quadrupole[QXY]*gqxz15+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqxz19))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz13+particleJ.quadrupole[QYY]*gqyz18+particleJ.quadrupole[QZZ]*gqyz20
              +2.0*(particleJ.quadrupole[QXY]*gqyz15+particleJ.quadrupole[QXZ]*gqyz16+particleJ.quadrupole[QYZ]*gqyz19)));

    dewkdz = particleI.charge*(particleJ.dipole[0]*gux4+particleJ.dipole[1]*guy4+particleJ.dipole[2]*guz4)
                      -particleJ.charge*(particleI.dipole[0]*gc7+particleI.dipole[1]*gc9+particleI.dipole[2]*gc10)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
              +2.0*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc13+particleI.quadrupole[QYY]*gc18+particleI.quadrupole[QZZ]*gc20
              +2.0*(particleI.quadrupole[QXY]*gc15+particleI.quadrupole[QXZ]*gc16+particleI.quadrupole[QYZ]*gc19))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
              +2.0*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux13+particleI.quadrupole[QYY]*gux18+particleI.quadrupole[QZZ]*gux20
              +2.0*(particleI.quadrupole[QXY]*gux15+particleI.quadrupole[QXZ]*gux16+particleI.quadrupole[QYZ]*gux19))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy13+particleI.quadrupole[QYY]*guy18+particleI.quadrupole[QZZ]*guy20
              +2.0*(particleI.quadrupole[QXY]*guy15+particleI.quadrupole[QXZ]*guy16+particleI.quadrupole[QYZ]*guy19))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz13+particleI.quadrupole[QYY]*guz18+particleI.quadrupole[QZZ]*guz20
              +2.0*(particleI.quadrupole[QXY]*guz15+particleI.quadrupole[QXZ]*guz16+particleI.quadrupole[QYZ]*guz19))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqyy13+particleJ.quadrupole[QZZ]*gqzz13
              +2.0*(particleJ.quadrupole[QXY]*gqxy13+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqyz13))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx18+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqzz18
              +2.0*(particleJ.quadrupole[QXY]*gqxy18+particleJ.quadrupole[QXZ]*gqxz18+particleJ.quadrupole[QYZ]*gqyz18))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx20+particleJ.quadrupole[QYY]*gqyy20+particleJ.quadrupole[QZZ]*gqzz20
              +2.0*(particleJ.quadrupole[QXY]*gqxy20+particleJ.quadrupole[QXZ]*gqxz20+particleJ.quadrupole[QYZ]*gqyz20))
       + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +2.0*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx16+particleJ.quadrupole[QYY]*gqyy16+particleJ.quadrupole[QZZ]*gqzz16
              +2.0*(particleJ.quadrupole[QXY]*gqxy16+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqyz16))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx19+particleJ.quadrupole[QYY]*gqyy19+particleJ.quadrupole[QZZ]*gqzz19
              +2.0*(particleJ.quadrupole[QXY]*gqxy19+particleJ.quadrupole[QXZ]*gqxz19+particleJ.quadrupole[QYZ]*gqyz19)));

    dedz = desymdz + 0.5*(dewidz + dewkdz);

    desymdr = particleI.charge*particleJ.charge*gc21
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux22+particleJ.dipole[1]*guy22+particleJ.dipole[2]*guz22)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux23+particleJ.dipole[1]*guy23+particleJ.dipole[2]*guz23)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux24+particleJ.dipole[1]*guy24+particleJ.dipole[2]*guz24));

    dewidr = particleI.charge*(particleJ.dipole[0]*gc22+particleJ.dipole[1]*gc23+particleJ.dipole[2]*gc24)
                      -particleJ.charge*(particleI.dipole[0]*gux21+particleI.dipole[1]*guy21+particleI.dipole[2]*guz21)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc25+particleJ.quadrupole[QYY]*gc28+particleJ.quadrupole[QZZ]*gc30
              +2.0*(particleJ.quadrupole[QXY]*gc26+particleJ.quadrupole[QXZ]*gc27+particleJ.quadrupole[QYZ]*gc29))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx21+particleI.quadrupole[QYY]*gqyy21+particleI.quadrupole[QZZ]*gqzz21
              +2.0*(particleI.quadrupole[QXY]*gqxy21+particleI.quadrupole[QXZ]*gqxz21+particleI.quadrupole[QYZ]*gqyz21))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux25+particleJ.quadrupole[QYY]*gux28+particleJ.quadrupole[QZZ]*gux30
              +2.0*(particleJ.quadrupole[QXY]*gux26+particleJ.quadrupole[QXZ]*gux27+particleJ.quadrupole[QYZ]*gux29))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy25+particleJ.quadrupole[QYY]*guy28+particleJ.quadrupole[QZZ]*guy30
              +2.0*(particleJ.quadrupole[QXY]*guy26+particleJ.quadrupole[QXZ]*guy27+particleJ.quadrupole[QYZ]*guy29))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz25+particleJ.quadrupole[QYY]*guz28+particleJ.quadrupole[QZZ]*guz30
              +2.0*(particleJ.quadrupole[QXY]*guz26+particleJ.quadrupole[QXZ]*guz27+particleJ.quadrupole[QYZ]*guz29))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx22+particleI.quadrupole[QYY]*gqyy22+particleI.quadrupole[QZZ]*gqzz22
              +2.0*(particleI.quadrupole[QXY]*gqxy22+particleI.quadrupole[QXZ]*gqxz22+particleI.quadrupole[QYZ]*gqyz22))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx23+particleI.quadrupole[QYY]*gqyy23+particleI.quadrupole[QZZ]*gqzz23
              +2.0*(particleI.quadrupole[QXY]*gqxy23+particleI.quadrupole[QXZ]*gqxz23+particleI.quadrupole[QYZ]*gqyz23))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx24+particleI.quadrupole[QYY]*gqyy24+particleI.quadrupole[QZZ]*gqzz24
              +2.0*(particleI.quadrupole[QXY]*gqxy24+particleI.quadrupole[QXZ]*gqxz24+particleI.quadrupole[QYZ]*gqyz24))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx25+particleJ.quadrupole[QYY]*gqxx28+particleJ.quadrupole[QZZ]*gqxx30
              +2.0*(particleJ.quadrupole[QXY]*gqxx26+particleJ.quadrupole[QXZ]*gqxx27+particleJ.quadrupole[QYZ]*gqxx29))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy25+particleJ.quadrupole[QYY]*gqyy28+particleJ.quadrupole[QZZ]*gqyy30
              +2.0*(particleJ.quadrupole[QXY]*gqyy26+particleJ.quadrupole[QXZ]*gqyy27+particleJ.quadrupole[QYZ]*gqyy29))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz25+particleJ.quadrupole[QYY]*gqzz28+particleJ.quadrupole[QZZ]*gqzz30
              +2.0*(particleJ.quadrupole[QXY]*gqzz26+particleJ.quadrupole[QXZ]*gqzz27+particleJ.quadrupole[QYZ]*gqzz29))
              + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy25+particleJ.quadrupole[QYY]*gqxy28+particleJ.quadrupole[QZZ]*gqxy30
              +2.0*(particleJ.quadrupole[QXY]*gqxy26+particleJ.quadrupole[QXZ]*gqxy27+particleJ.quadrupole[QYZ]*gqxy29))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz25+particleJ.quadrupole[QYY]*gqxz28+particleJ.quadrupole[QZZ]*gqxz30
              +2.0*(particleJ.quadrupole[QXY]*gqxz26+particleJ.quadrupole[QXZ]*gqxz27+particleJ.quadrupole[QYZ]*gqxz29))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz25+particleJ.quadrupole[QYY]*gqyz28+particleJ.quadrupole[QZZ]*gqyz30
              +2.0*(particleJ.quadrupole[QXY]*gqyz26+particleJ.quadrupole[QXZ]*gqyz27+particleJ.quadrupole[QYZ]*gqyz29)));

    dewkdr = particleI.charge*(particleJ.dipole[0]*gux21+particleJ.dipole[1]*guy21+particleJ.dipole[2]*guz21)
                      -particleJ.charge*(particleI.dipole[0]*gc22+particleI.dipole[1]*gc23+particleI.dipole[2]*gc24)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx21+particleJ.quadrupole[QYY]*gqyy21+particleJ.quadrupole[QZZ]*gqzz21
              +2.0*(particleJ.quadrupole[QXY]*gqxy21+particleJ.quadrupole[QXZ]*gqxz21+particleJ.quadrupole[QYZ]*gqyz21))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc25+particleI.quadrupole[QYY]*gc28+particleI.quadrupole[QZZ]*gc30
              +2.0*(particleI.quadrupole[QXY]*gc26+particleI.quadrupole[QXZ]*gc27+particleI.quadrupole[QYZ]*gc29))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx22+particleJ.quadrupole[QYY]*gqyy22+particleJ.quadrupole[QZZ]*gqzz22
              +2.0*(particleJ.quadrupole[QXY]*gqxy22+particleJ.quadrupole[QXZ]*gqxz22+particleJ.quadrupole[QYZ]*gqyz22))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx23+particleJ.quadrupole[QYY]*gqyy23+particleJ.quadrupole[QZZ]*gqzz23
              +2.0*(particleJ.quadrupole[QXY]*gqxy23+particleJ.quadrupole[QXZ]*gqxz23+particleJ.quadrupole[QYZ]*gqyz23))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx24+particleJ.quadrupole[QYY]*gqyy24+particleJ.quadrupole[QZZ]*gqzz24
              +2.0*(particleJ.quadrupole[QXY]*gqxy24+particleJ.quadrupole[QXZ]*gqxz24+particleJ.quadrupole[QYZ]*gqyz24))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux25+particleI.quadrupole[QYY]*gux28+particleI.quadrupole[QZZ]*gux30
              +2.0*(particleI.quadrupole[QXY]*gux26+particleI.quadrupole[QXZ]*gux27+particleI.quadrupole[QYZ]*gux29))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy25+particleI.quadrupole[QYY]*guy28+particleI.quadrupole[QZZ]*guy30
              +2.0*(particleI.quadrupole[QXY]*guy26+particleI.quadrupole[QXZ]*guy27+particleI.quadrupole[QYZ]*guy29))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz25+particleI.quadrupole[QYY]*guz28+particleI.quadrupole[QZZ]*guz30
              +2.0*(particleI.quadrupole[QXY]*guz26+particleI.quadrupole[QXZ]*guz27+particleI.quadrupole[QYZ]*guz29))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx25+particleJ.quadrupole[QYY]*gqyy25+particleJ.quadrupole[QZZ]*gqzz25
              +2.0*(particleJ.quadrupole[QXY]*gqxy25+particleJ.quadrupole[QXZ]*gqxz25+particleJ.quadrupole[QYZ]*gqyz25))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx28+particleJ.quadrupole[QYY]*gqyy28+particleJ.quadrupole[QZZ]*gqzz28
              +2.0*(particleJ.quadrupole[QXY]*gqxy28+particleJ.quadrupole[QXZ]*gqxz28+particleJ.quadrupole[QYZ]*gqyz28))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx30+particleJ.quadrupole[QYY]*gqyy30+particleJ.quadrupole[QZZ]*gqzz30
              +2.0*(particleJ.quadrupole[QXY]*gqxy30+particleJ.quadrupole[QXZ]*gqxz30+particleJ.quadrupole[QYZ]*gqyz30))
              + 2.0*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx26+particleJ.quadrupole[QYY]*gqyy26+particleJ.quadrupole[QZZ]*gqzz26
              +2.0*(particleJ.quadrupole[QXY]*gqxy26+particleJ.quadrupole[QXZ]*gqxz26+particleJ.quadrupole[QYZ]*gqyz26))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx27+particleJ.quadrupole[QYY]*gqyy27+particleJ.quadrupole[QZZ]*gqzz27
              +2.0*(particleJ.quadrupole[QXY]*gqxy27+particleJ.quadrupole[QXZ]*gqxz27+particleJ.quadrupole[QYZ]*gqyz27))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx29+particleJ.quadrupole[QYY]*gqyy29+particleJ.quadrupole[QZZ]*gqzz29
              +2.0*(particleJ.quadrupole[QXY]*gqxy29+particleJ.quadrupole[QXZ]*gqxz29+particleJ.quadrupole[QYZ]*gqyz29)));

    dsumdr = desymdr + 0.5*(dewidr + dewkdr);
    drbi = _bornRadii[jIndex]*dsumdr;
    drbk = _bornRadii[iIndex]*dsumdr;

    // torque on permanent dipoles due to permanent reaction field

    RealOpenMM trq1   = 0.0;
    RealOpenMM trq2   = 0.0;
    RealOpenMM trq3   = 0.0;

    RealOpenMM trq_k1 = 0.0;
    RealOpenMM trq_k2 = 0.0;
    RealOpenMM trq_k3 = 0.0;

    if (xr != 0.0 || yr != 0.0 || zr != 0.0)
    {

        RealOpenMM fid1 = particleJ.dipole[0]*gux2 + particleJ.dipole[1]*gux3 + particleJ.dipole[2]*gux4
                + 0.5*(particleJ.charge*gux1+particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
                      +2.0*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9)
                      +particleJ.charge*gc2+particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
                      +2.0*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2));

        RealOpenMM fid2 = particleJ.dipole[0]*guy2 + particleJ.dipole[1]*guy3 + particleJ.dipole[2]*guy4
                + 0.5*(particleJ.charge*guy1+particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
                      +2.0*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9)
                      +particleJ.charge*gc3+particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
                      +2.0*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3));

        RealOpenMM fid3 = particleJ.dipole[0]*guz2 + particleJ.dipole[1]*guz3 + particleJ.dipole[2]*guz4
                + 0.5*(particleJ.charge*guz1+particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
                      +2.0*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9)
                      +particleJ.charge*gc4+particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
                      +2.0*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4));

        RealOpenMM fkd1 = particleI.dipole[0]*gux2 + particleI.dipole[1]*gux3 + particleI.dipole[2]*gux4
                - 0.5*(particleI.charge*gux1+particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
                      +2.0*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9)
                      +particleI.charge*gc2+particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
                      +2.0*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2));

        RealOpenMM fkd2 = particleI.dipole[0]*guy2 + particleI.dipole[1]*guy3 + particleI.dipole[2]*guy4
                - 0.5*(particleI.charge*guy1+particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
                      +2.0*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9)
                      +particleI.charge*gc3+particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
                      +2.0*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3));

        RealOpenMM fkd3 = particleI.dipole[0]*guz2 + particleI.dipole[1]*guz3 + particleI.dipole[2]*guz4
                - 0.5*(particleI.charge*guz1+particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
                      +2.0*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9)
                      +particleI.charge*gc4+particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
                      +2.0*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4));

        trq1    = particleI.dipole[1]*fid3 - particleI.dipole[2]*fid2;
        trq2    = particleI.dipole[2]*fid1 - particleI.dipole[0]*fid3;
        trq3    = particleI.dipole[0]*fid2 - particleI.dipole[1]*fid1;

        trq_k1  = particleJ.dipole[1]*fkd3 - particleJ.dipole[2]*fkd2;
        trq_k2  = particleJ.dipole[2]*fkd1 - particleJ.dipole[0]*fkd3;
        trq_k3  = particleJ.dipole[0]*fkd2 - particleJ.dipole[1]*fkd1;

        // torque on quadrupoles due to permanent reaction field gradient

        RealOpenMM fidg11 =
                - 0.5*(particleJ.charge*gqxx1+particleJ.dipole[0]*gqxx2+particleJ.dipole[1]*gqxx3+particleJ.dipole[2]*gqxx4
                      +particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqxx8+particleJ.quadrupole[QZZ]*gqxx10
                      +2.0*(particleJ.quadrupole[QXY]*gqxx6+particleJ.quadrupole[QXZ]*gqxx7+particleJ.quadrupole[QYZ]*gqxx9)
                      +particleJ.charge*gc5+particleJ.dipole[0]*gux5+particleJ.dipole[1]*guy5+particleJ.dipole[2]*guz5
                      +particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
                      +2.0*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5));

        RealOpenMM fidg12 =
                - 0.5*(particleJ.charge*gqxy1+particleJ.dipole[0]*gqxy2+particleJ.dipole[1]*gqxy3+particleJ.dipole[2]*gqxy4
                      +particleJ.quadrupole[QXX]*gqxy5+particleJ.quadrupole[QYY]*gqxy8+particleJ.quadrupole[QZZ]*gqxy10
                      +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxy7+particleJ.quadrupole[QYZ]*gqxy9)
                      +particleJ.charge*gc6+particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6
                      +particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                      +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6));

        RealOpenMM fidg13 =
                - 0.5*(particleJ.charge*gqxz1+particleJ.dipole[0]*gqxz2+particleJ.dipole[1]*gqxz3+particleJ.dipole[2]*gqxz4
                      +particleJ.quadrupole[QXX]*gqxz5+particleJ.quadrupole[QYY]*gqxz8+particleJ.quadrupole[QZZ]*gqxz10
                      +2.0*(particleJ.quadrupole[QXY]*gqxz6+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqxz9)
                      +particleJ.charge*gc7+particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7
                      +particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                      +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7));

        RealOpenMM fidg22 =
                - 0.5*(particleJ.charge*gqyy1+particleJ.dipole[0]*gqyy2+particleJ.dipole[1]*gqyy3+particleJ.dipole[2]*gqyy4
                      +particleJ.quadrupole[QXX]*gqyy5+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqyy10
                      +2.0*(particleJ.quadrupole[QXY]*gqyy6+particleJ.quadrupole[QXZ]*gqyy7+particleJ.quadrupole[QYZ]*gqyy9)
                      +particleJ.charge*gc8+particleJ.dipole[0]*gux8+particleJ.dipole[1]*guy8+particleJ.dipole[2]*guz8
                      +particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
                      +2.0*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8));

        RealOpenMM fidg23 =
                - 0.5*(particleJ.charge*gqyz1+particleJ.dipole[0]*gqyz2+particleJ.dipole[1]*gqyz3+particleJ.dipole[2]*gqyz4
                      +particleJ.quadrupole[QXX]*gqyz5+particleJ.quadrupole[QYY]*gqyz8+particleJ.quadrupole[QZZ]*gqyz10
                      +2.0*(particleJ.quadrupole[QXY]*gqyz6+particleJ.quadrupole[QXZ]*gqyz7+particleJ.quadrupole[QYZ]*gqyz9)
                      +particleJ.charge*gc9+particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9
                      +particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                      +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9));

        RealOpenMM fidg33 =
                - 0.5*(particleJ.charge*gqzz1+particleJ.dipole[0]*gqzz2+particleJ.dipole[1]*gqzz3+particleJ.dipole[2]*gqzz4
                      +particleJ.quadrupole[QXX]*gqzz5+particleJ.quadrupole[QYY]*gqzz8+particleJ.quadrupole[QZZ]*gqzz10
                      +2.0*(particleJ.quadrupole[QXY]*gqzz6+particleJ.quadrupole[QXZ]*gqzz7+particleJ.quadrupole[QYZ]*gqzz9)
                      +particleJ.charge*gc10+particleJ.dipole[0]*gux10+particleJ.dipole[1]*guy10+particleJ.dipole[2]*guz10
                      +particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
                   +2.0*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10));

        RealOpenMM fidg21 = fidg12;
        RealOpenMM fidg31 = fidg13;
        RealOpenMM fidg32 = fidg23;

        RealOpenMM fkdg11 =
                - 0.5*(particleI.charge*gqxx1-particleI.dipole[0]*gqxx2-particleI.dipole[1]*gqxx3-particleI.dipole[2] *gqxx4
                      +particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqxx8+particleI.quadrupole[QZZ]*gqxx10
                      +2.0*(particleI.quadrupole[QXY]*gqxx6+particleI.quadrupole[QXZ]*gqxx7+particleI.quadrupole[QYZ]*gqxx9)
                      +particleI.charge*gc5-particleI.dipole[0]*gux5-particleI.dipole[1]*guy5-particleI.dipole[2]*guz5
                      +particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
                      +2.0*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5));

        RealOpenMM fkdg12 =
                - 0.5*(particleI.charge*gqxy1-particleI.dipole[0]*gqxy2-particleI.dipole[1]*gqxy3-particleI.dipole[2]*gqxy4
                      +particleI.quadrupole[QXX]*gqxy5+particleI.quadrupole[QYY]*gqxy8+particleI.quadrupole[QZZ]*gqxy10
                      +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxy7+particleI.quadrupole[QYZ]*gqxy9)
                      +particleI.charge*gc6-particleI.dipole[0]*gux6-particleI.dipole[1]*guy6-particleI.dipole[2]*guz6
                      +particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                      +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6));

        RealOpenMM fkdg13 =
                - 0.5*(particleI.charge*gqxz1-particleI.dipole[0]*gqxz2-particleI.dipole[1]*gqxz3-particleI.dipole[2]*gqxz4
                      +particleI.quadrupole[QXX]*gqxz5+particleI.quadrupole[QYY]*gqxz8+particleI.quadrupole[QZZ]*gqxz10
                      +2.0*(particleI.quadrupole[QXY]*gqxz6+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqxz9)
                      +particleI.charge*gc7-particleI.dipole[0]*gux7-particleI.dipole[1]*guy7-particleI.dipole[2]*guz7
                      +particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                      +2.0*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7));

        RealOpenMM fkdg22 =
                - 0.5*(particleI.charge*gqyy1-particleI.dipole[0]*gqyy2-particleI.dipole[1]*gqyy3-particleI.dipole[2]*gqyy4
                      +particleI.quadrupole[QXX]*gqyy5+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqyy10
                      +2.0*(particleI.quadrupole[QXY]*gqyy6+particleI.quadrupole[QXZ]*gqyy7+particleI.quadrupole[QYZ]*gqyy9)
                      +particleI.charge*gc8-particleI.dipole[0]*gux8-particleI.dipole[1]*guy8-particleI.dipole[2]*guz8
                      +particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
                      +2.0*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8));

        RealOpenMM fkdg23 =
                - 0.5*(particleI.charge*gqyz1-particleI.dipole[0]*gqyz2-particleI.dipole[1]*gqyz3-particleI.dipole[2]*gqyz4
                      +particleI.quadrupole[QXX]*gqyz5+particleI.quadrupole[QYY]*gqyz8+particleI.quadrupole[QZZ]*gqyz10
                      +2.0*(particleI.quadrupole[QXY]*gqyz6+particleI.quadrupole[QXZ]*gqyz7+particleI.quadrupole[QYZ]*gqyz9)
                      +particleI.charge*gc9-particleI.dipole[0]*gux9-particleI.dipole[1]*guy9-particleI.dipole[2]*guz9
                      +particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                      +2.0*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9));
        RealOpenMM fkdg33 =
                - 0.5*(particleI.charge*gqzz1-particleI.dipole[0]*gqzz2-particleI.dipole[1]*gqzz3-particleI.dipole[2]*gqzz4
                      +particleI.quadrupole[QXX]*gqzz5+particleI.quadrupole[QYY]*gqzz8+particleI.quadrupole[QZZ]*gqzz10
                      +2.0*(particleI.quadrupole[QXY]*gqzz6+particleI.quadrupole[QXZ]*gqzz7+particleI.quadrupole[QYZ]*gqzz9)
                      +particleI.charge*gc10-particleI.dipole[0]*gux10-particleI.dipole[1]*guy10-particleI.dipole[2]*guz10
                      +particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
                    +2.0*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10));

        RealOpenMM fkdg21 = fkdg12;
        RealOpenMM fkdg31 = fkdg13;
        RealOpenMM fkdg32 = fkdg23;

        trq1   += 2.0* (particleI.quadrupole[QXY]*fidg13+particleI.quadrupole[QYY]*fidg23+particleI.quadrupole[QYZ]*fidg33
                           -particleI.quadrupole[QXZ]*fidg12-particleI.quadrupole[QYZ]*fidg22-particleI.quadrupole[QZZ]*fidg32);

        trq2   += 2.0*(particleI.quadrupole[QXZ]*fidg11+particleI.quadrupole[QYZ]*fidg21+particleI.quadrupole[QZZ]*fidg31
                         -particleI.quadrupole[QXX]*fidg13-particleI.quadrupole[QXY]*fidg23-particleI.quadrupole[QXZ]*fidg33);

        trq3   += 2.0*(particleI.quadrupole[QXX]*fidg12+particleI.quadrupole[QXY]*fidg22+particleI.quadrupole[QXZ]*fidg32
                         -particleI.quadrupole[QXY]*fidg11-particleI.quadrupole[QYY]*fidg21-particleI.quadrupole[QYZ]*fidg31);

        trq_k1 += 2.0*
                          (particleJ.quadrupole[QXY]*fkdg13+particleJ.quadrupole[QYY]*fkdg23+particleJ.quadrupole[QYZ]*fkdg33
                          -particleJ.quadrupole[QXZ]*fkdg12-particleJ.quadrupole[QYZ]*fkdg22-particleJ.quadrupole[QZZ]*fkdg32);

        trq_k2 += 2.0*
                          (particleJ.quadrupole[QXZ]*fkdg11+particleJ.quadrupole[QYZ]*fkdg21+particleJ.quadrupole[QZZ]*fkdg31
                          -particleJ.quadrupole[QXX]*fkdg13-particleJ.quadrupole[QXY]*fkdg23-particleJ.quadrupole[QXZ]*fkdg33);

        trq_k3 += 2.0*
                          (particleJ.quadrupole[QXX]*fkdg12+particleJ.quadrupole[QXY]*fkdg22+particleJ.quadrupole[QXZ]*fkdg32
                          -particleJ.quadrupole[QXY]*fkdg11-particleJ.quadrupole[QYY]*fkdg21-particleJ.quadrupole[QYZ]*fkdg31);
    }

    // electrostatic solvation energy of the permanent multipoles in
    // the GK reaction potential of the induced dipoles

    esymi =              -particleI.dipole[0]*(_inducedDipoleS[jIndex][0]*gux2+_inducedDipoleS[jIndex][1]*guy2+_inducedDipoleS[jIndex][2]*guz2)
                        - particleI.dipole[1]*(_inducedDipoleS[jIndex][0]*gux3+_inducedDipoleS[jIndex][1]*guy3+_inducedDipoleS[jIndex][2]*guz3)
                        - particleI.dipole[2]*(_inducedDipoleS[jIndex][0]*gux4+_inducedDipoleS[jIndex][1]*guy4+_inducedDipoleS[jIndex][2]*guz4)
                        - particleJ.dipole[0]*(_inducedDipoleS[iIndex][0]*gux2+_inducedDipoleS[iIndex][1]*guy2+_inducedDipoleS[iIndex][2]*guz2)
                        - particleJ.dipole[1]*(_inducedDipoleS[iIndex][0]*gux3+_inducedDipoleS[iIndex][1]*guy3+_inducedDipoleS[iIndex][2]*guz3)
                        - particleJ.dipole[2]*(_inducedDipoleS[iIndex][0]*gux4+_inducedDipoleS[iIndex][1]*guy4+_inducedDipoleS[iIndex][2]*guz4);

    ewii = particleI.charge*(_inducedDipoleS[jIndex][0]*gc2+_inducedDipoleS[jIndex][1]*gc3+_inducedDipoleS[jIndex][2]*gc4)
                      - particleJ.charge*(_inducedDipoleS[iIndex][0]*gux1+_inducedDipoleS[iIndex][1]*guy1+_inducedDipoleS[iIndex][2]*guz1)
                      - _inducedDipoleS[iIndex][0]*(particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
                     +2.0*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9))
                      - _inducedDipoleS[iIndex][1]*(particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
                     +2.0*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9))
                      - _inducedDipoleS[iIndex][2]*(particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
                     +2.0*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9))
                      + _inducedDipoleS[jIndex][0]*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
                     +2.0*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
                      + _inducedDipoleS[jIndex][1]*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
                     +2.0*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
                      + _inducedDipoleS[jIndex][2]*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
                     +2.0*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4));

    ewki = particleI.charge*(_inducedDipoleS[jIndex][0]*gux1+_inducedDipoleS[jIndex][1]*guy1+_inducedDipoleS[jIndex][2]*guz1)
                      - particleJ.charge*(_inducedDipoleS[iIndex][0]*gc2+_inducedDipoleS[iIndex][1]*gc3+_inducedDipoleS[iIndex][2]*gc4)
                      - _inducedDipoleS[iIndex][0]*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
                     +2.0*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
                      - _inducedDipoleS[iIndex][1]*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
                     +2.0*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
                      - _inducedDipoleS[iIndex][2]*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
                     +2.0*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
                      + _inducedDipoleS[jIndex][0]*(particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
                     +2.0*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9))
                      + _inducedDipoleS[jIndex][1]*(particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
                     +2.0*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9))
                      + _inducedDipoleS[jIndex][2]*(particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
                     +2.0*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9));

    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    RealOpenMM dpsymdx =       - particleI.dipole[0]*(sxk*gux5+syk*guy5+szk*guz5)
                          - particleI.dipole[1]*(sxk*gux6+syk*guy6+szk*guz6)
                          - particleI.dipole[2]*(sxk*gux7+syk*guy7+szk*guz7)

                          - particleJ.dipole[0]*(sxi*gux5+syi*guy5+szi*guz5)
                          - particleJ.dipole[1]*(sxi*gux6+syi*guy6+szi*guz6)
                          - particleJ.dipole[2]*(sxi*gux7+syi*guy7+szi*guz7);

    dpwidx = particleI.charge*(sxk*gc5+syk*gc6+szk*gc7)
                        - particleJ.charge*(sxi*gux2+syi*guy2+szi*guz2)
                      - sxi*(particleJ.quadrupole[QXX]*gux11+particleJ.quadrupole[QYY]*gux14+particleJ.quadrupole[QZZ]*gux16
                     +2.0*(particleJ.quadrupole[QXY]*gux12+particleJ.quadrupole[QXZ]*gux13+particleJ.quadrupole[QYZ]*gux15))
                      - syi*(particleJ.quadrupole[QXX]*guy11+particleJ.quadrupole[QYY]*guy14+particleJ.quadrupole[QZZ]*guy16
                     +2.0*(particleJ.quadrupole[QXY]*guy12+particleJ.quadrupole[QXZ]*guy13+particleJ.quadrupole[QYZ]*guy15))
                      - szi*(particleJ.quadrupole[QXX]*guz11+particleJ.quadrupole[QYY]*guz14+particleJ.quadrupole[QZZ]*guz16
                     +2.0*(particleJ.quadrupole[QXY]*guz12+particleJ.quadrupole[QXZ]*guz13+particleJ.quadrupole[QYZ]*guz15))
                      + sxk*(particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
                     +2.0*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5))
                      + syk*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                     +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
                      + szk*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                     +2.0*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7));

    dpwkdx = particleI.charge*(sxk*gux2+syk*guy2+szk*guz2)
                        - particleJ.charge*(sxi*gc5+syi*gc6+szi*gc7)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
                     +2.0*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
                      - syi*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                     +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
                      - szi*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                     +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
                      + sxk*(particleI.quadrupole[QXX]*gux11+particleI.quadrupole[QYY]*gux14+particleI.quadrupole[QZZ]*gux16
                     +2.0*(particleI.quadrupole[QXY]*gux12+particleI.quadrupole[QXZ]*gux13+particleI.quadrupole[QYZ]*gux15))
                      + syk*(particleI.quadrupole[QXX]*guy11+particleI.quadrupole[QYY]*guy14+particleI.quadrupole[QZZ]*guy16
                     +2.0*(particleI.quadrupole[QXY]*guy12+particleI.quadrupole[QXZ]*guy13+particleI.quadrupole[QYZ]*guy15))
                      + szk*(particleI.quadrupole[QXX]*guz11+particleI.quadrupole[QYY]*guz14+particleI.quadrupole[QZZ]*guz16
                     +2.0*(particleI.quadrupole[QXY]*guz12+particleI.quadrupole[QXZ]*guz13+particleI.quadrupole[QYZ]*guz15));

    dpdx = 0.5*(dpsymdx + 0.5*(dpwidx + dpwkdx));

    dpsymdy = -particleI.dipole[0]*(sxk*gux6+syk*guy6+szk*guz6)
                          - particleI.dipole[1]*(sxk*gux8+syk*guy8+szk*guz8)
                          - particleI.dipole[2]*(sxk*gux9+syk*guy9+szk*guz9)
                          - particleJ.dipole[0]*(sxi*gux6+syi*guy6+szi*guz6)
                          - particleJ.dipole[1]*(sxi*gux8+syi*guy8+szi*guz8)
                          - particleJ.dipole[2]*(sxi*gux9+syi*guy9+szi*guz9);

    dpwidy = particleI.charge*(sxk*gc6+syk*gc8+szk*gc9)
                        - particleJ.charge*(sxi*gux3+syi*guy3+szi*guz3)
                         - sxi*(particleJ.quadrupole[QXX]*gux12+particleJ.quadrupole[QYY]*gux17+particleJ.quadrupole[QZZ]*gux19
                        +2.0*(particleJ.quadrupole[QXY]*gux14+particleJ.quadrupole[QXZ]*gux15+particleJ.quadrupole[QYZ]*gux18))
                         - syi*(particleJ.quadrupole[QXX]*guy12+particleJ.quadrupole[QYY]*guy17+particleJ.quadrupole[QZZ]*guy19
                        +2.0*(particleJ.quadrupole[QXY]*guy14+particleJ.quadrupole[QXZ]*guy15+particleJ.quadrupole[QYZ]*guy18))
                         - szi*(particleJ.quadrupole[QXX]*guz12+particleJ.quadrupole[QYY]*guz17+particleJ.quadrupole[QZZ]*guz19
                        +2.0*(particleJ.quadrupole[QXY]*guz14+particleJ.quadrupole[QXZ]*guz15+particleJ.quadrupole[QYZ]*guz18))
                         + sxk*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                        +2.0*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
                         + syk*(particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
                        +2.0*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8))
                         + szk*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                        +2.0*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9));

    dpwkdy = particleI.charge*(sxk*gux3+syk*guy3+szk*guz3)
                        - particleJ.charge*(sxi*gc6+syi*gc8+szi*gc9)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                     +2.0*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
                      - syi*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
                     +2.0*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
                      - szi*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                     +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
                      + sxk*(particleI.quadrupole[QXX]*gux12+particleI.quadrupole[QYY]*gux17+particleI.quadrupole[QZZ]*gux19
                     +2.0*(particleI.quadrupole[QXY]*gux14+particleI.quadrupole[QXZ]*gux15+particleI.quadrupole[QYZ]*gux18))
                      + syk*(particleI.quadrupole[QXX]*guy12+particleI.quadrupole[QYY]*guy17+particleI.quadrupole[QZZ]*guy19
                     +2.0*(particleI.quadrupole[QXY]*guy14+particleI.quadrupole[QXZ]*guy15+particleI.quadrupole[QYZ]*guy18))
                      + szk*(particleI.quadrupole[QXX]*guz12+particleI.quadrupole[QYY]*guz17+particleI.quadrupole[QZZ]*guz19
                     +2.0*(particleI.quadrupole[QXY]*guz14+particleI.quadrupole[QXZ]*guz15+particleI.quadrupole[QYZ]*guz18));

    dpdy    = 0.5*(dpsymdy + 0.5*(dpwidy + dpwkdy));

    dpsymdz = -particleI.dipole[0]*(sxk*gux7+syk*guy7+szk*guz7)
                          - particleI.dipole[1]*(sxk*gux9+syk*guy9+szk*guz9)
                          - particleI.dipole[2]*(sxk*gux10+syk*guy10+szk*guz10)
                          - particleJ.dipole[0]*(sxi*gux7+syi*guy7+szi*guz7)
                          - particleJ.dipole[1]*(sxi*gux9+syi*guy9+szi*guz9)
                          - particleJ.dipole[2]*(sxi*gux10+syi*guy10+szi*guz10);

    dpwidz = particleI.charge*(sxk*gc7+syk*gc9+szk*gc10)
                        - particleJ.charge*(sxi*gux4+syi*guy4+szi*guz4)
                      - sxi*(particleJ.quadrupole[QXX]*gux13+particleJ.quadrupole[QYY]*gux18+particleJ.quadrupole[QZZ]*gux20
                     +2.0*(particleJ.quadrupole[QXY]*gux15+particleJ.quadrupole[QXZ]*gux16+particleJ.quadrupole[QYZ]*gux19))
                      - syi*(particleJ.quadrupole[QXX]*guy13+particleJ.quadrupole[QYY]*guy18+particleJ.quadrupole[QZZ]*guy20
                     +2.0*(particleJ.quadrupole[QXY]*guy15+particleJ.quadrupole[QXZ]*guy16+particleJ.quadrupole[QYZ]*guy19))
                      - szi*(particleJ.quadrupole[QXX]*guz13+particleJ.quadrupole[QYY]*guz18+particleJ.quadrupole[QZZ]*guz20
                     +2.0*(particleJ.quadrupole[QXY]*guz15+particleJ.quadrupole[QXZ]*guz16+particleJ.quadrupole[QYZ]*guz19))
                      + sxk*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                     +2.0*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
                      + syk*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                     +2.0*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
                      + szk*(particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
                     +2.0*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10));

    dpwkdz = particleI.charge*(sxk*gux4+syk*guy4+szk*guz4)
                        - particleJ.charge*(sxi*gc7+syi*gc9+szi*gc10)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                     +2.0*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
                      - syi*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                     +2.0*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
                      - szi*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
                     +2.0*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
                      + sxk*(particleI.quadrupole[QXX]*gux13+particleI.quadrupole[QYY]*gux18+particleI.quadrupole[QZZ]*gux20
                     +2.0*(particleI.quadrupole[QXY]*gux15+particleI.quadrupole[QXZ]*gux16+particleI.quadrupole[QYZ]*gux19))
                      + syk*(particleI.quadrupole[QXX]*guy13+particleI.quadrupole[QYY]*guy18+particleI.quadrupole[QZZ]*guy20
                     +2.0*(particleI.quadrupole[QXY]*guy15+particleI.quadrupole[QXZ]*guy16+particleI.quadrupole[QYZ]*guy19))
                      + szk*(particleI.quadrupole[QXX]*guz13+particleI.quadrupole[QYY]*guz18+particleI.quadrupole[QZZ]*guz20
                     +2.0*(particleI.quadrupole[QXY]*guz15+particleI.quadrupole[QXZ]*guz16+particleI.quadrupole[QYZ]*guz19));

    dpdz = 0.5*(dpsymdz + 0.5*(dpwidz + dpwkdz));

    // effective radii chain rule terms for the
    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    dsymdr = -particleI.dipole[0]*(sxk*gux22+syk*guy22+szk*guz22)
                          - particleI.dipole[1]*(sxk*gux23+syk*guy23+szk*guz23)
                          - particleI.dipole[2]*(sxk*gux24+syk*guy24+szk*guz24)
                          - particleJ.dipole[0]*(sxi*gux22+syi*guy22+szi*guz22)
                          - particleJ.dipole[1]*(sxi*gux23+syi*guy23+szi*guz23)
                          - particleJ.dipole[2]*(sxi*gux24+syi*guy24+szi*guz24);

    dwipdr = particleI.charge*(sxk*gc22+syk*gc23+szk*gc24)
                         - particleJ.charge*(sxi*gux21+syi*guy21+szi*guz21)
                      - sxi*(particleJ.quadrupole[QXX]*gux25+particleJ.quadrupole[QYY]*gux28+particleJ.quadrupole[QZZ]*gux30
                     +2.0*(particleJ.quadrupole[QXY]*gux26+particleJ.quadrupole[QXZ]*gux27+particleJ.quadrupole[QYZ]*gux29))
                      - syi*(particleJ.quadrupole[QXX]*guy25+particleJ.quadrupole[QYY]*guy28+particleJ.quadrupole[QZZ]*guy30
                     +2.0*(particleJ.quadrupole[QXY]*guy26+particleJ.quadrupole[QXZ]*guy27+particleJ.quadrupole[QYZ]*guy29))
                      - szi*(particleJ.quadrupole[QXX]*guz25+particleJ.quadrupole[QYY]*guz28+particleJ.quadrupole[QZZ]*guz30
                     +2.0*(particleJ.quadrupole[QXY]*guz26+particleJ.quadrupole[QXZ]*guz27+particleJ.quadrupole[QYZ]*guz29))
                      + sxk*(particleI.quadrupole[QXX]*gqxx22+particleI.quadrupole[QYY]*gqyy22+particleI.quadrupole[QZZ]*gqzz22
                     +2.0*(particleI.quadrupole[QXY]*gqxy22+particleI.quadrupole[QXZ]*gqxz22+particleI.quadrupole[QYZ]*gqyz22))
                      + syk*(particleI.quadrupole[QXX]*gqxx23+particleI.quadrupole[QYY]*gqyy23+particleI.quadrupole[QZZ]*gqzz23
                     +2.0*(particleI.quadrupole[QXY]*gqxy23+particleI.quadrupole[QXZ]*gqxz23+particleI.quadrupole[QYZ]*gqyz23))
                      + szk*(particleI.quadrupole[QXX]*gqxx24+particleI.quadrupole[QYY]*gqyy24+particleI.quadrupole[QZZ]*gqzz24
                     +2.0*(particleI.quadrupole[QXY]*gqxy24+particleI.quadrupole[QXZ]*gqxz24+particleI.quadrupole[QYZ]*gqyz24));

    dwkpdr = particleI.charge*(sxk*gux21+syk*guy21+szk*guz21)
                         - particleJ.charge*(sxi*gc22+syi*gc23+szi*gc24)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx22+particleJ.quadrupole[QYY]*gqyy22+particleJ.quadrupole[QZZ]*gqzz22
                     +2.0*(particleJ.quadrupole[QXY]*gqxy22+particleJ.quadrupole[QXZ]*gqxz22+particleJ.quadrupole[QYZ]*gqyz22))
                      - syi*(particleJ.quadrupole[QXX]*gqxx23+particleJ.quadrupole[QYY]*gqyy23+particleJ.quadrupole[QZZ]*gqzz23
                     +2.0*(particleJ.quadrupole[QXY]*gqxy23+particleJ.quadrupole[QXZ]*gqxz23+particleJ.quadrupole[QYZ]*gqyz23))
                      - szi*(particleJ.quadrupole[QXX]*gqxx24+particleJ.quadrupole[QYY]*gqyy24+particleJ.quadrupole[QZZ]*gqzz24
                     +2.0*(particleJ.quadrupole[QXY]*gqxy24+particleJ.quadrupole[QXZ]*gqxz24+particleJ.quadrupole[QYZ]*gqyz24))
                      + sxk*(particleI.quadrupole[QXX]*gux25+particleI.quadrupole[QYY]*gux28+particleI.quadrupole[QZZ]*gux30
                     +2.0*(particleI.quadrupole[QXY]*gux26+particleI.quadrupole[QXZ]*gux27+particleI.quadrupole[QYZ]*gux29))
                      + syk*(particleI.quadrupole[QXX]*guy25+particleI.quadrupole[QYY]*guy28+particleI.quadrupole[QZZ]*guy30
                     +2.0*(particleI.quadrupole[QXY]*guy26+particleI.quadrupole[QXZ]*guy27+particleI.quadrupole[QYZ]*guy29))
                      + szk*(particleI.quadrupole[QXX]*guz25+particleI.quadrupole[QYY]*guz28+particleI.quadrupole[QZZ]*guz30
                     +2.0*(particleI.quadrupole[QXY]*guz26+particleI.quadrupole[QXZ]*guz27+particleI.quadrupole[QYZ]*guz29));

    dsumdr = dsymdr + 0.5*(dwipdr + dwkpdr);
    dpbi = 0.5*_bornRadii[jIndex]*dsumdr;
    dpbk = 0.5*_bornRadii[iIndex]*dsumdr;

    // mutual polarization electrostatic solvation free energy gradient

    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Mutual || getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {

        dpdx = dpdx - 0.5 *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux5+_inducedDipolePolarS[jIndex][1]*gux6+_inducedDipolePolarS[jIndex][2]*gux7)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy5+_inducedDipolePolarS[jIndex][1]*guy6+_inducedDipolePolarS[jIndex][2]*guy7)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz5+_inducedDipolePolarS[jIndex][1]*guz6+_inducedDipolePolarS[jIndex][2]*guz7)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux5+_inducedDipolePolarS[iIndex][1]*gux6+_inducedDipolePolarS[iIndex][2]*gux7)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy5+_inducedDipolePolarS[iIndex][1]*guy6+_inducedDipolePolarS[iIndex][2]*guy7)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz5+_inducedDipolePolarS[iIndex][1]*guz6+_inducedDipolePolarS[iIndex][2]*guz7));

        dpdy = dpdy - 0.5 *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux6+_inducedDipolePolarS[jIndex][1]*gux8+_inducedDipolePolarS[jIndex][2]*gux9)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy6+_inducedDipolePolarS[jIndex][1]*guy8+_inducedDipolePolarS[jIndex][2]*guy9)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz6+_inducedDipolePolarS[jIndex][1]*guz8+_inducedDipolePolarS[jIndex][2]*guz9)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux6+_inducedDipolePolarS[iIndex][1]*gux8+_inducedDipolePolarS[iIndex][2]*gux9)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy6+_inducedDipolePolarS[iIndex][1]*guy8+_inducedDipolePolarS[iIndex][2]*guy9)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz6+_inducedDipolePolarS[iIndex][1]*guz8+_inducedDipolePolarS[iIndex][2]*guz9));

        dpdz = dpdz - 0.5 *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux7+_inducedDipolePolarS[jIndex][1]*gux9+_inducedDipolePolarS[jIndex][2]*gux10)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy7+_inducedDipolePolarS[jIndex][1]*guy9+_inducedDipolePolarS[jIndex][2]*guy10)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz7+_inducedDipolePolarS[jIndex][1]*guz9+_inducedDipolePolarS[jIndex][2]*guz10)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux7+_inducedDipolePolarS[iIndex][1]*gux9+_inducedDipolePolarS[iIndex][2]*gux10)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy7+_inducedDipolePolarS[iIndex][1]*guy9+_inducedDipolePolarS[iIndex][2]*guy10)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz7+_inducedDipolePolarS[iIndex][1]*guz9+_inducedDipolePolarS[iIndex][2]*guz10));

        duvdr = _inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux22+_inducedDipolePolarS[jIndex][1]*gux23+_inducedDipolePolarS[jIndex][2]*gux24)
                            + _inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy22+_inducedDipolePolarS[jIndex][1]*guy23+_inducedDipolePolarS[jIndex][2]*guy24)
                            + _inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz22+_inducedDipolePolarS[jIndex][1]*guz23+_inducedDipolePolarS[jIndex][2]*guz24)
                            + _inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux22+_inducedDipolePolarS[iIndex][1]*gux23+_inducedDipolePolarS[iIndex][2]*gux24)
                            + _inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy22+_inducedDipolePolarS[iIndex][1]*guy23+_inducedDipolePolarS[iIndex][2]*guy24)
                            + _inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz22+_inducedDipolePolarS[iIndex][1]*guz23+_inducedDipolePolarS[iIndex][2]*guz24);

        dpbi = dpbi - 0.5*_bornRadii[jIndex]*duvdr;
        dpbk = dpbk - 0.5*_bornRadii[iIndex]*duvdr;
    }

    // torque due to induced reaction field on permanent dipoles

    RealOpenMM fid1 = 0.5*(sxk*gux2+syk*guy2+szk*guz2);
    RealOpenMM fid2 = 0.5*(sxk*gux3+syk*guy3+szk*guz3);
    RealOpenMM fid3 = 0.5*(sxk*gux4+syk*guy4+szk*guz4);
    RealOpenMM fkd1 = 0.5*(sxi*gux2+syi*guy2+szi*guz2);
    RealOpenMM fkd2 = 0.5*(sxi*gux3+syi*guy3+szi*guz3);
    RealOpenMM fkd3 = 0.5*(sxi*gux4+syi*guy4+szi*guz4);

    RealOpenMM trqi1   = particleI.dipole[1]*fid3 - particleI.dipole[2]*fid2;
    RealOpenMM trqi2   = particleI.dipole[2]*fid1 - particleI.dipole[0]*fid3;
    RealOpenMM trqi3   = particleI.dipole[0]*fid2 - particleI.dipole[1]*fid1;

    RealOpenMM trqi_k1 = particleJ.dipole[1]*fkd3 - particleJ.dipole[2]*fkd2;
    RealOpenMM trqi_k2 = particleJ.dipole[2]*fkd1 - particleJ.dipole[0]*fkd3;
    RealOpenMM trqi_k3 = particleJ.dipole[0]*fkd2 - particleJ.dipole[1]*fkd1;

    // torque due to induced reaction field gradient on quadrupoles;

    RealOpenMM fidg11 = -0.25 *
                              ((sxk*gqxx2+syk*gqxx3+szk*gqxx4)
                              + (sxk*gux5+syk*guy5+szk*guz5));

    RealOpenMM fidg12 = -0.25 *
                              ((sxk*gqxy2+syk*gqxy3+szk*gqxy4)
                              + (sxk*gux6+syk*guy6+szk*guz6));

    RealOpenMM fidg13 = -0.25 *
                              ((sxk*gqxz2+syk*gqxz3+szk*gqxz4)
                              + (sxk*gux7+syk*guy7+szk*guz7));

    RealOpenMM fidg22 = -0.25 *
                              ((sxk*gqyy2+syk*gqyy3+szk*gqyy4)
                              + (sxk*gux8+syk*guy8+szk*guz8));

    RealOpenMM fidg23 = -0.25 *
                              ((sxk*gqyz2+syk*gqyz3+szk*gqyz4)
                              + (sxk*gux9+syk*guy9+szk*guz9));

    RealOpenMM fidg33 = -0.25 *
                              ((sxk*gqzz2+syk*gqzz3+szk*gqzz4)
                              + (sxk*gux10+syk*guy10+szk*guz10));

    RealOpenMM fidg21 = fidg12;
    RealOpenMM fidg31 = fidg13;
    RealOpenMM fidg32 = fidg23;

    RealOpenMM fkdg11 = 0.25 *
                              ((sxi*gqxx2+syi*gqxx3+szi*gqxx4)
                              + (sxi*gux5+syi*guy5+szi*guz5));

    RealOpenMM fkdg12 = 0.25 *
                              ((sxi*gqxy2+syi*gqxy3+szi*gqxy4)
                              + (sxi*gux6+syi*guy6+szi*guz6));
    RealOpenMM fkdg13 = 0.25 *
                              ((sxi*gqxz2+syi*gqxz3+szi*gqxz4)
                              + (sxi*gux7+syi*guy7+szi*guz7));
    RealOpenMM fkdg22 = 0.25 *
                              ((sxi*gqyy2+syi*gqyy3+szi*gqyy4)
                              + (sxi*gux8+syi*guy8+szi*guz8));
    RealOpenMM fkdg23 = 0.25 *
                              ((sxi*gqyz2+syi*gqyz3+szi*gqyz4)
                              + (sxi*gux9+syi*guy9+szi*guz9));
    RealOpenMM fkdg33 = 0.25 *
                              ((sxi*gqzz2+syi*gqzz3+szi*gqzz4)
                              + (sxi*gux10+syi*guy10+szi*guz10));
    RealOpenMM fkdg21 = fkdg12;
    RealOpenMM fkdg31 = fkdg13;
    RealOpenMM fkdg32 = fkdg23;

    trqi1 += 2.0*(particleI.quadrupole[QXY]*fidg13+particleI.quadrupole[QYY]*fidg23+particleI.quadrupole[QYZ]*fidg33
                        -particleI.quadrupole[QXZ]*fidg12-particleI.quadrupole[QYZ]*fidg22-particleI.quadrupole[QZZ]*fidg32);
    trqi2 += 2.0*(particleI.quadrupole[QXZ]*fidg11+particleI.quadrupole[QYZ]*fidg21+particleI.quadrupole[QZZ]*fidg31
                        -particleI.quadrupole[QXX]*fidg13-particleI.quadrupole[QXY]*fidg23-particleI.quadrupole[QXZ]*fidg33);

    trqi3 += 2.0*(particleI.quadrupole[QXX]*fidg12+particleI.quadrupole[QXY]*fidg22+particleI.quadrupole[QXZ]*fidg32
                        -particleI.quadrupole[QXY]*fidg11-particleI.quadrupole[QYY]*fidg21-particleI.quadrupole[QYZ]*fidg31);

    trqi_k1 += 2.0*
                        (particleJ.quadrupole[QXY]*fkdg13+particleJ.quadrupole[QYY]*fkdg23+particleJ.quadrupole[QYZ]*fkdg33
                        -particleJ.quadrupole[QXZ]*fkdg12-particleJ.quadrupole[QYZ]*fkdg22-particleJ.quadrupole[QZZ]*fkdg32);

    trqi_k2 += 2.0*
                        (particleJ.quadrupole[QXZ]*fkdg11+particleJ.quadrupole[QYZ]*fkdg21+particleJ.quadrupole[QZZ]*fkdg31
                        -particleJ.quadrupole[QXX]*fkdg13-particleJ.quadrupole[QXY]*fkdg23-particleJ.quadrupole[QXZ]*fkdg33);

    trqi_k3 += 2.0*
                        (particleJ.quadrupole[QXX]*fkdg12+particleJ.quadrupole[QXY]*fkdg22+particleJ.quadrupole[QXZ]*fkdg32
                        -particleJ.quadrupole[QXY]*fkdg11-particleJ.quadrupole[QYY]*fkdg21-particleJ.quadrupole[QYZ]*fkdg31);

    // total permanent and induced energies for this interaction

    e                        = esym + 0.5*(ewi+ewk);
    ei                       = 0.5*(esymi + 0.5*(ewii+ewki));
    RealOpenMM energy        = e + ei;
    energy                   = iIndex == jIndex ? 0.5*energy : energy;

    forces[iIndex][0]       += (dedx + dpdx);
    forces[iIndex][1]       += (dedy + dpdy);
    forces[iIndex][2]       += (dedz + dpdz);

    forces[jIndex][0]       -= (dedx + dpdx);
    forces[jIndex][1]       -= (dedy + dpdy);
    forces[jIndex][2]       -= (dedz + dpdz);

    torques[iIndex][0]      += (trq1 + trqi1);
    torques[iIndex][1]      += (trq2 + trqi2);
    torques[iIndex][2]      += (trq3 + trqi3);

    dBorn[iIndex]           += drbi + dpbi;
    if (iIndex != jIndex) {

        torques[jIndex][0]  += (trq_k1 + trqi_k1);
        torques[jIndex][1]  += (trq_k2 + trqi_k2);
        torques[jIndex][2]  += (trq_k3 + trqi_k3);

        dBorn[jIndex]       += drbk + dpbk;
    }

    return (energy);
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateElectrostatic(const vector<MultipoleParticleData>& particleData,
                                                                                    vector<RealVec>& torques,
                                                                                    vector<RealVec>& forces)
{

    RealOpenMM energy = AmoebaReferenceMultipoleForce::calculateElectrostatic(particleData, torques, forces);

    vector<RealOpenMM> dBorn;
    initializeRealOpenMMVector(dBorn);

    // Kirkwood loop over particle pairs

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii; jj < particleData.size(); jj++) {
            energy += calculateKirkwoodPairIxn(particleData[ii], particleData[jj], forces, torques, dBorn);
        }
    }

    // cavity term

    if (getIncludeCavityTerm()) {
        energy += calculateCavityTermEnergyAndForces(dBorn);
    }

    // apply Born chain rule; skip diagonal terms since these make no contribution to forces

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii+1; jj < particleData.size(); jj++) {
            calculateGrycukChainRulePairIxn(particleData[ii], particleData[jj], dBorn, forces);
            calculateGrycukChainRulePairIxn(particleData[jj], particleData[ii], dBorn, forces);
        }
    }

    // correct vacuum to SCRF derivatives (ediff1 in TINKER)

    vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for (unsigned int kk = 0; kk < scaleFactors.size(); kk++) {
        scaleFactors[kk] = 1.0;
    }

    RealOpenMM eDiffEnergy = 0.0;
    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii+1; jj < particleData.size(); jj++) {

            if (jj <= _maxScaleIndex[ii]) {
                getMultipoleScaleFactors(ii, jj, scaleFactors);
            }

            eDiffEnergy += calculateKirkwoodEDiffPairIxn(particleData[ii], particleData[jj],
                                                         scaleFactors[P_SCALE], scaleFactors[D_SCALE], forces, torques);

            if (jj <= _maxScaleIndex[ii]) {
                for (unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++) {
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }
    energy += (_electric/_dielectric)*eDiffEnergy;

    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
        RealOpenMM prefac = (_electric/_dielectric);
        for (int i = 0; i < _numParticles; i++) {
            // Compute the µ(m) T µ(n) force contributions here
            for (int l = 0; l < _maxPTOrder-1; ++l) {
                for (int m = 0; m < _maxPTOrder-1-l; ++m) {
                    RealOpenMM p = _extPartCoefficients[l+m+1];
                    if(std::fabs(p) < 1e-6) continue;
                    forces[i][0] -= 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+0]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+4]);
                    forces[i][1] -= 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+1]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+5]);
                    forces[i][2] -= 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+4]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+5]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+2]);
                    forces[i][0] -= 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+0]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+4]);
                    forces[i][1] -= 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+1]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+5]);
                    forces[i][2] -= 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+4]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+5]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+2]);
                    forces[i][0] += 0.5*p*prefac*(_ptDipoleDS[l][i][0]*_ptDipoleFieldGradientPS[m][6*i+0]
                                                + _ptDipoleDS[l][i][1]*_ptDipoleFieldGradientPS[m][6*i+3]
                                                + _ptDipoleDS[l][i][2]*_ptDipoleFieldGradientPS[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipoleDS[l][i][0]*_ptDipoleFieldGradientPS[m][6*i+3]
                                                + _ptDipoleDS[l][i][1]*_ptDipoleFieldGradientPS[m][6*i+1]
                                                + _ptDipoleDS[l][i][2]*_ptDipoleFieldGradientPS[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipoleDS[l][i][0]*_ptDipoleFieldGradientPS[m][6*i+4]
                                                + _ptDipoleDS[l][i][1]*_ptDipoleFieldGradientPS[m][6*i+5]
                                                + _ptDipoleDS[l][i][2]*_ptDipoleFieldGradientPS[m][6*i+2]);
                    forces[i][0] += 0.5*p*prefac*(_ptDipolePS[l][i][0]*_ptDipoleFieldGradientDS[m][6*i+0]
                                                + _ptDipolePS[l][i][1]*_ptDipoleFieldGradientDS[m][6*i+3]
                                                + _ptDipolePS[l][i][2]*_ptDipoleFieldGradientDS[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipolePS[l][i][0]*_ptDipoleFieldGradientDS[m][6*i+3]
                                                + _ptDipolePS[l][i][1]*_ptDipoleFieldGradientDS[m][6*i+1]
                                                + _ptDipolePS[l][i][2]*_ptDipoleFieldGradientDS[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipolePS[l][i][0]*_ptDipoleFieldGradientDS[m][6*i+4]
                                                + _ptDipolePS[l][i][1]*_ptDipoleFieldGradientDS[m][6*i+5]
                                                + _ptDipolePS[l][i][2]*_ptDipoleFieldGradientDS[m][6*i+2]);
                }
            }
        }
    }

    return energy;
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateGrycukChainRulePairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                                                       const vector<RealOpenMM>& dBorn, vector<RealVec>& forces) const
{
    unsigned int iIndex                 = particleI.particleIndex;
    unsigned int jIndex                 = particleJ.particleIndex;

    RealOpenMM pi43                     = (4.0/3.0)*M_PI;

    RealOpenMM lik, uik;
    RealOpenMM lik4, uik4;
    RealOpenMM factor                   = -POW(M_PI, (1.0/3.0))*POW(6.0,(2.0/3.0))/9.0;
    RealOpenMM term                     = pi43/(_bornRadii[iIndex]*_bornRadii[iIndex]*_bornRadii[iIndex]);
          term                          = factor/POW(term, (4.0/3.0));

    RealVec deltaR                      = particleJ.position - particleI.position;

    RealOpenMM sk                       = _scaledRadii[jIndex];
    RealOpenMM sk2                      = sk*sk;
    RealOpenMM r2                       = deltaR.dot(deltaR);
    RealOpenMM r                        = SQRT(r2);
    RealOpenMM de                       = 0.0;

    if ((_atomicRadii[iIndex] + r) < sk) {
        RealOpenMM uik4;
        uik                = sk - r;
        uik4               = uik*uik;
        uik4               = uik4*uik4;
        de                 = -4.0*M_PI/uik4;
    }

    if ((_atomicRadii[iIndex] + r) < sk) {
        lik          = sk - r;
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += 0.25*M_PI*(sk2-4.0*sk*r+17.0*r2)/ (r2*lik4);
    } else if (r < (_atomicRadii[iIndex] +sk)) {
        lik          = _atomicRadii[iIndex];
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += 0.25*M_PI*(2.0*_atomicRadii[iIndex]*_atomicRadii[iIndex]-sk2-r2)/ (r2*lik4);
    } else {
        lik          = r - sk;
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += 0.25*M_PI*(sk2-4.0*sk*r+r2)/ (r2*lik4);
    }
    uik                = r + sk;
    uik4               = uik*uik;
    uik4               = uik4*uik4;

    de                -= 0.25*M_PI*(sk2+4.0*sk*r+r2)/ (r2*uik4);
    RealOpenMM dbr     = term*de/r;
          de           = dbr*dBorn[iIndex];

    forces[iIndex]    -= deltaR*de;
    forces[jIndex]    += deltaR*de;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateCavityTermEnergyAndForces(vector<RealOpenMM>& dBorn) const
{

    RealOpenMM energy       = 0.0;
    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        RealOpenMM r        = _atomicRadii[ii] + _dielectricOffset + _probeRadius;
        RealOpenMM ratio    = ((_atomicRadii[ii] + _dielectricOffset)/_bornRadii[ii]);
        RealOpenMM saTerm   = _surfaceAreaFactor*r*r*POW(ratio, 6.0);
        dBorn[ii]          += saTerm/_bornRadii[ii];
        energy             += saTerm;
    }
    return (energy/-6.0);
}

/*****************************************************************************

       ediff1 correct vacuum to SCRF derivatives

       calculates the energy and derivatives of polarizing
       the vacuum induced dipoles to their SCRF polarized values

******************************************************************************/

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateKirkwoodEDiffPairIxn(const MultipoleParticleData& particleI,
                                                                                            const MultipoleParticleData& particleJ,
                                                                                            RealOpenMM pscale, RealOpenMM dscale,
                                                                                            vector<RealVec>& forces,
                                                                                            vector<RealVec>& torques) const
{
    static const RealOpenMM uscale      = 1.0;

    RealOpenMM scale3,scale5;
    RealOpenMM scale7;
    RealOpenMM psc3,psc5,psc7;
    RealOpenMM dsc3,dsc5,dsc7;
    RealOpenMM scale3i,scale5i;
    RealOpenMM r,rr1,rr3;
    RealOpenMM rr5,rr7,rr9;
    RealOpenMM pgamma;

    unsigned int iIndex   = particleI.particleIndex;
    unsigned int jIndex   = particleJ.particleIndex;

    // deltaR

    RealVec deltaR        = particleJ.position - particleI.position;
    RealOpenMM r2         = deltaR.dot(deltaR);

    RealOpenMM xr         = deltaR[0];
    RealOpenMM yr         = deltaR[1];
    RealOpenMM zr         = deltaR[2];

    r                     = SQRT(r2);
    rr1                   = 1.0/r;
    rr3                   = rr1/r2;
    rr5                   = 3.0*rr3/r2;
    rr7                   = 5.0*rr5/r2;
    rr9                   = 7.0*rr7/r2;

    scale3                = 1.0;
    scale5                = 1.0;
    scale7                = 1.0;

    RealOpenMM ddsc3_1    = 0.0;
    RealOpenMM ddsc3_2    = 0.0;
    RealOpenMM ddsc3_3    = 0.0;

    RealOpenMM ddsc5_1    = 0.0;
    RealOpenMM ddsc5_2    = 0.0;
    RealOpenMM ddsc5_3    = 0.0;

    RealOpenMM ddsc7_1    = 0.0;
    RealOpenMM ddsc7_2    = 0.0;
    RealOpenMM ddsc7_3    = 0.0;

    // apply Thole polarization damping to scale factors

    RealOpenMM damp = particleI.dampingFactor*particleJ.dampingFactor;
    if (damp != 0.0) {
        pgamma          = particleJ.thole > particleI.thole ? particleI.thole : particleJ.thole;
        RealOpenMM ratio     = (r/damp);
        damp            = -pgamma*ratio*ratio*ratio;
        if (damp > -50.0) {
            RealOpenMM dampE  = EXP(damp);
            RealOpenMM damp2  = damp*damp;
            scale3       = 1.0 - dampE;
            scale5       = 1.0 - (1.0 - damp)*dampE;
            scale7       = 1.0 - (1.0 - damp + 0.6*damp2)*dampE;

            ddsc3_1     = -3.0*damp*EXP(damp)*xr/r2;
            ddsc3_2     = -3.0*damp*EXP(damp)*yr/r2;
            ddsc3_3     = -3.0*damp*EXP(damp)*zr/r2;

            ddsc5_1     = -damp*ddsc3_1;
            ddsc5_2     = -damp*ddsc3_2;
            ddsc5_3     = -damp*ddsc3_3;

            ddsc7_1     = (-0.2-0.6*damp)*ddsc5_1;
            ddsc7_2     = (-0.2-0.6*damp)*ddsc5_2;
            ddsc7_3     = (-0.2-0.6*damp)*ddsc5_3;
        }
    }

    scale3i             = scale3*uscale;
    scale5i             = scale5*uscale;

    dsc3                = scale3*dscale;
    dsc5                = scale5*dscale;
    dsc7                = scale7*dscale;

    psc3                = scale3*pscale;
    psc5                = scale5*pscale;
    psc7                = scale7*pscale;

    // construct auxiliary vectors for permanent terms

    RealOpenMM dixr1             = particleI.dipole[1]*zr - particleI.dipole[2]*yr;
    RealOpenMM dixr2             = particleI.dipole[2]*xr - particleI.dipole[0]*zr;
    RealOpenMM dixr3             = particleI.dipole[0]*yr - particleI.dipole[1]*xr;

    RealOpenMM dkxr1             = particleJ.dipole[1]*zr - particleJ.dipole[2]*yr;
    RealOpenMM dkxr2             = particleJ.dipole[2]*xr - particleJ.dipole[0]*zr;
    RealOpenMM dkxr3             = particleJ.dipole[0]*yr - particleJ.dipole[1]*xr;

    RealOpenMM qir1              = particleI.quadrupole[QXX]*xr + particleI.quadrupole[QXY]*yr + particleI.quadrupole[QXZ]*zr;
    RealOpenMM qir2              = particleI.quadrupole[QXY]*xr + particleI.quadrupole[QYY]*yr + particleI.quadrupole[QYZ]*zr;
    RealOpenMM qir3              = particleI.quadrupole[QXZ]*xr + particleI.quadrupole[QYZ]*yr + particleI.quadrupole[QZZ]*zr;

    RealOpenMM qkr1              = particleJ.quadrupole[QXX]*xr + particleJ.quadrupole[QXY]*yr + particleJ.quadrupole[QXZ]*zr;
    RealOpenMM qkr2              = particleJ.quadrupole[QXY]*xr + particleJ.quadrupole[QYY]*yr + particleJ.quadrupole[QYZ]*zr;
    RealOpenMM qkr3              = particleJ.quadrupole[QXZ]*xr + particleJ.quadrupole[QYZ]*yr + particleJ.quadrupole[QZZ]*zr;

    RealOpenMM rxqir1            = yr*qir3 - zr*qir2;
    RealOpenMM rxqir2            = zr*qir1 - xr*qir3;
    RealOpenMM rxqir3            = xr*qir2 - yr*qir1;

    RealOpenMM rxqkr1            = yr*qkr3 - zr*qkr2;
    RealOpenMM rxqkr2            = zr*qkr1 - xr*qkr3;
    RealOpenMM rxqkr3            = xr*qkr2 - yr*qkr1;

    // get intermediate variables for permanent energy terms

    RealOpenMM sc3               = particleI.dipole[0]*xr  + particleI.dipole[1]*yr  + particleI.dipole[2]*zr;
    RealOpenMM sc4               = particleJ.dipole[0]*xr  + particleJ.dipole[1]*yr  + particleJ.dipole[2]*zr;
    RealOpenMM sc5               = qir1*xr + qir2*yr + qir3*zr;
    RealOpenMM sc6               = qkr1*xr + qkr2*yr + qkr3*zr;

    // construct auxiliary vectors for induced terms

    RealOpenMM dixuk1            = particleI.dipole[1]*_inducedDipoleS[jIndex][2] - particleI.dipole[2]*_inducedDipoleS[jIndex][1];
    RealOpenMM dixuk2            = particleI.dipole[2]*_inducedDipoleS[jIndex][0] - particleI.dipole[0]*_inducedDipoleS[jIndex][2];
    RealOpenMM dixuk3            = particleI.dipole[0]*_inducedDipoleS[jIndex][1] - particleI.dipole[1]*_inducedDipoleS[jIndex][0];

    RealOpenMM dkxui1            = particleJ.dipole[1]*_inducedDipoleS[iIndex][2] - particleJ.dipole[2]*_inducedDipoleS[iIndex][1];
    RealOpenMM dkxui2            = particleJ.dipole[2]*_inducedDipoleS[iIndex][0] - particleJ.dipole[0]*_inducedDipoleS[iIndex][2];
    RealOpenMM dkxui3            = particleJ.dipole[0]*_inducedDipoleS[iIndex][1] - particleJ.dipole[1]*_inducedDipoleS[iIndex][0];

    RealOpenMM dixukp1           = particleI.dipole[1]*_inducedDipolePolarS[jIndex][2] - particleI.dipole[2]*_inducedDipolePolarS[jIndex][1];
    RealOpenMM dixukp2           = particleI.dipole[2]*_inducedDipolePolarS[jIndex][0] - particleI.dipole[0]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM dixukp3           = particleI.dipole[0]*_inducedDipolePolarS[jIndex][1] - particleI.dipole[1]*_inducedDipolePolarS[jIndex][0];

    RealOpenMM dkxuip1           = particleJ.dipole[1]*_inducedDipolePolarS[iIndex][2] - particleJ.dipole[2]*_inducedDipolePolarS[iIndex][1];
    RealOpenMM dkxuip2           = particleJ.dipole[2]*_inducedDipolePolarS[iIndex][0] - particleJ.dipole[0]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM dkxuip3           = particleJ.dipole[0]*_inducedDipolePolarS[iIndex][1] - particleJ.dipole[1]*_inducedDipolePolarS[iIndex][0];

    RealOpenMM qiuk1             = particleI.quadrupole[QXX]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipoleS[jIndex][2];
    RealOpenMM qiuk2             = particleI.quadrupole[QXY]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipoleS[jIndex][2];
    RealOpenMM qiuk3             = particleI.quadrupole[QXZ]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipoleS[jIndex][2];

    RealOpenMM qkui1             = particleJ.quadrupole[QXX]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipoleS[iIndex][2];
    RealOpenMM qkui2             = particleJ.quadrupole[QXY]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipoleS[iIndex][2];
    RealOpenMM qkui3             = particleJ.quadrupole[QXZ]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipoleS[iIndex][2];

    RealOpenMM qiukp1            = particleI.quadrupole[QXX]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM qiukp2            = particleI.quadrupole[QXY]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM qiukp3            = particleI.quadrupole[QXZ]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipolePolarS[jIndex][2];

    RealOpenMM qkuip1            = particleJ.quadrupole[QXX]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM qkuip2            = particleJ.quadrupole[QXY]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM qkuip3            = particleJ.quadrupole[QXZ]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipolePolarS[iIndex][2];

    RealOpenMM uixqkr1           = _inducedDipoleS[iIndex][1]*qkr3 - _inducedDipoleS[iIndex][2]*qkr2;
    RealOpenMM uixqkr2           = _inducedDipoleS[iIndex][2]*qkr1 - _inducedDipoleS[iIndex][0]*qkr3;
    RealOpenMM uixqkr3           = _inducedDipoleS[iIndex][0]*qkr2 - _inducedDipoleS[iIndex][1]*qkr1;

    RealOpenMM ukxqir1           = _inducedDipoleS[jIndex][1]*qir3 - _inducedDipoleS[jIndex][2]*qir2;
    RealOpenMM ukxqir2           = _inducedDipoleS[jIndex][2]*qir1 - _inducedDipoleS[jIndex][0]*qir3;
    RealOpenMM ukxqir3           = _inducedDipoleS[jIndex][0]*qir2 - _inducedDipoleS[jIndex][1]*qir1;

    RealOpenMM uixqkrp1          = _inducedDipolePolarS[iIndex][1]*qkr3 - _inducedDipolePolarS[iIndex][2]*qkr2;
    RealOpenMM uixqkrp2          = _inducedDipolePolarS[iIndex][2]*qkr1 - _inducedDipolePolarS[iIndex][0]*qkr3;
    RealOpenMM uixqkrp3          = _inducedDipolePolarS[iIndex][0]*qkr2 - _inducedDipolePolarS[iIndex][1]*qkr1;

    RealOpenMM ukxqirp1          = _inducedDipolePolarS[jIndex][1]*qir3 - _inducedDipolePolarS[jIndex][2]*qir2;
    RealOpenMM ukxqirp2          = _inducedDipolePolarS[jIndex][2]*qir1 - _inducedDipolePolarS[jIndex][0]*qir3;
    RealOpenMM ukxqirp3          = _inducedDipolePolarS[jIndex][0]*qir2 - _inducedDipolePolarS[jIndex][1]*qir1;

    RealOpenMM rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    RealOpenMM rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    RealOpenMM rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    RealOpenMM rxqkui1           = yr*qkui3 - zr*qkui2;
    RealOpenMM rxqkui2           = zr*qkui1 - xr*qkui3;
    RealOpenMM rxqkui3           = xr*qkui2 - yr*qkui1;

    RealOpenMM rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    RealOpenMM rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    RealOpenMM rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    RealOpenMM rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    RealOpenMM rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    RealOpenMM rxqkuip3          = xr*qkuip2 - yr*qkuip1;

    // get intermediate variables for induction energy terms

    RealOpenMM sci1              = _inducedDipoleS[iIndex][0]*particleJ.dipole[0] + _inducedDipoleS[iIndex][1]*particleJ.dipole[1] +
                              _inducedDipoleS[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipoleS[jIndex][0] +
                              particleI.dipole[1]*_inducedDipoleS[jIndex][1] + particleI.dipole[2]*_inducedDipoleS[jIndex][2];

    RealOpenMM sci3              = _inducedDipoleS[iIndex][0]*xr + _inducedDipoleS[iIndex][1]*yr + _inducedDipoleS[iIndex][2]*zr;
    RealOpenMM sci4              = _inducedDipoleS[jIndex][0]*xr + _inducedDipoleS[jIndex][1]*yr + _inducedDipoleS[jIndex][2]*zr;

    RealOpenMM sci7              = qir1*_inducedDipoleS[jIndex][0] + qir2*_inducedDipoleS[jIndex][1] + qir3*_inducedDipoleS[jIndex][2];
    RealOpenMM sci8              = qkr1*_inducedDipoleS[iIndex][0] + qkr2*_inducedDipoleS[iIndex][1] + qkr3*_inducedDipoleS[iIndex][2];

    RealOpenMM scip1             = _inducedDipolePolarS[iIndex][0]*particleJ.dipole[0] + _inducedDipolePolarS[iIndex][1]*particleJ.dipole[1] +
                              _inducedDipolePolarS[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipolePolarS[jIndex][0] +
                              particleI.dipole[1]*_inducedDipolePolarS[jIndex][1] + particleI.dipole[2]*_inducedDipolePolarS[jIndex][2];

    RealOpenMM scip2             = _inducedDipoleS[iIndex][0]*_inducedDipolePolarS[jIndex][0] + _inducedDipoleS[iIndex][1]*_inducedDipolePolarS[jIndex][1] +
                              _inducedDipoleS[iIndex][2]*_inducedDipolePolarS[jIndex][2] + _inducedDipolePolarS[iIndex][0]*_inducedDipoleS[jIndex][0] +
                              _inducedDipolePolarS[iIndex][1]*_inducedDipoleS[jIndex][1] + _inducedDipolePolarS[iIndex][2]*_inducedDipoleS[jIndex][2];

    RealOpenMM scip3             = _inducedDipolePolarS[iIndex][0]*xr + _inducedDipolePolarS[iIndex][1]*yr + _inducedDipolePolarS[iIndex][2]*zr;
    RealOpenMM scip4             = _inducedDipolePolarS[jIndex][0]*xr + _inducedDipolePolarS[jIndex][1]*yr + _inducedDipolePolarS[jIndex][2]*zr;

    RealOpenMM scip7             = qir1*_inducedDipolePolarS[jIndex][0] + qir2*_inducedDipolePolarS[jIndex][1] + qir3*_inducedDipolePolarS[jIndex][2];
    RealOpenMM scip8             = qkr1*_inducedDipolePolarS[iIndex][0] + qkr2*_inducedDipolePolarS[iIndex][1] + qkr3*_inducedDipolePolarS[iIndex][2];

    // calculate the gl functions for potential energy

    RealOpenMM gli1              = particleJ.charge*sci3 - particleI.charge*sci4;
    RealOpenMM gli2              = -sc3*sci4 - sci3*sc4;
    RealOpenMM gli3              = sci3*sc6 - sci4*sc5;
    RealOpenMM gli6              = sci1;
    RealOpenMM gli7              = 2.0*(sci7-sci8);
    RealOpenMM glip1             = particleJ.charge*scip3 - particleI.charge*scip4;
    RealOpenMM glip2             = -sc3*scip4 - scip3*sc4;
    RealOpenMM glip3             = scip3*sc6 - scip4*sc5;
    RealOpenMM glip6             = scip1;
    RealOpenMM glip7             = 2.0*(scip7-scip8);

    // get the permanent multipole and induced energies

    RealOpenMM energy            = 0.5*(rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms

    RealOpenMM gfi1              = 0.5*rr5*((gli1+gli6)*psc3 +
                                        (glip1+glip6)*dsc3+scip2*scale3i) +
                              0.5*rr7*((gli7+gli2)*psc5 +
                                        (glip7+glip2)*dsc5 -
                              (sci3*scip4+scip3*sci4)*scale5i) +
                              0.5*rr9*(gli3*psc7+glip3*dsc7);

    RealOpenMM gfi4              = 2.0*rr5;
    RealOpenMM gfi5              = rr7*(sci4*psc7+scip4*dsc7);
    RealOpenMM gfi6              = -rr7*(sci3*psc7+scip3*dsc7);

    // get the induced force

    RealOpenMM ftm2i1            = gfi1*xr + 0.5*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][0]*psc3+_inducedDipolePolarS[iIndex][0]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][0]*psc5+_inducedDipolePolarS[iIndex][0]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][0]*psc7+_inducedDipolePolarS[iIndex][0]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][0]*psc3+_inducedDipolePolarS[jIndex][0]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][0]*psc5+_inducedDipolePolarS[jIndex][0]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][0]*psc7+_inducedDipolePolarS[jIndex][0]*dsc7))*0.5
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][0]+scip4*_inducedDipoleS[iIndex][0]
                           + sci3*_inducedDipolePolarS[jIndex][0]+scip3*_inducedDipoleS[jIndex][0])*0.5
                           + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[0]
                           + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[0]
                           + 0.5*gfi4*((qkui1-qiuk1)*psc5
                           + (qkuip1-qiukp1)*dsc5)
                           + gfi5*qir1 + gfi6*qkr1;

    RealOpenMM ftm2i2            = gfi1*yr + 0.5*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][1]*psc3+_inducedDipolePolarS[iIndex][1]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][1]*psc5+_inducedDipolePolarS[iIndex][1]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][1]*psc7+_inducedDipolePolarS[iIndex][1]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][1]*psc3+_inducedDipolePolarS[jIndex][1]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][1]*psc5+_inducedDipolePolarS[jIndex][1]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][1]*psc7+_inducedDipolePolarS[jIndex][1]*dsc7))*0.5
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][1]+scip4*_inducedDipoleS[iIndex][1]
                           + sci3*_inducedDipolePolarS[jIndex][1]+scip3*_inducedDipoleS[jIndex][1])*0.5
                           + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[1]
                           + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[1]
                           + 0.5*gfi4*((qkui2-qiuk2)*psc5
                           + (qkuip2-qiukp2)*dsc5)
                           + gfi5*qir2 + gfi6*qkr2;

    RealOpenMM ftm2i3            = gfi1*zr  + 0.5*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][2]*psc3+_inducedDipolePolarS[iIndex][2]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][2]*psc5+_inducedDipolePolarS[iIndex][2]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][2]*psc7+_inducedDipolePolarS[iIndex][2]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][2]*psc3+_inducedDipolePolarS[jIndex][2]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][2]*psc5+_inducedDipolePolarS[jIndex][2]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][2]*psc7+_inducedDipolePolarS[jIndex][2]*dsc7))*0.5
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][2]+scip4*_inducedDipoleS[iIndex][2]
                           + sci3*_inducedDipolePolarS[jIndex][2]+scip3*_inducedDipoleS[jIndex][2])*0.5
                           + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[2]
                           + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[2]
                           + 0.5*gfi4*((qkui3-qiuk3)*psc5
                           + (qkuip3-qiukp3)*dsc5)
                           + gfi5*qir3 + gfi6*qkr3;

    // intermediate values needed for partially excluded interactions

    RealOpenMM fridmp1           = 0.5*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_1
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_1
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    RealOpenMM fridmp2           = 0.5*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_2
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_2
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    RealOpenMM fridmp3           = 0.5*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_3
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_3
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms

    RealOpenMM findmp1           = 0.5*uscale*(scip2*rr3*ddsc3_1
                          - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    RealOpenMM findmp2           = 0.5*uscale*(scip2*rr3*ddsc3_2
                          - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    RealOpenMM findmp3           = 0.5*uscale*(scip2*rr3*ddsc3_3
                          - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions

    ftm2i1           -= fridmp1 + findmp1;
    ftm2i2           -= fridmp2 + findmp2;
    ftm2i3           -= fridmp3 + findmp3;

    // correction to convert mutual to direct polarization force

    if (getPolarizationType() != AmoebaReferenceMultipoleForce::Mutual) {
        RealOpenMM gfd      = 0.5*(rr5*scip2*scale3i - rr7*(scip3*sci4+sci3*scip4)*scale5i);
        RealOpenMM fdir1    = gfd*xr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][0]+scip4*_inducedDipoleS[iIndex][0] + sci3*_inducedDipolePolarS[jIndex][0]+scip3*_inducedDipoleS[jIndex][0]);
        RealOpenMM fdir2    = gfd*yr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][1]+scip4*_inducedDipoleS[iIndex][1] + sci3*_inducedDipolePolarS[jIndex][1]+scip3*_inducedDipoleS[jIndex][1]);
        RealOpenMM fdir3    = gfd*zr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][2]+scip4*_inducedDipoleS[iIndex][2] + sci3*_inducedDipolePolarS[jIndex][2]+scip3*_inducedDipoleS[jIndex][2]);
        ftm2i1        -= fdir1 - findmp1;
        ftm2i2        -= fdir2 - findmp2;
        ftm2i3        -= fdir3 - findmp3;

    }

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and j

    RealOpenMM gti2              = 0.5*(sci4*psc5+scip4*dsc5)*rr5;
    RealOpenMM gti3              = 0.5*(sci3*psc5+scip3*dsc5)*rr5;
    RealOpenMM gti4              = gfi4;
    RealOpenMM gti5              = gfi5;
    RealOpenMM gti6              = gfi6;

    // calculate the induced torque components

    RealOpenMM ttm2i1            = -rr3*(dixuk1*psc3+dixukp1*dsc3)*0.5
                          + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                          +(ukxqirp1+rxqiukp1)*dsc5)*0.5 - gti5*rxqir1;

    RealOpenMM ttm2i2            = -rr3*(dixuk2*psc3+dixukp2*dsc3)*0.5
                          + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                          +(ukxqirp2+rxqiukp2)*dsc5)*0.5 - gti5*rxqir2;

    RealOpenMM ttm2i3            = -rr3*(dixuk3*psc3+dixukp3*dsc3)*0.5
                          + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                          +(ukxqirp3+rxqiukp3)*dsc5)*0.5 - gti5*rxqir3;

    RealOpenMM ttm3i1            = -rr3*(dkxui1*psc3+dkxuip1*dsc3)*0.5
                          + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                          +(uixqkrp1+rxqkuip1)*dsc5)*0.5 - gti6*rxqkr1;

    RealOpenMM ttm3i2            = -rr3*(dkxui2*psc3+dkxuip2*dsc3)*0.5
                          + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                          +(uixqkrp2+rxqkuip2)*dsc5)*0.5 - gti6*rxqkr2;

    RealOpenMM ttm3i3            = -rr3*(dkxui3*psc3+dkxuip3*dsc3)*0.5
                          + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                          +(uixqkrp3+rxqkuip3)*dsc5)*0.5 - gti6*rxqkr3;

    // update force and torque

    RealVec force, torqueI, torqueJ;
    force[0]      = -ftm2i1;
    force[1]      = -ftm2i2;
    force[2]      = -ftm2i3;

    torqueI[0]    = ttm2i1;
    torqueI[1]    = ttm2i2;
    torqueI[2]    = ttm2i3;

    torqueJ[0]    = ttm3i1;
    torqueJ[1]    = ttm3i2;
    torqueJ[2]    = ttm3i3;

    // construct auxiliary vectors for induced terms

    dixuk1            = particleI.dipole[1]*_inducedDipole[jIndex][2] - particleI.dipole[2]*_inducedDipole[jIndex][1];
    dixuk2            = particleI.dipole[2]*_inducedDipole[jIndex][0] - particleI.dipole[0]*_inducedDipole[jIndex][2];
    dixuk3            = particleI.dipole[0]*_inducedDipole[jIndex][1] - particleI.dipole[1]*_inducedDipole[jIndex][0];

    dkxui1            = particleJ.dipole[1]*_inducedDipole[iIndex][2] - particleJ.dipole[2]*_inducedDipole[iIndex][1];
    dkxui2            = particleJ.dipole[2]*_inducedDipole[iIndex][0] - particleJ.dipole[0]*_inducedDipole[iIndex][2];
    dkxui3            = particleJ.dipole[0]*_inducedDipole[iIndex][1] - particleJ.dipole[1]*_inducedDipole[iIndex][0];

    dixukp1           = particleI.dipole[1]*_inducedDipolePolar[jIndex][2] - particleI.dipole[2]*_inducedDipolePolar[jIndex][1];
    dixukp2           = particleI.dipole[2]*_inducedDipolePolar[jIndex][0] - particleI.dipole[0]*_inducedDipolePolar[jIndex][2];
    dixukp3           = particleI.dipole[0]*_inducedDipolePolar[jIndex][1] - particleI.dipole[1]*_inducedDipolePolar[jIndex][0];

    dkxuip1           = particleJ.dipole[1]*_inducedDipolePolar[iIndex][2] - particleJ.dipole[2]*_inducedDipolePolar[iIndex][1];
    dkxuip2           = particleJ.dipole[2]*_inducedDipolePolar[iIndex][0] - particleJ.dipole[0]*_inducedDipolePolar[iIndex][2];
    dkxuip3           = particleJ.dipole[0]*_inducedDipolePolar[iIndex][1] - particleJ.dipole[1]*_inducedDipolePolar[iIndex][0];

    qiuk1             = particleI.quadrupole[QXX]*_inducedDipole[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipole[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipole[jIndex][2];
    qiuk2             = particleI.quadrupole[QXY]*_inducedDipole[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipole[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipole[jIndex][2];
    qiuk3             = particleI.quadrupole[QXZ]*_inducedDipole[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipole[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipole[jIndex][2];

    qkui1             = particleJ.quadrupole[QXX]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipole[iIndex][2];
    qkui2             = particleJ.quadrupole[QXY]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipole[iIndex][2];
    qkui3             = particleJ.quadrupole[QXZ]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipole[iIndex][2];

    qiukp1            = particleI.quadrupole[QXX]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipolePolar[jIndex][2];
    qiukp2            = particleI.quadrupole[QXY]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipolePolar[jIndex][2];
    qiukp3            = particleI.quadrupole[QXZ]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipolePolar[jIndex][2];

    qkuip1            = particleJ.quadrupole[QXX]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipolePolar[iIndex][2];
    qkuip2            = particleJ.quadrupole[QXY]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipolePolar[iIndex][2];
    qkuip3            = particleJ.quadrupole[QXZ]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipolePolar[iIndex][2];

    uixqkr1           = _inducedDipole[iIndex][1]*qkr3 - _inducedDipole[iIndex][2]*qkr2;
    uixqkr2           = _inducedDipole[iIndex][2]*qkr1 - _inducedDipole[iIndex][0]*qkr3;
    uixqkr3           = _inducedDipole[iIndex][0]*qkr2 - _inducedDipole[iIndex][1]*qkr1;

    ukxqir1           = _inducedDipole[jIndex][1]*qir3 - _inducedDipole[jIndex][2]*qir2;
    ukxqir2           = _inducedDipole[jIndex][2]*qir1 - _inducedDipole[jIndex][0]*qir3;
    ukxqir3           = _inducedDipole[jIndex][0]*qir2 - _inducedDipole[jIndex][1]*qir1;

    uixqkrp1          = _inducedDipolePolar[iIndex][1]*qkr3 - _inducedDipolePolar[iIndex][2]*qkr2;
    uixqkrp2          = _inducedDipolePolar[iIndex][2]*qkr1 - _inducedDipolePolar[iIndex][0]*qkr3;
    uixqkrp3          = _inducedDipolePolar[iIndex][0]*qkr2 - _inducedDipolePolar[iIndex][1]*qkr1;

    ukxqirp1          = _inducedDipolePolar[jIndex][1]*qir3 - _inducedDipolePolar[jIndex][2]*qir2;
    ukxqirp2          = _inducedDipolePolar[jIndex][2]*qir1 - _inducedDipolePolar[jIndex][0]*qir3;
    ukxqirp3          = _inducedDipolePolar[jIndex][0]*qir2 - _inducedDipolePolar[jIndex][1]*qir1;

    rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    rxqkui1           = yr*qkui3 - zr*qkui2;
    rxqkui2           = zr*qkui1 - xr*qkui3;
    rxqkui3           = xr*qkui2 - yr*qkui1;

    rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    rxqkuip3          = xr*qkuip2 - yr*qkuip1;

    // get intermediate variables for induction energy terms

    sci1              = _inducedDipole[iIndex][0]*particleJ.dipole[0] + _inducedDipole[iIndex][1]*particleJ.dipole[1]
                          + _inducedDipole[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipole[jIndex][0]
                          + particleI.dipole[1]*_inducedDipole[jIndex][1] + particleI.dipole[2]*_inducedDipole[jIndex][2];

    sci3              = _inducedDipole[iIndex][0]*xr + _inducedDipole[iIndex][1]*yr + _inducedDipole[iIndex][2]*zr;
    sci4              = _inducedDipole[jIndex][0]*xr + _inducedDipole[jIndex][1]*yr + _inducedDipole[jIndex][2]*zr;

    sci7              = qir1*_inducedDipole[jIndex][0] + qir2*_inducedDipole[jIndex][1] + qir3*_inducedDipole[jIndex][2];
    sci8              = qkr1*_inducedDipole[iIndex][0] + qkr2*_inducedDipole[iIndex][1] + qkr3*_inducedDipole[iIndex][2];

    scip1             = _inducedDipolePolar[iIndex][0]*particleJ.dipole[0] + _inducedDipolePolar[iIndex][1]*particleJ.dipole[1] + _inducedDipolePolar[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipolePolar[jIndex][0] + particleI.dipole[1]*_inducedDipolePolar[jIndex][1] + particleI.dipole[2]*_inducedDipolePolar[jIndex][2];
    scip2             = _inducedDipole[iIndex][0]*_inducedDipolePolar[jIndex][0]+_inducedDipole[iIndex][1]*_inducedDipolePolar[jIndex][1] + _inducedDipole[iIndex][2]*_inducedDipolePolar[jIndex][2]+_inducedDipolePolar[iIndex][0]*_inducedDipole[jIndex][0] + _inducedDipolePolar[iIndex][1]*_inducedDipole[jIndex][1]+_inducedDipolePolar[iIndex][2]*_inducedDipole[jIndex][2];

    scip3             = _inducedDipolePolar[iIndex][0]*xr + _inducedDipolePolar[iIndex][1]*yr + _inducedDipolePolar[iIndex][2]*zr;
    scip4             = _inducedDipolePolar[jIndex][0]*xr + _inducedDipolePolar[jIndex][1]*yr + _inducedDipolePolar[jIndex][2]*zr;

    scip7             = qir1*_inducedDipolePolar[jIndex][0] + qir2*_inducedDipolePolar[jIndex][1] + qir3*_inducedDipolePolar[jIndex][2];
    scip8             = qkr1*_inducedDipolePolar[iIndex][0] + qkr2*_inducedDipolePolar[iIndex][1] + qkr3*_inducedDipolePolar[iIndex][2];

    // calculate the gl functions for potential energy

    gli1              = particleJ.charge*sci3 - particleI.charge*sci4;
    gli2              = -sc3*sci4 - sci3*sc4;
    gli3              = sci3*sc6 - sci4*sc5;
    gli6              = sci1;
    gli7              = 2.0*(sci7-sci8);

    glip1             = particleJ.charge*scip3 - particleI.charge*scip4;
    glip2             = -sc3*scip4 - scip3*sc4;
    glip3             = scip3*sc6 - scip4*sc5;
    glip6             = scip1;
    glip7             = 2.0*(scip7-scip8);

    // get the permanent multipole and induced energies

    energy           += -0.5*(rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms

    gfi1             = 0.5*rr5*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3+scip2*scale3i)
                         + 0.5*rr7*((gli7+gli2)*psc5 +(glip7+glip2)*dsc5
                         -(sci3*scip4+scip3*sci4)*scale5i)
                         + 0.5*rr9*(gli3*psc7+glip3*dsc7);

    gfi4             = 2.0*rr5;
    gfi5             = rr7*(sci4*psc7+scip4*dsc7);
    gfi6             = -rr7*(sci3*psc7+scip3*dsc7);

    // get the induced force

    ftm2i1           = gfi1*xr + 0.5*
                        (- rr3*particleJ.charge*(_inducedDipole[iIndex][0]*psc3+_inducedDipolePolar[iIndex][0]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][0]*psc5+_inducedDipolePolar[iIndex][0]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][0]*psc7+_inducedDipolePolar[iIndex][0]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][0]*psc3+_inducedDipolePolar[jIndex][0]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][0]*psc5+_inducedDipolePolar[jIndex][0]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][0]*psc7+_inducedDipolePolar[jIndex][0]*dsc7))*0.5
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0]
                         + sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0])*0.5
                         + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[0]
                         + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[0]
                         + 0.5*gfi4*((qkui1-qiuk1)*psc5
                         + (qkuip1-qiukp1)*dsc5)
                         + gfi5*qir1 + gfi6*qkr1;

    ftm2i2           = gfi1*yr + 0.5*
                        (- rr3*particleJ.charge*(_inducedDipole[iIndex][1]*psc3+_inducedDipolePolar[iIndex][1]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][1]*psc5+_inducedDipolePolar[iIndex][1]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][1]*psc7+_inducedDipolePolar[iIndex][1]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][1]*psc3+_inducedDipolePolar[jIndex][1]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][1]*psc5+_inducedDipolePolar[jIndex][1]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][1]*psc7+_inducedDipolePolar[jIndex][1]*dsc7))*0.5
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1]
                         + sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1])*0.5
                         + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[1]
                         + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[1]
                         + 0.5*gfi4*((qkui2-qiuk2)*psc5
                         + (qkuip2-qiukp2)*dsc5)
                         + gfi5*qir2 + gfi6*qkr2;

    ftm2i3           = gfi1*zr  + 0.5*
                         (- rr3*particleJ.charge*(_inducedDipole[iIndex][2]*psc3+_inducedDipolePolar[iIndex][2]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][2]*psc5+_inducedDipolePolar[iIndex][2]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][2]*psc7+_inducedDipolePolar[iIndex][2]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][2]*psc3+_inducedDipolePolar[jIndex][2]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][2]*psc5+_inducedDipolePolar[jIndex][2]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][2]*psc7+_inducedDipolePolar[jIndex][2]*dsc7))*0.5
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2]
                         + sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2])*0.5
                         + 0.5*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[2]
                         + 0.5*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[2]
                         + 0.5*gfi4*((qkui3-qiuk3)*psc5
                         + (qkuip3-qiukp3)*dsc5)
                         + gfi5*qir3 + gfi6*qkr3;

    // intermediate values needed for partially excluded interactions

    fridmp1          = 0.5*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_1
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_1
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    fridmp2          = 0.5*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_2
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_2
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    fridmp3          = 0.5*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_3
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_3
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms

    findmp1          = 0.5*uscale*(scip2*rr3*ddsc3_1
                         - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    findmp2          = 0.5*uscale*(scip2*rr3*ddsc3_2
                         - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    findmp3          = 0.5*uscale*(scip2*rr3*ddsc3_3
                         - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions

    ftm2i1           = ftm2i1 - fridmp1 - findmp1;
    ftm2i2           = ftm2i2 - fridmp2 - findmp2;
    ftm2i3           = ftm2i3 - fridmp3 - findmp3;

    // correction to convert mutual to direct polarization force

    if (getPolarizationType() != AmoebaReferenceMultipoleForce::Mutual) {

        RealOpenMM gfd    = 0.5*(rr5*scip2*scale3i- rr7*(scip3*sci4+sci3*scip4)*scale5i);
        RealOpenMM fdir1  = gfd*xr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0] + sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0]);
        RealOpenMM fdir2  = gfd*yr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1] + sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1]);
        RealOpenMM fdir3  = gfd*zr + 0.5*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2] + sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2]);
        ftm2i1      -= fdir1 - findmp1;
        ftm2i2      -= fdir2 - findmp2;
        ftm2i3      -= fdir3 - findmp3;
    }

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and j

    gti2             = 0.5*(sci4*psc5+scip4*dsc5)*rr5;
    gti3             = 0.5*(sci3*psc5+scip3*dsc5)*rr5;
    gti4             = gfi4;
    gti5             = gfi5;
    gti6             = gfi6;

    // calculate the induced torque components

    torqueI[0]      -= -rr3*(dixuk1*psc3+dixukp1*dsc3)*0.5
                         + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                         +(ukxqirp1+rxqiukp1)*dsc5)*0.5 - gti5*rxqir1;

    torqueI[1]      -= -rr3*(dixuk2*psc3+dixukp2*dsc3)*0.5
                         + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                         +(ukxqirp2+rxqiukp2)*dsc5)*0.5 - gti5*rxqir2;

    torqueI[2]      -= -rr3*(dixuk3*psc3+dixukp3*dsc3)*0.5
                         + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                         +(ukxqirp3+rxqiukp3)*dsc5)*0.5 - gti5*rxqir3;

    torqueJ[0]      -= -rr3*(dkxui1*psc3+dkxuip1*dsc3)*0.5
                         + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                         +(uixqkrp1+rxqkuip1)*dsc5)*0.5 - gti6*rxqkr1;

    torqueJ[1]      -= -rr3*(dkxui2*psc3+dkxuip2*dsc3)*0.5
                         + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                         +(uixqkrp2+rxqkuip2)*dsc5)*0.5 - gti6*rxqkr2;

    torqueJ[2]      -= -rr3*(dkxui3*psc3+dkxuip3*dsc3)*0.5
                         + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                         +(uixqkrp3+rxqkuip3)*dsc5)*0.5 - gti6*rxqkr3;

    // update force and torque

    force[0]       += ftm2i1;
    force[1]       += ftm2i2;
    force[2]       += ftm2i3;

    force          *= (_electric/_dielectric);
    forces[iIndex] += force;
    forces[jIndex] -= force;

    torqueI         *= (_electric/_dielectric);
    torqueJ         *= (_electric/_dielectric);
    torques[iIndex] += torqueI;
    torques[jIndex] += torqueJ;

    return energy;
}

const int AmoebaReferencePmeMultipoleForce::AMOEBA_PME_ORDER = 5;

const RealOpenMM AmoebaReferencePmeMultipoleForce::SQRT_PI = 1.77245385091;

AmoebaReferencePmeMultipoleForce::AmoebaReferencePmeMultipoleForce() :
               AmoebaReferenceMultipoleForce(PME),
               _cutoffDistance(1.0), _cutoffDistanceSquared(1.0),
               _pmeGridSize(0), _totalGridSize(0), _alphaEwald(0.0)
{

    _fftplan = NULL;
    _pmeGrid = NULL;
    _pmeGridDimensions = IntVec(-1, -1, -1);
}

AmoebaReferencePmeMultipoleForce::~AmoebaReferencePmeMultipoleForce()
{
    if (_fftplan) {
        fftpack_destroy(_fftplan);
    }
    if (_pmeGrid) {
        delete _pmeGrid;
    }
};

RealOpenMM AmoebaReferencePmeMultipoleForce::getCutoffDistance() const
{
     return _cutoffDistance;
};

void AmoebaReferencePmeMultipoleForce::setCutoffDistance(RealOpenMM cutoffDistance)
{
     _cutoffDistance        = cutoffDistance;
     _cutoffDistanceSquared = cutoffDistance*cutoffDistance;
};

RealOpenMM AmoebaReferencePmeMultipoleForce::getAlphaEwald() const
{
     return _alphaEwald;
};

void AmoebaReferencePmeMultipoleForce::setAlphaEwald(RealOpenMM alphaEwald)
{
     _alphaEwald = alphaEwald;
};

void AmoebaReferencePmeMultipoleForce::getPmeGridDimensions(vector<int>& pmeGridDimensions) const
{

    pmeGridDimensions.resize(3);

    pmeGridDimensions[0] = _pmeGridDimensions[0];
    pmeGridDimensions[1] = _pmeGridDimensions[1];
    pmeGridDimensions[2] = _pmeGridDimensions[2];
};

void AmoebaReferencePmeMultipoleForce::setPmeGridDimensions(vector<int>& pmeGridDimensions)
{

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

void AmoebaReferencePmeMultipoleForce::setPeriodicBoxSize(OpenMM::RealVec* vectors)
{

    if (vectors[0][0] == 0.0 || vectors[1][1] == 0.0 || vectors[2][2] == 0.0) {
        std::stringstream message;
        message << "Box size of zero is invalid.";
        throw OpenMMException(message.str());
    }

    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
    RealOpenMM determinant = vectors[0][0]*vectors[1][1]*vectors[2][2];
    assert(determinant > 0);
    RealOpenMM scale = 1.0/determinant;
    _recipBoxVectors[0] = RealVec(vectors[1][1]*vectors[2][2], 0, 0)*scale;
    _recipBoxVectors[1] = RealVec(-vectors[1][0]*vectors[2][2], vectors[0][0]*vectors[2][2], 0)*scale;
    _recipBoxVectors[2] = RealVec(vectors[1][0]*vectors[2][1]-vectors[1][1]*vectors[2][0], -vectors[0][0]*vectors[2][1], vectors[0][0]*vectors[1][1])*scale;
};

int compareInt2(const int2& v1, const int2& v2)
{
    return v1[1] < v2[1];
}

void AmoebaReferencePmeMultipoleForce::resizePmeArrays()
{

    _totalGridSize = _pmeGridDimensions[0]*_pmeGridDimensions[1]*_pmeGridDimensions[2];
    if (_pmeGridSize < _totalGridSize) {
        if (_pmeGrid) {
            delete _pmeGrid;
        }
        _pmeGrid      = new t_complex[_totalGridSize];
        _pmeGridSize  = _totalGridSize;
    }

    for (unsigned int ii = 0; ii < 3; ii++) {
       _pmeBsplineModuli[ii].resize(_pmeGridDimensions[ii]);
       _thetai[ii].resize(AMOEBA_PME_ORDER*_numParticles);
    }

    _iGrid.resize(_numParticles);
    _phi.resize(20*_numParticles);
    _phid.resize(10*_numParticles);
    _phip.resize(10*_numParticles);
    _phidp.resize(20*_numParticles);
}

void AmoebaReferencePmeMultipoleForce::initializePmeGrid()
{
    if (_pmeGrid == NULL)
        return;

    for (int jj = 0; jj < _totalGridSize; jj++) {
        _pmeGrid[jj].re = _pmeGrid[jj].im = 0.0;
    }
}

void AmoebaReferencePmeMultipoleForce::getPeriodicDelta(RealVec& deltaR) const
{
    deltaR -= _periodicBoxVectors[2]*FLOOR(deltaR[2]*_recipBoxVectors[2][2]+0.5);
    deltaR -= _periodicBoxVectors[1]*FLOOR(deltaR[1]*_recipBoxVectors[1][1]+0.5);
    deltaR -= _periodicBoxVectors[0]*FLOOR(deltaR[0]*_recipBoxVectors[0][0]+0.5);
}

void AmoebaReferencePmeMultipoleForce::getDampedInverseDistances(const MultipoleParticleData& particleI,
                                                                 const MultipoleParticleData& particleJ,
                                                                 RealOpenMM dscale, RealOpenMM pscale, RealOpenMM r,
                                                                 RealVec& dampedDInverseDistances,
                                                                 RealVec& dampedPInverseDistances) const
{

    RealVec scaleFactor    = RealVec(1.0, 1.0, 1.0);
    RealOpenMM damp        = particleI.dampingFactor*particleJ.dampingFactor;
    if (damp != 0.0) {

        RealOpenMM ratio   = (r/damp);
                   ratio   = ratio*ratio*ratio;

        RealOpenMM pgamma  = particleI.thole < particleJ.thole ? particleI.thole : particleJ.thole;
                   damp    = -pgamma*ratio;

        if (damp > -50.0) {
            RealOpenMM expdamp = EXP(damp);
            scaleFactor[0]     = 1.0 - expdamp;
            scaleFactor[1]     = 1.0 - expdamp*(1.0-damp);
            scaleFactor[2]     = 1.0 - expdamp*(1.0-damp+(0.6f*damp*damp));
        }
    }
    RealVec dampedDScale       = scaleFactor*dscale;

    RealOpenMM r2              = r*r;
    RealOpenMM r3              = r*r2;
    RealOpenMM r5              = r3*r2;
    RealOpenMM r7              = r5*r2;

    dampedDInverseDistances[0] =      (1.0-dampedDScale[0])/r3;
    dampedDInverseDistances[1] =  3.0*(1.0-dampedDScale[1])/r5;
    dampedDInverseDistances[2] = 15.0*(1.0-dampedDScale[2])/r7;
    if (pscale == dscale) {
        dampedPInverseDistances = dampedDInverseDistances;
    } else {
        RealVec dampedPScale       = scaleFactor*pscale;
        dampedPInverseDistances[0] =      (1.0-dampedPScale[0])/r3;
        dampedPInverseDistances[1] =  3.0*(1.0-dampedPScale[1])/r5;
        dampedPInverseDistances[2] = 15.0*(1.0-dampedPScale[2])/r7;
    }
}

void AmoebaReferencePmeMultipoleForce::initializeBSplineModuli()
{

    // Initialize the b-spline moduli.

    int maxSize = -1;
    for (unsigned int ii = 0; ii < 3; ii++) {
       _pmeBsplineModuli[ii].resize(_pmeGridDimensions[ii]);
        maxSize = maxSize  > _pmeGridDimensions[ii] ? maxSize : _pmeGridDimensions[ii];
    }

    RealOpenMM array[AMOEBA_PME_ORDER];
    RealOpenMM x = 0.0;
    array[0]     = 1.0 - x;
    array[1]     = x;
    for (int k = 2; k < AMOEBA_PME_ORDER; k++) {
        RealOpenMM denom = 1.0/k;
        array[k] = x*array[k-1]*denom;
        for (int i = 1; i < k; i++) {
            array[k-i] = ((x+i)*array[k-i-1] + ((k-i+1)-x)*array[k-i])*denom;
        }
        array[0] = (1.0-x)*array[0]*denom;
    }

    vector<RealOpenMM> bsarray(maxSize+1, 0.0);
    for (int i = 2; i <= AMOEBA_PME_ORDER+1; i++) {
        bsarray[i] = array[i-2];
    }
    for (int dim = 0; dim < 3; dim++) {

        int size = _pmeGridDimensions[dim];

        // get the modulus of the discrete Fourier transform

        RealOpenMM factor = 2.0*M_PI/size;
        for (int i = 0; i < size; i++) {
            RealOpenMM sum1 = 0.0;
            RealOpenMM sum2 = 0.0;
            for (int j = 1; j <= size; j++) {
                RealOpenMM arg = factor*i*(j-1);
                sum1          += bsarray[j]*COS(arg);
                sum2          += bsarray[j]*SIN(arg);
            }
            _pmeBsplineModuli[dim][i] = (sum1*sum1 + sum2*sum2);
        }

        // fix for exponential Euler spline interpolation failure

        RealOpenMM eps = 1.0e-7;
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
            RealOpenMM zeta;
            if (k == 0) {
                zeta = 1.0;
            } else {
                RealOpenMM sum1 = 1.0;
                RealOpenMM sum2 = 1.0;
                factor          = M_PI*k/size;
                for (int j = 1; j <= jcut; j++) {
                    RealOpenMM arg = factor/(factor+M_PI*j);
                    sum1           = sum1 + POW(arg,   AMOEBA_PME_ORDER);
                    sum2           = sum2 + POW(arg, 2*AMOEBA_PME_ORDER);
                }
                for (int j = 1; j <= jcut; j++) {
                    RealOpenMM arg  = factor/(factor-M_PI*j);
                    sum1           += POW(arg,   AMOEBA_PME_ORDER);
                    sum2           += POW(arg, 2*AMOEBA_PME_ORDER);
                }
                zeta = sum2/sum1;
            }
            _pmeBsplineModuli[dim][i-1] = _pmeBsplineModuli[dim][i-1]*(zeta*zeta);
        }
    }
}

void AmoebaReferencePmeMultipoleForce::calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI,
                                                                           const MultipoleParticleData& particleJ,
                                                                           RealOpenMM dscale, RealOpenMM pscale)
{

    unsigned int iIndex    = particleI.particleIndex;
    unsigned int jIndex    = particleJ.particleIndex;

    // compute the real space portion of the Ewald summation

    if (particleI.particleIndex == particleJ.particleIndex)
        return;

    RealVec deltaR    = particleJ.position - particleI.position;
    getPeriodicDelta(deltaR);
    RealOpenMM r2     = deltaR.dot(deltaR);

    if (r2 > _cutoffDistanceSquared)
        return;

    RealOpenMM r           = SQRT(r2);

    // calculate the error function damping terms

    RealOpenMM ralpha      = _alphaEwald*r;

    RealOpenMM bn0         = erfc(ralpha)/r;
    RealOpenMM alsq2       = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM alsq2n      = 1.0/(SQRT_PI*_alphaEwald);
    RealOpenMM exp2a       = EXP(-(ralpha*ralpha));
    alsq2n                *= alsq2;
    RealOpenMM bn1         = (bn0+alsq2n*exp2a)/r2;

    alsq2n                *= alsq2;
    RealOpenMM bn2         = (3.0*bn1+alsq2n*exp2a)/r2;

    alsq2n                *= alsq2;
    RealOpenMM bn3         = (5.0*bn2+alsq2n*exp2a)/r2;

    // compute the error function scaled and unscaled terms

    RealVec dampedDInverseDistances;
    RealVec dampedPInverseDistances;
    getDampedInverseDistances(particleI, particleJ, dscale, pscale, r, dampedDInverseDistances, dampedPInverseDistances);

    RealOpenMM drr3        = dampedDInverseDistances[0];
    RealOpenMM drr5        = dampedDInverseDistances[1];
    RealOpenMM drr7        = dampedDInverseDistances[2];

    RealOpenMM prr3        = dampedPInverseDistances[0];
    RealOpenMM prr5        = dampedPInverseDistances[1];
    RealOpenMM prr7        = dampedPInverseDistances[2];

    RealOpenMM dir         = particleI.dipole.dot(deltaR);

    RealVec qxI            = RealVec(particleI.quadrupole[QXX], particleI.quadrupole[QXY], particleI.quadrupole[QXZ]);
    RealVec qyI            = RealVec(particleI.quadrupole[QXY], particleI.quadrupole[QYY], particleI.quadrupole[QYZ]);
    RealVec qzI            = RealVec(particleI.quadrupole[QXZ], particleI.quadrupole[QYZ], particleI.quadrupole[QZZ]);

    RealVec qi             = RealVec(qxI.dot(deltaR), qyI.dot(deltaR), qzI.dot(deltaR));
    RealOpenMM qir         = qi.dot(deltaR);

    RealOpenMM djr         = particleJ.dipole.dot(deltaR);

    RealVec qxJ            = RealVec(particleJ.quadrupole[QXX], particleJ.quadrupole[QXY], particleJ.quadrupole[QXZ]);
    RealVec qyJ            = RealVec(particleJ.quadrupole[QXY], particleJ.quadrupole[QYY], particleJ.quadrupole[QYZ]);
    RealVec qzJ            = RealVec(particleJ.quadrupole[QXZ], particleJ.quadrupole[QYZ], particleJ.quadrupole[QZZ]);

    RealVec qj             = RealVec(qxJ.dot(deltaR), qyJ.dot(deltaR), qzJ.dot(deltaR));
    RealOpenMM qjr         = qj.dot(deltaR);

    RealVec fim            = qj*(2.0*bn2)  - particleJ.dipole*bn1  - deltaR*(bn1*particleJ.charge - bn2*djr+bn3*qjr);
    RealVec fjm            = qi*(-2.0*bn2)  - particleI.dipole*bn1  + deltaR*(bn1*particleI.charge + bn2*dir+bn3*qir);

    RealVec fid            = qj*(2.0*drr5) - particleJ.dipole*drr3 - deltaR*(drr3*particleJ.charge - drr5*djr+drr7*qjr);
    RealVec fjd            = qi*(-2.0*drr5) - particleI.dipole*drr3 + deltaR*(drr3*particleI.charge + drr5*dir+drr7*qir);

    RealVec fip            = qj*(2.0*prr5) - particleJ.dipole*prr3 - deltaR*(prr3*particleJ.charge - prr5*djr+prr7*qjr);
    RealVec fjp            = qi*(-2.0*prr5) - particleI.dipole*prr3 + deltaR*(prr3*particleI.charge + prr5*dir+prr7*qir);

    // increment the field at each site due to this interaction


    _fixedMultipoleField[iIndex]      += fim - fid;
    _fixedMultipoleField[jIndex]      += fjm - fjd;

    _fixedMultipoleFieldPolar[iIndex] += fim - fip;
    _fixedMultipoleFieldPolar[jIndex] += fjm - fjp;
}

void AmoebaReferencePmeMultipoleForce::calculateFixedMultipoleField(const vector<MultipoleParticleData>& particleData)
{

    // first calculate reciprocal space fixed multipole fields

    resizePmeArrays();
    computeAmoebaBsplines(particleData);
    initializePmeGrid();
    spreadFixedMultipolesOntoGrid(particleData);
    fftpack_exec_3d(_fftplan, FFTPACK_FORWARD, _pmeGrid, _pmeGrid);
    performAmoebaReciprocalConvolution();
    fftpack_exec_3d(_fftplan, FFTPACK_BACKWARD, _pmeGrid, _pmeGrid);
    computeFixedPotentialFromGrid();
    recordFixedMultipoleField();

    // include self-energy portion of the multipole field
    // and initialize _fixedMultipoleFieldPolar to _fixedMultipoleField

    RealOpenMM term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (unsigned int jj = 0; jj < _numParticles; jj++) {
        RealVec selfEnergy             = particleData[jj].dipole*term;
        _fixedMultipoleField[jj]      += selfEnergy;
        _fixedMultipoleFieldPolar[jj]  = _fixedMultipoleField[jj];
    }

    // include direct space fixed multipole fields

    this->AmoebaReferenceMultipoleForce::calculateFixedMultipoleField(particleData);
}

#define ARRAY(x,y) array[(x)-1+((y)-1)*AMOEBA_PME_ORDER]

/**
 * This is called from computeBsplines().  It calculates the spline coefficients for a single atom along a single axis.
 */
void AmoebaReferencePmeMultipoleForce::computeBSplinePoint(vector<RealOpenMM4>& thetai, RealOpenMM w)
{

    RealOpenMM array[AMOEBA_PME_ORDER*AMOEBA_PME_ORDER];

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
        RealOpenMM denom = 1.0 / k;
        ARRAY(i,i) = denom * w * ARRAY(k,k);
        for (int j = 1; j <= i-2; j++) {
            ARRAY(i,i-j) = denom * ((w+j)*ARRAY(k,i-j-1)+(i-j-w)*ARRAY(k,i-j));
        }
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

    for (int i = 1; i <= AMOEBA_PME_ORDER; i++) {
        thetai[i-1] = RealOpenMM4(ARRAY(AMOEBA_PME_ORDER,i), ARRAY(AMOEBA_PME_ORDER-1,i), ARRAY(AMOEBA_PME_ORDER-2,i), ARRAY(AMOEBA_PME_ORDER-3,i));
    }
}

/**
 * Compute b-spline coefficients.
 */
void AmoebaReferencePmeMultipoleForce::computeAmoebaBsplines(const vector<MultipoleParticleData>& particleData)
{
    //  get the B-spline coefficients for each multipole site

    for (unsigned int ii = 0; ii < _numParticles; ii++) {
        RealVec position  = particleData[ii].position;
        getPeriodicDelta(position);
        IntVec igrid;
        for (unsigned int jj = 0; jj < 3; jj++) {

            RealOpenMM w  = position[0]*_recipBoxVectors[0][jj]+position[1]*_recipBoxVectors[1][jj]+position[2]*_recipBoxVectors[2][jj];
            RealOpenMM fr = _pmeGridDimensions[jj]*(w-(int)(w+0.5)+0.5);
            int ifr       = static_cast<int>(floor(fr));
            w             = fr - ifr;
            igrid[jj]     = ifr - AMOEBA_PME_ORDER + 1;
            igrid[jj]    += igrid[jj] < 0 ? _pmeGridDimensions[jj] : 0;
            vector<RealOpenMM4> thetaiTemp(AMOEBA_PME_ORDER);
            computeBSplinePoint(thetaiTemp, w);
            for (unsigned int kk = 0; kk < AMOEBA_PME_ORDER; kk++) {
                _thetai[jj][ii*AMOEBA_PME_ORDER+kk] = thetaiTemp[kk];
            }
        }

        // Record the grid point.

        _iGrid[ii]               = igrid;
    }
}

void AmoebaReferencePmeMultipoleForce::transformMultipolesToFractionalCoordinates(const vector<MultipoleParticleData>& particleData) {
    // Build matrices for transforming the dipoles and quadrupoles.

    RealVec a[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a[j][i] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    int index1[] = {0, 0, 0, 1, 1, 2};
    int index2[] = {0, 1, 2, 1, 2, 2};
    RealOpenMM b[6][6];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            b[i][j] = a[index1[i]][index1[j]]*a[index2[i]][index2[j]];
            if (index1[i] != index2[i])
                b[i][j] += a[index1[i]][index2[j]]*a[index2[i]][index1[j]];
        }
    }

    // Transform the multipoles.

    _transformed.resize(particleData.size());
    double quadScale[] = {1, 2, 2, 1, 2, 1};
    for (int i = 0; i < (int) particleData.size(); i++) {
        _transformed[i].charge = particleData[i].charge;
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

void AmoebaReferencePmeMultipoleForce::transformPotentialToCartesianCoordinates(const vector<RealOpenMM>& fphi, vector<RealOpenMM>& cphi) const {
    // Build a matrix for transforming the potential.

    RealVec a[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    int index1[] = {0, 1, 2, 0, 0, 1};
    int index2[] = {0, 1, 2, 1, 2, 2};
    RealOpenMM b[6][6];
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

void AmoebaReferencePmeMultipoleForce::spreadFixedMultipolesOntoGrid(const vector<MultipoleParticleData>& particleData)
{

    transformMultipolesToFractionalCoordinates(particleData);

    // Clear the grid.

    for (int gridIndex = 0; gridIndex < _totalGridSize; gridIndex++)
        _pmeGrid[gridIndex] = t_complex(0, 0);

    // Loop over atoms and spread them on the grid.

    for (int atomIndex = 0; atomIndex < _numParticles; atomIndex++) {
        RealOpenMM atomCharge       = _transformed[atomIndex].charge;
        RealVec atomDipole          = RealVec(_transformed[atomIndex].dipole[0],
                                              _transformed[atomIndex].dipole[1],
                                              _transformed[atomIndex].dipole[2]);

        RealOpenMM atomQuadrupoleXX = _transformed[atomIndex].quadrupole[QXX];
        RealOpenMM atomQuadrupoleXY = _transformed[atomIndex].quadrupole[QXY];
        RealOpenMM atomQuadrupoleXZ = _transformed[atomIndex].quadrupole[QXZ];
        RealOpenMM atomQuadrupoleYY = _transformed[atomIndex].quadrupole[QYY];
        RealOpenMM atomQuadrupoleYZ = _transformed[atomIndex].quadrupole[QYZ];
        RealOpenMM atomQuadrupoleZZ = _transformed[atomIndex].quadrupole[QZZ];
        IntVec& gridPoint = _iGrid[atomIndex];
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
            int x = (gridPoint[0]+ix) % _pmeGridDimensions[0];
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int y = (gridPoint[1]+iy) % _pmeGridDimensions[1];
                for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
                    int z = (gridPoint[2]+iz) % _pmeGridDimensions[2];
                    RealOpenMM4 t = _thetai[0][atomIndex*AMOEBA_PME_ORDER+ix];
                    RealOpenMM4 u = _thetai[1][atomIndex*AMOEBA_PME_ORDER+iy];
                    RealOpenMM4 v = _thetai[2][atomIndex*AMOEBA_PME_ORDER+iz];
                    RealOpenMM term0 = atomCharge*u[0]*v[0] + atomDipole[1]*u[1]*v[0] + atomDipole[2]*u[0]*v[1] + atomQuadrupoleYY*u[2]*v[0] + atomQuadrupoleZZ*u[0]*v[2] + atomQuadrupoleYZ*u[1]*v[1];
                    RealOpenMM term1 = atomDipole[0]*u[0]*v[0] + atomQuadrupoleXY*u[1]*v[0] + atomQuadrupoleXZ*u[0]*v[1];
                    RealOpenMM term2 = atomQuadrupoleXX * u[0] * v[0];
                    t_complex& gridValue = _pmeGrid[x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z];
                    gridValue.re += term0*t[0] + term1*t[1] + term2*t[2];
                }
            }
        }
    }
}

void AmoebaReferencePmeMultipoleForce::performAmoebaReciprocalConvolution()
{

    RealOpenMM expFactor   = (M_PI*M_PI)/(_alphaEwald*_alphaEwald);
    RealOpenMM scaleFactor = 1.0/(M_PI*_periodicBoxVectors[0][0]*_periodicBoxVectors[1][1]*_periodicBoxVectors[2][2]);

    for (int index = 0; index < _totalGridSize; index++)
    {
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

        RealOpenMM mhx = mx*_recipBoxVectors[0][0];
        RealOpenMM mhy = mx*_recipBoxVectors[1][0]+my*_recipBoxVectors[1][1];
        RealOpenMM mhz = mx*_recipBoxVectors[2][0]+my*_recipBoxVectors[2][1]+mz*_recipBoxVectors[2][2];

        RealOpenMM bx = _pmeBsplineModuli[0][kx];
        RealOpenMM by = _pmeBsplineModuli[1][ky];
        RealOpenMM bz = _pmeBsplineModuli[2][kz];

        RealOpenMM m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        RealOpenMM denom = m2*bx*by*bz;
        RealOpenMM eterm = scaleFactor*EXP(-expFactor*m2)/denom;

        _pmeGrid[index].re *= eterm;
        _pmeGrid[index].im *= eterm;
    }
}

void AmoebaReferencePmeMultipoleForce::computeFixedPotentialFromGrid()
{
    // extract the permanent multipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        IntVec gridPoint = _iGrid[m];
        RealOpenMM tuv000 = 0.0;
        RealOpenMM tuv001 = 0.0;
        RealOpenMM tuv010 = 0.0;
        RealOpenMM tuv100 = 0.0;
        RealOpenMM tuv200 = 0.0;
        RealOpenMM tuv020 = 0.0;
        RealOpenMM tuv002 = 0.0;
        RealOpenMM tuv110 = 0.0;
        RealOpenMM tuv101 = 0.0;
        RealOpenMM tuv011 = 0.0;
        RealOpenMM tuv300 = 0.0;
        RealOpenMM tuv030 = 0.0;
        RealOpenMM tuv003 = 0.0;
        RealOpenMM tuv210 = 0.0;
        RealOpenMM tuv201 = 0.0;
        RealOpenMM tuv120 = 0.0;
        RealOpenMM tuv021 = 0.0;
        RealOpenMM tuv102 = 0.0;
        RealOpenMM tuv012 = 0.0;
        RealOpenMM tuv111 = 0.0;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            RealOpenMM4 v = _thetai[2][m*AMOEBA_PME_ORDER+iz];
            RealOpenMM tu00 = 0.0;
            RealOpenMM tu10 = 0.0;
            RealOpenMM tu01 = 0.0;
            RealOpenMM tu20 = 0.0;
            RealOpenMM tu11 = 0.0;
            RealOpenMM tu02 = 0.0;
            RealOpenMM tu30 = 0.0;
            RealOpenMM tu21 = 0.0;
            RealOpenMM tu12 = 0.0;
            RealOpenMM tu03 = 0.0;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                RealOpenMM4 u = _thetai[1][m*AMOEBA_PME_ORDER+iy];
                RealOpenMM4 t = RealOpenMM4(0.0, 0.0, 0.0, 0.0);
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    RealOpenMM tq = _pmeGrid[gridIndex].re;
                    RealOpenMM4 tadd = _thetai[0][m*AMOEBA_PME_ORDER+ix];
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

void AmoebaReferencePmeMultipoleForce::spreadInducedDipolesOnGrid(const vector<RealVec>& inputInducedDipole,
                                                                  const vector<RealVec>& inputInducedDipolePolar) {
    // Create the matrix to convert from Cartesian to fractional coordinates.

    RealVec cartToFrac[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cartToFrac[j][i] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];

    // Clear the grid.

    for (int gridIndex = 0; gridIndex < _totalGridSize; gridIndex++)
        _pmeGrid[gridIndex] = t_complex(0, 0);

    // Loop over atoms and spread them on the grid.

    for (int atomIndex = 0; atomIndex < _numParticles; atomIndex++) {
        RealVec inducedDipole = RealVec(inputInducedDipole[atomIndex][0]*cartToFrac[0][0] + inputInducedDipole[atomIndex][1]*cartToFrac[0][1] + inputInducedDipole[atomIndex][2]*cartToFrac[0][2],
                                        inputInducedDipole[atomIndex][0]*cartToFrac[1][0] + inputInducedDipole[atomIndex][1]*cartToFrac[1][1] + inputInducedDipole[atomIndex][2]*cartToFrac[1][2],
                                        inputInducedDipole[atomIndex][0]*cartToFrac[2][0] + inputInducedDipole[atomIndex][1]*cartToFrac[2][1] + inputInducedDipole[atomIndex][2]*cartToFrac[2][2]);
        RealVec inducedDipolePolar = RealVec(inputInducedDipolePolar[atomIndex][0]*cartToFrac[0][0] + inputInducedDipolePolar[atomIndex][1]*cartToFrac[0][1] + inputInducedDipolePolar[atomIndex][2]*cartToFrac[0][2],
                                             inputInducedDipolePolar[atomIndex][0]*cartToFrac[1][0] + inputInducedDipolePolar[atomIndex][1]*cartToFrac[1][1] + inputInducedDipolePolar[atomIndex][2]*cartToFrac[1][2],
                                             inputInducedDipolePolar[atomIndex][0]*cartToFrac[2][0] + inputInducedDipolePolar[atomIndex][1]*cartToFrac[2][1] + inputInducedDipolePolar[atomIndex][2]*cartToFrac[2][2]);
        IntVec& gridPoint = _iGrid[atomIndex];
        for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
            int x = (gridPoint[0]+ix) % _pmeGridDimensions[0];
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int y = (gridPoint[1]+iy) % _pmeGridDimensions[1];
                for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
                    int z = (gridPoint[2]+iz) % _pmeGridDimensions[2];

                    RealOpenMM4 t = _thetai[0][atomIndex*AMOEBA_PME_ORDER+ix];
                    RealOpenMM4 u = _thetai[1][atomIndex*AMOEBA_PME_ORDER+iy];
                    RealOpenMM4 v = _thetai[2][atomIndex*AMOEBA_PME_ORDER+iz];

                    RealOpenMM term01 = inducedDipole[1]*u[1]*v[0] + inducedDipole[2]*u[0]*v[1];
                    RealOpenMM term11 = inducedDipole[0]*u[0]*v[0];
                    RealOpenMM term02 = inducedDipolePolar[1]*u[1]*v[0] + inducedDipolePolar[2]*u[0]*v[1];
                    RealOpenMM term12 = inducedDipolePolar[0]*u[0]*v[0];

                    t_complex& gridValue = _pmeGrid[x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z];
                    gridValue.re += term01*t[0] + term11*t[1];
                    gridValue.im += term02*t[0] + term12*t[1];
                }
            }
        }
    }
}

void AmoebaReferencePmeMultipoleForce::computeInducedPotentialFromGrid()
{
    // extract the induced dipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        IntVec gridPoint = _iGrid[m];
        RealOpenMM tuv100_1 = 0.0;
        RealOpenMM tuv010_1 = 0.0;
        RealOpenMM tuv001_1 = 0.0;
        RealOpenMM tuv200_1 = 0.0;
        RealOpenMM tuv020_1 = 0.0;
        RealOpenMM tuv002_1 = 0.0;
        RealOpenMM tuv110_1 = 0.0;
        RealOpenMM tuv101_1 = 0.0;
        RealOpenMM tuv011_1 = 0.0;
        RealOpenMM tuv100_2 = 0.0;
        RealOpenMM tuv010_2 = 0.0;
        RealOpenMM tuv001_2 = 0.0;
        RealOpenMM tuv200_2 = 0.0;
        RealOpenMM tuv020_2 = 0.0;
        RealOpenMM tuv002_2 = 0.0;
        RealOpenMM tuv110_2 = 0.0;
        RealOpenMM tuv101_2 = 0.0;
        RealOpenMM tuv011_2 = 0.0;
        RealOpenMM tuv000 = 0.0;
        RealOpenMM tuv001 = 0.0;
        RealOpenMM tuv010 = 0.0;
        RealOpenMM tuv100 = 0.0;
        RealOpenMM tuv200 = 0.0;
        RealOpenMM tuv020 = 0.0;
        RealOpenMM tuv002 = 0.0;
        RealOpenMM tuv110 = 0.0;
        RealOpenMM tuv101 = 0.0;
        RealOpenMM tuv011 = 0.0;
        RealOpenMM tuv300 = 0.0;
        RealOpenMM tuv030 = 0.0;
        RealOpenMM tuv003 = 0.0;
        RealOpenMM tuv210 = 0.0;
        RealOpenMM tuv201 = 0.0;
        RealOpenMM tuv120 = 0.0;
        RealOpenMM tuv021 = 0.0;
        RealOpenMM tuv102 = 0.0;
        RealOpenMM tuv012 = 0.0;
        RealOpenMM tuv111 = 0.0;
        for (int iz = 0; iz < AMOEBA_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            RealOpenMM4 v = _thetai[2][m*AMOEBA_PME_ORDER+iz];
            RealOpenMM tu00_1 = 0.0;
            RealOpenMM tu01_1 = 0.0;
            RealOpenMM tu10_1 = 0.0;
            RealOpenMM tu20_1 = 0.0;
            RealOpenMM tu11_1 = 0.0;
            RealOpenMM tu02_1 = 0.0;
            RealOpenMM tu00_2 = 0.0;
            RealOpenMM tu01_2 = 0.0;
            RealOpenMM tu10_2 = 0.0;
            RealOpenMM tu20_2 = 0.0;
            RealOpenMM tu11_2 = 0.0;
            RealOpenMM tu02_2 = 0.0;
            RealOpenMM tu00 = 0.0;
            RealOpenMM tu10 = 0.0;
            RealOpenMM tu01 = 0.0;
            RealOpenMM tu20 = 0.0;
            RealOpenMM tu11 = 0.0;
            RealOpenMM tu02 = 0.0;
            RealOpenMM tu30 = 0.0;
            RealOpenMM tu21 = 0.0;
            RealOpenMM tu12 = 0.0;
            RealOpenMM tu03 = 0.0;
            for (int iy = 0; iy < AMOEBA_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                RealOpenMM4 u = _thetai[1][m*AMOEBA_PME_ORDER+iy];
                RealOpenMM t0_1 = 0.0;
                RealOpenMM t1_1 = 0.0;
                RealOpenMM t2_1 = 0.0;
                RealOpenMM t0_2 = 0.0;
                RealOpenMM t1_2 = 0.0;
                RealOpenMM t2_2 = 0.0;
                RealOpenMM t3 = 0.0;
                for (int ix = 0; ix < AMOEBA_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    t_complex tq = _pmeGrid[gridIndex];
                    RealOpenMM4 tadd = _thetai[0][m*AMOEBA_PME_ORDER+ix];
                    t0_1 += tq.re*tadd[0];
                    t1_1 += tq.re*tadd[1];
                    t2_1 += tq.re*tadd[2];
                    t0_2 += tq.im*tadd[0];
                    t1_2 += tq.im*tadd[1];
                    t2_2 += tq.im*tadd[2];
                    t3 += (tq.re+tq.im)*tadd[3];
                }
                tu00_1 += t0_1*u[0];
                tu10_1 += t1_1*u[0];
                tu01_1 += t0_1*u[1];
                tu20_1 += t2_1*u[0];
                tu11_1 += t1_1*u[1];
                tu02_1 += t0_1*u[2];
                tu00_2 += t0_2*u[0];
                tu10_2 += t1_2*u[0];
                tu01_2 += t0_2*u[1];
                tu20_2 += t2_2*u[0];
                tu11_2 += t1_2*u[1];
                tu02_2 += t0_2*u[2];
                RealOpenMM t0 = t0_1 + t0_2;
                RealOpenMM t1 = t1_1 + t1_2;
                RealOpenMM t2 = t2_1 + t2_2;
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
            tuv100_1 += tu10_1*v[0];
            tuv010_1 += tu01_1*v[0];
            tuv001_1 += tu00_1*v[1];
            tuv200_1 += tu20_1*v[0];
            tuv020_1 += tu02_1*v[0];
            tuv002_1 += tu00_1*v[2];
            tuv110_1 += tu11_1*v[0];
            tuv101_1 += tu10_1*v[1];
            tuv011_1 += tu01_1*v[1];
            tuv100_2 += tu10_2*v[0];
            tuv010_2 += tu01_2*v[0];
            tuv001_2 += tu00_2*v[1];
            tuv200_2 += tu20_2*v[0];
            tuv020_2 += tu02_2*v[0];
            tuv002_2 += tu00_2*v[2];
            tuv110_2 += tu11_2*v[0];
            tuv101_2 += tu10_2*v[1];
            tuv011_2 += tu01_2*v[1];
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
        _phid[10*m]   = 0.0;
        _phid[10*m+1] = tuv100_1;
        _phid[10*m+2] = tuv010_1;
        _phid[10*m+3] = tuv001_1;
        _phid[10*m+4] = tuv200_1;
        _phid[10*m+5] = tuv020_1;
        _phid[10*m+6] = tuv002_1;
        _phid[10*m+7] = tuv110_1;
        _phid[10*m+8] = tuv101_1;
        _phid[10*m+9] = tuv011_1;

        _phip[10*m]   = 0.0;
        _phip[10*m+1] = tuv100_2;
        _phip[10*m+2] = tuv010_2;
        _phip[10*m+3] = tuv001_2;
        _phip[10*m+4] = tuv200_2;
        _phip[10*m+5] = tuv020_2;
        _phip[10*m+6] = tuv002_2;
        _phip[10*m+7] = tuv110_2;
        _phip[10*m+8] = tuv101_2;
        _phip[10*m+9] = tuv011_2;

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

RealOpenMM AmoebaReferencePmeMultipoleForce::computeReciprocalSpaceFixedMultipoleForceAndEnergy(const vector<MultipoleParticleData>& particleData,
                                                                                                vector<RealVec>& forces, vector<RealVec>& torques) const
{
    RealOpenMM multipole[10];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    vector<RealOpenMM> cphi(10*_numParticles);
    transformPotentialToCartesianCoordinates(_phi, cphi);
    RealVec fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    RealOpenMM energy = 0.0;
    for (int i = 0; i < _numParticles; i++) {

        // Compute the torque.

        multipole[0] = particleData[i].charge;

        multipole[1] =  particleData[i].dipole[0];
        multipole[2] =  particleData[i].dipole[1];
        multipole[3] =  particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];

        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        const RealOpenMM* phi = &cphi[10*i];
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

        RealVec f = RealVec(0.0, 0.0, 0.0);
        for (int k = 0; k < 10; k++) {
            energy += multipole[k]*_phi[20*i+k];
            f[0]   += multipole[k]*_phi[20*i+deriv1[k]];
            f[1]   += multipole[k]*_phi[20*i+deriv2[k]];
            f[2]   += multipole[k]*_phi[20*i+deriv3[k]];
        }
        f              *= (_electric);
        forces[i]      -= RealVec(f[0]*fracToCart[0][0] + f[1]*fracToCart[0][1] + f[2]*fracToCart[0][2],
                                  f[0]*fracToCart[1][0] + f[1]*fracToCart[1][1] + f[2]*fracToCart[1][2],
                                  f[0]*fracToCart[2][0] + f[1]*fracToCart[2][1] + f[2]*fracToCart[2][2]);

    }
    return (0.5*_electric*energy);
}

/**
 * Compute the forces due to the reciprocal space PME calculation for induced dipoles.
 */
RealOpenMM AmoebaReferencePmeMultipoleForce::computeReciprocalSpaceInducedDipoleForceAndEnergy(AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                                               const vector<MultipoleParticleData>& particleData,
                                                                                               vector<RealVec>& forces, vector<RealVec>& torques) const
{

    RealOpenMM multipole[10];
    RealOpenMM inducedDipole[3];
    RealOpenMM inducedDipolePolar[3];
    RealOpenMM scales[3];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    vector<RealOpenMM> cphi(10*_numParticles);
    transformPotentialToCartesianCoordinates(_phidp, cphi);
    RealVec cartToFrac[3], fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cartToFrac[j][i] = fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    RealOpenMM energy = 0.0;
    for (int i = 0; i < _numParticles; i++) {

        // Compute the torque.

        unsigned int iIndex = particleData[i].particleIndex;

        multipole[0] = particleData[i].charge;

        multipole[1] = particleData[i].dipole[0];
        multipole[2] = particleData[i].dipole[1];
        multipole[3] = particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];
        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        const RealOpenMM* phi = &cphi[10*i];
        torques[iIndex][0] += 0.5*_electric*(multipole[3]*phi[2] - multipole[2]*phi[3]
                      + 2.0*(multipole[6]-multipole[5])*phi[9]
                      + multipole[8]*phi[7] + multipole[9]*phi[5]
                      - multipole[7]*phi[8] - multipole[9]*phi[6]);

        torques[iIndex][1] += 0.5*_electric*(multipole[1]*phi[3] - multipole[3]*phi[1]
                      + 2.0*(multipole[4]-multipole[6])*phi[8]
                      + multipole[7]*phi[9] + multipole[8]*phi[6]
                      - multipole[8]*phi[4] - multipole[9]*phi[7]);

        torques[iIndex][2] += 0.5*_electric*(multipole[2]*phi[1] - multipole[1]*phi[2]
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
        inducedDipolePolar[0] = _inducedDipolePolar[i][0]*cartToFrac[0][0] + _inducedDipolePolar[i][1]*cartToFrac[0][1] + _inducedDipolePolar[i][2]*cartToFrac[0][2];
        inducedDipolePolar[1] = _inducedDipolePolar[i][0]*cartToFrac[1][0] + _inducedDipolePolar[i][1]*cartToFrac[1][1] + _inducedDipolePolar[i][2]*cartToFrac[1][2];
        inducedDipolePolar[2] = _inducedDipolePolar[i][0]*cartToFrac[2][0] + _inducedDipolePolar[i][1]*cartToFrac[2][1] + _inducedDipolePolar[i][2]*cartToFrac[2][2];

        energy += (inducedDipole[0]+inducedDipolePolar[0])*_phi[20*i+1];
        energy += (inducedDipole[1]+inducedDipolePolar[1])*_phi[20*i+2];
        energy += (inducedDipole[2]+inducedDipolePolar[2])*_phi[20*i+3];

        RealVec f        = RealVec(0.0, 0.0, 0.0);

        for (int k = 0; k < 3; k++) {

            int j1 = deriv1[k+1];
            int j2 = deriv2[k+1];
            int j3 = deriv3[k+1];

            f[0] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j1];
            f[1] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j2];
            f[2] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j3];
 
            if (polarizationType == AmoebaReferenceMultipoleForce::Mutual) {
                f[0] += (inducedDipole[k]*_phip[10*i+j1] + inducedDipolePolar[k]*_phid[10*i+j1]);
                f[1] += (inducedDipole[k]*_phip[10*i+j2] + inducedDipolePolar[k]*_phid[10*i+j2]);
                f[2] += (inducedDipole[k]*_phip[10*i+j3] + inducedDipolePolar[k]*_phid[10*i+j3]);
            }
        }

        for (int k = 0; k < 10; k++) {
            f[0] += multipole[k]*_phidp[20*i+deriv1[k]];
            f[1] += multipole[k]*_phidp[20*i+deriv2[k]];
            f[2] += multipole[k]*_phidp[20*i+deriv3[k]];
        }

        f              *= (0.5*_electric);
        forces[iIndex] -= RealVec(f[0]*fracToCart[0][0] + f[1]*fracToCart[0][1] + f[2]*fracToCart[0][2],
                                  f[0]*fracToCart[1][0] + f[1]*fracToCart[1][1] + f[2]*fracToCart[1][2],
                                  f[0]*fracToCart[2][0] + f[1]*fracToCart[2][1] + f[2]*fracToCart[2][2]);
    }
    return (0.25*_electric*energy);
}

void AmoebaReferencePmeMultipoleForce::recordFixedMultipoleField()
{
    RealVec fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    for (int i = 0; i < _numParticles; i++) {
        _fixedMultipoleField[i][0] = -(_phi[20*i+1]*fracToCart[0][0] + _phi[20*i+2]*fracToCart[0][1] + _phi[20*i+3]*fracToCart[0][2]);
        _fixedMultipoleField[i][1] = -(_phi[20*i+1]*fracToCart[1][0] + _phi[20*i+2]*fracToCart[1][1] + _phi[20*i+3]*fracToCart[1][2]);
        _fixedMultipoleField[i][2] = -(_phi[20*i+1]*fracToCart[2][0] + _phi[20*i+2]*fracToCart[2][1] + _phi[20*i+3]*fracToCart[2][2]);
    }
}

void AmoebaReferencePmeMultipoleForce::initializeInducedDipoles(vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    this->AmoebaReferenceMultipoleForce::initializeInducedDipoles(updateInducedDipoleFields);
    calculateReciprocalSpaceInducedDipoleField(updateInducedDipoleFields);
}

void AmoebaReferencePmeMultipoleForce::recordInducedDipoleField(vector<RealVec>& field, vector<RealVec>& fieldPolar)
{
    RealVec fracToCart[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];
    for (int i = 0; i < _numParticles; i++) {

        field[i][0] -= _phid[10*i+1]*fracToCart[0][0] + _phid[10*i+2]*fracToCart[0][1] + _phid[10*i+3]*fracToCart[0][2];
        field[i][1] -= _phid[10*i+1]*fracToCart[1][0] + _phid[10*i+2]*fracToCart[1][1] + _phid[10*i+3]*fracToCart[1][2];
        field[i][2] -= _phid[10*i+1]*fracToCart[2][0] + _phid[10*i+2]*fracToCart[2][1] + _phid[10*i+3]*fracToCart[2][2];

        fieldPolar[i][0] -= _phip[10*i+1]*fracToCart[0][0] + _phip[10*i+2]*fracToCart[0][1] + _phip[10*i+3]*fracToCart[0][2];
        fieldPolar[i][1] -= _phip[10*i+1]*fracToCart[1][0] + _phip[10*i+2]*fracToCart[1][1] + _phip[10*i+3]*fracToCart[1][2];
        fieldPolar[i][2] -= _phip[10*i+1]*fracToCart[2][0] + _phip[10*i+2]*fracToCart[2][1] + _phip[10*i+3]*fracToCart[2][2];
    }
}

void AmoebaReferencePmeMultipoleForce::calculateReciprocalSpaceInducedDipoleField(vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{
    // Perform PME for the induced dipoles.

    initializePmeGrid();
    spreadInducedDipolesOnGrid(*updateInducedDipoleFields[0].inducedDipoles, *updateInducedDipoleFields[1].inducedDipoles);
    fftpack_exec_3d(_fftplan, FFTPACK_FORWARD, _pmeGrid, _pmeGrid);
    performAmoebaReciprocalConvolution();
    fftpack_exec_3d(_fftplan, FFTPACK_BACKWARD, _pmeGrid, _pmeGrid);
    computeInducedPotentialFromGrid();
    recordInducedDipoleField(updateInducedDipoleFields[0].inducedDipoleField, updateInducedDipoleFields[1].inducedDipoleField);
}

void AmoebaReferencePmeMultipoleForce::calculateInducedDipoleFields(const vector<MultipoleParticleData>& particleData,
                                                                     vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{
    // Initialize the fields to zero.

    RealVec zeroVec(0.0, 0.0, 0.0);
    for (unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++)
        std::fill(updateInducedDipoleFields[ii].inducedDipoleField.begin(), updateInducedDipoleFields[ii].inducedDipoleField.end(), zeroVec);

    // Add fields from direct space interactions.

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii + 1; jj < particleData.size(); jj++) {
            calculateDirectInducedDipolePairIxns(particleData[ii], particleData[jj], updateInducedDipoleFields);
        }
    }

    // reciprocal space ixns

    calculateReciprocalSpaceInducedDipoleField(updateInducedDipoleFields);

    if(getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
        // While we have the reciprocal space (fractional coordinate) field gradient available, add it to the real space
        // terms computed above, after transforming to Cartesian coordinates.  This allows real and reciprocal space
        // dipole response force contributions to be computed together.
        RealVec fracToCart[3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                fracToCart[i][j] = _pmeGridDimensions[j]*_recipBoxVectors[i][j];


        for (int i = 0; i < _numParticles; i++) {
            RealOpenMM EmatD[3][3] = {
                { _phid[10*i+4], _phid[10*i+7], _phid[10*i+8] },
                { _phid[10*i+7], _phid[10*i+5], _phid[10*i+9] },
                { _phid[10*i+8], _phid[10*i+9], _phid[10*i+6] }
            };

            RealOpenMM Exx = 0.0, Eyy = 0.0, Ezz = 0.0, Exy = 0.0, Exz = 0.0, Eyz = 0.0;
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    Exx += fracToCart[0][k] * EmatD[k][l] * fracToCart[0][l];
                    Eyy += fracToCart[1][k] * EmatD[k][l] * fracToCart[1][l];
                    Ezz += fracToCart[2][k] * EmatD[k][l] * fracToCart[2][l];
                    Exy += fracToCart[0][k] * EmatD[k][l] * fracToCart[1][l];
                    Exz += fracToCart[0][k] * EmatD[k][l] * fracToCart[2][l];
                    Eyz += fracToCart[1][k] * EmatD[k][l] * fracToCart[2][l];
                }
            }
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][0] -= Exx;
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][1] -= Eyy;
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][2] -= Ezz;
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][3] -= Exy;
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][4] -= Exz;
            updateInducedDipoleFields[0].inducedDipoleFieldGradient[i][5] -= Eyz;

            RealOpenMM EmatP[3][3] = {
                { _phip[10*i+4], _phip[10*i+7], _phip[10*i+8] },
                { _phip[10*i+7], _phip[10*i+5], _phip[10*i+9] },
                { _phip[10*i+8], _phip[10*i+9], _phip[10*i+6] }
            };

            Exx = 0.0; Eyy = 0.0; Ezz = 0.0; Exy = 0.0; Exz = 0.0; Eyz = 0.0;
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    Exx += fracToCart[0][k] * EmatP[k][l] * fracToCart[0][l];
                    Eyy += fracToCart[1][k] * EmatP[k][l] * fracToCart[1][l];
                    Ezz += fracToCart[2][k] * EmatP[k][l] * fracToCart[2][l];
                    Exy += fracToCart[0][k] * EmatP[k][l] * fracToCart[1][l];
                    Exz += fracToCart[0][k] * EmatP[k][l] * fracToCart[2][l];
                    Eyz += fracToCart[1][k] * EmatP[k][l] * fracToCart[2][l];
                }
            }

            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][0] -= Exx;
            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][1] -= Eyy;
            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][2] -= Ezz;
            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][3] -= Exy;
            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][4] -= Exz;
            updateInducedDipoleFields[1].inducedDipoleFieldGradient[i][5] -= Eyz;
        }
    }

    // self ixn

    RealOpenMM term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++) {
        vector<RealVec>& inducedDipoles = *updateInducedDipoleFields[ii].inducedDipoles;
        vector<RealVec>& field          = updateInducedDipoleFields[ii].inducedDipoleField;
        for (unsigned int jj = 0; jj < particleData.size(); jj++) {
            field[jj] += inducedDipoles[jj]*term;
        }
    }
}

void AmoebaReferencePmeMultipoleForce::calculateDirectInducedDipolePairIxn(unsigned int iIndex, unsigned int jIndex,
                                                                           RealOpenMM preFactor1, RealOpenMM preFactor2,
                                                                           const RealVec& delta,
                                                                           const vector<RealVec>& inducedDipole,
                                                                           vector<RealVec>& field) const
{

    // field at i due induced dipole at j

    RealOpenMM dur  = inducedDipole[jIndex].dot(delta);
    field[iIndex]  += delta*(dur*preFactor2) + inducedDipole[jIndex]*preFactor1;

    // field at j due induced dipole at i

               dur  = inducedDipole[iIndex].dot(delta);
    field[jIndex]  += delta*(dur*preFactor2) + inducedDipole[iIndex]*preFactor1;
}

void AmoebaReferencePmeMultipoleForce::calculateDirectInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                                                            const MultipoleParticleData& particleJ,
                                                                            vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    // compute the real space portion of the Ewald summation

    RealOpenMM uscale = 1.0;
    RealVec deltaR    = particleJ.position - particleI.position;

    // periodic boundary conditions

    getPeriodicDelta(deltaR);
    RealOpenMM r2     = deltaR.dot(deltaR);

    if (r2 > _cutoffDistanceSquared)
        return;

    RealOpenMM r           = SQRT(r2);

    // calculate the error function damping terms

    RealOpenMM ralpha      = _alphaEwald*r;

    RealOpenMM bn0         = erfc(ralpha)/r;
    RealOpenMM alsq2       = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM alsq2n      = 1.0/(SQRT_PI*_alphaEwald);
    RealOpenMM exp2a       = EXP(-(ralpha*ralpha));
    alsq2n                *= alsq2;
    RealOpenMM bn1         = (bn0+alsq2n*exp2a)/r2;

    alsq2n                *= alsq2;
    RealOpenMM bn2         = (3.0*bn1+alsq2n*exp2a)/r2;

    alsq2n                *= alsq2;
    RealOpenMM bn3         = (5.0*bn2+alsq2n*exp2a)/r2;

    // compute the error function scaled and unscaled terms

    RealOpenMM scale3      = 1.0;
    RealOpenMM scale5      = 1.0;
    RealOpenMM scale7      = 1.0;
    RealOpenMM damp        = particleI.dampingFactor*particleJ.dampingFactor;
    if (damp != 0.0) {

        RealOpenMM ratio  = (r/damp);
              ratio       = ratio*ratio*ratio;
        RealOpenMM pgamma = particleI.thole < particleJ.thole ? particleI.thole : particleJ.thole;
              damp        = -pgamma*ratio;

        if (damp > -50.0) {
            RealOpenMM expdamp = expf(damp);
            scale3        = 1.0 - expdamp;
            scale5        = 1.0 - expdamp*(1.0-damp);
            scale7        = 1.0 - (1.0 - damp + (0.6*damp*damp))*expdamp;
        }
    }
    RealOpenMM dsc3        = uscale*scale3;
    RealOpenMM dsc5        = uscale*scale5;
    RealOpenMM dsc7        = uscale*scale7;

    RealOpenMM r3          = (r*r2);
    RealOpenMM r5          = (r3*r2);
    RealOpenMM r7          = (r5*r2);
    RealOpenMM rr3         = (1.0-dsc3)/r3;
    RealOpenMM rr5         = 3.0*(1.0-dsc5)/r5;
    RealOpenMM rr7         = 15.0*(1.0-dsc7)/r7;

    RealOpenMM preFactor1  = rr3 - bn1;
    RealOpenMM preFactor2  = bn2 - rr5;
    RealOpenMM preFactor3  = bn3 - rr7;

    for (unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++) {
        calculateDirectInducedDipolePairIxn(particleI.particleIndex, particleJ.particleIndex, preFactor1, preFactor2, deltaR,
                                            *updateInducedDipoleFields[ii].inducedDipoles,
                                            updateInducedDipoleFields[ii].inducedDipoleField);
        if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
            // Compute and store the field gradient for later use.
            RealOpenMM dx = deltaR[0];
            RealOpenMM dy = deltaR[1];
            RealOpenMM dz = deltaR[2];

            OpenMM::RealVec &dipolesI = (*updateInducedDipoleFields[ii].inducedDipoles)[particleI.particleIndex];
            RealOpenMM xDipole = dipolesI[0];
            RealOpenMM yDipole = dipolesI[1];
            RealOpenMM zDipole = dipolesI[2];
            RealOpenMM muDotR = xDipole*dx + yDipole*dy + zDipole*dz;
            RealOpenMM Exx = muDotR*dx*dx*preFactor3 - (2.0*xDipole*dx + muDotR)*preFactor2;
            RealOpenMM Eyy = muDotR*dy*dy*preFactor3 - (2.0*yDipole*dy + muDotR)*preFactor2;
            RealOpenMM Ezz = muDotR*dz*dz*preFactor3 - (2.0*zDipole*dz + muDotR)*preFactor2;
            RealOpenMM Exy = muDotR*dx*dy*preFactor3 - (xDipole*dy + yDipole*dx)*preFactor2;
            RealOpenMM Exz = muDotR*dx*dz*preFactor3 - (xDipole*dz + zDipole*dx)*preFactor2;
            RealOpenMM Eyz = muDotR*dy*dz*preFactor3 - (yDipole*dz + zDipole*dy)*preFactor2;

            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][0] -= Exx;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][1] -= Eyy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][2] -= Ezz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][3] -= Exy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][4] -= Exz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleJ.particleIndex][5] -= Eyz;

            OpenMM::RealVec &dipolesJ = (*updateInducedDipoleFields[ii].inducedDipoles)[particleJ.particleIndex];
            xDipole = dipolesJ[0];
            yDipole = dipolesJ[1];
            zDipole = dipolesJ[2];
            muDotR = xDipole*dx + yDipole*dy + zDipole*dz;
            Exx = muDotR*dx*dx*preFactor3 - (2.0*xDipole*dx + muDotR)*preFactor2;
            Eyy = muDotR*dy*dy*preFactor3 - (2.0*yDipole*dy + muDotR)*preFactor2;
            Ezz = muDotR*dz*dz*preFactor3 - (2.0*zDipole*dz + muDotR)*preFactor2;
            Exy = muDotR*dx*dy*preFactor3 - (xDipole*dy + yDipole*dx)*preFactor2;
            Exz = muDotR*dx*dz*preFactor3 - (xDipole*dz + zDipole*dx)*preFactor2;
            Eyz = muDotR*dy*dz*preFactor3 - (yDipole*dz + zDipole*dy)*preFactor2;

            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][0] += Exx;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][1] += Eyy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][2] += Ezz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][3] += Exy;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][4] += Exz;
            updateInducedDipoleFields[ii].inducedDipoleFieldGradient[particleI.particleIndex][5] += Eyz;
        }
    }
}

RealOpenMM AmoebaReferencePmeMultipoleForce::calculatePmeSelfEnergy(const vector<MultipoleParticleData>& particleData) const
{
    RealOpenMM cii = 0.0;
    RealOpenMM dii = 0.0;
    RealOpenMM qii = 0.0;
    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        const MultipoleParticleData& particleI = particleData[ii];

        cii += particleI.charge*particleI.charge;

        RealVec dipole(particleI.sphericalDipole[1], particleI.sphericalDipole[2], particleI.sphericalDipole[0]);
        dii += dipole.dot(dipole + (_inducedDipole[ii]+_inducedDipolePolar[ii])*0.5);

        qii += (particleI.sphericalQuadrupole[0]*particleI.sphericalQuadrupole[0]
               +particleI.sphericalQuadrupole[1]*particleI.sphericalQuadrupole[1]
               +particleI.sphericalQuadrupole[2]*particleI.sphericalQuadrupole[2]
               +particleI.sphericalQuadrupole[3]*particleI.sphericalQuadrupole[3]
               +particleI.sphericalQuadrupole[4]*particleI.sphericalQuadrupole[4]);
    }
    RealOpenMM prefac = -_alphaEwald * _electric / (_dielectric*SQRT_PI);
    RealOpenMM a2 = _alphaEwald * _alphaEwald;
    RealOpenMM a4 = a2*a2;
    RealOpenMM energy = prefac*(cii + twoThirds*a2*dii + fourOverFifteen*a4*qii);
    return energy;
}

void AmoebaReferencePmeMultipoleForce::calculatePmeSelfTorque(const vector<MultipoleParticleData>& particleData,
                                                              vector<RealVec>& torques) const
{
    RealOpenMM term = (2.0/3.0)*(_electric/_dielectric)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        const MultipoleParticleData& particleI = particleData[ii];
        RealVec ui                             = (_inducedDipole[ii] + _inducedDipolePolar[ii]);
        RealVec dipole(particleI.sphericalDipole[1], particleI.sphericalDipole[2], particleI.sphericalDipole[0]);
        RealVec torque                         = dipole.cross(ui);
                torque                        *= term;
       torques[ii]                            += torque;
    }
}

RealOpenMM AmoebaReferencePmeMultipoleForce::calculatePmeDirectElectrostaticPairIxn(const MultipoleParticleData& particleI,
                                                                                    const MultipoleParticleData& particleJ,
                                                                                    const vector<RealOpenMM>& scalingFactors,
                                                                                    vector<RealVec>& forces,
                                                                                    vector<RealVec>& torques) const
{

    unsigned int iIndex = particleI.particleIndex;
    unsigned int jIndex = particleJ.particleIndex;

    RealOpenMM energy;
    RealVec deltaR   = particleJ.position - particleI.position;
    getPeriodicDelta(deltaR);
    RealOpenMM r2    = deltaR.dot(deltaR);

    if (r2 > _cutoffDistanceSquared)return 0.0;

    RealOpenMM r      = SQRT(r2);

    // Start by constructing rotation matrices to put dipoles and
    // quadrupoles into the QI frame, from the lab frame.
    RealOpenMM qiRotationMatrix1[3][3];
    formQIRotationMatrix(particleI.position, particleJ.position, deltaR, r, qiRotationMatrix1);
    RealOpenMM qiRotationMatrix2[5][5];
    buildSphericalQuadrupoleRotationMatrix(qiRotationMatrix1, qiRotationMatrix2);
    // The force rotation matrix rotates the QI forces into the lab
    // frame, and makes sure the result is in {x,y,z} ordering. Its
    // transpose is used to rotate the induced dipoles to the QI frame.
    RealOpenMM forceRotationMatrix[3][3];
    forceRotationMatrix[0][0] = qiRotationMatrix1[1][1];
    forceRotationMatrix[0][1] = qiRotationMatrix1[2][1];
    forceRotationMatrix[0][2] = qiRotationMatrix1[0][1];
    forceRotationMatrix[1][0] = qiRotationMatrix1[1][2];
    forceRotationMatrix[1][1] = qiRotationMatrix1[2][2];
    forceRotationMatrix[1][2] = qiRotationMatrix1[0][2];
    forceRotationMatrix[2][0] = qiRotationMatrix1[1][0];
    forceRotationMatrix[2][1] = qiRotationMatrix1[2][0];
    forceRotationMatrix[2][2] = qiRotationMatrix1[0][0];
    // For efficiency, we go ahead and cache that transposed version
    // now, because we need to do 4 rotations in total (I,J, and p,d).
    // We also fold in the factor of 0.5 needed to average the p and d
    // components.
    RealOpenMM inducedDipoleRotationMatrix[3][3];
    inducedDipoleRotationMatrix[0][0] = 0.5*qiRotationMatrix1[0][1];
    inducedDipoleRotationMatrix[0][1] = 0.5*qiRotationMatrix1[0][2];
    inducedDipoleRotationMatrix[0][2] = 0.5*qiRotationMatrix1[0][0];
    inducedDipoleRotationMatrix[1][0] = 0.5*qiRotationMatrix1[1][1];
    inducedDipoleRotationMatrix[1][1] = 0.5*qiRotationMatrix1[1][2];
    inducedDipoleRotationMatrix[1][2] = 0.5*qiRotationMatrix1[1][0];
    inducedDipoleRotationMatrix[2][0] = 0.5*qiRotationMatrix1[2][1];
    inducedDipoleRotationMatrix[2][1] = 0.5*qiRotationMatrix1[2][2];
    inducedDipoleRotationMatrix[2][2] = 0.5*qiRotationMatrix1[2][0];

    // Rotate the induced dipoles to the QI frame.
    RealOpenMM qiUindI[3], qiUindJ[3], qiUinpI[3], qiUinpJ[3];
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM valIP = 0.0;
        RealOpenMM valID = 0.0;
        RealOpenMM valJP = 0.0;
        RealOpenMM valJD = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            valIP += inducedDipoleRotationMatrix[ii][jj] * _inducedDipolePolar[iIndex][jj];
            valID += inducedDipoleRotationMatrix[ii][jj] * _inducedDipole[iIndex][jj];
            valJP += inducedDipoleRotationMatrix[ii][jj] * _inducedDipolePolar[jIndex][jj];
            valJD += inducedDipoleRotationMatrix[ii][jj] * _inducedDipole[jIndex][jj];
        }
        qiUindI[ii] = valID;
        qiUinpI[ii] = valIP;
        qiUindJ[ii] = valJD;
        qiUinpJ[ii] = valJP;
    }

    // The Qtilde intermediates (QI frame multipoles) for atoms I and J
    RealOpenMM qiQI[9], qiQJ[9];
    // Rotate the permanent multipoles to the QI frame.
    qiQI[0] = particleI.charge;
    qiQJ[0] = particleJ.charge;
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM valI = 0.0;
        RealOpenMM valJ = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            valI += qiRotationMatrix1[ii][jj] * particleI.sphericalDipole[jj];
            valJ += qiRotationMatrix1[ii][jj] * particleJ.sphericalDipole[jj];
        }
        qiQI[ii+1] = valI;
        qiQJ[ii+1] = valJ;
    }
    for (int ii = 0; ii < 5; ii++) {
        RealOpenMM valI = 0.0;
        RealOpenMM valJ = 0.0;
        for (int jj = 0; jj < 5; jj++) {
            valI += qiRotationMatrix2[ii][jj] * particleI.sphericalQuadrupole[jj];
            valJ += qiRotationMatrix2[ii][jj] * particleJ.sphericalQuadrupole[jj];
        }
        qiQI[ii+4] = valI;
        qiQJ[ii+4] = valJ;
    }

    // The Qtilde{x,y,z} torque intermediates for atoms I and J, which are used to obtain the torques on the permanent moments.
    RealOpenMM qiQIX[9] = {0.0, qiQI[3], 0.0, -qiQI[1], sqrtThree*qiQI[6], qiQI[8], -sqrtThree*qiQI[4] - qiQI[7], qiQI[6], -qiQI[5]};
    RealOpenMM qiQIY[9] = {0.0, -qiQI[2], qiQI[1], 0.0, -sqrtThree*qiQI[5], sqrtThree*qiQI[4] - qiQI[7], -qiQI[8], qiQI[5], qiQI[6]};
    RealOpenMM qiQIZ[9] = {0.0, 0.0, -qiQI[3], qiQI[2], 0.0, -qiQI[6], qiQI[5], -2.0*qiQI[8], 2.0*qiQI[7]};
    RealOpenMM qiQJX[9] = {0.0, qiQJ[3], 0.0, -qiQJ[1], sqrtThree*qiQJ[6], qiQJ[8], -sqrtThree*qiQJ[4] - qiQJ[7], qiQJ[6], -qiQJ[5]};
    RealOpenMM qiQJY[9] = {0.0, -qiQJ[2], qiQJ[1], 0.0, -sqrtThree*qiQJ[5], sqrtThree*qiQJ[4] - qiQJ[7], -qiQJ[8], qiQJ[5], qiQJ[6]};
    RealOpenMM qiQJZ[9] = {0.0, 0.0, -qiQJ[3], qiQJ[2], 0.0, -qiQJ[6], qiQJ[5], -2.0*qiQJ[8], 2.0*qiQJ[7]};

    // The field derivatives at I due to permanent and induced moments on J, and vice-versa.
    // Also, their derivatives w.r.t. R, which are needed for force calculations
    RealOpenMM Vij[9], Vji[9], VjiR[9], VijR[9];
    // The field derivatives at I due to only permanent moments on J, and vice-versa.
    RealOpenMM Vijp[3], Vijd[3], Vjip[3], Vjid[3];
    RealOpenMM rInvVec[7], alphaRVec[8], bVec[5];

    RealOpenMM prefac = (_electric/_dielectric);
    RealOpenMM rInv = 1.0 / r;

    // The rInvVec array is defined such that the ith element is R^-i, with the
    // dieleectric constant folded in, to avoid conversions later.
    rInvVec[1] = prefac * rInv;
    for (int i = 2; i < 7; ++i)
        rInvVec[i] = rInvVec[i-1] * rInv;

    // The alpharVec array is defined such that the ith element is (alpha R)^i,
    // where kappa (alpha in OpenMM parlance) is the Ewald attenuation parameter.
    alphaRVec[1] = _alphaEwald * r;
    for (int i = 2; i < 8; ++i)
        alphaRVec[i] = alphaRVec[i-1] * alphaRVec[1];

    RealOpenMM erfAlphaR = erf(alphaRVec[1]);
    RealOpenMM X = 2.0*EXP(-alphaRVec[2])/SQRT_PI;
    RealOpenMM mScale = scalingFactors[M_SCALE];
    RealOpenMM dScale = scalingFactors[D_SCALE];
    RealOpenMM pScale = scalingFactors[P_SCALE];
    RealOpenMM uScale = scalingFactors[U_SCALE];

    int doubleFactorial = 1, facCount = 1;
    RealOpenMM tmp = alphaRVec[1];
    bVec[1] = -erfAlphaR;
    for (int i=2; i < 5; ++i) {
        bVec[i] = bVec[i-1] + tmp * X / (RealOpenMM)(doubleFactorial);
        facCount = facCount + 2;
        doubleFactorial = doubleFactorial * facCount;
        tmp *= 2.0 * alphaRVec[2];
    }

    RealOpenMM dmp = particleI.dampingFactor*particleJ.dampingFactor;
    RealOpenMM a = particleI.thole < particleJ.thole ? particleI.thole : particleJ.thole;
    RealOpenMM u = std::abs(dmp) > 1.0E-5 ? r/dmp : 1E10;
    RealOpenMM au3 = a*u*u*u;
    RealOpenMM expau3 = au3 < 50.0 ? EXP(-au3) : 0.0;
    RealOpenMM a2u6 = au3*au3;
    RealOpenMM a3u9 = a2u6*au3;
    // Thole damping factors for energies
    RealOpenMM thole_c  = 1.0 - expau3;
    RealOpenMM thole_d0 = 1.0 - expau3*(1.0 + 1.5*au3);
    RealOpenMM thole_d1 = 1.0 - expau3;
    RealOpenMM thole_q0 = 1.0 - expau3*(1.0 + au3 + a2u6);
    RealOpenMM thole_q1 = 1.0 - expau3*(1.0 + au3);
    // Thole damping factors for derivatives
    RealOpenMM dthole_c  = 1.0 - expau3*(1.0 + 1.5*au3);
    RealOpenMM dthole_d0 = 1.0 - expau3*(1.0 + au3 + 1.5*a2u6);
    RealOpenMM dthole_d1 = 1.0 - expau3*(1.0 + au3);
    RealOpenMM dthole_q0 = 1.0 - expau3*(1.0 + au3 + 0.25*a2u6 + 0.75*a3u9);
    RealOpenMM dthole_q1 = 1.0 - expau3*(1.0 + au3 + 0.75*a2u6);

    // Now we compute the (attenuated) Coulomb operator and its derivatives, contracted with
    // permanent moments and induced dipoles.  Note that the coefficient of the permanent force
    // terms is half of the expected value; this is because we compute the interaction of I with
    // the sum of induced and permanent moments on J, as well as the interaction of J with I's
    // permanent and induced moments; doing so double counts the permanent-permanent interaction.
    RealOpenMM ePermCoef, dPermCoef, eUIndCoef, dUIndCoef, eUInpCoef, dUInpCoef;

    // C-C terms (m=0)
    ePermCoef = rInvVec[1]*(mScale + bVec[2] - alphaRVec[1]*X);
    dPermCoef = -0.5*(mScale + bVec[2])*rInvVec[2];
    Vij[0]  = ePermCoef*qiQJ[0];
    Vji[0]  = ePermCoef*qiQI[0];
    VijR[0] = dPermCoef*qiQJ[0];
    VjiR[0] = dPermCoef*qiQI[0];

    // C-D and C-Uind terms (m=0)
    ePermCoef = rInvVec[2]*(mScale + bVec[2]);
    eUIndCoef = rInvVec[2]*(pScale*thole_c + bVec[2]);
    eUInpCoef = rInvVec[2]*(dScale*thole_c + bVec[2]);
    dPermCoef = -rInvVec[3]*(mScale + bVec[2] + alphaRVec[3]*X);
    dUIndCoef = -2.0*rInvVec[3]*(pScale*dthole_c + bVec[2] + alphaRVec[3]*X);
    dUInpCoef = -2.0*rInvVec[3]*(dScale*dthole_c + bVec[2] + alphaRVec[3]*X);
    Vij[0]  += -(ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0]);
    Vji[1]   = -(ePermCoef*qiQI[0]);
    VijR[0] += -(dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0]);
    VjiR[1]  = -(dPermCoef*qiQI[0]);
    Vjip[0]  = -(eUInpCoef*qiQI[0]);
    Vjid[0]  = -(eUIndCoef*qiQI[0]);
    // D-C and Uind-C terms (m=0)
    Vij[1]   = ePermCoef*qiQJ[0];
    Vji[0]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1]  = dPermCoef*qiQJ[0];
    VjiR[0] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0]  = eUInpCoef*qiQJ[0];
    Vijd[0]  = eUIndCoef*qiQJ[0];

    // D-D and D-Uind terms (m=0)
    ePermCoef = -twoThirds*rInvVec[3]*(3.0*(mScale + bVec[3]) + alphaRVec[3]*X);
    eUIndCoef = -twoThirds*rInvVec[3]*(3.0*(pScale*thole_d0 + bVec[3]) + alphaRVec[3]*X);
    eUInpCoef = -twoThirds*rInvVec[3]*(3.0*(dScale*thole_d0 + bVec[3]) + alphaRVec[3]*X);
    dPermCoef = rInvVec[4]*(3.0*(mScale + bVec[3]) + 2.*alphaRVec[5]*X);
    dUIndCoef = rInvVec[4]*(6.0*(pScale*dthole_d0 + bVec[3]) + 4.0*alphaRVec[5]*X);
    dUInpCoef = rInvVec[4]*(6.0*(dScale*dthole_d0 + bVec[3]) + 4.0*alphaRVec[5]*X);
    Vij[1]  += ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0];
    Vji[1]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1] += dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0];
    VjiR[1] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0] += eUInpCoef*qiQJ[1];
    Vijd[0] += eUIndCoef*qiQJ[1];
    Vjip[0] += eUInpCoef*qiQI[1];
    Vjid[0] += eUIndCoef*qiQI[1];
    // D-D and D-Uind terms (m=1)
    ePermCoef = rInvVec[3]*(mScale + bVec[3] - twoThirds*alphaRVec[3]*X);
    eUIndCoef = rInvVec[3]*(pScale*thole_d1 + bVec[3] - twoThirds*alphaRVec[3]*X);
    eUInpCoef = rInvVec[3]*(dScale*thole_d1 + bVec[3] - twoThirds*alphaRVec[3]*X);
    dPermCoef = -1.5*rInvVec[4]*(mScale + bVec[3]);
    dUIndCoef = -3.0*rInvVec[4]*(pScale*dthole_d1 + bVec[3]);
    dUInpCoef = -3.0*rInvVec[4]*(dScale*dthole_d1 + bVec[3]);
    Vij[2]  = ePermCoef*qiQJ[2] + eUIndCoef*qiUindJ[1] + eUInpCoef*qiUinpJ[1];
    Vji[2]  = ePermCoef*qiQI[2] + eUIndCoef*qiUindI[1] + eUInpCoef*qiUinpI[1];
    VijR[2] = dPermCoef*qiQJ[2] + dUIndCoef*qiUindJ[1] + dUInpCoef*qiUinpJ[1];
    VjiR[2] = dPermCoef*qiQI[2] + dUIndCoef*qiUindI[1] + dUInpCoef*qiUinpI[1];
    Vij[3]  = ePermCoef*qiQJ[3] + eUIndCoef*qiUindJ[2] + eUInpCoef*qiUinpJ[2];
    Vji[3]  = ePermCoef*qiQI[3] + eUIndCoef*qiUindI[2] + eUInpCoef*qiUinpI[2];
    VijR[3] = dPermCoef*qiQJ[3] + dUIndCoef*qiUindJ[2] + dUInpCoef*qiUinpJ[2];
    VjiR[3] = dPermCoef*qiQI[3] + dUIndCoef*qiUindI[2] + dUInpCoef*qiUinpI[2];
    Vijp[1] = eUInpCoef*qiQJ[2];
    Vijd[1] = eUIndCoef*qiQJ[2];
    Vjip[1] = eUInpCoef*qiQI[2];
    Vjid[1] = eUIndCoef*qiQI[2];
    Vijp[2] = eUInpCoef*qiQJ[3];
    Vijd[2] = eUIndCoef*qiQJ[3];
    Vjip[2] = eUInpCoef*qiQI[3];
    Vjid[2] = eUIndCoef*qiQI[3];

    // C-Q terms (m=0)
    ePermCoef = (mScale + bVec[3])*rInvVec[3];
    dPermCoef = -oneThird*rInvVec[4]*(4.5*(mScale + bVec[3]) + 2.0*alphaRVec[5]*X);
    Vij[0]  += ePermCoef*qiQJ[4];
    Vji[4]   = ePermCoef*qiQI[0];
    VijR[0] += dPermCoef*qiQJ[4];
    VjiR[4]  = dPermCoef*qiQI[0];
    // Q-C terms (m=0)
    Vij[4]   = ePermCoef*qiQJ[0];
    Vji[0]  += ePermCoef*qiQI[4];
    VijR[4]  = dPermCoef*qiQJ[0];
    VjiR[0] += dPermCoef*qiQI[4];

    // D-Q and Uind-Q terms (m=0)
    ePermCoef = rInvVec[4]*(3.0*(mScale + bVec[3]) + fourThirds*alphaRVec[5]*X);
    eUIndCoef = rInvVec[4]*(3.0*(pScale*thole_q0 + bVec[3]) + fourThirds*alphaRVec[5]*X);
    eUInpCoef = rInvVec[4]*(3.0*(dScale*thole_q0 + bVec[3]) + fourThirds*alphaRVec[5]*X);
    dPermCoef = -fourThirds*rInvVec[5]*(4.5*(mScale + bVec[3]) + (1.0 + alphaRVec[2])*alphaRVec[5]*X);
    dUIndCoef = -fourThirds*rInvVec[5]*(9.0*(pScale*dthole_q0 + bVec[3]) + 2.0*(1.0 + alphaRVec[2])*alphaRVec[5]*X);
    dUInpCoef = -fourThirds*rInvVec[5]*(9.0*(dScale*dthole_q0 + bVec[3]) + 2.0*(1.0 + alphaRVec[2])*alphaRVec[5]*X);
    Vij[1]  += ePermCoef*qiQJ[4];
    Vji[4]  += ePermCoef*qiQI[1] + eUIndCoef*qiUindI[0] + eUInpCoef*qiUinpI[0];
    VijR[1] += dPermCoef*qiQJ[4];
    VjiR[4] += dPermCoef*qiQI[1] + dUIndCoef*qiUindI[0] + dUInpCoef*qiUinpI[0];
    Vijp[0] += eUInpCoef*qiQJ[4];
    Vijd[0] += eUIndCoef*qiQJ[4];
    // Q-D and Q-Uind terms (m=0)
    Vij[4]  += -(ePermCoef*qiQJ[1] + eUIndCoef*qiUindJ[0] + eUInpCoef*qiUinpJ[0]);
    Vji[1]  += -(ePermCoef*qiQI[4]);
    VijR[4] += -(dPermCoef*qiQJ[1] + dUIndCoef*qiUindJ[0] + dUInpCoef*qiUinpJ[0]);
    VjiR[1] += -(dPermCoef*qiQI[4]);
    Vjip[0] += -(eUInpCoef*qiQI[4]);
    Vjid[0] += -(eUIndCoef*qiQI[4]);

    // D-Q and Uind-Q terms (m=1)
    ePermCoef = -sqrtThree*rInvVec[4]*(mScale + bVec[3]);
    eUIndCoef = -sqrtThree*rInvVec[4]*(pScale*thole_q1 + bVec[3]);
    eUInpCoef = -sqrtThree*rInvVec[4]*(dScale*thole_q1 + bVec[3]);
    dPermCoef = fourSqrtOneThird*rInvVec[5]*(1.5*(mScale + bVec[3]) + 0.5*alphaRVec[5]*X);
    dUIndCoef = fourSqrtOneThird*rInvVec[5]*(3.0*(pScale*dthole_q1 + bVec[3]) + alphaRVec[5]*X);
    dUInpCoef = fourSqrtOneThird*rInvVec[5]*(3.0*(dScale*dthole_q1 + bVec[3]) + alphaRVec[5]*X);
    Vij[2]  += ePermCoef*qiQJ[5];
    Vji[5]   = ePermCoef*qiQI[2] + eUIndCoef*qiUindI[1] + eUInpCoef*qiUinpI[1];
    VijR[2] += dPermCoef*qiQJ[5];
    VjiR[5]  = dPermCoef*qiQI[2] + dUIndCoef*qiUindI[1] + dUInpCoef*qiUinpI[1];
    Vij[3]  += ePermCoef*qiQJ[6];
    Vji[6]   = ePermCoef*qiQI[3] + eUIndCoef*qiUindI[2] + eUInpCoef*qiUinpI[2];
    VijR[3] += dPermCoef*qiQJ[6];
    VjiR[6]  = dPermCoef*qiQI[3] + dUIndCoef*qiUindI[2] + dUInpCoef*qiUinpI[2];
    Vijp[1] += eUInpCoef*qiQJ[5];
    Vijd[1] += eUIndCoef*qiQJ[5];
    Vijp[2] += eUInpCoef*qiQJ[6];
    Vijd[2] += eUIndCoef*qiQJ[6];
    // D-Q and Uind-Q terms (m=1)
    Vij[5]   = -(ePermCoef*qiQJ[2] + eUIndCoef*qiUindJ[1] + eUInpCoef*qiUinpJ[1]);
    Vji[2]  += -(ePermCoef*qiQI[5]);
    VijR[5]  = -(dPermCoef*qiQJ[2] + dUIndCoef*qiUindJ[1] + dUInpCoef*qiUinpJ[1]);
    VjiR[2] += -(dPermCoef*qiQI[5]);
    Vij[6]   = -(ePermCoef*qiQJ[3] + eUIndCoef*qiUindJ[2] + eUInpCoef*qiUinpJ[2]);
    Vji[3]  += -(ePermCoef*qiQI[6]);
    VijR[6]  = -(dPermCoef*qiQJ[3] + dUIndCoef*qiUindJ[2] + dUInpCoef*qiUinpJ[2]);
    VjiR[3] += -(dPermCoef*qiQI[6]);
    Vjip[1] += -(eUInpCoef*qiQI[5]);
    Vjid[1] += -(eUIndCoef*qiQI[5]);
    Vjip[2] += -(eUInpCoef*qiQI[6]);
    Vjid[2] += -(eUIndCoef*qiQI[6]);

    // Q-Q terms (m=0)
    ePermCoef = rInvVec[5]*(6.0*(mScale + bVec[4]) + fourOverFortyFive*(-3.0 + 10.0*alphaRVec[2])*alphaRVec[5]*X);
    dPermCoef = -oneNinth*rInvVec[6]*(135.0*(mScale + bVec[4]) + 4.0*(1.0 + 2.0*alphaRVec[2])*alphaRVec[7]*X);
    Vij[4]  += ePermCoef*qiQJ[4];
    Vji[4]  += ePermCoef*qiQI[4];
    VijR[4] += dPermCoef*qiQJ[4];
    VjiR[4] += dPermCoef*qiQI[4];
    // Q-Q terms (m=1)
    ePermCoef = -fourOverFifteen*rInvVec[5]*(15.0*(mScale + bVec[4]) + alphaRVec[5]*X);
    dPermCoef = rInvVec[6]*(10.0*(mScale + bVec[4]) + fourThirds*alphaRVec[7]*X);
    Vij[5]  += ePermCoef*qiQJ[5];
    Vji[5]  += ePermCoef*qiQI[5];
    VijR[5] += dPermCoef*qiQJ[5];
    VjiR[5] += dPermCoef*qiQI[5];
    Vij[6]  += ePermCoef*qiQJ[6];
    Vji[6]  += ePermCoef*qiQI[6];
    VijR[6] += dPermCoef*qiQJ[6];
    VjiR[6] += dPermCoef*qiQI[6];
    // Q-Q terms (m=2)
    ePermCoef = rInvVec[5]*(mScale + bVec[4] - fourOverFifteen*alphaRVec[5]*X);
    dPermCoef = -2.5*(mScale + bVec[4])*rInvVec[6];
    Vij[7]  = ePermCoef*qiQJ[7];
    Vji[7]  = ePermCoef*qiQI[7];
    VijR[7] = dPermCoef*qiQJ[7];
    VjiR[7] = dPermCoef*qiQI[7];
    Vij[8]  = ePermCoef*qiQJ[8];
    Vji[8]  = ePermCoef*qiQI[8];
    VijR[8] = dPermCoef*qiQJ[8];
    VjiR[8] = dPermCoef*qiQI[8];

    // Evaluate the energies, forces and torques due to permanent+induced moments
    // interacting with just the permanent moments.
    energy = 0.5*(qiQI[0]*Vij[0] + qiQJ[0]*Vji[0]);
    RealOpenMM fIZ = qiQI[0]*VijR[0];
    RealOpenMM fJZ = qiQJ[0]*VjiR[0];
    RealOpenMM EIX = 0.0, EIY = 0.0, EIZ = 0.0, EJX = 0.0, EJY = 0.0, EJZ = 0.0;
    for (int i = 1; i < 9; ++i) {
        energy += 0.5*(qiQI[i]*Vij[i] + qiQJ[i]*Vji[i]);
        fIZ += qiQI[i]*VijR[i];
        fJZ += qiQJ[i]*VjiR[i];
        EIX += qiQIX[i]*Vij[i];
        EIY += qiQIY[i]*Vij[i];
        EIZ += qiQIZ[i]*Vij[i];
        EJX += qiQJX[i]*Vji[i];
        EJY += qiQJY[i]*Vji[i];
        EJZ += qiQJZ[i]*Vji[i];
    }

    // Define the torque intermediates for the induced dipoles. These are simply the induced dipole torque
    // intermediates dotted with the field due to permanent moments only, at each center. We inline the
    // induced dipole torque intermediates here, for simplicity. N.B. There are no torques on the dipoles
    // themselves, so we accumulate the torque intermediates into separate variables to allow them to be
    // used only in the force calculation.
    //
    // The torque about the x axis (needed to obtain the y force on the induced dipoles, below)
    //    qiUindIx[0] = qiQUindI[2];    qiUindIx[1] = 0;    qiUindIx[2] = -qiQUindI[0]
    RealOpenMM iEIX = qiUinpI[2]*Vijp[0] + qiUindI[2]*Vijd[0] - qiUinpI[0]*Vijp[2] - qiUindI[0]*Vijd[2];
    RealOpenMM iEJX = qiUinpJ[2]*Vjip[0] + qiUindJ[2]*Vjid[0] - qiUinpJ[0]*Vjip[2] - qiUindJ[0]*Vjid[2];
    // The torque about the y axis (needed to obtain the x force on the induced dipoles, below)
    //    qiUindIy[0] = -qiQUindI[1];   qiUindIy[1] = qiQUindI[0];    qiUindIy[2] = 0
    RealOpenMM iEIY = qiUinpI[0]*Vijp[1] + qiUindI[0]*Vijd[1] - qiUinpI[1]*Vijp[0] - qiUindI[1]*Vijd[0];
    RealOpenMM iEJY = qiUinpJ[0]*Vjip[1] + qiUindJ[0]*Vjid[1] - qiUinpJ[1]*Vjip[0] - qiUindJ[1]*Vjid[0];

    // Add in the induced-induced terms, if needed.
    if(getPolarizationType() == AmoebaReferenceMultipoleForce::Mutual) {
        // Uind-Uind terms (m=0)
        RealOpenMM eCoef = -fourThirds*rInvVec[3]*(3.0*(uScale*thole_d0 + bVec[3]) + alphaRVec[3]*X);
        RealOpenMM dCoef = rInvVec[4]*(6.0*(uScale*dthole_d0 + bVec[3]) + 4.0*alphaRVec[5]*X);
        iEIX += eCoef*(qiUinpI[2]*qiUindJ[0] + qiUindI[2]*qiUinpJ[0]);
        iEJX += eCoef*(qiUinpJ[2]*qiUindI[0] + qiUindJ[2]*qiUinpI[0]);
        iEIY -= eCoef*(qiUinpI[1]*qiUindJ[0] + qiUindI[1]*qiUinpJ[0]);
        iEJY -= eCoef*(qiUinpJ[1]*qiUindI[0] + qiUindJ[1]*qiUinpI[0]);
        fIZ += dCoef*(qiUinpI[0]*qiUindJ[0] + qiUindI[0]*qiUinpJ[0]);
        fIZ += dCoef*(qiUinpJ[0]*qiUindI[0] + qiUindJ[0]*qiUinpI[0]);
        // Uind-Uind terms (m=1)
        eCoef = 2.0*rInvVec[3]*(uScale*thole_d1 + bVec[3] - twoThirds*alphaRVec[3]*X);
        dCoef = -3.0*rInvVec[4]*(uScale*dthole_d1 + bVec[3]);
        iEIX -= eCoef*(qiUinpI[0]*qiUindJ[2] + qiUindI[0]*qiUinpJ[2]);
        iEJX -= eCoef*(qiUinpJ[0]*qiUindI[2] + qiUindJ[0]*qiUinpI[2]);
        iEIY += eCoef*(qiUinpI[0]*qiUindJ[1] + qiUindI[0]*qiUinpJ[1]);
        iEJY += eCoef*(qiUinpJ[0]*qiUindI[1] + qiUindJ[0]*qiUinpI[1]);
        fIZ += dCoef*(qiUinpI[1]*qiUindJ[1] + qiUindI[1]*qiUinpJ[1] + qiUinpI[2]*qiUindJ[2] + qiUindI[2]*qiUinpJ[2]);
        fIZ += dCoef*(qiUinpJ[1]*qiUindI[1] + qiUindJ[1]*qiUinpI[1] + qiUinpJ[2]*qiUindI[2] + qiUindJ[2]*qiUinpI[2]);
    }

    // The quasi-internal frame forces and torques.  Note that the induced torque intermediates are
    // used in the force expression, but not in the torques; the induced dipoles are isotropic.
    RealOpenMM qiForce[3] = {rInv*(EIY+EJY+iEIY+iEJY), -rInv*(EIX+EJX+iEIX+iEJX), -(fJZ+fIZ)};
    RealOpenMM qiTorqueI[3] = {-EIX, -EIY, -EIZ};
    RealOpenMM qiTorqueJ[3] = {-EJX, -EJY, -EJZ};

    // Rotate the forces and torques back to the lab frame
    for (int ii = 0; ii < 3; ii++) {
        RealOpenMM forceVal = 0.0;
        RealOpenMM torqueIVal = 0.0;
        RealOpenMM torqueJVal = 0.0;
        for (int jj = 0; jj < 3; jj++) {
            forceVal   += forceRotationMatrix[ii][jj] * qiForce[jj];
            torqueIVal += forceRotationMatrix[ii][jj] * qiTorqueI[jj];
            torqueJVal += forceRotationMatrix[ii][jj] * qiTorqueJ[jj];
        }
        torques[iIndex][ii] += torqueIVal;
        torques[jIndex][ii] += torqueJVal;
        forces[iIndex][ii]  -= forceVal;
        forces[jIndex][ii]  += forceVal;
    }
    return energy;

}

RealOpenMM AmoebaReferencePmeMultipoleForce::calculateElectrostatic(const vector<MultipoleParticleData>& particleData,
                                                                    vector<RealVec>& torques, vector<RealVec>& forces)
{
    RealOpenMM energy = 0.0;
    vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for (unsigned int kk = 0; kk < scaleFactors.size(); kk++) {
        scaleFactors[kk] = 1.0;
    }

    // loop over particle pairs for direct space interactions

    for (unsigned int ii = 0; ii < particleData.size(); ii++) {
        for (unsigned int jj = ii+1; jj < particleData.size(); jj++) {

            if (jj <= _maxScaleIndex[ii]) {
                getMultipoleScaleFactors(ii, jj, scaleFactors);
            }

            energy += calculatePmeDirectElectrostaticPairIxn(particleData[ii], particleData[jj], scaleFactors, forces, torques);

            if (jj <= _maxScaleIndex[ii]) {
                for (unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++) {
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }

    // The polarization energy
    calculatePmeSelfTorque(particleData, torques);
    energy += computeReciprocalSpaceInducedDipoleForceAndEnergy(getPolarizationType(), particleData, forces, torques);
    energy += computeReciprocalSpaceFixedMultipoleForceAndEnergy(particleData, forces, torques);
    energy += calculatePmeSelfEnergy(particleData);

    // Now that both the direct and reciprocal space contributions have been added, we can compute the dipole
    // response contributions to the forces, if we're using the extrapolated polarization algorithm.
    if (getPolarizationType() == AmoebaReferenceMultipoleForce::Extrapolated) {
        RealOpenMM prefac = (_electric/_dielectric);
        for (int i = 0; i < _numParticles; i++) {
            // Compute the µ(m) T µ(n) force contributions here
            for (int l = 0; l < _maxPTOrder-1; ++l) {
                for (int m = 0; m < _maxPTOrder-1-l; ++m) {
                    RealOpenMM p = _extPartCoefficients[l+m+1];
                    if(std::fabs(p) < 1e-6) continue;
                    forces[i][0] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+0]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+3]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+1]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipoleD[l][i][0]*_ptDipoleFieldGradientP[m][6*i+4]
                                                + _ptDipoleD[l][i][1]*_ptDipoleFieldGradientP[m][6*i+5]
                                                + _ptDipoleD[l][i][2]*_ptDipoleFieldGradientP[m][6*i+2]);
                    forces[i][0] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+0]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+4]);
                    forces[i][1] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+3]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+1]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+5]);
                    forces[i][2] += 0.5*p*prefac*(_ptDipoleP[l][i][0]*_ptDipoleFieldGradientD[m][6*i+4]
                                                + _ptDipoleP[l][i][1]*_ptDipoleFieldGradientD[m][6*i+5]
                                                + _ptDipoleP[l][i][2]*_ptDipoleFieldGradientD[m][6*i+2]);
                }
            }
        }
    }
    return energy;
}
