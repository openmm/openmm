/* Portions copyright (c) 2006 Stanford University and Simbios.
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include <cmath>

using std::vector;
using namespace OpenMM;

AmoebaReferenceGeneralizedKirkwoodForce::AmoebaReferenceGeneralizedKirkwoodForce() : _numParticles(0),
                                                                                     _includeCavityTerm(1),
                                                                                     _directPolarization(0),
                                                                                     _soluteDielectric(1.0),
                                                                                     _solventDielectric(78.3),
                                                                                     _dielectricOffset(0.009),
                                                                                     _probeRadius(0.14),
                                                                                     _surfaceAreaFactor(0.0054),
                                                                                     _tanhRescaling(false),
                                                                                     _beta0(0.9563),
                                                                                     _beta1(0.2578),
                                                                                     _beta2(0.0810),
                                                                                     _descreenOffset(0.0) {
}

void AmoebaReferenceGeneralizedKirkwoodForce::setNumParticles(int numParticles) {
    _numParticles = numParticles;
}

int AmoebaReferenceGeneralizedKirkwoodForce::getNumParticles() const {
    return _numParticles;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setIncludeCavityTerm(int includeCavityTerm) {
    _includeCavityTerm = includeCavityTerm;
}

int AmoebaReferenceGeneralizedKirkwoodForce::getIncludeCavityTerm() const {
    return _includeCavityTerm;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setTanhRescaling(bool tanhRescaling) {
    _tanhRescaling = tanhRescaling;
}

bool AmoebaReferenceGeneralizedKirkwoodForce::getTanhRescaling() const {
    return _tanhRescaling;
}

/**
 * Get Tanh parameters beta0, beta1 and beta2.
 */
void AmoebaReferenceGeneralizedKirkwoodForce::getTanhParameters(double& b0, double& b1, double& b2) const {
    b0 = _beta0;
    b1 = _beta1;
    b2 = _beta2;
}

/**
 * Set the flag signaling whether the solute integral is rescaled by a Tanh function
 * to account for interstitial spaces.
*/
void AmoebaReferenceGeneralizedKirkwoodForce::setTanhParameters(double b0, double b1, double b2) {
    _beta0 = b0;
    _beta1 = b1;
    _beta2 = b2;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setDescreenOffset(double descreenOffset) {
    _descreenOffset  = descreenOffset;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getDescreenOffset() const {
    return _descreenOffset;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setDirectPolarization(int directPolarization) {
    _directPolarization = directPolarization;
}

int AmoebaReferenceGeneralizedKirkwoodForce::getDirectPolarization() const {
    return _directPolarization;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setSoluteDielectric(double soluteDielectric) {
    _soluteDielectric  = soluteDielectric;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getSoluteDielectric() const {
    return _soluteDielectric;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setSolventDielectric(double solventDielectric) {
    _solventDielectric  = solventDielectric;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getSolventDielectric() const {
    return _solventDielectric;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setDielectricOffset(double dielectricOffset) {
    _dielectricOffset  = dielectricOffset;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getDielectricOffset() const {
    return _dielectricOffset;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setProbeRadius(double probeRadius) {
    _probeRadius  = probeRadius;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getProbeRadius() const {
    return _probeRadius;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setSurfaceAreaFactor(double surfaceAreaFactor) {
    _surfaceAreaFactor  = surfaceAreaFactor;
}

double AmoebaReferenceGeneralizedKirkwoodForce::getSurfaceAreaFactor() const {
    return _surfaceAreaFactor;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setAtomicRadii(const vector<double>& atomicRadii) {
    _atomicRadii.resize(atomicRadii.size());
    copy(atomicRadii.begin(), atomicRadii.end(), _atomicRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getAtomicRadii(vector<double>& atomicRadii) const {
    atomicRadii.resize(_atomicRadii.size());
    copy(_atomicRadii.begin(), _atomicRadii.end(), atomicRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setScaleFactors(const vector<double>& scaleFactors) {
    _scaleFactors.resize(scaleFactors.size());
    copy(scaleFactors.begin(), scaleFactors.end(), _scaleFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getScaleFactors(vector<double>& scaleFactors) const {
    scaleFactors.resize(_scaleFactors.size());
    copy(_scaleFactors.begin(), _scaleFactors.end(), scaleFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setCharges(const vector<double>& charges) {
    _charges.resize(charges.size());
    copy(charges.begin(), charges.end(), _charges.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setDescreenRadii(const vector<double> &descreenRadii) {
    _descreenRadii.resize(descreenRadii.size());
    copy(descreenRadii.begin(), descreenRadii.end(), _descreenRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getDescreenRadii(vector<double> &descreenRadii) const {
    descreenRadii.resize(_atomicRadii.size());
    copy(_descreenRadii.begin(), _descreenRadii.end(), descreenRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setNeckFactors(const vector<double> &neckFactors) {
    _neckFactors.resize(neckFactors.size());
    copy(neckFactors.begin(), neckFactors.end(), _neckFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getNeckFactors(vector<double> &neckFactors) const {
    neckFactors.resize(_neckFactors.size());
    copy(_neckFactors.begin(), _neckFactors.end(), neckFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getGrycukBornRadii(vector<double>& bornRadii) const {
    bornRadii.resize(_bornRadii.size());
    copy(_bornRadii.begin(), _bornRadii.end(), bornRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getSoluteIntegral(vector<double> &soluteIntegral) const {
    soluteIntegral.resize(_soluteIntegral.size());
    copy(_soluteIntegral.begin(), _soluteIntegral.end(), soluteIntegral.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::calculateGrycukBornRadii(const vector<Vec3> &particlePositions) {

    // Set the radius to 30 Angstroms (3 nm) if either the base radius is zero, or the
    // descreening integral is negative.
    const double MAX_RADIUS = 3.0;
    const double RECIP_MAX_RADIUS3 = pow(MAX_RADIUS, -3.0);
    const double PI4_3 = 4.0 / 3.0 * M_PI;
    const double INVERSE_PI4_3 = 1.0 / PI4_3;
    const double ONE_THIRD = 1.0 / 3.0;

    _bornRadii.resize(_numParticles);
    _soluteIntegral.resize(_numParticles);

    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        if (_atomicRadii[ii] <= 0.0) {
            _bornRadii[ii] = MAX_RADIUS;
            continue;
        }

        double integralStartI = max(_atomicRadii[ii], _descreenRadii[ii]) + _descreenOffset;

        double bornSum = 0.0;
        double neckSum = 0.0;
        for (unsigned int jj = 0; jj < _numParticles; jj++) {

            double sk = _descreenRadii[jj] * _scaleFactors[jj];

            if (ii == jj || integralStartI <= 0.0 || sk <= 0.0) continue;

            double xr = particlePositions[jj][0] - particlePositions[ii][0];
            double yr = particlePositions[jj][1] - particlePositions[ii][1];
            double zr = particlePositions[jj][2] - particlePositions[ii][2];

            double r2 = xr * xr + yr * yr + zr * zr;
            double r = sqrt(r2);

            // If atom ii engulfs the descreening atom, then continue.
            if (integralStartI > r + sk) continue;

            double sk2 = sk * sk;

            if ((integralStartI + r) < sk) {
                double lik = integralStartI;
                double uik = sk - r;
                double lik3 = lik * lik * lik;
                double uik3 = uik * uik * uik;
                bornSum -= (1.0 / uik3 - 1.0 / lik3);
            }

            double uik = r + sk;
            double lik;
            if ((integralStartI + r) < sk)
                lik = sk - r;
            else if (r < (integralStartI + sk))
                lik = integralStartI;
            else
                lik = r - sk;

            double l2 = lik * lik;
            double l4 = l2 * l2;
            double lr = lik * r;
            double l4r = l4 * r;

            double u2 = uik * uik;
            double u4 = u2 * u2;
            double ur = uik * r;
            double u4r = u4 * r;

            double term =
                    (3.0 * (r2 - sk2) + 6.0 * u2 - 8.0 * ur) / u4r
                    - (3.0 * (r2 - sk2) + 6.0 * l2 - 8.0 * lr) / l4r;
            bornSum += term / 16.0;

            double mixedNeckScale = 0.5 * (_neckFactors[ii] + _neckFactors[jj]);
            if (mixedNeckScale > 0.0) {
                neckSum += AmoebaGeneralizedKirkwoodForceImpl::neckDescreen(r, integralStartI, _descreenRadii[jj], mixedNeckScale);
            }
        }

        bornSum = bornSum * PI4_3 + neckSum;
        _soluteIntegral[ii] = bornSum;

        double baseRadiusI3 = _atomicRadii[ii] * _atomicRadii[ii] * _atomicRadii[ii];

        if (_tanhRescaling) {
            // Set up tanh function components
            double rhoi3Psi = baseRadiusI3 * _soluteIntegral[ii];
            double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
            double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
            // If the output of the tanh function is 1.0, then the Born radius will be MaxBornRadius
            double tanh_constant = PI4_3 * ((1.0 / baseRadiusI3) - RECIP_MAX_RADIUS3);
            bornSum = tanh_constant * tanh(_beta0 * rhoi3Psi - _beta1 * rhoi6Psi2 + _beta2 * rhoi9Psi3);
        }

        bornSum = PI4_3 / baseRadiusI3 - bornSum;

        _bornRadii[ii] = (bornSum <= 0.0) ? MAX_RADIUS : pow(INVERSE_PI4_3 * bornSum, -ONE_THIRD);

        // Born radius must be at least as large as its base radius.
        if (_bornRadii[ii] < _atomicRadii[ii]) {
            _bornRadii[ii] = _atomicRadii[ii];
        }

        // Maximum Born radius is 30.0 Angstroms.
        if (_bornRadii[ii] > MAX_RADIUS) {
            _bornRadii[ii] = MAX_RADIUS;
        }
    }

}
