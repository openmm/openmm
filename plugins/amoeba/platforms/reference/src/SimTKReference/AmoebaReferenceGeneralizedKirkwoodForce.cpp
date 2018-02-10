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

#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
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
                                                                                      _surfaceAreaFactor(0.0054) {

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

void AmoebaReferenceGeneralizedKirkwoodForce::getGrycukBornRadii(vector<double>& bornRadii) const {
    bornRadii.resize(_bornRadii.size());
    copy(_bornRadii.begin(), _bornRadii.end(), bornRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::calculateGrycukBornRadii(const vector<Vec3>& particlePositions) {

    const double bigRadius = 1000.0;

    _bornRadii.resize(_numParticles);
    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        if (_atomicRadii[ii] <= 0.0) {
            _bornRadii[ii] = bigRadius;
            continue;
        }

        double bornSum = 0.0;
        for (unsigned int jj = 0; jj < _numParticles; jj++) {

            if (ii == jj || _atomicRadii[jj] < 0.0)continue;
          
            double xr       = particlePositions[jj][0] - particlePositions[ii][0];
            double yr       = particlePositions[jj][1] - particlePositions[ii][1];
            double zr       = particlePositions[jj][2] - particlePositions[ii][2];

            double r2       = xr*xr + yr*yr + zr*zr;
            double r        = sqrt(r2);

            double sk       = _atomicRadii[jj]*_scaleFactors[jj];

            // If atom ii engulfs the descreening atom, then continue.
            if (_atomicRadii[ii] > r + sk) continue;

            double sk2      = sk*sk;

            if ((_atomicRadii[ii] + r) < sk) {
                double lik       = _atomicRadii[ii];
                double uik       = sk - r;  
                double lik3      = lik*lik*lik;
                double uik3      = uik*uik*uik;
                bornSum             -= (1.0/uik3 - 1.0/lik3);
            }   
        
            double uik = r + sk; 
            double lik;
            if ((_atomicRadii[ii] + r) < sk) {
                lik = sk - r;  
            } else if (r < (_atomicRadii[ii] + sk)) {
                lik = _atomicRadii[ii];
            } else {
                lik = r - sk; 
            }   
        
            double l2          = lik*lik; 
            double l4          = l2*l2;
            double lr          = lik*r;
            double l4r         = l4*r;
        
            double u2          = uik*uik;
            double u4          = u2*u2;
            double ur          = uik*r;
            double u4r         = u4*r;
        
            double term        = (3.0*(r2-sk2) + 6.0*u2 - 8.0*ur)/u4r - (3.0*(r2-sk2) + 6.0*l2 - 8.0*lr)/l4r;
            bornSum           += term/16.0;
        
        }
        bornSum        = 1.0/(_atomicRadii[ii]*_atomicRadii[ii]*_atomicRadii[ii]) - bornSum;
        _bornRadii[ii] = (bornSum <= 0.0) ? bigRadius : pow(bornSum, -1.0/3.0);
    }

    return;
}
