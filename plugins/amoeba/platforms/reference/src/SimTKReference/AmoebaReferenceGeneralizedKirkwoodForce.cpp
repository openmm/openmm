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

void AmoebaReferenceGeneralizedKirkwoodForce::setSoluteDielectric(RealOpenMM soluteDielectric) {
    _soluteDielectric  = soluteDielectric;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodForce::getSoluteDielectric() const {
    return _soluteDielectric;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setSolventDielectric(RealOpenMM solventDielectric) {
    _solventDielectric  = solventDielectric;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodForce::getSolventDielectric() const {
    return _solventDielectric;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setDielectricOffset(RealOpenMM dielectricOffset) {
    _dielectricOffset  = dielectricOffset;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodForce::getDielectricOffset() const {
    return _dielectricOffset;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setProbeRadius(RealOpenMM probeRadius) {
    _probeRadius  = probeRadius;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodForce::getProbeRadius() const {
    return _probeRadius;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setSurfaceAreaFactor(RealOpenMM surfaceAreaFactor) {
    _surfaceAreaFactor  = surfaceAreaFactor;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodForce::getSurfaceAreaFactor() const {
    return _surfaceAreaFactor;
}

void AmoebaReferenceGeneralizedKirkwoodForce::setAtomicRadii(const vector<RealOpenMM>& atomicRadii) {
    _atomicRadii.resize(atomicRadii.size());
    copy(atomicRadii.begin(), atomicRadii.end(), _atomicRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getAtomicRadii(vector<RealOpenMM>& atomicRadii) const {
    atomicRadii.resize(_atomicRadii.size());
    copy(_atomicRadii.begin(), _atomicRadii.end(), atomicRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setScaleFactors(const vector<RealOpenMM>& scaleFactors) {
    _scaleFactors.resize(scaleFactors.size());
    copy(scaleFactors.begin(), scaleFactors.end(), _scaleFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getScaleFactors(vector<RealOpenMM>& scaleFactors) const {
    scaleFactors.resize(_scaleFactors.size());
    copy(_scaleFactors.begin(), _scaleFactors.end(), scaleFactors.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::setCharges(const vector<RealOpenMM>& charges) {
    _charges.resize(charges.size());
    copy(charges.begin(), charges.end(), _charges.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::getGrycukBornRadii(vector<RealOpenMM>& bornRadii) const {
    bornRadii.resize(_bornRadii.size());
    copy(_bornRadii.begin(), _bornRadii.end(), bornRadii.begin());
}

void AmoebaReferenceGeneralizedKirkwoodForce::calculateGrycukBornRadii(const vector<RealVec>& particlePositions) {

    const RealOpenMM zero      = 0.0;
    const RealOpenMM one       = 1.0;
    const RealOpenMM three     = 3.0;
    const RealOpenMM six       = 6.0;
    const RealOpenMM eight     = 8.0;
    const RealOpenMM sixteen   = 16.0;
    const RealOpenMM oneThird  = 1.0/3.0;
    const RealOpenMM bigRadius = 1000.0;

    _bornRadii.resize(_numParticles);
    for (unsigned int ii = 0; ii < _numParticles; ii++) {

        if (_atomicRadii[ii] <= zero) {
            _bornRadii[ii] = bigRadius;
            continue;
        }

        RealOpenMM bornSum = zero;
        for (unsigned int jj = 0; jj < _numParticles; jj++) {

            if (ii == jj || _atomicRadii[jj] < zero)continue;
          
            RealOpenMM xr       = particlePositions[jj][0] - particlePositions[ii][0];
            RealOpenMM yr       = particlePositions[jj][1] - particlePositions[ii][1];
            RealOpenMM zr       = particlePositions[jj][2] - particlePositions[ii][2];

            RealOpenMM r2       = xr*xr + yr*yr + zr*zr;
            RealOpenMM r        = SQRT(r2);

            RealOpenMM sk       = _atomicRadii[jj]*_scaleFactors[jj];
            RealOpenMM sk2      = sk*sk;

            if ((_atomicRadii[ii] + r) < sk) {
                RealOpenMM lik       = _atomicRadii[ii];
                RealOpenMM uik       = sk - r;  
                RealOpenMM lik3      = lik*lik*lik;
                RealOpenMM uik3      = uik*uik*uik;
                bornSum             -= (one/uik3 - one/lik3);
            }   
        
            RealOpenMM uik = r + sk; 
            RealOpenMM lik;
            if ((_atomicRadii[ii] + r) < sk) {
                lik = sk - r;  
            } else if (r < (_atomicRadii[ii] + sk)) {
                lik = _atomicRadii[ii];
            } else {
                lik = r - sk; 
            }   
        
            RealOpenMM l2          = lik*lik; 
            RealOpenMM l4          = l2*l2;
            RealOpenMM lr          = lik*r;
            RealOpenMM l4r         = l4*r;
        
            RealOpenMM u2          = uik*uik;
            RealOpenMM u4          = u2*u2;
            RealOpenMM ur          = uik*r;
            RealOpenMM u4r         = u4*r;
        
            RealOpenMM term        = (three*(r2-sk2) + six*u2 - eight*ur)/u4r - (three*(r2-sk2) + six*l2 - eight*lr)/l4r;
            bornSum               += term/sixteen;
        
        }
        bornSum        = one/(_atomicRadii[ii]*_atomicRadii[ii]*_atomicRadii[ii]) - bornSum;
        _bornRadii[ii] = (bornSum <= zero) ? bigRadius : POW(bornSum, -oneThird);
    }

    return;
}
