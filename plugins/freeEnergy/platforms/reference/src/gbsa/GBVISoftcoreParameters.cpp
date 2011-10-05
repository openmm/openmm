/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
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

#include <math.h>
#include <iostream>
#include <sstream>
#include <string.h>

#include "openmm/OpenMMException.h"
#include "GBVISoftcoreParameters.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

/**---------------------------------------------------------------------------------------

   GBVISoftcoreParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters::GBVISoftcoreParameters( int numberOfAtoms ) : _numberOfAtoms(numberOfAtoms), _soluteDielectric(1.0), _solventDielectric(78.3),
                                                                      _electricConstant(-0.5*ONE_4PI_EPS0), _quinticLowerLimitFactor(0.8), _bornRadiusScalingSoftcoreMethod(NoScaling),
                                                                      _cutoff(false), _periodic(false) {

    // ---------------------------------------------------------------------------------------

    _atomicRadii.resize( numberOfAtoms );
    _scaledRadii.resize( numberOfAtoms );
    _gammaParameters.resize( numberOfAtoms );
    _bornRadiusScaleFactors.resize( numberOfAtoms );
 
    setQuinticUpperBornRadiusLimit( static_cast<RealOpenMM>(5.0) );

}

/**---------------------------------------------------------------------------------------

   GBVISoftcoreParameters destructor

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters::~GBVISoftcoreParameters( ){
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::getNumberOfAtoms( void ) const {
    return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get electric constant

   @return electric constant

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getElectricConstant( void ) const {
    return _electricConstant;
}

/**---------------------------------------------------------------------------------------

   Get solvent dielectric

   @return solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getSolventDielectric( void ) const {
    return _solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric

   @param solventDielectric solvent dielectric

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setSolventDielectric( RealOpenMM solventDielectric ){
    _solventDielectric = solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Get solute dielectric

   @return soluteDielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getSoluteDielectric( void ) const {
    return _soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solute dielectric

   @param soluteDielectric solute dielectric

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setSoluteDielectric( RealOpenMM soluteDielectric ){
    _soluteDielectric = soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Get the quintic spline lower limit factor

   @return quintic spline lower limit factor

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticLowerLimitFactor( void ) const {
    return _quinticLowerLimitFactor;
}

/**---------------------------------------------------------------------------------------

   Set the quintic spline lower limit factor

   @param quintic spline lower limit factor

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setQuinticLowerLimitFactor( RealOpenMM quinticLowerLimitFactor ){
    _quinticLowerLimitFactor = quinticLowerLimitFactor;
}

/**---------------------------------------------------------------------------------------

   Get the quintic spline upper limit 

   @return quintic spline upper limit

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticUpperBornRadiusLimit( void ) const {
    return _quinticUpperBornRadiusLimit;
}

/**---------------------------------------------------------------------------------------

   Set the quintic spline upper limit

   @param quintic spline upper limit

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setQuinticUpperBornRadiusLimit( RealOpenMM quinticUpperBornRadiusLimit ){
    _quinticUpperBornRadiusLimit  = quinticUpperBornRadiusLimit;
    _quinticUpperSplineLimit      = POW( _quinticUpperBornRadiusLimit, static_cast<RealOpenMM>(-3.0) ); 
}

/**---------------------------------------------------------------------------------------

   Get the quintic upper spline limit

   @return the quintic upper spline limit

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticUpperSplineLimit( void ) const {
    return _quinticUpperSplineLimit;
}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVISoftcoreParameters::getAtomicRadii( void ) const {
    return _atomicRadii;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii vector

   @param atomicRadii vector of atomic radii

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

   // ---------------------------------------------------------------------------------------

    if( atomicRadii.size() == _atomicRadii.size() ){
        for( unsigned int ii = 0; ii < atomicRadii.size(); ii++ ){
            _atomicRadii[ii] = atomicRadii[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVISoftcoreParameters: input size for atomic radii does not agree w/ current size: input=";
        msg << atomicRadii.size();
        msg << " current size=" << _atomicRadii.size();
        throw OpenMM::OpenMMException(msg.str());
    }
}

/**---------------------------------------------------------------------------------------

   Return scaled radii

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVISoftcoreParameters::getScaledRadii( void ) const {
    return _scaledRadii;
}

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setScaledRadii( const RealOpenMMVector& scaledRadii ){

   // ---------------------------------------------------------------------------------------

    if( scaledRadii.size() == _scaledRadii.size() ){
        for( unsigned int ii = 0; ii < scaledRadii.size(); ii++ ){
            _scaledRadii[ii] = scaledRadii[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVISoftcoreParameters: input size for scaled radii does not agree w/ current size: input=";
        msg << scaledRadii.size();
        msg << " current size=" << _scaledRadii.size();
        throw OpenMM::OpenMMException(msg.str());
    }
}

/**---------------------------------------------------------------------------------------

   Return gamma parameters

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVISoftcoreParameters::getGammaParameters( void ) const {
    return _gammaParameters;
}

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gammas

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setGammaParameters( const RealOpenMMVector& gammas ){

   // ---------------------------------------------------------------------------------------

    if( gammas.size() == _gammaParameters.size() ){
        for( unsigned int ii = 0; ii < gammas.size(); ii++ ){
            _gammaParameters[ii] = gammas[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVISoftcoreParameters: input size for gammas does not agree w/ current size: input=";
        msg << gammas.size();
        msg << " current size=" << _gammaParameters.size();
        throw OpenMM::OpenMMException(msg.str());
    }
}

/**---------------------------------------------------------------------------------------

   Return BornRadiusScaleFactors

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVISoftcoreParameters::getBornRadiusScaleFactors( void ) const {
    return _bornRadiusScaleFactors;
}

/**---------------------------------------------------------------------------------------

   Set bornRadiusScaleFactors parameters

   @param bornRadiusScaleFactors  bornRadiusScaleFactors

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setBornRadiusScaleFactors( const RealOpenMMVector& bornRadiusScaleFactors ){

   // ---------------------------------------------------------------------------------------

    if( bornRadiusScaleFactors.size() == _bornRadiusScaleFactors.size() ){
        for( int ii = 0; ii < (int) bornRadiusScaleFactors.size(); ii++ ){
            _bornRadiusScaleFactors[ii] = bornRadiusScaleFactors[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVISoftcoreParameters: input size for bornRadiusScaleFactors does not agree w/ current size: input=";
        msg << bornRadiusScaleFactors.size();
        msg << " current size=" << _bornRadiusScaleFactors.size();
        throw OpenMM::OpenMMException(msg.str());
    }
}

/**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance

     --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setUseCutoff( RealOpenMM distance ) {
     _cutoff         = true;
     _cutoffDistance = distance;
}

/**---------------------------------------------------------------------------------------

     Get whether to use a cutoff.

     --------------------------------------------------------------------------------------- */

bool GBVISoftcoreParameters::getUseCutoff() {
     return _cutoff;
}

/**---------------------------------------------------------------------------------------

     Get the cutoff distance.

     --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getCutoffDistance() {
     return _cutoffDistance;
}

/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setPeriodic( RealOpenMM* boxSize ) {

    assert(_cutoff);

    assert(boxSize[0] >= 2.0*_cutoffDistance);
    assert(boxSize[1] >= 2.0*_cutoffDistance);
    assert(boxSize[2] >= 2.0*_cutoffDistance);

    _periodic           = true;
    _periodicBoxSize[0] = boxSize[0];
    _periodicBoxSize[1] = boxSize[1];
    _periodicBoxSize[2] = boxSize[2];
}

/**---------------------------------------------------------------------------------------

     Get whether to use periodic boundary conditions.

     --------------------------------------------------------------------------------------- */

bool GBVISoftcoreParameters::getPeriodic() {
    return _periodic;
}

/**---------------------------------------------------------------------------------------

     Get the periodic box dimension

     --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVISoftcoreParameters::getPeriodicBox() {
    return _periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

   Get tau prefactor

   @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getTau( void ) const {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero = 0.0;
    static const RealOpenMM one  = 1.0;
 
    // ---------------------------------------------------------------------------------------
 
    RealOpenMM tau;
    if( getSoluteDielectric() != zero && getSolventDielectric() != zero ){
       tau = (one/getSoluteDielectric()) - (one/getSolventDielectric());
    } else {
       tau = zero;
    }   
 
    return tau;
}

/**---------------------------------------------------------------------------------------

   Get Born radii switching function method

   @return method

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters::BornRadiusScalingSoftcoreMethod GBVISoftcoreParameters::getBornRadiusScalingSoftcoreMethod( void ) const {
    return _bornRadiusScalingSoftcoreMethod;
}

/**---------------------------------------------------------------------------------------

   Set Born radii switching function method

   @param method

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setBornRadiusScalingSoftcoreMethod( BornRadiusScalingSoftcoreMethod bornRadiusScalingSoftcoreMethod ){
    _bornRadiusScalingSoftcoreMethod = bornRadiusScalingSoftcoreMethod;
}
