
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
#include <sstream>
#include <string.h>

#include "openmm/OpenMMException.h"
#include "GBVIParameters.h"
#include "SimTKOpenMMCommon.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

    GBVIParameters constructor 

    @param numberOfAtoms       number of atoms

    --------------------------------------------------------------------------------------- */

GBVIParameters::GBVIParameters( int numberOfAtoms ) : _numberOfAtoms(numberOfAtoms), 
                                                      _soluteDielectric(1.0),
                                                      _solventDielectric(78.3),
                                                      _electricConstant(-0.5*ONE_4PI_EPS0),
                                                      _cutoff(false),
                                                      _periodic(false),
                                                      _bornRadiusScalingMethod(0),
                                                      _quinticLowerLimitFactor(0.8),
                                                      _quinticUpperBornRadiusLimit(5.0) {

    _atomicRadii.resize( numberOfAtoms );
    _scaledRadii.resize( numberOfAtoms );
    _gammaParameters.resize( numberOfAtoms ); 

}

/**---------------------------------------------------------------------------------------

    GBVIParameters destructor 

    --------------------------------------------------------------------------------------- */

GBVIParameters::~GBVIParameters( ){
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int GBVIParameters::getNumberOfAtoms( void ) const {
    return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get electric constant

   @return electric constant

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getElectricConstant( void ) const {
    return _electricConstant;
}

/**---------------------------------------------------------------------------------------

   Get solvent dielectric

   @return solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getSolventDielectric( void ) const {
    return _solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric

   @param solventDielectric solvent dielectric

   --------------------------------------------------------------------------------------- */

void GBVIParameters::setSolventDielectric( RealOpenMM solventDielectric ){
    _solventDielectric = solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Get solute dielectric

   @return soluteDielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getSoluteDielectric( void ) const {
    return _soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solute dielectric

   @param soluteDielectric solute dielectric

   --------------------------------------------------------------------------------------- */

void GBVIParameters::setSoluteDielectric( RealOpenMM soluteDielectric ){
    _soluteDielectric = soluteDielectric;
}

/**---------------------------------------------------------------------------------------

    Get AtomicRadii array

    @return array of atomic radii

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVIParameters::getAtomicRadii( void ) const {
    return _atomicRadii;
}

/**---------------------------------------------------------------------------------------

    Set AtomicRadii array

    @param atomicRadii vector of atomic radii

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

    if( atomicRadii.size() == _atomicRadii.size() ){
        for( unsigned int ii = 0; ii < atomicRadii.size(); ii++ ){
            _atomicRadii[ii] = atomicRadii[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "GBVIParameters: input size for atomic radii does not agree w/ current size: input=";
        msg << atomicRadii.size();
        msg << " current size=" << _atomicRadii.size();
        throw OpenMM::OpenMMException(msg.str());
    }   

}

/**---------------------------------------------------------------------------------------

    Return scaled radii
    @return array 

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVIParameters::getScaledRadii( void ) const {
    return _scaledRadii;
}

/**---------------------------------------------------------------------------------------

    Set scaled radii

    @param scaledRadii  scaledRadii

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setScaledRadii( const RealOpenMMVector& scaledRadii ){

    if( scaledRadii.size() == _scaledRadii.size() ){
        for( unsigned int ii = 0; ii < scaledRadii.size(); ii++ ){
            _scaledRadii[ii] = scaledRadii[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVIParameters: input size for scaled radii does not agree w/ current size: input=";
        msg << scaledRadii.size();
        msg << " current size=" << _scaledRadii.size();
        throw OpenMM::OpenMMException(msg.str());
    }

}

/**---------------------------------------------------------------------------------------

    Return gamma parameters
    If not previously set, allocate space

    @return array 

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& GBVIParameters::getGammaParameters( void ) const {
    return _gammaParameters;
}

/**---------------------------------------------------------------------------------------

    Set gamma parameters

    @param gammas  gammas

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setGammaParameters( const RealOpenMMVector& gammas ){

    if( gammas.size() == _gammaParameters.size() ){
        for( unsigned int ii = 0; ii < gammas.size(); ii++ ){
            _gammaParameters[ii] = gammas[ii];
        }
    } else {
        std::stringstream msg;
        msg << "GBVIParameters: input size for gammas does not agree w/ current size: input=";
        msg << gammas.size();
        msg << " current size=" << _gammaParameters.size();
        throw OpenMM::OpenMMException(msg.str());
    }

}

/**---------------------------------------------------------------------------------------

      Set the force to use a cutoff.

      @param distance            the cutoff distance

      --------------------------------------------------------------------------------------- */

void GBVIParameters::setUseCutoff( RealOpenMM distance ) {

     _cutoff          = true;
     _cutoffDistance = distance;
}

/**---------------------------------------------------------------------------------------

      Get whether to use a cutoff.

      --------------------------------------------------------------------------------------- */

bool GBVIParameters::getUseCutoff() {
     return _cutoff;
}

/**---------------------------------------------------------------------------------------

      Get the cutoff distance.

      --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getCutoffDistance() {
     return _cutoffDistance;
}

/**---------------------------------------------------------------------------------------

      Set the force to use periodic boundary conditions.  This requires that a cutoff has
      also been set, and the smallest side of the periodic box is at least twice the cutoff
      distance.

      @param boxSize             the X, Y, and Z widths of the periodic box

      --------------------------------------------------------------------------------------- */

void GBVIParameters::setPeriodic( RealVec& boxSize ) {

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

bool GBVIParameters::getPeriodic() {
     return _periodic;
}

/**---------------------------------------------------------------------------------------

      Get the periodic box dimension

      --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVIParameters::getPeriodicBox() {
     return _periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

    Get tau prefactor

    @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric

    --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getTau( void ) const {

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

    Get bornRadiusScalingMethod

    @return bornRadiusScalingMethod

    --------------------------------------------------------------------------------------- */

int GBVIParameters::getBornRadiusScalingMethod( void ) const {
    return _bornRadiusScalingMethod;
}

/**---------------------------------------------------------------------------------------

    Set bornRadiusScalingMethod 

    @param bornRadiusScalingMethod bornRadiusScalingMethod

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setBornRadiusScalingMethod( int bornRadiusScalingMethod ){
    _bornRadiusScalingMethod    = bornRadiusScalingMethod;
}

/**---------------------------------------------------------------------------------------

    Get quinticLowerLimitFactor

    @return quinticLowerLimitFactor

    --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getQuinticLowerLimitFactor( void ) const {
    return _quinticLowerLimitFactor;
}

/**---------------------------------------------------------------------------------------

    Set quinticLowerLimitFactor 

    @param quinticLowerLimitFactor quinticLowerLimitFactor

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setQuinticLowerLimitFactor( RealOpenMM quinticLowerLimitFactor ){
    _quinticLowerLimitFactor    = quinticLowerLimitFactor;
}

/**---------------------------------------------------------------------------------------

    Get quinticUpperBornRadiusLimit

    @return quinticUpperBornRadiusLimit

    --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getQuinticUpperBornRadiusLimit( void ) const {
    return _quinticUpperBornRadiusLimit;
}

/**---------------------------------------------------------------------------------------

    Set quinticUpperBornRadiusLimit 

    @param quinticUpperBornRadiusLimit quinticUpperBornRadiusLimit

    --------------------------------------------------------------------------------------- */

void GBVIParameters::setQuinticUpperBornRadiusLimit( RealOpenMM quinticUpperBornRadiusLimit ){
    _quinticUpperBornRadiusLimit    = quinticUpperBornRadiusLimit;
}
