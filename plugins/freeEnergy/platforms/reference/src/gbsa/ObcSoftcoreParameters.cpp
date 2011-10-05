
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
#include <string.h>
#include <sstream>

#include "openmm/OpenMMException.h"
#include "ObcSoftcoreParameters.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"

/**---------------------------------------------------------------------------------------

    ObcSoftcoreParameters constructor  

    @param numberOfAtoms       number of atoms
    @param obcType             OBC type (Eq. 7 or 8 in paper)

    --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::ObcSoftcoreParameters( int numberOfAtoms, ObcSoftcoreParameters::ObcType obcType ) :
                                              _numberOfAtoms(numberOfAtoms),                                          
                                              _obcType(obcType),
                                              _dielectricOffset(0.009),
                                              _nonPolarPreFactor(2.25936),
                                              _soluteDielectric(1.0), 
                                              _solventDielectric(78.3),
                                               _probeRadius(0.14),
                                              _electricConstant(-0.5*ONE_4PI_EPS0),
                                              _pi4Asolv( 28.3919551),
                                              _cutoff(false),
                                              _periodic(false) {

    _atomicRadii.resize( numberOfAtoms );
    _scaledRadiusFactors.resize( numberOfAtoms );
    _nonPolarScaleFactors.resize( numberOfAtoms );

    setObcTypeParameters( obcType );

}

/**---------------------------------------------------------------------------------------

    ObcSoftcoreParameters destructor  

    --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::~ObcSoftcoreParameters( ){
}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::getNumberOfAtoms( void ) const {
    return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

    Get OBC type

    @return OBC type

    --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::ObcType ObcSoftcoreParameters::getObcType( void ) const {
    return _obcType;
}

/**---------------------------------------------------------------------------------------

    Set OBC type specific parameters

    @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)

    --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setObcTypeParameters( ObcSoftcoreParameters::ObcType obcType ){

    if( obcType == ObcTypeI ){
       _alphaObc   = 0.8f;
       _betaObc    = 0.0f;
       _gammaObc   = 2.91f;
    } else {
       _alphaObc   = 1.0f;
       _betaObc    = 0.8f;
       _gammaObc   = 4.85f;
    }
    _obcType = obcType;
}

/**---------------------------------------------------------------------------------------

    Get dielectricOffset

    @return _dielectricOffset

    --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getDielectricOffset( void ) const {
    return _dielectricOffset;
}

/**---------------------------------------------------------------------------------------

    Get alpha OBC (Eqs. 6 & 7) in Proteins paper

    @return alphaObc

    --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getAlphaObc( void ) const {
    return _alphaObc;
}

/**---------------------------------------------------------------------------------------

    Get beta OBC (Eqs. 6 & 7) in Proteins paper

    @return betaObc

    --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getBetaObc( void ) const {
    return _betaObc;
}

/**---------------------------------------------------------------------------------------

    Get gamma OBC (Eqs. 6 & 7) in Proteins paper

    @return gammaObc

    --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getGammaObc( void ) const {
    return _gammaObc;
}

/**---------------------------------------------------------------------------------------

   Get solvent dielectric

   @return solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getSolventDielectric( void ) const {
    return _solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric

   @param solventDielectric solvent dielectric

   --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setSolventDielectric( RealOpenMM solventDielectric ){
    _solventDielectric = solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Get solute dielectric

   @return soluteDielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getSoluteDielectric( void ) const {
    return _soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solute dielectric

   @param soluteDielectric solute dielectric

   --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setSoluteDielectric( RealOpenMM soluteDielectric ){
    _soluteDielectric = soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Get electric constant 

   @return electricConstant

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getElectricConstant( void ) const {
    return _electricConstant;
}

/**---------------------------------------------------------------------------------------

   Get probe radius 

   @return probeRadius

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getProbeRadius( void ) const {
    return _probeRadius;
}

/**---------------------------------------------------------------------------------------

   Set probe radius  

   @param probeRadius   probe radius

   --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setProbeRadius( RealOpenMM probeRadius ){
    _probeRadius = probeRadius;
}

/**---------------------------------------------------------------------------------------

   Get pi*4*Asolv:  used in ACE approximation for nonpolar term  
         ((RealOpenMM) M_PI)*4.0f*0.0049*1000.0; (Still) 
         ((RealOpenMM) M_PI)*4.0f*0.0054*1000.0; (OBC) 

   @return pi4Asolv

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getPi4Asolv( void ) const {
    return _pi4Asolv;
}

/**---------------------------------------------------------------------------------------

    Get AtomicRadii array

    @return array of atomic radii

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& ObcSoftcoreParameters::getAtomicRadii( void ) const {
    return _atomicRadii;
}

/**---------------------------------------------------------------------------------------

    Set AtomicRadii array

    @param atomicRadii vector of atomic radii

    --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

    if( atomicRadii.size() == _atomicRadii.size() ){
        for( unsigned int ii = 0; ii < atomicRadii.size(); ii++ ){
            _atomicRadii[ii] = atomicRadii[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcSoftcoreParameters: input size for atomic radii does not agree w/ current size: input=";
        msg << atomicRadii.size();
        msg << " current size=" << _atomicRadii.size();
        throw OpenMM::OpenMMException(msg.str());
    }   
}

/**---------------------------------------------------------------------------------------

    Return OBC scale factors

    @return array 

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& ObcSoftcoreParameters::getScaledRadiusFactors( void ) const {
    return _scaledRadiusFactors;
}

/**---------------------------------------------------------------------------------------

    Set OBC scale factors

    @param scaledRadiusFactors  scaledRadiusFactors

    --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors ){

    if( scaledRadiusFactors.size() == _scaledRadiusFactors.size() ){
        for( unsigned int ii = 0; ii < scaledRadiusFactors.size(); ii++ ){
            _scaledRadiusFactors[ii] = scaledRadiusFactors[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcSoftcoreParameters: input size for scaled radius factors does not agree w/ current size: input=";
        msg << scaledRadiusFactors.size();
        msg << " current size=" << _scaledRadiusFactors.size();
        throw OpenMM::OpenMMException(msg.str());
    }   
}

/**---------------------------------------------------------------------------------------

      Set the force to use a cutoff.

      @param distance            the cutoff distance

      --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setUseCutoff( RealOpenMM distance ) {
     _cutoff         = true;
     _cutoffDistance = distance;
}

/**---------------------------------------------------------------------------------------

      Get whether to use a cutoff.

      --------------------------------------------------------------------------------------- */

bool ObcSoftcoreParameters::getUseCutoff() {
     return _cutoff;
}

/**---------------------------------------------------------------------------------------

      Get the cutoff distance.

      --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getCutoffDistance() {
     return _cutoffDistance;
}

/**---------------------------------------------------------------------------------------

      Set the force to use periodic boundary conditions.  This requires that a cutoff has
      also been set, and the smallest side of the periodic box is at least twice the cutoff
      distance.

      @param boxSize             the X, Y, and Z widths of the periodic box

      --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setPeriodic( RealOpenMM* boxSize ) {

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

bool ObcSoftcoreParameters::getPeriodic() {
     return _periodic;
}

/**---------------------------------------------------------------------------------------

      Get the periodic box dimension

      --------------------------------------------------------------------------------------- */

const RealOpenMM* ObcSoftcoreParameters::getPeriodicBox() {
     return _periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

    Return non-polar scale factors

    @return array 

    --------------------------------------------------------------------------------------- */

const RealOpenMMVector& ObcSoftcoreParameters::getNonPolarScaleFactors( void ) const {
    return _nonPolarScaleFactors;
}

/**---------------------------------------------------------------------------------------

    Set non-polar scale factors

    @param nonPolarScaleFactors  nonPolarScaleFactors

    --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setNonPolarScaleFactors( const RealOpenMMVector& nonPolarScaleFactors ){

    if( nonPolarScaleFactors.size() == _nonPolarScaleFactors.size() ){
        for( unsigned int ii = 0; ii < nonPolarScaleFactors.size(); ii++ ){
            _nonPolarScaleFactors[ii] = nonPolarScaleFactors[ii];
        }   
    } else {
        std::stringstream msg;
        msg << "ObcSoftcoreParameters: input size for non-polar scale factors does not agree w/ current size: input=";
        msg << nonPolarScaleFactors.size();
        msg << " current size=" << _nonPolarScaleFactors.size();
        throw OpenMM::OpenMMException(msg.str());
    }   
}

/**---------------------------------------------------------------------------------------

    Set OBC scale factors

    @param nonPolarScaleFactors  nonPolarScaleFactors

    --------------------------------------------------------------------------------------- */

void ObcSoftcoreParameters::setNonPolarPrefactor( RealOpenMM nonPolarPreFactor ){
    _nonPolarPreFactor = nonPolarPreFactor;
}

