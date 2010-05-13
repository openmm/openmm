
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

#include "ObcSoftcoreParameters.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

// #define UseGromacsMalloc 1

#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

const std::string ObcSoftcoreParameters::ParameterFileName = std::string( "params.agb" );

/**---------------------------------------------------------------------------------------

   ObcSoftcoreParameters:

		Calculates for each atom

			(1) the van der Waal radii
         (2) volume
         (3) fixed terms in Obc equation gPol
         (4) list of atoms that should be excluded in calculating
				 force -- nonbonded atoms (1-2, and 1-3 atoms)

	Implementation:

		Slightly different sequence of calls when running on CPU vs GPU.
		Difference arise because the CPU-side data arrays for the Brook
		streams are allocated by the BrookStreamWrapper objects. These
		arrays are then used by ObcSoftcoreParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 ObcSoftcoreParameters* obcParameters = new ObcSoftcoreParameters( numberOfAtoms, log );
          obcParameters->initializeParameters( top );

		Gpu:

			obcParameters   = new ObcSoftcoreParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that ObcSoftcoreParameters destructor does not free arrays 
			
			obcParameters->setVdwRadii(  getBrookStreamWrapperAtIndex( GpuObc::obcVdwRadii  )->getData() );
			obcParameters->setVolume(    getBrookStreamWrapperAtIndex( GpuObc::obcVolume    )->getData() );
			obcParameters->setGPolFixed( getBrookStreamWrapperAtIndex( GpuObc::obcGpolFixed )->getData() );
			obcParameters->setBornRadii( getBrookStreamWrapperAtIndex( GpuObc::obcBornRadii )->getData() );
			
			obcParameters->setFreeArrays( false );
			
			obcParameters->initializeParameters( top );
 

   Issues:

		Tinker's atom radii are used. 
      The logic for mapping the Gromacs atom names to Tinker type may be incomplete;
      only tested for generic proteins
		see mapGmxAtomNameToTinkerAtomNumber()

   --------------------------------------------------------------------------------------- */


/**---------------------------------------------------------------------------------------

   ObcSoftcoreParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms
   @param obcType             OBC type (Eq. 7 or 8 in paper)

   --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::ObcSoftcoreParameters( int numberOfAtoms, ObcSoftcoreParameters::ObcType obcType ) : ImplicitSolventParameters( numberOfAtoms ), cutoff(false), periodic(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::ObcSoftcoreParameters";
   
   // ---------------------------------------------------------------------------------------

   _obcType                      = obcType;
   _dielectricOffset             = 0.009f;
   _ownScaledRadiusFactors       = 0;
   _scaledRadiusFactors          = NULL;

   _ownNonPolarScaleFactors      = 0;
   _nonPolarScaleFactors         = NULL;


   _nonPolarPreFactor            = (RealOpenMM) 2.25936;

   setObcTypeParameters( obcType );

}

/**---------------------------------------------------------------------------------------

   ObcSoftcoreParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::~ObcSoftcoreParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::~ObcSoftcoreParameters";
   
   // ---------------------------------------------------------------------------------------

   // in GPU runs, arrays may be 'owned' by BrookStreamWrapper -- hence they should not
   // be freed here, i.e., _freeArrays should be 'false'   

#ifdef UseGromacsMalloc

/*
   if( _freeArrays ){

      if( _vdwRadii != NULL ){
         save_free( "_vdwRadii", __FILE__, __LINE__, _vdwRadii );
      }
   
   } */

#else

   if( _ownScaledRadiusFactors ){
      delete[] _scaledRadiusFactors;
   }

   if( _ownNonPolarScaleFactors ){
      delete[] _nonPolarScaleFactors;
   }

/*
   if( getFreeArrays() ){

   } */

#endif

}

/**---------------------------------------------------------------------------------------

   Get OBC type

   @return OBC type

   --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters::ObcType ObcSoftcoreParameters::getObcType( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getObcType:";

   // ---------------------------------------------------------------------------------------

   return _obcType;
}

/**---------------------------------------------------------------------------------------

   Set OBC type specific parameters

   @param obcType OBC type (ObcTypeI or ObcTypeII -- Eq. 7 or 8)

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setObcTypeParameters( ObcSoftcoreParameters::ObcType obcType ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::setObcTypeParameters:";

   // ---------------------------------------------------------------------------------------

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

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get dielectricOffset

   @return _dielectricOffset

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getDielectricOffset( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getDielectricOffset:";

   // ---------------------------------------------------------------------------------------

   return _dielectricOffset;
}

/**---------------------------------------------------------------------------------------

   Get alpha OBC (Eqs. 6 & 7) in Proteins paper

   @return alphaObc

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getAlphaObc( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getAlphaObc:";

   // ---------------------------------------------------------------------------------------

   return _alphaObc;
}

/**---------------------------------------------------------------------------------------

   Get beta OBC (Eqs. 6 & 7) in Proteins paper

   @return betaObc

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getBetaObc( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getBetaObc:";

   // ---------------------------------------------------------------------------------------

   return _betaObc;
}

/**---------------------------------------------------------------------------------------

   Get gamma OBC (Eqs. 6 & 7) in Proteins paper

   @return gammaObc

   --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getGammaObc( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getGammaObc:";

   // ---------------------------------------------------------------------------------------

   return _gammaObc;
}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

RealOpenMM* ObcSoftcoreParameters::getAtomicRadii( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   RealOpenMM* atomicRadii = ImplicitSolventParameters::getAtomicRadii();

   // if dielectric offset applied, then unapply

   return atomicRadii;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii array of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setAtomicRadii( RealOpenMM* atomicRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii vector of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nObcSoftcoreParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Return OBC scale factors
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMM* ObcSoftcoreParameters::getScaledRadiusFactors( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadiusFactors == NULL ){
      ObcSoftcoreParameters* localThis = const_cast<ObcSoftcoreParameters* const>(this);
      localThis->_scaledRadiusFactors    = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownScaledRadiusFactors = true;
      memset( _scaledRadiusFactors, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }   
   return _scaledRadiusFactors;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether scale factors array should be deleted

   @param ownScaledRadiusFactors flag indicating whether scale factors 
                                 array should be deleted

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setOwnScaleFactors( int ownScaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownScaledRadiusFactors = ownScaledRadiusFactors;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set OBC scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setScaledRadiusFactors( RealOpenMM* scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadiusFactors && _scaledRadiusFactors != scaledRadiusFactors ){
      delete[] _scaledRadiusFactors;
      _ownScaledRadiusFactors = false;
   }

   _scaledRadiusFactors = scaledRadiusFactors;

   return 0;

}

#if RealOpenMMType == 0

/**---------------------------------------------------------------------------------------

   Set OBC scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setScaledRadiusFactors( float* scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadiusFactors == NULL ){
      _scaledRadiusFactors    = new RealOpenMM[getNumberOfAtoms()];
      _ownScaledRadiusFactors = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _scaledRadiusFactors[ii] = (RealOpenMM) scaledRadiusFactors[ii];
   }

   return 0;

}

#endif

/**---------------------------------------------------------------------------------------

   Set OBC scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadiusFactors && _scaledRadiusFactors != NULL ){
      delete[] _scaledRadiusFactors;
   }
   _ownScaledRadiusFactors = true;
   _scaledRadiusFactors    = new RealOpenMM[getNumberOfAtoms()];
   for( int ii = 0; ii < (int) scaledRadiusFactors.size(); ii++ ){
      _scaledRadiusFactors[ii] = scaledRadiusFactors[ii];
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Map Gmx atom name to Tinker atom number (Simbios)

   @param atomName            atom name (CA, HA, ...); upper and lower case should both work
   @param log                 if set, then print error messages to log file

   @return Tinker atom number if atom name is valid; else return -1

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   static int mapCreated = 0;
   static int atomNameMap[26];
        
   // ---------------------------------------------------------------------------------------

   // set up atomNameMap array on first call to this method

   // atomNameMap[ii] = Tinker atom number
   // where ii = (the ASCII index - 65) of the first character in the
   // input atom name; name may be lower case

   if( !mapCreated ){

      mapCreated = 1;

      for( int ii = 0; ii < 26; ii++ ){
         atomNameMap[ii] = -1;
      }

      // H
      atomNameMap[7]  = 1;

      // C
      atomNameMap[2]  = 6;

      // N
      atomNameMap[13] = 7;

      // O
      atomNameMap[14] = 8;

      // S
      atomNameMap[18] = 16;
   }

   // map first letter in atom name to Tinker atom number

   int firstAsciiValue = ((int) atomName[0]) - 65;

   // check for lower case

   if( firstAsciiValue > 25 ){
      firstAsciiValue -= 32;
   }

   // validate

   if( firstAsciiValue < 0 || firstAsciiValue > 25 ){ 
      if( log != NULL ){
         (void) fprintf( log, "Atom name=<%s> unrecognized.", atomName );
      }
      (void) fprintf( stderr, "Atom name=<%s> unrecognized.", atomName );
      return -1;
   }
   return atomNameMap[firstAsciiValue];
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state
   
   @param title               title (optional)
      
   @return string
      
   --------------------------------------------------------------------------------------- */

std::string ObcSoftcoreParameters::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nObcSoftcoreParameters::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << ImplicitSolventParameters::getStateString( title );

   std::string tab           = getStringTab();

   if( getObcType() == ObcTypeI ){
      message << tab << "OBC type:    Type I";
   } else {
      message << tab << "OBC type:    Type II";
   }
   message << tab << "Alpha:                " << getAlphaObc();
   message << tab << "Beta:                 " << getBetaObc();
   message << tab << "Gamma:                " << getGammaObc();

   return message.str();


}

/**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setUseCutoff( RealOpenMM distance ) {

    cutoff = true;
    cutoffDistance = distance;
    return 0;
}

/**---------------------------------------------------------------------------------------

     Get whether to use a cutoff.

     --------------------------------------------------------------------------------------- */

bool ObcSoftcoreParameters::getUseCutoff() {
    return cutoff;
}

/**---------------------------------------------------------------------------------------

     Get the cutoff distance.

     --------------------------------------------------------------------------------------- */

RealOpenMM ObcSoftcoreParameters::getCutoffDistance() {
    return cutoffDistance;
}

/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setPeriodic( RealOpenMM* boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
    return 0;
}

/**---------------------------------------------------------------------------------------

     Get whether to use periodic boundary conditions.

     --------------------------------------------------------------------------------------- */

bool ObcSoftcoreParameters::getPeriodic() {
    return periodic;
}

/**---------------------------------------------------------------------------------------

     Get the periodic box dimension

     --------------------------------------------------------------------------------------- */

const RealOpenMM* ObcSoftcoreParameters::getPeriodicBox() {
    return periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

   Return non-polar scale factors
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMM* ObcSoftcoreParameters::getNonPolarScaleFactors( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "ObcSoftcoreParameters::getNonPolarScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _nonPolarScaleFactors == NULL ){
      ObcSoftcoreParameters* localThis     = const_cast<ObcSoftcoreParameters* const>(this);
      localThis->_nonPolarScaleFactors     = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownNonPolarScaleFactors  = true;
      memset( _nonPolarScaleFactors, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }   
   return _nonPolarScaleFactors;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether scale factors array should be deleted

   @param ownNonPolarScaleFactors flag indicating whether scale factors 
                                 array should be deleted

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setOwnNonPolarScaleFactors( int ownNonPolarScaleFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownNonPolarScaleFactors = ownNonPolarScaleFactors;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set non-polar scale factors

   @param nonPolarScaleFactors  nonPolarScaleFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setNonPolarScaleFactors( const RealOpenMMVector& nonPolarScaleFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setNonPolarScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownNonPolarScaleFactors ){
      delete[] _nonPolarScaleFactors;
   }
   _ownNonPolarScaleFactors = true;
   _nonPolarScaleFactors    = new RealOpenMM[nonPolarScaleFactors.size()];

   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _nonPolarScaleFactors[ii] = nonPolarScaleFactors[ii];
   }

   return 0;

}

#if RealOpenMMType == 0

/**---------------------------------------------------------------------------------------

   Set non-polar scale factors

   @param nonPolarScaleFactors  nonPolarScaleFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setNonPolarScaleFactors( float* nonPolarScaleFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setNonPolarScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _nonPolarScaleFactors == NULL ){
      _nonPolarScaleFactors    = new RealOpenMM[getNumberOfAtoms()];
      _ownNonPolarScaleFactors = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _nonPolarScaleFactors[ii] = (RealOpenMM) nonPolarScaleFactors[ii];
   }

   return 0;

}

#endif

/**---------------------------------------------------------------------------------------

   Set OBC scale factors

   @param nonPolarScaleFactors  nonPolarScaleFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int ObcSoftcoreParameters::setNonPolarPrefactor( RealOpenMM nonPolarPreFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setNonPolarScaleFactors";

   // ---------------------------------------------------------------------------------------

   _nonPolarPreFactor = nonPolarPreFactor;
   return 0;
}

