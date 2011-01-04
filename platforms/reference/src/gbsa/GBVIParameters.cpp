
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

#include "GBVIParameters.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

using std::vector;
using OpenMM::RealVec;

const std::string GBVIParameters::ParameterFileName = std::string( "params.agb" );

/**---------------------------------------------------------------------------------------

   GBVIParameters:

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
		arrays are then used by GBVIParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 GBVIParameters* gb_VIParameters = new GBVIParameters( numberOfAtoms, log );
          gb_VIParameters->initializeParameters( top );

		Gpu:

			gb_VIParameters   = new GBVIParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that GBVIParameters destructor does not free arrays 
			
			gb_VIParameters->setVdwRadii(  getBrookStreamWrapperAtIndex( GpuObc::gb_VIVdwRadii  )->getData() );
			gb_VIParameters->setVolume(    getBrookStreamWrapperAtIndex( GpuObc::gb_VIVolume    )->getData() );
			gb_VIParameters->setGPolFixed( getBrookStreamWrapperAtIndex( GpuObc::gb_VIGpolFixed )->getData() );
			gb_VIParameters->setBornRadii( getBrookStreamWrapperAtIndex( GpuObc::gb_VIBornRadii )->getData() );
			
			gb_VIParameters->setFreeArrays( false );
			
			gb_VIParameters->initializeParameters( top );
 

   Issues:

		Tinker's atom radii are used. 
      The logic for mapping the Gromacs atom names to Tinker type may be incomplete;
      only tested for generic proteins
		see mapGmxAtomNameToTinkerAtomNumber()

   --------------------------------------------------------------------------------------- */


/**---------------------------------------------------------------------------------------

   GBVIParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms

   --------------------------------------------------------------------------------------- */

GBVIParameters::GBVIParameters( int numberOfAtoms ) : ImplicitSolventParameters( numberOfAtoms ), cutoff(false), periodic(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::GBVIParameters";
   
   // ---------------------------------------------------------------------------------------

   _ownScaledRadii       = 0;
   _scaledRadii          = NULL;
   _ownGammaParameters   = 0;
   _gammaParameters      = NULL;

}

/**---------------------------------------------------------------------------------------

   GBVIParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

GBVIParameters::~GBVIParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::~GBVIParameters";
   
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

   if( _ownScaledRadii ){
      delete[] _scaledRadii;
   }
   delete[] _gammaParameters;
/*
   if( getFreeArrays() ){

   } */

#endif

}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

RealOpenMM* GBVIParameters::getAtomicRadii( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   RealOpenMM* atomicRadii = ImplicitSolventParameters::getAtomicRadii();

   return atomicRadii;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii array of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setAtomicRadii( RealOpenMM* atomicRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii vector of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGBVIParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Return scaled radii
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVIParameters::getScaledRadii( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::getScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadii == NULL ){
      GBVIParameters* localThis  = const_cast<GBVIParameters* const>(this);
      localThis->_scaledRadii    = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownScaledRadii = true;
      memset( _scaledRadii, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }   
   return _scaledRadii;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether scale factors array should be deleted

   @param ownScaledRadii flag indicating whether scale factors 
                                 array should be deleted

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setOwnScaledRadii( int ownScaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownScaledRadii = ownScaledRadii;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setScaledRadii( RealOpenMM* scaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadii && _scaledRadii != scaledRadii ){
      delete[] _scaledRadii;
      _ownScaledRadii = false;
   }

   _scaledRadii = scaledRadii;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setScaledRadii( const RealOpenMMVector& scaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadii && _scaledRadii != NULL ){
      delete[] _scaledRadii;
   }
   _ownScaledRadii = true;
   _scaledRadii    = new RealOpenMM[getNumberOfAtoms()];
   for( int ii = 0; ii < (int) scaledRadii.size(); ii++ ){
      _scaledRadii[ii] = scaledRadii[ii];
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Return gamma parameters
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

RealOpenMM* GBVIParameters::getGammaParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::getGammaParameters";

   // ---------------------------------------------------------------------------------------

   if( _gammaParameters == NULL ){
      GBVIParameters* localThis = const_cast<GBVIParameters* const>(this);
      localThis->_gammaParameters    = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownGammaParameters = true;
      memset( _gammaParameters, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }   
   return _gammaParameters;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether scale factors array should be deleted

   @param ownGammaParameters   flag indicating whether gamma parameter
                               array should be deleted

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setOwnGammaParameters( int ownGammaParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownGammaParameters = ownGammaParameters;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gamma parameters

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setGammaParameters( RealOpenMM* gammas ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setGammas";

   // ---------------------------------------------------------------------------------------

   if( _ownGammaParameters && _gammaParameters != gammas ){
      delete[] _gammaParameters;
      _ownGammaParameters = false;
   }

   _gammaParameters = gammas;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gammas

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVIParameters::setGammaParameters( const RealOpenMMVector& gammas ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setGammas";

   // ---------------------------------------------------------------------------------------

   if( _ownGammaParameters && _gammaParameters != NULL ){
      delete[] _gammaParameters;
   }
   _ownGammaParameters = true;

   _gammaParameters    = new RealOpenMM[getNumberOfAtoms()];
   for( int ii = 0; ii < (int) gammas.size(); ii++ ){
      _gammaParameters[ii] = gammas[ii];
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state
   
   @param title               title (optional)
      
   @return string
      
   --------------------------------------------------------------------------------------- */

std::string GBVIParameters::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << ImplicitSolventParameters::getStateString( title );

   std::string tab           = getStringTab();

   return message.str();

}

/**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int GBVIParameters::setUseCutoff( RealOpenMM distance ) {

    cutoff = true;
    cutoffDistance = distance;
    return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

     Get whether to use a cutoff.

     --------------------------------------------------------------------------------------- */

bool GBVIParameters::getUseCutoff() {
    return cutoff;
}

/**---------------------------------------------------------------------------------------

     Get the cutoff distance.

     --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getCutoffDistance() {
    return cutoffDistance;
}

/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int GBVIParameters::setPeriodic( RealVec& boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
    return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

     Get whether to use periodic boundary conditions.

     --------------------------------------------------------------------------------------- */

bool GBVIParameters::getPeriodic() {
    return periodic;
}

/**---------------------------------------------------------------------------------------

     Get the periodic box dimension

     --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVIParameters::getPeriodicBox() {
    return periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

   Get tau prefactor

   @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVIParameters::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVIParameters::getTau:";

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
