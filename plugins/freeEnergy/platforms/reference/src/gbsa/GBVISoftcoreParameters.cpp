
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

#include "GBVISoftcoreParameters.h"
#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

// #define UseGromacsMalloc 1

#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

/**---------------------------------------------------------------------------------------

   GBVISoftcoreParameters:

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
		arrays are then used by GBVISoftcoreParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 GBVISoftcoreParameters* gb_VIParameters = new GBVISoftcoreParameters( numberOfAtoms, log );
          gb_VIParameters->initializeParameters( top );

		Gpu:

			gb_VIParameters   = new GBVISoftcoreParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that GBVISoftcoreParameters destructor does not free arrays 
			
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

   GBVISoftcoreParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters::GBVISoftcoreParameters( int numberOfAtoms ) : ImplicitSolventParameters( numberOfAtoms ), cutoff(false), periodic(false) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::GBVISoftcoreParameters";
   
   // ---------------------------------------------------------------------------------------

   _ownScaledRadii                      = 0;
   _scaledRadii                         = NULL;
   _ownGammaParameters                  = 0;
   _gammaParameters                     = NULL;
   _ownBornRadiusScaleFactors           = 0;
   _bornRadiusScaleFactors              = NULL;

   _bornRadiusScalingSoftcoreMethod     = NoScaling;
   _quinticLowerLimitFactor             = static_cast<RealOpenMM>(0.8);
   setQuinticUpperBornRadiusLimit( static_cast<RealOpenMM>(5.0) );

}

/**---------------------------------------------------------------------------------------

   GBVISoftcoreParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters::~GBVISoftcoreParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::~GBVISoftcoreParameters";
   
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
   delete[] _bornRadiusScaleFactors;
/*
   if( getFreeArrays() ){

   } */

#endif

}

/**---------------------------------------------------------------------------------------

   Get the quintic spline lower limit factor

   @return quintic spline lower limit factor

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticLowerLimitFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GBVISoftcoreParameters::getQuinticLowerLimitFactor:";

   // ---------------------------------------------------------------------------------------

   return _quinticLowerLimitFactor;
}

/**---------------------------------------------------------------------------------------

   Set the quintic spline lower limit factor

   @param quintic spline lower limit factor
   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setQuinticLowerLimitFactor( RealOpenMM quinticLowerLimitFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GBVISoftcoreParameters::setQuinticLowerLimitFactor:";

   // ---------------------------------------------------------------------------------------

   _quinticLowerLimitFactor = quinticLowerLimitFactor;
   return 0;
}

/**---------------------------------------------------------------------------------------

   Get the quintic spline upper limit 

   @return quintic spline upper limit

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticUpperBornRadiusLimit( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GBVISoftcoreParameters::getQuinticUpperBornRadiusLimit:";

   // ---------------------------------------------------------------------------------------

   return _quinticUpperBornRadiusLimit;
}

/**---------------------------------------------------------------------------------------

   Set the quintic spline upper limit

   @param quintic spline upper limit

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setQuinticUpperBornRadiusLimit( RealOpenMM quinticUpperBornRadiusLimit ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GBVISoftcoreParameters::setQuinticUpperBornRadiusLimit:";

   // ---------------------------------------------------------------------------------------

   _quinticUpperBornRadiusLimit  = quinticUpperBornRadiusLimit;
   _quinticUpperSplineLimit      = POW( _quinticUpperBornRadiusLimit, static_cast<RealOpenMM>(-3.0) ); 
   return 0;
}

/**---------------------------------------------------------------------------------------

   Get the quintic upper spline limit

   @return the quintic upper spline limit

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getQuinticUpperSplineLimit( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GBVISoftcoreParameters::getQuinticUpperSplineLimit:";

   // ---------------------------------------------------------------------------------------

   return _quinticUpperSplineLimit;
}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

RealOpenMM* GBVISoftcoreParameters::getAtomicRadii( void ) const {

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

int GBVISoftcoreParameters::setAtomicRadii( RealOpenMM* atomicRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii vector of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGBVISoftcoreParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii );
}

/**---------------------------------------------------------------------------------------

   Return scaled radii
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVISoftcoreParameters::getScaledRadii( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::getScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadii == NULL ){
      GBVISoftcoreParameters* localThis  = const_cast<GBVISoftcoreParameters* const>(this);
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

int GBVISoftcoreParameters::setOwnScaledRadii( int ownScaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownScaledRadii = ownScaledRadii;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setScaledRadii( RealOpenMM* scaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadii && _scaledRadii != scaledRadii ){
      delete[] _scaledRadii;
      _ownScaledRadii = false;
   }

   _scaledRadii = scaledRadii;

   return 0;

}

#if RealOpenMMType == 0

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setScaledRadii( float* scaledRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setScaledRadii";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadii == NULL ){
      _scaledRadii    = new RealOpenMM[getNumberOfAtoms()];
      _ownScaledRadii = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _scaledRadii[ii] = (RealOpenMM) scaledRadii[ii];
   }

   return 0;

}

#endif

/**---------------------------------------------------------------------------------------

   Set scaled radii

   @param scaledRadii  scaledRadii

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setScaledRadii( const RealOpenMMVector& scaledRadii ){

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

   return 0;
}

/**---------------------------------------------------------------------------------------

   Return gamma parameters
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

RealOpenMM* GBVISoftcoreParameters::getGammaParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::getGammaParameters";

   // ---------------------------------------------------------------------------------------

   if( _gammaParameters == NULL ){
      GBVISoftcoreParameters* localThis = const_cast<GBVISoftcoreParameters* const>(this);
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

int GBVISoftcoreParameters::setOwnGammaParameters( int ownGammaParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownGammaParameters = ownGammaParameters;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gamma parameters

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setGammaParameters( RealOpenMM* gammas ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setGammas";

   // ---------------------------------------------------------------------------------------

   if( _ownGammaParameters && _gammaParameters != gammas ){
      delete[] _gammaParameters;
      _ownGammaParameters = false;
   }

   _gammaParameters = gammas;

   return 0;

}

#if RealOpenMMType == 0

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gammas

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setGammaParameters( float* gammas ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setGammas";

   // ---------------------------------------------------------------------------------------

   if( _gammaParameters == NULL ){
      _gammaParameters    = new RealOpenMM[getNumberOfAtoms()];
      _ownGammaParameters = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _gammaParameters[ii] = (RealOpenMM) gammas[ii];
   }

   return 0;

}

#endif

/**---------------------------------------------------------------------------------------

   Set gamma parameters

   @param gammas  gammas

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setGammaParameters( const RealOpenMMVector& gammas ){

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

   return 0;
}

/**---------------------------------------------------------------------------------------

   Return BornRadiusScaleFactors
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

RealOpenMM* GBVISoftcoreParameters::getBornRadiusScaleFactors( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::getBornRadiusScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _bornRadiusScaleFactors == NULL ){
      GBVISoftcoreParameters* localThis = const_cast<GBVISoftcoreParameters* const>(this);
      localThis->_bornRadiusScaleFactors    = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownBornRadiusScaleFactors = true;
      memset( _bornRadiusScaleFactors, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }   
   return _bornRadiusScaleFactors;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether scale factors array should be deleted

   @param ownBornRadiusScaleFactors   flag indicating whether Born radius scale factors 
                                       array should be deleted

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setOwnBornRadiusScaleFactors( int ownBornRadiusScaleFactorsParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownBornRadiusScaleFactors = ownBornRadiusScaleFactorsParameters;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set BornRadiusScaleFactors 

   @param bornRadiusScaleFactors Born radius scale factors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setBornRadiusScaleFactors( RealOpenMM* bornRadiusScaleFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setBornRadiusScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownBornRadiusScaleFactors && _bornRadiusScaleFactors != bornRadiusScaleFactors ){
      delete[] _bornRadiusScaleFactors;
      _ownBornRadiusScaleFactors = false;
   }

   _bornRadiusScaleFactors = bornRadiusScaleFactors;

   return 0;

}

#if RealOpenMMType == 0

/**---------------------------------------------------------------------------------------

   Set bornRadiusScaleFactors parameters

   @param bornRadiusScaleFactorss  bornRadiusScaleFactorss

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setBornRadiusScaleFactors( float* bornRadiusScaleFactorss ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::setBornRadiusScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _bornRadiusScaleFactors == NULL ){
      _bornRadiusScaleFactors    = new RealOpenMM[getNumberOfAtoms()];
      _ownBornRadiusScaleFactors = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _bornRadiusScaleFactors[ii] = (RealOpenMM) bornRadiusScaleFactorss[ii];
   }

   return 0;

}

#endif

/**---------------------------------------------------------------------------------------

   Set bornRadiusScaleFactors parameters

   @param bornRadiusScaleFactors  bornRadiusScaleFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setBornRadiusScaleFactors( const RealOpenMMVector& bornRadiusScaleFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setBornRadiusScaleFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownBornRadiusScaleFactors && _bornRadiusScaleFactors != NULL ){
      delete[] _bornRadiusScaleFactors;
   }
   _ownBornRadiusScaleFactors = true;

   _bornRadiusScaleFactors    = new RealOpenMM[getNumberOfAtoms()];
   for( int ii = 0; ii < (int) bornRadiusScaleFactors.size(); ii++ ){
      _bornRadiusScaleFactors[ii] = bornRadiusScaleFactors[ii];
   }

   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state
   
   @param title               title (optional)
      
   @return string
      
   --------------------------------------------------------------------------------------- */

std::string GBVISoftcoreParameters::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << ImplicitSolventParameters::getStateString( title );

   std::string tab           = getStringTab();

   return message.str();

}

/**---------------------------------------------------------------------------------------
            
   Return zero value if all parameters set; else return nonzero
         
   @return ready status
            
   --------------------------------------------------------------------------------------- */
 
int GBVISoftcoreParameters::isNotReady( void ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGBVISoftcoreParameters::isNotReady";

   // ---------------------------------------------------------------------------------------

   int isReady = ImplicitSolventParameters::isNotReady();

   int errors  = 0;
   std::stringstream message;
   message << methodName;

   const RealOpenMM* scaledRadii = getScaledRadii();
   if( scaledRadii == NULL || scaledRadii[0] <= 0.0 ){
      errors++;
      message << "\n   scaledRadii are not set";
   }

   const RealOpenMM* gamma = getGammaParameters();
   if( gamma == NULL ){
      errors++;
      message << "\n   gamma parameters are not set";
   }

   // check scale factors are in correct units

   RealOpenMM average, stdDev, maxValue, minValue;
   int minIndex, maxIndex;
   SimTKOpenMMUtilities::getArrayStatistics( getNumberOfAtoms(), scaledRadii, &average,
                                             &stdDev, &minValue, &minIndex,
                                             &maxValue, &maxIndex );

   if( average < 0.3 || average > 2.0 || minValue < 0.1 ){
      errors++;
      message << "\n   scale factors for atomic radii appear not to be set correctly -- radii should be in nm";
      message << "\n   average radius=" << average << " min radius=" << minValue << " at atom index=" << minIndex;
   }   

   if( errors ){
      message << std::endl;
      (void) fprintf( stderr, "%s", message.str().c_str() );
   }

   errors += isReady;

   return errors;
}

/**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setUseCutoff( RealOpenMM distance ) {

    cutoff = true;
    cutoffDistance = distance;
    return 0;
}

/**---------------------------------------------------------------------------------------

     Get whether to use a cutoff.

     --------------------------------------------------------------------------------------- */

bool GBVISoftcoreParameters::getUseCutoff() {
    return cutoff;
}

/**---------------------------------------------------------------------------------------

     Get the cutoff distance.

     --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getCutoffDistance() {
    return cutoffDistance;
}

/**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

int GBVISoftcoreParameters::setPeriodic( RealOpenMM* boxSize ) {

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

bool GBVISoftcoreParameters::getPeriodic() {
    return periodic;
}

/**---------------------------------------------------------------------------------------

     Get the periodic box dimension

     --------------------------------------------------------------------------------------- */

const RealOpenMM* GBVISoftcoreParameters::getPeriodicBox() {
    return periodicBoxSize;
}

/**---------------------------------------------------------------------------------------

   Get tau prefactor

   @return (1/e1 - 1/e0), where e1 = solute dielectric, e0 = solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM GBVISoftcoreParameters::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGBVISoftcoreParameters::getTau:";

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

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = " GBVISoftcoreParameters::getBornRadiusScalingSoftcoreMethod:";

   // ---------------------------------------------------------------------------------------

   return _bornRadiusScalingSoftcoreMethod;

}

/**---------------------------------------------------------------------------------------

   Set Born radii switching function method

   @param method

   --------------------------------------------------------------------------------------- */

void GBVISoftcoreParameters::setBornRadiusScalingSoftcoreMethod( BornRadiusScalingSoftcoreMethod bornRadiusScalingSoftcoreMethod ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = " GBVISoftcoreParameters::setBornRadiusScalingSoftcoreMethod:";

   // ---------------------------------------------------------------------------------------

   _bornRadiusScalingSoftcoreMethod = bornRadiusScalingSoftcoreMethod;

}
