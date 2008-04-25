
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

#include <math.h>
#include <sstream>

#include "GrycukParameters.h"
#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"

// #define UseGromacsMalloc 1

#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

const std::string GrycukParameters::ParameterFileName = std::string( "params.agb" );

/**---------------------------------------------------------------------------------------

   GrycukParameters:

		Calculates for each atom

			(1) the van der Waal radii
         (2) volume
         (3) fixed terms in Grycuk equation gPol
         (4) list of atoms that should be excluded in calculating
				 force -- nonbonded atoms (1-2, and 1-3 atoms)

	Implementation:

		Slightly different sequence of calls when running on CPU vs GPU.
		Difference arise because the CPU-side data arrays for the Brook
		streams are allocated by the BrookStreamWrapper objects. These
		arrays are then used by GrycukParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 GrycukParameters* grycukParameters = new GrycukParameters( numberOfAtoms, log );
          grycukParameters->initializeParameters( top );

		Gpu:

			grycukParameters   = new GrycukParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that GrycukParameters destructor does not free arrays 
			
			grycukParameters->setVdwRadii(  getBrookStreamWrapperAtIndex( GpuGrycuk::grycukVdwRadii  )->getData() );
			grycukParameters->setVolume(    getBrookStreamWrapperAtIndex( GpuGrycuk::grycukVolume    )->getData() );
			grycukParameters->setGPolFixed( getBrookStreamWrapperAtIndex( GpuGrycuk::grycukGpolFixed )->getData() );
			grycukParameters->setBornRadii( getBrookStreamWrapperAtIndex( GpuGrycuk::grycukBornRadii )->getData() );
			
			grycukParameters->setFreeArrays( false );
			
			grycukParameters->initializeParameters( top );
 

   Issues:

		Tinker's atom radii are used. 
      The logic for mapping the Gromacs atom names to Tinker type may be incomplete;
      only tested for generic proteins
		see mapGmxAtomNameToTinkerAtomNumber()

   --------------------------------------------------------------------------------------- */


/**---------------------------------------------------------------------------------------

   GrycukParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms

   --------------------------------------------------------------------------------------- */

GrycukParameters::GrycukParameters( int numberOfAtoms ) : ImplicitSolventParameters( numberOfAtoms ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGrycukParameters::GrycukParameters";
   
   // ---------------------------------------------------------------------------------------

   _ownScaledRadiusFactors       = 0;
   _scaledRadiusFactors          = NULL;

}

/**---------------------------------------------------------------------------------------

   GrycukParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

GrycukParameters::~GrycukParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGrycukParameters::~GrycukParameters";
   
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
/*
   if( getFreeArrays() ){

   } */

#endif

}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

RealOpenMM* GrycukParameters::getAtomicRadii( void ) const {

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
   @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                  SimTKOpenMMCommon::MdUnits 

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GrycukParameters::setAtomicRadii( RealOpenMM* atomicRadii, int units ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGrycukParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii, units );
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii vector of atomic radii
   @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                  SimTKOpenMMCommon::MdUnits 

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int GrycukParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii, int units ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGrycukParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   return ImplicitSolventParameters::setAtomicRadii( atomicRadii, units );
}

/**---------------------------------------------------------------------------------------

   Return Grycuk scale factors
   If not previously set, allocate space

   @return array 

   --------------------------------------------------------------------------------------- */

const RealOpenMM* GrycukParameters::getScaledRadiusFactors( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::getScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadiusFactors == NULL ){
      GrycukParameters* localThis = const_cast<GrycukParameters* const>(this);
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

int GrycukParameters::setOwnScaleFactors( int ownScaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::setOwnScaleFactors";

   // ---------------------------------------------------------------------------------------

   _ownScaledRadiusFactors = ownScaledRadiusFactors;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set Grycuk scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GrycukParameters::setScaledRadiusFactors( RealOpenMM* scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadiusFactors && _scaledRadiusFactors != scaledRadiusFactors ){
      delete[] _scaledRadiusFactors;
      _ownScaledRadiusFactors = false;
   }

   _scaledRadiusFactors = scaledRadiusFactors;

   return SimTKOpenMMCommon::DefaultReturn;

}

#if RealOpenMMType == 2

/**---------------------------------------------------------------------------------------

   Set Grycuk scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GrycukParameters::setScaledRadiusFactors( float* scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _scaledRadiusFactors == NULL ){
      _scaledRadiusFactors    = new RealOpenMM[getNumberOfAtoms()];
      _ownScaledRadiusFactors = true;
   }   
   for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
      _scaledRadiusFactors[ii] = (RealOpenMM) scaledRadiusFactors[ii];
   }

   return SimTKOpenMMCommon::DefaultReturn;

}

#endif

/**---------------------------------------------------------------------------------------

   Set Grycuk scale factors

   @param scaledRadiusFactors  scaledRadiusFactors

   @return SimTKOpenMMCommon::DefaultReturn always

   --------------------------------------------------------------------------------------- */

int GrycukParameters::setScaledRadiusFactors( const RealOpenMMVector& scaledRadiusFactors ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGrycuk::setScaledRadiusFactors";

   // ---------------------------------------------------------------------------------------

   if( _ownScaledRadiusFactors && _scaledRadiusFactors != NULL ){
      delete[] _scaledRadiusFactors;
   }
   _ownScaledRadiusFactors = true;
   _scaledRadiusFactors    = new RealOpenMM[getNumberOfAtoms()];
   for( int ii = 0; ii < (int) scaledRadiusFactors.size(); ii++ ){
      _scaledRadiusFactors[ii] = scaledRadiusFactors[ii];
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Map Gmx atom name to Tinker atom number (Simbios)

   @param atomName            atom name (CA, HA, ...); upper and lower case should both work
   @param log                 if set, then print error messages to log file

   @return Tinker atom number if atom name is valid; else return -1

   --------------------------------------------------------------------------------------- */

int GrycukParameters::mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const {

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

std::string GrycukParameters::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGrycukParameters::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << ImplicitSolventParameters::getStateString( title );

   // std::string tab           = getStringTab();

   return message.str();


}

/**---------------------------------------------------------------------------------------
            
   Return zero value if all parameters set; else return nonzero
         
   @return ready status
            
   --------------------------------------------------------------------------------------- */
 
int GrycukParameters::isNotReady( void ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGrycukParameters::isNotReady";

   // ---------------------------------------------------------------------------------------

   int isReady = ImplicitSolventParameters::isNotReady();

   int errors  = 0;
   std::stringstream message;
   message << methodName;

   const RealOpenMM* scaledRadiusFactors = getScaledRadiusFactors();
   if( scaledRadiusFactors == NULL || scaledRadiusFactors[0] <= 0.0 ){
      errors++;
      message << "\n   scaledRadiusFactors is not set";
   }

   // check scale factors are in correct units

   RealOpenMM average, stdDev, maxValue, minValue;
   int minIndex, maxIndex;
   SimTKOpenMMUtilities::getArrayStatistics( getNumberOfAtoms(), scaledRadiusFactors, &average,
                                             &stdDev, &minValue, &minIndex,
                                             &maxValue, &maxIndex );

   if( average < 0.3 || average > 2.0 || minValue < 0.1 ){
      errors++;
      message << "\n   scale factors for atomic radii appear not to be set correctly -- radii should be in Angstroms";
      message << "\n   average radius=" << average << " min radius=" << minValue << " at atom index=" << minIndex;
   }   

   if( errors ){
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   errors += isReady;

   return errors;
}

