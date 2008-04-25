
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

// generic c includes

#include <math.h>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMGromacsUtilities.h"
#include "GbsaAtomParameter.h"
#include "GbsaParameters.h"

#define UseGromacsMalloc 1

#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif


/**---------------------------------------------------------------------------------------

   GbsaParameters:

		Calculates for each atom

			(1) the van der Waal radii
         (2) volume
         (3) fixed terms in Gbsa equation gPol
         (4) list of atoms that should be excluded in calculating
				 force -- nonbonded atoms (1-2, and 1-3 atoms)

	Implementation:

		Slightly different sequence of calls when running on CPU vs GPU.
		Difference arise because the CPU-side data arrays for the Brook
		streams are allocated by the BrookStreamWrapper objects. These
		arrays are then used by GbsaParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 GbsaParameters* gbsaParameters = new GbsaParameters( numberOfAtoms, log );
          gbsaParameters->initializeParameters( top );

		Gpu:

			gbsaParameters   = new GbsaParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that GbsaParameters destructor does not free arrays 
			
			gbsaParameters->setVdwRadii(  getBrookStreamWrapperAtIndex( GpuGbsa::gbsaVdwRadii  )->getData() );
			gbsaParameters->setVolume(    getBrookStreamWrapperAtIndex( GpuGbsa::gbsaVolume    )->getData() );
			gbsaParameters->setGPolFixed( getBrookStreamWrapperAtIndex( GpuGbsa::gbsaGpolFixed )->getData() );
			gbsaParameters->setBornRadii( getBrookStreamWrapperAtIndex( GpuGbsa::gbsaBornRadii )->getData() );
			
			gbsaParameters->setFreeArrays( false );
			
			gbsaParameters->initializeParameters( top );
 

   Issues:

		Tinker's atom radii are used. 
      The logic for mapping the Gromacs atom names to Tinker type may be incomplete;
      only tested for generic proteins
		see mapGmxAtomNameToTinkerAtomNumber()

   --------------------------------------------------------------------------------------- */


/**---------------------------------------------------------------------------------------

   GbsaParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms
   @param log                 file descriptor (can be null)

   --------------------------------------------------------------------------------------- */

GbsaParameters::GbsaParameters( int numberOfAtoms, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::GbsaParameters";
   
   // ---------------------------------------------------------------------------------------

   _numberOfAtoms          = numberOfAtoms;
   _log                    = log;

   _gbsaBondsArray         = NULL;
   _exclusionWorkArray     = NULL;

   _vdwRadii               = NULL;
   _volume                 = NULL;
   _gPolFixed              = NULL;
   _bornRadii              = NULL;

   // see comments in ~GbsaParameters for explanation

   _freeArrays             = true;

}

/**---------------------------------------------------------------------------------------

   GbsaParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

GbsaParameters::~GbsaParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::~GbsaParameters";
   
   // ---------------------------------------------------------------------------------------

   // in GPU runs, arrays may be 'owned' by BrookStreamWrapper -- hence they should not
   // be freed here, i.e., _freeArrays should be 'false'   

#ifdef UseGromacsMalloc

   if( _freeArrays ){

      if( _vdwRadii != NULL ){
         save_free( "_vdwRadii", __FILE__, __LINE__, _vdwRadii );
      }
   
      if( _volume != NULL ){
         save_free( "_volume", __FILE__, __LINE__, _volume );
      }
   
      if( _gPolFixed != NULL ){
         save_free( "_gPolFixed", __FILE__, __LINE__, _gPolFixed );
      }
   
      if( _bornRadii != NULL ){
         save_free( "_bornRadii", __FILE__, __LINE__, _bornRadii );
      }

   }

   if( _exclusionWorkArray != NULL ){
      save_free( "_exclusionWorkArray", __FILE__, __LINE__, _exclusionWorkArray );
   }

#else

   if( _freeArrays ){

      if( _vdwRadii != NULL ){
         delete[] _vdwRadii;
      }
   
      if( _volume != NULL ){
         delete[] _volume;
      }
   
      if( _gPolFixed != NULL ){
         delete[] _gPolFixed;
      }
   
      if( _bornRadii != NULL ){
         delete[] _bornRadii;
      }

   }

   if( _exclusionWorkArray != NULL ){
      delete[] _exclusionWorkArray;
   }

#endif

   if( _gbsaBondsArray != NULL ){
      freeBondArray( _numberOfAtoms, _gbsaBondsArray );
   }

}

/**---------------------------------------------------------------------------------------

   Initialize Gbsa Parameters (Simbios) 

   --------------------------------------------------------------------------------------- */

void GbsaParameters::initializeConstants( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::initializeConstants:";

   // ---------------------------------------------------------------------------------------

   // GBSA parameters from '97 paper
   
   _innerDielectric         =    1.0f;
   _solventDielectric       =   78.3f;
   _kcalA_To_kJNm           =    4.184f*0.1f;
   _phi                     =   -0.09f;
   _P1                      =    0.073f;
   _P2                      =    0.921f;
   _P3                      =    6.211f;
   _P4                      =   15.236f;
   _P4_2                    =    0.5f*sqrt( _P4 );
   _P5                      =    1.254f;
   _probeRadius             =    0.14f;
   _electricConstant        = -166.02691f;
   
   // ---------------------------------------------------------------------------------------
   
   _resetPreFactor();
   
   _P5Inverse               = 1.0f/_P5;
   _piP5                    = ((Real) M_PI)*_P5;
   _pi4Asolv                = ((Real) M_PI)*4.0f*0.0049f*1000.0f;

}

/**---------------------------------------------------------------------------------------

   Reset prefactor (Simbios) 

	called when _electricConstant, _innerDielectric, or _solventDielectric are modified

   --------------------------------------------------------------------------------------- */

void GbsaParameters::_resetPreFactor( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::_resetPreFactor:";

   // ---------------------------------------------------------------------------------------

   if( _innerDielectric != 0.0f && _solventDielectric != 0.0f ){
      _preFactor = 2.0f*_electricConstant*( (1.0f/_innerDielectric) - (1.0f/_solventDielectric) );
   } else {
      _preFactor = 0.0;
   }
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric (Simbios) 

   @param solventDielectric			solvent dielectric

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::setSolventDielectric( Real solventDielectric ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::setSolventDielectric:";

   // ---------------------------------------------------------------------------------------

   _solventDielectric = solventDielectric;
   _resetPreFactor();

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set inner dielectric (Simbios) 

   @param innerDielectric:			inner dielectric

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::setInnerDielectric( Real innerDielectric ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::setInnerDielectric:";

   // ---------------------------------------------------------------------------------------

   _innerDielectric = innerDielectric;
   _resetPreFactor();

   return 0;

}

/**---------------------------------------------------------------------------------------

   Set electric constant (Simbios) 

   @param electricConstant			electric constant

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::setElectricConstant( Real electricConstant ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::setElectricConstant:";

   // ---------------------------------------------------------------------------------------

   _electricConstant = electricConstant;
   _resetPreFactor();

   return 0;
}

/**---------------------------------------------------------------------------------------

   Initialize Gbsa Parameters (Simbios) 

   @param top                 GMX topology data struct -- used to get harmonic angle

   @return 0 if no errors

   --------------------------------------------------------------------------------------- */

int GbsaParameters::initializeParameters( const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   bool printParameters          = false;
   bool debugOn                  = false;
   // bool debugOn                  = false;
   static const char* methodName = "\nGbsaParameters::initializeParameters:";

   // ---------------------------------------------------------------------------------------
   
   FILE* log = getLog();

   initializeConstants();

   // find stretch bond entries via F_BONDS & F_SHAKE (for H's) entries
   // Settle 'bonds' handled separtely

   gbsaBonds** gbsaBondsArray = (gbsaBonds**) allocateBondsArray( getNumberOfAtoms() );
   int errors                 = addStretchBonds( getNumberOfAtoms(), gbsaBondsArray, top, F_BONDS,  3, 1, log );
                   errors    += addStretchBonds( getNumberOfAtoms(), gbsaBondsArray, top, F_SHAKE,  3, 1, log );

                   errors    += addSettleStretchBonds( getNumberOfAtoms(), gbsaBondsArray, top, F_SETTLE, 2, 1, log );

   if( errors > 0 ){
      exit(-1);
   } else if( log != NULL ){
      (void) fprintf( log, "\nNo 1-2 bond errors." );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   // find angle bond entries based on stretch bond entries

   errors = addAngleBonds( getNumberOfAtoms(), gbsaBondsArray, top, F_ANGLES, log );
   setGbsaBondsArray( gbsaBondsArray );

   if( errors > 0 ){
      printBondsArray( getNumberOfAtoms(), top, gbsaBondsArray, "Bonds", log );
      exit(-1);
   } else if( log != NULL ){
      (void) fprintf( log, "\nNo 1-3 bond errors." );
      (void) fflush( log );
   }

   // set exclusion array

   setGbsaExclusions( top, debugOn, log );
   if( log != NULL ){
      (void) fprintf( log, "%s exclusions set", methodName );
      (void) fflush( log );
   }

   // print bonds and report if errors

   if( debugOn && log ){
      printBondsArray( getNumberOfAtoms(), top, gbsaBondsArray, "Bonds", log );
   }

   // ---------------------------------------------------------------------------------------

   // allocate memory for work arrays, if not already allocated

   // when Gbsa is run on GPU, the memory for these arrays may
   // have been allocated when the BrookStreamWrappers 
   // were created (GpuGbsa constructor) and hence do not
   // need to be allocated here

   unsigned int arraySz = getNumberOfAtoms()*sizeof( Real );

#ifdef UseGromacsMalloc
   if( _vdwRadii == NULL ){
      _vdwRadii = (Real*) save_malloc( "_vdwRadii", __FILE__, __LINE__, arraySz );
   }

   if( _volume == NULL ){
      _volume = (Real*) save_malloc( "_volume", __FILE__, __LINE__, arraySz );
   }

   if( _gPolFixed == NULL ){
      _gPolFixed = (Real*) save_malloc( "_gPolFixed", __FILE__, __LINE__, arraySz );
   }

   if( _bornRadii == NULL ){
      _bornRadii = (Real*) save_malloc( "_bornRadii", __FILE__, __LINE__, arraySz );
   }

#else

   if( _vdwRadii == NULL ){
      _vdwRadii  = new Real[getNumberOfAtoms()];
   }

   if( _volume == NULL ){
      _volume    = new Real[getNumberOfAtoms()];
   }

   if( _gPolFixed == NULL ){
      _gPolFixed = new Real[getNumberOfAtoms()];
   }

   if( _bornRadii == NULL ){
      _bornRadii = new Real[getNumberOfAtoms()];
   }

#endif

   memset( _vdwRadii, 0, arraySz );
   memset( _volume, 0, arraySz );
   memset( _gPolFixed, 0, arraySz );
   memset( _bornRadii, 0, arraySz );

   // compute atomic volumes (eqs 6 & 7) in J. Phys. Chem. A V101 No 16, p. 3005
   // In Tinker: ksolv.f ~ line 340 (solvtyp.eq.'STILL')

   // get Macromodel vdW radii and volumes (excluding overlap)
   // and calculate fixed GBSA terms

   std::string parameterFileName = "Params_macromodel_agb.txt";
fprintf( log, "\nCalling  getMacroModelAtomicRadii %s", parameterFileName.c_str() );
   getMacroModelAtomicRadii( getNumberOfAtoms(), parameterFileName, top, top->atoms.atomname, _vdwRadii, log );
/*
   getMacroModelAtomicRadii( getNumberOfAtoms(), gbsaBondsArray, top->atoms.atomname, _bornRadii, log );
char*** names = top->atoms.atomname;
for( int ii = 0; ii < getNumberOfAtoms(); ii++ ){
   Real diff = abs( _bornRadii[ii] - _vdwRadii[ii] );
   fprintf( log, "\n%d %s %.3f %3f %.3f", ii, (*(names[ii])), diff, _vdwRadii[ii], _bornRadii[ii] );
}
fprintf( log, "\nDone getMacroModelAtomicRadii %s", parameterFileName.c_str() );
fflush( log );
exit(0);
*/
   getMacroModelAtomicVolumes( getNumberOfAtoms(), top, gbsaBondsArray, _vdwRadii, _volume, log );
   getFixedGBSA_GPol( getNumberOfAtoms(), top, gbsaBondsArray, _vdwRadii, _volume, _gPolFixed, log );

   // ---------------------------------------------------------------------------------------

   // print parameters

   if( printParameters && log != NULL ){
      printParameterInfo( getNumberOfAtoms(), top,
                          "GmxParameters.txt", log );
   }

   // ---------------------------------------------------------------------------------------

   if( log != NULL ){
      (void) fprintf( log, "%s done", methodName );
      (void) fflush( log );
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Set Gbsa Exclusions (Simbios) 

   @param top                 GMX topology data struct -- used to get harmonic angle
   @param printOn             print diagnostics 
   @param log                 file descriptor 

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::setGbsaExclusions( const t_topology* top, bool printOn, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const int  collectionBinMax = 500;
   int collectionBin[collectionBinMax];

   bool saveExclusionIndices = true;
   // static const char* methodName = "\nGbsaParameters::setGbsaExclusions";
   
   // ---------------------------------------------------------------------------------------

   int numberOfAtoms          = getNumberOfAtoms();
   gbsaBonds** gbsaBondsArray = getGbsaBondsArray();

   // excludedAtoms: used to mark which atoms are to excluded

#ifdef UseGromacsMalloc
   int* excludedAtoms = (int*) save_malloc( "excludedAtoms", __FILE__, __LINE__, (numberOfAtoms + 2)*sizeof( int ) );
#else
   int* excludedAtoms = new int[numberOfAtoms+2];
#endif

   for( int ii = 0; ii < numberOfAtoms + 2; ii++ ){
      excludedAtoms[ii]   = -1;
   }

   // gather indices to be excluded for atom ii and then
   // set values in excludeIndices array
   // excludedAtoms is used to insure each atom is only 
   // excluded once

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      
      if( gbsaBondsArray[ii] != NULL ){

         int numberOfIndicesToExclude =
            collectGbsaBondIndices( gbsaBondsArray[ii], collectionBinMax, collectionBin, ii, log );

         gbsaBondsArray[ii]->minExclusionIndex = numberOfAtoms + 1;
         gbsaBondsArray[ii]->maxExclusionIndex = -1;

// SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( ii, top, MAX_BUFFER_SZ, buffer, numberOfAtoms, ATOM_ID_STRING_TAB );
// (void) fprintf( log, "\n[%d %s] CalcExcnt=%d ", ii, buffer, numberOfIndicesToExclude );

         for( int jj = 0; jj < numberOfIndicesToExclude; jj++ ){
            int excludeIndex = collectionBin[jj];

            if( excludedAtoms[excludeIndex] != ii ){

                excludedAtoms[excludeIndex] = ii;

               if( saveExclusionIndices && gbsaBondsArray[ii]->numberOfExclusions < MAX_BOND_EXCLUSIONS ){
                  gbsaBondsArray[ii]->exclusions[gbsaBondsArray[ii]->numberOfExclusions++] = excludeIndex;
// (void) fprintf( log, "[%d %d] ", excludeIndex, gbsaBondsArray[ii]->numberOfExclusions );
               }

               if( gbsaBondsArray[ii]->minExclusionIndex > excludeIndex ){
                  gbsaBondsArray[ii]->minExclusionIndex = excludeIndex;
               }
               if( gbsaBondsArray[ii]->maxExclusionIndex < excludeIndex ){
                  gbsaBondsArray[ii]->maxExclusionIndex = excludeIndex;
               }

            }
         }

      }

   }

   // free excludedAtoms array

#ifdef UseGromacsMalloc
   save_free( "excludedAtoms", __FILE__, __LINE__, excludedAtoms );
#else
   delete[] excludedAtoms;
#endif

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Collect Gbsa bond indices for a particular atom (Simbios) 

   @param gbsaBond            gbsaBond 
   @param collectionBin       array of indices to be excluded
   @param excludeIndex        exclude index

   @return number of collected indices (includes excludeIndex at end of array)

   --------------------------------------------------------------------------------------- */

int GbsaParameters::collectGbsaBondIndices( gbsaBonds* gbsaBond, int collectionBinMax,
                                            int* collectionBin, int excludeIndex, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const int  gbsaStretchBondArrayMax = 500;
   gbsaStretchBond* gbsaStretchBondArray[gbsaStretchBondArrayMax];

   static const char* methodName = "\nGbsaParameters:;collectGbsaBondIndices";
   
   // ---------------------------------------------------------------------------------------
 
   // first collect all stretch bonds from the stretchBonds and angleBonds data structs

   int numberOfStretchBonds = 0;
   gbsaStretchBond* stretchBonds = gbsaBond->stretchBonds;
   while( stretchBonds != NULL && numberOfStretchBonds < gbsaStretchBondArrayMax ){
      gbsaStretchBondArray[numberOfStretchBonds++] = stretchBonds; 
      stretchBonds = stretchBonds->nextBond;
   }

/*
   (void) fprintf( log, "\n   excludeIndex=%d str=%d", excludeIndex, numberOfStretchBonds );
   (void) fflush( log );
*/

   gbsaAngleBond* angleBonds = gbsaBond->angleBonds;
   while( angleBonds != NULL && numberOfStretchBonds < (gbsaStretchBondArrayMax-1) ){

      if( angleBonds->stretchBondI != NULL ){
         gbsaStretchBondArray[numberOfStretchBonds++] = angleBonds->stretchBondI;
      }
      if( angleBonds->stretchBondJ != NULL ){
         gbsaStretchBondArray[numberOfStretchBonds++] = angleBonds->stretchBondJ;
      }
      angleBonds = angleBonds->nextBond;
   }

   if( numberOfStretchBonds >= gbsaStretchBondArrayMax ){
      (void) fprintf( log ? log : stderr, "%s gbsaStretchBondArrayMax=%d too small", methodName, gbsaStretchBondArrayMax );
      return -1;
   }

/*
   (void) fprintf( log, " str2=%d", numberOfStretchBonds );
   (void) fflush( log );
*/

   // ---------------------------------------------------------------------------------------

   // now collect all indices including the exclude index

   int collectionBinSize = 0;
   for( int ii = 0; ii < numberOfStretchBonds && collectionBinSize < (collectionBinMax+2); ii++ ){

/*
   (void) fprintf( log, " %s", (gbsaStretchBondArray[ii] == NULL) ? "NULL" : "Ok" );
   (void) fflush( log );
*/

      if( gbsaStretchBondArray[ii]->atomI != excludeIndex ){
         collectionBin[collectionBinSize++] = gbsaStretchBondArray[ii]->atomI;
      }

      if( gbsaStretchBondArray[ii]->atomJ != excludeIndex ){
         collectionBin[collectionBinSize++] = gbsaStretchBondArray[ii]->atomJ;
      }

/*
   (void) fprintf( log, " collect=[%d %d %d]", collectionBinSize,  gbsaStretchBondArray[ii]->atomI,  gbsaStretchBondArray[ii]->atomJ );
   (void) fflush( log );
*/

   }
   if( collectionBinSize < collectionBinMax ){
      collectionBin[collectionBinSize++] = excludeIndex;
   } else {
      (void) fprintf( log ? log : stderr, "%s collectionBinMax=%d is too small.", methodName, collectionBinMax );
   } 
 
   return collectionBinSize;
}

/**---------------------------------------------------------------------------------------

   Allocate memory for bond array (Simbios) 

   @param maxAtoms max number of atoms

   array entries are initialized to zero

   free array[0] and then array when done

   @returns ptr to allocated array or NULL if out of memory

   --------------------------------------------------------------------------------------- */

gbsaBonds** GbsaParameters::allocateBondsArray( int maxAtoms ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGbsaParameters::allocateBondsArray";

   // ---------------------------------------------------------------------------------------

#ifdef UseGromacsMalloc
   gbsaBonds** gbsaBondsArray      = (gbsaBonds**) save_malloc( "gbsaBondsArray", __FILE__, __LINE__, maxAtoms*sizeof( gbsaBonds* ) );
   gbsaBonds*  gbsaBondsBlock      = (gbsaBonds* ) save_malloc( "gbsaBondsBlock", __FILE__, __LINE__, maxAtoms*sizeof( gbsaBonds  ) );
#else
   gbsaBonds** gbsaBondsArray      = new gbsaBonds*[maxAtoms];
   gbsaBonds*  gbsaBondsBlock      = new gbsaBonds[maxAtoms];
#endif

   if( gbsaBondsArray == NULL || gbsaBondsBlock == NULL ){
      (void) fprintf( stderr, "%s apparently out of memory: maxAtoms=%d", methodName, maxAtoms );
      (void) fflush( stderr );
      return NULL;
   }
   memset( gbsaBondsBlock, 0, sizeof( gbsaBonds )*maxAtoms );

   for( int ii = 0; ii < maxAtoms; ii++ ){
      gbsaBondsArray[ii] = gbsaBondsBlock++;
   }

   return gbsaBondsArray;
}

/**---------------------------------------------------------------------------------------

   Deallocate memory for bond array (Simbios) 

   @param maxAtoms  max number of atoms
   @param gbsaBondsArray   array to be freed

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::freeBondArray( int maxAtoms, gbsaBonds** gbsaBondsArray ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::freeBondArray";

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < maxAtoms; ii++ ){

      if( gbsaBondsArray[ii] != NULL ){

         gbsaStretchBond* nextStretch;
         gbsaStretchBond* stretchBonds = gbsaBondsArray[ii]->stretchBonds;
         while( stretchBonds != NULL ){
            nextStretch = stretchBonds->nextBond;

#ifdef UseGromacsMalloc
            save_free( "stretchBonds", __FILE__, __LINE__, stretchBonds );
#else
            delete stretchBonds;
#endif

            stretchBonds = nextStretch;
         }

         gbsaAngleBond* angleBonds = gbsaBondsArray[ii]->angleBonds;
         gbsaAngleBond* nextAngle;
         while( angleBonds != NULL ){

            nextAngle = angleBonds->nextBond;

#ifdef UseGromacsMalloc
            save_free( "angleBonds", __FILE__, __LINE__, angleBonds );
#else
            delete angleBonds;
#endif
            angleBonds = nextAngle;
         }
      }
   }

#ifdef UseGromacsMalloc
   save_free( "gbsaBondsArray[0]", __FILE__, __LINE__, gbsaBondsArray[0] );
   save_free( "gbsaBondsArray", __FILE__, __LINE__, gbsaBondsArray );
#else
   delete[] gbsaBondsArray[0];
   delete[] gbsaBondsArray;
#endif

   return 0;
}

/**---------------------------------------------------------------------------------------

   Print bond array (Simbios) 

   @param numberOfAtoms       number of atoms
   @param top                 GMX topology data struct
   @param gbsaBonds           bond array
   @param title               title string (optional)
   @param log                 file descriptor 

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::printBondsArray( int numberOfAtoms, const t_topology* top,
                                     gbsaBonds** gbsaBondsArray, const char* title, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::printBondsArray";

   // ---------------------------------------------------------------------------------------

   if( !log ){
      return 0;
   }

   // print title (if set) and bonds

   if( title != NULL ){
      (void) fprintf( log, "\n%s", title );
   }

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      if( gbsaBondsArray[ii] != NULL ){
         printBond( numberOfAtoms, top, ii, gbsaBondsArray[ii], log );
      }
   }

   (void) fflush( log );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Print bond (Simbios) 

   @param top                 GMX topology data struct -- used to get harmonic angle
   @param atomIndex           atom index
   @param gbsaBond            bond data struct
   @param log                 file descriptor 

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::printBond( int numberOfAtoms, const t_topology* top, int atomIndex,
                               const gbsaBonds* gbsaBond, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   int stretchCount;
   int angleCount;
   char buffer[4][MAX_BUFFER_SZ];
   // static const char* methodName = "\nGbsaParameters::printBond";

   // ---------------------------------------------------------------------------------------

   // print bond

   if( log == NULL ){
      return 0;
   }

   if( gbsaBond != NULL ){

      getBondCounts( gbsaBond, &stretchCount, &angleCount );
      SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atomIndex, top, MAX_BUFFER_SZ, buffer[0], numberOfAtoms, ATOM_ID_STRING_TAB );

      (void) fprintf( log, "\n%d %s StretchCount=%d AngleCount=%d ExclRange=[%d %d] ExclWidth=%d Save=[%d %d] width=%d", 
                      atomIndex, buffer[0], stretchCount, angleCount,
                      gbsaBond->minExclusionIndex,  gbsaBond->maxExclusionIndex,
                      gbsaBond->maxExclusionIndex - gbsaBond->minExclusionIndex,
                      gbsaBond->startExclusionIndex,
                      gbsaBond->stopExclusionIndex,
                      gbsaBond->stopExclusionIndex - gbsaBond->startExclusionIndex
                    );
      if( gbsaBond->numberOfExclusions ){
          (void) fprintf( log, "\nExclusions: %d [ ", gbsaBond->numberOfExclusions ); 
          for( int ii = 0; ii < gbsaBond->numberOfExclusions; ii++ ){
             (void) fprintf( log, "%d ", gbsaBond->exclusions[ii] ); 
          }
          (void) fprintf( log, "]", gbsaBond->numberOfExclusions ); 
      }

      gbsaStretchBond* stretchBonds = gbsaBond->stretchBonds;
      int count = 0;
      while( stretchBonds != NULL ){
         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBonds->atomJ, top, MAX_BUFFER_SZ, buffer[0], 
                                                        numberOfAtoms, ATOM_ID_STRING_TAB );
         (void) fprintf( log, "\n   %d %s %.6f", count++, buffer[0], stretchBonds->bondLength );
         stretchBonds = stretchBonds->nextBond;
      }

      gbsaAngleBond* angleBonds = gbsaBond->angleBonds;
      count                    = 0;
      while( angleBonds != NULL ){

         gbsaStretchBond* stretchBondI = angleBonds->stretchBondI;
         gbsaStretchBond* stretchBondJ = angleBonds->stretchBondJ;

         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBondI->atomI, top, MAX_BUFFER_SZ, buffer[0],
                                                        numberOfAtoms, ATOM_ID_STRING_TAB );

         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBondI->atomJ, top, MAX_BUFFER_SZ, buffer[1],
                                                        numberOfAtoms, ATOM_ID_STRING_TAB );

         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBondJ->atomI, top, MAX_BUFFER_SZ, buffer[2],
                                                        numberOfAtoms, ATOM_ID_STRING_TAB );

         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBondJ->atomJ, top, MAX_BUFFER_SZ, buffer[3],
                                                        numberOfAtoms, ATOM_ID_STRING_TAB );

         (void) fprintf( log, "\n   %d [ %s %s %s %s] [%.4f %.4f] pvt=%d vol=%d Ang=%.2f Ord=[%d %d %d]",
                         count++, buffer[0], buffer[1], buffer[2], buffer[3],
                         stretchBondI->bondLength,
                         stretchBondJ->bondLength,
                         angleBonds->pivotAtomIndex,
                         angleBonds->volumeIndex,
                         angleBonds->harmonicAngleWidth,
                         angleBonds->orderedIndices[0],
                         angleBonds->orderedIndices[1],
                         angleBonds->orderedIndices[2] 
                       );
         angleBonds = angleBonds->nextBond;
      }

   }

   (void) fflush( log );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Print parameter info (Simbios) 

   @param numberOfAtoms				number of atoms
   @param top							GMX topology data struct -- used to get harmonic angle
   @param parameterInfoFileName	parameterInfo FileName
   @param log							file descriptor 

   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::printParameterInfo( int numberOfAtoms, const t_topology* top,
                                        const char* parameterInfoFileName, FILE* log ){

   // ---------------------------------------------------------------------------------------

   char buffer[4][MAX_BUFFER_SZ];
   static const char* methodName = "\nGbsaParameters::printParameterInfo";

   // ---------------------------------------------------------------------------------------

   gbsaBonds** gbsaBondsArray = getGbsaBondsArray();
   FILE* parameterInfoFile = NULL;
#ifdef WIN32
   fopen_s( &parameterInfoFile, parameterInfoFileName, "w" );
#else
   parameterInfoFile = fopen( parameterInfoFileName, "w" );
#endif

   if( parameterInfoFile != NULL ){
      if( log != NULL ){
         (void) fprintf( log, "%s opened file=<%s>.", methodName, parameterInfoFileName );
         (void) fflush( log );
      }
   } else {
      if( log != NULL ){
         (void) fprintf( log, "%s could not open file=<%s> -- abort output.",
                         methodName, parameterInfoFile );
         (void) fflush( log );
      }
      return -1;
   }

   (void) fprintf( parameterInfoFile, "# %d atoms\n", numberOfAtoms );

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( ii, top, MAX_BUFFER_SZ, buffer[0], numberOfAtoms, ATOM_ID_STRING_TAB );
      (void) fprintf( parameterInfoFile, "%d %s\n", ii, buffer[0] );
   }
   (void) fflush( parameterInfoFile );

   (void) fprintf( parameterInfoFile, "StretchBond\n" );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      if( gbsaBondsArray[ii] != NULL ){
         gbsaStretchBond* stretchBonds = gbsaBondsArray[ii]->stretchBonds;
         while( stretchBonds != NULL ){

            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBonds->atomI, top, MAX_BUFFER_SZ, buffer[0],
                                                           numberOfAtoms, ATOM_ID_STRING_TAB );

            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( stretchBonds->atomJ, top, MAX_BUFFER_SZ, buffer[1],
                                                           numberOfAtoms, ATOM_ID_STRING_TAB );

            (void) fprintf( parameterInfoFile, "%d %d %.6e %s %s\n",
                            stretchBonds->atomI, stretchBonds->atomJ, stretchBonds->bondLength,
                             buffer[0], buffer[1] );
            stretchBonds = stretchBonds->nextBond;
         }
      }
   }
   (void) fflush( parameterInfoFile );

   (void) fprintf( parameterInfoFile, "AngleBond\n" );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      if( gbsaBondsArray[ii] != NULL ){
         gbsaAngleBond* angleBonds = gbsaBondsArray[ii]->angleBonds;
         while( angleBonds != NULL ){

            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( angleBonds->orderedIndices[0], top, MAX_BUFFER_SZ, buffer[0],
                                                           numberOfAtoms, ATOM_ID_STRING_TAB );

            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( angleBonds->orderedIndices[1], top, MAX_BUFFER_SZ, buffer[1],
                                                           numberOfAtoms, ATOM_ID_STRING_TAB );

            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( angleBonds->orderedIndices[2], top, MAX_BUFFER_SZ, buffer[2],
                                                           numberOfAtoms, ATOM_ID_STRING_TAB );

            (void) fprintf( parameterInfoFile, "%d %d %d %.6e %s %s %s\n",
                            angleBonds->orderedIndices[0],
                            angleBonds->orderedIndices[1],
                            angleBonds->orderedIndices[2],
                            angleBonds->harmonicAngleWidth,  buffer[0],  buffer[1],  buffer[2] );
            angleBonds = angleBonds->nextBond;
         }
      }
   }

   (void) fclose( parameterInfoFile );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Add angle bonds using stretch bond lists (Simbios) 
   Also load in angle width

   @param maxAtoms            max number of atoms
   @param gbsaBondsArray      array of gbsaBonds
   @param top                 GMX topology data struct -- used to get harmonic angle
   @param idefArrayIndex      parameter index for GMX iparams data struct
   @param log                 if set, then print error messages to log file

   @return 0 if no errors

   --------------------------------------------------------------------------------------- */

int GbsaParameters::addAngleBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                                   const t_topology* top, int idefArrayIndex, FILE* log  ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGbsaParameters::addAngleBonds";

   // ---------------------------------------------------------------------------------------

   // loop over atoms
   //   for each atomI loop over stretch bonds
   //      for each atomJ in stretch bond, check if atomI/atomJ != ii
   //         create new angle bond w/ atomJ as the pivot atom

   for( int ii = 0; ii < maxAtoms; ii++ ){

      gbsaStretchBond* stretchBonds = gbsaBondsArray[ii]->stretchBonds;
      while( stretchBonds != NULL ){

         int atomJ = stretchBonds->atomJ;
         gbsaStretchBond* stretchBondsAtomJ = gbsaBondsArray[atomJ]->stretchBonds;

         while( stretchBondsAtomJ != NULL ){
            if( stretchBondsAtomJ->atomI != ii && stretchBondsAtomJ->atomJ != ii ){

               addStretchBondsToAngleBondList( stretchBonds, stretchBondsAtomJ,
                                               atomJ, gbsaBondsArray, ii, log );

               if( stretchBondsAtomJ->atomI != atomJ ){
                  addStretchBondsToAngleBondList( stretchBonds, stretchBondsAtomJ, 
                                                  atomJ, gbsaBondsArray, stretchBondsAtomJ->atomI, log );
               }
               if( stretchBondsAtomJ->atomJ != atomJ ){
                  addStretchBondsToAngleBondList( stretchBonds, stretchBondsAtomJ, 
                                                  atomJ, gbsaBondsArray, stretchBondsAtomJ->atomJ, log );
               }
            }
            stretchBondsAtomJ = stretchBondsAtomJ->nextBond;
         }
         stretchBonds = stretchBonds->nextBond;
      }
   }

   // printBondsArray( maxAtoms, top, gbsaBondsArray, "PreAngleWidth", log );

   // add angle 'width'
   // match atom indices w/ GMX indices and then use GMX harmonic.rA value
   // to set angle width

   t_iatom* atoms      = top->idef.il[idefArrayIndex].iatoms;
   t_iparams* params   = top->idef.iparams;
   int offset          = 4;

   // loop over GMX angles
   // set angle width in each of the 3 atoms
   // print warning message if a bond was not found

   for( int ii = 0; ii < top->idef.il[idefArrayIndex].nr; ii += offset ){

      int type          = atoms[ii];
      Real angleWidth   = params[ type ].harmonic.rA*SimTKOpenMMCommon::DegreeToRadians ;
      for( int jj = 0; jj < 3; jj++ ){
         if( jj != 1 ){
            gbsaAngleBond* angleBond = findAngleBond( gbsaBondsArray[atoms[ii + jj + 1]],
                                                      atoms[ii + 1], 
                                                      atoms[ii + 2],
                                                      atoms[ii + 3] );
            if( angleBond != NULL && angleBond->harmonicAngleWidth <= 0.0f ){
               angleBond->harmonicAngleWidth = angleWidth;
            } else if( log != NULL ){
               char buffer[3][MAX_BUFFER_SZ];
               for( int kk = 0; kk < 3; kk++ ){
                  SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atoms[ii + kk + 1], top, MAX_BUFFER_SZ, buffer[kk],
                                                                 maxAtoms, ATOM_ID_STRING_TAB );
               }
               (void) fprintf( log, "%s Warning No angle bond for [%s %s %s] [%d %d %d] in atm=%d", 
                               methodName, buffer[0], buffer[1], buffer[2], atoms[ii + 1],
                               atoms[ii + 2], atoms[ii + 3], atoms[ii + jj + 1] );
            }
         }
      }
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Add stretch bonds to angle data struct (Simbios) 

   @param stretchBondI        stretch bond I
   @param stretchBondJ        stretch bond J
   @param pivotAtomIndex      index of shared atom
   @param gbsaBonds           gbsaBond
   @param saveIndex           index to save????
   @param log                 if set, then print error messages to log file

   @return allocated gbsaAngleBond 

   --------------------------------------------------------------------------------------- */

gbsaAngleBond* GbsaParameters::addStretchBondsToAngleBondList( gbsaStretchBond* stretchBondI, 
                                                               gbsaStretchBond* stretchBondJ,
                                                               int pivotAtomIndex, gbsaBonds** gbsaBonds,
                                                               int saveIndex, FILE* log ){

   // ---------------------------------------------------------------------------------------

   int orderedIndices[3];

   static const char* methodName = "\nGbsaParameters::addStretchBondsToAngleBondList";

   // ---------------------------------------------------------------------------------------

   int duplicateCount =
        sortAtomAngleIndices( stretchBondI->atomI, stretchBondI->atomJ,
                              stretchBondJ->atomI, stretchBondJ->atomJ,
                              orderedIndices );

   if( duplicateCount != 1 && log ){
      (void) fprintf( log, "%s duplicate count from sortAtomAngleIndices is not 1 (=%d) ind=[%d %d %d %d]",
                      methodName, duplicateCount, stretchBondI->atomI, stretchBondI->atomJ,
                      stretchBondJ->atomI, stretchBondJ->atomJ );   
   }

   gbsaAngleBond* newAngleBond;
   if( findAngleBondGivenOrderedArray( gbsaBonds[saveIndex], orderedIndices ) == NULL ){

#ifdef UseGromacsMalloc
      newAngleBond                       = (gbsaAngleBond*) save_malloc( "gbsaAngleBond", __FILE__, __LINE__, sizeof( gbsaAngleBond ) );
#else
      newAngleBond                       = new gbsaAngleBond;
#endif

      memset( newAngleBond, 0, sizeof( gbsaAngleBond ) );

      // gbsaAngleBond* newAngleBond;
      // snew( newAngleBond, 1 );

      newAngleBond->stretchBondI         = stretchBondI;
      newAngleBond->stretchBondJ         = stretchBondJ;
      newAngleBond->pivotAtomIndex       = pivotAtomIndex;
      memcpy( newAngleBond->orderedIndices, orderedIndices, sizeof(int)*3 );
      int volumeIndex = -1;
      for( int ii = 0; ii < 3 && volumeIndex < 0; ii++ ){
         if( orderedIndices[ii] != saveIndex &&  orderedIndices[ii] != pivotAtomIndex ){
            volumeIndex =  orderedIndices[ii];
         }
      }
      newAngleBond->volumeIndex = volumeIndex;

      // push current entry onto list

      newAngleBond->nextBond             = gbsaBonds[saveIndex]->angleBonds;
      gbsaBonds[saveIndex]->angleBonds    = newAngleBond;
      gbsaBonds[saveIndex]->angleBondCount++;

   } else {
      newAngleBond = NULL;
   }

   return newAngleBond;

}

/**---------------------------------------------------------------------------------------

   Find angle bond with atom indices that match input angle indices (Simbios) 

   @param gbsaBond            gbsaBond
   @param atomI               index of atomI
   @param atomJ               index of atomJ
   @param atomK               index of atomK

   @return allocated gbsaAngleBond 

   --------------------------------------------------------------------------------------- */

gbsaAngleBond* GbsaParameters::findAngleBond( const gbsaBonds* gbsaBond, int atomI,
                                              int atomJ, int atomK ) const {

   // ---------------------------------------------------------------------------------------

   int orderedAtomIndices[3];

   // ---------------------------------------------------------------------------------------

   // sort atom indices and them compare to all available bonds
   // return bond when all indices match or NULL if no match 

   sortAtomAngleIndices( atomI, atomI, atomJ, atomK, orderedAtomIndices );

   return findAngleBondGivenOrderedArray( gbsaBond, orderedAtomIndices );

}

/**---------------------------------------------------------------------------------------

   Find angle bond with atom indices that match input ordered angle indices (Simbios) 

   @param gbsaBond            gbsaBond
   @param orderedAtomIndices  array of ordered indices

   @return gbsaAngleBond if found

   --------------------------------------------------------------------------------------- */

gbsaAngleBond* GbsaParameters::findAngleBondGivenOrderedArray( const gbsaBonds* gbsaBond,
                                                               int* orderedAtomIndices ) const {

   // ---------------------------------------------------------------------------------------

   // ---------------------------------------------------------------------------------------

   // sort atom indices and them compare to all available bonds
   // return bond when all indices match or NULL if no match 

   gbsaAngleBond* angleBond = gbsaBond->angleBonds;
   while( angleBond != NULL ){
      if( angleBond->orderedIndices[0] == orderedAtomIndices[0] &&
          angleBond->orderedIndices[1] == orderedAtomIndices[1] &&
          angleBond->orderedIndices[2] == orderedAtomIndices[2] ){
         return angleBond;
      }
      angleBond = angleBond->nextBond;
   }

   return NULL;

}

/**---------------------------------------------------------------------------------------

   Sort atom indices

   @param atomI               index of atomI
   @param atomJ               index of atomJ
   @param atomK               index of atomK
   @param atomL               index of atomL
   @param orderedIndices      output array of ordered indices assumed to be of size 3

   @return count (should be 3, if 4, then all the atom indices were distinct)

   --------------------------------------------------------------------------------------- */

int GbsaParameters::sortAtomAngleIndices( int atomI, int atomJ, 
                                          int atomK, int atomL,
                                          int* orderedIndices ) const { 

   // ---------------------------------------------------------------------------------------

   int orderedAtomIndices[4];

   // ---------------------------------------------------------------------------------------

   // sort atom indices and then remove index with more than 1 entry
   // input should be of the form: (i, j, k, k ), i.e., one duplicate

   orderedAtomIndices[0] = atomI;
   orderedAtomIndices[1] = atomJ;
   orderedAtomIndices[2] = atomK;
   orderedAtomIndices[3] = atomL;

   size_t numberToSort   = 4;
   qsort( orderedAtomIndices, numberToSort, sizeof( int ), &integerComparison );

   orderedIndices[0]  = orderedAtomIndices[0];
   int count          = 1;
   int duplicateCount = 0;
   for( int ii = 1; ii < 4; ii++ ){
      if( orderedAtomIndices[ii] != orderedAtomIndices[ii-1] && count < 3 ){
         orderedIndices[count++] = orderedAtomIndices[ii];
      } else {
         duplicateCount++;
      }
   }

   return duplicateCount;

}

/**---------------------------------------------------------------------------------------

   Get bond counts (Simbios) 

   @param gbsaBond             gbsaBonds to check
   @param stretchCount        stretch count on return
   @param angleCount          angle count on return
 
   @return 0

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getBondCounts( const gbsaBonds* gbsaBond, int* stretchCount, int* angleCount ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\ngetBondCounts";

   // ---------------------------------------------------------------------------------------

   *stretchCount = getStretchBondCount( gbsaBond );
   *angleCount   = getAngleBondCount( gbsaBond );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get stretch bond count (Simbios) 

   @param gbsaBond              gbsaBonds to check

   @return stretch bond count

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getStretchBondCount( const gbsaBonds* gbsaBond ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GbsaParameters::getStretchBondCount";

   // ---------------------------------------------------------------------------------------

   gbsaStretchBond* stretchBonds = gbsaBond->stretchBonds;
   int stretchCount              = 0;
   while( stretchBonds != NULL ){
      stretchCount++;
      stretchBonds = stretchBonds->nextBond;
   }

   return stretchCount;
}

/**---------------------------------------------------------------------------------------

   Get angle bond count (Simbios) 

   @param gbsaBond              gbsaBonds to check

   @return angle bond count

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getAngleBondCount( const gbsaBonds* gbsaBond ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GbsaParameters::getAngleBondCount";

   // ---------------------------------------------------------------------------------------

   gbsaAngleBond* angleBonds = gbsaBond->angleBonds;
   int angleCount             = 0;
   while( angleBonds != NULL ){
      angleCount++;
      angleBonds = angleBonds->nextBond;
   }

   return angleCount;
}

/**---------------------------------------------------------------------------------------

   Add stretch (1-2) bonds (Simbios) 

   @param maxAtoms            max number of atoms
   @param gbsaBondsArray      array of gbsaBonds data structs
   @param top                 Gromacs t_topolgy struct
   @param idefArrayIndex      index to bond parameters (F_BONDS, F_SHAKE, ... )
   @param offset              number of entries for each bond block
   @param atomIndexOffset     offset into block for atom indices
   @param log                 if set, then print error messages to log file

   @return 0 if no errors or
   return x, where x is the number of errors encountered

   --------------------------------------------------------------------------------------- */

int GbsaParameters::addStretchBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                                     const t_topology* top, int idefArrayIndex, 
                                     int offset, int atomIndexOffset, FILE* log  ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GbsaParameters::addStretchBonds";

   // ---------------------------------------------------------------------------------------

   // load 1-2 bonds

   t_iatom* atoms      = top->idef.il[idefArrayIndex].iatoms;
   t_iparams* params   = top->idef.iparams;
   int errors          = 0;

   for( int ii = 0; ii < top->idef.il[idefArrayIndex].nr; ii += offset ){

      int type  = (int) atoms[ ii ];
      int atomI = (int) atoms[ ii + atomIndexOffset ];
      int atomJ = (int) atoms[ ii + atomIndexOffset + 1 ];

      // validate indices

      if( atomI >= maxAtoms || atomI < 0 ){
         if( log ){
            (void) fprintf( log, "\nAtom indexI=%d (ii=%d) too large: max=%d", atomI, ii, maxAtoms );
         }
         errors++;
         atomI = -1;
      }

      if( atomJ >= maxAtoms || atomJ < 0 ){
         if( log ){
            (void) fprintf( log, "\nAtom indexJ=%d (ii=%d) too large: max=%d", atomJ, ii, maxAtoms );
         }
         errors++;
         atomJ = -1;
      }

      // add new stretch bond, if atomJ and/or atomI are not in stretch bond lists

      Real bondLength = params[type].harmonic.rA; 

      if( atomI >= 0 && atomJ >= 0 ){

         if( !isAtomInStretchBondList( atomJ, gbsaBondsArray[atomI], log ) ){
            addStretchBond( atomJ, atomI, gbsaBondsArray[atomI], bondLength, log );
         }

         if( !isAtomInStretchBondList( atomI, gbsaBondsArray[atomJ], log ) ){
            addStretchBond( atomI, atomJ, gbsaBondsArray[atomJ], bondLength, log );
         }
      }

   }

   return errors;
}

/**---------------------------------------------------------------------------------------

   Add SETTLE stretch (1-2) bonds (Simbios) 

   @param maxAtoms            max number of atoms
   @param gbsaBondsArray      array of gbsaBonds data structs
   @param top                 Gromacs t_topolgy struct
   @param idefArrayIndex      index to bond parameters (F_BONDS, F_SHAKE, ... )
   @param offset              number of entries for each bond block
   @param atomIndexOffset     offset into block for atom indices
   @param log                 if set, then print error messages to log file

   @return 0 if no errors or
   return x, where x is the number of errors encountered

   --------------------------------------------------------------------------------------- */

int GbsaParameters::addSettleStretchBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                                           const t_topology* top, int idefArrayIndex, 
                                           int offset, int atomIndexOffset, FILE* log  ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "GbsaParameters::addStretchBonds";

   // ---------------------------------------------------------------------------------------

   // load bonds via SETTLE parameter

   t_iatom* atoms      = top->idef.il[idefArrayIndex].iatoms;
   t_iparams* params   = top->idef.iparams;
   int errors          = 0;

// fprintf( log, "\naddBond idx=%d nr=%d off=%d %d", idefArrayIndex, top->idef.il[idefArrayIndex].nr, offset,
//        atomIndexOffset );
// fflush( log );

   for( int ii = 0; ii < top->idef.il[idefArrayIndex].nr; ii += offset ){

      int type            = (int) atoms[ ii ];
      int atomI           = (int) atoms[ ii + atomIndexOffset ];
      Real OHBondLength   = params[type].harmonic.rA; 
      Real HHBondLength   = params[type].harmonic.krA;
      Real bondAngle      = 1.0 - ( (HHBondLength*HHBondLength)/ (2.0*OHBondLength*OHBondLength));
            bondAngle     = ACOS( bondAngle );

//fprintf( log, "\n%d type=%d [%d %d %d] OH=%.3f HH=%.4f angle=%.3f (deg)",
//         ii, type, atomI, atomI+1, atomI+2, OHBondLength, HHBondLength, bondAngle*180.0f/M_PI );
//fflush( log );

      // validate indices

      if( (atomI + 2) >= maxAtoms || atomI < 0 ){
         if( log ){
            (void) fprintf( log, "\nAtom indexI=%d (ii=%d) too large: max=%d", atomI, ii, maxAtoms );
            (void) fflush( log );
         }
         errors++;
         atomI = -1;
      }

      // add new stretch bond, if atomJ and/or atomI are not in stretch bond lists


      if( atomI >= 0 ){
         gbsaStretchBond* oxygenStretchBonds[2] = { NULL , NULL };
         int hits                               = 0;
         for( int jj = 1; jj <= 2; jj++ ){

            if( !isAtomInStretchBondList( atomI+jj, gbsaBondsArray[atomI], log ) ){
               hits++;
               oxygenStretchBonds[jj-1] = addStretchBond( atomI+jj, atomI, gbsaBondsArray[atomI], OHBondLength, log );
            }

            if( !isAtomInStretchBondList( atomI, gbsaBondsArray[atomI+jj], log ) ){
               addStretchBond( atomI, atomI+jj, gbsaBondsArray[atomI+jj], OHBondLength, log );
            }
         }

         // SETTLE values

         // HOH angle ~ 104.5 deg
         // OH distance=0.09572 
         // HH distance=0.15139

         // HH distance= sqrt( (1-cos(1.82))*2*OH distance*OH distance)

         if( hits == 2 && oxygenStretchBonds[0] && oxygenStretchBonds[1] ){
            for( int jj = 0; jj < 3; jj++ ){
               gbsaAngleBond* gbsaAngleBond = addStretchBondsToAngleBondList( oxygenStretchBonds[0], oxygenStretchBonds[1], 
                                                                              atomI, gbsaBondsArray, atomI + jj, log );
               gbsaAngleBond->harmonicAngleWidth = bondAngle;
            }
         }
      }

   }

   return errors;
}

/**---------------------------------------------------------------------------------------

   Check if atom is in bond list (as atom j) (Simbios)

   @param atomIndex           index of atom to be searched for
   @param bond                gbsaBond data struct
   @param log                 if set, then print error messages to log file

   @return 0 if atom is in StretchBondList; 1 otherwise

   --------------------------------------------------------------------------------------- */

int GbsaParameters::isAtomInStretchBondList( int atomIndex, const gbsaBonds* bond, FILE* log ) const {

   // check if atomIndex is in stretch bond list

   gbsaStretchBond* stretchBond = bond->stretchBonds;
   while( stretchBond != NULL ){ 
      if( stretchBond->atomJ == atomIndex ){
         return 1;
      }
      stretchBond = stretchBond->nextBond;
   }

   return 0;
}

/**---------------------------------------------------------------------------------------

   Add a stretch bonds (Simbios) 

   @param atomIndexI          index of atom I
   @param atomIndexJ          index of atom J
   @param bonds                bond data struct
   @param bondLength          bond length
   @param log                 if set, then print error messages to log file

   @return gbsaStretchBond

   --------------------------------------------------------------------------------------- */

gbsaStretchBond* GbsaParameters::addStretchBond( int atomIndexJ, int atomIndexI,
                                                 gbsaBonds* bonds, Real bondLength, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // ---------------------------------------------------------------------------------------

   // allocate memory, set fields and add to stretchBonds list

#ifdef UseGromacsMalloc
   gbsaStretchBond* newStretchBond       = (gbsaStretchBond*) save_malloc( "gbsaStretchBond", __FILE__, __LINE__, sizeof( gbsaStretchBond ) );
#else
   gbsaStretchBond* newStretchBond       = new gbsaStretchBond;
#endif

   newStretchBond->atomI                = atomIndexI;
   newStretchBond->atomJ                = atomIndexJ;
   newStretchBond->bondLength           = bondLength;

   // push current entry onto list

   newStretchBond->nextBond             = bonds->stretchBonds;
   bonds->stretchBonds                  = newStretchBond;
   bonds->stretchBondCount++;

   return newStretchBond;
}

/**---------------------------------------------------------------------------------------

   Assign standard radii for GB/SA methods other than ACE;
   taken from Macromodel and OPLS-AA, except for hydrogens (Simbios)

   Logic follows that in Tinker's ksolv.f

   Currently only works for standard amino acid atoms
   If invalid atom name is encountered, a message is printed to log file and the
   radius for that atom is set to 1.0f

   @param numberOfAtoms       number of atoms
   @param gbsaBondsArray      array of gbsaBonds
   @param atomNames           array of atom names from GMX top data struct
   @param radii               array to store Macromodel radii for each atom
   @param log                 if set, then print error messages to log file

   @return 0 always

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getMacroModelAtomicRadii( int numberOfAtoms, gbsaBonds** gbsaBondsArray,
                                              char*** atomNames, Real* radii, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGbsaParameters::getMacroModelAtomicRadii";

   // ---------------------------------------------------------------------------------------

   // loop over atoms
   // get Tinker atom number from atom name
   // using atom number and bonds12 array, set atom radius

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      int tinkerAtomNumber = mapGmxAtomNameToTinkerAtomNumber( (*(atomNames[atomI])), log );
      int stretchBondCount = gbsaBondsArray[atomI]->stretchBondCount;

      Real radius;
      int bondedAtom;
      int bondedTinkerAtomNumber;
      switch( tinkerAtomNumber ){

         // H
         case  1:

            // NH?
            bondedAtom             = gbsaBondsArray[atomI]->stretchBonds->atomJ;
            bondedTinkerAtomNumber = mapGmxAtomNameToTinkerAtomNumber( (*(atomNames[bondedAtom])), log );

/*
(void) fprintf( log, "\n    H j=%d Name=<%s>", bondedAtom, (*(atomNames[bondedAtom])) );
(void) fprintf( log, " tinkerNo=%d", bondedTinkerAtomNumber );
(void) fflush( log );
*/
            if( bondedTinkerAtomNumber == 7 ){
               radius = 1.15f;
            } else if( bondedTinkerAtomNumber == 8 ){
               radius = 1.05f;
            } else {
               radius = 1.25f;
            }
            break;

         // ?
         case  3:

            radius = 1.432f;
            break;

         // C
         case  6:

            // C-terminal

            if( stretchBondCount == 3 ){
               radius = 1.875f;
            } else if( stretchBondCount == 2 ){
               radius = 1.825f;
            } else {
               radius = 1.90f;
            }
            break;

         // N
         case 7:

            if( stretchBondCount == 4 ){
               radius = 1.625f;
            } else if( stretchBondCount == 1 ){
               radius = 1.60f;
            } else {
               radius = 1.7063f;
            }
            break;

         // O
         case 8:

            if( stretchBondCount == 1 ){
               radius = 1.48f;
            } else {
               radius = 1.535f;
            }
            break;

         case 9:
            radius = 1.47f;
            break;
         case 10:
            radius = 1.39f;
            break;
         case 11:
            radius = 1.992f;
            break;
         case 12:
            radius = 1.70f;
            break;
         case 14:
            radius = 1.80f;
            break;
         case 15:
            radius = 1.87f;
            break;
         case 16:
            radius = 1.775f;
            break;
         case 17:
            radius = 1.735f;
            break;
         case 18:
            radius = 1.70f;
            break;
         case 19:
            radius = 2.123f;
            break;
         case 20:
            radius = 1.817f;
            break;
         case 35:
            radius = 1.90f;
            break;
         case 36:
            radius = 1.812f;
            break;
         case 37:
            radius = 2.26f;
            break;
         case 53:
            radius = 2.10f;
            break;
         case 54:
            radius = 1.967f;
            break;
         case 55:
            radius = 2.507f;
            break;
         case 56:
            radius = 2.188f;
            break;
         default:
            radius = 1.0f;
            (void) fprintf( log ? log : stderr, "\nWarning: %s Tinker atom number=%d unrecognized -- name=<%s>.",
                            methodName, tinkerAtomNumber, (*(atomNames[bondedAtom])) );
           break;
      }

      // convert from A to nm

      radii[atomI] = 0.1f*radius;
   }
      
   return 0;
}

/**---------------------------------------------------------------------------------------

   Map Gmx atom name to Tinker atom number (Simbios)

   @param atomName            atom name (CA, HA, ...); upper and lower case should both work
   @param log                 if set, then print error messages to log file

   @return Tinker atom number if atom name is valid; else return -1

   --------------------------------------------------------------------------------------- */

int GbsaParameters::mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const {

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

   Get atomic volumes minus subvolumes that lie inside directly bonded atoms
   Eqs 6 & 7 in J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   Use Macromodel radii

   @param numberOfAtoms       number of atoms
   @param top                 GMX t_topology struct
   @param gbsaBondsArray      array of gbsaBonds
   @param vdwRadius           array of Macromodel radii for each atom
   @param volume              output array of volumes
   @param log                 if set, then print error messages to log file

   @return 0 always; atomic volumes in array 'volume'

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getMacroModelAtomicVolumes( int numberOfAtoms, const t_topology* top,
                                                gbsaBonds** gbsaBondsArray, Real* vdwRadius,
                                                Real* volume, FILE* log ){

   // ---------------------------------------------------------------------------------------

   char buffer[MAX_BUFFER_SZ];

   bool debugOn                  = false;
   // static const char* methodName = "\nGbsaParameters::getMacroModelAtomicVolumes";

   // ---------------------------------------------------------------------------------------

   debugOn = debugOn && log != NULL;

   Real pi_3 = ((Real) M_PI)/3.0f;
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // include pi/3 after finished summing over h_ij

      Real volume_I = 4.0*POW( vdwRadius[atomI], 3.0 );

      Real atomI_R  = vdwRadius[atomI];
      Real atomI_R2 = atomI_R*atomI_R;

      // loop over directly bonded neighbors

      gbsaStretchBond* stretchBond = gbsaBondsArray[atomI]->stretchBonds;
      while( stretchBond != NULL ){ 

         int atomJ           =  stretchBond->atomJ;
         Real bondLength     =  stretchBond->bondLength*1.01f;
         Real atomJ_R        = vdwRadius[atomJ];
         Real ratio          = (atomJ_R*atomJ_R - atomI_R2 - bondLength*bondLength)/(2.0f*atomI_R*bondLength);
         Real h_IJ           = atomI_R*( 1.0f + ratio );
         volume_I           -= h_IJ*h_IJ*( 3.0f*atomI_R - h_IJ ); 

         if( debugOn ){
            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atomJ, top, MAX_BUFFER_SZ, buffer, numberOfAtoms, ATOM_ID_STRING_TAB );
            (void) fprintf( log, "\n   atmJ=%d %s bndL=%.3f jR=%.3f h=%.3f offset=%.3f",
                            atomJ, buffer, bondLength, atomJ_R, h_IJ,  (h_IJ*h_IJ*( 3.0f*atomI_R - h_IJ )) );  
         }

         stretchBond = stretchBond->nextBond;
      }
      volume[atomI]     = volume_I*pi_3;

      if( debugOn ){
         SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atomI, top, MAX_BUFFER_SZ, buffer, numberOfAtoms, ATOM_ID_STRING_TAB );
         (void) fprintf( log, "\n%s %.3f %.3f %d  Bonded: ", buffer,
                         10.0f*vdwRadius[atomI], 1000.0f*volume[atomI], gbsaBondsArray[atomI]->stretchBondCount );
         stretchBond = gbsaBondsArray[atomI]->stretchBonds;
         while( stretchBond != NULL ){ 
            (void) fprintf( log, " %d", stretchBond->atomJ );
            stretchBond = stretchBond->nextBond;
         }
      }
   }

   return 0;

}

/**---------------------------------------------------------------------------------------

   Calculate first two terms in Gpol equation (nearest and next-nearest neighbors) 
   Eqs 6 & 7 in J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   Use Macromodel radii

   @param numberOfAtoms       number of atoms
   @param top                 GMX t_topology struct
   @param gbsaBondsArray      array of gbsa bonds
   @param vdwRadius           array of Macromodel radii for each atom
   @param volume              array of atomic volumes
   @param gPolFixed           output array of fixed GBSA values for GPol
   @param log                 if set, then print error messages to log file

   @return 0 always; fixed value for G_pol in array 'gPolFixed'

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getFixedGBSA_GPol( int numberOfAtoms, const  t_topology* top,
                                       gbsaBonds** gbsaBondsArray, Real* vdwRadius,
                                       Real* volume, Real* gPolFixed, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGbsaParameters::getFixedGBSA_GPol";

   // ---------------------------------------------------------------------------------------

   // GBSA parameters from '97 paper

   Real phi               =  getPhi();
   Real P1                =  getP1();
   Real P2                =  getP2();
   Real P3                =  getP3();
   Real electricConstant  = getElectricConstant();

   // ---------------------------------------------------------------------------------------

   // 1st term
 
   // convert to nm

   Real P1_Phi = 0.1f*(P1 + phi);
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      gPolFixed[atomI] = electricConstant/(vdwRadius[atomI]+P1_Phi);
   }

   // stretch term (1-2 bonds) 

   char idBuffer[MAX_BUFFER_SZ];
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // loop over directly bonded neighbors

      Real volumeR4Sum = 0.0f;
      gbsaStretchBond* stretchBond = gbsaBondsArray[atomI]->stretchBonds;
      while( stretchBond != NULL ){ 
         int atomJ           = stretchBond->atomJ;
         Real volumeJ_R      = volume[atomJ];
         Real bondLength     = stretchBond->bondLength;
         Real r4             = POW( bondLength, 4.0 ); 
         if( r4 > 0.0 ){
            volumeR4Sum     += volumeJ_R/r4;
         } else if( log != NULL ){
            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atomI, top, MAX_BUFFER_SZ, idBuffer, numberOfAtoms, ATOM_ID_STRING_TAB );
            (void) fprintf( log, "%s %s 1-2 Bonded atom atomIndex=%d has bond length=0.",
                            methodName, idBuffer, atomJ ); 
         }
         stretchBond = stretchBond->nextBond;
         
      }
      gPolFixed[atomI]    += P2*volumeR4Sum;
   }

   // bend term (1-3 bonds) 

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // loop over angles

      Real volumeR4Sum = 0.0f;
      gbsaAngleBond* angleBond = gbsaBondsArray[atomI]->angleBonds;
      while( angleBond != NULL ){ 

         gbsaStretchBond* gbsaStretchBond1 = angleBond->stretchBondI;
         gbsaStretchBond* gbsaStretchBond2 = angleBond->stretchBondJ;

         Real bondLength1                 = gbsaStretchBond1->bondLength;
         Real bondLength2                 = gbsaStretchBond2->bondLength;

         Real cosOfAngle                  = cos( angleBond->harmonicAngleWidth );

         Real distance2                   = bondLength1*bondLength1 +
                                             bondLength2*bondLength2 -
                                             2.0f*bondLength1*bondLength2*cosOfAngle;

         Real r4                          = distance2*distance2;

         if( r4 > 0.0 ){
            volumeR4Sum += volume[angleBond->volumeIndex]/r4;
         } else if( log != NULL ){
            SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( atomI, top, MAX_BUFFER_SZ, idBuffer, numberOfAtoms, ATOM_ID_STRING_TAB );
            (void) fprintf( log, "%s %s angle-bonded atom has distance 0.",
                            methodName, idBuffer ); 
         }
         angleBond = angleBond->nextBond;
      }

      gPolFixed[atomI]    += P3*volumeR4Sum;
   }

   return 0;

}

/**---------------------------------------------------------------------------------------

   Get exclusions for specific atom (Simbios)

   @param numberOfAtoms       number of atoms
   @param gbsaBondsArray      array of gbsa bonds
   @param atomI               atom index of atom for which exclusions are to be set
   @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
   @param                      value may be null on input in which space is allocated
   @param previousIndex       previousIndex -- if < 0, then iniitialize all entries to
   @param                      0
   @param log                 if set, then print error messages to log file

   @return 0

   Abort if exclusionWorkArray == NULL

   --------------------------------------------------------------------------------------- */

int GbsaParameters::getExclusionsForAtom( int numberOfAtoms, const gbsaBonds** gbsaBondsArray,
                                          int atomI, int* exclusionWorkArray, int previousIndex,
                                          FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nGbsaParameters::getExclusionsForAtom";

   // ---------------------------------------------------------------------------------------

   // allocate memory if exclusionWorkArray is null
 
   if( exclusionWorkArray == NULL ){
      (void) fprintf( log ? log : stderr, "%s error exclusionWorkArray is NULL.", methodName );
      (void) fflush( log ? log : stderr );
      exit(-1);
      return -1;
   }

   // set exclusion factors

   if( previousIndex < 0 ){
      if( previousIndex == -1 ){
         for( int ii = 0; ii < numberOfAtoms; ii++ ){
            exclusionWorkArray[ii] = -1;
         }
      } else {
         memset( exclusionWorkArray, 0, numberOfAtoms*sizeof( int ) );
      }
   } else {
      setExclusionValue( gbsaBondsArray[previousIndex], -1, exclusionWorkArray, log );
   }

   setExclusionValue( gbsaBondsArray[atomI], atomI, exclusionWorkArray, log );

   return 0;

}

/**---------------------------------------------------------------------------------------

   Set exclusions for specific atom (Simbios)

   @param gbsaBonds           gbsa bond
   @param setValue            set value
   @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
   @param                      value may be null on input in which space is allocated
   @param log                 if set, then print error messages to log file

   @return array of Born radii

   --------------------------------------------------------------------------------------- */

int GbsaParameters::setExclusionValue( const gbsaBonds* gbsaBonds, int setValue,
                                       int* exclusionWorkArray, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::setExclusionsValue";

   // ---------------------------------------------------------------------------------------

   gbsaStretchBond* stretchBond = gbsaBonds->stretchBonds;
   while( stretchBond != NULL ){
      exclusionWorkArray[stretchBond->atomI] = setValue;
      exclusionWorkArray[stretchBond->atomJ] = setValue;
      stretchBond = stretchBond->nextBond;
   }

   gbsaAngleBond* angleBond = gbsaBonds->angleBonds;
   while( angleBond != NULL ){ 
      exclusionWorkArray[angleBond->stretchBondI->atomI] = setValue;
      exclusionWorkArray[angleBond->stretchBondI->atomJ] = setValue;
      exclusionWorkArray[angleBond->stretchBondJ->atomI] = setValue;
      exclusionWorkArray[angleBond->stretchBondJ->atomJ] = setValue;
      angleBond = angleBond->nextBond;
   }

   return 0;

}

/**---------------------------------------------------------------------------------------

   Get array for accessing exclusion entries (Simbios)
   (work array)

   @return array of size _numberOfAtoms (int)

   --------------------------------------------------------------------------------------- */

int* GbsaParameters::getExclusionWorkArray( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::getExclusionWorkArray";

   // ---------------------------------------------------------------------------------------

   if( _exclusionWorkArray == NULL && _numberOfAtoms > 0 ){
#ifdef UseGromacsMalloc
      _exclusionWorkArray = (int*) save_malloc( "_exclusionWorkArray", __FILE__, __LINE__, sizeof( int )*_numberOfAtoms );
#else
      _exclusionWorkArray = new int[_numberOfAtoms];
#endif

      memset( _exclusionWorkArray, 0, _numberOfAtoms*sizeof( int ) );
   }
   return _exclusionWorkArray;

}

/**---------------------------------------------------------------------------------------
      
   Print state to log file (Simbios)
   
   @param title               title (optional)
   @param log                 print state to log file
      
   @return 0 always
      
   --------------------------------------------------------------------------------------- */

int GbsaParameters::logState( const char* title, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nGbsaParameters::logState";

   // ---------------------------------------------------------------------------------------

   if( log == NULL ){
      return 0;
   }
   if( title ){
      (void) fprintf( log, "%s", title );
   }

   (void) fprintf( log, "\nGbsaParameter info:\n" );
   (void) fprintf( log, "\n   Solvent dielectric:      %12.4f", _solventDielectric );
   (void) fprintf( log, "\n     Inner dielectric:      %12.4f", _innerDielectric ); 
   (void) fprintf( log, "\n    Electric constant:      %12.4f", _electricConstant ); 
   (void) fprintf( log, "\n       Gpol prefactor:      %12.4f", _preFactor ); 
   (void) fprintf( log, "\n         Probe radius:      %12.4f", _probeRadius );
   (void) fprintf( log, "\n                  Phi:      %12.4f", _phi );
   (void) fprintf( log, "\n                   P1:      %12.4f", _P1 );
   (void) fprintf( log, "\n                   P2:      %12.4f", _P2 );
   (void) fprintf( log, "\n                   P3:      %12.4f", _P3 );
   (void) fprintf( log, "\n                   P4:      %12.4f", _P4 );
   (void) fprintf( log, "\n                   P5:      %12.4f", _P5 );
   (void) fprintf( log, "\n      Number of atoms:      %12d",   _numberOfAtoms );
   (void) fprintf( log, "\n          Free arrays:      %12d",   _freeArrays );

   return 0;


}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

void* GbsaParameters::operator new( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

/*
   (void) fprintf( stdout, "\nGbsaParameters new called -- size=%u", size );
   (void) fflush( stdout );
*/

   return ptr;
}

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

void GbsaParameters::operator delete( void *ptr ){

/*
   (void) fprintf( stdout, "\nGbsaParameters delete called." );
   (void) fflush( stdout );
*/

   sfree( ptr ); 
}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

void* GbsaParameters::operator new[]( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

/*
   (void) fprintf( stdout, "\nGbsaParameters new[] called -- size=%u", size );
   (void) fflush( stdout );
*/

   return ptr;
}

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

void GbsaParameters::operator delete[]( void *ptr ){

/*
   (void) fprintf( stdout, "\nGbsaParameters delete[] called." );
   (void) fflush( stdout );
*/

   sfree( ptr ); 
}

/**---------------------------------------------------------------------------------------

   Get Vdw radii from parameter file (Simbios) 
      
   @param numberOfAtoms       number of atoms
   @param parameterFileName   parameterFileName
   @param top                 Gromacs topology data struct
   @param atomNames           array of atom names from GMX top data struct
   @param radii               array to store Macromodel radii for each atom
   @param log                 if set, then print error messages to log file
      
   @return 0 always
      
   --------------------------------------------------------------------------------------- */

int GbsaParameters::getMacroModelAtomicRadii( int numberOfAtoms, const std::string parameterFileName,
                                              const t_topology* top, char*** atomNames,
                                              Real* radii, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nGbsaParameters::getMacroModelAtomicRadii";

   // ---------------------------------------------------------------------------------------

   StringVector* fileContents         = readFile( parameterFileName ); 
   GbsaAtomParameterMap* parameterMap = parseParameterFile( *fileContents );
   delete fileContents;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      char* atomTypeC = *(top->atoms.atomtype[ii]);   
      std::string atomType = atomTypeC;
      GbsaAtomParameterMapCI entry = parameterMap->find( atomType ); 
      if( entry != parameterMap->end() ){
         GbsaAtomParameter* gbsaAtomParameter = (*entry).second;
         radii[ii]  = (Real) gbsaAtomParameter->getVdwRadius();
         radii[ii] *= 0.1;

         std::stringstream message;
         message << " type found for atom=<" << (*(atomNames[ii])) << "> type=<" << atomType << "> r=" << radii[ii];
         (void) fprintf( log, "%s\n" , message.str().c_str() );

      } else {
         radii[ii] = (Real) 0.0001f;
         std::stringstream message;
         message << methodName;
         message << " no type found for atom=<" << (*(atomNames[ii])) << "> type=<" << atomType << ">";
         (void) fprintf( log, "%s\n" , message.str().c_str() );
      }
   }
   
   delete parameterMap;

   return 0;
}

/**---------------------------------------------------------------------------------------

   Read file (Simbios) 

   @param fileName       file name

   @return vector of strings, one string per line; return NULL if file is not found

   --------------------------------------------------------------------------------------- */

StringVector* GbsaParameters::readFile( const std::string& fileName ){

   // ---------------------------------------------------------------------------------------

   static int bufferSize = 2048;
   static char buffer[2048];
   static const std::string methodName = "\nGbsaParameters::readParameterFile";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* file = NULL;
#ifdef WIN32
   fopen_s( &file, fileName.c_str(), "r" );
#else
   file = fopen( fileName.c_str(), "r" );
#endif

   if( file != NULL ){
      std::stringstream message;
      message << methodName.c_str() << " opened file=<" <<  fileName.c_str() << ">.";
      // AmoebaLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName.c_str() << " could not open file=<" <<  fileName.c_str() << ">.";
      (void) printf( "\n%s\n", message.str().c_str() ); 
      (void) fflush( stdout );
      // AmoebaLog::printMessage( message );
      exit(-1);
      return NULL;
   }

   // loop over all lines in file, checking for end-of-file

   int lineNumber             = 0;
   StringVector* fileContents = new StringVector(); 

   while( !feof( file ) ){ 

      // read next line

      int bufferLen;
      lineNumber++;
      if( fgets( buffer, bufferSize, file ) != NULL ){
         bufferLen = (int) strlen( buffer );
         if( bufferLen > 0 ){
            buffer[bufferLen-1] = '\0';
         }
         // (void) fprintf( log, "%s", buffer );
         // (void) fflush( log );
         fileContents->push_back( buffer );
      }
   }

   // done

   (void) fclose( file );

   if( file ){
      //std::stringstream message;
      //message << methodName.c_str() << " read " << lineNumber << " lines from file=<" << fileName->c_str() << ">.";
      // AmoebaLog::printMessage( message );
   }

   return fileContents;
}

/**---------------------------------------------------------------------------------------

   Parse parameter file (Simbios) 

   @param fileContents  file contents

   @return vector of strings, one string per line

   --------------------------------------------------------------------------------------- */

GbsaAtomParameterMap* GbsaParameters::parseParameterFile( const StringVector& fileContents ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nGbsaParameters::parseParameterFile";
   static const std::string delimiter  = " ";

   // ---------------------------------------------------------------------------------------

   GbsaAtomParameterMap* map = new GbsaAtomParameterMap();
   for( StringVectorCI ii = fileContents.begin(); ii != fileContents.end(); ii++ ){

      std::string line = *ii;
      StringVector tokenVector;
      tokenizeString( line, tokenVector, delimiter, 0 );
      if( tokenVector[0] != "#" && tokenVector[0] != "@" ){
         GbsaAtomParameter* gbsaAtomParameter = new GbsaAtomParameter( tokenVector );      
         (*map)[tokenVector[0]] = gbsaAtomParameter;
      }
   }

   return map;
}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param line                 string to tokenize
   @param tokenVector          upon return vector of tokens
   @param delimiter            token delimter
   @param clearTokenVector     if true, clear tokenVector

   @return 1

   --------------------------------------------------------------------------------------- */

int GbsaParameters::tokenizeString( const std::string& line, StringVector& tokenVector,
                                    const std::string& delimiter, int clearTokenVector ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nGbsaParameters::tokenizeString";

   // ---------------------------------------------------------------------------------------

   char *ptr_c;
   if( clearTokenVector ){
      tokenVector.clear();
   }   
   char* lineBuffer = (char*) malloc( 2048*sizeof( char) );
#ifdef WIN32
   (void) sprintf_s( lineBuffer, 2048, "%s", line.c_str() );
#else
   (void) sprintf( lineBuffer, "%s", line.c_str() );
#endif

   while( (ptr_c = strsep( &lineBuffer, delimiter.c_str() )) != NULL ){
      if( *ptr_c ){
         tokenVector.push_back( std::string( ptr_c ) );
      }   
   }
   free( lineBuffer );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse parameter file lines

   Should be moved to Utilities file

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args; if return value equals maxTokens, then more tokens than allocated

   --------------------------------------------------------------------------------------- */

char* GbsaParameters::strsep( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nTinkerParameterSet::strsep"

   char *s;
   const char *spanp;
   int c, sc;
   char *tok;

   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }

   for( tok = s;; ){
      c = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}

/**---------------------------------------------------------------------------------------

   Qsort/heapsort integer comparison (Simbios) 

   @param a                    first value to compare
   @param b                    second value to compare

   @return 1, 0, -1 based on comparison

   --------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b){

   int diff = *(int *) a - *(int *) b;

   if( diff < 0 ){ 
     return -1;
   } else if( diff > 0 ){ 
     return 1;
   } else {
     return 0;
   }
}
