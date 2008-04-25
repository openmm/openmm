
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

#include <string.h>
#include <math.h>

#include "UtilitiesSimTk.h"
#include "CpuGbsa.h"

#define UseGromacsMalloc 1

// Replacement new/delete w/ Gromac's smalloc() and sfree()

#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

// static data member created by call to cpuSetGbsaParameters() 
// stores parameter settings, ...
// used to calculate GBSA forces/energy

CpuGbsa* CpuGbsa::_cpuGbsa = NULL;

/**---------------------------------------------------------------------------------------

   CpuGbsa constructor

   @param numberOfAtoms       number of atoms
   @param top                 Gromac's topology struct
   @param log                 if set, then print error messages to log file
   
   --------------------------------------------------------------------------------------- */

CpuGbsa::CpuGbsa( int numberOfAtoms, const t_topology* top, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::CpuGbsa";

   // ---------------------------------------------------------------------------------------

   initializeDataMembers( );

   _gbsaParameters = new GbsaParameters( numberOfAtoms, log );
   _gbsaParameters->initializeParameters( top );
}

/**---------------------------------------------------------------------------------------

   CpuGbsa constructor

   gbsaParameters      gbsaParameters object
   
   --------------------------------------------------------------------------------------- */

CpuGbsa::CpuGbsa( GbsaParameters* gbsaParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::CpuGbsa(2)";

   // ---------------------------------------------------------------------------------------

   initializeDataMembers( );
   _gbsaParameters = gbsaParameters;

}

/**---------------------------------------------------------------------------------------

   CpuGbsa destructor

   --------------------------------------------------------------------------------------- */

CpuGbsa::~CpuGbsa( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::~CpuGbsa";

   // ---------------------------------------------------------------------------------------

   if( _gbsaParameters != NULL ){
      delete _gbsaParameters;
   }

   if( _gbsaBornForce != NULL ){
#ifdef UseGromacsMalloc
      save_free( "_gbsaBornForce", __FILE__, __LINE__, _gbsaBornForce );
#else
      delete[] _gbsaBornForce;
#endif
   }

   if( _gbsaBornRadiiTemp != NULL ){
#ifdef UseGromacsMalloc
      save_free( "_gbsaBornRadiiTemp", __FILE__, __LINE__, _gbsaBornRadiiTemp );
#else
      delete[] _gbsaBornRadiiTemp;
#endif
   }

}

/**---------------------------------------------------------------------------------------

   Delete static _cpuGbsa object if set

   @return 0 if _cpuGbsa was set; otherwise return -1

   --------------------------------------------------------------------------------------- */

int CpuGbsa::deleteCpuGbsa( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::deleteCpuGbsa";

   // ---------------------------------------------------------------------------------------

   if( _cpuGbsa != NULL ){
      delete _cpuGbsa;
      UtilitiesSimTk::getAtomIdStringGivenAtomIndex( -1, NULL, NULL, 0, 0 );
      _cpuGbsa = NULL;
      return 0;
   } else {
      return -1;
   }
}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuGbsa::initializeDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _gbsaBornForce            = NULL;
   _gbsaBornRadiiTemp        = NULL;
 
   _includeAceApproximation  = false;

   _forceConversionFactor    = 1.0f;
   _energyConversionFactor   = 1.0f;

   _gbsaEnergy               = 0.0f;
}

/**---------------------------------------------------------------------------------------

   Return _gbsaBornForce, a work array
   On first call, memory for array is allocated

   --------------------------------------------------------------------------------------- */

float* CpuGbsa::getGbsaBornForce( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::getGbsaBornForce";

   // ---------------------------------------------------------------------------------------

   if( _gbsaBornForce == NULL ){
#ifdef UseGromacsMalloc
      // smalloc( (void*) _gbsaBornForce, sizeof( float )*_gbsaParameters->getNumberOfAtoms() );
      _gbsaBornForce = (float*) save_malloc( "_gbsaBornForce", __FILE__, __LINE__, sizeof( float )*_gbsaParameters->getNumberOfAtoms() );
#else
      _gbsaBornForce = new float[_gbsaParameters->getNumberOfAtoms()];
#endif
   }
   return _gbsaBornForce;

}

/**---------------------------------------------------------------------------------------

   Return _gbsaBornRadiiTemp, a work array
   On first call, memory for array is allocated

   --------------------------------------------------------------------------------------- */

float* CpuGbsa::getGbsaBornRadiiTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::getGbsaBornRadiiTemp";

   // ---------------------------------------------------------------------------------------

   if( _gbsaBornRadiiTemp == NULL ){
#ifdef UseGromacsMalloc
      _gbsaBornRadiiTemp = (float*) save_malloc( "_gbsaBornRadiiTemp", __FILE__, __LINE__, sizeof( float )*_gbsaParameters->getNumberOfAtoms() );
#else
      _gbsaBornRadiiTemp = new float[_gbsaParameters->getNumberOfAtoms()];
#endif
   }
   return _gbsaBornRadiiTemp;
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii

   @return array of Born radii

   --------------------------------------------------------------------------------------- */

int CpuGbsa::computeBornRadii( const rvec *atomCoordinates, float* bornRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   GbsaParameters* gbsaParameters = getGbsaParameters();

   int numberOfAtoms              = gbsaParameters->getNumberOfAtoms();
   FILE* log                      = gbsaParameters->getLog();
   float* vdwRadius               = gbsaParameters->getVdwRadii();
   float* volume                  = gbsaParameters->getVolume();
   float* gPolFixed               = gbsaParameters->getGPolFixed();
   gbsaBonds** gbsaBondsArray     = gbsaParameters->getGbsaBondsArray();

   // ---------------------------------------------------------------------------------------

   float P4                      = gbsaParameters->getP4();
   float P4_2                    = gbsaParameters->getP4_2();
   float invertGpol              = gbsaParameters->getElectricConstant();
   float P5Inverse               = gbsaParameters->getP5Inverse();
   float piP5                    = gbsaParameters->getPiP5();

   // ---------------------------------------------------------------------------------------

   // exclusion work array
 
   int* exclusionWorkArray = gbsaParameters->getExclusionWorkArray();

   // intialize Born radii w/ fixed portion of Gpol

   memcpy( bornRadii, gPolFixed, sizeof( float )*numberOfAtoms );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
     
      // set exclusion factors

      getExclusionsForAtom( numberOfAtoms, gbsaBondsArray,
                            atomI, exclusionWorkArray, atomI - 1, log );

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( exclusionWorkArray[atomJ] != atomI ){

            float deltaX = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            float deltaY = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            float deltaZ = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
            float r2     = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            float r4     = r2*r2;
            float vdWSum = vdwRadius[atomI] +  vdwRadius[atomJ];
            float ratio  = r2/(vdWSum*vdWSum);

            float ccf;
            if( ratio > P5Inverse ){
               ccf         = P4;
            } else {
               float term  = P4_2*(1.0f-cos(ratio*piP5));
               ccf         = term*term;
            }

            bornRadii[atomI] +=  (ccf*volume[atomJ])/r4;

         }
      }
      if( bornRadii[atomI] != 0.0 ){
         bornRadii[atomI] = invertGpol/bornRadii[atomI];
      } else {
         bornRadii[atomI] = 0.1f;
      }
   }

   return 0;

}

/**---------------------------------------------------------------------------------------

   Get exclusions for specific atom (Simbios)

   @param numberOfAtoms       number of atoms
   @param gbsaBondsArray      array of gpu bonds
   @param atomI               atom index of atom for which exclusions are to be set
   @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
                              value may be null on input in which space is allocated
   @param previousIndex       previousIndex -- if < 0, then iniitialize all entries to 0
   @param log                 if set, then print error messages to log file

   @return 0

   --------------------------------------------------------------------------------------- */

int CpuGbsa::getExclusionsForAtom( int numberOfAtoms, gbsaBonds** gbsaBondsArray,
                                   int atomI, int* exclusionWorkArray, int previousIndex,
                                   FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuGbsa::getExclusionsForAtom";

   // ---------------------------------------------------------------------------------------

   // abort if exclusionWorkArray is null
 
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

   @param gbsaBonds           gpu bond
   @param setValue            set value
   @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded \n
                              value may be null on input in which space is allocated \n
   @param log                 if set, then print error messages to log file

   @return 0

   --------------------------------------------------------------------------------------- */

int CpuGbsa::setExclusionValue( gbsaBonds* gbsaBonds, int setValue, int* exclusionWorkArray, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::setExclusionsValue";

   // ---------------------------------------------------------------------------------------

   gbsaStretchBond* stretchBond = gbsaBonds->stretchBonds;
   while( stretchBond != NULL ){
      exclusionWorkArray[stretchBond->atomI] = setValue;
      exclusionWorkArray[stretchBond->atomJ] = setValue;
      stretchBond                            = stretchBond->nextBond;
   }

   gbsaAngleBond* angleBond = gbsaBonds->angleBonds;
   while( angleBond != NULL ){ 
      exclusionWorkArray[angleBond->stretchBondI->atomI] = setValue;
      exclusionWorkArray[angleBond->stretchBondI->atomJ] = setValue;
      exclusionWorkArray[angleBond->stretchBondJ->atomI] = setValue;
      exclusionWorkArray[angleBond->stretchBondJ->atomJ] = setValue;
      angleBond                                          = angleBond->nextBond;
   }

   return 0;

}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
   Use static member of class, _cpuGbsa, which should have been initialized
   prior to the call here

   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces (output)
   @param log                 if set, then print error messages to log file

   @return force array

   --------------------------------------------------------------------------------------- */

float** CpuGbsa::computeGbsaForces( const rvec *atomCoordinates, const float* partialCharges,
                                    float** forces, FILE* log ){
   // ---------------------------------------------------------------------------------------

   int printSampleOutput         = 0;
   static int callId             = 0;
   static const char* methodName = "\nCpuGbsa::computeGbsaForces";

   // ---------------------------------------------------------------------------------------

   // check intialized correctly

   if( _cpuGbsa == NULL ){
      (void) fprintf( log ? log : stderr, "%s CpuGbsa has not been initialized!", methodName );
      return forces;
   }
   callId++;

   GbsaParameters* gbsaParameters = _cpuGbsa->getGbsaParameters();

   // check to see if Born radii have been previously calculated
   // if not, then calculate;
   // logic here assumes that the radii are intitialized to zero 
   // and then once computed, always greater than zero.

   float* bornRadii = gbsaParameters->getBornRadii();
   if( bornRadii[0] < 0.0001f ){

      if( log && gbsaParameters->getLog() == NULL ){
         gbsaParameters->setLog( log );
      }

      _cpuGbsa->computeBornRadii( atomCoordinates, bornRadii );

      // diagnostics

      if( printSampleOutput && log ){
         int numberOfAtoms = gbsaParameters->getNumberOfAtoms();
         (void) fprintf( log, "%s initialize Born radii for %d atoms on call=%d", 
                         methodName, numberOfAtoms, callId  );
         for( int ii = 0; ii < printSampleOutput && ii < numberOfAtoms; ii++ ){
            (void) fprintf( log, "\n%5d %12.4f", ii, bornRadii[ii] );
         }
         int startIndex = gbsaParameters->getNumberOfAtoms() - printSampleOutput > 0 ?
                          gbsaParameters->getNumberOfAtoms() - printSampleOutput : 0;
         for( int ii = startIndex; ii < numberOfAtoms; ii++ ){
            (void) fprintf( log, "\n%5d %12.4f", ii, bornRadii[ii] );
         }
      }
   }

   forces = _cpuGbsa->computeBornEnergyForces( gbsaParameters->getBornRadii(), atomCoordinates,
                                               partialCharges, forces );

   // diagnostics

   if( printSampleOutput && log ){

      int numberOfAtoms = gbsaParameters->getNumberOfAtoms();
      (void) fprintf( log, "%s call=%d E=%.4f forces %d atoms", methodName, callId, _cpuGbsa->_gbsaEnergy, numberOfAtoms  );
      for( int ii = 0; ii < printSampleOutput && ii < numberOfAtoms; ii++ ){
         (void) fprintf( log, "\n%5d [%12.4f %12.4f %12.4f] %12.4f", ii, forces[ii][0], forces[ii][1], forces[ii][2], bornRadii[ii] );
      }
      int startIndex = gbsaParameters->getNumberOfAtoms() - printSampleOutput > 0 ?
                       gbsaParameters->getNumberOfAtoms() - printSampleOutput : 0;
      for( int ii = startIndex; ii < numberOfAtoms; ii++ ){
         (void) fprintf( log, "\n%5d [%12.4f %12.4f %12.4f] %12.4f", ii, forces[ii][0], forces[ii][1], forces[ii][2], bornRadii[ii] );
      }
   }

   return forces;
}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   @param bornRadii           Born radii -- optional; if NULL, then GbsaParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces
   @param log                 if set, then print error messages to log file

   @return force array

   The array bornRadii is also updated and the gbsaEnergy

   --------------------------------------------------------------------------------------- */

float** CpuGbsa::computeBornEnergyForces( float* bornRadii, const rvec *atomCoordinates,
                                          const float* partialCharges, float** forces ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::computeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   GbsaParameters* gbsaParameters = getGbsaParameters();
   int numberOfAtoms              = gbsaParameters->getNumberOfAtoms();
   FILE* log                      = gbsaParameters->getLog();
   float* vdwRadius               = gbsaParameters->getVdwRadii();
   float* volume                  = gbsaParameters->getVolume();
   float* gPolFixed               = gbsaParameters->getGPolFixed();
   if( bornRadii == NULL ){
      bornRadii                   = gbsaParameters->getBornRadii();
   }
   gbsaBonds** gbsaBondsArray     = gbsaParameters->getGbsaBondsArray( );

   // ---------------------------------------------------------------------------------------

   // constants

   float preFactor               = gbsaParameters->getPreFactor();
   float P5Inverse               = gbsaParameters->getP5Inverse();
   float P4                      = gbsaParameters->getP4();
   float piP5                    = gbsaParameters->getPiP5();
   float electricConstant        = gbsaParameters->getElectricConstant();

   // ---------------------------------------------------------------------------------------

   // get exclusion work array
 
   int* exclusionWorkArray       = gbsaParameters->getExclusionWorkArray();

   // set energy to zero

   _gbsaEnergy                   = 0.0f;

   // ---------------------------------------------------------------------------------------

   // allocate force array if not allocated and initialize values to 0

   if( !forces ){
      int numberForceBlocks         = 3;
      // forces                        = UtilitiesSimTk::allocate2DFloatArray( numberForceBlocks, numberOfAtoms, forces, false, 0.0f, log );
      forces                        = UtilitiesSimTk::allocate2DFloatArray( numberOfAtoms, numberForceBlocks, forces, false, 0.0f, log );
   }
   unsigned int arraySzInBytes      = sizeof( float )*numberOfAtoms;

/*
   for( int ii = 0; ii < 3; ii++ ){
      memset( forces[ii], 0, arraySzInBytes );
   }
*/

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      memset( forces[ii], 0, 3*sizeof( float ));
   }

   float* bornForces                = getGbsaBornForce();

   // N*( 8 + pow) ACE
    
   if( includeAceApproximation() ){

      // compute the nonpolar solvation via ACE approximation

      float probeRadius       = gbsaParameters->getProbeRadius();
      float surfaceAreaFactor = gbsaParameters->getPi4Asolv();

      // 1 + 1 + pow + 3 + 1 + 2 FLOP

      for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
         if( bornRadii[atomI] > 0.0f ){
            float r            = vdwRadius[atomI] + probeRadius;
            float ratio6       = powf( vdwRadius[atomI]/bornRadii[atomI], 6.0f );
            float saTerm       = surfaceAreaFactor*r*r*ratio6;
            _gbsaEnergy       += saTerm;
            bornForces[atomI]  = -6.0f*saTerm/bornRadii[atomI]; 
         }
      }

   } else {
      memset( bornForces, 0, arraySzInBytes );
   }

   // intialize Born radii w/ fixed portion of Gpol

   float* bornRadiiTemp             = getGbsaBornRadiiTemp();
   memcpy( bornRadiiTemp, gPolFixed, arraySzInBytes );

   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      float partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         // 3 FLOP

         float deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
         float deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
         float deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
         // 5 FLOP

         float r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;

         // 3 FLOP

         float alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         float D_ij               = r2/(4.0f*alpha2_ij);

         // exp + 2 + sqrt FLOP 

         float expTerm            = exp( -D_ij );
         float denominator2       = r2 + alpha2_ij*expTerm; 
         float denominator        = sqrt( denominator2 ); 
         
         // 6 FLOP

         float Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
  
         // dGpol/dr              = -1/2*(Gpol/denominator2)*(2r - r/2*exp() )
         float dGpol_dr           = -Gpol*( 1.0f - 0.25f*expTerm )/denominator2;  

         // 5 FLOP

         float dGpol_dalpha2_ij   = -0.5f*Gpol*expTerm*( 1.0f + D_ij )/denominator2;

         // 11 FLOP

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;
/*
             forces[0][atomI]  -= deltaX;
             forces[1][atomI]  -= deltaY;
             forces[2][atomI]  -= deltaZ;

             forces[0][atomJ]  += deltaX;
             forces[1][atomJ]  += deltaY;
             forces[2][atomJ]  += deltaZ;
*/
             forces[atomI][0]  -= deltaX;
             forces[atomI][1]  -= deltaY;
             forces[atomI][2]  -= deltaZ;

             forces[atomJ][0]  += deltaX;
             forces[atomJ][1]  += deltaY;
             forces[atomJ][2]  += deltaZ;

         } else {
            Gpol *= 0.5;
         }

         // 3 FLOP

         _gbsaEnergy       += Gpol;
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

      }
   }

   _gbsaEnergy *= getEnergyConversionFactor();

   // second main loop (born1 in Tinker esolv1.f)

   // chain terms from Born radii

   // precompute gpi terms and perform force conversion if necessary

   // 6 FLOP

   float forceFactor    = getForceConversionFactor();
   float constantFactor = P4/electricConstant;
   if( fabs(forceFactor - 1.0f) > 1.0e-04 ){
      constantFactor *= forceFactor;
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         bornForces[ii] *= constantFactor*bornRadii[ii]*bornRadii[ii];
         forces[ii][0]  *= forceFactor;
         forces[ii][1]  *= forceFactor;
         forces[ii][2]  *= forceFactor;
      }
   } else {
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         bornForces[ii] *= constantFactor*bornRadii[ii]*bornRadii[ii];
      }
   }

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // set exclusion factors

      getExclusionsForAtom( numberOfAtoms, gbsaBondsArray,
                            atomI, exclusionWorkArray, atomI - 1, log );

      for( int atomJ = atomI + 1; atomJ < numberOfAtoms; atomJ++ ){

         if( exclusionWorkArray[atomJ] != atomI ){

            // 3 FLOP

            float deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            float deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            float deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            // 7 FLOP

            float r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            float r4                 = r2*r2;
            float r6                 = r4*r2;

            // 3 FLOP

            float vdwSum             = vdwRadius[atomI] + vdwRadius[atomJ];
            float ratio              = r2/(vdwSum*vdwSum);

            float ccf;
            float dccf;
            float term;
            float bornFactor;

            // 0 or 8 + cos/sin

            if( ratio > P5Inverse ){
               ccf                   = 4.0f;
               dccf                  = 0.0f;
               term                  = 2.0f;
               bornFactor            = P4;
            } else {
               float theta           = ratio*piP5;
               term                  = 1.0f - cos(theta);
               ccf                   = term*term;
               dccf                  = term*sin(theta)*piP5*ratio;
               bornFactor            = 0.25f*P4*ccf;
            }

            // 9 FLOP

            float dE                 = (bornForces[atomI]*volume[atomJ] + bornForces[atomJ]*volume[atomI])*(ccf-dccf)/r6;
   
            deltaX                  *= dE;
            deltaY                  *= dE;
            deltaZ                  *= dE;

/*
            forces[0][atomI]        -= deltaX;
            forces[1][atomI]        -= deltaY;
            forces[2][atomI]        -= deltaZ;

            forces[0][atomJ]        += deltaX;
            forces[1][atomJ]        += deltaY;
            forces[2][atomJ]        += deltaZ;
*/
            // 12 FLOP

             forces[atomI][0]       -= deltaX;
             forces[atomI][1]       -= deltaY;
             forces[atomI][2]       -= deltaZ;

             forces[atomJ][0]       += deltaX;
             forces[atomJ][1]       += deltaY;
             forces[atomJ][2]       += deltaZ;

            bornRadiiTemp[atomI]    +=  (bornFactor*volume[atomJ])/r4;
            bornRadiiTemp[atomJ]    +=  (bornFactor*volume[atomI])/r4;
         }
      }
   }

   // finish up new Born radii and then copy 
   // updated values into 'permanent' array

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      if( bornRadiiTemp[atomI] != 0.0 ){
         bornRadiiTemp[atomI] = electricConstant/bornRadiiTemp[atomI];
      } else {
         bornRadiiTemp[atomI] = 0.1f;
      }
   }
   memcpy( bornRadii, bornRadiiTemp, arraySzInBytes );

   return forces;

}

/**---------------------------------------------------------------------------------------

   Write Born energy and forces (Simbios)

   @param top                 GMX t_topology struct
   @param atomCoordinates     atomic coordinates
   @param float forces        forces
   @param gbsaResultsFileName output file name
   @param log                 if set, then print error messages to log file

   @return 0 always

   --------------------------------------------------------------------------------------- */

int CpuGbsa::writeBornEnergyForces( const t_topology* top, const rvec *atomCoordinates,
                                    const float* partialCharges, float** forces,
                                    const char* gbsaResultsFileName, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuGbsa::writeBornEnergyForces";

   char buffer[MAX_BUFFER_SZ];

   // ---------------------------------------------------------------------------------------

   GbsaParameters* gbsaParameters = getGbsaParameters();
   int numberOfAtoms              = gbsaParameters->getNumberOfAtoms();
   float* vdwRadius               = gbsaParameters->getVdwRadii();
   float* volume                  = gbsaParameters->getVolume();
   float* gPolFixed               = gbsaParameters->getGPolFixed();
   float* bornRadii               = gbsaParameters->getBornRadii();

   // ---------------------------------------------------------------------------------------

   // char* gbsaResultsFileName = "111UBA.gbsa";
   
   FILE* gbsaResultsFile = fopen( gbsaResultsFileName, "w" );
   if( gbsaResultsFile != NULL ){
      if( log ){
         (void) fprintf( log, "%s Opened file=<%s>.", methodName, gbsaResultsFileName );
      }
   } else {
      if( log ){
         (void) fprintf( log, "%s could not open file=<%s> -- abort output.",
                         methodName, gbsaResultsFileName );
      }
      return -1;
   }

   (void) fprintf( gbsaResultsFile, "\n# %d atoms", numberOfAtoms );
   (void) fprintf( gbsaResultsFile, "\n#  Rdw Vol Gpol_fx Radii Q   Forces" );

   float forceConversion  = 0.01f;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      UtilitiesSimTk::getAtomIdStringGivenAtomIndex( ii, top, buffer, numberOfAtoms, ATOM_ID_STRING_TAB );
      (void) fprintf( gbsaResultsFile, "\n%d %.5e %.5e %.5e ",ii,
                      10.0f*vdwRadius[ii], 1000.0f*volume[ii],
                      0.1f*gPolFixed[ii] );
      if( forces != NULL && atomCoordinates != NULL ){
         (void) fprintf( gbsaResultsFile, "%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e ",
                         10.0f*bornRadii[ii], partialCharges[ii],
                         10.0f*atomCoordinates[ii][0],
                         10.0f*atomCoordinates[ii][1], 
                         10.0f*atomCoordinates[ii][2],
                         forceConversion*forces[0][ii],
                         forceConversion*forces[1][ii],
                         forceConversion*forces[2][ii]
                       );
      }
      (void) fprintf( gbsaResultsFile, "%s", buffer );
   }
   (void) fclose( gbsaResultsFile );

   return 0;

}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   @param top                 GMX t_topology struct
   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      atom partial charges
   @param energy              energy
   @param debugOn             enable debug
   @param debugAtomI          debug flag for atomI (outer loop atom)
   @param debugAtomJ          debug flag for atomJ (inner loop atom)
   @param debugReport         flag signalling what quantities are to be saved for comparisons
   @param saveLoopForces      flag signalling whether intermediate results for each loop
                              are to be retained
   @param printOn             print flag 
   @param unsetDebugValue     ?
   @param numberOfDebugStreams 
                              number of debug streams
   @param debugStreams        array of debug streams
   @param forces              force array (may be null)
   @param log                 if set, then print error messages to log file

   @return 0 always; fixed value for G_pol in array 'gPolFixed'

   --------------------------------------------------------------------------------------- */

float** CpuGbsa::computeBornEnergyForcesDebug( const t_topology* top,
                                               float* bornRadii, const rvec *atomCoordinates,
                                               const float* partialCharges, float* energy,
                                               bool debugOn, int debugAtomI, int debugAtomJ,
                                               int debugReport, bool* saveLoopForces,
                                               bool printOn, float unsetDebugValue, int numberOfDebugStreams,
                                               float** debugStreams, float** forces, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuGbsa::computeBornEnergyForces";

   char buffer[MAX_BUFFER_SZ];
   static int debugBufferAllocated = 0;
   static char* debugBuffer[MAX_DEBUG_FIELDS];

   static int floatValueBlockAllocated = 0;
   static float* floatValues[2] = { NULL, NULL };

   // debug stuff

   // int atomTrack               = 438;  // N 27 ALA
   int atomTrack               = -442;    // CB 27 ALA

   // ---------------------------------------------------------------------------------------

   // prototype not working for debugBuffer[MAX_DEBUG_FIELDS][MAX_BUFFER_SZ] -- silly workaround

   if( atomTrack > 0 && !debugBufferAllocated ){
      debugBufferAllocated  = 1;
      for( int ii = 0; ii < MAX_DEBUG_FIELDS; ii++ ){
#ifdef UseGromacsMalloc
         debugBuffer[ii] = (char*) save_malloc( "debugBuffer", __FILE__, __LINE__, sizeof( char )*MAX_BUFFER_SZ );
#else
         debugBuffer[ii] = new char[MAX_BUFFER_SZ];
#endif
      }
   }

   // ---------------------------------------------------------------------------------------

   GbsaParameters* gbsaParameters = getGbsaParameters();

   int numberOfAtoms              = gbsaParameters->getNumberOfAtoms();
   // FILE* log                      = gbsaParameters->getLog();
   float* vdwRadius               = gbsaParameters->getVdwRadii();
   float* volume                  = gbsaParameters->getVolume();
   // float* gPolFixed               = gbsaParameters->getGPolFixed();
   if( bornRadii == NULL ){
      bornRadii               = gbsaParameters->getBornRadii();
   }
   gbsaBonds** gbsaBondsArray     = gbsaParameters->getGbsaBondsArray( );

   // ---------------------------------------------------------------------------------------

   // constants

   float preFactor               = gbsaParameters->getPreFactor();
   float P5Inverse               = gbsaParameters->getP5Inverse();
   float P4                      = gbsaParameters->getP4();
   float P5                      = gbsaParameters->getP5();
   float piP5                    = gbsaParameters->getPiP5();
   float electricConstant        = gbsaParameters->getElectricConstant();

   // ---------------------------------------------------------------------------------------

   // exclusion work array
 
   int* exclusionWorkArray = gbsaParameters->getExclusionWorkArray();

   // prototype not working for debugBuffer[MAX_DEBUG_FIELDS][MAX_BUFFER_SZ] -- silly workaround
   if( !floatValueBlockAllocated ){
      floatValueBlockAllocated = 1;
      for( int ii = 0; ii < 2; ii++ ){
#ifdef UseGromacsMalloc
         floatValues[ii] = (float*) save_malloc( "floatValues", __FILE__, __LINE__, sizeof( float )*10 );
#else
         floatValues[ii] = new float[10];
#endif
      }
   }

   // ---------------------------------------------------------------------------------------

   // allocate force array if not allocated and initialize value to 0

   int numberForceBlocks         = 4;
   int numberForceBlockIndex[3];
   for( int ii = 0; ii < 3; ii++ ){
      if( saveLoopForces[ii] ){
         numberForceBlockIndex[ii]  = numberForceBlocks;
         numberForceBlocks         += 4;
      } else {
         numberForceBlockIndex[ii]  = -1;
      }
      if( log ){
         (void) fprintf( log, "\nForceBlocks: %d idx=%d total=%d", ii, numberForceBlockIndex[ii], numberForceBlocks );
      }
   }
   forces                        = UtilitiesSimTk::allocate2DFloatArray( numberForceBlocks, numberOfAtoms,
                                                                         forces, true, 0.0f, log );
   float* bornForces             = forces[3];
   const int jUnroll             = 4;

   // ---------------------------------------------------------------------------------------

   int numberOfDebugFields     = 5;
   float debugFields[30];

   FILE* debugFile             = NULL;
   char* debugFileName         = "111UBQ.dbg0";
   float energyConversion      = 0.2388f;
   float forceConversion       = energyConversion*0.1f;

   if( atomTrack > -1 ){
      UtilitiesSimTk::getAtomIdStringGivenAtomIndex( atomTrack, top, buffer, numberOfAtoms, ATOM_ID_STRING_TAB );
      (void) sprintf( debugBuffer[0], "# Inner 0 debugAtom=%d %s [%.4f %.4f %.4f] add 1 to get Tinker id if files synched",
                      atomTrack, buffer, atomCoordinates[atomTrack][0], atomCoordinates[atomTrack][1],
                      atomCoordinates[atomTrack][2] );
      debugFile = UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, debugBuffer[0],
                                  debugFileName, WriteDebugFile, debugFile, log );
   }

   // ---------------------------------------------------------------------------------------

/*
     if (use_gbsa .and. solvtyp.ne.'ONION') then
         term = 4.0d0 * pi
         do i = 1, n
            ai = asolv(i)
            ri = rsolv(i)
            rb = rborn(i)
            if (rb .ne. 0.0d0) then
               e = ai * term * (ri+probe)**2 * (ri/rb)**6
               es = es + e
               drb(i) = drb(i) - 6.0d0*e/rb
            end if 
         end do
     end if
*/

   // compute the nonpolar solvation via ACE approximation
/*
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      float r            = vdwRadius[atomI] + probe;
      float ratio6       = powf( vdwRadius[atomI]/bornRadii[atomI], 6.0 );
      float saTerm       = pi4Asolv*r*r*ratio6;
      // *energy           += saTerm;
      bornForces[atomI]  = -saTerm/bornRadii[atomI]; 
   }
*/

   // ---------------------------------------------------------------------------------------

   // intialize debug stream

   if( debugOn ){
      UtilitiesSimTk::initialize2DFloatArray( numberOfDebugStreams, numberOfAtoms, debugStreams, unsetDebugValue, log );
   }

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // gpu-cpu comparison

      if( debugOn && debugReport == 20 ){

         float* floatPtr  = (float*) &(debugStreams[0][jUnroll*atomI]);
         *floatPtr        = atomCoordinates[atomI][0]; floatPtr++;
         *floatPtr        = atomCoordinates[atomI][1]; floatPtr++;
         *floatPtr        = atomCoordinates[atomI][2]; floatPtr++;
         *floatPtr        = partialCharges[atomI];

         floatPtr         = (float*) &(debugStreams[1][jUnroll*atomI]);
         *floatPtr        = bornRadii[atomI]; floatPtr++;
         *floatPtr        = vdwRadius[atomI]; floatPtr++;
         // *floatPtr      = bornForces[atomI]; floatPtr++;
         *floatPtr        = 0.0f; floatPtr++;
         *floatPtr        = 0.0f;

      }
                
      // ---------------------------------------------------------------------------------------

      // first main loop

      float partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         float deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
         float deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
         float deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
 
         float r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;

         float alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         float D_ij               = r2/(4.0f*alpha2_ij);

         float expTerm            = exp( -D_ij );
         float denominator2       = r2 + alpha2_ij*expTerm; 
         float denominator        = sqrt( denominator2 ); 
         
         float Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
  
         // dGpol/dr              = -1/2*(Gpol/denominator2)*(2r - r/2*exp() )
         float dGpol_dr           = -Gpol*( 1.0f - 0.25f*expTerm )/denominator2;  

         float dGpol_dalpha2_ij   = -0.5f*Gpol*expTerm*( 1.0f + D_ij )/denominator2;

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;

             forces[0][atomI]  -= deltaX;
             forces[1][atomI]  -= deltaY;
             forces[2][atomI]  -= deltaZ;

             forces[0][atomJ]  += deltaX;
             forces[1][atomJ]  += deltaY;
             forces[2][atomJ]  += deltaZ;

         } else {
            Gpol *= 0.5;
         }

         *energy           += Gpol;
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

/*
fprintf( log, "\nBFF %d %d r2=%.3f Gpol=%.3e denominator2=%.3e D_ij=%.3e expTerm=%.3e dGpol_dalpha2_ij=%.3e"
              " Br=%.3e Bf=%.3e, preFactor=%.3e pQI=%.4f pQJ=%.3e", atomI, atomJ, r2, Gpol,
              denominator2, D_ij, expTerm, dGpol_dalpha2_ij, bornRadii[atomJ], bornForces[atomI],
              preFactor, partialCharges[atomI], partialCharges[atomJ] );
*/

         // ---------------------------------------------------------------------------------------

         // collect terms for comparison w/ gpu

         if( debugOn && debugReport > 20 && debugReport < 30 ){

            bool useUpper = false;
            for( int kk = 0; kk < jUnroll; kk++ ){
               if( debugReport == 21 ){
                  floatValues[0][kk] = r2;
                  floatValues[1][kk] = bornRadii[debugAtomJ + kk];
               } else if( debugReport == 22 ){
                  floatValues[0][kk] = dGpol_dr;
                  floatValues[1][kk] = Gpol;
               } else if( debugReport == 23 ){
                  floatValues[0][kk] = deltaX;
                  floatValues[1][kk] = deltaY;
               }
            }
            UtilitiesSimTk::storeInGpuFormat( atomI, atomJ, debugAtomJ,
                                              jUnroll, 2, debugStreams, floatValues,
                                              unsetDebugValue, useUpper, log );
         }

         // ---------------------------------------------------------------------------------------

         // Tinker comparisons

         if( atomTrack > -1 && (atomI == atomTrack || atomJ == atomTrack) ){

               int atomToPrint;
               if( atomI == atomTrack ){
                  atomToPrint = atomJ;
               } else {
                  atomToPrint = atomI;
               }
               if( atomI == atomJ ){
                  deltaX = 0.0f;
                  deltaY = 0.0f;
                  deltaZ = 0.0f;
               }
         
               numberOfDebugFields                     = 0;
               debugFields[numberOfDebugFields++]      = (float) atomToPrint;
         
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][0];
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][1];
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][2];
         
               debugFields[numberOfDebugFields++]      = -forceConversion*deltaX;
               debugFields[numberOfDebugFields++]      = -forceConversion*deltaY;
               debugFields[numberOfDebugFields++]      = -forceConversion*deltaZ;
         
               debugFields[numberOfDebugFields++]      = 10.0f*bornRadii[atomToPrint];
               debugFields[numberOfDebugFields++]      = D_ij;
         
               debugFields[numberOfDebugFields++]      = 10.0f*sqrt( r2 );
               debugFields[numberOfDebugFields++]      = 100.0f*alpha2_ij;
               debugFields[numberOfDebugFields++]      = 10.0f*denominator;
               debugFields[numberOfDebugFields++]      = energyConversion*Gpol;
               debugFields[numberOfDebugFields++]      = 0.1f*forceConversion*dGpol_dr;
               debugFields[numberOfDebugFields++]      = 0.1f*forceConversion*dGpol_dalpha2_ij;
               debugFields[numberOfDebugFields++]      = expTerm;
               debugFields[numberOfDebugFields++]      = (-dGpol_dalpha2_ij*denominator2)/(Gpol*expTerm);

               debugFields[numberOfDebugFields++]      = forceConversion*dGpol_dalpha2_ij*bornRadii[atomToPrint];
               debugFields[numberOfDebugFields++]      = partialCharges[atomToPrint];
        
               UtilitiesSimTk::getAtomIdStringGivenAtomIndex( atomToPrint, top, debugBuffer[0], numberOfAtoms, ATOM_ID_STRING_TAB );

               debugFile = UtilitiesSimTk::writeDebugFile( numberOfDebugFields, debugFields, 1, debugBuffer, NULL,
                                           NULL, WriteDebugFile, debugFile, log );
         }

   // ---------------------------------------------------------------------------------------

      }

   }

   // ---------------------------------------------------------------------------------------

   // close debug file

   if( atomTrack > -1 ){
      UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, NULL, debugFileName, CloseDebugFile, debugFile, log );
   }

   // ---------------------------------------------------------------------------------------

   if( atomTrack > -1 ){
      writeBornEnergyForces( top, atomCoordinates, partialCharges, forces, "111UBQ.gbsa0", log );
   }

   // ---------------------------------------------------------------------------------------

   if( atomTrack > -1 ){

      debugFile                   = NULL;
      debugFileName               = "111UBQ.dbg1";

      (void) sprintf( debugBuffer[0], "# debugAtom=%d E=%.5e %s", atomTrack, *energy, buffer );
      debugFile  = UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, debugBuffer[0],
                                   debugFileName, WriteDebugFile, debugFile, log );

      int stretchCount;
      int angleCount;
      for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

         // sum no. exclusion atoms

         gbsaParameters->getBondCounts( gbsaBondsArray[atomI], &stretchCount, &angleCount );
         getExclusionsForAtom( numberOfAtoms, gbsaBondsArray,
                               atomI, exclusionWorkArray, atomI-1, log );

         int sum  =  0;
         for( int ii = 0; ii < numberOfAtoms; ii++ ){
            if( exclusionWorkArray[ii] == atomI ){
               sum += 1;
            }
         }
         UtilitiesSimTk::getAtomIdStringGivenAtomIndex( atomI, top, debugBuffer[0], numberOfAtoms, ATOM_ID_STRING_TAB );

         numberOfDebugFields                     = 0;
   
         debugFields[numberOfDebugFields++]      = (float) atomI;

         debugFields[numberOfDebugFields++]      = (float) stretchCount;
         debugFields[numberOfDebugFields++]      = (float) angleCount;
         debugFields[numberOfDebugFields++]      = (float) sum;
   
         debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomI][0];
         debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomI][1];
         debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomI][2];

         debugFields[numberOfDebugFields++]      = 10.0f*bornRadii[atomI];
         debugFields[numberOfDebugFields++]      = forceConversion*bornForces[atomI];

         debugFile = UtilitiesSimTk::writeDebugFile( numberOfDebugFields, debugFields, 1, debugBuffer, NULL,
                                     NULL, WriteDebugFile, debugFile, log );
      }

      // close file

      UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, NULL,
                      debugFileName, CloseDebugFile, debugFile, log );
   }

   // ---------------------------------------------------------------------------------------

   if( atomTrack > -1 ){

      debugFile                   = NULL;
      debugFileName               = "111UBQ.dbg2";

      (void) sprintf( debugBuffer[0], "# debugAtom=%d E=%.5e %s", atomTrack, *energy, buffer );
      debugFile  = UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, debugBuffer[0],
                                   debugFileName, WriteDebugFile, debugFile, log );

   }

   // ---------------------------------------------------------------------------------------

// if( atomTrack > -1 ){
//   return forces;
// }

   // save forces

   if( numberForceBlockIndex[0] >= 0 && saveLoopForces[0] ){
      int index = numberForceBlockIndex[0];
      unsigned int totalSz = sizeof( float )*numberOfAtoms;
      for( int jj = 0; jj < 4; jj++ ){
         memcpy( forces[index+jj], forces[jj], totalSz );
      }
   }

   if( printOn && log ){
      (void) fprintf( log, "%s CPU bornForces to 1", methodName );
      for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
         (void) fprintf( log, "\n%d [%12.4e %12.4e %12.4e %12.4e]",
                         atomI, forces[0][atomI], forces[1][atomI], forces[2][atomI], bornForces[atomI] );
   // bornForces[atomI] = 1.0f;
      }
   }

/*

// zero forces -- debugging

(void) fprintf( log, "%s CPU zeroed forces", methodName );
initialize2DFloatArray( 3, numberOfAtoms, forces, 0.0f, log );
(void) fflush( log );
*/

   // ---------------------------------------------------------------------------------------

   // second main loop (born1 in Tinker esolv1.f)

   // chain terms from Born radii

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // i-specific factors

      float gpi =  (bornForces[atomI]*P4*bornRadii[atomI]*bornRadii[atomI])/electricConstant;
      if( debugOn && debugReport == 30 ){

         float* floatPtr = (float*) &(debugStreams[0][jUnroll*atomI]);
         *floatPtr       = atomCoordinates[atomI][0]; floatPtr++;
         *floatPtr       = atomCoordinates[atomI][1]; floatPtr++;
         *floatPtr       = atomCoordinates[atomI][2]; floatPtr++;
         *floatPtr       = gpi;

         floatPtr         = (float*) &(debugStreams[1][jUnroll*atomI]);
         *floatPtr      = vdwRadius[atomI]; floatPtr++;
         *floatPtr      = volume[atomI]; floatPtr++;

         int direction    = ( atomI % 2 ) ? -1 : 1;
         *floatPtr      = (float)
                            ( (gbsaBondsArray[atomI]->startExclusionIndex < gbsaBondsArray[atomI+direction]->startExclusionIndex) ?
                               gbsaBondsArray[atomI]->startExclusionIndex : gbsaBondsArray[atomI+direction]->startExclusionIndex );
         floatPtr++;

         *floatPtr      = (float) 
                            ( (gbsaBondsArray[atomI]->stopExclusionIndex > gbsaBondsArray[atomI+direction]->stopExclusionIndex ) ?
                               gbsaBondsArray[atomI]->stopExclusionIndex : gbsaBondsArray[atomI+direction]->stopExclusionIndex );
      }

      // set exclusion factors

      getExclusionsForAtom( numberOfAtoms, gbsaBondsArray,
                            atomI, exclusionWorkArray, atomI - 1, log );

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( exclusionWorkArray[atomJ] != atomI ){

            float deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            float deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            float deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            float r2                 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
            float r6                 = r2*r2*r2;

            float vdwSum             = vdwRadius[atomI] + vdwRadius[atomJ];
            float ratio              = r2/(vdwSum*vdwSum);

            float ccf;
            float dccf;
            float term;
            if( ratio > P5Inverse ){
               ccf              = 4.0f;
               dccf             = 0.0f;
               term             = 2.0f;
            } else {
               float theta      = ratio*piP5;
               term             = 1.0f - cos(theta);
               ccf              = term*term;
               dccf             = term*sin(theta)*piP5*ratio;
            }

            bool oldWay         = false;
            float dE;
            if( oldWay ){
               dE               = (gpi*volume[atomJ])*(ccf-dccf)/r6;
            } else {
               float gpj        =  (bornForces[atomJ]*P4*bornRadii[atomJ]*bornRadii[atomJ])/electricConstant;
               dE               = (gpi*volume[atomJ] + gpj*volume[atomI])*(ccf-dccf)/r6;
            }
   
            deltaX             *= dE;
            deltaY             *= dE;
            deltaZ             *= dE;

            forces[0][atomI]   -= deltaX;
            forces[1][atomI]   -= deltaY;
            forces[2][atomI]   -= deltaZ;

            if( debugOn && debugReport != 30 ){
               if( oldWay ){
                  forces[0][atomJ]   += deltaX;
                  forces[1][atomJ]   += deltaY;
                  forces[2][atomJ]   += deltaZ;
               }
            }
   
           // save forces (Loop 2 only)

           if( numberForceBlockIndex[1] >= 0 && saveLoopForces[1] ){
               int index                = numberForceBlockIndex[1];
               forces[index][atomI]    -= deltaX;
               forces[index+1][atomI]  -= deltaY;
               forces[index+2][atomI]  -= deltaZ;
           }

           // save 'nonexcluded' forces from Loop 2

           if(   numberForceBlockIndex[2] >= 0 && saveLoopForces[2] &&
               ( atomJ >=  gbsaBondsArray[atomI]->startExclusionIndex &&
                 atomJ <   gbsaBondsArray[atomI]->stopExclusionIndex ) ){
               int index                = numberForceBlockIndex[2];
               forces[index][atomI]    -= deltaX;
               forces[index+1][atomI]  -= deltaY;
               forces[index+2][atomI]  -= deltaZ;
               // QQ1
               if( atomI == debugAtomI ){
/*
                  (void) fprintf( log, "\n%d %d Add=[%.4f %.4f %.4f] Ttl=[%.4f %.4f %.4f] [%d %d]",
                                  atomI, atomJ, -deltaX, -deltaY, -deltaZ,
                                  forces[index][atomI], forces[index+1][atomI], forces[index+2][atomI],
                                  gbsaBondsArray[atomI]->startExclusionIndex,
                                  gbsaBondsArray[atomI]->stopExclusionIndex );
*/
                  if( log ){
                     (void) fprintf( log, "\n%d %d Add=[%.4f %.4f %.4f] de=%.4f ratio=%.4f term=%.4f"
                                     " dccf=%.4f vdWS=%.4f r2=%.4f TtlF=[%.4f %.4f %.4f] [%d %d]",
                                     atomI, atomJ, -deltaX, -deltaY, -deltaZ,
                                     dE, P5*ratio, term, dccf, vdwSum, r2,
                                     forces[index][atomI], forces[index+1][atomI], forces[index+2][atomI],
                                     gbsaBondsArray[atomI]->startExclusionIndex,
                                     gbsaBondsArray[atomI]->stopExclusionIndex );
                  }
               }
           }

            // ---------------------------------------------------------------------------------------

            // collect terms for comparison w/ gpu
   
            if( debugOn && debugReport >= 31 && debugReport < 34 &&

                atomJ >= debugAtomJ && atomJ < (debugAtomJ+jUnroll) ){
   
               bool useUpper = true;
               if( debugReport == 31 ){
                  // useUpper = false; 
                  for( int kk = 0; kk < jUnroll; kk++ ){
                     floatValues[0][kk] = r2;
                     floatValues[1][kk] = dE;
                  }
// fprintf( log, "\n%d %d r2=%.3f dE=%.3f %d Br=%.3e Bf=%.3e", atomI, atomJ, r2, dE, debugReport, bornRadii[atomJ], bornForces[atomJ] );

               } else if( debugReport == 32 ){

                  for( int kk = 0; kk < jUnroll; kk++ ){
                     floatValues[0][kk] = ccf;
                     floatValues[1][kk] = dccf;
                  }

               } else if( debugReport == 33 ){

                  for( int kk = 0; kk < jUnroll; kk++ ){
                     if( debugAtomJ+kk < numberOfAtoms ){
                        floatValues[0][kk] = volume[debugAtomJ+kk];
                        floatValues[1][kk] = gpi;
                     } else {
                        floatValues[0][kk] = volume[debugAtomJ];
                        // floatValues[0][kk] = (float) atomJ;
                        floatValues[1][kk] = (bornRadii[debugAtomJ]*bornRadii[debugAtomJ])/electricConstant;
                     }
                  }
               }
               UtilitiesSimTk::storeInGpuFormat( atomI, atomJ, debugAtomJ,
                                                 jUnroll, 2, debugStreams, floatValues,
                                                 unsetDebugValue, useUpper, log );
            }
   
            // ---------------------------------------------------------------------------------------

            // Tinker comparison

            if( atomTrack > -1 && atomI == atomTrack ){
   
               int atomToPrint                         = atomJ;
            
               numberOfDebugFields                     = 0;
               debugFields[numberOfDebugFields++]      = (float) atomToPrint;
            
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][0];
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][1];
               debugFields[numberOfDebugFields++]      = 10.0f*atomCoordinates[atomToPrint][2];
            
               // debugFields[numberOfDebugFields++]      = -forceConversion*deltaX;
               // debugFields[numberOfDebugFields++]      = -forceConversion*deltaY;
               // debugFields[numberOfDebugFields++]      = -forceConversion*deltaZ;
            
               debugFields[numberOfDebugFields++]      = -0.01f*deltaX;
               debugFields[numberOfDebugFields++]      = -0.01f*deltaY;
               debugFields[numberOfDebugFields++]      = -0.01f*deltaZ;
            
               debugFields[numberOfDebugFields++]      = ratio;
               debugFields[numberOfDebugFields++]      = 10.0f*sqrt( r2 );
               debugFields[numberOfDebugFields++]      = ccf/4.0f;
               debugFields[numberOfDebugFields++]      = dccf;

               // debugFields[numberOfDebugFields++]      = energyConversion*dE;
               debugFields[numberOfDebugFields++]      = -0.001f*dE;

               debugFields[numberOfDebugFields++]      = P5Inverse;
               debugFields[numberOfDebugFields++]      = 1000.0f*volume[atomJ];

               // float gpi =  (bornForces[atomI]*P4*bornRadii[atomI]*bornRadii[atomI])/electricConstant;
               float localGpi                          = gpi;
               if( bornForces[atomI] > 0.0f ){
                  localGpi /= bornForces[atomI];
               }
               localGpi    /= P4;
               debugFields[numberOfDebugFields++]      = -localGpi*100.0f;

               debugFields[numberOfDebugFields++]      = 10.0f*bornRadii[atomI];
               debugFields[numberOfDebugFields++]      = -localGpi/(bornRadii[atomI]*bornRadii[atomI]);
            
               UtilitiesSimTk::getAtomIdStringGivenAtomIndex( atomToPrint, top, debugBuffer[0],
                                                              numberOfAtoms, ATOM_ID_STRING_TAB );
   
               debugFile = UtilitiesSimTk::writeDebugFile( numberOfDebugFields, debugFields, 1, debugBuffer, NULL,
                                           NULL, WriteDebugFile, debugFile, log );
            }

   // ---------------------------------------------------------------------------------------

         }
      }
   }

   // close debug file

   if( debugFile ){
      UtilitiesSimTk::writeDebugFile( 0, NULL, 0, NULL, NULL,
                                      debugFileName, CloseDebugFile, debugFile, log );
   }

   // ---------------------------------------------------------------------------------------

   if( atomTrack > -1 ){
      writeBornEnergyForces( top, atomCoordinates, partialCharges, forces, "111UBQ.gbsa1", log );
   }

   // ---------------------------------------------------------------------------------------

   return forces;

}

/**---------------------------------------------------------------------------------------
      
   Print state to log file (Simbios)
   
   @param title               title (optional)
   @param log                 print state to log file
      
   @return 0 always
      
   --------------------------------------------------------------------------------------- */

int CpuGbsa::logState( const char* title, FILE* log ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGbsa::logState";

   // ---------------------------------------------------------------------------------------

   if( log == NULL ){
      return 0;
   }
   if( title ){
      (void) fprintf( log, "%s", title );
   }

   (void) fprintf( log, "\nGbsa info:" );
   (void) fprintf( log, "\nForce conversion=%12.4f Energy conversion=%12.4f",
                   _forceConversionFactor, _energyConversionFactor ); 

   (void) fprintf( log, "\nInclude ACE approximation=%s", _includeAceApproximation ? "Y" : "N" ); 

   // log parameters

   getGbsaParameters()->logState( NULL, log );

   return 0;
}

/**---------------------------------------------------------------------------------------

   Write Tinker xyz file (Simbios)

   @param numberOfAtoms        number of atoms
   @param atomCoordinates      atom coordinates
   @param header               header
   @param xyzFileName          output file name
   @param gbsaBondsArray       bond array -- used to print 1-2 bonds
   @param top                  Gromac's topology struct
   @param log                  if set, then print error messages to log file

   @return 0 unless error detected

   Currently no attempt is made to get the atom name/type to accurately 
   reflect the Tinker names/types. Rather method is used to output atoms
   in Gromacs order and then reorder those in a corresponding xyz file
   w/ the correct atom names/types so that they match the Gromacs order
   This makes it easier to compare results between Gromacs and Tinker

   --------------------------------------------------------------------------------------- */

int CpuGbsa::writeXyzFile( int numberOfAtoms, const rvec* atomCoordinates, 
                           const char* header, const char* xyzFileName,
                           gbsaBonds** gbsaBondsArray, const t_topology* top, FILE* log ){

   // ---------------------------------------------------------------------------------------

   char atomName[32];
   static const char* methodName = "\nCpuGbsa::writeXyzFile";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* xyzFile = fopen( xyzFileName, "w" );
   if( xyzFile != NULL ){
      if( log != NULL ){
         (void) fprintf( log, "%s Opened file=<%s>.", methodName, xyzFileName );
      }
   } else {
      if( log != NULL ){
         (void) fprintf( log, "%s could not open file=<%s> -- abort output.",
                         methodName, xyzFileName );
         (void) fflush( log );
      }
      return -1;
   }

   // first line

   (void) fprintf( xyzFile, "%d", numberOfAtoms );
   if( header != NULL ){
      (void) fprintf( xyzFile, " %s", header );
   } else {
      (void) fprintf( xyzFile, " Unknown" );
   }
   (void) fprintf( xyzFile, "\n" );

/*
  1232  CHROMOSOMAL PROTEIN                     02-JAN-87   1UBQ 
     1  N3    27.340000   24.430000    2.614000   472     2     5     6     7    
     2  CT    26.266000   25.413000    2.842000   473     1     3     8     9    
     3  C     26.913000   26.639000    3.531000   474     2     4    20   
*/

   int type = 1;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      // truncate atoms names to at most 3 chars

      UtilitiesSimTk::getAtomNameGivenAtomIndex( ii, atomName, top );
      if( strlen( atomName ) > 3 ){
         atomName[3] = '\0';
      }

      // scale coordinates by 10 for Angstrom -> nanometer conversion

      (void) fprintf( xyzFile, "%6d  %-3s %16.9f %16.9f %16.9f %6d ", ii+1, atomName,
                      10.0f*atomCoordinates[ii][0], 10.0f*atomCoordinates[ii][1], 10.0f*atomCoordinates[ii][2],
                      type );

      // include 1-2 bonds

      gbsaStretchBond* stretchBond = gbsaBondsArray[ii]->stretchBonds;
      while( stretchBond != NULL ){
         if( stretchBond->atomI != ii ){
            (void) fprintf( xyzFile, "%6d ", stretchBond->atomI + 1 );
         }
         if( stretchBond->atomJ != ii ){
            (void) fprintf( xyzFile, "%6d ", stretchBond->atomJ + 1 );
         }
         stretchBond = stretchBond->nextBond;
      }
      (void) fprintf( xyzFile, "\n" );
   }
   (void) fflush( xyzFile );
   (void) fclose( xyzFile );

   if( log != NULL ){
      (void) fprintf( log, "%s Closed file=<%s> atoms=%d.", methodName, xyzFileName, numberOfAtoms );
   }

   return 0;
}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

void* CpuGbsa::operator new( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

/*
   (void) fprintf( stdout, "\nCpuGbsa new called -- size=%u", size );
   (void) fflush( stdout );
*/

   return ptr;
}

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

void CpuGbsa::operator delete( void *ptr ){

/*
   (void) fprintf( stdout, "\nCpuGbsa delete called." );
   (void) fflush( stdout );
*/

   sfree( ptr ); 
}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

void* CpuGbsa::operator new[]( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

/*
   (void) fprintf( stdout, "\nCpuGbsa new[] called -- size=%u", size );
   (void) fflush( stdout );
*/

   return ptr;
}

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

void CpuGbsa::operator delete[]( void *ptr ){

/*
   (void) fprintf( stdout, "\nCpuGbsa delete[] called." );
   (void) fflush( stdout );
*/

   sfree( ptr ); 
}
