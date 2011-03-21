
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

#include <string.h>
#include <stdlib.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "CpuObc.h"
#include "../SimTKReference/ReferenceForce.h"
#include <cmath>
#include <cstdio>

using namespace OpenMM;
using namespace std;

/**---------------------------------------------------------------------------------------

   CpuObc constructor

   obcParameters      obcParameters object
   
   --------------------------------------------------------------------------------------- */

CpuObc::CpuObc( ImplicitSolventParameters* obcParameters ) : CpuImplicitSolvent( obcParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::CpuObc";

   // ---------------------------------------------------------------------------------------

   _initializeObcDataMembers( );

   _obcParameters = static_cast<ObcParameters*> (obcParameters);

}

/**---------------------------------------------------------------------------------------

   CpuObc destructor

   --------------------------------------------------------------------------------------- */

CpuObc::~CpuObc( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::~CpuObc";

   // ---------------------------------------------------------------------------------------
}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuObc::_initializeObcDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _obcParameters = NULL;
}

/**---------------------------------------------------------------------------------------

   Get ObcParameters reference

   @return ObcParameters reference

   --------------------------------------------------------------------------------------- */

ObcParameters* CpuObc::getObcParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcParameters";

   // ---------------------------------------------------------------------------------------

   return _obcParameters;
}

/**---------------------------------------------------------------------------------------

   Set ObcParameters reference

   @param ObcParameters reference

   --------------------------------------------------------------------------------------- */

void CpuObc::setObcParameters(  ObcParameters* obcParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::setObcParameters";

   // ---------------------------------------------------------------------------------------

   _obcParameters = obcParameters;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuObc::getObcChain( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcChain";

   // ---------------------------------------------------------------------------------------

   if( _obcChain.size() == 0 ){
      _obcChain.resize(_obcParameters->getNumberOfAtoms());
   }
   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()

   @return array

   --------------------------------------------------------------------------------------- */

const vector<RealOpenMM>& CpuObc::getObcChainConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getObcChain";

   // ---------------------------------------------------------------------------------------

   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain temp work array of size=_obcParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuObc::getObcChainTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::getImplicitSolventObcChainTemp";

   // ---------------------------------------------------------------------------------------

   return _obcChainTemp;
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on papers:

      J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
      Proteins: Structure, Function, and Bioinformatcis 55:383-394 (2004) (OBC paper)

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii

   --------------------------------------------------------------------------------------- */

void CpuObc::computeBornRadii( vector<RealVec>& atomCoordinates, vector<RealOpenMM>& bornRadii ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;

   static const char* methodName = "\nCpuObc::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   ObcParameters* obcParameters          = getObcParameters();

   int numberOfAtoms                     = obcParameters->getNumberOfAtoms();
   RealOpenMM* atomicRadii               = obcParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();
   vector<RealOpenMM>& obcChain          = getObcChain();

   RealOpenMM dielectricOffset           = obcParameters->getDielectricOffset();
   RealOpenMM alphaObc                   = obcParameters->getAlphaObc();
   RealOpenMM betaObc                    = obcParameters->getBetaObc();
   RealOpenMM gammaObc                   = obcParameters->getGammaObc();

   // ---------------------------------------------------------------------------------------

   // calculate Born radii

//FILE* logFile = SimTKOpenMMLog::getSimTKOpenMMLogFile( );
//FILE* logFile = NULL;
//FILE* logFile = fopen( "bR", "w" );

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
     
      RealOpenMM radiusI         = atomicRadii[atomI];
      RealOpenMM offsetRadiusI   = radiusI - dielectricOffset;

      RealOpenMM radiusIInverse  = one/offsetRadiusI;
      RealOpenMM sum             = zero;

      // HCT code

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
            if (_obcParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR );
            else
                ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
            RealOpenMM r               = deltaR[ReferenceForce::RIndex];
            if (_obcParameters->getUseCutoff() && r > _obcParameters->getCutoffDistance())
                continue;
            RealOpenMM offsetRadiusJ   = atomicRadii[atomJ] - dielectricOffset; 
            RealOpenMM scaledRadiusJ   = offsetRadiusJ*scaledRadiusFactor[atomJ];
            RealOpenMM rScaledRadiusJ  = r + scaledRadiusJ;

            if( offsetRadiusI < rScaledRadiusJ ){
               RealOpenMM rInverse = one/r;
               RealOpenMM l_ij     = offsetRadiusI > FABS( r - scaledRadiusJ ) ? offsetRadiusI : FABS( r - scaledRadiusJ );
                          l_ij     = one/l_ij;

               RealOpenMM u_ij     = one/rScaledRadiusJ;

               RealOpenMM l_ij2    = l_ij*l_ij;
               RealOpenMM u_ij2    = u_ij*u_ij;
 
               RealOpenMM ratio    = LN( (u_ij/l_ij) );
               RealOpenMM term     = l_ij - u_ij + fourth*r*(u_ij2 - l_ij2)  + ( half*rInverse*ratio) + (fourth*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2 - u_ij2);

               // this case (atom i completely inside atom j) is not considered in the original paper
               // Jay Ponder and the authors of Tinker recognized this and
               // worked out the details

               if( offsetRadiusI < (scaledRadiusJ - r) ){
                  term += two*( radiusIInverse - l_ij);
               }
               sum += term;

/*
if( logFile && atomI == 0 ){
   (void) fprintf( logFile, "\nRR %d %d r=%.4f rads[%.6f %.6f] scl=[%.3f %.3f] sum=%12.6e %12.6e %12.6e %12.6e",
                   atomI, atomJ, r, offsetRadiusI, offsetRadiusJ, scaledRadiusFactor[atomI], scaledRadiusFactor[atomJ], 0.5f*sum,
                   l_ij, u_ij, term );
}
*/

            }
         }
      }
 
      // OBC-specific code (Eqs. 6-8 in paper)

      sum                  *= half*offsetRadiusI;
      RealOpenMM sum2       = sum*sum;
      RealOpenMM sum3       = sum*sum2;
      RealOpenMM tanhSum    = TANH( alphaObc*sum - betaObc*sum2 + gammaObc*sum3 );
      
      bornRadii[atomI]      = one/( one/offsetRadiusI - tanhSum/radiusI ); 
 
      obcChain[atomI]       = offsetRadiusI*( alphaObc - two*betaObc*sum + three*gammaObc*sum2 );
      obcChain[atomI]       = (one - tanhSum*tanhSum)*obcChain[atomI]/radiusI;

/*
if( logFile ){
   (void) fprintf( logFile, "\nRRQ %d sum %12.6e tanhS %12.6e radI %.5f %.5f born %18.10e obc %12.6e",
                   atomI, sum, tanhSum, radiusI, offsetRadiusI, bornRadii[atomI], obcChain[atomI] );
}
*/

   }
/*
if( logFile ){
      (void) fclose( logFile );
} */
}

/**---------------------------------------------------------------------------------------

   Get Obc Born energy and forces

   @param bornRadii           Born radii -- optional; if NULL, then ObcParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   The array bornRadii is also updated and the obcEnergy

   --------------------------------------------------------------------------------------- */

void CpuObc::computeBornEnergyForces( vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                     const RealOpenMM* partialCharges, vector<RealVec>& inputForces ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObc::computeBornEnergyForces";

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM four    = (RealOpenMM) 4.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;
   static const RealOpenMM eighth  = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   const ObcParameters* obcParameters = getObcParameters();
   const int numberOfAtoms            = obcParameters->getNumberOfAtoms();

   // ---------------------------------------------------------------------------------------

   // constants

   const RealOpenMM preFactor           = obcParameters->getPreFactor();
   const RealOpenMM dielectricOffset    = obcParameters->getDielectricOffset();

   // ---------------------------------------------------------------------------------------

#if 0
{
   RealOpenMM* atomicRadii               = obcParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();
   RealOpenMM* obcChain                  = getObcChain();
   FILE* logFile = fopen( "bornParameters", "w" );
   (void) fprintf( logFile, "%5d dielOff=%.4e rad::hct::q::bR::Chain::coords\n", numberOfAtoms, dielectricOffset );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "%5d %10.5f %10.5f %10.5f %14.7e %14.7e %14.7e %14.7e %14.7e\n", atomI,
                      atomicRadii[atomI], scaledRadiusFactor[atomI], partialCharges[atomI], bornRadii[atomI], obcChain[atomI],
                      atomCoordinates[atomI][0], atomCoordinates[atomI][1], atomCoordinates[atomI][2] );
   }
   (void) fclose( logFile );
}
#endif

   // set energy/forces to zero

   RealOpenMM obcEnergy                 = zero;

   RealOpenMM** forces  = (RealOpenMM**) malloc( sizeof( RealOpenMM* )*numberOfAtoms );
   RealOpenMM*  block   = (RealOpenMM*)  malloc( sizeof( RealOpenMM )*numberOfAtoms*3 );
	memset( block, 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
	RealOpenMM* blockPtr = block;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      forces[ii] = blockPtr;
		blockPtr  += 3;
   }

   vector<RealOpenMM>& bornForces = getBornForce();
   bornForces.assign(numberOfAtoms, 0.0);

   // ---------------------------------------------------------------------------------------

   // N*( 8 + pow) ACE
   // compute the nonpolar solvation via ACE approximation
    
   if( includeAceApproximation() ){
      computeAceNonPolarForce( obcParameters, bornRadii, &obcEnergy, bornForces );
   }
 
   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      RealOpenMM partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
         if (_obcParameters->getPeriodic())
             ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR );
         else
             ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
         if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _obcParameters->getCutoffDistance())
             continue;
         RealOpenMM r2                 = deltaR[ReferenceForce::R2Index];
         RealOpenMM deltaX             = deltaR[ReferenceForce::XIndex];
         RealOpenMM deltaY             = deltaR[ReferenceForce::YIndex];
         RealOpenMM deltaZ             = deltaR[ReferenceForce::ZIndex];

         // 3 FLOP

         RealOpenMM alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         RealOpenMM D_ij               = r2/(four*alpha2_ij);

         // exp + 2 + sqrt FLOP 

         RealOpenMM expTerm            = EXP( -D_ij );
         RealOpenMM denominator2       = r2 + alpha2_ij*expTerm; 
         RealOpenMM denominator        = SQRT( denominator2 ); 
         
         // 6 FLOP

         RealOpenMM Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
         RealOpenMM dGpol_dr           = -Gpol*( one - fourth*expTerm )/denominator2;  

         // 5 FLOP

         RealOpenMM dGpol_dalpha2_ij   = -half*Gpol*expTerm*( one + D_ij )/denominator2;

         // 11 FLOP

         if( atomI != atomJ ){

             bornForces[atomJ] += dGpol_dalpha2_ij*bornRadii[atomI];

             deltaX            *= dGpol_dr;
             deltaY            *= dGpol_dr;
             deltaZ            *= dGpol_dr;

             forces[atomI][0]  += deltaX;
             forces[atomI][1]  += deltaY;
             forces[atomI][2]  += deltaZ;

             forces[atomJ][0]  -= deltaX;
             forces[atomJ][1]  -= deltaY;
             forces[atomJ][2]  -= deltaZ;

         } else {
            Gpol *= half;
         }

         // 3 FLOP

         obcEnergy         += Gpol;
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

      }
   }


   //obcEnergy *= getEnergyConversionFactor();

   // ---------------------------------------------------------------------------------------

   // second main loop

   // initialize Born radii & ObcChain temp arrays -- contain values
   // used in next iteration

   vector<RealOpenMM>& bornRadiiTemp     = getBornRadiiTemp();
   bornRadiiTemp.assign(numberOfAtoms, 0.0);

   vector<RealOpenMM>& obcChainTemp      = getObcChainTemp();
   obcChainTemp.assign(numberOfAtoms, 0.0);

   vector<RealOpenMM>& obcChain          = getObcChain();
   const RealOpenMM* atomicRadii         = obcParameters->getAtomicRadii();

   const RealOpenMM alphaObc             = obcParameters->getAlphaObc();
   const RealOpenMM betaObc              = obcParameters->getBetaObc();
   const RealOpenMM gammaObc             = obcParameters->getGammaObc();
   const RealOpenMM* scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();

    // compute factor that depends only on the outer loop index

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      bornForces[atomI] *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
   }

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // radius w/ dielectric offset applied

      RealOpenMM radiusI        = atomicRadii[atomI];
      RealOpenMM offsetRadiusI  = radiusI - dielectricOffset;

      // used to compute Born radius for next iteration

      RealOpenMM bornSum        = zero;

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
            if (_obcParameters->getPeriodic())
               ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR );
            else 
               ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
            if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _obcParameters->getCutoffDistance())
                   continue;
   
            RealOpenMM deltaX             = deltaR[ReferenceForce::XIndex];
            RealOpenMM deltaY             = deltaR[ReferenceForce::YIndex];
            RealOpenMM deltaZ             = deltaR[ReferenceForce::ZIndex];
            RealOpenMM r                  = deltaR[ReferenceForce::RIndex];
 
            // radius w/ dielectric offset applied

            RealOpenMM offsetRadiusJ      = atomicRadii[atomJ] - dielectricOffset;

            RealOpenMM scaledRadiusJ      = offsetRadiusJ*scaledRadiusFactor[atomJ];
            RealOpenMM scaledRadiusJ2     = scaledRadiusJ*scaledRadiusJ;
            RealOpenMM rScaledRadiusJ     = r + scaledRadiusJ;

            // dL/dr & dU/dr are zero (this can be shown analytically)
            // removed from calculation

            if( offsetRadiusI < rScaledRadiusJ ){

               RealOpenMM l_ij          = offsetRadiusI > FABS( r - scaledRadiusJ ) ? offsetRadiusI : FABS( r - scaledRadiusJ );
                    l_ij                = one/l_ij;

               RealOpenMM u_ij          = one/rScaledRadiusJ;

               RealOpenMM l_ij2         = l_ij*l_ij;

               RealOpenMM u_ij2         = u_ij*u_ij;
 
               RealOpenMM rInverse      = one/r;
               RealOpenMM r2Inverse     = rInverse*rInverse;

               RealOpenMM t3            = eighth*(one + scaledRadiusJ2*r2Inverse)*(l_ij2 - u_ij2) + fourth*LN( u_ij/l_ij )*r2Inverse;

               RealOpenMM de            = bornForces[atomI]*t3*rInverse;

               deltaX                  *= de;
               deltaY                  *= de;
               deltaZ                  *= de;
   
               forces[atomI][0]        -= deltaX;
               forces[atomI][1]        -= deltaY;
               forces[atomI][2]        -= deltaZ;
  
               forces[atomJ][0]        += deltaX;
               forces[atomJ][1]        += deltaY;
               forces[atomJ][2]        += deltaZ;
 
               // Born radius term

               RealOpenMM term          =  l_ij - u_ij  + fourth*r*(u_ij2 - l_ij2) + (half*rInverse)*LN(u_ij/l_ij)   +
                                           (fourth*scaledRadiusJ*scaledRadiusJ*rInverse)*(l_ij2-u_ij2);

               if( offsetRadiusI < (scaledRadiusJ - r) ){
                  term += two*( (one/offsetRadiusI) - l_ij);
               }
               bornSum += term; 
            }
         }
      }

      // OBC-specific code (Eqs. 6-8 in paper)

      bornSum                   *= half*offsetRadiusI;
      RealOpenMM sum2            = bornSum*bornSum;
      RealOpenMM sum3            = bornSum*sum2;
      RealOpenMM tanhSum         = TANH( alphaObc*bornSum - betaObc*sum2 + gammaObc*sum3 );
      
      bornRadiiTemp[atomI]       = one/( one/offsetRadiusI - tanhSum/radiusI );
 
      obcChainTemp[atomI]        = offsetRadiusI*( alphaObc - two*betaObc*bornSum + three*gammaObc*sum2 );
      obcChainTemp[atomI]        = (one - tanhSum*tanhSum)*obcChainTemp[atomI]/radiusI;
   }

   // cal to Joule conversion

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      inputForces[atomI][0] += forces[atomI][0];
      inputForces[atomI][1] += forces[atomI][1];
      inputForces[atomI][2] += forces[atomI][2];
   }
   setEnergy( obcEnergy );

   // copy new Born radii and obcChain values into permanent array

   bornRadii = bornRadiiTemp;
   obcChain = obcChainTemp;

    free( (char*) block );
    free( (char*) forces );
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      
   @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuObc::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << CpuImplicitSolvent::getStateString( title );

   return message.str();
}
