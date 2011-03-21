
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
#include <sstream>
#include <stdlib.h>

//#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "CpuObcSoftcore.h"
#include "../SimTKReference/ReferenceForce.h"
#include <cmath>
#include <cstdio>
#include <vector>

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   CpuObcSoftcore constructor

   obcSoftcoreParameters      obcSoftcoreParameters object
   
   --------------------------------------------------------------------------------------- */

CpuObcSoftcore::CpuObcSoftcore( ImplicitSolventParameters* obcSoftcoreParameters ) : CpuImplicitSolvent( obcSoftcoreParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::CpuObcSoftcore";

   // ---------------------------------------------------------------------------------------

   _initializeObcDataMembers( );

   _obcSoftcoreParameters = static_cast<ObcSoftcoreParameters*> (obcSoftcoreParameters);

}

/**---------------------------------------------------------------------------------------

   CpuObcSoftcore destructor

   --------------------------------------------------------------------------------------- */

CpuObcSoftcore::~CpuObcSoftcore( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::~CpuObcSoftcore";

   // ---------------------------------------------------------------------------------------
}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuObcSoftcore::_initializeObcDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _obcSoftcoreParameters = NULL;
}

/**---------------------------------------------------------------------------------------

   Get ObcSoftcoreParameters reference

   @return ObcSoftcoreParameters reference

   --------------------------------------------------------------------------------------- */

ObcSoftcoreParameters* CpuObcSoftcore::getObcSoftcoreParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::getObcSoftcoreParameters";

   // ---------------------------------------------------------------------------------------

   return _obcSoftcoreParameters;
}

/**---------------------------------------------------------------------------------------

   Set ObcSoftcoreParameters reference

   @param ObcSoftcoreParameters reference

   --------------------------------------------------------------------------------------- */

void CpuObcSoftcore::setObcSoftcoreParameters(  ObcSoftcoreParameters* obcSoftcoreParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::setObcSoftcoreParameters";

   // ---------------------------------------------------------------------------------------

   _obcSoftcoreParameters = obcSoftcoreParameters;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcSoftcoreParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuObcSoftcore::getObcChain( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::getObcChain";

   // ---------------------------------------------------------------------------------------

   if( _obcChain.size() == 0 ){
      _obcChain.resize(_obcSoftcoreParameters->getNumberOfAtoms());
   }
   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcSoftcoreParameters->getNumberOfAtoms()

   @return array

   --------------------------------------------------------------------------------------- */

const vector<RealOpenMM>& CpuObcSoftcore::getObcChainConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::getObcChain";

   // ---------------------------------------------------------------------------------------

   return _obcChain;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain temp work array of size=_obcSoftcoreParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuObcSoftcore::getObcChainTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::getImplicitSolventObcChainTemp";

   // ---------------------------------------------------------------------------------------

   if( _obcChainTemp.size() == 0 ){
      _obcChainTemp.resize(_obcSoftcoreParameters->getNumberOfAtoms());
   }
   return _obcChainTemp;
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on papers:

      J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
      Proteins: Structure, Function, and Bioinformatcis 55:383-394 (2004) (OBC paper)

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii

   --------------------------------------------------------------------------------------- */

void CpuObcSoftcore::computeBornRadii( vector<RealVec>& atomCoordinates, RealOpenMM* bornRadii ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;

   static const char* methodName   = "\nCpuObcSoftcore::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   ObcSoftcoreParameters* obcSoftcoreParameters             = getObcSoftcoreParameters();

   int numberOfAtoms                        = obcSoftcoreParameters->getNumberOfAtoms();
   RealOpenMM* atomicRadii                  = obcSoftcoreParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor     = obcSoftcoreParameters->getScaledRadiusFactors();
   vector<RealOpenMM>& obcChain             = getObcChain();

   const RealOpenMM* nonPolarScaleFactors   = obcSoftcoreParameters->getNonPolarScaleFactors();
   RealOpenMM dielectricOffset              = obcSoftcoreParameters->getDielectricOffset();
   RealOpenMM alphaObc                      = obcSoftcoreParameters->getAlphaObc();
   RealOpenMM betaObc                       = obcSoftcoreParameters->getBetaObc();
   RealOpenMM gammaObc                      = obcSoftcoreParameters->getGammaObc();

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
            if (_obcSoftcoreParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcSoftcoreParameters->getPeriodicBox(), deltaR );
            else
                ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
            RealOpenMM r               = deltaR[ReferenceForce::RIndex];
            if (_obcSoftcoreParameters->getUseCutoff() && r > _obcSoftcoreParameters->getCutoffDistance())
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
               sum += nonPolarScaleFactors[atomJ]*term;

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

#if 0
if( logFile && atomI >= 0 ){
   (void) fprintf( logFile, "\nRRQ %d sum %12.6e tanhS %12.6e radI %.5f %.5f born %18.10e obc %12.6e",
                   atomI, sum, tanhSum, radiusI, offsetRadiusI, bornRadii[atomI], obcChain[atomI] );
}
#endif

   }

#if 0
if( logFile ){
      (void) fclose( logFile );
}
#endif
}

/**---------------------------------------------------------------------------------------

   Get nonpolar solvation force constribution via ACE approximation

   @param obcSoftcoreParameters     parameters
   @param vdwRadii                  Vdw radii
   @param bornRadii                 Born radii
   @param energy                    energy (output): value is incremented from input value 
   @param forces                    forces: values are incremented from input values

   --------------------------------------------------------------------------------------- */

void CpuObcSoftcore::computeAceNonPolarForce( const ObcSoftcoreParameters* obcSoftcoreParameters,
                                             const vector<RealOpenMM>& bornRadii, RealOpenMM* energy,
                                             vector<RealOpenMM>& forces ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::computeAceNonPolarForce";

   static const RealOpenMM minusSix = -6.0;

   // ---------------------------------------------------------------------------------------

   // compute the nonpolar solvation via ACE approximation

   const RealOpenMM probeRadius           = obcSoftcoreParameters->getProbeRadius();
   const RealOpenMM surfaceAreaFactor     = obcSoftcoreParameters->getPi4Asolv();

   const RealOpenMM* atomicRadii          = obcSoftcoreParameters->getAtomicRadii();
   const RealOpenMM* nonPolarScaleFactors = obcSoftcoreParameters->getNonPolarScaleFactors();
   int numberOfAtoms                      = obcSoftcoreParameters->getNumberOfAtoms();

   // 1 + 1 + pow + 3 + 1 + 2 FLOP

   // the original ACE equation is based on Eq.2 of

   // M. Schaefer, C. Bartels and M. Karplus, "Solution Conformations
   // and Thermodynamics of Structured Peptides: Molecular Dynamics
   // Simulation with an Implicit Solvation Model", J. Mol. Biol.,
   // 284, 835-848 (1998)  (ACE Method)

   // The original equation includes the factor (atomicRadii[atomI]/bornRadii[atomI]) to the first power,
   // whereas here the ratio is raised to the sixth power: (atomicRadii[atomI]/bornRadii[atomI])**6

   // This modification was made by Jay Ponder who observed it gave better correlations w/
   // observed values. He did not think it was important enough to write up, so there is
   // no paper to cite.

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      if( bornRadii[atomI] > 0.0 ){
         RealOpenMM r            = atomicRadii[atomI] + probeRadius;
         RealOpenMM ratio6       = POW( atomicRadii[atomI]/bornRadii[atomI], (RealOpenMM) 6.0 );
         RealOpenMM saTerm       = nonPolarScaleFactors[atomI]*surfaceAreaFactor*r*r*ratio6;
         *energy                += saTerm;
         forces[atomI]          += minusSix*saTerm/bornRadii[atomI]; 
      }
   }
}

/**---------------------------------------------------------------------------------------

   Get Obc Born energy and forces

   @param bornRadii           Born radii -- optional; if NULL, then ObcSoftcoreParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   The array bornRadii is also updated and the obcEnergy

   --------------------------------------------------------------------------------------- */

void CpuObcSoftcore::computeBornEnergyForces( vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                             const RealOpenMM* partialCharges, vector<RealVec>& inputForces ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuObcSoftcore::computeBornEnergyForces";

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM four    = (RealOpenMM) 4.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;
   static const RealOpenMM eighth  = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   const ObcSoftcoreParameters* obcSoftcoreParameters = getObcSoftcoreParameters();
   const int numberOfAtoms            = obcSoftcoreParameters->getNumberOfAtoms();

   // ---------------------------------------------------------------------------------------

   // constants

   const RealOpenMM preFactor           = obcSoftcoreParameters->getPreFactor();
   const RealOpenMM dielectricOffset    = obcSoftcoreParameters->getDielectricOffset();

   // ---------------------------------------------------------------------------------------

#if 0
{
   RealOpenMM* atomicRadii               = obcSoftcoreParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = obcSoftcoreParameters->getScaledRadiusFactors();
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
      computeAceNonPolarForce( obcSoftcoreParameters, bornRadii, &obcEnergy, bornForces );
   }
 
   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      RealOpenMM partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
         if (_obcSoftcoreParameters->getPeriodic())
             ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcSoftcoreParameters->getPeriodicBox(), deltaR );
         else
             ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
         if (_obcSoftcoreParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _obcSoftcoreParameters->getCutoffDistance())
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
   const RealOpenMM* atomicRadii         = obcSoftcoreParameters->getAtomicRadii();

   const RealOpenMM alphaObc             = obcSoftcoreParameters->getAlphaObc();
   const RealOpenMM betaObc              = obcSoftcoreParameters->getBetaObc();
   const RealOpenMM gammaObc             = obcSoftcoreParameters->getGammaObc();
   const RealOpenMM* scaledRadiusFactor  = obcSoftcoreParameters->getScaledRadiusFactors();

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
            if (_obcSoftcoreParameters->getPeriodic())
               ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcSoftcoreParameters->getPeriodicBox(), deltaR );
            else 
               ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
            if (_obcSoftcoreParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _obcSoftcoreParameters->getCutoffDistance())
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

#if 0
{
   RealOpenMM* atomicRadii               = obcSoftcoreParameters->getAtomicRadii();
   const RealOpenMM* scaledRadiusFactor  = obcSoftcoreParameters->getScaledRadiusFactors();
   RealOpenMM* obcChain                  = getObcChain();
   //FILE* logFile = fopen( "bornParameters", "w" );
   FILE* logFile = stderr;
   (void) fprintf( logFile, "%5d dielOff=%.4e rad::hct::q::bR::Chain::bF::f::coords\n", numberOfAtoms, dielectricOffset );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "%5d %10.5f %10.5f q=%10.5f b[%14.7e %14.7e %14.7e] f[%14.7e %14.7e %14.7e] x[%14.7e %14.7e %14.7e]\n", atomI,
                      atomicRadii[atomI], scaledRadiusFactor[atomI], partialCharges[atomI], bornRadii[atomI], obcChain[atomI],
                      conversion*bornForces[atomI], 
                      conversion*forces[atomI][0], conversion*forces[atomI][1], conversion*forces[atomI][2],
                        atomCoordinates[atomI][0],   atomCoordinates[atomI][1],   atomCoordinates[atomI][2] );
   }
   if( logFile != stderr || logFile != stdout ){
      (void) fclose( logFile );
   }
}
#endif

   // copy new Born radii and obcChain values into permanent array

   bornRadii = bornRadiiTemp;
   obcChain = obcChainTemp;

    free( (char*) block );
    free( (char*) forces );
}
