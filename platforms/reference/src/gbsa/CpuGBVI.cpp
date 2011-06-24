
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

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "CpuGBVI.h"
#include "../SimTKReference/ReferenceForce.h"
#include <math.h>

using namespace std;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   CpuGBVI constructor

   gbviParameters      gbviParameters object
   
   --------------------------------------------------------------------------------------- */

CpuGBVI::CpuGBVI( ImplicitSolventParameters* gbviParameters ) : CpuImplicitSolvent( gbviParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::CpuGBVI";

   // ---------------------------------------------------------------------------------------

   _initializeGBVIDataMembers( );

   _gbviParameters = static_cast<GBVIParameters*> (gbviParameters);

}

/**---------------------------------------------------------------------------------------

   CpuGBVI destructor

   --------------------------------------------------------------------------------------- */

CpuGBVI::~CpuGBVI( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::~CpuGBVI";

   // ---------------------------------------------------------------------------------------

   //if( _gbviParameters != NULL ){
     // delete _gbviParameters;
   //}

}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuGBVI::_initializeGBVIDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _gbviParameters  = NULL;
}

/**---------------------------------------------------------------------------------------

   Get GBVIParameters reference

   @return GBVIParameters reference

   --------------------------------------------------------------------------------------- */

GBVIParameters* CpuGBVI::getGBVIParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::getGBVIParameters";

   // ---------------------------------------------------------------------------------------

   return _gbviParameters;
}

/**---------------------------------------------------------------------------------------

   Set GBVIParameters reference

   @param GBVIParameters reference

   --------------------------------------------------------------------------------------- */

void CpuGBVI::setGBVIParameters(  GBVIParameters* gbviParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::setGBVIParameters";

   // ---------------------------------------------------------------------------------------

   _gbviParameters = gbviParameters;
   if( _switchDeriviative.size() <= 0 ){
      _switchDeriviative.resize(_gbviParameters->getNumberOfAtoms());
   }
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

std::vector<RealOpenMM>& CpuGBVI::getSwitchDeriviative( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::getSwitchDeriviative";

   // ---------------------------------------------------------------------------------------

   if( _switchDeriviative.size() <= 0 ){
      _switchDeriviative.resize(_gbviParameters->getNumberOfAtoms());
   }
   return _switchDeriviative;
}

/**---------------------------------------------------------------------------------------

   Return switching function derivative

   @return array

   --------------------------------------------------------------------------------------- */

/*
std::vector<RealOpenMM> CpuGBVI::getSwitchDeriviativeConst( void ) const {
   return _switchDeriviative;
}

*/

/**---------------------------------------------------------------------------------------

   Compute quintic spline value and associated derviative

   @param x                   value to compute spline at
   @param rl                  lower cutoff value
   @param ru                  upper cutoff value
   @param outValue            value of spline at x
   @param outDerivative       value of derivative of spline at x

   --------------------------------------------------------------------------------------- */

#define GBVIDebug 0

void CpuGBVI::quinticSpline( RealOpenMM x, RealOpenMM rl, RealOpenMM ru,
                             RealOpenMM* outValue, RealOpenMM* outDerivative ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM one           = (RealOpenMM)   1.0;
   static const RealOpenMM minusSix      = (RealOpenMM)  -6.0;
   static const RealOpenMM minusTen      = (RealOpenMM) -10.0;
   static const RealOpenMM minusThirty   = (RealOpenMM) -30.0;
   static const RealOpenMM fifteen       = (RealOpenMM)  15.0;
   static const RealOpenMM sixty         = (RealOpenMM)  60.0;

   // static const char* methodName         = "CpuGBVI::quinticSpline";

   // ---------------------------------------------------------------------------------------

   RealOpenMM numerator    = x  - rl;
   RealOpenMM denominator  = ru - rl;
   RealOpenMM ratio        = numerator/denominator;
   RealOpenMM ratio2       = ratio*ratio;
   RealOpenMM ratio3       = ratio2*ratio;

   *outValue               = one + ratio3*(minusTen + fifteen*ratio + minusSix*ratio2);
   *outDerivative          = ratio2*(minusThirty + sixty*ratio + minusThirty*ratio2)/denominator;
}

/**---------------------------------------------------------------------------------------

   Compute Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])
   and quintic splice switching function

   @param atomicRadius3       atomic radius cubed
   @param bornSum             Born sum (volume integral)
   @param gbviParameters      Gbvi parameters (parameters used in spline
                              QuinticLowerLimitFactor & QuinticUpperBornRadiusLimit)
   @param bornRadius          output Born radius
   @param switchDeriviative   output switching function deriviative

   --------------------------------------------------------------------------------------- */

#define GBVIDebug 0

void CpuGBVI::computeBornRadiiUsingQuinticSpline( RealOpenMM atomicRadius3, RealOpenMM bornSum,
                                                  GBVIParameters* gbviParameters, 
                                                  RealOpenMM* bornRadius, RealOpenMM* switchDeriviative ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero          = (RealOpenMM)  0.0;
   static const RealOpenMM one           = (RealOpenMM)  1.0;
   static const RealOpenMM minusOne      = (RealOpenMM) -1.0;
   static const RealOpenMM minusThree    = (RealOpenMM) -3.0;
   static const RealOpenMM oneEighth     = (RealOpenMM)  0.125;
   static const RealOpenMM minusOneThird = (RealOpenMM) (-1.0/3.0);
   static const RealOpenMM three         = (RealOpenMM)  3.0;

   static const char* methodName         = "CpuGBVI::computeBornRadiiUsingQuinticSpline";

   // ---------------------------------------------------------------------------------------

#if( GBVIDebug == 1 )
   FILE* logFile                         = stderr;
#endif

   // R                = [ S(V)*(A - V) ]**(-1/3)

   // S(V)             = 1                                 V < L
   // S(V)             = qSpline + U/(A-V)                 L < V < A
   // S(V)             = U/(A-V)                           U < V 

   // dR/dr            = (-1/3)*[ S(V)*(A - V) ]**(-4/3)*[ d{ S(V)*(A-V) }/dr

   // d{ S(V)*(A-V) }/dr   = (dV/dr)*[ (A-V)*dS/dV - S(V) ]

   //  (A - V)*dS/dV - S(V)  = 0 - 1                             V < L

   //  (A - V)*dS/dV - S(V)  = (A-V)*d(qSpline) + (A-V)*U/(A-V)**2 - qSpline - U/(A-V) 

	//                        = (A-V)*d(qSpline) - qSpline        L < V < A**(-3)

   //  (A - V)*dS/dV - S(V)  = (A-V)*U*/(A-V)**2 - U/(A-V) = 0   U < V

   RealOpenMM splineL          = gbviParameters->getQuinticLowerLimitFactor()*atomicRadius3;
   RealOpenMM sum;
   if( bornSum > splineL ){
      if( bornSum < atomicRadius3 ){
         RealOpenMM splineValue, splineDerivative;
         quinticSpline( bornSum, splineL, atomicRadius3, &splineValue, &splineDerivative ); 
         sum                 = (atomicRadius3 - bornSum)*splineValue + gbviParameters->getQuinticUpperBornRadiusLimit();
         *switchDeriviative  = splineValue - (atomicRadius3 - bornSum)*splineDerivative;
#if( GBVIDebug == 1 )
      (void) fprintf( logFile, " Qv=%14.6e splnDrvtv=%14.6e spline[%10.3e %10.3e] ", splineValue, splineDerivative,
                      splineL, gbviParameters->getQuinticUpperBornRadiusLimit() );
#endif
      } else {   
         sum                = gbviParameters->getQuinticUpperBornRadiusLimit();
         *switchDeriviative = zero;
      }
   } else {
      sum                = atomicRadius3 - bornSum; 
      *switchDeriviative = one;
   }
   *bornRadius = POW( sum, minusOneThird );
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii
   @param chain               not used here

   --------------------------------------------------------------------------------------- */

#define GBVIDebug 0

void CpuGBVI::computeBornRadii( vector<RealVec>& atomCoordinates, vector<RealOpenMM>& bornRadii ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero          = (RealOpenMM) 0.0;
   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM minusThree    = (RealOpenMM) -3.0;
   static const RealOpenMM oneEighth     = (RealOpenMM) 0.125;
   static const RealOpenMM minusOneThird = (RealOpenMM) (-1.0/3.0);
   static const RealOpenMM three         = (RealOpenMM) 3.0;

   static const char* methodName         = "CpuGBVI::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   GBVIParameters* gbviParameters                   = getGBVIParameters();
   int numberOfAtoms                                = gbviParameters->getNumberOfAtoms();
   RealOpenMM* atomicRadii                          = gbviParameters->getAtomicRadii();
   const RealOpenMM* scaledRadii                    = gbviParameters->getScaledRadii();

   std::vector<RealOpenMM>& switchDeriviatives      = getSwitchDeriviative();

   // ---------------------------------------------------------------------------------------

#if( GBVIDebug == 1 )
   FILE* logFile                         = stderr;
   (void) fprintf( logFile, "\n%s scalingMethdd=%d\n", methodName,  _gbviParameters->getBornRadiusScalingMethod() );
   (void) fprintf( logFile, "Lower=%14.6e Upper=%10.3e\n", gbviParameters->getQuinticLowerLimitFactor(), gbviParameters->getQuinticUpperBornRadiusLimit() );
#endif

   // calculate Born radii

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
     
      RealOpenMM radiusI         = atomicRadii[atomI];
      RealOpenMM sum             = zero;

      // sum over volumes

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
            if (_gbviParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _gbviParameters->getPeriodicBox(), deltaR );
            else
                ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );

            RealOpenMM r               = deltaR[ReferenceForce::RIndex];

            if (_gbviParameters->getUseCutoff() && r > _gbviParameters->getCutoffDistance())
                continue;

            sum  += CpuGBVI::getVolume( r, radiusI, scaledRadii[atomJ] );

#if( GBVIDebug == 1 )
if( atomI == 1568 || atomJ == 1568 ){
            (void) fprintf( logFile, "%d addJ=%d scR=%14.6e %14.6e sum=%14.6e rI=%14.6e r=%14.6e S-R=%14.6e\n",
                            atomI, atomJ, scaledRadii[atomJ], getVolume( r, radiusI, scaledRadii[atomJ] ), sum,
                            radiusI, r, (scaledRadii[atomJ]-radiusI) );
}
#endif
         }
      }

#if( GBVIDebug == 1 )
      (void) fprintf( logFile, "%d Born radius sum=%14.6e %14.6e %14.6e q=%d ", 
                      atomI, sum, POW( radiusI, minusThree ), (POW( radiusI, minusThree ) - sum), GBVIParameters::QuinticSpline );
#endif

      RealOpenMM atomicRadius3 = POW( radiusI, minusThree );
      if( _gbviParameters->getBornRadiusScalingMethod() == GBVIParameters::QuinticSpline ){
         sum                        = atomicRadius3 - sum; 
         bornRadii[atomI]           = POW( sum, minusOneThird );
         switchDeriviatives[atomI]  = one; 
      } else {
         RealOpenMM bornRadius, switchDeriviative;
         computeBornRadiiUsingQuinticSpline( atomicRadius3, sum, gbviParameters, 
                                             &bornRadius, &switchDeriviative );
         bornRadii[atomI]           = bornRadius;
         switchDeriviatives[atomI]  = switchDeriviative;
      }    

#if( GBVIDebug == 1 )
      (void) fprintf( logFile, " br=%14.6e swDeriv=%14.6e\n", atomI, bornRadii[atomI], switchDeriviatives[atomI] );
#endif

   }

#if( GBVIDebug == 2 )
   std::string fileName  = "br.txt";
   FILE* brFile          = NULL;
#ifdef WIN32
   fopen_s( &brFile, fileName.c_str(), "w" );
#else
   brFile = fopen( fileName.c_str(), "w" );
#endif

   (void) fprintf( brFile, "%5d\n", numberOfAtoms );
   for( unsigned int ii = 0; ii < numberOfAtoms; ii++ ){
       (void) fprintf( brFile, "%6u %14.6e %9.4f %9.4f %14.6e %14.6e %14.6e\n",
                       ii, bornRadii[ii], atomicRadii[ii], scaledRadii[ii], 
                       atomCoordinates[ii][0], atomCoordinates[ii][1], atomCoordinates[ii][2] );
   }
   (void) fclose( brFile );
#endif

}

#undef GBVIDebug

/**---------------------------------------------------------------------------------------

   Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return volume

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::getVolume( RealOpenMM r, RealOpenMM R, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::getVolume";

   static const RealOpenMM zero         = (RealOpenMM)  0.0;
   static const RealOpenMM minusThree   = (RealOpenMM) -3.0;

   RealOpenMM              diff         = (S - R);
   if( FABS( diff ) < r ){

      RealOpenMM lowerBound = (R > (r - S)) ? R : (r - S);
      return (CpuGBVI::getL( r, (r + S),    S ) -
              CpuGBVI::getL( r, lowerBound, S ));

   } else if( r <= diff ){

      return CpuGBVI::getL( r, (r + S), S ) -
             CpuGBVI::getL( r, (r - S), S ) + 
             POW( R, minusThree );

   } else {
      return zero;
   }
}

/**---------------------------------------------------------------------------------------

   Get L (used in analytical solution for volume integrals) 

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return L value (Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::getL( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::getL";

   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM threeHalves   = (RealOpenMM) 1.5;
   static const RealOpenMM third         = (RealOpenMM) (1.0/3.0);
   static const RealOpenMM fourth        = (RealOpenMM) 0.25;
   static const RealOpenMM eighth        = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   RealOpenMM rInv   = one/r;

   RealOpenMM xInv   = one/x;
   RealOpenMM xInv2  = xInv*xInv;
   RealOpenMM xInv3  = xInv2*xInv;

   RealOpenMM diff2  = (r + S)*(r - S);

   return (threeHalves*xInv)*( (xInv*fourth*rInv) - (xInv2*third) + (diff2*xInv3*eighth*rInv) );
}

/**---------------------------------------------------------------------------------------

   Get partial derivative of L wrt r

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::dL_dr( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::dL_dr";

   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM threeHalves   = (RealOpenMM) 1.5;
   static const RealOpenMM threeEights   = (RealOpenMM) 0.375;
   static const RealOpenMM third         = (RealOpenMM) (1.0/3.0);
   static const RealOpenMM fourth        = (RealOpenMM) 0.25;
   static const RealOpenMM eighth        = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   RealOpenMM rInv   = one/r;
   RealOpenMM rInv2  = rInv*rInv;

   RealOpenMM xInv   = one/x;
   RealOpenMM xInv2  = xInv*xInv;
   RealOpenMM xInv3  = xInv2*xInv;

   RealOpenMM diff2  = (r + S)*(r - S);

   return ( (-threeHalves*xInv2*rInv2)*( fourth + eighth*diff2*xInv2 ) + threeEights*xInv3*xInv );
}

/**---------------------------------------------------------------------------------------

   Get partial derivative of L wrt x

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::dL_dx( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::dL_dx";

   static const RealOpenMM one           = (RealOpenMM)  1.0;
   static const RealOpenMM half          = (RealOpenMM)  0.5;
   static const RealOpenMM threeHalvesM  = (RealOpenMM) -1.5;
   static const RealOpenMM third         = (RealOpenMM)  (1.0/3.0);

   // ---------------------------------------------------------------------------------------

   RealOpenMM rInv   = one/r;

   RealOpenMM xInv   = one/x;
   RealOpenMM xInv2  = xInv*xInv;
   RealOpenMM xInv3  = xInv2*xInv;

   RealOpenMM diff   = (r + S)*(r - S);

   return (threeHalvesM*xInv3)*( (half*rInv) - xInv + (half*diff*xInv2*rInv) );
}

/**---------------------------------------------------------------------------------------

   Sgb function

   @param t                   r*r*G_i*G_j

   @return Sgb (top of p. 1694 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::Sgb( RealOpenMM t ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::Sgb";

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;

   // ---------------------------------------------------------------------------------------

   return ( (t != zero) ? one/SQRT( (one + (fourth*EXP( -t ))/t) ) : zero);
}

#define GBVIDebug 0

/**---------------------------------------------------------------------------------------

   Get GB/VI energy

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVI::computeBornEnergy( const vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                       const RealOpenMM* partialCharges ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "CpuGBVI::computeBornEnergy";

   static const RealOpenMM zero          = (RealOpenMM) 0.0;
   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM two           = (RealOpenMM) 2.0;
   static const RealOpenMM three         = (RealOpenMM) 3.0;
   static const RealOpenMM four          = (RealOpenMM) 4.0;
   static const RealOpenMM half          = (RealOpenMM) 0.5;
   static const RealOpenMM fourth        = (RealOpenMM) 0.25;
   static const RealOpenMM eighth        = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   const GBVIParameters* gbviParameters = getGBVIParameters();
   const RealOpenMM preFactor           = gbviParameters->getElectricConstant();
   const int numberOfAtoms              = gbviParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii        = gbviParameters->getAtomicRadii();
   const RealOpenMM* gammaParameters    = gbviParameters->getGammaParameters();

#if( GBVIDebug == 1 )
   FILE* logFile                        = stderr;
   (void) fprintf( logFile, "\n%s\n", methodName );
   (void) fflush( logFile );
#endif

   // ---------------------------------------------------------------------------------------

   // Eq.2 of Labute paper [JCC 29 p. 1693-1698 2008]
   // to minimze roundoff error sum cavityEnergy separately since in general much
   // smaller than other contributions

   RealOpenMM energy                 = zero;
   RealOpenMM cavityEnergy           = zero;

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      RealOpenMM partialChargeI   = partialCharges[atomI];

      // self-energy term

      RealOpenMM  atomIEnergy     = half*partialChargeI/bornRadii[atomI];

      // cavity term

      RealOpenMM ratio            = (atomicRadii[atomI]/bornRadii[atomI]);
      cavityEnergy               += gammaParameters[atomI]*ratio*ratio*ratio;

/*
RealOpenMM e1 = partialChargeI*partialCharges[atomI]/bornRadii[atomI];
RealOpenMM e2 = gammaParameters[atomI]*ratio*ratio*ratio;
(void) fprintf( stderr, "E %d self=%.4e gamma=%.4e e=%.4e\n", atomI, e1, e2, energy );
*/

      for( int atomJ = atomI + 1; atomJ < numberOfAtoms; atomJ++ ){

         RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
         if (_gbviParameters->getPeriodic())
             ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _gbviParameters->getPeriodicBox(), deltaR );
         else
             ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
         if (_gbviParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _gbviParameters->getCutoffDistance())
             continue;

         RealOpenMM r2                 = deltaR[ReferenceForce::R2Index];
         RealOpenMM t                  = fourth*r2/(bornRadii[atomI]*bornRadii[atomJ]);         
         atomIEnergy                  += partialCharges[atomJ]*Sgb( t )/deltaR[ReferenceForce::RIndex];

/*
RealOpenMM e3 =  -partialChargeI*partialCharges[atomJ]*Sgb( t )/deltaR[ReferenceForce::RIndex];
partialCharges[atomJ]*Sgb( t )/deltaR[ReferenceForce::RIndex];
(void) fprintf( stderr, "E %d %d e3=%.4e r2=%4e t=%.3e sgb=%.4e e=%.5e\n", atomI, atomJ, e3, r2, t, Sgb( t ), energy );
*/
      }

      energy += two*partialChargeI*atomIEnergy;
   }
   energy *= preFactor;
   energy -= cavityEnergy;

#if( GBVIDebug == 1 )
   (void) fprintf( logFile, "ElectricConstant=%.4e Tau=%.4e e=%.5e eOut=%.5e\n", preFactor, gbviParameters->getTau(), energy, gbviParameters->getTau()*energy );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "bR %d bR=%16.8e\n", atomI, bornRadii[atomI] );
   }
   (void) fflush( logFile );
#endif

   RealOpenMM conversion = (RealOpenMM)(gbviParameters->getTau());  
   return (conversion*energy);
 
}

#undef GBVIDebug

#define GBVIDebug 0

/**---------------------------------------------------------------------------------------

   Get GB/VI forces

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   --------------------------------------------------------------------------------------- */


void CpuGBVI::computeBornForces( const std::vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                const RealOpenMM* partialCharges, std::vector<OpenMM::RealVec>& inputForces){

   // ---------------------------------------------------------------------------------------

   static const char* methodName        = "CpuGBVI::computeBornEnergyForces";

   static const RealOpenMM zero         = (RealOpenMM) 0.0;
   static const RealOpenMM one          = (RealOpenMM) 1.0;
   static const RealOpenMM two          = (RealOpenMM) 2.0;
   static const RealOpenMM three        = (RealOpenMM) 3.0;
   static const RealOpenMM four         = (RealOpenMM) 4.0;
   static const RealOpenMM half         = (RealOpenMM) 0.5;
   static const RealOpenMM oneThird     = (RealOpenMM) (1.0/3.0);
   static const RealOpenMM fourth       = (RealOpenMM) 0.25;
   static const RealOpenMM eighth       = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

#if( GBVIDebug == 1 )
   FILE* logFile                        = stderr;
   (void) fprintf( logFile, "\n%s\n", methodName );
   (void) fflush( logFile );
#endif

   const GBVIParameters* gbviParameters = getGBVIParameters();
   const int numberOfAtoms              = gbviParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii        = gbviParameters->getAtomicRadii();
   const RealOpenMM* gammaParameters    = gbviParameters->getGammaParameters();

   // ---------------------------------------------------------------------------------------

   // constants

   const RealOpenMM preFactor           = two*gbviParameters->getElectricConstant();

   // ---------------------------------------------------------------------------------------

   // set energy/forces to zero

   RealOpenMM** forces  = new RealOpenMM*[numberOfAtoms];
   RealOpenMM*  block   = new RealOpenMM[numberOfAtoms*3];
	memset( block, 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
	RealOpenMM* blockPtr = block;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      forces[ii] = blockPtr;
		blockPtr  += 3;
   }

   vector<RealOpenMM>& bornForces = getBornForce();
   bornForces.assign(numberOfAtoms, 0.0);

   // ---------------------------------------------------------------------------------------

   // first main loop

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      // partial of polar term wrt Born radius
      // and (dGpol/dr)(dr/dx)

      RealOpenMM partialChargeI = preFactor*partialCharges[atomI];
      for( int atomJ = atomI; atomJ < numberOfAtoms; atomJ++ ){

         RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
         if (_gbviParameters->getPeriodic())
             ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _gbviParameters->getPeriodicBox(), deltaR );
         else
             ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
         if (_gbviParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _gbviParameters->getCutoffDistance())
             continue;

         RealOpenMM r2                 = deltaR[ReferenceForce::R2Index];

         RealOpenMM deltaX             = deltaR[ReferenceForce::XIndex];
         RealOpenMM deltaY             = deltaR[ReferenceForce::YIndex];
         RealOpenMM deltaZ             = deltaR[ReferenceForce::ZIndex];

         RealOpenMM alpha2_ij          = bornRadii[atomI]*bornRadii[atomJ];
         RealOpenMM D_ij               = r2/(four*alpha2_ij);

         RealOpenMM expTerm            = EXP( -D_ij );
         RealOpenMM denominator2       = r2 + alpha2_ij*expTerm; 
         RealOpenMM denominator        = SQRT( denominator2 ); 
         
         RealOpenMM Gpol               = (partialChargeI*partialCharges[atomJ])/denominator; 
         RealOpenMM dGpol_dr           = -Gpol*( one - fourth*expTerm )/denominator2;  

         RealOpenMM dGpol_dalpha2_ij   = -half*Gpol*expTerm*( one + D_ij )/denominator2;

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

         }

         // 3 FLOP

#if 0
if( atomI == 0 ){
   (void) fprintf( logFile, "bFCalc: %6d %6d %14.6e %14.6e   %14.6e %14.6e\n", atomI, atomJ, dGpol_dalpha2_ij, bornRadii[atomJ], bornForces[atomI], bornRadii[atomI] );
}
#endif
         bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

      }
   }

#if( GBVIDebug == 1 )
{
   double stupidFactor                   = three;
   RealOpenMM conversion                 = (RealOpenMM)(gbviParameters->getTau());  
   int maxPrint                          = 10;
   const RealOpenMM* scaledRadii         = gbviParameters->getScaledRadii();

   (void) fprintf( logFile, "Conversion=%14.6e %14.6e*%14.6e (tau)\n", conversion, 1, gbviParameters->getTau() );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      RealOpenMM R        = atomicRadii[atomI];

      // partial of cavity term wrt Born radius
     
      RealOpenMM  ratio   = (atomicRadii[atomI]/bornRadii[atomI]);
      bornForces[atomI]  += (stupidFactor*gammaParameters[atomI]*ratio*ratio*ratio)/bornRadii[atomI]; 

      RealOpenMM b2       = bornRadii[atomI]*bornRadii[atomI];
      double xx           = bornForces[atomI]*oneThird*b2*b2;

      // xx*conversion should agree w/ values pulled out of kReduceGBVIBornForces_kernel in kForces.cu

      (void) fprintf( logFile, "F1 %6d r/sclR[%14.6e %14.6e] bR=%14.6e bF=%14.6e %14.6e f[%14.6e %14.6e %14.6e](cnvrtd)"
                               " x[%14.6e %14.6e %14.6e]\n",
                      atomI, atomicRadii[atomI], scaledRadii[atomI], bornRadii[atomI], bornForces[atomI], xx*conversion,
//                      forces[atomI][0], forces[atomI][1], forces[atomI][2],
                      conversion*forces[atomI][0], conversion*forces[atomI][1],  conversion*forces[atomI][2],
                      atomCoordinates[atomI][0], atomCoordinates[atomI][1], atomCoordinates[atomI][2] );
      if( atomI == maxPrint ){
         atomI = numberOfAtoms - maxPrint;
         if( atomI < maxPrint )atomI = maxPrint;
      }
   }
   (void) fflush( logFile );
}
#endif

   // ---------------------------------------------------------------------------------------

   // second main loop: (dGpol/dBornRadius)(dBornRadius/dr)(dr/dx)

   // dGpol/dBornRadius) = bornForces[]
   // dBornRadius/dr     = (1/3)*(bR**4)*(dV/dr)

#if 0
   (void) fprintf( logFile, "Clearing forces before loop2 periodic=%d cutoff=%d cutoffR=%14.7e\n",
                   _gbviParameters->getPeriodic(), _gbviParameters->getUseCutoff(),  _gbviParameters->getCutoffDistance() );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      forces[atomI][0]  = zero;
      forces[atomI][1]  = zero;
      forces[atomI][2]  = zero;
   } 
   (void) fflush( logFile );
#endif

   const RealOpenMM* scaledRadii                 = gbviParameters->getScaledRadii();
   std::vector<RealOpenMM>& switchDeriviative    = getSwitchDeriviative();
   RealOpenMM stupidFactor                       = three;
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
      RealOpenMM R        = atomicRadii[atomI];

      // partial of cavity term wrt Born radius
     
      RealOpenMM  ratio   = (atomicRadii[atomI]/bornRadii[atomI]);
      bornForces[atomI]  += (stupidFactor*gammaParameters[atomI]*ratio*ratio*ratio)/bornRadii[atomI]; 

      RealOpenMM b2       = bornRadii[atomI]*bornRadii[atomI];
      bornForces[atomI]  *= switchDeriviative[atomI]*oneThird*b2*b2;

      for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

         if( atomJ != atomI ){

            RealOpenMM deltaX             = atomCoordinates[atomJ][0] - atomCoordinates[atomI][0];
            RealOpenMM deltaY             = atomCoordinates[atomJ][1] - atomCoordinates[atomI][1];
            RealOpenMM deltaZ             = atomCoordinates[atomJ][2] - atomCoordinates[atomI][2];
    
            RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
            if (_gbviParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _gbviParameters->getPeriodicBox(), deltaR );
            else
                ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
            if (_gbviParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > _gbviParameters->getCutoffDistance())
                continue;
   
            RealOpenMM r2                 = deltaR[ReferenceForce::R2Index];
                       deltaX             = deltaR[ReferenceForce::XIndex];
                       deltaY             = deltaR[ReferenceForce::YIndex];
                       deltaZ             = deltaR[ReferenceForce::ZIndex];

            RealOpenMM r                  = SQRT( r2 );
 
            RealOpenMM S                  = scaledRadii[atomJ];
            RealOpenMM diff               = (S - R);

            RealOpenMM de                 = zero;

            // find dRb/dr, where Rb is the Born radius

            de = CpuGBVI::dL_dr( r, r+S, S ) + CpuGBVI::dL_dx( r, r+S, S );   
            if( FABS( diff ) < r ){
               if( R > (r - S) ){
                  de -= CpuGBVI::dL_dr( r, R, S );  
               } else {
                  de -= ( CpuGBVI::dL_dr( r, (r-S), S ) + CpuGBVI::dL_dx( r, (r-S), S ) );
               }
            } else if( r < (S - R) ){
               de -= ( CpuGBVI::dL_dr( r, r-S, S ) + CpuGBVI::dL_dx( r, r-S, S ) );   
            }

#if 0
   RealOpenMM delta = (RealOpenMM) 1.0e-02;
   (void) fprintf( stderr, "\n" );
   for( int kk = 0; kk < 5; kk++ ){
      RealOpenMM V1    = CpuGBVI::getVolume( r, R, S );
      RealOpenMM V2    = CpuGBVI::getVolume( r+delta, R, S );
      RealOpenMM df    = (V2-V1)/delta;
      (void) fprintf( stderr, "df %d %d [%14.6e %14.6e] V[%14.6e %14.6e] %.2e\n", atomI, atomJ, de, df, V2, V1, delta );
      delta *= (RealOpenMM) 0.1;
   }

   double deltaD = 1.0e-02;
   double ded = CpuGBVI::dL_drD( (double) r, r+S, S ) + CpuGBVI::dL_dxD( r, r+S, S ) - ( CpuGBVI::dL_drD( r, (r-S), S ) + CpuGBVI::dL_dxD( r, (r-S), S ) );
   for( int kk = 0; kk < 5; kk++ ){
      double V1    = CpuGBVI::getVolumeD( r, R, S );
      double V2    = CpuGBVI::getVolumeD( r+deltaD, R, S );
      double df    = (V2-V1)/deltaD;
      (void) fprintf( stderr, "df %d %d [%14.6e %14.6e] V[%14.6e %14.6e] %.2e\n", atomI, atomJ, ded, df, V2, V1, deltaD );
      deltaD *= 0.1;
   }
#endif


             // de = (dG/dRb)(dRb/dr)

            de                      *= bornForces[atomI]/r;

            deltaX                  *= de;
            deltaY                  *= de;
            deltaZ                  *= de;
   
            forces[atomI][0]        += deltaX;
            forces[atomI][1]        += deltaY;
            forces[atomI][2]        += deltaZ;
  
            forces[atomJ][0]        -= deltaX;
            forces[atomJ][1]        -= deltaY;
            forces[atomJ][2]        -= deltaZ;

         }
      }

   }

#if( GBVIDebug == 1 )
   (void) fprintf( logFile, "Atom      BornRadii      BornForce                                         Forces\n" );
   double forceSum[3] = { 0.0, 0.0, 0.0 };
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      forceSum[0] += forces[atomI][0];
      forceSum[1] += forces[atomI][1];
      forceSum[2] += forces[atomI][2];
      (void) fprintf( logFile, "%4d %14.6e %14.6e [%14.6e %14.6e %14.6e]\n", atomI, bornRadii[atomI], bornForces[atomI],  forces[atomI][0], forces[atomI][1], forces[atomI][2] );
   }   
   (void) fprintf( logFile, "F sum=[%14.6e %14.6e %14.6e]\n", forceSum[0], forceSum[1], forceSum[2] );
   (void) fflush( logFile );
#endif

   // convert from cal to Joule & apply prefactor tau = (1/diel_solute - 1/diel_solvent)

   RealOpenMM conversion = (RealOpenMM)(gbviParameters->getTau());  
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      inputForces[atomI][0] += conversion*forces[atomI][0];
      inputForces[atomI][1] += conversion*forces[atomI][1];
      inputForces[atomI][2] += conversion*forces[atomI][2];
   }

#if( GBVIDebug == 1 )
{
   (void) fprintf( logFile, "\nPost conversion\n" );
   (void) fprintf( logFile, "Atom      BornRadii      BornForce                                         Forces\n" );
   int maxPrint = 10;
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "%4d %14.6e %14.6e [%14.6e %14.6e %14.6e]\n", atomI, bornRadii[atomI], conversion*bornForces[atomI],
                      inputForces[atomI][0], inputForces[atomI][1], inputForces[atomI][2] );
      if( atomI == maxPrint ){
         atomI = numberOfAtoms - maxPrint;
         if( atomI < maxPrint )atomI = numberOfAtoms;
      }
   }
   (void) fflush( logFile );
}
#endif
#undef GBVIDebug

   delete[] forces;
   delete[] block;

}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      
   @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuGBVI::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << CpuImplicitSolvent::getStateString( title );

   return message.str();
}

/**---------------------------------------------------------------------------------------

   Write Born energy and forces (Simbios)

   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces
   @param resultsFileName     output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuGBVI::writeBornEnergyForces( vector<RealVec>& atomCoordinates,
                                   const RealOpenMM* partialCharges, vector<RealVec>& forces,
                                   const std::string& resultsFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nCpuGBVI::writeBornEnergyForces";

   // ---------------------------------------------------------------------------------------
/*
   ImplicitSolventParameters* implicitSolventParameters = getImplicitSolventParameters();
   const GBVIParameters* gbviParameters                   = static_cast<const GBVIParameters*>(implicitSolventParameters);
   

   int numberOfAtoms                    = gbviParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii        = gbviParameters->getAtomicRadii();
   const RealOpenMM* bornRadii          = getBornRadiiConst();
   const RealOpenMM* scaledRadii        = gbviParameters->getScaledRadiusFactors();
   const RealOpenMM* gbviChain           = getObcChainConst();
   const RealOpenMM  energy             = getEnergy();

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* implicitSolventResultsFile = NULL;
#ifdef WIN32
   fopen_s( &implicitSolventResultsFile, resultsFileName.c_str(), "w" );
#else
   implicitSolventResultsFile = fopen( resultsFileName.c_str(), "w" );
#endif

   // diganostics

   std::stringstream message;
   message << methodName;
   if( implicitSolventResultsFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << resultsFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << resultsFileName << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // header

   (void) fprintf( implicitSolventResultsFile, "# %d atoms E=%.7e   format: coords(3) bornRadii(input) q atomicRadii scaleFactors forces gbviChain\n",
                   numberOfAtoms, energy );

   RealOpenMM forceConversion  = (RealOpenMM) 1.0;
   RealOpenMM lengthConversion = (RealOpenMM) 1.0;

   // output

   if( forces != NULL && atomCoordinates != NULL && partialCharges != NULL && atomicRadii != NULL ){
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
            (void) fprintf( implicitSolventResultsFile, "%.7e %.7e %.7e %.7e %.5f %.5f %.5f %.7e %.7e %.7e %.7e\n",
                            lengthConversion*atomCoordinates[ii][0],
                            lengthConversion*atomCoordinates[ii][1], 
                            lengthConversion*atomCoordinates[ii][2],
                           (bornRadii != NULL ? lengthConversion*bornRadii[ii] : 0.0),
                            partialCharges[ii], lengthConversion*atomicRadii[ii], scaledRadii[ii],
                            forceConversion*forces[ii][0],
                            forceConversion*forces[ii][1],
                            forceConversion*forces[ii][2],
                            forceConversion*gbviChain[ii]
                          );
      }
   }
   (void) fclose( implicitSolventResultsFile );

*/
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get Obc Born energy and forces -- used debugging

   @param bornRadii           Born radii -- optional; if NULL, then GBVIParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   The array bornRadii is also updated and the obcEnergy

   --------------------------------------------------------------------------------------- */

void CpuGBVI::computeBornEnergyForces( RealOpenMM* bornRadii, vector<RealVec>& atomCoordinates,
                                      const RealOpenMM* partialCharges, vector<RealVec>& forces ){
 
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVI::computeBornEnergyForcesPrint";

}

/**---------------------------------------------------------------------------------------

   Use double precision 

   Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return volume

   --------------------------------------------------------------------------------------- */

double CpuGBVI::getVolumeD( double r, double R, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::getVolume";

   static const double zero         = 0.0;
   static const double minusThree   = -3.0;

   double              diff    = (S - R);
   if( fabs( diff ) < r ){

      double lowerBound = (R > (r - S)) ? R : (r - S);

      return (CpuGBVI::getLD( r, (r + S),    S ) -
              CpuGBVI::getLD( r, lowerBound, S ));

   } else if( r < diff ){

      return CpuGBVI::getLD( r, (r + S), S ) -
             CpuGBVI::getLD( r, (r - S), S ) + 
             pow( R, minusThree );

   } else {
      return zero;
   }
}

/**---------------------------------------------------------------------------------------

   Use double precision 

   Get L (used in analytical solution for volume integrals) 

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return L value (Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

double CpuGBVI::getLD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::getL";

   static const double one           =  1.0;
   static const double threeHalves   =  1.5;
   static const double third         =  1.0/3.0;
   static const double fourth        =  0.25;
   static const double eighth        =  0.125;

   // ---------------------------------------------------------------------------------------

   double rInv   = one/r;

   double xInv   = one/x;
   double xInv2  = xInv*xInv;
   double xInv3  = xInv2*xInv;

   double diff2  = (r + S)*(r - S);

   return (threeHalves*xInv)*( (xInv*fourth*rInv) - (xInv2*third) + (diff2*xInv3*eighth*rInv) );
}

/**---------------------------------------------------------------------------------------

   Use double precision 

   Get partial derivative of L wrt r

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

double CpuGBVI::dL_drD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::dL_dr";

   static const double one           =  1.0;
   static const double threeHalves   =  1.5;
   static const double threeEights   =  0.375;
   static const double third         =  1.0/3.0;
   static const double fourth        =  0.25;
   static const double eighth        =  0.125;

   // ---------------------------------------------------------------------------------------

   double rInv   = one/r;
   double rInv2  = rInv*rInv;

   double xInv   = one/x;
   double xInv2  = xInv*xInv;
   double xInv3  = xInv2*xInv;

   double diff2  = (r + S)*(r - S);

   return ( (-threeHalves*xInv2*rInv2)*( fourth + eighth*diff2*xInv2 ) + threeEights*xInv3*xInv );
}

/**---------------------------------------------------------------------------------------

   Use double precision 

   Get partial derivative of L wrt x

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

double CpuGBVI::dL_dxD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVI::dL_dx";

   static const double one           =   1.0;
   static const double half          =   0.5;
   static const double threeHalvesM  =  -1.5;
   static const double third         =   1.0/3.0;

   // ---------------------------------------------------------------------------------------

   double rInv   = one/r;

   double xInv   = one/x;
   double xInv2  = xInv*xInv;
   double xInv3  = xInv2*xInv;

   double diff   = (r + S)*(r - S);

   return (threeHalvesM*xInv3)*( (half*rInv) - xInv + (half*diff*xInv2*rInv) );
}
