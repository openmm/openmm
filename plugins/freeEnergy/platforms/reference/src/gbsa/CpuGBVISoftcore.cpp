
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
#include "CpuGBVISoftcore.h"
#include "../SimTKReference/ReferenceForce.h"
#include <math.h>

/**---------------------------------------------------------------------------------------

   CpuGBVISoftcore constructor

   gbviParameters      gbviParameters object
   
   --------------------------------------------------------------------------------------- */

CpuGBVISoftcore::CpuGBVISoftcore( ImplicitSolventParameters* gbviParameters ) : CpuImplicitSolvent( gbviParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::CpuGBVISoftcore";

   // ---------------------------------------------------------------------------------------

   _initializeGBVISoftcoreDataMembers( );

   _gbviParameters = static_cast<GBVISoftcoreParameters*> (gbviParameters);

}

/**---------------------------------------------------------------------------------------

   CpuGBVISoftcore destructor

   --------------------------------------------------------------------------------------- */

CpuGBVISoftcore::~CpuGBVISoftcore( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::~CpuGBVISoftcore";

   // ---------------------------------------------------------------------------------------

   delete _switchDeriviative;
   //if( _gbviParameters != NULL ){
     // delete _gbviParameters;
   //}

}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::_initializeGBVISoftcoreDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _gbviParameters    = NULL;
   _switchDeriviative = NULL;
}

/**---------------------------------------------------------------------------------------

   Get GBVISoftcoreParameters reference

   @return GBVISoftcoreParameters reference

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters* CpuGBVISoftcore::getGBVISoftcoreParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::getGBVISoftcoreParameters";

   // ---------------------------------------------------------------------------------------

   return _gbviParameters;
}

/**---------------------------------------------------------------------------------------

   Set GBVISoftcoreParameters reference

   @param GBVISoftcoreParameters reference

   @return 0;

   --------------------------------------------------------------------------------------- */

int CpuGBVISoftcore::setGBVISoftcoreParameters( GBVISoftcoreParameters* gbviParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::setGBVISoftcoreParameters";

   // ---------------------------------------------------------------------------------------

   _gbviParameters = gbviParameters;
   return 0;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcSoftcoreParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMM* CpuGBVISoftcore::getSwitchDeriviative( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::getSwitchDeriviative";

   // ---------------------------------------------------------------------------------------

   if( _switchDeriviative == NULL && _gbviParameters != NULL ){
      _switchDeriviative = new RealOpenMM[_gbviParameters->getNumberOfAtoms()];
   }
   return _switchDeriviative;
}

/**---------------------------------------------------------------------------------------

   Return switching function derivative

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMM* CpuGBVISoftcore::getSwitchDeriviativeConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::getSwitchDeriviative";

   // ---------------------------------------------------------------------------------------

   return _switchDeriviative;
}

/**---------------------------------------------------------------------------------------

   Compute quintic spline value and associated derviative

   @param x                   value to compute spline at
   @param rl                  lower cutoff value
   @param ru                  upper cutoff value
   @param outValue            value of spline at x
   @param outDerivative       value of derivative of spline at x

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

#define GBVISoftcoreDebug 0

int CpuGBVISoftcore::quinticSpline( RealOpenMM x, RealOpenMM rl, RealOpenMM ru,
                                    RealOpenMM* outValue, RealOpenMM* outDerivative ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM one           = (RealOpenMM)   1.0;
   static const RealOpenMM minusSix      = (RealOpenMM)  -6.0;
   static const RealOpenMM minusTen      = (RealOpenMM) -10.0;
   static const RealOpenMM minusThirty   = (RealOpenMM) -30.0;
   static const RealOpenMM fifteen       = (RealOpenMM)  15.0;
   static const RealOpenMM sixty         = (RealOpenMM)  60.0;

   // static const char* methodName         = "CpuGBVISoftcore::quinticSpline";

   // ---------------------------------------------------------------------------------------

   RealOpenMM numerator    = x  - rl;
   RealOpenMM denominator  = ru - rl;
   RealOpenMM ratio        = numerator/denominator;
   RealOpenMM ratio2       = ratio*ratio;
   RealOpenMM ratio3       = ratio2*ratio;

   *outValue               = one + ratio3*(minusTen + fifteen*ratio + minusSix*ratio2);
   *outDerivative          = ratio2*(minusThirty + sixty*ratio + minusThirty*ratio2)/denominator;

   return 0;
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

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

#define GBVISoftcoreDebug 0

int CpuGBVISoftcore::computeBornRadiiUsingQuinticSpline( RealOpenMM atomicRadius3, RealOpenMM bornSum,
                                                         GBVISoftcoreParameters* gbviParameters, 
                                                         RealOpenMM* bornRadius, RealOpenMM* switchDeriviative ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero          = (RealOpenMM)  0.0;
   static const RealOpenMM one           = (RealOpenMM)  1.0;
   static const RealOpenMM minusOne      = (RealOpenMM) -1.0;
   static const RealOpenMM minusThree    = (RealOpenMM) -3.0;
   static const RealOpenMM oneEighth     = (RealOpenMM)  0.125;
   static const RealOpenMM minusOneThird = (RealOpenMM) (-1.0/3.0);
   static const RealOpenMM three         = (RealOpenMM)  3.0;

   static const char* methodName         = "CpuGBVISoftcore::computeBornRadiiUsingQuinticSpline";

   // ---------------------------------------------------------------------------------------

#if( GBVISoftcoreDebug == 1 )
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
         sum                 = (atomicRadius3 - bornSum)*splineValue + gbviParameters->getQuinticUpperSplineLimit();
         *switchDeriviative  = splineValue - (atomicRadius3 - bornSum)*splineDerivative;
#if( GBVISoftcoreDebug == 1 )
      (void) fprintf( logFile, " Qv=%14.6e splnDrvtv=%14.6e spline[%10.3e %10.3e] ", splineValue, splineDerivative,
                      splineL, gbviParameters->getQuinticUpperSplineLimit() );
#endif
      } else {   
         sum                = gbviParameters->getQuinticUpperSplineLimit();
         *switchDeriviative = zero;
      }
   } else {
      sum                = atomicRadius3 - bornSum; 
      *switchDeriviative = one;
   }
   *bornRadius = POW( sum, minusOneThird );
  
   return 0; 
}

#undef GBVISoftcoreDebug

/**---------------------------------------------------------------------------------------

   Get Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii
   @param chain               not used here

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

#define GBVISoftcoreDebug 0

int CpuGBVISoftcore::computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii, RealOpenMM* switchDeriviative ){

   // ---------------------------------------------------------------------------------------

   static const RealOpenMM zero          = (RealOpenMM) 0.0;
   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM minusThree    = (RealOpenMM) -3.0;
   static const RealOpenMM oneEighth     = (RealOpenMM) 0.125;
   static const RealOpenMM minusOneThird = (RealOpenMM) (-1.0/3.0);
   static const RealOpenMM three         = (RealOpenMM) 3.0;

   static const char* methodName         = "CpuGBVISoftcore::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   GBVISoftcoreParameters* gbviParameters           = getGBVISoftcoreParameters();
   int numberOfAtoms                                = gbviParameters->getNumberOfAtoms();
   RealOpenMM* atomicRadii                          = gbviParameters->getAtomicRadii();
   const RealOpenMM* scaledRadii                    = gbviParameters->getScaledRadii();
   const RealOpenMM* bornRadiusScaleFactors         = gbviParameters->getBornRadiusScaleFactors();

   if( switchDeriviative == NULL ){
      switchDeriviative = getSwitchDeriviative();
   }

   // ---------------------------------------------------------------------------------------

#if( GBVISoftcoreDebug == 1 )
   FILE* logFile                         = stderr;
   (void) fprintf( logFile, "\n%s\n", methodName );
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

            sum  += bornRadiusScaleFactors[atomJ]*CpuGBVISoftcore::getVolume( r, radiusI, scaledRadii[atomJ] );

#if( GBVISoftcoreDebug == -1 )
if( atomI == 0 || atomI == 1 ){
            (void) fprintf( logFile, "%d addJ=%d scR=%14.6e %14.6e sum=%14.6e rI=%14.6e r=%14.6e S-R=%14.6e\n",
                            atomI, atomJ, scaledRadii[atomJ], getVolume( r, radiusI, scaledRadii[atomJ] ), sum,
                            radiusI, r, (scaledRadii[atomJ]-radiusI) );
}
#endif
         }
      }

#if( GBVISoftcoreDebug == 1 )
      (void) fprintf( logFile, "%6d BornSum=%14.6e r=%14.6e r3=%14.6e (r3-sum)=%14.6e method=%d ",
                       atomI, sum, radiusI, POW( radiusI, minusThree ), (POW( radiusI, minusThree ) - sum), _gbviParameters->getBornRadiusScalingSoftcoreMethod() );
#endif

      RealOpenMM atomicRadius3 = POW( radiusI, minusThree );
      if( _gbviParameters->getBornRadiusScalingSoftcoreMethod() == GBVISoftcoreParameters::NoScaling ){
         sum                       = atomicRadius3 - sum;
         bornRadii[atomI]          = POW( sum, minusOneThird );
         switchDeriviative[atomI]  = one;
      } else if( _gbviParameters->getBornRadiusScalingSoftcoreMethod() == GBVISoftcoreParameters::QuinticSpline ){
         computeBornRadiiUsingQuinticSpline( atomicRadius3, sum, gbviParameters, 
                                             bornRadii + atomI, switchDeriviative + atomI );
      }

#if( GBVISoftcoreDebug == 1 )
      (void) fprintf( logFile, "br=%14.6e swDrvtv=%14.6e %s\n", bornRadii[atomI], switchDeriviative[atomI],  (fabs( switchDeriviative[atomI] - 1.0 ) > 1.0e-05 ? "SWWWWW" : "")  );
#endif

   }

   return 0;

}

#undef GBVISoftcoreDebug

/**---------------------------------------------------------------------------------------

   Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return volume

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVISoftcore::getVolume( RealOpenMM r, RealOpenMM R, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::getVolume";

   static const RealOpenMM zero         = (RealOpenMM)  0.0;
   static const RealOpenMM minusThree   = (RealOpenMM) -3.0;

   RealOpenMM              diff         = (S - R);
   if( FABS( diff ) < r ){

      RealOpenMM lowerBound = (R > (r - S)) ? R : (r - S);

      return (CpuGBVISoftcore::getL( r, (r + S),    S ) -
              CpuGBVISoftcore::getL( r, lowerBound, S ));

   } else if( r <= diff ){

      return CpuGBVISoftcore::getL( r, (r + S), S ) -
             CpuGBVISoftcore::getL( r, (r - S), S ) + 
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

RealOpenMM CpuGBVISoftcore::getL( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::getL";

   static const RealOpenMM one           = (RealOpenMM) 1.0;
   static const RealOpenMM threeHalves   = (RealOpenMM) 1.5;
   static const RealOpenMM third         = (RealOpenMM) (1.0/3.0);
   static const RealOpenMM fourth        = (RealOpenMM) 0.25;
   static const RealOpenMM eighth        = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   RealOpenMM rInv   = one/r;

   RealOpenMM xInv   = one/x;
   RealOpenMM xInv2  = xInv*xInv;

   RealOpenMM diff2  = (r + S)*(r - S);

   return (threeHalves*xInv2)*( (fourth*rInv) - (xInv*third) + (diff2*xInv2*eighth*rInv) );
}

/**---------------------------------------------------------------------------------------

   Get partial derivative of L wrt r

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return partial derivative based on Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVISoftcore::dL_dr( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuGBVISoftcore::dL_dr";

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

RealOpenMM CpuGBVISoftcore::dL_dx( RealOpenMM r, RealOpenMM x, RealOpenMM S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::dL_dx";

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

RealOpenMM CpuGBVISoftcore::Sgb( RealOpenMM t ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::Sgb";

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;

   // ---------------------------------------------------------------------------------------

   return ( (t != zero) ? one/SQRT( (one + (fourth*EXP( -t ))/t) ) : zero);
}

#define GBVISoftcoreDebug 0

/**---------------------------------------------------------------------------------------

   Get GB/VI energy

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVISoftcore::computeBornEnergy( const RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                               const RealOpenMM* partialCharges ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "CpuGBVISoftcore::computeBornEnergy";

   static const RealOpenMM zero    = (RealOpenMM) 0.0;
   static const RealOpenMM one     = (RealOpenMM) 1.0;
   static const RealOpenMM two     = (RealOpenMM) 2.0;
   static const RealOpenMM three   = (RealOpenMM) 3.0;
   static const RealOpenMM four    = (RealOpenMM) 4.0;
   static const RealOpenMM half    = (RealOpenMM) 0.5;
   static const RealOpenMM fourth  = (RealOpenMM) 0.25;
   static const RealOpenMM eighth  = (RealOpenMM) 0.125;

   // ---------------------------------------------------------------------------------------

   const GBVISoftcoreParameters* gbviParameters = getGBVISoftcoreParameters();
   const RealOpenMM preFactor           = gbviParameters->getElectricConstant();
   const int numberOfAtoms              = gbviParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii        = gbviParameters->getAtomicRadii();
   const RealOpenMM* gammaParameters    = gbviParameters->getGammaParameters();

   if( bornRadii == NULL ){
      bornRadii   = getBornRadii();
   }

#if( GBVISoftcoreDebug == 1 )
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
RealOpenMM e3 = -partialChargeI2*partialCharges[atomJ]*Sgb( t )/deltaR[ReferenceForce::RIndex];
(void) fprintf( stderr, "E %d %d e3=%.4e r2=%4e t=%.3e sgb=%.4e e=%.5e\n", atomI, atomJ, e3, r2, t, Sgb( t ), energy );
*/
      }

      energy += two*partialChargeI*atomIEnergy;
   }
   energy *= preFactor;
   energy -= cavityEnergy;

#if( GBVISoftcoreDebug == 1 )
   (void) fprintf( logFile, "ElectricConstant=%.4e Tau=%.4e e=%.5e eOut=%.5e\n", preFactor, gbviParameters->getTau(), energy, gbviParameters->getTau()*energy );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "bR %d bR=%16.8e\n", atomI, bornRadii[atomI] );
   }
   (void) fflush( logFile );
#endif

   RealOpenMM conversion = (RealOpenMM)(gbviParameters->getTau());  
   return (conversion*energy);
 
}

#undef GBVISoftcoreDebug

#define GBVISoftcoreDebug 0

/**---------------------------------------------------------------------------------------

   Get GB/VI forces

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return 0;

   --------------------------------------------------------------------------------------- */


int CpuGBVISoftcore::computeBornForces( const RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                        const RealOpenMM* partialCharges, RealOpenMM** inputForces ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName        = "CpuGBVISoftcore::computeBornForces";

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

#if( GBVISoftcoreDebug == 1 || GBVISoftcoreDebug == 2 )
   FILE* logFile                        = stderr;
   (void) fprintf( logFile, "\n%s\n", methodName );
   (void) fflush( logFile );
#endif

   const GBVISoftcoreParameters* gbviParameters = getGBVISoftcoreParameters();
   const int numberOfAtoms              = gbviParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii        = gbviParameters->getAtomicRadii();
   const RealOpenMM* gammaParameters    = gbviParameters->getGammaParameters();

   if( bornRadii == NULL ){
      bornRadii   = getBornRadii();
   }

   // ---------------------------------------------------------------------------------------

   // constants

   const RealOpenMM preFactor           = two*gbviParameters->getElectricConstant();

   // ---------------------------------------------------------------------------------------

   // set energy/forces to zero

   const unsigned int arraySzInBytes    = sizeof( RealOpenMM )*numberOfAtoms;

   RealOpenMM** forces  = new RealOpenMM*[numberOfAtoms];
   RealOpenMM*  block   = new RealOpenMM[numberOfAtoms*3];
	memset( block, 0, sizeof( RealOpenMM )*numberOfAtoms*3 );
	RealOpenMM* blockPtr = block;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      forces[ii] = blockPtr;
		blockPtr  += 3;
   }

   RealOpenMM* bornForces = getBornForce();
   memset( bornForces, 0, arraySzInBytes );

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

#if( GBVISoftcoreDebug == 1 )
{
   double stupidFactor                   = three;
   RealOpenMM conversion                 = (RealOpenMM)(gbviParameters->getTau());  
   int maxPrint                          = 20;
   const RealOpenMM* scaledRadii         = gbviParameters->getScaledRadii();
   RealOpenMM* switchDeriviative         = getSwitchDeriviative();

   (void) fprintf( logFile, "F1: Conversion=%14.6e %14.6e*%14.6e (tau)\n", conversion, 1, gbviParameters->getTau() );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      RealOpenMM R        = atomicRadii[atomI];
      RealOpenMM  ratio   = (atomicRadii[atomI]/bornRadii[atomI]);
      RealOpenMM  bF      = bornForces[atomI]  + (stupidFactor*gammaParameters[atomI]*ratio*ratio*ratio)/bornRadii[atomI]; 

      RealOpenMM b2       = bornRadii[atomI]*bornRadii[atomI];
      double xx           = switchDeriviative[atomI]*bF*oneThird*b2*b2;

      // xx*conversion should agree w/ values pulled out of kReduceGBVISoftcoreBornForces_kernel in kForces.cu
      (void) fprintf( logFile, "F1 %6d r/sclR[%14.6e %14.6e] bR=%14.6e bF=%14.6e sw=%14.6e f[%14.6e %14.6e %14.6e](cnvrtd)"
                               " x[%14.6e %14.6e %14.6e]\n",
                      atomI, atomicRadii[atomI], scaledRadii[atomI], bornRadii[atomI], xx*conversion, switchDeriviative[atomI],
                      conversion*forces[atomI][0], conversion*forces[atomI][1],  conversion*forces[atomI][2],
                      atomCoordinates[atomI][0], atomCoordinates[atomI][1], atomCoordinates[atomI][2] );
      if( atomI == maxPrint ){
         atomI = numberOfAtoms - maxPrint;
         if( atomI < maxPrint )atomI = maxPrint;
      }
   }
   (void) fflush( logFile );
   int clearForces = 0;
   if( clearForces ){
      (void) fprintf( logFile, "Forces cleared after loop 1\n" );
      for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
          forces[atomI][0] = 0.0f;
          forces[atomI][1] = 0.0f;
          forces[atomI][2] = 0.0f;
      }
   }
}
#endif

#if( GBVISoftcoreDebug == 2 )
{
   double stupidFactor                   = three;
   RealOpenMM conversion                 = (RealOpenMM)(gbviParameters->getTau());  
   int maxPrint                          = 1000000;
   const RealOpenMM* scaledRadii         = gbviParameters->getScaledRadii();
   RealOpenMM* switchDeriviative         = getSwitchDeriviative();

   (void) fprintf( logFile, "F1: Conversion=%14.6e %14.6e*%14.6e (tau)\n", conversion, 1, gbviParameters->getTau() );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      RealOpenMM R        = atomicRadii[atomI];
      RealOpenMM  ratio   = (atomicRadii[atomI]/bornRadii[atomI]);
      RealOpenMM  bF      = bornForces[atomI]  + (stupidFactor*gammaParameters[atomI]*ratio*ratio*ratio)/bornRadii[atomI]; 

      RealOpenMM b2       = bornRadii[atomI]*bornRadii[atomI];
      double xx           = switchDeriviative[atomI]*bF*oneThird*b2*b2;

      // xx*conversion should agree w/ values pulled out of kReduceGBVISoftcoreBornForces_kernel in kForces.cu
/*
      (void) fprintf( logFile, "F1 %6d r/sclR[%14.6e %14.6e] bR=%14.6e sw=%14.6e bF=%14.6e %14.6e f[%14.6e %14.6e %14.6e](cnvrtd)"
                               " x[%14.6e %14.6e %14.6e]\n",
                      atomI, atomicRadii[atomI], scaledRadii[atomI], bornRadii[atomI], bF, switchDeriviative[atomI], xx*conversion,
                      conversion*forces[atomI][0], conversion*forces[atomI][1],  conversion*forces[atomI][2],
                      atomCoordinates[atomI][0], atomCoordinates[atomI][1], atomCoordinates[atomI][2] );
*/
      (void) fprintf( logFile, "%6d %14.6e %14.6e %14.6e %14.6e %14.6e    %14.6e %14.6e %14.6e    %14.6e %14.6e %14.6e   %14.6e\n",
                      atomI, atomicRadii[atomI], scaledRadii[atomI], bornRadii[atomI], bF, xx*conversion,
                      conversion*forces[atomI][0], conversion*forces[atomI][1],  conversion*forces[atomI][2],
                      atomCoordinates[atomI][0], atomCoordinates[atomI][1], atomCoordinates[atomI][2], switchDeriviative[atomI] );
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

   const RealOpenMM* scaledRadii                    = gbviParameters->getScaledRadii();
   RealOpenMM* switchDeriviative                    = getSwitchDeriviative();
   RealOpenMM stupidFactor                          = three;
   const RealOpenMM* bornRadiusScaleFactors         = gbviParameters->getBornRadiusScaleFactors();
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

            // find dRb/dr, where Rb is the Born radius

            RealOpenMM de = CpuGBVISoftcore::dL_dr( r, r+S, S ) + CpuGBVISoftcore::dL_dx( r, r+S, S );   
            if( FABS( diff ) < r ){
               if( R > (r - S) ){
                  de -= CpuGBVISoftcore::dL_dr( r, R, S );  
               } else {
                  de -= ( CpuGBVISoftcore::dL_dr( r, (r-S), S ) + CpuGBVISoftcore::dL_dx( r, (r-S), S ) );
               }
            } else if( r < (S - R) ){
               de -= ( CpuGBVISoftcore::dL_dr( r, r-S, S ) + CpuGBVISoftcore::dL_dx( r, r-S, S ) );   
            }


#if 0
   for( int kk = 0; kk < 5; kk++ ){
      RealOpenMM V1    = CpuGBVISoftcore::getVolume( r, R, S );
      RealOpenMM V2    = CpuGBVISoftcore::getVolume( r+delta, R, S );
      RealOpenMM df    = (V2-V1)/delta;
      (void) fprintf( stderr, "df %d %d [%14.6e %14.6e] V[%14.6e %14.6e] %.2e\n", atomI, atomJ, de, df, V2, V1, delta );
      delta *= (RealOpenMM) 0.1;
   }

   double deltaD = 1.0e-02;
   double ded = CpuGBVISoftcore::dL_drD( (double) r, r+S, S ) + CpuGBVISoftcore::dL_dxD( r, r+S, S ) - ( CpuGBVISoftcore::dL_drD( r, (r-S), S ) + CpuGBVISoftcore::dL_dxD( r, (r-S), S ) );
   for( int kk = 0; kk < 5; kk++ ){
      double V1    = CpuGBVISoftcore::getVolumeD( r, R, S );
      double V2    = CpuGBVISoftcore::getVolumeD( r+deltaD, R, S );
      double df    = (V2-V1)/deltaD;
      (void) fprintf( stderr, "df %d %d [%14.6e %14.6e] V[%14.6e %14.6e] %.2e\n", atomI, atomJ, ded, df, V2, V1, deltaD );
      deltaD *= 0.1;
   }
#endif

             // de = (dG/dRb)(dRb/dr)

            de                      *= bornRadiusScaleFactors[atomJ]*bornForces[atomI]/r;

            deltaX                  *= de;
            deltaY                  *= de;
            deltaZ                  *= de;
   
            forces[atomI][0]        += deltaX;
            forces[atomI][1]        += deltaY;
            forces[atomI][2]        += deltaZ;
  
            forces[atomJ][0]        -= deltaX;
            forces[atomJ][1]        -= deltaY;
            forces[atomJ][2]        -= deltaZ;

#if 0
   if( atomI == 2613 ){
   (void) fprintf( stderr, "AtomJ %5d r=%14.7e de=%14.7e bfI=%14.7e finalDe=%14.7e [%14.7e %14.7e %14.7e]\n",
                   atomJ, r, de, bornForces[atomI], (de*bornForces[atomI]/r), 
                   forces[atomI][0], forces[atomI][1], forces[atomI][2]  );
   } else if( atomJ == 2613 ){
   (void) fprintf( stderr, "AtomI %5d r=%14.7e de=%14.7e bfI=%14.7e finalDe=%14.7e [%14.7e %14.7e %14.7e]\n",
                   atomI, r, de, bornForces[atomI], (de*bornForces[atomI]/r), 
                   forces[atomJ][0], forces[atomJ][1], forces[atomJ][2]  );
   }
#endif
         }
      }

   }

#if( GBVISoftcoreDebug == 9 )
{
   (void) fprintf( logFile, "\nPre conversion\n" );
   (void) fprintf( logFile, "Atom        ScaledRadii    BornRadii      BornForce      SwitchDrv                                         Forces\n" );
   double forceSum[3] = { 0.0, 0.0, 0.0 };
   RealOpenMM conversion = (RealOpenMM)(gbviParameters->getTau());  
   const RealOpenMM* scaledRadii                    = gbviParameters->getScaledRadii();
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      forceSum[0] += forces[atomI][0];
      forceSum[1] += forces[atomI][1];
      forceSum[2] += forces[atomI][2];
      (void) fprintf( logFile, "%6d %14.6e %14.6e %14.6e %14.6e [%14.6e %14.6e %14.6e]\n",
                      atomI, scaledRadii[atomI], bornRadii[atomI], conversion*bornForces[atomI], switchDeriviative[atomI],
                      conversion*forces[atomI][0], conversion*forces[atomI][1], conversion*forces[atomI][2] );
   }   
   (void) fprintf( logFile, "F sum=[%14.6e %14.6e %14.6e]\n", forceSum[0], forceSum[1], forceSum[2] );
   (void) fflush( logFile );
}
#endif

   // convert from cal to Joule & apply prefactor tau = (1/diel_solute - 1/diel_solvent)

   RealOpenMM conversion = (RealOpenMM)(gbviParameters->getTau());  
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      inputForces[atomI][0] += conversion*forces[atomI][0];
      inputForces[atomI][1] += conversion*forces[atomI][1];
      inputForces[atomI][2] += conversion*forces[atomI][2];
   }

#if( GBVISoftcoreDebug == 1 )
{
   (void) fprintf( logFile, "\nPost conversion\n" );
   (void) fprintf( logFile, "Atom        BornRadii      BornForce      SwitchDrv                                         Forces\n" );
   int maxPrint = 20;
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "%6d %14.6e %14.6e %14.6e 2[%14.6e %14.6e %14.6e] ttlF[%14.6e %14.6e %14.6e] %s\n",
                      atomI, bornRadii[atomI], conversion*bornForces[atomI], switchDeriviative[atomI],
                      conversion*forces[atomI][0], conversion*forces[atomI][1], conversion*forces[atomI][2], 
                      inputForces[atomI][0], inputForces[atomI][1], inputForces[atomI][2], 
                      (fabs( switchDeriviative[atomI] - 1.0 ) > 1.0e-05 ? "SWWWWW" : "") );
      if( atomI == maxPrint ){
         atomI = numberOfAtoms - maxPrint;
         if( atomI < maxPrint )atomI = numberOfAtoms;
      }
   }
   (void) fflush( logFile );
}
#endif
#if( GBVISoftcoreDebug == 2 )
{
   (void) fprintf( logFile, "\nAtom        BornRadii      BornForce      SwitchDrv                                         Forces  Post conversion\n" );
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      (void) fprintf( logFile, "%6d %14.6e %14.6e    %14.6e %14.6e %14.6e    %14.6e\n", atomI, bornRadii[atomI], conversion*bornForces[atomI],
                      inputForces[atomI][0], inputForces[atomI][1], inputForces[atomI][2], switchDeriviative[atomI] );
   }
   (void) fflush( logFile );
}
#endif
#undef GBVISoftcoreDebug

   delete[] forces;
   delete[] block;

   return 0;

}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      
   @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuGBVISoftcore::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << CpuImplicitSolvent::getStateString( title );

   return message.str();
}

/**---------------------------------------------------------------------------------------

   Use double precision 

   Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return volume

   --------------------------------------------------------------------------------------- */

double CpuGBVISoftcore::getVolumeD( double r, double R, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::getVolume";

   static const double zero         = 0.0;
   static const double minusThree   = -3.0;

   double              diff    = (S - R);
   if( fabs( diff ) < r ){

      double lowerBound = (R > (r - S)) ? R : (r - S);

      return (CpuGBVISoftcore::getLD( r, (r + S),    S ) -
              CpuGBVISoftcore::getLD( r, lowerBound, S ));

   } else if( r < diff ){

      return CpuGBVISoftcore::getLD( r, (r + S), S ) -
             CpuGBVISoftcore::getLD( r, (r - S), S ) + 
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

double CpuGBVISoftcore::getLD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::getL";

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

double CpuGBVISoftcore::dL_drD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::dL_dr";

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

double CpuGBVISoftcore::dL_dxD( double r, double x, double S ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "CpuGBVISoftcore::dL_dx";

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
