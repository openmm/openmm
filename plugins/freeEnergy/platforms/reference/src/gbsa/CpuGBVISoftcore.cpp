
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
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <vector>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKReference/ReferenceForce.h"
#include "CpuGBVISoftcore.h"

using std::vector;
using OpenMM::RealVec;

/**---------------------------------------------------------------------------------------

   CpuGBVISoftcore constructor

   gbviParameters      gbviParameters object
   
   --------------------------------------------------------------------------------------- */

CpuGBVISoftcore::CpuGBVISoftcore( GBVISoftcoreParameters* gbviParameters ) : _gbviParameters(gbviParameters) {
    _switchDeriviative.resize(_gbviParameters->getNumberOfAtoms());

}

/**---------------------------------------------------------------------------------------

   CpuGBVISoftcore destructor

   --------------------------------------------------------------------------------------- */

CpuGBVISoftcore::~CpuGBVISoftcore( ){
}

/**---------------------------------------------------------------------------------------

   Get GBVISoftcoreParameters reference

   @return GBVISoftcoreParameters reference

   --------------------------------------------------------------------------------------- */

GBVISoftcoreParameters* CpuGBVISoftcore::getGBVISoftcoreParameters( void ) const {
    return _gbviParameters;
}

/**---------------------------------------------------------------------------------------

   Set GBVISoftcoreParameters reference

   @param GBVISoftcoreParameters reference

   --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::setGBVISoftcoreParameters( GBVISoftcoreParameters* gbviParameters ){
    _gbviParameters = gbviParameters;
}

/**---------------------------------------------------------------------------------------

   Return OBC chain derivative: size = _obcSoftcoreParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMMVector& CpuGBVISoftcore::getSwitchDeriviative( void ){
    return _switchDeriviative;
}

/**---------------------------------------------------------------------------------------

   Compute quintic spline value and associated derviative

   @param x                   value to compute spline at
   @param rl                  lower cutoff value
   @param ru                  upper cutoff value
   @param outValue            value of spline at x
   @param outDerivative       value of derivative of spline at x

   --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::quinticSpline( RealOpenMM x, RealOpenMM rl, RealOpenMM ru,
                                     RealOpenMM* outValue, RealOpenMM* outDerivative ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM one           = static_cast<RealOpenMM>(   1.0 );
    static const RealOpenMM minusSix      = static_cast<RealOpenMM>(  -6.0 );
    static const RealOpenMM minusTen      = static_cast<RealOpenMM>( -10.0 );
    static const RealOpenMM minusThirty   = static_cast<RealOpenMM>( -30.0 );
    static const RealOpenMM fifteen       = static_cast<RealOpenMM>(  15.0 );
    static const RealOpenMM sixty         = static_cast<RealOpenMM>(  60.0 );

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

void CpuGBVISoftcore::computeBornRadiiUsingQuinticSpline( RealOpenMM atomicRadius3, RealOpenMM bornSum,
                                                          GBVISoftcoreParameters* gbviParameters, 
                                                          RealOpenMM& bornRadius, RealOpenMM* switchDeriviative ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = static_cast<RealOpenMM>(  0.0 );
    static const RealOpenMM one           = static_cast<RealOpenMM>(  1.0 );
    static const RealOpenMM minusOne      = static_cast<RealOpenMM>( -1.0 );
    static const RealOpenMM minusThree    = static_cast<RealOpenMM>( -3.0 );
    static const RealOpenMM oneEighth     = static_cast<RealOpenMM>(  0.125 );
    static const RealOpenMM minusOneThird = static_cast<RealOpenMM>( (-1.0/3.0) );
    static const RealOpenMM three         = static_cast<RealOpenMM>(  3.0 );

    // ---------------------------------------------------------------------------------------

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
        } else {   
            sum                = gbviParameters->getQuinticUpperSplineLimit();
            *switchDeriviative = zero;
        }
    } else {
        sum                = atomicRadius3 - bornSum; 
        *switchDeriviative = one;
    }
    bornRadius = POW( sum, minusOneThird );
}

/**---------------------------------------------------------------------------------------

   Get Born radii based on Eq. 3 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii
   @param chain               not used here

   --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::computeBornRadii( std::vector<OpenMM::RealVec>& atomCoordinates,
                                        RealOpenMMVector& bornRadii ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one           = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM minusThree    = static_cast<RealOpenMM>( -3.0 );
    static const RealOpenMM oneEighth     = static_cast<RealOpenMM>( 0.125 );
    static const RealOpenMM minusOneThird = static_cast<RealOpenMM>( (-1.0/3.0) );
    static const RealOpenMM three         = static_cast<RealOpenMM>( 3.0 );

    // ---------------------------------------------------------------------------------------

    GBVISoftcoreParameters* gbviParameters          = getGBVISoftcoreParameters();
    int numberOfAtoms                               = gbviParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii             = gbviParameters->getAtomicRadii();
    const RealOpenMMVector& scaledRadii             = gbviParameters->getScaledRadii();
    const RealOpenMMVector& bornRadiusScaleFactors  = gbviParameters->getBornRadiusScaleFactors();

    RealOpenMMVector& switchDeriviative             = getSwitchDeriviative();

    // ---------------------------------------------------------------------------------------

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

            }
        }

        RealOpenMM atomicRadius3 = POW( radiusI, minusThree );
        if( _gbviParameters->getBornRadiusScalingSoftcoreMethod() == GBVISoftcoreParameters::NoScaling ){
            sum                       = atomicRadius3 - sum;
            bornRadii[atomI]          = POW( sum, minusOneThird );
            switchDeriviative[atomI]  = one;
        } else if( _gbviParameters->getBornRadiusScalingSoftcoreMethod() == GBVISoftcoreParameters::QuinticSpline ){
            RealOpenMM switchDeriviativeValue;
            computeBornRadiiUsingQuinticSpline( atomicRadius3, sum, gbviParameters, 
                                                bornRadii[atomI], &switchDeriviativeValue );
            switchDeriviative[atomI] = switchDeriviativeValue;

        }

    }

}

/**---------------------------------------------------------------------------------------

   Get volume Eq. 4 of Labute paper [JCC 29 p. 1693-1698 2008])

   @param r                   distance between atoms i & j
   @param R                   atomic radius
   @param S                   scaled atomic radius

   @return volume

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVISoftcore::getVolume( RealOpenMM r, RealOpenMM R, RealOpenMM S ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero         = static_cast<RealOpenMM>(  0.0 );
    static const RealOpenMM minusThree   = static_cast<RealOpenMM>( -3.0 );

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

    static const RealOpenMM one           = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM threeHalves   = static_cast<RealOpenMM>( 1.5 );
    static const RealOpenMM third         = static_cast<RealOpenMM>( (1.0/3.0) );
    static const RealOpenMM fourth        = static_cast<RealOpenMM>( 0.25 );
    static const RealOpenMM eighth        = static_cast<RealOpenMM>( 0.125 );

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

    static const RealOpenMM one           = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM threeHalves   = static_cast<RealOpenMM>( 1.5 );
    static const RealOpenMM threeEights   = static_cast<RealOpenMM>( 0.375 );
    static const RealOpenMM third         = static_cast<RealOpenMM>( (1.0/3.0) );
    static const RealOpenMM fourth        = static_cast<RealOpenMM>( 0.25 );
    static const RealOpenMM eighth        = static_cast<RealOpenMM>( 0.125 );

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

    static const RealOpenMM one           = static_cast<RealOpenMM>(  1.0 );
    static const RealOpenMM half          = static_cast<RealOpenMM>(  0.5 );
    static const RealOpenMM threeHalvesM  = static_cast<RealOpenMM>( -1.5 );
    static const RealOpenMM third         = static_cast<RealOpenMM>(  (1.0/3.0) );

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

    static const RealOpenMM zero    = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one     = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM fourth  = static_cast<RealOpenMM>( 0.25 );

    // ---------------------------------------------------------------------------------------

    return ( (t != zero) ? one/SQRT( (one + (fourth*EXP( -t ))/t) ) : zero);
}

/**---------------------------------------------------------------------------------------

   Get GB/VI energy

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuGBVISoftcore::computeBornEnergy( const vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                               const RealOpenMM* partialCharges ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero    = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one     = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM two     = static_cast<RealOpenMM>( 2.0 );
    static const RealOpenMM three   = static_cast<RealOpenMM>( 3.0 );
    static const RealOpenMM four    = static_cast<RealOpenMM>( 4.0 );
    static const RealOpenMM half    = static_cast<RealOpenMM>( 0.5 );
    static const RealOpenMM fourth  = static_cast<RealOpenMM>( 0.25 );
    static const RealOpenMM eighth  = static_cast<RealOpenMM>( 0.125 );

    // ---------------------------------------------------------------------------------------

    const GBVISoftcoreParameters* gbviParameters       = getGBVISoftcoreParameters();
    const RealOpenMM preFactor                         = gbviParameters->getElectricConstant();
    const int numberOfAtoms                            = gbviParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii                = gbviParameters->getAtomicRadii();
    const RealOpenMMVector& gammaParameters            = gbviParameters->getGammaParameters();

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
        }

        energy += two*partialChargeI*atomIEnergy;
    }
    energy *= preFactor;
    energy -= cavityEnergy;

    RealOpenMM conversion = static_cast<RealOpenMM>(gbviParameters->getTau());  
    return (conversion*energy);
 
}

/**---------------------------------------------------------------------------------------

   Get GB/VI forces

   @param bornRadii           Born radii
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::computeBornForces( const vector<RealOpenMM>& bornRadii, vector<RealVec>& atomCoordinates,
                                         const RealOpenMM* partialCharges, vector<RealVec>& inputForces ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero         = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one          = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM two          = static_cast<RealOpenMM>( 2.0 );
    static const RealOpenMM three        = static_cast<RealOpenMM>( 3.0 );
    static const RealOpenMM four         = static_cast<RealOpenMM>( 4.0 );
    static const RealOpenMM half         = static_cast<RealOpenMM>( 0.5 );
    static const RealOpenMM oneThird     = static_cast<RealOpenMM>( (1.0/3.0) );
    static const RealOpenMM fourth       = static_cast<RealOpenMM>( 0.25 );
    static const RealOpenMM eighth       = static_cast<RealOpenMM>( 0.125 );

    // ---------------------------------------------------------------------------------------

    const GBVISoftcoreParameters* gbviParameters       = getGBVISoftcoreParameters();
    const int numberOfAtoms                            = gbviParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii                = gbviParameters->getAtomicRadii();
    const RealOpenMMVector& gammaParameters            = gbviParameters->getGammaParameters();
    const RealOpenMM preFactor                         = two*gbviParameters->getElectricConstant();

    // ---------------------------------------------------------------------------------------

    // set energy/forces to zero

    RealOpenMMVector bornForces( numberOfAtoms, 0.0 );
    std::vector<RealVec> forces( numberOfAtoms, RealVec(0.0, 0.0, 0.0) );

    // ---------------------------------------------------------------------------------------

    // first main loop

#undef TARGET
#define TARGET -1926
RealOpenMMVector bornForcesT( numberOfAtoms, 0.0 );

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

if( atomJ == TARGET ){
bornForcesT[atomI] = dGpol_dalpha2_ij*bornRadii[atomI];
}

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
           bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];
if( atomI == TARGET ){
bornForcesT[atomJ] = dGpol_dalpha2_ij*bornRadii[atomJ];
}

       }
    }


    // second main loop: (dGpol/dBornRadius)(dBornRadius/dr)(dr/dx)

    // dGpol/dBornRadius) = bornForces[]
    // dBornRadius/dr     = (1/3)*(bR**4)*(dV/dr)

    const RealOpenMMVector& scaledRadii                 = gbviParameters->getScaledRadii();
    RealOpenMMVector& switchDeriviative                 = getSwitchDeriviative();
    const RealOpenMMVector& bornRadiusScaleFactors      = gbviParameters->getBornRadiusScaleFactors();

if( 0 ){
double bF = 0.0;
int count = 0;
fprintf( stderr, "PdeRef %d\n", TARGET );
RealOpenMM tau      = static_cast<RealOpenMM>(gbviParameters->getTau());
for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
    bF += bornForcesT[atomI];
    if( fabs( bornForcesT[atomI] ) > 0.0 ){
        count++;
        fprintf( stderr, "%6d %6d bFAdd=%15.7e bR=%15.7e\n", atomI, count, tau*bornForcesT[atomI], bornRadii[atomI] );
    }
}
RealOpenMM  ratio   = (atomicRadii[TARGET]/bornRadii[TARGET]);
RealOpenMM    bF1   = bF + (three*gammaParameters[TARGET]*ratio*ratio*ratio)/bornRadii[TARGET];
RealOpenMM b2       = bornRadii[TARGET]*bornRadii[TARGET];
bF1                *= switchDeriviative[TARGET]*oneThird*b2*b2;
fprintf( stderr, "sumbF Ref %6d %15.7e %15.7e %15.7e %15.7e\n", TARGET, bF, tau*bF1, bornRadii[TARGET], switchDeriviative[TARGET] );
}

    // ---------------------------------------------------------------------------------------

    for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
        RealOpenMM R        = atomicRadii[atomI];

        // partial of cavity term wrt Born radius
       
        RealOpenMM  ratio   = (atomicRadii[atomI]/bornRadii[atomI]);
        bornForces[atomI]  += (three*gammaParameters[atomI]*ratio*ratio*ratio)/bornRadii[atomI]; 

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

          }
       }

    }

    //printGbvi( atomCoordinates, partialCharges, bornRadii,bornForces, forces, "Reference: GBVI softcore post loop2", stderr );

    // apply prefactor tau = (1/diel_solute - 1/diel_solvent)

    RealOpenMM conversion = static_cast<RealOpenMM>(gbviParameters->getTau());  
    for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
        inputForces[atomI][0] += conversion*forces[atomI][0];
        inputForces[atomI][1] += conversion*forces[atomI][1];
        inputForces[atomI][2] += conversion*forces[atomI][2];
    }

}

/**---------------------------------------------------------------------------------------

    Print GB/VI parameters, radii, forces, ...

    @param atomCoordinates     atomic coordinates
    @param partialCharges      partial charges
    @param bornRadii           Born radii (may be empty)
    @param bornForces          Born forces (may be empty)
    @param forces              forces (may be empty)
    @param idString            id string (who is calling)
    @param log                 log file

    --------------------------------------------------------------------------------------- */

void CpuGBVISoftcore::printGbvi( const std::vector<OpenMM::RealVec>& atomCoordinates,
                                 const RealOpenMM* partialCharges,
                                 const RealOpenMMVector& bornRadii,
                                 const RealOpenMMVector& bornForces,
                                 const std::vector<OpenMM::RealVec>& forces,
                                 const std::string& idString, FILE* log ){

    // ---------------------------------------------------------------------------------------

    const GBVISoftcoreParameters* gbviParameters    = getGBVISoftcoreParameters();
    const int numberOfAtoms                         = gbviParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii             = gbviParameters->getAtomicRadii();
    const RealOpenMMVector& gammaParameters         = gbviParameters->getGammaParameters();
    const RealOpenMMVector& bornRadiusScaleFactors  = gbviParameters->getBornRadiusScaleFactors();
    const RealOpenMM preFactor                      = 2.0*gbviParameters->getElectricConstant();

    // ---------------------------------------------------------------------------------------

    const RealOpenMMVector& scaledRadii       = gbviParameters->getScaledRadii();
    const RealOpenMMVector& switchDeriviative = getSwitchDeriviative();
    RealOpenMM tau                            = static_cast<RealOpenMM>(gbviParameters->getTau());  

    (void) fprintf( log, "Reference Gbvi     %s atoms=%d\n", idString.c_str(), numberOfAtoms );
    (void) fprintf( log, "    tau            %15.7e\n", tau ); 
    (void) fprintf( log, "    scaleMethod    %d (QuinticEnum=%d)\n", 
                    gbviParameters->getBornRadiusScalingSoftcoreMethod(), GBVISoftcoreParameters::QuinticSpline );
    (void) fprintf( log, "    preFactor      %15.7e)\n", preFactor );
 
    for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
        (void) fprintf( log, "%6d r=%15.7e rSc=%15.7e swd=%15.7e lmda=%4.2f tau*gam=%15.7e q=%15.7e", atomI,
                        atomicRadii[atomI],    
                        scaledRadii[atomI],    
                        switchDeriviative[atomI],    
                        bornRadiusScaleFactors[atomI],
                        tau*gammaParameters[atomI],    
                        partialCharges[atomI] );
        if( bornRadii.size() > atomI ){
             (void) fprintf( log, " bR=%15.7e", bornRadii[atomI] );
        }
        if( bornForces.size() > atomI ){
             (void) fprintf( log, " tau*bF=%15.7e", tau*bornForces[atomI] );    
        }
        (void) fprintf( log, "\n" );
    }   

    return;

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
