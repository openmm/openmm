
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
#include <cmath>
#include <cstdio>

#include "SimTKOpenMMCommon.h"
#include "ReferenceForce.h"
#include "CpuObc.h"

using namespace OpenMM;
using namespace std;

/**---------------------------------------------------------------------------------------

    CpuObc constructor

    obcParameters      obcParameters object
    
    --------------------------------------------------------------------------------------- */

CpuObc::CpuObc( ObcParameters* obcParameters ) : _obcParameters(obcParameters), _includeAceApproximation(1) {
    _obcChain.resize(_obcParameters->getNumberOfAtoms());
}

/**---------------------------------------------------------------------------------------

    CpuObc destructor

    --------------------------------------------------------------------------------------- */

CpuObc::~CpuObc( ){
}

/**---------------------------------------------------------------------------------------

    Get ObcParameters reference

    @return ObcParameters reference

    --------------------------------------------------------------------------------------- */

ObcParameters* CpuObc::getObcParameters( void ) const {
    return _obcParameters;
}

/**---------------------------------------------------------------------------------------

    Set ObcParameters reference

    @param ObcParameters reference

    --------------------------------------------------------------------------------------- */

void CpuObc::setObcParameters(  ObcParameters* obcParameters ){
    _obcParameters = obcParameters;
}

/**---------------------------------------------------------------------------------------

   Return flag signalling whether AceApproximation for nonpolar term is to be included

   @return flag

   --------------------------------------------------------------------------------------- */

int CpuObc::includeAceApproximation( void ) const {
    return _includeAceApproximation;
}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether AceApproximation is to be included

   @param includeAceApproximation new includeAceApproximation value

   --------------------------------------------------------------------------------------- */

void CpuObc::setIncludeAceApproximation( int includeAceApproximation ){
    _includeAceApproximation = includeAceApproximation;
}

/**---------------------------------------------------------------------------------------

    Return OBC chain derivative: size = _obcParameters->getNumberOfAtoms()

    @return array

    --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuObc::getObcChain( void ){
    return _obcChain;
}

/**---------------------------------------------------------------------------------------

    Get Born radii based on papers:

       J. Phys. Chem. 1996 100, 19824-19839 (HCT paper)
       Proteins: Structure, Function, and Bioinformatcis 55:383-394 (2004) (OBC paper)

    @param atomCoordinates     atomic coordinates
    @param bornRadii           output array of Born radii

    --------------------------------------------------------------------------------------- */

void CpuObc::computeBornRadii( const vector<RealVec>& atomCoordinates, vector<RealOpenMM>& bornRadii ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero    = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one     = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM two     = static_cast<RealOpenMM>( 2.0 );
    static const RealOpenMM three   = static_cast<RealOpenMM>( 3.0 );
    static const RealOpenMM half    = static_cast<RealOpenMM>( 0.5 );
    static const RealOpenMM fourth  = static_cast<RealOpenMM>( 0.25 );

    // ---------------------------------------------------------------------------------------

    ObcParameters* obcParameters                = getObcParameters();

    int numberOfAtoms                           = obcParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii         = obcParameters->getAtomicRadii();
    const RealOpenMMVector& scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();
    RealOpenMMVector& obcChain                  = getObcChain();

    RealOpenMM dielectricOffset                 = obcParameters->getDielectricOffset();
    RealOpenMM alphaObc                         = obcParameters->getAlphaObc();
    RealOpenMM betaObc                          = obcParameters->getBetaObc();
    RealOpenMM gammaObc                         = obcParameters->getGammaObc();

    // ---------------------------------------------------------------------------------------

    // calculate Born radii

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

    }
}

/**---------------------------------------------------------------------------------------

    Get nonpolar solvation force constribution via ACE approximation

    @param obcParameters parameters
    @param bornRadii                 Born radii
    @param energy                    energy (output): value is incremented from input value 
    @param forces                    forces: values are incremented from input values

    --------------------------------------------------------------------------------------- */

void CpuObc::computeAceNonPolarForce( const ObcParameters* obcParameters,
                                      const RealOpenMMVector& bornRadii, 
                                      RealOpenMM* energy,
                                      RealOpenMMVector& forces ) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero     = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM minusSix = -6.0;
    static const RealOpenMM six      = static_cast<RealOpenMM>( 6.0 );

    // ---------------------------------------------------------------------------------------

    // compute the nonpolar solvation via ACE approximation

    const RealOpenMM probeRadius          = obcParameters->getProbeRadius();
    const RealOpenMM surfaceAreaFactor    = obcParameters->getPi4Asolv();

    const RealOpenMMVector& atomicRadii   = obcParameters->getAtomicRadii();
    int numberOfAtoms                     = obcParameters->getNumberOfAtoms();

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
        if( bornRadii[atomI] > zero ){
            RealOpenMM r            = atomicRadii[atomI] + probeRadius;
            RealOpenMM ratio6       = POW( atomicRadii[atomI]/bornRadii[atomI], six );
            RealOpenMM saTerm       = surfaceAreaFactor*r*r*ratio6;
            *energy                += saTerm;
            forces[atomI]          += minusSix*saTerm/bornRadii[atomI]; 
        }
    }
}

/**---------------------------------------------------------------------------------------

    Get Obc Born energy and forces

    @param atomCoordinates     atomic coordinates
    @param partialCharges      partial charges
    @param forces              forces

    The array bornRadii is also updated and the obcEnergy

    --------------------------------------------------------------------------------------- */

RealOpenMM CpuObc::computeBornEnergyForces( const vector<RealVec>& atomCoordinates,
                                            const RealOpenMMVector& partialCharges, vector<RealVec>& inputForces ){

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero    = static_cast<RealOpenMM>( 0.0 );
    static const RealOpenMM one     = static_cast<RealOpenMM>( 1.0 );
    static const RealOpenMM two     = static_cast<RealOpenMM>( 2.0 );
    static const RealOpenMM three   = static_cast<RealOpenMM>( 3.0 );
    static const RealOpenMM four    = static_cast<RealOpenMM>( 4.0 );
    static const RealOpenMM half    = static_cast<RealOpenMM>( 0.5 );
    static const RealOpenMM fourth  = static_cast<RealOpenMM>( 0.25 );
    static const RealOpenMM eighth  = static_cast<RealOpenMM>( 0.125 );

    // constants

    const int numberOfAtoms = _obcParameters->getNumberOfAtoms();
    const RealOpenMM dielectricOffset = _obcParameters->getDielectricOffset();
    const RealOpenMM cutoffDistance = _obcParameters->getCutoffDistance();
    const RealOpenMM soluteDielectric = _obcParameters->getSoluteDielectric();
    const RealOpenMM solventDielectric = _obcParameters->getSolventDielectric();
    RealOpenMM preFactor;
    if (soluteDielectric != zero && solventDielectric != zero)
        preFactor = two*_obcParameters->getElectricConstant()*((one/soluteDielectric) - (one/solventDielectric));
    else
        preFactor = zero;

    // ---------------------------------------------------------------------------------------

    // compute Born radii

    RealOpenMMVector bornRadii( numberOfAtoms );
    computeBornRadii( atomCoordinates, bornRadii );

    // set energy/forces to zero

    RealOpenMM obcEnergy                 = zero;
    RealOpenMMVector bornForces( numberOfAtoms, 0.0 );

    // ---------------------------------------------------------------------------------------

    // compute the nonpolar solvation via ACE approximation
     
    if( includeAceApproximation() ){
       computeAceNonPolarForce( _obcParameters, bornRadii, &obcEnergy, bornForces );
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
          if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > cutoffDistance)
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
          
          RealOpenMM energy = Gpol;

          if( atomI != atomJ ){

              if (_obcParameters->getUseCutoff())
                  energy -= partialChargeI*partialCharges[atomJ]/cutoffDistance;
              
              bornForces[atomJ]        += dGpol_dalpha2_ij*bornRadii[atomI];

              deltaX                   *= dGpol_dr;
              deltaY                   *= dGpol_dr;
              deltaZ                   *= dGpol_dr;

              inputForces[atomI][0]    += deltaX;
              inputForces[atomI][1]    += deltaY;
              inputForces[atomI][2]    += deltaZ;

              inputForces[atomJ][0]    -= deltaX;
              inputForces[atomJ][1]    -= deltaY;
              inputForces[atomJ][2]    -= deltaZ;

          } else {
             energy *= half;
          }

          obcEnergy         += energy;
          bornForces[atomI] += dGpol_dalpha2_ij*bornRadii[atomJ];

       }
    }

    // ---------------------------------------------------------------------------------------

    // second main loop

    const RealOpenMMVector& obcChain            = getObcChain();
    const RealOpenMMVector& atomicRadii         = _obcParameters->getAtomicRadii();

    const RealOpenMM alphaObc                   = _obcParameters->getAlphaObc();
    const RealOpenMM betaObc                    = _obcParameters->getBetaObc();
    const RealOpenMM gammaObc                   = _obcParameters->getGammaObc();
    const RealOpenMMVector& scaledRadiusFactor  = _obcParameters->getScaledRadiusFactors();

    // compute factor that depends only on the outer loop index

    for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
       bornForces[atomI] *= bornRadii[atomI]*bornRadii[atomI]*obcChain[atomI];      
    }

    for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
 
       // radius w/ dielectric offset applied

       RealOpenMM radiusI        = atomicRadii[atomI];
       RealOpenMM offsetRadiusI  = radiusI - dielectricOffset;

       for( int atomJ = 0; atomJ < numberOfAtoms; atomJ++ ){

          if( atomJ != atomI ){

             RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
             if (_obcParameters->getPeriodic())
                ReferenceForce::getDeltaRPeriodic( atomCoordinates[atomI], atomCoordinates[atomJ], _obcParameters->getPeriodicBox(), deltaR );
             else 
                ReferenceForce::getDeltaR( atomCoordinates[atomI], atomCoordinates[atomJ], deltaR );
             if (_obcParameters->getUseCutoff() && deltaR[ReferenceForce::RIndex] > cutoffDistance)
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
    
                inputForces[atomI][0]   -= deltaX;
                inputForces[atomI][1]   -= deltaY;
                inputForces[atomI][2]   -= deltaZ;
  
                inputForces[atomJ][0]   += deltaX;
                inputForces[atomJ][1]   += deltaY;
                inputForces[atomJ][2]   += deltaZ;
 
             }
          }
       }

    }

    //printObc( atomCoordinates, partialCharges, bornRadii, bornForces, inputForces, "Obc Post loop2", stderr );

    return obcEnergy;
}

/**---------------------------------------------------------------------------------------

    Print Obc parameters, radii, forces, ...

    @param atomCoordinates     atomic coordinates
    @param partialCharges      partial charges
    @param bornRadii           Born radii (may be empty)
    @param bornForces          Born forces (may be empty)
    @param forces              forces (may be empty)
    @param idString            id string (who is calling)
    @param log                 log file

    --------------------------------------------------------------------------------------- */

void CpuObc::printObc( const std::vector<OpenMM::RealVec>& atomCoordinates,
                       const RealOpenMMVector& partialCharges,
                       const RealOpenMMVector& bornRadii,
                       const RealOpenMMVector& bornForces,
                       const std::vector<OpenMM::RealVec>& forces,
                       const std::string& idString, FILE* log ){

    // ---------------------------------------------------------------------------------------

    const ObcParameters* obcParameters          = getObcParameters();
    const int numberOfAtoms                     = obcParameters->getNumberOfAtoms();
    const RealOpenMMVector& atomicRadii         = obcParameters->getAtomicRadii();
    const RealOpenMM preFactor                  = 2.0*obcParameters->getElectricConstant();
    const RealOpenMMVector& obcChain            = getObcChain();
    const RealOpenMMVector& scaledRadiusFactor  = obcParameters->getScaledRadiusFactors();

    const RealOpenMM alphaObc                   = obcParameters->getAlphaObc();
    const RealOpenMM betaObc                    = obcParameters->getBetaObc();
    const RealOpenMM gammaObc                   = obcParameters->getGammaObc();

    const int comparisonFormat                  = 1;

    // ---------------------------------------------------------------------------------------

    (void) fprintf( log, "Reference Obc      %s atoms=%d\n", idString.c_str(), numberOfAtoms );
    if( comparisonFormat ){    
        (void) fprintf( log, "Reference Obc  %s atoms=%d Chain/Radii/Force\n", idString.c_str(), numberOfAtoms );
        for( unsigned int atomI = 0; atomI < static_cast<unsigned int>(numberOfAtoms); atomI++ ){
            (void) fprintf( log, "%6d ", atomI );
            if( obcChain.size() > atomI ){
                 (void) fprintf( log, " %15.7e", obcChain[atomI] );
            }
            if( bornRadii.size() > atomI ){
                 (void) fprintf( log, " %15.7e", bornRadii[atomI] );
            }
            if( bornForces.size() > atomI ){
                 (void) fprintf( log, " %15.7e", bornForces[atomI] );    
            }
            (void) fprintf( log, " %15.7e %6.3f", atomicRadii[atomI], partialCharges[atomI] );
            (void) fprintf( log, "\n" );
        }   
    } else { 
        (void) fprintf( log, "Reference Obc      %s atoms=%d\n", idString.c_str(), numberOfAtoms );
        (void) fprintf( log, "    preFactor      %15.7e\n", preFactor );
        (void) fprintf( log, "    alpha          %15.7e\n", alphaObc);
        (void) fprintf( log, "    beta           %15.7e\n", betaObc);
        (void) fprintf( log, "    gamma          %15.7e\n", gammaObc );
     
        for( unsigned int atomI = 0; atomI < static_cast<unsigned int>(numberOfAtoms); atomI++ ){
            (void) fprintf( log, "%6d r=%15.7e q=%6.3f", atomI,
                            atomicRadii[atomI], partialCharges[atomI] );
            if( obcChain.size() > atomI ){
                 (void) fprintf( log, " bChn=%15.7e", obcChain[atomI] );
            }
            if( bornRadii.size() > atomI ){
                 (void) fprintf( log, " bR=%15.7e", bornRadii[atomI] );
            }
            if( bornForces.size() > atomI ){
                 (void) fprintf( log, " bF=%15.7e", bornForces[atomI] );    
            }
            (void) fprintf( log, "\n" );
        }   
    }   
    
    return;

}

