
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

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceVdwForce.h"
#include <algorithm>
#include <cctype>

using std::vector;
using OpenMM::RealVec;

AmoebaReferenceVdwForce::AmoebaReferenceVdwForce( ) : _nonbondedMethod(NoCutoff) {

    _cutoff  = 1.0e+10;
    setSigmaCombiningRule( "ARITHMETIC" );
    setEpsilonCombiningRule( "GEOMETRIC" );
}

AmoebaReferenceVdwForce::AmoebaReferenceVdwForce( const std::string& sigmaCombiningRule, const std::string& epsilonCombiningRule ) : _nonbondedMethod(NoCutoff) {

    setSigmaCombiningRule( sigmaCombiningRule );
    setEpsilonCombiningRule( epsilonCombiningRule );
}

AmoebaReferenceVdwForce::NonbondedMethod AmoebaReferenceVdwForce::getNonbondedMethod( void ) const {
    return _nonbondedMethod;
}

void AmoebaReferenceVdwForce::setNonbondedMethod( AmoebaReferenceVdwForce::NonbondedMethod nonbondedMethod ){
    _nonbondedMethod = nonbondedMethod;
}

void AmoebaReferenceVdwForce::setCutoff( double cutoff ){
    _cutoff = cutoff;
}

double AmoebaReferenceVdwForce::getCutoff( void ) const {
    return _cutoff;
}

void AmoebaReferenceVdwForce::setSigmaCombiningRule( const std::string& sigmaCombiningRule ){

    _sigmaCombiningRule = sigmaCombiningRule;

    // convert to upper case and set combining function

    std::transform( _sigmaCombiningRule.begin(), _sigmaCombiningRule.end(), _sigmaCombiningRule.begin(),  (int(*)(int)) std::toupper);
    if( _sigmaCombiningRule == "GEOMETRIC" ){
        _combineSigmas = &AmoebaReferenceVdwForce::geometricSigmaCombiningRule;
    } else if( _sigmaCombiningRule == "CUBIC-MEAN" ){
        _combineSigmas = &AmoebaReferenceVdwForce::cubicMeanSigmaCombiningRule;
    } else {
        _combineSigmas = &AmoebaReferenceVdwForce::arithmeticSigmaCombiningRule;
    }
}

std::string AmoebaReferenceVdwForce::getSigmaCombiningRule( void ) const {
    return _sigmaCombiningRule;
}

RealOpenMM AmoebaReferenceVdwForce::arithmeticSigmaCombiningRule( RealOpenMM sigmaI, RealOpenMM sigmaJ ) const {
    return (sigmaI + sigmaJ);
}

RealOpenMM AmoebaReferenceVdwForce::geometricSigmaCombiningRule( RealOpenMM sigmaI, RealOpenMM sigmaJ ) const {
    return 2.0*SQRT(sigmaI*sigmaJ);
}

RealOpenMM AmoebaReferenceVdwForce::cubicMeanSigmaCombiningRule( RealOpenMM sigmaI, RealOpenMM sigmaJ ) const {

    const RealOpenMM zero = 0.0;

    RealOpenMM sigmaI2    = sigmaI*sigmaI;
    RealOpenMM sigmaJ2    = sigmaJ*sigmaJ;

    return sigmaI != zero && sigmaJ != 0.0 ? 2.0*(sigmaI2*sigmaI + sigmaJ2*sigmaJ)/(sigmaI2 + sigmaJ2) : zero;
}

void AmoebaReferenceVdwForce::setEpsilonCombiningRule( const std::string& epsilonCombiningRule ){

    _epsilonCombiningRule = epsilonCombiningRule;
    std::transform( _epsilonCombiningRule.begin(), _epsilonCombiningRule.end(), _epsilonCombiningRule.begin(),  (int(*)(int)) std::toupper);

    // convert to upper case and set combining function

    if( _epsilonCombiningRule == "ARITHMETIC" ){
         _combineEpsilons = &AmoebaReferenceVdwForce::arithmeticEpsilonCombiningRule;
    } else if( _epsilonCombiningRule == "HARMONIC" ){
         _combineEpsilons = &AmoebaReferenceVdwForce::harmonicEpsilonCombiningRule;
    } else if( _epsilonCombiningRule == "HHG" ){
         _combineEpsilons = &AmoebaReferenceVdwForce::hhgEpsilonCombiningRule;
    } else {
         _combineEpsilons = &AmoebaReferenceVdwForce::geometricEpsilonCombiningRule;
    }
}

std::string AmoebaReferenceVdwForce::getEpsilonCombiningRule( void ) const {
    return _epsilonCombiningRule;
}

RealOpenMM AmoebaReferenceVdwForce::arithmeticEpsilonCombiningRule( RealOpenMM epsilonI, RealOpenMM epsilonJ ) const {
    return 0.5*(epsilonI + epsilonJ);
}

RealOpenMM AmoebaReferenceVdwForce::geometricEpsilonCombiningRule( RealOpenMM epsilonI, RealOpenMM epsilonJ ) const {
    return SQRT(epsilonI*epsilonJ);
}

RealOpenMM AmoebaReferenceVdwForce::harmonicEpsilonCombiningRule( RealOpenMM epsilonI, RealOpenMM epsilonJ ) const {
    return (epsilonI != 0.0 && epsilonJ != 0.0) ? 2.0*(epsilonI*epsilonJ)/(epsilonI + epsilonJ) : 0.0;
}

RealOpenMM AmoebaReferenceVdwForce::hhgEpsilonCombiningRule( RealOpenMM epsilonI, RealOpenMM epsilonJ ) const {
    RealOpenMM denominator = SQRT(epsilonI) + SQRT(epsilonJ);
    return (epsilonI != 0.0 && epsilonJ != 0.0) ? 4.0*(epsilonI*epsilonJ)/(denominator*denominator) : 0.0;
}

void AmoebaReferenceVdwForce::addReducedForce( unsigned int particleI, unsigned int particleIV,
                                               RealOpenMM reduction, RealOpenMM sign,
                                               Vec3& force, vector<RealVec>& forces ) const {

// ---------------------------------------------------------------------------------------

    static const RealOpenMM one           = 1.0;

   // ---------------------------------------------------------------------------------------

    forces[particleI][0]  += sign*force[0]*reduction;
    forces[particleI][1]  += sign*force[1]*reduction;
    forces[particleI][2]  += sign*force[2]*reduction;

    forces[particleIV][0] += sign*force[0]*(one - reduction);
    forces[particleIV][1] += sign*force[1]*(one - reduction);
    forces[particleIV][2] += sign*force[2]*(one - reduction);
}

RealOpenMM AmoebaReferenceVdwForce::calculatePairIxn( RealOpenMM combindedSigma, RealOpenMM combindedEpsilon,
                                                      const Vec3& particleIPosition,
                                                      const Vec3& particleJPosition,
                                                      Vec3& force ) const {

   // ---------------------------------------------------------------------------------------

    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM seven         = 7.0;

    static const RealOpenMM dhal          = 0.07;
    static const RealOpenMM ghal          = 0.12;

   // ---------------------------------------------------------------------------------------

    RealOpenMM xr           = particleIPosition[0] - particleJPosition[0];
    RealOpenMM yr           = particleIPosition[1] - particleJPosition[1];
    RealOpenMM zr           = particleIPosition[2] - particleJPosition[2];
    RealOpenMM r_ij_2       = xr*xr + yr*yr + zr*zr;
                 
    RealOpenMM sigma_7      = combindedSigma*combindedSigma*combindedSigma;
               sigma_7      = sigma_7*sigma_7*combindedSigma;

    RealOpenMM r_ij         = SQRT(r_ij_2);
    RealOpenMM r_ij_6       = r_ij_2*r_ij_2*r_ij_2;
    RealOpenMM r_ij_7       = r_ij_6*r_ij;

    RealOpenMM rho          = r_ij_7 + ghal*sigma_7;

    RealOpenMM tau          = (dhal + one)/(r_ij + dhal*combindedSigma);
    RealOpenMM tau_7        = tau*tau*tau;
               tau_7        = tau_7*tau_7*tau;

    RealOpenMM dtau         = tau/(dhal + one);

    RealOpenMM ratio        = (sigma_7/rho);
    RealOpenMM gtau         = combindedEpsilon*tau_7*r_ij_6*(ghal+one)*ratio*ratio;

    RealOpenMM energy       = combindedEpsilon*tau_7*sigma_7*( (ghal+one)*sigma_7/rho - two);

    RealOpenMM dEdR         = -seven*(dtau*energy + gtau);
               dEdR        /= r_ij;

    force[0]                = dEdR*xr;
    force[1]                = dEdR*yr;
    force[2]                = dEdR*zr;

   return energy;

}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergyNoCutoff( int numParticles,
                                                                     const vector<RealVec>& particlePositions,
                                                                     const std::vector<int>& indexIVs, 
                                                                     const std::vector<RealOpenMM>& sigmas,
                                                                     const std::vector<RealOpenMM>& epsilons,
                                                                     const std::vector<RealOpenMM>& reductions,
                                                                     const std::vector< std::vector<int> >& allExclusions,
                                                                     vector<RealVec>& forces ) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;

    // ---------------------------------------------------------------------------------------

    // set reduced coordinates

    std::vector<Vec3> reducedPositions;
    reducedPositions.resize(numParticles);
    for( unsigned int ii = 0; ii <  static_cast<unsigned int>(numParticles); ii++ ){
        if( reductions[ii] != zero ){
            int reductionIndex     = indexIVs[ii];
            reducedPositions[ii]   = Vec3( reductions[ii]*( particlePositions[ii][0] - particlePositions[reductionIndex][0] ) + particlePositions[reductionIndex][0], 
                                           reductions[ii]*( particlePositions[ii][1] - particlePositions[reductionIndex][1] ) + particlePositions[reductionIndex][1], 
                                           reductions[ii]*( particlePositions[ii][2] - particlePositions[reductionIndex][2] ) + particlePositions[reductionIndex][2] ); 
        } else {
            reducedPositions[ii] = Vec3( particlePositions[ii][0], particlePositions[ii][1], particlePositions[ii][2] ); 
        }
    }
 
    // loop over all ixns
    //    (1) initialize exclusion vector

    RealOpenMM energy = zero;
    std::vector<unsigned int> exclusions(numParticles, 0);
    for( unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++ ){
 
        RealOpenMM sigmaI      = sigmas[ii];
        RealOpenMM epsilonI    = epsilons[ii];
        for( unsigned int jj = 0; jj < allExclusions[ii].size(); jj++ ){
            exclusions[allExclusions[ii][jj]] = 1;
        }

        for( unsigned int jj = ii+1; jj < static_cast<unsigned int>(numParticles); jj++ ){
            if( exclusions[jj] == 0 ){

                RealOpenMM combindedSigma   = (this->*_combineSigmas)(sigmaI, sigmas[jj] );
                RealOpenMM combindedEpsilon = (this->*_combineEpsilons)(epsilonI, epsilons[jj] );

                Vec3 force;
                energy                     += calculatePairIxn( combindedSigma, combindedEpsilon,
                                                                reducedPositions[ii], reducedPositions[jj],
                                                                force );
                
                if( indexIVs[ii] == ii ){
                    forces[ii][0] -= force[0];
                    forces[ii][1] -= force[1];
                    forces[ii][2] -= force[2];
                } else {
                    addReducedForce( ii, indexIVs[ii], reductions[ii], -one, force, forces );
                }
                if( indexIVs[jj] == jj ){
                    forces[jj][0] += force[0];
                    forces[jj][1] += force[1];
                    forces[jj][2] += force[2];
                } else {
                    addReducedForce( jj, indexIVs[jj], reductions[jj], one, force, forces );
                }

            }
        }

        for( unsigned int jj = 0; jj < allExclusions[ii].size(); jj++ ){
            exclusions[allExclusions[ii][jj]] = 0;
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergyApplyCutoff( int numParticles,
                                                                        const vector<RealVec>& particlePositions,
                                                                        const std::vector<int>& indexIVs, 
                                                                        const std::vector<RealOpenMM>& sigmas,
                                                                        const std::vector<RealOpenMM>& epsilons,
                                                                        const std::vector<RealOpenMM>& reductions,
                                                                        const std::vector< std::vector<int> >& allExclusions,
                                                                        vector<RealVec>& forces ) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;

    // ---------------------------------------------------------------------------------------

    // set reduced coordinates

    std::vector<Vec3> reducedPositions;
    reducedPositions.resize(numParticles);
    for( unsigned int ii = 0; ii <  static_cast<unsigned int>(numParticles); ii++ ){
        if( reductions[ii] != zero ){
            int reductionIndex     = indexIVs[ii];
            reducedPositions[ii]   = Vec3( reductions[ii]*( particlePositions[ii][0] - particlePositions[reductionIndex][0] ) + particlePositions[reductionIndex][0], 
                                           reductions[ii]*( particlePositions[ii][1] - particlePositions[reductionIndex][1] ) + particlePositions[reductionIndex][1], 
                                           reductions[ii]*( particlePositions[ii][2] - particlePositions[reductionIndex][2] ) + particlePositions[reductionIndex][2] ); 
        } else {
            reducedPositions[ii] = Vec3( particlePositions[ii][0], particlePositions[ii][1], particlePositions[ii][2] ); 
        }
    }
 
    // loop over all ixns
    //    (1) initialize exclusion vector

    RealOpenMM energy = zero;
    std::vector<unsigned int> exclusions(numParticles, 0);
    for( unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++ ){
 
        RealOpenMM sigmaI      = sigmas[ii];
        RealOpenMM epsilonI    = epsilons[ii];
        for( unsigned int jj = 0; jj < allExclusions[ii].size(); jj++ ){
            exclusions[allExclusions[ii][jj]] = 1;
        }

        for( unsigned int jj = ii+1; jj < static_cast<unsigned int>(numParticles); jj++ ){
            if( exclusions[jj] == 0 ){

                RealOpenMM combindedSigma   = (this->*_combineSigmas)(sigmaI, sigmas[jj] );
                RealOpenMM combindedEpsilon = (this->*_combineEpsilons)(epsilonI, epsilons[jj] );

                Vec3 force;
                energy                     += calculatePairIxn( combindedSigma, combindedEpsilon,
                                                                reducedPositions[ii], reducedPositions[jj],
                                                                force );
                
                if( indexIVs[ii] == ii ){
                    forces[ii][0] -= force[0];
                    forces[ii][1] -= force[1];
                    forces[ii][2] -= force[2];
                } else {
                    addReducedForce( ii, indexIVs[ii], reductions[ii], -one, force, forces );
                }
                if( indexIVs[jj] == jj ){
                    forces[jj][0] += force[0];
                    forces[jj][1] += force[1];
                    forces[jj][2] += force[2];
                } else {
                    addReducedForce( jj, indexIVs[jj], reductions[jj], one, force, forces );
                }

            }
        }

        for( unsigned int jj = 0; jj < allExclusions[ii].size(); jj++ ){
            exclusions[allExclusions[ii][jj]] = 0;
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergy( int numParticles,
                                                             const vector<RealVec>& particlePositions,
                                                             const std::vector<int>& indexIVs, 
                                                             const std::vector<RealOpenMM>& sigmas,
                                                             const std::vector<RealOpenMM>& epsilons,
                                                             const std::vector<RealOpenMM>& reductions,
                                                             const std::vector< std::vector<int> >& allExclusions,
                                                             vector<RealVec>& forces ) const {

    
    if( getNonbondedMethod() == NoCutoff ){
        return calculateForceAndEnergyNoCutoff( numParticles, particlePositions, 
                                                indexIVs, sigmas, epsilons,
                                                reductions, allExclusions, forces );
    } else {
        return calculateForceAndEnergyApplyCutoff( numParticles, particlePositions, 
                                                   indexIVs, sigmas, epsilons,
                                                   reductions, allExclusions, forces );
    }

}    
/*
    double cutoff = force.getCutoff();
    double taperCutoff = cutoff*0.9;
    replacements["CUTOFF_DISTANCE"] = cu.doubleToString(force.getCutoff());
    replacements["TAPER_CUTOFF"] = cu.doubleToString(taperCutoff);
    replacements["TAPER_C3"] = cu.doubleToString(10/pow(taperCutoff-cutoff, 3.0));
    replacements["TAPER_C4"] = cu.doubleToString(15/pow(taperCutoff-cutoff, 4.0));
    replacements["TAPER_C5"] = cu.doubleToString(6/pow(taperCutoff-cutoff, 5.0));

*/
