
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
#include "ReferenceForce.h"
#include <algorithm>
#include <cctype>

using std::vector;
using namespace OpenMM;

AmoebaReferenceVdwForce::AmoebaReferenceVdwForce() : _nonbondedMethod(NoCutoff), _cutoff(1.0e+10), _taperCutoffFactor(0.9) {

    setTaperCoefficients(_cutoff);
    setSigmaCombiningRule("ARITHMETIC");
    setEpsilonCombiningRule("GEOMETRIC");
}


AmoebaReferenceVdwForce::AmoebaReferenceVdwForce(const std::string& sigmaCombiningRule, const std::string& epsilonCombiningRule) : _nonbondedMethod(NoCutoff), _cutoff(1.0e+10), _taperCutoffFactor(0.9) {

    setTaperCoefficients(_cutoff);
    setSigmaCombiningRule(sigmaCombiningRule);
    setEpsilonCombiningRule(epsilonCombiningRule);
}

AmoebaReferenceVdwForce::NonbondedMethod AmoebaReferenceVdwForce::getNonbondedMethod() const {
    return _nonbondedMethod;
}

void AmoebaReferenceVdwForce::setNonbondedMethod(AmoebaReferenceVdwForce::NonbondedMethod nonbondedMethod) {
    _nonbondedMethod = nonbondedMethod;
}

void AmoebaReferenceVdwForce::setTaperCoefficients(double cutoff) {
    _taperCutoff = cutoff*_taperCutoffFactor;
    if (_taperCutoff != cutoff) {
        _taperCoefficients[C3] = 10.0/pow(_taperCutoff - cutoff, 3.0);
        _taperCoefficients[C4] = 15.0/pow(_taperCutoff - cutoff, 4.0);
        _taperCoefficients[C5] =  6.0/pow(_taperCutoff - cutoff, 5.0);
    } else {
        _taperCoefficients[C3] = 0.0;
        _taperCoefficients[C4] = 0.0;
        _taperCoefficients[C5] = 0.0;
    }
}

void AmoebaReferenceVdwForce::setCutoff(double cutoff) {
    _cutoff  = cutoff;
    setTaperCoefficients(_cutoff);
}

double AmoebaReferenceVdwForce::getCutoff() const {
    return _cutoff;
}

void AmoebaReferenceVdwForce::setPeriodicBox(OpenMM::RealVec* vectors) {
    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
}

void AmoebaReferenceVdwForce::setSigmaCombiningRule(const std::string& sigmaCombiningRule) {

    _sigmaCombiningRule = sigmaCombiningRule;

    // convert to upper case and set combining function

    std::transform(_sigmaCombiningRule.begin(), _sigmaCombiningRule.end(), _sigmaCombiningRule.begin(),  (int(*)(int)) std::toupper);
    if (_sigmaCombiningRule == "GEOMETRIC") {
        _combineSigmas = &AmoebaReferenceVdwForce::geometricSigmaCombiningRule;
    } else if (_sigmaCombiningRule == "CUBIC-MEAN") {
        _combineSigmas = &AmoebaReferenceVdwForce::cubicMeanSigmaCombiningRule;
    } else {
        _combineSigmas = &AmoebaReferenceVdwForce::arithmeticSigmaCombiningRule;
    }
}

std::string AmoebaReferenceVdwForce::getSigmaCombiningRule() const {
    return _sigmaCombiningRule;
}

RealOpenMM AmoebaReferenceVdwForce::arithmeticSigmaCombiningRule(RealOpenMM sigmaI, RealOpenMM sigmaJ) const {
    return (sigmaI + sigmaJ);
}

RealOpenMM AmoebaReferenceVdwForce::geometricSigmaCombiningRule(RealOpenMM sigmaI, RealOpenMM sigmaJ) const {
    return 2.0*SQRT(sigmaI*sigmaJ);
}

RealOpenMM AmoebaReferenceVdwForce::cubicMeanSigmaCombiningRule(RealOpenMM sigmaI, RealOpenMM sigmaJ) const {

    const RealOpenMM zero = 0.0;

    RealOpenMM sigmaI2    = sigmaI*sigmaI;
    RealOpenMM sigmaJ2    = sigmaJ*sigmaJ;

    return sigmaI != zero && sigmaJ != 0.0 ? 2.0*(sigmaI2*sigmaI + sigmaJ2*sigmaJ)/(sigmaI2 + sigmaJ2) : zero;
}

void AmoebaReferenceVdwForce::setEpsilonCombiningRule(const std::string& epsilonCombiningRule) {

    _epsilonCombiningRule = epsilonCombiningRule;
    std::transform(_epsilonCombiningRule.begin(), _epsilonCombiningRule.end(), _epsilonCombiningRule.begin(),  (int(*)(int)) std::toupper);

    // convert to upper case and set combining function

    if (_epsilonCombiningRule == "ARITHMETIC") {
         _combineEpsilons = &AmoebaReferenceVdwForce::arithmeticEpsilonCombiningRule;
    } else if (_epsilonCombiningRule == "HARMONIC") {
         _combineEpsilons = &AmoebaReferenceVdwForce::harmonicEpsilonCombiningRule;
    } else if (_epsilonCombiningRule == "HHG") {
         _combineEpsilons = &AmoebaReferenceVdwForce::hhgEpsilonCombiningRule;
    } else {
         _combineEpsilons = &AmoebaReferenceVdwForce::geometricEpsilonCombiningRule;
    }
}

std::string AmoebaReferenceVdwForce::getEpsilonCombiningRule() const {
    return _epsilonCombiningRule;
}

RealOpenMM AmoebaReferenceVdwForce::arithmeticEpsilonCombiningRule(RealOpenMM epsilonI, RealOpenMM epsilonJ) const {
    return 0.5*(epsilonI + epsilonJ);
}

RealOpenMM AmoebaReferenceVdwForce::geometricEpsilonCombiningRule(RealOpenMM epsilonI, RealOpenMM epsilonJ) const {
    return SQRT(epsilonI*epsilonJ);
}

RealOpenMM AmoebaReferenceVdwForce::harmonicEpsilonCombiningRule(RealOpenMM epsilonI, RealOpenMM epsilonJ) const {
    return (epsilonI != 0.0 && epsilonJ != 0.0) ? 2.0*(epsilonI*epsilonJ)/(epsilonI + epsilonJ) : 0.0;
}

RealOpenMM AmoebaReferenceVdwForce::hhgEpsilonCombiningRule(RealOpenMM epsilonI, RealOpenMM epsilonJ) const {
    RealOpenMM denominator = SQRT(epsilonI) + SQRT(epsilonJ);
    return (epsilonI != 0.0 && epsilonJ != 0.0) ? 4.0*(epsilonI*epsilonJ)/(denominator*denominator) : 0.0;
}

void AmoebaReferenceVdwForce::addReducedForce(unsigned int particleI, unsigned int particleIV,
                                              RealOpenMM reduction, RealOpenMM sign,
                                              Vec3& force, vector<RealVec>& forces) const {

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

RealOpenMM AmoebaReferenceVdwForce::calculatePairIxn(RealOpenMM combinedSigma, RealOpenMM combinedEpsilon,
                                                     const Vec3& particleIPosition,
                                                     const Vec3& particleJPosition,
                                                     Vec3& force) const {

   // ---------------------------------------------------------------------------------------

    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;
    static const RealOpenMM seven         = 7.0;

    static const RealOpenMM dhal          = 0.07;
    static const RealOpenMM ghal          = 0.12;

   // ---------------------------------------------------------------------------------------

    // get deltaR, R2, and R between 2 atoms

    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (_nonbondedMethod == CutoffPeriodic)
        ReferenceForce::getDeltaRPeriodic(particleJPosition, particleIPosition, _periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(particleJPosition, particleIPosition, deltaR);

    RealOpenMM r_ij_2       = deltaR[ReferenceForce::R2Index];
    RealOpenMM r_ij         = deltaR[ReferenceForce::RIndex];
    RealOpenMM sigma_7      = combinedSigma*combinedSigma*combinedSigma;
               sigma_7      = sigma_7*sigma_7*combinedSigma;

    RealOpenMM r_ij_6       = r_ij_2*r_ij_2*r_ij_2;
    RealOpenMM r_ij_7       = r_ij_6*r_ij;

    RealOpenMM rho          = r_ij_7 + ghal*sigma_7;

    RealOpenMM tau          = (dhal + one)/(r_ij + dhal*combinedSigma);
    RealOpenMM tau_7        = tau*tau*tau;
               tau_7        = tau_7*tau_7*tau;

    RealOpenMM dtau         = tau/(dhal + one);

    RealOpenMM ratio        = (sigma_7/rho);
    RealOpenMM gtau         = combinedEpsilon*tau_7*r_ij_6*(ghal+one)*ratio*ratio;

    RealOpenMM energy       = combinedEpsilon*tau_7*sigma_7*((ghal+one)*sigma_7/rho - two);

    RealOpenMM dEdR         = -seven*(dtau*energy + gtau);

    // tapering

    if ((_nonbondedMethod == CutoffNonPeriodic || _nonbondedMethod == CutoffPeriodic) && r_ij > _taperCutoff) {
        RealOpenMM delta    = r_ij - _taperCutoff;
        RealOpenMM taper    = 1.0 + delta*delta*delta*(_taperCoefficients[C3] + delta*(_taperCoefficients[C4] + delta*_taperCoefficients[C5]));
        RealOpenMM dtaper   = delta*delta*(3.0*_taperCoefficients[C3] + delta*(4.0*_taperCoefficients[C4] + delta*5.0*_taperCoefficients[C5]));
        dEdR                = energy*dtaper + dEdR*taper;
        energy             *= taper;
    }

    dEdR                   /= r_ij;

    force[0]                = dEdR*deltaR[0];
    force[1]                = dEdR*deltaR[1];
    force[2]                = dEdR*deltaR[2];

    return energy;

}

void AmoebaReferenceVdwForce::setReducedPositions(int numParticles,
                                                  const vector<RealVec>& particlePositions,
                                                  const std::vector<int>& indexIVs, 
                                                  const std::vector<RealOpenMM>& reductions,
                                                  std::vector<Vec3>& reducedPositions) const {

    static const RealOpenMM zero          = 0.0;

    reducedPositions.resize(numParticles);
    for (unsigned int ii = 0; ii <  static_cast<unsigned int>(numParticles); ii++) {
        if (reductions[ii] != zero) {
            int reductionIndex     = indexIVs[ii];
            reducedPositions[ii]   = Vec3(reductions[ii]*(particlePositions[ii][0] - particlePositions[reductionIndex][0]) + particlePositions[reductionIndex][0], 
                                          reductions[ii]*(particlePositions[ii][1] - particlePositions[reductionIndex][1]) + particlePositions[reductionIndex][1], 
                                          reductions[ii]*(particlePositions[ii][2] - particlePositions[reductionIndex][2]) + particlePositions[reductionIndex][2]); 
        } else {
            reducedPositions[ii]   = Vec3(particlePositions[ii][0], particlePositions[ii][1], particlePositions[ii][2]); 
        }
    }
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergy(int numParticles,
                                                            const vector<RealVec>& particlePositions,
                                                            const std::vector<int>& indexIVs, 
                                                            const std::vector<RealOpenMM>& sigmas,
                                                            const std::vector<RealOpenMM>& epsilons,
                                                            const std::vector<RealOpenMM>& reductions,
                                                            const std::vector< std::set<int> >& allExclusions,
                                                            vector<RealVec>& forces) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;

    // ---------------------------------------------------------------------------------------

    // set reduced coordinates

    std::vector<Vec3> reducedPositions;
    setReducedPositions(numParticles, particlePositions, indexIVs, reductions, reducedPositions);

    // loop over all particle pairs

    //    (1) initialize exclusion vector
    //    (2) calculate pair ixn, if not excluded
    //    (3) accumulate forces: if particle is a site where interaction position != particle position,
    //        then call addReducedForce() to apportion force to particle and its covalent partner
    //        based on reduction factor
    //    (4) reset exclusion vector

    RealOpenMM energy = zero;
    std::vector<unsigned int> exclusions(numParticles, 0);
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++) {
 
        RealOpenMM sigmaI      = sigmas[ii];
        RealOpenMM epsilonI    = epsilons[ii];
        for (std::set<int>::const_iterator jj = allExclusions[ii].begin(); jj != allExclusions[ii].end(); jj++) {
            exclusions[*jj] = 1;
        }

        for (unsigned int jj = ii+1; jj < static_cast<unsigned int>(numParticles); jj++) {
            if (exclusions[jj] == 0) {

                RealOpenMM combinedSigma   = (this->*_combineSigmas)(sigmaI, sigmas[jj]);
                RealOpenMM combinedEpsilon = (this->*_combineEpsilons)(epsilonI, epsilons[jj]);

                Vec3 force;
                energy                     += calculatePairIxn(combinedSigma, combinedEpsilon,
                                                               reducedPositions[ii], reducedPositions[jj],
                                                               force);
                
                if (indexIVs[ii] == ii) {
                    forces[ii][0] -= force[0];
                    forces[ii][1] -= force[1];
                    forces[ii][2] -= force[2];
                } else {
                    addReducedForce(ii, indexIVs[ii], reductions[ii], -one, force, forces);
                }
                if (indexIVs[jj] == jj) {
                    forces[jj][0] += force[0];
                    forces[jj][1] += force[1];
                    forces[jj][2] += force[2];
                } else {
                    addReducedForce(jj, indexIVs[jj], reductions[jj], one, force, forces);
                }

            }
        }

        for (std::set<int>::const_iterator jj = allExclusions[ii].begin(); jj != allExclusions[ii].end(); jj++) {
            exclusions[*jj] = 0;
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergy(int numParticles,
                                                            const vector<RealVec>& particlePositions,
                                                            const std::vector<int>& indexIVs, 
                                                            const std::vector<RealOpenMM>& sigmas,
                                                            const std::vector<RealOpenMM>& epsilons,
                                                            const std::vector<RealOpenMM>& reductions,
                                                            const NeighborList& neighborList,
                                                            vector<RealVec>& forces) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero          = 0.0;
    static const RealOpenMM one           = 1.0;
    static const RealOpenMM two           = 2.0;

    // ---------------------------------------------------------------------------------------

    // set reduced coordinates

    std::vector<Vec3> reducedPositions;
    setReducedPositions(numParticles, particlePositions, indexIVs, reductions, reducedPositions);
 
    // loop over neighbor list
    //    (1) calculate pair vdw ixn
    //    (2) accumulate forces: if particle is a site where interaction position != particle position,
    //        then call addReducedForce() to apportion force to particle and its covalent partner
    //        based on reduction factor

    RealOpenMM energy = zero;
    for (unsigned int ii = 0; ii < neighborList.size(); ii++) {

        OpenMM::AtomPair pair       = neighborList[ii];
        int siteI                   = pair.first;
        int siteJ                   = pair.second;

        RealOpenMM combinedSigma   = (this->*_combineSigmas)(sigmas[siteI], sigmas[siteJ]);
        RealOpenMM combinedEpsilon = (this->*_combineEpsilons)(epsilons[siteI], epsilons[siteJ]);

        Vec3 force;
        energy                     += calculatePairIxn(combinedSigma, combinedEpsilon,
                                                       reducedPositions[siteI], reducedPositions[siteJ], force);
                
        if (indexIVs[siteI] == siteI) {
            forces[siteI][0] -= force[0];
            forces[siteI][1] -= force[1];
            forces[siteI][2] -= force[2];
        } else {
            addReducedForce(siteI, indexIVs[siteI], reductions[siteI], -one, force, forces);
        }
        if (indexIVs[siteJ] == siteJ) {
            forces[siteJ][0] += force[0];
            forces[siteJ][1] += force[1];
            forces[siteJ][2] += force[2];
        } else {
            addReducedForce(siteJ, indexIVs[siteJ], reductions[siteJ], one, force, forces);
        }

    }

    return energy;
}
