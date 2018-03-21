
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

AmoebaReferenceVdwForce::AmoebaReferenceVdwForce()
    : _nonbondedMethod(NoCutoff)
    , _functionalForm(BUFFERED_14_7)
    , _coupleMethod(Decouple)
    , _cutoff(1.0e+10)
    , _taperCutoffFactor(0.9) {
    setTaperCoefficients(_cutoff);
}

AmoebaReferenceVdwForce::NonbondedMethod AmoebaReferenceVdwForce::getNonbondedMethod() const {
    return _nonbondedMethod;
}

void AmoebaReferenceVdwForce::setNonbondedMethod(AmoebaReferenceVdwForce::NonbondedMethod nonbondedMethod) {
    _nonbondedMethod = nonbondedMethod;
}

void AmoebaReferenceVdwForce::setTaperCoefficients(double cutoff) {
    _taperCutoff = cutoff * _taperCutoffFactor;
    if (_taperCutoff != cutoff) {
        _taperCoefficients[C3] = 10.0 / pow(_taperCutoff - cutoff, 3.0);
        _taperCoefficients[C4] = 15.0 / pow(_taperCutoff - cutoff, 4.0);
        _taperCoefficients[C5] = 6.0 / pow(_taperCutoff - cutoff, 5.0);
    }
    else {
        _taperCoefficients[C3] = 0.0;
        _taperCoefficients[C4] = 0.0;
        _taperCoefficients[C5] = 0.0;
    }
}

void AmoebaReferenceVdwForce::setCutoff(double cutoff) {
    _cutoff = cutoff;
    setTaperCoefficients(_cutoff);
}

double AmoebaReferenceVdwForce::getCutoff() const {
    return _cutoff;
}

int AmoebaReferenceVdwForce::getCoupleMethod() const {
    return _coupleMethod;
}

void AmoebaReferenceVdwForce::setCoupleMethod(int method) {
    _coupleMethod = Decouple;
    if (method != Decouple) {
        _coupleMethod = Annihilate;
    }
}

void AmoebaReferenceVdwForce::setFunctionalForm(FunctionalForm functionalForm) {
    _functionalForm = functionalForm;
}

AmoebaReferenceVdwForce::FunctionalForm AmoebaReferenceVdwForce::getFunctionalForm() const {
    return _functionalForm;
}

void AmoebaReferenceVdwForce::setPeriodicBox(OpenMM::RealVec* vectors) {
    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
}

void AmoebaReferenceVdwForce::addReducedForce(unsigned int particleI, unsigned int particleIV,
    RealOpenMM reduction, RealOpenMM sign,
    Vec3& force, vector<RealVec>& forces) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM one = 1.0;

    // ---------------------------------------------------------------------------------------

    forces[particleI][0] += sign * force[0] * reduction;
    forces[particleI][1] += sign * force[1] * reduction;
    forces[particleI][2] += sign * force[2] * reduction;

    forces[particleIV][0] += sign * force[0] * (one - reduction);
    forces[particleIV][1] += sign * force[1] * (one - reduction);
    forces[particleIV][2] += sign * force[2] * (one - reduction);
}

RealOpenMM AmoebaReferenceVdwForce::calculatePairIxn(RealOpenMM combinedSigma,
    RealOpenMM combinedEpsilon,
    RealOpenMM combinedLambda,
    const Vec3& particleIPosition,
    const Vec3& particleJPosition,
    Vec3& force) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM one = 1.0;
    static const RealOpenMM two = 2.0;
    static const RealOpenMM seven = 7.0;

    static const RealOpenMM dhal = 0.07;
    static const RealOpenMM ghal = 0.12;

    // ---------------------------------------------------------------------------------------

    // get deltaR, R2, and R between 2 atoms

    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (_nonbondedMethod == CutoffPeriodic)
        ReferenceForce::getDeltaRPeriodic(particleJPosition, particleIPosition, _periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(particleJPosition, particleIPosition, deltaR);

    RealOpenMM energy = 0.0;
    RealOpenMM dEdR = 0.0;
    RealOpenMM r_ij_2 = deltaR[ReferenceForce::R2Index];
    RealOpenMM r_ij = deltaR[ReferenceForce::RIndex];
    if (_functionalForm == BUFFERED_14_7) {
        RealOpenMM comblambda2 = combinedLambda * combinedLambda;
        combinedEpsilon = combinedEpsilon * comblambda2 * comblambda2 * combinedLambda;
        RealOpenMM rho = r_ij / combinedSigma;
        RealOpenMM rho2 = rho * rho;
        RealOpenMM rho6 = rho2 * rho2 * rho2;
        RealOpenMM rhoplus = rho + 0.07;
        RealOpenMM rhodec2 = rhoplus * rhoplus;
        RealOpenMM rhodec = rhodec2 * rhodec2 * rhodec2;
        RealOpenMM scal = 0.7 * (1 - combinedLambda) * (1 - combinedLambda);
        RealOpenMM s1 = 1 / (scal + rhodec * rhoplus);
        RealOpenMM s2 = 1 / (scal + rho6 * rho + 0.12);
        RealOpenMM point72 = 1.07 * 1.07;
        RealOpenMM t1 = 1.07 * point72 * point72 * point72 * s1;
        RealOpenMM t2 = 1.12 * s2;
        RealOpenMM t2min = t2 - 2;
        RealOpenMM dt1 = -7.0 * rhodec * t1 * s1;
        RealOpenMM dt2 = -7.0 * rho6 * t2 * s2;
        energy = combinedEpsilon * t1 * (t2min);
        dEdR = combinedEpsilon * (dt1 * t2min + t1 * dt2) / combinedSigma;
    } else if (_functionalForm == LENNARD_JONES) {
        RealOpenMM pp1 = combinedSigma / r_ij;
        RealOpenMM pp2 = pp1 * pp1;
        RealOpenMM pp3 = pp2 * pp1;
        RealOpenMM pp6 = pp3 * pp3;
        RealOpenMM pp12 = pp6 * pp6;
        energy = combinedEpsilon * (pp12 - ((RealOpenMM)2.0) * pp6);
        dEdR = -combinedEpsilon * (pp12 - pp6) * ((RealOpenMM)12.0) / r_ij;
    }
    // tapering

    if ((_nonbondedMethod == CutoffNonPeriodic || _nonbondedMethod == CutoffPeriodic) && r_ij > _taperCutoff) {
        RealOpenMM delta = r_ij - _taperCutoff;
        RealOpenMM taper = 1.0 + delta * delta * delta * (_taperCoefficients[C3] + delta * (_taperCoefficients[C4] + delta * _taperCoefficients[C5]));
        RealOpenMM dtaper = delta * delta * (3.0 * _taperCoefficients[C3] + delta * (4.0 * _taperCoefficients[C4] + delta * 5.0 * _taperCoefficients[C5]));
        dEdR = energy * dtaper + dEdR * taper;
        energy *= taper;
    }

    dEdR /= r_ij;

    force[0] = dEdR * deltaR[0];
    force[1] = dEdR * deltaR[1];
    force[2] = dEdR * deltaR[2];

    return energy;
}

void AmoebaReferenceVdwForce::setReducedPositions(int numParticles,
    const vector<RealVec>& particlePositions,
    const std::vector<int>& indexIVs,
    const std::vector<RealOpenMM>& reductions,
    std::vector<Vec3>& reducedPositions) const {
    static const RealOpenMM zero = 0.0;

    reducedPositions.resize(numParticles);
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(numParticles); ii++) {
        if (reductions[ii] != zero) {
            int reductionIndex = indexIVs[ii];
            reducedPositions[ii] = Vec3(reductions[ii] * (particlePositions[ii][0] - particlePositions[reductionIndex][0]) + particlePositions[reductionIndex][0],
                reductions[ii] * (particlePositions[ii][1] - particlePositions[reductionIndex][1]) + particlePositions[reductionIndex][1],
                reductions[ii] * (particlePositions[ii][2] - particlePositions[reductionIndex][2]) + particlePositions[reductionIndex][2]);
        }
        else {
            reducedPositions[ii] = Vec3(particlePositions[ii][0], particlePositions[ii][1], particlePositions[ii][2]);
        }
    }
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergy(int numParticles, int numVdwprTypes,
    const std::vector<int>& vdwprTypes,
    const std::vector<OpenMM::RealVec>& particlePositions,
    const std::vector<int>& indexIVs,
    const std::vector<RealOpenMM>& combinedSigmas, const std::vector<RealOpenMM>& combinedEpsilons,
    const std::vector<RealOpenMM>& reductions, const std::vector<RealOpenMM>& lambdas,
    const std::vector<std::set<int> >& vdwExclusions,
    std::vector<OpenMM::RealVec>& forces) const {
    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero = 0.0;
    static const RealOpenMM one = 1.0;
    static const RealOpenMM two = 2.0;

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
        int typeI = vdwprTypes[ii];
        for (std::set<int>::const_iterator jj = vdwExclusions[ii].begin(); jj != vdwExclusions[ii].end(); jj++) {
            exclusions[*jj] = 1;
        }

        for (unsigned int jj = ii + 1; jj < static_cast<unsigned int>(numParticles); jj++) {
            if (exclusions[jj] == 0) {
                int typeJ = vdwprTypes[jj];
                int k = typeI * numVdwprTypes + typeJ;
                RealOpenMM combinedSigma = combinedSigmas[k];
                RealOpenMM combinedEpsilon = combinedEpsilons[k];
                RealOpenMM combinedLambda = 1.0;
                if (_coupleMethod == Decouple) {
                    combinedLambda = (lambdas[ii] == lambdas[jj] ? 1.0f : (lambdas[ii] < lambdas[jj] ? lambdas[ii] : lambdas[jj]));
                } else if (_coupleMethod == Annihilate) {
                    combinedLambda = (lambdas[ii] < lambdas[jj] ? lambdas[ii] : lambdas[jj]);
                }

                Vec3 force;

                energy += calculatePairIxn(combinedSigma, combinedEpsilon, combinedLambda,
                    reducedPositions[ii], reducedPositions[jj],
                    force);

                if (indexIVs[ii] == ii) {
                    forces[ii][0] -= force[0];
                    forces[ii][1] -= force[1];
                    forces[ii][2] -= force[2];
                }
                else {
                    addReducedForce(ii, indexIVs[ii], reductions[ii], -one, force, forces);
                }
                if (indexIVs[jj] == jj) {
                    forces[jj][0] += force[0];
                    forces[jj][1] += force[1];
                    forces[jj][2] += force[2];
                }
                else {
                    addReducedForce(jj, indexIVs[jj], reductions[jj], one, force, forces);
                }
            }
        }

        for (std::set<int>::const_iterator jj = vdwExclusions[ii].begin(); jj != vdwExclusions[ii].end(); jj++) {
            exclusions[*jj] = 0;
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceVdwForce::calculateForceAndEnergy(int numParticles, int numVdwprTypes,
    const std::vector<int>& vdwprTypes,
    const std::vector<OpenMM::RealVec>& particlePositions,
    const std::vector<int>& indexIVs,
    const std::vector<RealOpenMM>& combinedSigmas, const std::vector<RealOpenMM>& combinedEpsilons,
    const std::vector<RealOpenMM>& reductions, const std::vector<RealOpenMM>& lambdas,
    const NeighborList& neighborList,
    std::vector<OpenMM::RealVec>& forces) const {
    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero = 0.0;
    static const RealOpenMM one = 1.0;
    static const RealOpenMM two = 2.0;

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

        OpenMM::AtomPair pair = neighborList[ii];
        int siteI = pair.first;
        int siteJ = pair.second;
        int typeI = vdwprTypes[siteI];
        int typeJ = vdwprTypes[siteJ];
        int k = typeI * numVdwprTypes + typeJ;
        RealOpenMM combinedSigma = combinedSigmas[k];
        RealOpenMM combinedEpsilon = combinedEpsilons[k];
        RealOpenMM combinedLambda = 1.0f;
        if (_coupleMethod == Decouple) {
            combinedLambda = (lambdas[siteI] == lambdas[siteJ] ? 1.0f : (lambdas[siteI] < lambdas[siteJ] ? lambdas[siteI] : lambdas[siteJ]));
        } else if (_coupleMethod == Annihilate) {
            combinedLambda = (lambdas[siteI] == lambdas[siteJ] ? lambdas[siteI] : lambdas[siteJ]);
        }
        Vec3 force;
        energy += calculatePairIxn(combinedSigma, combinedEpsilon, combinedLambda,
            reducedPositions[siteI], reducedPositions[siteJ], force);

        if (indexIVs[siteI] == siteI) {
            forces[siteI][0] -= force[0];
            forces[siteI][1] -= force[1];
            forces[siteI][2] -= force[2];
        }
        else {
            addReducedForce(siteI, indexIVs[siteI], reductions[siteI], -one, force, forces);
        }
        if (indexIVs[siteJ] == siteJ) {
            forces[siteJ][0] += force[0];
            forces[siteJ][1] += force[1];
            forces[siteJ][2] += force[2];
        }
        else {
            addReducedForce(siteJ, indexIVs[siteJ], reductions[siteJ], one, force, forces);
        }
    }

    return energy;
}
