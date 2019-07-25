
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
#include "AmoebaReferenceCTForce.h"
#include "ReferenceForce.h"
#include <algorithm>
#include <cctype>

using std::vector;
using namespace OpenMM;

AmoebaReferenceCTForce::AmoebaReferenceCTForce()
    : _nonbondedMethod(NoCutoff)
    , _cutoff(1.0e+10)
    , _taperCutoffFactor(0.9) {
    setTaperCoefficients(_cutoff);
}

AmoebaReferenceCTForce::NonbondedMethod AmoebaReferenceCTForce::getNonbondedMethod() const {
    return _nonbondedMethod;
}

void AmoebaReferenceCTForce::setNonbondedMethod(AmoebaReferenceCTForce::NonbondedMethod nonbondedMethod) {
    _nonbondedMethod = nonbondedMethod;
}

void AmoebaReferenceCTForce::setTaperCoefficients(double cutoff) {
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

void AmoebaReferenceCTForce::setCutoff(double cutoff) {
    _cutoff = cutoff;
    setTaperCoefficients(_cutoff);
}

double AmoebaReferenceCTForce::getCutoff() const {
    return _cutoff;
}


void AmoebaReferenceCTForce::setPeriodicBox(OpenMM::RealVec* vectors) {
    _periodicBoxVectors[0] = vectors[0];
    _periodicBoxVectors[1] = vectors[1];
    _periodicBoxVectors[2] = vectors[2];
}


RealOpenMM AmoebaReferenceCTForce::calculatePairIxn(RealOpenMM combinedApre,
    RealOpenMM combinedBexp,
    const Vec3& particleIPosition,
    const Vec3& particleJPosition,
    Vec3& force) const {

    // get deltaR, R2, and R between 2 atoms

    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (_nonbondedMethod == CutoffPeriodic)
        ReferenceForce::getDeltaRPeriodic(particleJPosition, particleIPosition, _periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(particleJPosition, particleIPosition, deltaR);

    RealOpenMM energy = 0.0;
    RealOpenMM dEdR = 0.0;
    RealOpenMM r_ij = deltaR[ReferenceForce::RIndex];

    energy = -combinedApre*1000.0 *EXP(-combinedBexp*r_ij); 
    dEdR   =  combinedBexp*combinedApre*1000.0 *EXP(-combinedBexp*r_ij); 

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

RealOpenMM AmoebaReferenceCTForce::calculateForceAndEnergy(int numParticles, int numCTprTypes,
    const std::vector<int>& CTprTypes,
    const std::vector<OpenMM::RealVec>& particlePositions,
    const std::vector<RealOpenMM>& combinedApres, const std::vector<RealOpenMM>& combinedBexps,
    const std::vector<std::set<int> >& CTExclusions,
    std::vector<OpenMM::RealVec>& forces) const {

    // ---------------------------------------------------------------------------------------
    static const RealOpenMM zero = 0.0;
    static const RealOpenMM one = 1.0;
    static const RealOpenMM two = 2.0;

    // ---------------------------------------------------------------------------------------

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
        int typeI = CTprTypes[ii];
        for (std::set<int>::const_iterator jj = CTExclusions[ii].begin(); jj != CTExclusions[ii].end(); jj++) {
            exclusions[*jj] = 1;
        }

        for (unsigned int jj = ii + 1; jj < static_cast<unsigned int>(numParticles); jj++) {
            if (exclusions[jj] == 0) {
                int typeJ = CTprTypes[jj];
                int k = typeI * numCTprTypes + typeJ;
                RealOpenMM combinedApre = combinedApres[k];
                RealOpenMM combinedBexp = combinedBexps[k];

                Vec3 force;

                energy += calculatePairIxn(combinedApre, combinedBexp, 
                    particlePositions[ii], particlePositions[jj],
                    force);

                forces[ii][0] -= force[0];
                forces[ii][1] -= force[1];
                forces[ii][2] -= force[2];

                forces[jj][0] += force[0];
                forces[jj][1] += force[1];
                forces[jj][2] += force[2];

            }
        }

        for (std::set<int>::const_iterator jj = CTExclusions[ii].begin(); jj != CTExclusions[ii].end(); jj++) {
            exclusions[*jj] = 0;
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceCTForce::calculateForceAndEnergy(int numParticles, int numCTprTypes,
    const std::vector<int>& CTprTypes,
    const std::vector<OpenMM::RealVec>& particlePositions,
    const std::vector<RealOpenMM>& combinedApres, const std::vector<RealOpenMM>& combinedBexps,
    const NeighborList& neighborList,
    std::vector<OpenMM::RealVec>& forces) const {
    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero = 0.0;
    static const RealOpenMM one = 1.0;
    static const RealOpenMM two = 2.0;

    // ---------------------------------------------------------------------------------------

    // loop over neighbor list
    //    (1) calculate pair CT ixn
    //    (2) accumulate forces: if particle is a site where interaction position != particle position,
    //        then call addReducedForce() to apportion force to particle and its covalent partner
    //        based on reduction factor

    RealOpenMM energy = zero;
    for (unsigned int ii = 0; ii < neighborList.size(); ii++) {

        OpenMM::AtomPair pair = neighborList[ii];
        int siteI = pair.first;
        int siteJ = pair.second;
        int typeI = CTprTypes[siteI];
        int typeJ = CTprTypes[siteJ];
        int k = typeI * numCTprTypes + typeJ;
        RealOpenMM combinedApre = combinedApres[k];
        RealOpenMM combinedBexp = combinedBexps[k];
        Vec3 force;
        energy += calculatePairIxn(combinedApre, combinedBexp, 
            particlePositions[siteI], particlePositions[siteJ], force);

        forces[siteI][0] -= force[0];
        forces[siteI][1] -= force[1];
        forces[siteI][2] -= force[2];

        forces[siteJ][0] += force[0];
        forces[siteJ][1] += force[1];
        forces[siteJ][2] += force[2];

    }

    return energy;
}
