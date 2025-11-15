/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceForce.h"
#include "ReferenceLCPOIxn.h"
#include "ReferenceNeighborList.h"
#include "SimTKOpenMMRealType.h"

using namespace OpenMM;
using namespace std;

ReferenceLCPOIxn::ReferenceLCPOIxn(const vector<int>& activeParticles, const vector<int>& activeParticlesInv, const vector<array<double, 4> >& parameters, double cutoff, bool usePeriodic) :
        activeParticles(activeParticles), activeParticlesInv(activeParticlesInv), parameters(parameters), cutoff(cutoff), usePeriodic(usePeriodic) {
    numParticles = activeParticlesInv.size();
    numActiveParticles = activeParticles.size();
}

double ReferenceLCPOIxn::execute(const Vec3* boxVectors, const std::vector<Vec3>& posData, std::vector<Vec3>& forceData, bool includeForces, bool includeEnergy) {
    if (!numActiveParticles) {
        return 0.0;
    }

    // We want a neighbor list in a form allowing us to query neighbors of each
    // particle.  Here, we use a NeighborList and then process its output:
    // inefficient, but simple, so suitable for this reference implementation.

    NeighborList neighborList;
    computeNeighborListVoxelHash(neighborList, numParticles, posData, vector<set<int> >(numParticles), boxVectors, usePeriodic, cutoff, 0.0, false);
    vector<map<int, NeighborInfo> > neighbors(numActiveParticles);
    for (auto atomPair : neighborList) {
        // Only include particles participating in the LCPO interaction.

        int i = activeParticlesInv[atomPair.first];
        int j = activeParticlesInv[atomPair.second];
        if (i == -1 || j == -1) {
            continue;
        }

        double iRadius = parameters[i][RadiusIndex];
        double jRadius = parameters[j][RadiusIndex];

        // Only include particles close enough to each other.

        double deltaR[ReferenceForce::LastDeltaRIndex];
        if (usePeriodic) {
            ReferenceForce::getDeltaRPeriodic(posData[atomPair.first], posData[atomPair.second], boxVectors, deltaR);
        }
        else {
            ReferenceForce::getDeltaR(posData[atomPair.first], posData[atomPair.second], deltaR);
        }
        double r = deltaR[ReferenceForce::RIndex];
        if (r >= iRadius + jRadius) {
            continue;
        }
        double rRecip = 1.0 / r;

        // Precompute and store the buried areas and their derivatives.

        double iRadiusPi = PI_M * iRadius;
        double jRadiusPi = PI_M * jRadius;
        double deltaRadiusR = (iRadius * iRadius - jRadius * jRadius) * rRecip;
        double deltaRadiusRSq = deltaRadiusR * rRecip;
        Vec3 direction = Vec3(deltaR[ReferenceForce::XIndex], deltaR[ReferenceForce::YIndex], deltaR[ReferenceForce::ZIndex]) * rRecip;
        neighbors[i][j] = {iRadiusPi * (2.0 * iRadius - r - deltaRadiusR), iRadiusPi * (deltaRadiusRSq - 1.0) * direction};
        neighbors[j][i] = {jRadiusPi * (2.0 * jRadius - r + deltaRadiusR), jRadiusPi * (deltaRadiusRSq + 1.0) * direction};
    }

    double energy = 0.0;

    // Compute LCPO two- and three-body energy and forces.

    for (int i = 0; i < numActiveParticles; i++) {
        double p2 = parameters[i][P2Index];
        double p3 = parameters[i][P3Index];
        double p4 = parameters[i][P4Index];
        double term2 = 0.0;
        double term3 = 0.0;
        double term4 = 0.0;
        Vec3 iForce;

        for (auto jNeighbor : neighbors[i]) {
            int j = jNeighbor.first;
            double Aij = jNeighbor.second.Aij;
            Vec3 dAij = jNeighbor.second.dAij;
            Vec3 jForce;

            // Two-body term.

            term2 += Aij;
            if (includeForces) {
                Vec3 ijForce2Body = p2 * dAij;
                iForce += ijForce2Body;
                jForce -= ijForce2Body;
            }

            // Three-body term: includes all pairs (j, k) of neighbors of i that
            // are also neighbors of each other.

            for (auto kNeighbor : neighbors[j]) {
                int k = kNeighbor.first;
                double Ajk = kNeighbor.second.Aij;
                Vec3 dAjk = kNeighbor.second.dAij;

                if (neighbors[i].find(k) == neighbors[i].end()) {
                    continue;
                }

                term3 += Ajk;
                term4 += Aij * Ajk;
                if (includeForces) {
                    Vec3 jkForce3Body = (p3 + p4 * Aij) * dAjk;
                    Vec3 ijForce3Body = p4 * dAij * Ajk;
                    iForce += ijForce3Body;
                    jForce += jkForce3Body - ijForce3Body;
                    forceData[activeParticles[k]] -= jkForce3Body;
                }
            }

            if (includeForces) {
                forceData[activeParticles[j]] += jForce;
            }
        }

        energy += p2 * term2 + p3 * term3 + p4 * term4;
        if (includeForces) {
            forceData[activeParticles[i]] += iForce;
        }
    }

    return energy;
}
