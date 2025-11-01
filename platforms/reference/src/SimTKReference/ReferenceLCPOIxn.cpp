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

ReferenceLCPOIxn::ReferenceLCPOIxn(const vector<int>& indices, const vector<int>& particles, const vector<array<double, 4> >& parameters, double cutoff, bool usePeriodic) :
        indices(indices), particles(particles), parameters(parameters), cutoff(cutoff), usePeriodic(usePeriodic) {
    numParticles = indices.size();
    numActiveParticles = particles.size();
}

double ReferenceLCPOIxn::execute(const Vec3* boxVectors, const std::vector<Vec3>& posData, std::vector<Vec3>& forceData, bool includeForces, bool includeEnergy) {
    // We want a neighbor list in a form allowing us to query neighbors of each
    // particle.  Here, we use a NeighborList and then process its output: very
    // inefficient, but simple, so suitable for this reference implementation.

    NeighborList neighborList;
    computeNeighborListVoxelHash(neighborList, numParticles, posData, vector<set<int> >(numParticles), boxVectors, usePeriodic, cutoff, 0.0, false);
    vector<map<int, double> > neighbors(numActiveParticles);
    for (auto pair : neighborList) {
        // Only include particles participating in the LCPO interaction.

        int i = indices[pair.first];
        int j = indices[pair.second];
        if (i == -1 || j == -1) {
            continue;
        }

        double iRadius = parameters[i][RadiusIndex];
        double jRadius = parameters[j][RadiusIndex];

        // Only include particles close enough to each other.

        double deltaR[ReferenceForce::LastDeltaRIndex];
        if (usePeriodic) {
            ReferenceForce::getDeltaRPeriodic(posData[pair.first], posData[pair.second], boxVectors, deltaR);
        }
        else {
            ReferenceForce::getDeltaR(posData[pair.first], posData[pair.second], deltaR);
        }
        double r = deltaR[ReferenceForce::RIndex];
        if (r >= iRadius + jRadius) {
            continue;
        }

        // Precompute and store the buried areas.

        neighbors[i][j] = 2.0 * PI_M * iRadius * (iRadius - r / 2.0 - (iRadius * iRadius - jRadius * jRadius) / (2.0 * r));
        neighbors[j][i] = 2.0 * PI_M * jRadius * (jRadius - r / 2.0 - (jRadius * jRadius - iRadius * iRadius) / (2.0 * r));
    }

    double energy = 0.0;

    // Compute LCPO two- and three-body energy and forces.

    for (int i = 0; i < numActiveParticles; i++) {
        double term2 = 0.0;
        double term3 = 0.0;
        double term4 = 0.0;

        for (auto jNeighbor : neighbors[i]) {
            int j = jNeighbor.first;
            double Aij = jNeighbor.second;

            // Two-body term.

            term2 += Aij;

            // Three-body term: includes all pairs (j, k) of neighbors of i that
            // are also neighbors of each other.

            for (auto kNeighbor : neighbors[j]) {
                int k = kNeighbor.first;
                double Ajk = kNeighbor.second;

                if (neighbors[i].find(k) == neighbors[i].end()) {
                    continue;
                }

                term3 += Ajk;
                term4 += Aij * Ajk;
            }
        }

        energy += parameters[i][P2Index] * term2 + parameters[i][P3Index] * term3 + parameters[i][P4Index] * term4;
    }

    return energy;
}