/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "ReferenceConstraints.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/OpenMMException.h"
#include <map>
#include <utility>
#include <vector>

using namespace OpenMM;
using namespace std;

ReferenceConstraints::ReferenceConstraints(const System& system, RealOpenMM tolerance) : ccma(NULL), settle(NULL) {
    this->tolerance = tolerance;
    int numParticles = system.getNumParticles();
    vector<RealOpenMM> masses(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = (RealOpenMM) system.getParticleMass(i);

    // Record the set of constraints and how many constraints each atom is involved in.

    vector<int> atom1;
    vector<int> atom2;
    vector<double> distance;
    vector<int> constraintCount(numParticles, 0);
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int p1, p2;
        double d;
        system.getConstraintParameters(i, p1, p2, d);
        if (masses[p1] != 0 || masses[p2] != 0) {
            atom1.push_back(p1);
            atom2.push_back(p2);
            distance.push_back(d);
            constraintCount[p1]++;
            constraintCount[p2]++;
        }
    }

    // Identify clusters of three atoms that can be treated with SETTLE.  First, for every
    // atom that might be part of such a cluster, make a list of the two other atoms it is
    // connected to.

    vector<map<int, float> > settleConstraints(numParticles);
    for (int i = 0; i < (int)atom1.size(); i++) {
        if (constraintCount[atom1[i]] == 2 && constraintCount[atom2[i]] == 2) {
            settleConstraints[atom1[i]][atom2[i]] = (float) distance[i];
            settleConstraints[atom2[i]][atom1[i]] = (float) distance[i];
        }
    }

    // Now remove the ones that don't actually form closed loops of three atoms.

    vector<int> settleClusters;
    for (int i = 0; i < (int)settleConstraints.size(); i++) {
        if (settleConstraints[i].size() == 2) {
            int partner1 = settleConstraints[i].begin()->first;
            int partner2 = (++settleConstraints[i].begin())->first;
            if (settleConstraints[partner1].size() != 2 || settleConstraints[partner2].size() != 2 ||
                    settleConstraints[partner1].find(partner2) == settleConstraints[partner1].end())
                settleConstraints[i].clear();
            else if (i < partner1 && i < partner2)
                settleClusters.push_back(i);
        }
        else
            settleConstraints[i].clear();
    }

    // Record the SETTLE clusters.

    vector<bool> isSettleAtom(numParticles, false);
    int numSETTLE = settleClusters.size();
    if (numSETTLE > 0) {
        vector<int> atom1(numSETTLE);
        vector<int> atom2(numSETTLE);
        vector<int> atom3(numSETTLE);
        vector<RealOpenMM> distance1(numSETTLE);
        vector<RealOpenMM> distance2(numSETTLE);
        for (int i = 0; i < numSETTLE; i++) {
            int p1 = settleClusters[i];
            int p2 = settleConstraints[p1].begin()->first;
            int p3 = (++settleConstraints[p1].begin())->first;
            float dist12 = settleConstraints[p1].find(p2)->second;
            float dist13 = settleConstraints[p1].find(p3)->second;
            float dist23 = settleConstraints[p2].find(p3)->second;
            if (dist12 == dist13) {
                // p1 is the central atom
                atom1[i] = p1;
                atom2[i] = p2;
                atom3[i] = p3;
                distance1[i] = dist12;
                distance2[i] = dist23;
            }
            else if (dist12 == dist23) {
                // p2 is the central atom
                atom1[i] = p2;
                atom2[i] = p1;
                atom3[i] = p3;
                distance1[i] = dist12;
                distance2[i] = dist13;
            }
            else if (dist13 == dist23) {
                // p3 is the central atom
                atom1[i] = p3;
                atom2[i] = p1;
                atom3[i] = p2;
                distance1[i] = dist13;
                distance2[i] = dist12;
            }
            else
                throw OpenMMException("Two of the three distances constrained with SETTLE must be the same.");
            isSettleAtom[p1] = true;
            isSettleAtom[p2] = true;
            isSettleAtom[p3] = true;
        }
        settle = new ReferenceSETTLEAlgorithm(atom1, atom2, atom3, distance1, distance2, masses, tolerance);
    }

    // All other constraints are handled with CCMA.

    vector<int> ccmaConstraints;
    for (unsigned i = 0; i < atom1.size(); i++)
        if (!isSettleAtom[atom1[i]])
            ccmaConstraints.push_back(i);
    int numCCMA = ccmaConstraints.size();
    if (numCCMA > 0) {
        // Record particles and distances for CCMA.
        
        vector<pair<int, int> > ccmaIndices(numCCMA);
        vector<RealOpenMM> ccmaDistance(numCCMA);
        for (int i = 0; i < numCCMA; i++) {
            int index = ccmaConstraints[i];
            ccmaIndices[i] = make_pair(atom1[index], atom2[index]);
            ccmaDistance[i] = distance[index];
        }

        // Look up angles for CCMA.
        
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        for (int i = 0; i < system.getNumForces(); i++) {
            const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
            if (force != NULL) {
                for (int j = 0; j < force->getNumAngles(); j++) {
                    int atom1, atom2, atom3;
                    double angle, k;
                    force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                    angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, (RealOpenMM) angle));
                }
            }
        }
        
        // Create the CCMA object.
        
        ccma = new ReferenceCCMAAlgorithm(numParticles, numCCMA, ccmaIndices, ccmaDistance, masses, angles, tolerance);
    }
}

ReferenceConstraints::~ReferenceConstraints() {
    if (ccma != NULL)
        delete ccma;
    if (settle != NULL)
        delete settle;
}

RealOpenMM ReferenceConstraints::getTolerance() const {
    return tolerance;
}

void ReferenceConstraints::setTolerance(RealOpenMM tolerance) {
    this->tolerance = tolerance;
    if (ccma != NULL)
        ccma->setTolerance(tolerance);
    if (settle != NULL)
        settle->setTolerance(tolerance);
}

int ReferenceConstraints::apply(int numberOfAtoms, vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& atomCoordinatesP, vector<RealOpenMM>& inverseMasses) {
    if (ccma != NULL)
        ccma->apply(numberOfAtoms, atomCoordinates, atomCoordinatesP, inverseMasses);
    if (settle != NULL)
        settle->apply(numberOfAtoms, atomCoordinates, atomCoordinatesP, inverseMasses);
    return SimTKOpenMMCommon::DefaultReturn;
}

int ReferenceConstraints::applyToVelocities(int numberOfAtoms, vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& velocities, vector<RealOpenMM>& inverseMasses) {
    if (ccma != NULL)
        ccma->applyToVelocities(numberOfAtoms, atomCoordinates, velocities, inverseMasses);
    if (settle != NULL)
        settle->applyToVelocities(numberOfAtoms, atomCoordinates, velocities, inverseMasses);
    return SimTKOpenMMCommon::DefaultReturn;
}
