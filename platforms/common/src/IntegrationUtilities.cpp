/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/common/IntegrationUtilities.h"
#include "openmm/common/ComputeContext.h"
#include "CommonKernelSources.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/VirtualSite.h"
#include "quern.h"
#include "ReferenceCCMAAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>

using namespace OpenMM;
using namespace std;

struct IntegrationUtilities::ShakeCluster {
    int centralID;
    int peripheralID[3];
    int size;
    bool valid;
    double distance;
    double centralInvMass, peripheralInvMass;
    ShakeCluster() : valid(true) {
    }
    ShakeCluster(int centralID, double invMass) : centralID(centralID), centralInvMass(invMass), size(0), valid(true) {
    }
    void addAtom(int id, double dist, double invMass) {
        if (size == 3 || (size > 0 && abs(dist-distance)/distance > 1e-8) || (size > 0 && abs(invMass-peripheralInvMass)/peripheralInvMass > 1e-8))
            valid = false;
        else {
            peripheralID[size++] = id;
            distance = dist;
            peripheralInvMass = invMass;
        }
    }
    void markInvalid(map<int, ShakeCluster>& allClusters, vector<bool>& invalidForShake)
    {
        valid = false;
        invalidForShake[centralID] = true;
        for (int i = 0; i < size; i++) {
            invalidForShake[peripheralID[i]] = true;
            map<int, ShakeCluster>::iterator otherCluster = allClusters.find(peripheralID[i]);
            if (otherCluster != allClusters.end() && otherCluster->second.valid)
                otherCluster->second.markInvalid(allClusters, invalidForShake);
        }
    }
};

struct IntegrationUtilities::ConstraintOrderer : public binary_function<int, int, bool> {
    const vector<int>& atom1;
    const vector<int>& atom2;
    const vector<int>& constraints;
    ConstraintOrderer(const vector<int>& atom1, const vector<int>& atom2, const vector<int>& constraints) : atom1(atom1), atom2(atom2), constraints(constraints) {
    }
    bool operator()(int x, int y) {
        int ix = constraints[x];
        int iy = constraints[y];
        if (atom1[ix] != atom1[iy])
            return atom1[ix] < atom1[iy];
        return atom2[ix] < atom2[iy];
    }
};

IntegrationUtilities::IntegrationUtilities(ComputeContext& context, const System& system) : context(context),
        randomPos(0), hasOverlappingVsites(false) {
    // Create workspace arrays.

    lastStepSize = mm_double2(0.0, 0.0);
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
        posDelta.initialize<mm_double4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<mm_double4> deltas(posDelta.getSize(), mm_double4(0.0, 0.0, 0.0, 0.0));
        posDelta.upload(deltas);
        stepSize.initialize<mm_double2>(context, 1, "stepSize");
        stepSize.upload(&lastStepSize);
    }
    else {
        posDelta.initialize<mm_float4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<mm_float4> deltas(posDelta.getSize(), mm_float4(0.0f, 0.0f, 0.0f, 0.0f));
        posDelta.upload(deltas);
        stepSize.initialize<mm_float2>(context, 1, "stepSize");
        mm_float2 lastStepSizeFloat = mm_float2(0.0f, 0.0f);
        stepSize.upload(&lastStepSizeFloat);
    }

    // Record the set of constraints and how many constraints each atom is involved in.

    vector<int> atom1;
    vector<int> atom2;
    vector<double> distance;
    vector<int> constraintCount(context.getNumAtoms(), 0);
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int p1, p2;
        double d;
        system.getConstraintParameters(i, p1, p2, d);
        if (system.getParticleMass(p1) != 0 || system.getParticleMass(p2) != 0) {
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

    int numAtoms = system.getNumParticles();
    vector<map<int, float> > settleConstraints(numAtoms);
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

    vector<bool> isShakeAtom(numAtoms, false);
    if (settleClusters.size() > 0) {
        vector<mm_int4> atoms;
        vector<mm_float2> params;
        for (int i = 0; i < (int) settleClusters.size(); i++) {
            int atom1 = settleClusters[i];
            int atom2 = settleConstraints[atom1].begin()->first;
            int atom3 = (++settleConstraints[atom1].begin())->first;
            float dist12 = settleConstraints[atom1].find(atom2)->second;
            float dist13 = settleConstraints[atom1].find(atom3)->second;
            float dist23 = settleConstraints[atom2].find(atom3)->second;
            if (dist12 == dist13) {
                // atom1 is the central atom
                atoms.push_back(mm_int4(atom1, atom2, atom3, 0));
                params.push_back(mm_float2(dist12, dist23));
            }
            else if (dist12 == dist23) {
                // atom2 is the central atom
                atoms.push_back(mm_int4(atom2, atom1, atom3, 0));
                params.push_back(mm_float2(dist12, dist13));
            }
            else if (dist13 == dist23) {
                // atom3 is the central atom
                atoms.push_back(mm_int4(atom3, atom1, atom2, 0));
                params.push_back(mm_float2(dist13, dist12));
            }
            else
                continue; // We can't handle this with SETTLE
            isShakeAtom[atom1] = true;
            isShakeAtom[atom2] = true;
            isShakeAtom[atom3] = true;
        }
        if (atoms.size() > 0) {
            settleAtoms.initialize<mm_int4>(context, atoms.size(), "settleAtoms");
            settleParams.initialize<mm_float2>(context, params.size(), "settleParams");
            settleAtoms.upload(atoms);
            settleParams.upload(params);
        }
    }

    // Find clusters consisting of a central atom with up to three peripheral atoms.

    map<int, ShakeCluster> clusters;
    vector<bool> invalidForShake(numAtoms, false);
    for (int i = 0; i < (int) atom1.size(); i++) {
        if (isShakeAtom[atom1[i]])
            continue; // This is being taken care of with SETTLE.

        // Determine which is the central atom.

        bool firstIsCentral;
        if (constraintCount[atom1[i]] > 1)
            firstIsCentral = true;
        else if (constraintCount[atom2[i]] > 1)
            firstIsCentral = false;
        else if (atom1[i] < atom2[i])
            firstIsCentral = true;
        else
            firstIsCentral = false;
        int centralID, peripheralID;
        if (firstIsCentral) {
            centralID = atom1[i];
            peripheralID = atom2[i];
        }
        else {
            centralID = atom2[i];
            peripheralID = atom1[i];
        }

        // Add it to the cluster.

        if (clusters.find(centralID) == clusters.end()) {
            clusters[centralID] = ShakeCluster(centralID, 1.0/system.getParticleMass(centralID));
        }
        ShakeCluster& cluster = clusters[centralID];
        cluster.addAtom(peripheralID, distance[i], 1.0/system.getParticleMass(peripheralID));
        if (constraintCount[peripheralID] != 1 || invalidForShake[atom1[i]] || invalidForShake[atom2[i]]) {
            cluster.markInvalid(clusters, invalidForShake);
            map<int, ShakeCluster>::iterator otherCluster = clusters.find(peripheralID);
            if (otherCluster != clusters.end() && otherCluster->second.valid)
                otherCluster->second.markInvalid(clusters, invalidForShake);
        }
    }
    int validShakeClusters = 0;
    for (map<int, ShakeCluster>::iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
        ShakeCluster& cluster = iter->second;
        if (cluster.valid) {
            cluster.valid = !invalidForShake[cluster.centralID] && cluster.size == constraintCount[cluster.centralID];
            for (int i = 0; i < cluster.size; i++)
                if (invalidForShake[cluster.peripheralID[i]])
                    cluster.valid = false;
            if (cluster.valid)
                ++validShakeClusters;
        }
    }

    // Record the SHAKE clusters.

    if (validShakeClusters > 0) {
        vector<mm_int4> atoms;
        vector<mm_float4> params;
        int index = 0;
        for (map<int, ShakeCluster>::const_iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
            const ShakeCluster& cluster = iter->second;
            if (!cluster.valid)
                continue;
            atoms.push_back(mm_int4(cluster.centralID, cluster.peripheralID[0], (cluster.size > 1 ? cluster.peripheralID[1] : -1), (cluster.size > 2 ? cluster.peripheralID[2] : -1)));
            params.push_back(mm_float4((float) cluster.centralInvMass, (float) (0.5/(cluster.centralInvMass+cluster.peripheralInvMass)), (float) (cluster.distance*cluster.distance), (float) cluster.peripheralInvMass));
            isShakeAtom[cluster.centralID] = true;
            isShakeAtom[cluster.peripheralID[0]] = true;
            if (cluster.size > 1)
                isShakeAtom[cluster.peripheralID[1]] = true;
            if (cluster.size > 2)
                isShakeAtom[cluster.peripheralID[2]] = true;
            ++index;
        }
        shakeAtoms.initialize<mm_int4>(context, atoms.size(), "shakeAtoms");
        shakeParams.initialize<mm_float4>(context, params.size(), "shakeParams");
        shakeAtoms.upload(atoms);
        shakeParams.upload(params);
    }

    // Find connected constraints for CCMA.

    vector<int> ccmaConstraints;
    for (unsigned i = 0; i < atom1.size(); i++)
        if (!isShakeAtom[atom1[i]])
            ccmaConstraints.push_back(i);

    // Record the connections between constraints.

    int numCCMA = (int) ccmaConstraints.size();
    if (numCCMA > 0) {
        // Record information needed by ReferenceCCMAAlgorithm.
        
        vector<pair<int, int> > refIndices(numCCMA);
        vector<double> refDistance(numCCMA);
        for (int i = 0; i < numCCMA; i++) {
            int index = ccmaConstraints[i];
            refIndices[i] = make_pair(atom1[index], atom2[index]);
            refDistance[i] = distance[index];
        }
        vector<double> refMasses(numAtoms);
        for (int i = 0; i < numAtoms; ++i)
            refMasses[i] = system.getParticleMass(i);

        // Look up angles for CCMA.
        
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        for (int i = 0; i < system.getNumForces(); i++) {
            const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
            if (force != NULL) {
                for (int j = 0; j < force->getNumAngles(); j++) {
                    int atom1, atom2, atom3;
                    double angle, k;
                    force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                    angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, angle));
                }
            }
        }
        
        // Create a ReferenceCCMAAlgorithm.  It will build and invert the constraint matrix for us.
        
        ReferenceCCMAAlgorithm ccma(numAtoms, numCCMA, refIndices, refDistance, refMasses, angles, 0.1);
        vector<vector<pair<int, double> > > matrix = ccma.getMatrix();
        int maxRowElements = 0;
        for (unsigned i = 0; i < matrix.size(); i++)
            maxRowElements = max(maxRowElements, (int) matrix[i].size());
        maxRowElements++;

        // Build the list of constraints for each atom.

        vector<vector<int> > atomConstraints(context.getNumAtoms());
        for (int i = 0; i < numCCMA; i++) {
            atomConstraints[atom1[ccmaConstraints[i]]].push_back(i);
            atomConstraints[atom2[ccmaConstraints[i]]].push_back(i);
        }
        int maxAtomConstraints = 0;
        for (unsigned i = 0; i < atomConstraints.size(); i++)
            maxAtomConstraints = max(maxAtomConstraints, (int) atomConstraints[i].size());

        // Sort the constraints.

        vector<int> constraintOrder(numCCMA);
        for (int i = 0; i < numCCMA; ++i)
            constraintOrder[i] = i;
        sort(constraintOrder.begin(), constraintOrder.end(), ConstraintOrderer(atom1, atom2, ccmaConstraints));
        vector<int> inverseOrder(numCCMA);
        for (int i = 0; i < numCCMA; ++i)
            inverseOrder[constraintOrder[i]] = i;
        for (int i = 0; i < (int)matrix.size(); ++i)
            for (int j = 0; j < (int)matrix[i].size(); ++j)
                matrix[i][j].first = inverseOrder[matrix[i][j].first];

        // Record the CCMA data structures.

        ccmaAtoms.initialize<mm_int2>(context, numCCMA, "CcmaAtoms");
        ccmaAtomConstraints.initialize<int>(context, numAtoms*maxAtomConstraints, "CcmaAtomConstraints");
        ccmaNumAtomConstraints.initialize<int>(context, numAtoms, "CcmaAtomConstraintsIndex");
        ccmaConstraintMatrixColumn.initialize<int>(context, numCCMA*maxRowElements, "ConstraintMatrixColumn");
        ccmaConverged.initialize<int>(context, 2, "ccmaConverged");
        vector<mm_int2> atomsVec(ccmaAtoms.getSize());
        vector<int> atomConstraintsVec(ccmaAtomConstraints.getSize());
        vector<int> numAtomConstraintsVec(ccmaNumAtomConstraints.getSize());
        vector<int> constraintMatrixColumnVec(ccmaConstraintMatrixColumn.getSize());
        int elementSize = (context.getUseDoublePrecision() || context.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
        ccmaDistance.initialize(context, numCCMA, 4*elementSize, "CcmaDistance");
        ccmaDelta1.initialize(context, numCCMA, elementSize, "CcmaDelta1");
        ccmaDelta2.initialize(context, numCCMA, elementSize, "CcmaDelta2");
        ccmaReducedMass.initialize(context, numCCMA, elementSize, "CcmaReducedMass");
        ccmaConstraintMatrixValue.initialize(context, numCCMA*maxRowElements, elementSize, "ConstraintMatrixValue");
        vector<mm_double4> distanceVec(ccmaDistance.getSize());
        vector<double> reducedMassVec(ccmaReducedMass.getSize());
        vector<double> constraintMatrixValueVec(ccmaConstraintMatrixValue.getSize());
        for (int i = 0; i < numCCMA; i++) {
            int index = constraintOrder[i];
            int c = ccmaConstraints[index];
            atomsVec[i].x = atom1[c];
            atomsVec[i].y = atom2[c];
            distanceVec[i].w = distance[c];
            reducedMassVec[i] = (0.5/(1.0/system.getParticleMass(atom1[c])+1.0/system.getParticleMass(atom2[c])));
            for (unsigned int j = 0; j < matrix[index].size(); j++) {
                constraintMatrixColumnVec[i+j*numCCMA] = matrix[index][j].first;
                constraintMatrixValueVec[i+j*numCCMA] = matrix[index][j].second;
            }
            constraintMatrixColumnVec[i+matrix[index].size()*numCCMA] = numCCMA;
        }
        ccmaDistance.upload(distanceVec, true);
        ccmaReducedMass.upload(reducedMassVec, true);
        ccmaConstraintMatrixValue.upload(constraintMatrixValueVec, true);
        for (unsigned int i = 0; i < atomConstraints.size(); i++) {
            numAtomConstraintsVec[i] = atomConstraints[i].size();
            for (unsigned int j = 0; j < atomConstraints[i].size(); j++) {
                bool forward = (atom1[ccmaConstraints[atomConstraints[i][j]]] == i);
                atomConstraintsVec[i+j*numAtoms] = (forward ? inverseOrder[atomConstraints[i][j]]+1 : -inverseOrder[atomConstraints[i][j]]-1);
            }
        }
        ccmaAtoms.upload(atomsVec);
        ccmaAtomConstraints.upload(atomConstraintsVec);
        ccmaNumAtomConstraints.upload(numAtomConstraintsVec);
        ccmaConstraintMatrixColumn.upload(constraintMatrixColumnVec);
    }
    
    // Build the list of virtual sites.
    
    vector<mm_int4> vsite2AvgAtomVec;
    vector<mm_double2> vsite2AvgWeightVec;
    vector<mm_int4> vsite3AvgAtomVec;
    vector<mm_double4> vsite3AvgWeightVec;
    vector<mm_int4> vsiteOutOfPlaneAtomVec;
    vector<mm_double4> vsiteOutOfPlaneWeightVec;
    vector<int> vsiteLocalCoordsIndexVec;
    vector<int> vsiteLocalCoordsAtomVec;
    vector<int> vsiteLocalCoordsStartVec;
    vector<double> vsiteLocalCoordsWeightVec;
    vector<mm_double4> vsiteLocalCoordsPosVec;
    for (int i = 0; i < numAtoms; i++) {
        if (system.isVirtualSite(i)) {
            if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                vsite2AvgAtomVec.push_back(mm_int4(i, site.getParticle(0), site.getParticle(1), 0));
                vsite2AvgWeightVec.push_back(mm_double2(site.getWeight(0), site.getWeight(1)));
            }
            else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                vsite3AvgAtomVec.push_back(mm_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsite3AvgWeightVec.push_back(mm_double4(site.getWeight(0), site.getWeight(1), site.getWeight(2), 0.0));
            }
            else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                vsiteOutOfPlaneAtomVec.push_back(mm_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsiteOutOfPlaneWeightVec.push_back(mm_double4(site.getWeight12(), site.getWeight13(), site.getWeightCross(), 0.0));
            }
            else if (dynamic_cast<const LocalCoordinatesSite*>(&system.getVirtualSite(i)) != NULL) {
                // A local coordinates site.
                
                const LocalCoordinatesSite& site = dynamic_cast<const LocalCoordinatesSite&>(system.getVirtualSite(i));
                int numParticles = site.getNumParticles();
                vector<double> origin, x, y;
                site.getOriginWeights(origin);
                site.getXWeights(x);
                site.getYWeights(y);
                vsiteLocalCoordsIndexVec.push_back(i);
                vsiteLocalCoordsStartVec.push_back(vsiteLocalCoordsAtomVec.size());
                for (int j = 0; j < numParticles; j++) {
                    vsiteLocalCoordsAtomVec.push_back(site.getParticle(j));
                    vsiteLocalCoordsWeightVec.push_back(origin[j]);
                    vsiteLocalCoordsWeightVec.push_back(x[j]);
                    vsiteLocalCoordsWeightVec.push_back(y[j]);
                }
                Vec3 pos = site.getLocalPosition();
                vsiteLocalCoordsPosVec.push_back(mm_double4(pos[0], pos[1], pos[2], 0.0));
            }
        }
    }
    vsiteLocalCoordsStartVec.push_back(vsiteLocalCoordsAtomVec.size());
    int num2Avg = vsite2AvgAtomVec.size();
    int num3Avg = vsite3AvgAtomVec.size();
    int numOutOfPlane = vsiteOutOfPlaneAtomVec.size();
    int numLocalCoords = vsiteLocalCoordsPosVec.size();
    numVsites = num2Avg+num3Avg+numOutOfPlane+numLocalCoords;
    vsite2AvgAtoms.initialize<mm_int4>(context, max(1, num2Avg), "vsite2AvgAtoms");
    vsite3AvgAtoms.initialize<mm_int4>(context, max(1, num3Avg), "vsite3AvgAtoms");
    vsiteOutOfPlaneAtoms.initialize<mm_int4>(context, max(1, numOutOfPlane), "vsiteOutOfPlaneAtoms");
    vsiteLocalCoordsIndex.initialize<int>(context, max(1, (int) vsiteLocalCoordsIndexVec.size()), "vsiteLocalCoordsIndex");
    vsiteLocalCoordsAtoms.initialize<int>(context, max(1, (int) vsiteLocalCoordsAtomVec.size()), "vsiteLocalCoordsAtoms");
    vsiteLocalCoordsStartIndex.initialize<int>(context, max(1, (int) vsiteLocalCoordsStartVec.size()), "vsiteLocalCoordsStartIndex");
    if (num2Avg > 0)
        vsite2AvgAtoms.upload(vsite2AvgAtomVec);
    if (num3Avg > 0)
        vsite3AvgAtoms.upload(vsite3AvgAtomVec);
    if (numOutOfPlane > 0)
        vsiteOutOfPlaneAtoms.upload(vsiteOutOfPlaneAtomVec);
    if (numLocalCoords > 0) {
        vsiteLocalCoordsIndex.upload(vsiteLocalCoordsIndexVec);
        vsiteLocalCoordsAtoms.upload(vsiteLocalCoordsAtomVec);
        vsiteLocalCoordsStartIndex.upload(vsiteLocalCoordsStartVec);
    }
    int elementSize = (context.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    vsite2AvgWeights.initialize(context, max(1, num2Avg), 2*elementSize, "vsite2AvgWeights");
    vsite3AvgWeights.initialize(context, max(1, num3Avg), 4*elementSize, "vsite3AvgWeights");
    vsiteOutOfPlaneWeights.initialize(context, max(1, numOutOfPlane), 4*elementSize, "vsiteOutOfPlaneWeights");
    vsiteLocalCoordsWeights.initialize(context, max(1, (int) vsiteLocalCoordsWeightVec.size()), elementSize, "vsiteLocalCoordsWeights");
    vsiteLocalCoordsPos.initialize(context, max(1, (int) vsiteLocalCoordsPosVec.size()), 4*elementSize, "vsiteLocalCoordsPos");
    if (num2Avg > 0)
        vsite2AvgWeights.upload(vsite2AvgWeightVec, true);
    if (num3Avg > 0)
        vsite3AvgWeights.upload(vsite3AvgWeightVec, true);
    if (numOutOfPlane > 0)
        vsiteOutOfPlaneWeights.upload(vsiteOutOfPlaneWeightVec, true);
    if (numLocalCoords > 0) {
        vsiteLocalCoordsWeights.upload(vsiteLocalCoordsWeightVec, true);
        vsiteLocalCoordsPos.upload(vsiteLocalCoordsPosVec, true);
    }

    // If multiple virtual sites depend on the same particle, make sure the force distribution
    // can be done safely.
    
    vector<int> atomCounts(numAtoms, 0);
    for (int i = 0; i < numAtoms; i++)
        if (system.isVirtualSite(i))
            for (int j = 0; j < system.getVirtualSite(i).getNumParticles(); j++)
                atomCounts[system.getVirtualSite(i).getParticle(j)]++;
    for (int i = 0; i < numAtoms; i++)
        if (atomCounts[i] > 1)
            hasOverlappingVsites = true;
    if (hasOverlappingVsites && !context.getSupports64BitGlobalAtomics())
        throw OpenMMException("This device does not support 64 bit atomics.  Cannot have multiple virtual sites that depend on the same atom.");

    // Create the kernels used by this class.

    map<string, string> defines;
    defines["NUM_CCMA_CONSTRAINTS"] = context.intToString(numCCMA);
    defines["NUM_ATOMS"] = context.intToString(numAtoms);
    defines["NUM_2_AVERAGE"] = context.intToString(num2Avg);
    defines["NUM_3_AVERAGE"] = context.intToString(num3Avg);
    defines["NUM_OUT_OF_PLANE"] = context.intToString(numOutOfPlane);
    defines["NUM_LOCAL_COORDS"] = context.intToString(numLocalCoords);
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    if (hasOverlappingVsites)
        defines["HAS_OVERLAPPING_VSITES"] = "1";
    ComputeProgram program = context.compileProgram(CommonKernelSources::integrationUtilities, defines);
    settlePosKernel = program->createKernel("applySettleToPositions");
    settleVelKernel = program->createKernel("applySettleToVelocities");
    shakePosKernel = program->createKernel("applyShakeToPositions");
    shakeVelKernel = program->createKernel("applyShakeToVelocities");
    ccmaDirectionsKernel = program->createKernel("computeCCMAConstraintDirections");
    ccmaPosForceKernel = program->createKernel("computeCCMAPositionConstraintForce");
    ccmaVelForceKernel = program->createKernel("computeCCMAVelocityConstraintForce");
    ccmaMultiplyKernel = program->createKernel("multiplyByCCMAConstraintMatrix");
    ccmaUpdateKernel = program->createKernel("updateCCMAAtomPositions");
    vsitePositionKernel = program->createKernel("computeVirtualSites");
    vsiteForceKernel = program->createKernel("distributeVirtualSiteForces");
    vsiteSaveForcesKernel = program->createKernel("saveDistributedForces");
    randomKernel = program->createKernel("generateRandomNumbers");
    timeShiftKernel = program->createKernel("timeShiftVelocities");

    // Set arguments for virtual site kernels.

    vsitePositionKernel->addArg(context.getPosq());
    if (context.getUseMixedPrecision())
        vsitePositionKernel->addArg(context.getPosqCorrection());
    else
        vsitePositionKernel->addArg(NULL);
    vsitePositionKernel->addArg(vsite2AvgAtoms);
    vsitePositionKernel->addArg(vsite2AvgWeights);
    vsitePositionKernel->addArg(vsite3AvgAtoms);
    vsitePositionKernel->addArg(vsite3AvgWeights);
    vsitePositionKernel->addArg(vsiteOutOfPlaneAtoms);
    vsitePositionKernel->addArg(vsiteOutOfPlaneWeights);
    vsitePositionKernel->addArg(vsiteLocalCoordsIndex);
    vsitePositionKernel->addArg(vsiteLocalCoordsAtoms);
    vsitePositionKernel->addArg(vsiteLocalCoordsWeights);
    vsitePositionKernel->addArg(vsiteLocalCoordsPos);
    vsitePositionKernel->addArg(vsiteLocalCoordsStartIndex);
    vsiteForceKernel->addArg(context.getPosq());
    if (context.getUseMixedPrecision())
        vsiteForceKernel->addArg(context.getPosqCorrection());
    else
        vsiteForceKernel->addArg(NULL);
    vsiteForceKernel->addArg(); // Skip argument 2: the force array hasn't been created yet.
    vsiteForceKernel->addArg(vsite2AvgAtoms);
    vsiteForceKernel->addArg(vsite2AvgWeights);
    vsiteForceKernel->addArg(vsite3AvgAtoms);
    vsiteForceKernel->addArg(vsite3AvgWeights);
    vsiteForceKernel->addArg(vsiteOutOfPlaneAtoms);
    vsiteForceKernel->addArg(vsiteOutOfPlaneWeights);
    vsiteForceKernel->addArg(vsiteLocalCoordsIndex);
    vsiteForceKernel->addArg(vsiteLocalCoordsAtoms);
    vsiteForceKernel->addArg(vsiteLocalCoordsWeights);
    vsiteForceKernel->addArg(vsiteLocalCoordsPos);
    vsiteForceKernel->addArg(vsiteLocalCoordsStartIndex);
    for (int i = 0; i < 3; i++)
        vsiteSaveForcesKernel->addArg();

    // Set arguments for constraint kernels.

    if (settleAtoms.isInitialized()) {
        settlePosKernel->addArg(settleAtoms.getSize());
        settlePosKernel->addArg();
        settlePosKernel->addArg(context.getPosq());
        settlePosKernel->addArg(posDelta);
        settlePosKernel->addArg(context.getVelm());
        settlePosKernel->addArg(settleAtoms);
        settlePosKernel->addArg(settleParams);
        if (context.getUseMixedPrecision())
            settlePosKernel->addArg(context.getPosqCorrection());
        settleVelKernel->addArg(settleAtoms.getSize());
        settleVelKernel->addArg();
        settleVelKernel->addArg(context.getPosq());
        settleVelKernel->addArg(posDelta);
        settleVelKernel->addArg(context.getVelm());
        settleVelKernel->addArg(settleAtoms);
        settleVelKernel->addArg(settleParams);
        if (context.getUseMixedPrecision())
            settleVelKernel->addArg(context.getPosqCorrection());
    }
    if (shakeAtoms.isInitialized()) {
        shakePosKernel->addArg(shakeAtoms.getSize());
        shakePosKernel->addArg();
        shakePosKernel->addArg(context.getPosq());
        shakePosKernel->addArg(posDelta);
        shakePosKernel->addArg(shakeAtoms);
        shakePosKernel->addArg(shakeParams);
        if (context.getUseMixedPrecision())
            shakePosKernel->addArg(context.getPosqCorrection());
        shakeVelKernel->addArg(shakeAtoms.getSize());
        shakeVelKernel->addArg();
        shakeVelKernel->addArg(context.getPosq());
        shakeVelKernel->addArg(context.getVelm());
        shakeVelKernel->addArg(shakeAtoms);
        shakeVelKernel->addArg(shakeParams);
        if (context.getUseMixedPrecision())
            shakeVelKernel->addArg(context.getPosqCorrection());
    }
    if (ccmaAtoms.isInitialized()) {
        ccmaDirectionsKernel->addArg(ccmaAtoms);
        ccmaDirectionsKernel->addArg(ccmaDistance);
        ccmaDirectionsKernel->addArg(context.getPosq());
        ccmaDirectionsKernel->addArg(ccmaConverged);
        if (context.getUseMixedPrecision())
            ccmaDirectionsKernel->addArg(context.getPosqCorrection());
        ccmaPosForceKernel->addArg(ccmaAtoms);
        ccmaPosForceKernel->addArg(ccmaDistance);
        ccmaPosForceKernel->addArg(posDelta);
        ccmaPosForceKernel->addArg(ccmaReducedMass);
        ccmaPosForceKernel->addArg(ccmaDelta1);
        ccmaPosForceKernel->addArg(ccmaConverged);
        ccmaPosForceKernel->addArg();
        ccmaPosForceKernel->addArg();
        ccmaPosForceKernel->addArg();
        ccmaVelForceKernel->addArg(ccmaAtoms);
        ccmaVelForceKernel->addArg(ccmaDistance);
        ccmaVelForceKernel->addArg(context.getVelm());
        ccmaVelForceKernel->addArg(ccmaReducedMass);
        ccmaVelForceKernel->addArg(ccmaDelta1);
        ccmaVelForceKernel->addArg(ccmaConverged);
        ccmaVelForceKernel->addArg();
        ccmaVelForceKernel->addArg();
        ccmaVelForceKernel->addArg();
        ccmaMultiplyKernel->addArg(ccmaDelta1);
        ccmaMultiplyKernel->addArg(ccmaDelta2);
        ccmaMultiplyKernel->addArg(ccmaConstraintMatrixColumn);
        ccmaMultiplyKernel->addArg(ccmaConstraintMatrixValue);
        ccmaMultiplyKernel->addArg(ccmaConverged);
        ccmaMultiplyKernel->addArg();
        ccmaUpdateKernel->addArg(ccmaNumAtomConstraints);
        ccmaUpdateKernel->addArg(ccmaAtomConstraints);
        ccmaUpdateKernel->addArg(ccmaDistance);
        ccmaUpdateKernel->addArg();
        ccmaUpdateKernel->addArg(context.getVelm());
        ccmaUpdateKernel->addArg(ccmaDelta1);
        ccmaUpdateKernel->addArg(ccmaDelta2);
        ccmaUpdateKernel->addArg(ccmaConverged);
        ccmaUpdateKernel->addArg();
    }

    // Arguments for time shift kernel will be set later.
    
    for (int i = 0; i < 3; i++)
        timeShiftKernel->addArg();
}

void IntegrationUtilities::setNextStepSize(double size) {
    if (size != lastStepSize.x || size != lastStepSize.y) {
        lastStepSize = mm_double2(size, size);
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            stepSize.upload(&lastStepSize);
        else {
            mm_float2 lastStepSizeFloat = mm_float2((float) size, (float) size);
            stepSize.upload(&lastStepSizeFloat);
        }
    }
}

double IntegrationUtilities::getLastStepSize() {
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
        stepSize.download(&lastStepSize);
    else {
        mm_float2 lastStepSizeFloat;
        stepSize.download(&lastStepSizeFloat);
        lastStepSize = mm_double2(lastStepSizeFloat.x, lastStepSizeFloat.y);
    }
    return lastStepSize.y;
}

void IntegrationUtilities::applyConstraints(double tol) {
    applyConstraintsImpl(false, tol);
}

void IntegrationUtilities::applyVelocityConstraints(double tol) {
    applyConstraintsImpl(true, tol);
}

void IntegrationUtilities::computeVirtualSites() {
    if (numVsites > 0)
        vsitePositionKernel->execute(numVsites);
}

void IntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random.isInitialized()) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("IntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random.initialize<mm_float4>(context, 4*context.getPaddedNumAtoms(), "random");
    randomSeed.initialize<mm_int4>(context, context.getNumThreadBlocks()*64, "randomSeed");
    randomPos = random.getSize();
    randomKernel->addArg(random.getSize());
    randomKernel->addArg(random);
    randomKernel->addArg(randomSeed);

    // Use a quick and dirty RNG to pick seeds for the real random number generator.

    vector<mm_int4> seed(randomSeed.getSize());
    unsigned int r = randomNumberSeed;
    if (r == 0)
        r = (unsigned int) osrngseed(); // A seed of 0 means use a unique one
    for (int i = 0; i < randomSeed.getSize(); i++) {
        seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
    }
    randomSeed.upload(seed);
}

int IntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random.getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    if (numValues > random.getSize()) {
        random.resize(numValues);
        randomKernel->setArg(0, numValues);
    }
    randomKernel->execute(random.getSize(), 64);
    randomPos = numValues;
    return 0;
}

void IntegrationUtilities::createCheckpoint(ostream& stream) {
    int numChains = noseHooverChainState.size();
    bool useDouble = context.getUseDoublePrecision() || context.getUseMixedPrecision();
    stream.write((char*) &numChains, sizeof(int));
    for (auto &chainState: noseHooverChainState){
        int chainID = chainState.first;
        int chainLength = chainState.second.getSize();
        stream.write((char*) &chainID, sizeof(int));
        stream.write((char*) &chainLength, sizeof(int));
        if (useDouble) {
            vector<mm_double2> stateVec;
            chainState.second.download(stateVec);
            stream.write((char*) stateVec.data(), sizeof(mm_double2)*chainLength);
        }
        else {
            vector<mm_float2> stateVec;
            chainState.second.download(stateVec);
            stream.write((char*) stateVec.data(), sizeof(mm_float2)*chainLength);
        }
    }
    if (!random.isInitialized())
        return;
    stream.write((char*) &randomPos, sizeof(int));
    vector<mm_float4> randomVec;
    random.download(randomVec);
    stream.write((char*) &randomVec[0], sizeof(mm_float4)*random.getSize());
    vector<mm_int4> randomSeedVec;
    randomSeed.download(randomSeedVec);
    stream.write((char*) &randomSeedVec[0], sizeof(mm_int4)*randomSeed.getSize());
}

void IntegrationUtilities::loadCheckpoint(istream& stream) {
    int numChains;
    bool useDouble = context.getUseDoublePrecision() || context.getUseMixedPrecision();
    stream.read((char*) &numChains, sizeof(int));
    noseHooverChainState.clear();
    for (int i = 0; i < numChains; i++) {
        int chainID, chainLength;
        stream.read((char*) &chainID, sizeof(int));
        stream.read((char*) &chainLength, sizeof(int));
        if (useDouble) {
            noseHooverChainState[chainID] = ComputeArray();
            noseHooverChainState[chainID].initialize<mm_double2>(context, chainLength, "chainState" + to_string(chainID));
            vector<mm_double2> stateVec(chainLength);
            stream.read((char*) &stateVec[0], sizeof(mm_double2)*chainLength);
            noseHooverChainState[chainID].upload(stateVec);
        }
        else {
            noseHooverChainState[chainID] = ComputeArray();
            noseHooverChainState[chainID].initialize<mm_float2>(context, chainLength, "chainState" + to_string(chainID));
            vector<mm_float2> stateVec(chainLength);
            stream.read((char*) &stateVec[0], sizeof(mm_float2)*chainLength);
            noseHooverChainState[chainID].upload(stateVec);
        }
    }
    if (!random.isInitialized())
        return;
    stream.read((char*) &randomPos, sizeof(int));
    vector<mm_float4> randomVec(random.getSize());
    stream.read((char*) &randomVec[0], sizeof(mm_float4)*random.getSize());
    random.upload(randomVec);
    vector<mm_int4> randomSeedVec(randomSeed.getSize());
    stream.read((char*) &randomSeedVec[0], sizeof(mm_int4)*randomSeed.getSize());
    randomSeed.upload(randomSeedVec);
}

double IntegrationUtilities::computeKineticEnergy(double timeShift) {
    int numParticles = context.getNumAtoms();
    if (timeShift != 0) {
        // Copy the velocities into the posDelta array while we temporarily modify them.

        context.getVelm().copyTo(posDelta);

        // Apply the time shift.

        timeShiftKernel->setArg(0, context.getVelm());
        timeShiftKernel->setArg(1, context.getLongForceBuffer());
        if (context.getUseDoublePrecision())
            timeShiftKernel->setArg(2, timeShift);
        else
            timeShiftKernel->setArg(2, (float) timeShift);
        timeShiftKernel->execute(numParticles);
        applyConstraintsImpl(true, 1e-4);
    }
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
        vector<mm_double4> velm;
        context.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++) {
            mm_double4 v = velm[i];
            if (v.w != 0)
                energy += (v.x*v.x+v.y*v.y+v.z*v.z)/v.w;
        }
    }
    else {
        vector<mm_float4> velm;
        context.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++) {
            mm_float4 v = velm[i];
            if (v.w != 0)
                energy += (v.x*v.x+v.y*v.y+v.z*v.z)/v.w;
        }
    }
    
    // Restore the velocities.
    
    if (timeShift != 0)
        posDelta.copyTo(context.getVelm());
    return 0.5*energy;
}
