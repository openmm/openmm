/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2018 Stanford University and the Authors.      *
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

#include "CudaIntegrationUtilities.h"
#include "CudaArray.h"
#include "CudaKernelSources.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/VirtualSite.h"
#include "quern.h"
#include "CudaExpressionUtilities.h"
#include "ReferenceCCMAAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) CHECK_RESULT2(result, errorMessage);
#define CHECK_RESULT2(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<context.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

struct CudaIntegrationUtilities::ShakeCluster {
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

struct CudaIntegrationUtilities::ConstraintOrderer : public binary_function<int, int, bool> {
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

CudaIntegrationUtilities::CudaIntegrationUtilities(CudaContext& context, const System& system) : context(context),
        randomPos(0), ccmaConvergedMemory(NULL) {
    // Create workspace arrays.

    lastStepSize = make_double2(0.0, 0.0);
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
        posDelta.initialize<double4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<double4> deltas(posDelta.getSize(), make_double4(0.0, 0.0, 0.0, 0.0));
        posDelta.upload(deltas);
        stepSize.initialize<double2>(context, 1, "stepSize");
        stepSize.upload(&lastStepSize);
    }
    else {
        posDelta.initialize<float4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<float4> deltas(posDelta.getSize(), make_float4(0.0f, 0.0f, 0.0f, 0.0f));
        posDelta.upload(deltas);
        stepSize.initialize<float2>(context, 1, "stepSize");
        float2 lastStepSizeFloat = make_float2(0.0f, 0.0f);
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
        vector<int4> atoms;
        vector<float2> params;
        for (int i = 0; i < (int) settleClusters.size(); i++) {
            int atom1 = settleClusters[i];
            int atom2 = settleConstraints[atom1].begin()->first;
            int atom3 = (++settleConstraints[atom1].begin())->first;
            float dist12 = settleConstraints[atom1].find(atom2)->second;
            float dist13 = settleConstraints[atom1].find(atom3)->second;
            float dist23 = settleConstraints[atom2].find(atom3)->second;
            if (dist12 == dist13) {
                // atom1 is the central atom
                atoms.push_back(make_int4(atom1, atom2, atom3, 0));
                params.push_back(make_float2(dist12, dist23));
            }
            else if (dist12 == dist23) {
                // atom2 is the central atom
                atoms.push_back(make_int4(atom2, atom1, atom3, 0));
                params.push_back(make_float2(dist12, dist13));
            }
            else if (dist13 == dist23) {
                // atom3 is the central atom
                atoms.push_back(make_int4(atom3, atom1, atom2, 0));
                params.push_back(make_float2(dist13, dist12));
            }
            else
                continue; // We can't handle this with SETTLE
            isShakeAtom[atom1] = true;
            isShakeAtom[atom2] = true;
            isShakeAtom[atom3] = true;
        }
        if (atoms.size() > 0) {
            settleAtoms.initialize<int4>(context, atoms.size(), "settleAtoms");
            settleParams.initialize<float2>(context, params.size(), "settleParams");
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
        vector<int4> atoms;
        vector<float4> params;
        int index = 0;
        for (map<int, ShakeCluster>::const_iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
            const ShakeCluster& cluster = iter->second;
            if (!cluster.valid)
                continue;
            atoms.push_back(make_int4(cluster.centralID, cluster.peripheralID[0], (cluster.size > 1 ? cluster.peripheralID[1] : -1), (cluster.size > 2 ? cluster.peripheralID[2] : -1)));
            params.push_back(make_float4((float) cluster.centralInvMass, (float) (0.5/(cluster.centralInvMass+cluster.peripheralInvMass)), (float) (cluster.distance*cluster.distance), (float) cluster.peripheralInvMass));
            isShakeAtom[cluster.centralID] = true;
            isShakeAtom[cluster.peripheralID[0]] = true;
            if (cluster.size > 1)
                isShakeAtom[cluster.peripheralID[1]] = true;
            if (cluster.size > 2)
                isShakeAtom[cluster.peripheralID[2]] = true;
            ++index;
        }
        shakeAtoms.initialize<int4>(context, atoms.size(), "shakeAtoms");
        shakeParams.initialize<float4>(context, params.size(), "shakeParams");
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

        ccmaAtoms.initialize<int2>(context, numCCMA, "CcmaAtoms");
        ccmaAtomConstraints.initialize<int>(context, numAtoms*maxAtomConstraints, "CcmaAtomConstraints");
        ccmaNumAtomConstraints.initialize<int>(context, numAtoms, "CcmaAtomConstraintsIndex");
        ccmaConstraintMatrixColumn.initialize<int>(context, numCCMA*maxRowElements, "ConstraintMatrixColumn");
        ccmaConverged.initialize<int>(context, 2, "ccmaConverged");
        CHECK_RESULT2(cuMemHostAlloc((void**) &ccmaConvergedMemory, sizeof(int), CU_MEMHOSTALLOC_DEVICEMAP), "Error allocating pinned memory");
        CHECK_RESULT2(cuMemHostGetDevicePointer(&ccmaConvergedDeviceMemory, ccmaConvergedMemory, 0), "Error getting device address for pinned memory");
        vector<int2> atomsVec(ccmaAtoms.getSize());
        vector<int> atomConstraintsVec(ccmaAtomConstraints.getSize());
        vector<int> numAtomConstraintsVec(ccmaNumAtomConstraints.getSize());
        vector<int> constraintMatrixColumnVec(ccmaConstraintMatrixColumn.getSize());
        int elementSize = (context.getUseDoublePrecision() || context.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
        ccmaDistance.initialize(context, numCCMA, 4*elementSize, "CcmaDistance");
        ccmaDelta1.initialize(context, numCCMA, elementSize, "CcmaDelta1");
        ccmaDelta2.initialize(context, numCCMA, elementSize, "CcmaDelta2");
        ccmaReducedMass.initialize(context, numCCMA, elementSize, "CcmaReducedMass");
        ccmaConstraintMatrixValue.initialize(context, numCCMA*maxRowElements, elementSize, "ConstraintMatrixValue");
        vector<double4> distanceVec(ccmaDistance.getSize());
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
    
    vector<int4> vsite2AvgAtomVec;
    vector<double2> vsite2AvgWeightVec;
    vector<int4> vsite3AvgAtomVec;
    vector<double4> vsite3AvgWeightVec;
    vector<int4> vsiteOutOfPlaneAtomVec;
    vector<double4> vsiteOutOfPlaneWeightVec;
    vector<int> vsiteLocalCoordsIndexVec;
    vector<int> vsiteLocalCoordsAtomVec;
    vector<int> vsiteLocalCoordsStartVec;
    vector<double> vsiteLocalCoordsWeightVec;
    vector<double4> vsiteLocalCoordsPosVec;
    for (int i = 0; i < numAtoms; i++) {
        if (system.isVirtualSite(i)) {
            if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                vsite2AvgAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), 0));
                vsite2AvgWeightVec.push_back(make_double2(site.getWeight(0), site.getWeight(1)));
            }
            else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                vsite3AvgAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsite3AvgWeightVec.push_back(make_double4(site.getWeight(0), site.getWeight(1), site.getWeight(2), 0.0));
            }
            else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                vsiteOutOfPlaneAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsiteOutOfPlaneWeightVec.push_back(make_double4(site.getWeight12(), site.getWeight13(), site.getWeightCross(), 0.0));
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
                vsiteLocalCoordsPosVec.push_back(make_double4(pos[0], pos[1], pos[2], 0.0));
            }
        }
    }
    vsiteLocalCoordsStartVec.push_back(vsiteLocalCoordsAtomVec.size());
    int num2Avg = vsite2AvgAtomVec.size();
    int num3Avg = vsite3AvgAtomVec.size();
    int numOutOfPlane = vsiteOutOfPlaneAtomVec.size();
    int numLocalCoords = vsiteLocalCoordsPosVec.size();
    vsite2AvgAtoms.initialize<int4>(context, max(1, num2Avg), "vsite2AvgAtoms");
    vsite3AvgAtoms.initialize<int4>(context, max(1, num3Avg), "vsite3AvgAtoms");
    vsiteOutOfPlaneAtoms.initialize<int4>(context, max(1, numOutOfPlane), "vsiteOutOfPlaneAtoms");
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

    // Create the kernels used by this class.

    map<string, string> defines;
    defines["NUM_CCMA_CONSTRAINTS"] = context.intToString(numCCMA);
    defines["NUM_ATOMS"] = context.intToString(numAtoms);
    defines["NUM_2_AVERAGE"] = context.intToString(num2Avg);
    defines["NUM_3_AVERAGE"] = context.intToString(num3Avg);
    defines["NUM_OUT_OF_PLANE"] = context.intToString(numOutOfPlane);
    defines["NUM_LOCAL_COORDS"] = context.intToString(numLocalCoords);
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    CUmodule module = context.createModule(CudaKernelSources::vectorOps+CudaKernelSources::integrationUtilities, defines);
    settlePosKernel = context.getKernel(module, "applySettleToPositions");
    settleVelKernel = context.getKernel(module, "applySettleToVelocities");
    shakePosKernel = context.getKernel(module, "applyShakeToPositions");
    shakeVelKernel = context.getKernel(module, "applyShakeToVelocities");
    ccmaDirectionsKernel = context.getKernel(module, "computeCCMAConstraintDirections");
    ccmaPosForceKernel = context.getKernel(module, "computeCCMAPositionConstraintForce");
    ccmaVelForceKernel = context.getKernel(module, "computeCCMAVelocityConstraintForce");
    ccmaMultiplyKernel = context.getKernel(module, "multiplyByCCMAConstraintMatrix");
    ccmaUpdateKernel = context.getKernel(module, "updateCCMAAtomPositions");
    CHECK_RESULT2(cuEventCreate(&ccmaEvent, CU_EVENT_DISABLE_TIMING), "Error creating event for CCMA");
    vsitePositionKernel = context.getKernel(module, "computeVirtualSites");
    vsiteForceKernel = context.getKernel(module, "distributeVirtualSiteForces");
    numVsites = num2Avg+num3Avg+numOutOfPlane+numLocalCoords;
    randomKernel = context.getKernel(module, "generateRandomNumbers");
    timeShiftKernel = context.getKernel(module, "timeShiftVelocities");
}

CudaIntegrationUtilities::~CudaIntegrationUtilities() {
    context.setAsCurrent();
    if (ccmaConvergedMemory != NULL)
        cuMemFreeHost(ccmaConvergedMemory);
}

void CudaIntegrationUtilities::setNextStepSize(double size) {
    if (size != lastStepSize.x || size != lastStepSize.y) {
        lastStepSize = make_double2(size, size);
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            stepSize.upload(&lastStepSize);
        else {
            float2 lastStepSizeFloat = make_float2((float) size, (float) size);
            stepSize.upload(&lastStepSizeFloat);
        }
    }
}

double CudaIntegrationUtilities::getLastStepSize() {
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
        stepSize.download(&lastStepSize);
    else {
        float2 lastStepSizeFloat;
        stepSize.download(&lastStepSizeFloat);
        lastStepSize = make_double2(lastStepSizeFloat.x, lastStepSizeFloat.y);
    }
    return lastStepSize.y;
}

void CudaIntegrationUtilities::applyConstraints(double tol) {
    applyConstraints(false, tol);
}

void CudaIntegrationUtilities::applyVelocityConstraints(double tol) {
    applyConstraints(true, tol);
}

void CudaIntegrationUtilities::applyConstraints(bool constrainVelocities, double tol) {
    CUfunction settleKernel, shakeKernel, ccmaForceKernel;
    if (constrainVelocities) {
        settleKernel = settleVelKernel;
        shakeKernel = shakeVelKernel;
        ccmaForceKernel = ccmaVelForceKernel;
    }
    else {
        settleKernel = settlePosKernel;
        shakeKernel = shakePosKernel;
        ccmaForceKernel = ccmaPosForceKernel;
    }
    float floatTol = (float) tol;
    void* tolPointer = (context.getUseDoublePrecision() || context.getUseMixedPrecision() ? (void*) &tol : (void*) &floatTol);
    CUdeviceptr posCorrection = (context.getUseMixedPrecision() ? context.getPosqCorrection().getDevicePointer() : 0);
    if (settleAtoms.isInitialized()) {
        int numClusters = settleAtoms.getSize();
        void* args[] = {&numClusters, tolPointer, &context.getPosq().getDevicePointer(), &posCorrection,
                &posDelta.getDevicePointer(), &context.getVelm().getDevicePointer(),
                &settleAtoms.getDevicePointer(), &settleParams.getDevicePointer()};
        context.executeKernel(settleKernel, args, settleAtoms.getSize());
    }
    if (shakeAtoms.isInitialized()) {
        int numClusters = shakeAtoms.getSize();
        void* args[] = {&numClusters, tolPointer, &context.getPosq().getDevicePointer(), &posCorrection,
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta.getDevicePointer(),
                &shakeAtoms.getDevicePointer(), &shakeParams.getDevicePointer()};
        context.executeKernel(shakeKernel, args, shakeAtoms.getSize());
    }
    if (ccmaAtoms.isInitialized()) {
        void* directionsArgs[] = {&ccmaAtoms.getDevicePointer(), &ccmaDistance.getDevicePointer(), &context.getPosq().getDevicePointer(), &posCorrection, &ccmaConverged.getDevicePointer()};
        context.executeKernel(ccmaDirectionsKernel, directionsArgs, ccmaAtoms.getSize());
        int i;
        void* forceArgs[] = {&ccmaAtoms.getDevicePointer(), &ccmaDistance.getDevicePointer(),
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta.getDevicePointer(),
                &ccmaReducedMass.getDevicePointer(), &ccmaDelta1.getDevicePointer(), &ccmaConverged.getDevicePointer(),
                &ccmaConvergedDeviceMemory, tolPointer, &i};
        void* multiplyArgs[] = {&ccmaDelta1.getDevicePointer(), &ccmaDelta2.getDevicePointer(),
                &ccmaConstraintMatrixColumn.getDevicePointer(), &ccmaConstraintMatrixValue.getDevicePointer(), &ccmaConverged.getDevicePointer(), &i};
        void* updateArgs[] = {&ccmaNumAtomConstraints.getDevicePointer(), &ccmaAtomConstraints.getDevicePointer(), &ccmaDistance.getDevicePointer(),
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta.getDevicePointer(),
                &context.getVelm().getDevicePointer(), &ccmaDelta1.getDevicePointer(), &ccmaDelta2.getDevicePointer(),
                &ccmaConverged.getDevicePointer(), &i};
        const int checkInterval = 4;
        ccmaConvergedMemory[0] = 0;
        for (i = 0; i < 150; i++) {
            context.executeKernel(ccmaForceKernel, forceArgs, ccmaAtoms.getSize());
            if ((i+1)%checkInterval == 0)
                CHECK_RESULT2(cuEventRecord(ccmaEvent, 0), "Error recording event for CCMA");
            context.executeKernel(ccmaMultiplyKernel, multiplyArgs, ccmaAtoms.getSize());
            context.executeKernel(ccmaUpdateKernel, updateArgs, context.getNumAtoms());
            if ((i+1)%checkInterval == 0) {
                CHECK_RESULT2(cuEventSynchronize(ccmaEvent), "Error synchronizing on event for CCMA");
                if (ccmaConvergedMemory[0])
                    break;
            }
        }
    }
}

void CudaIntegrationUtilities::computeVirtualSites() {
    if (numVsites > 0) {
        CUdeviceptr posCorrection = (context.getUseMixedPrecision() ? context.getPosqCorrection().getDevicePointer() : 0);
        void* args[] = {&context.getPosq().getDevicePointer(), &posCorrection, &vsite2AvgAtoms.getDevicePointer(), &vsite2AvgWeights.getDevicePointer(),
                &vsite3AvgAtoms.getDevicePointer(), &vsite3AvgWeights.getDevicePointer(),
                &vsiteOutOfPlaneAtoms.getDevicePointer(), &vsiteOutOfPlaneWeights.getDevicePointer(),
                &vsiteLocalCoordsIndex.getDevicePointer(), &vsiteLocalCoordsAtoms.getDevicePointer(),
                &vsiteLocalCoordsWeights.getDevicePointer(), &vsiteLocalCoordsPos.getDevicePointer(),
                &vsiteLocalCoordsStartIndex.getDevicePointer()};
        context.executeKernel(vsitePositionKernel, args, numVsites);
    }
}

void CudaIntegrationUtilities::distributeForcesFromVirtualSites() {
    if (numVsites > 0) {
        CUdeviceptr posCorrection = (context.getUseMixedPrecision() ? context.getPosqCorrection().getDevicePointer() : 0);
        void* args[] = {&context.getPosq().getDevicePointer(), &posCorrection, &context.getForce().getDevicePointer(),
                &vsite2AvgAtoms.getDevicePointer(), &vsite2AvgWeights.getDevicePointer(),
                &vsite3AvgAtoms.getDevicePointer(), &vsite3AvgWeights.getDevicePointer(),
                &vsiteOutOfPlaneAtoms.getDevicePointer(), &vsiteOutOfPlaneWeights.getDevicePointer(),
                &vsiteLocalCoordsIndex.getDevicePointer(), &vsiteLocalCoordsAtoms.getDevicePointer(),
                &vsiteLocalCoordsWeights.getDevicePointer(), &vsiteLocalCoordsPos.getDevicePointer(),
                &vsiteLocalCoordsStartIndex.getDevicePointer()};
        context.executeKernel(vsiteForceKernel, args, numVsites);
    }
}

void CudaIntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random.isInitialized()) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("CudaIntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random.initialize<float4>(context, 4*context.getPaddedNumAtoms(), "random");
    randomSeed.initialize<int4>(context, context.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
    randomPos = random.getSize();

    // Use a quick and dirty RNG to pick seeds for the real random number generator.

    vector<int4> seed(randomSeed.getSize());
    unsigned int r = randomNumberSeed;
    if (r == 0) r = (unsigned int) osrngseed();
    for (int i = 0; i < randomSeed.getSize(); i++) {
        seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
    }
    randomSeed.upload(seed);
}

int CudaIntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random.getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    if (numValues > random.getSize())
        random.resize(numValues);
    int size = random.getSize();
    void* args[] = {&size, &random.getDevicePointer(), &randomSeed.getDevicePointer()};
    context.executeKernel(randomKernel, args, random.getSize());
    randomPos = numValues;
    return 0;
}

void CudaIntegrationUtilities::createCheckpoint(ostream& stream) {
    if (!random.isInitialized()) 
        return;
    stream.write((char*) &randomPos, sizeof(int));
    vector<float4> randomVec;
    random.download(randomVec);
    stream.write((char*) &randomVec[0], sizeof(float4)*random.getSize());
    vector<int4> randomSeedVec;
    randomSeed.download(randomSeedVec);
    stream.write((char*) &randomSeedVec[0], sizeof(int4)*randomSeed.getSize());
}

void CudaIntegrationUtilities::loadCheckpoint(istream& stream) {
    if (!random.isInitialized()) 
        return;
    stream.read((char*) &randomPos, sizeof(int));
    vector<float4> randomVec(random.getSize());
    stream.read((char*) &randomVec[0], sizeof(float4)*random.getSize());
    random.upload(randomVec);
    vector<int4> randomSeedVec(randomSeed.getSize());
    stream.read((char*) &randomSeedVec[0], sizeof(int4)*randomSeed.getSize());
    randomSeed.upload(randomSeedVec);
}

double CudaIntegrationUtilities::computeKineticEnergy(double timeShift) {
    int numParticles = context.getNumAtoms();
    if (timeShift != 0) {
        float timeShiftFloat = (float) timeShift;
        void* timeShiftPtr = (context.getUseDoublePrecision() ? (void*) &timeShift : (void*) &timeShiftFloat);

        // Copy the velocities into the posDelta array while we temporarily modify them.

        context.getVelm().copyTo(posDelta);

        // Apply the time shift.

        void* args[] = {&context.getVelm().getDevicePointer(), &context.getForce().getDevicePointer(), timeShiftPtr};
        context.executeKernel(timeShiftKernel, args, numParticles);
        applyConstraints(true, 1e-4);
    }
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
        vector<double4> velm;
        context.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++) {
            double4 v = velm[i];
            if (v.w != 0)
                energy += (v.x*v.x+v.y*v.y+v.z*v.z)/v.w;
        }
    }
    else {
        vector<float4> velm;
        context.getVelm().download(velm);
        for (int i = 0; i < numParticles; i++) {
            float4 v = velm[i];
            if (v.w != 0)
                energy += (v.x*v.x+v.y*v.y+v.z*v.z)/v.w;
        }
    }
    
    // Restore the velocities.
    
    if (timeShift != 0)
        posDelta.copyTo(context.getVelm());
    return 0.5*energy;
}
