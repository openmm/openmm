/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2014 Stanford University and the Authors.      *
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
#include "openmm/HarmonicAngleForce.h"
#include "openmm/VirtualSite.h"
#include "quern.h"
#include "CudaExpressionUtilities.h"
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
        posDelta(NULL), settleAtoms(NULL), settleParams(NULL), shakeAtoms(NULL), shakeParams(NULL),
        random(NULL), randomSeed(NULL), randomPos(0), stepSize(NULL), ccmaAtoms(NULL), ccmaDistance(NULL),
        ccmaReducedMass(NULL), ccmaAtomConstraints(NULL), ccmaNumAtomConstraints(NULL), ccmaConstraintMatrixColumn(NULL),
        ccmaConstraintMatrixValue(NULL), ccmaDelta1(NULL), ccmaDelta2(NULL), ccmaConverged(NULL), ccmaConvergedMemory(NULL),
        vsite2AvgAtoms(NULL), vsite2AvgWeights(NULL), vsite3AvgAtoms(NULL), vsite3AvgWeights(NULL),
        vsiteOutOfPlaneAtoms(NULL), vsiteOutOfPlaneWeights(NULL), vsiteLocalCoordsAtoms(NULL), vsiteLocalCoordsParams(NULL) {
    // Create workspace arrays.

    if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
        posDelta = CudaArray::create<double4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<double4> deltas(posDelta->getSize(), make_double4(0.0, 0.0, 0.0, 0.0));
        posDelta->upload(deltas);
        stepSize = CudaArray::create<double2>(context, 1, "stepSize");
        vector<double2> step(1, make_double2(0.0, 0.0));
        stepSize->upload(step);
    }
    else {
        posDelta = CudaArray::create<float4>(context, context.getPaddedNumAtoms(), "posDelta");
        vector<float4> deltas(posDelta->getSize(), make_float4(0.0f, 0.0f, 0.0f, 0.0f));
        posDelta->upload(deltas);
        stepSize = CudaArray::create<float2>(context, 1, "stepSize");
        vector<float2> step(1, make_float2(0.0f, 0.0f));
        stepSize->upload(step);
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
                throw OpenMMException("Two of the three distances constrained with SETTLE must be the same.");
            isShakeAtom[atom1] = true;
            isShakeAtom[atom2] = true;
            isShakeAtom[atom3] = true;
        }
        settleAtoms = CudaArray::create<int4>(context, atoms.size(), "settleAtoms");
        settleParams = CudaArray::create<float2>(context, params.size(), "settleParams");
        settleAtoms->upload(atoms);
        settleParams->upload(params);
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
        shakeAtoms = CudaArray::create<int4>(context, atoms.size(), "shakeAtoms");
        shakeParams = CudaArray::create<float4>(context, params.size(), "shakeParams");
        shakeAtoms->upload(atoms);
        shakeParams->upload(params);
    }

    // Find connected constraints for CCMA.

    vector<int> ccmaConstraints;
    for (unsigned i = 0; i < atom1.size(); i++)
        if (!isShakeAtom[atom1[i]])
            ccmaConstraints.push_back(i);

    // Record the connections between constraints.

    int numCCMA = (int) ccmaConstraints.size();
    if (numCCMA > 0) {
        vector<vector<int> > atomConstraints(context.getNumAtoms());
        for (int i = 0; i < numCCMA; i++) {
            atomConstraints[atom1[ccmaConstraints[i]]].push_back(i);
            atomConstraints[atom2[ccmaConstraints[i]]].push_back(i);
        }
        vector<vector<int> > linkedConstraints(numCCMA);
        for (unsigned atom = 0; atom < atomConstraints.size(); atom++) {
            for (unsigned i = 0; i < atomConstraints[atom].size(); i++)
                for (unsigned j = 0; j < i; j++) {
                    int c1 = atomConstraints[atom][i];
                    int c2 = atomConstraints[atom][j];
                    linkedConstraints[c1].push_back(c2);
                    linkedConstraints[c2].push_back(c1);
                }
        }
        int maxLinks = 0;
        for (unsigned i = 0; i < linkedConstraints.size(); i++)
            maxLinks = max(maxLinks, (int) linkedConstraints[i].size());
        int maxAtomConstraints = 0;
        for (unsigned i = 0; i < atomConstraints.size(); i++)
            maxAtomConstraints = max(maxAtomConstraints, (int) atomConstraints[i].size());

        // Compute the constraint coupling matrix

        vector<vector<int> > atomAngles(numAtoms);
        HarmonicAngleForce const* angleForce = NULL;
        for (int i = 0; i < system.getNumForces() && angleForce == NULL; i++)
            angleForce = dynamic_cast<HarmonicAngleForce const*>(&system.getForce(i));
        if (angleForce != NULL)
            for (int i = 0; i < angleForce->getNumAngles(); i++) {
                int particle1, particle2, particle3;
                double angle, k;
                angleForce->getAngleParameters(i, particle1, particle2, particle3, angle, k);
                atomAngles[particle2].push_back(i);
            }
        vector<vector<pair<int, double> > > matrix(numCCMA);
        for (int j = 0; j < numCCMA; j++) {
            for (int k = 0; k < numCCMA; k++) {
                if (j == k) {
                    matrix[j].push_back(pair<int, double>(j, 1.0));
                    continue;
                }
                double scale;
                int cj = ccmaConstraints[j];
                int ck = ccmaConstraints[k];
                int atomj0 = atom1[cj];
                int atomj1 = atom2[cj];
                int atomk0 = atom1[ck];
                int atomk1 = atom2[ck];
                int atoma, atomb, atomc;
                double imj0 = 1.0/system.getParticleMass(atomj0);
                double imj1 = 1.0/system.getParticleMass(atomj1);
                if (atomj0 == atomk0) {
                    atoma = atomj1;
                    atomb = atomj0;
                    atomc = atomk1;
                    scale = imj0/(imj0+imj1);
                }
                else if (atomj1 == atomk1) {
                    atoma = atomj0;
                    atomb = atomj1;
                    atomc = atomk0;
                    scale = imj1/(imj0+imj1);
                }
                else if (atomj0 == atomk1) {
                    atoma = atomj1;
                    atomb = atomj0;
                    atomc = atomk0;
                    scale = imj0/(imj0+imj1);
                }
                else if (atomj1 == atomk0) {
                    atoma = atomj0;
                    atomb = atomj1;
                    atomc = atomk1;
                    scale = imj1/(imj0+imj1);
                }
                else
                    continue; // These constraints are not connected.

                // Look for a third constraint forming a triangle with these two.

                bool foundConstraint = false;
                for (int m = 0; m < numCCMA; m++) {
                    int other = ccmaConstraints[m];
                    if ((atom1[other] == atoma && atom2[other] == atomc) || (atom1[other] == atomc && atom2[other] == atoma)) {
                        double d1 = distance[cj];
                        double d2 = distance[ck];
                        double d3 = distance[other];
                        matrix[j].push_back(pair<int, double>(k, scale*(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2)));
                        foundConstraint = true;
                        break;
                    }
                }
                if (!foundConstraint && angleForce != NULL) {
                    // We didn't find one, so look for an angle force field term.

                    const vector<int>& angleCandidates = atomAngles[atomb];
                    for (vector<int>::const_iterator iter = angleCandidates.begin(); iter != angleCandidates.end(); iter++) {
                        int particle1, particle2, particle3;
                        double angle, ka;
                        angleForce->getAngleParameters(*iter, particle1, particle2, particle3, angle, ka);
                        if ((particle1 == atoma && particle3 == atomc) || (particle3 == atoma && particle1 == atomc)) {
                            matrix[j].push_back(pair<int, double>(k, scale*cos(angle)));
                            break;
                        }
                    }
                }
            }
        }

        // Invert it using QR.

        vector<int> matrixRowStart;
        vector<int> matrixColIndex;
        vector<double> matrixValue;
        for (int i = 0; i < numCCMA; i++) {
            matrixRowStart.push_back(matrixValue.size());
            for (int j = 0; j < (int) matrix[i].size(); j++) {
                pair<int, double> element = matrix[i][j];
                matrixColIndex.push_back(element.first);
                matrixValue.push_back(element.second);
            }
        }
        matrixRowStart.push_back(matrixValue.size());
        int *qRowStart, *qColIndex, *rRowStart, *rColIndex;
        double *qValue, *rValue;
        int result = QUERN_compute_qr(numCCMA, numCCMA, &matrixRowStart[0], &matrixColIndex[0], &matrixValue[0], NULL,
                &qRowStart, &qColIndex, &qValue, &rRowStart, &rColIndex, &rValue);
        vector<double> rhs(numCCMA);
        matrix.clear();
        matrix.resize(numCCMA);
        for (int i = 0; i < numCCMA; i++) {
            // Extract column i of the inverse matrix.

            for (int j = 0; j < numCCMA; j++)
                rhs[j] = (i == j ? 1.0 : 0.0);
            result = QUERN_multiply_with_q_transpose(numCCMA, qRowStart, qColIndex, qValue, &rhs[0]);
            result = QUERN_solve_with_r(numCCMA, rRowStart, rColIndex, rValue, &rhs[0], &rhs[0]);
            for (int j = 0; j < numCCMA; j++) {
                double value = rhs[j]*distance[ccmaConstraints[i]]/distance[ccmaConstraints[j]];
                if (abs(value) > 0.1)
                    matrix[j].push_back(pair<int, double>(i, value));
            }
        }
        QUERN_free_result(qRowStart, qColIndex, qValue);
        QUERN_free_result(rRowStart, rColIndex, rValue);
        int maxRowElements = 0;
        for (unsigned i = 0; i < matrix.size(); i++)
            maxRowElements = max(maxRowElements, (int) matrix[i].size());
        maxRowElements++;

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

        ccmaAtoms = CudaArray::create<int2>(context, numCCMA, "CcmaAtoms");
        ccmaAtomConstraints = CudaArray::create<int>(context, numAtoms*maxAtomConstraints, "CcmaAtomConstraints");
        ccmaNumAtomConstraints = CudaArray::create<int>(context, numAtoms, "CcmaAtomConstraintsIndex");
        ccmaConstraintMatrixColumn = CudaArray::create<int>(context, numCCMA*maxRowElements, "ConstraintMatrixColumn");
        ccmaConverged = CudaArray::create<int>(context, 2, "ccmaConverged");
        CHECK_RESULT2(cuMemHostAlloc((void**) &ccmaConvergedMemory, sizeof(int), CU_MEMHOSTALLOC_DEVICEMAP), "Error allocating pinned memory");
        CHECK_RESULT2(cuMemHostGetDevicePointer(&ccmaConvergedDeviceMemory, ccmaConvergedMemory, 0), "Error getting device address for pinned memory");
        vector<int2> atomsVec(ccmaAtoms->getSize());
        vector<int> atomConstraintsVec(ccmaAtomConstraints->getSize());
        vector<int> numAtomConstraintsVec(ccmaNumAtomConstraints->getSize());
        vector<int> constraintMatrixColumnVec(ccmaConstraintMatrixColumn->getSize());
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision()) {
            ccmaDistance = CudaArray::create<double4>(context, numCCMA, "CcmaDistance");
            ccmaDelta1 = CudaArray::create<double>(context, numCCMA, "CcmaDelta1");
            ccmaDelta2 = CudaArray::create<double>(context, numCCMA, "CcmaDelta2");
            ccmaReducedMass = CudaArray::create<double>(context, numCCMA, "CcmaReducedMass");
            ccmaConstraintMatrixValue = CudaArray::create<double>(context, numCCMA*maxRowElements, "ConstraintMatrixValue");
            vector<double4> distanceVec(ccmaDistance->getSize());
            vector<double> reducedMassVec(ccmaReducedMass->getSize());
            vector<double> constraintMatrixValueVec(ccmaConstraintMatrixValue->getSize());
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
            ccmaDistance->upload(distanceVec);
            ccmaReducedMass->upload(reducedMassVec);
            ccmaConstraintMatrixValue->upload(constraintMatrixValueVec);
        }
        else {
            ccmaDistance = CudaArray::create<float4>(context, numCCMA, "CcmaDistance");
            ccmaDelta1 = CudaArray::create<float>(context, numCCMA, "CcmaDelta1");
            ccmaDelta2 = CudaArray::create<float>(context, numCCMA, "CcmaDelta2");
            ccmaReducedMass = CudaArray::create<float>(context, numCCMA, "CcmaReducedMass");
            ccmaConstraintMatrixValue = CudaArray::create<float>(context, numCCMA*maxRowElements, "ConstraintMatrixValue");
            vector<float4> distanceVec(ccmaDistance->getSize());
            vector<float> reducedMassVec(ccmaReducedMass->getSize());
            vector<float> constraintMatrixValueVec(ccmaConstraintMatrixValue->getSize());
            for (int i = 0; i < numCCMA; i++) {
                int index = constraintOrder[i];
                int c = ccmaConstraints[index];
                atomsVec[i].x = atom1[c];
                atomsVec[i].y = atom2[c];
                distanceVec[i].w = (float) distance[c];
                reducedMassVec[i] = (float) (0.5/(1.0/system.getParticleMass(atom1[c])+1.0/system.getParticleMass(atom2[c])));
                for (unsigned int j = 0; j < matrix[index].size(); j++) {
                    constraintMatrixColumnVec[i+j*numCCMA] = matrix[index][j].first;
                    constraintMatrixValueVec[i+j*numCCMA] = (float) matrix[index][j].second;
                }
                constraintMatrixColumnVec[i+matrix[index].size()*numCCMA] = numCCMA;
            }
            ccmaDistance->upload(distanceVec);
            ccmaReducedMass->upload(reducedMassVec);
            ccmaConstraintMatrixValue->upload(constraintMatrixValueVec);
        }
        for (unsigned int i = 0; i < atomConstraints.size(); i++) {
            numAtomConstraintsVec[i] = atomConstraints[i].size();
            for (unsigned int j = 0; j < atomConstraints[i].size(); j++) {
                bool forward = (atom1[ccmaConstraints[atomConstraints[i][j]]] == i);
                atomConstraintsVec[i+j*numAtoms] = (forward ? inverseOrder[atomConstraints[i][j]]+1 : -inverseOrder[atomConstraints[i][j]]-1);
            }
        }
        ccmaAtoms->upload(atomsVec);
        ccmaAtomConstraints->upload(atomConstraintsVec);
        ccmaNumAtomConstraints->upload(numAtomConstraintsVec);
        ccmaConstraintMatrixColumn->upload(constraintMatrixColumnVec);
    }
    
    // Build the list of virtual sites.
    
    vector<int4> vsite2AvgAtomVec;
    vector<double2> vsite2AvgWeightVec;
    vector<int4> vsite3AvgAtomVec;
    vector<double4> vsite3AvgWeightVec;
    vector<int4> vsiteOutOfPlaneAtomVec;
    vector<double4> vsiteOutOfPlaneWeightVec;
    vector<int4> vsiteLocalCoordsAtomVec;
    vector<double> vsiteLocalCoordsParamVec;
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
                // An out of plane site.
                
                const LocalCoordinatesSite& site = dynamic_cast<const LocalCoordinatesSite&>(system.getVirtualSite(i));
                vsiteLocalCoordsAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                Vec3 origin = site.getOriginWeights();
                Vec3 x = site.getXWeights();
                Vec3 y = site.getYWeights();
                Vec3 pos = site.getLocalPosition();
                vsiteLocalCoordsParamVec.push_back(origin[0]);
                vsiteLocalCoordsParamVec.push_back(origin[1]);
                vsiteLocalCoordsParamVec.push_back(origin[2]);
                vsiteLocalCoordsParamVec.push_back(x[0]);
                vsiteLocalCoordsParamVec.push_back(x[1]);
                vsiteLocalCoordsParamVec.push_back(x[2]);
                vsiteLocalCoordsParamVec.push_back(y[0]);
                vsiteLocalCoordsParamVec.push_back(y[1]);
                vsiteLocalCoordsParamVec.push_back(y[2]);
                vsiteLocalCoordsParamVec.push_back(pos[0]);
                vsiteLocalCoordsParamVec.push_back(pos[1]);
                vsiteLocalCoordsParamVec.push_back(pos[2]);
            }
        }
    }
    int num2Avg = vsite2AvgAtomVec.size();
    int num3Avg = vsite3AvgAtomVec.size();
    int numOutOfPlane = vsiteOutOfPlaneAtomVec.size();
    int numLocalCoords = vsiteLocalCoordsAtomVec.size();
    vsite2AvgAtoms = CudaArray::create<int4>(context, max(1, num2Avg), "vsite2AvgAtoms");
    vsite3AvgAtoms = CudaArray::create<int4>(context, max(1, num3Avg), "vsite3AvgAtoms");
    vsiteOutOfPlaneAtoms = CudaArray::create<int4>(context, max(1, numOutOfPlane), "vsiteOutOfPlaneAtoms");
    vsiteLocalCoordsAtoms = CudaArray::create<int4>(context, max(1, numLocalCoords), "vsiteLocalCoordinatesAtoms");
    if (num2Avg > 0)
        vsite2AvgAtoms->upload(vsite2AvgAtomVec);
    if (num3Avg > 0)
        vsite3AvgAtoms->upload(vsite3AvgAtomVec);
    if (numOutOfPlane > 0)
        vsiteOutOfPlaneAtoms->upload(vsiteOutOfPlaneAtomVec);
    if (numLocalCoords > 0)
        vsiteLocalCoordsAtoms->upload(vsiteLocalCoordsAtomVec);
    if (context.getUseDoublePrecision()) {
        vsite2AvgWeights = CudaArray::create<double2>(context, max(1, num2Avg), "vsite2AvgWeights");
        vsite3AvgWeights = CudaArray::create<double4>(context, max(1, num3Avg), "vsite3AvgWeights");
        vsiteOutOfPlaneWeights = CudaArray::create<double4>(context, max(1, numOutOfPlane), "vsiteOutOfPlaneWeights");
        vsiteLocalCoordsParams = CudaArray::create<double>(context, max(1, 12*numLocalCoords), "vsiteLocalCoordinatesParams");
        if (num2Avg > 0)
            vsite2AvgWeights->upload(vsite2AvgWeightVec);
        if (num3Avg > 0)
            vsite3AvgWeights->upload(vsite3AvgWeightVec);
        if (numOutOfPlane > 0)
            vsiteOutOfPlaneWeights->upload(vsiteOutOfPlaneWeightVec);
        if (numLocalCoords > 0)
            vsiteLocalCoordsParams->upload(vsiteLocalCoordsParamVec);
    }
    else {
        vsite2AvgWeights = CudaArray::create<float2>(context, max(1, num2Avg), "vsite2AvgWeights");
        vsite3AvgWeights = CudaArray::create<float4>(context, max(1, num3Avg), "vsite3AvgWeights");
        vsiteOutOfPlaneWeights = CudaArray::create<float4>(context, max(1, numOutOfPlane), "vsiteOutOfPlaneWeights");
        vsiteLocalCoordsParams = CudaArray::create<float>(context, max(1, 12*numLocalCoords), "vsiteLocalCoordinatesParams");
        if (num2Avg > 0) {
            vector<float2> floatWeights(num2Avg);
            for (int i = 0; i < num2Avg; i++)
                floatWeights[i] = make_float2((float) vsite2AvgWeightVec[i].x, (float) vsite2AvgWeightVec[i].y);
            vsite2AvgWeights->upload(floatWeights);
        }
        if (num3Avg > 0) {
            vector<float4> floatWeights(num3Avg);
            for (int i = 0; i < num3Avg; i++)
                floatWeights[i] = make_float4((float) vsite3AvgWeightVec[i].x, (float) vsite3AvgWeightVec[i].y, (float) vsite3AvgWeightVec[i].z, 0.0f);
            vsite3AvgWeights->upload(floatWeights);
        }
        if (numOutOfPlane > 0) {
            vector<float4> floatWeights(numOutOfPlane);
            for (int i = 0; i < numOutOfPlane; i++)
                floatWeights[i] = make_float4((float) vsiteOutOfPlaneWeightVec[i].x, (float) vsiteOutOfPlaneWeightVec[i].y, (float) vsiteOutOfPlaneWeightVec[i].z, 0.0f);
            vsiteOutOfPlaneWeights->upload(floatWeights);
        }
        if (numLocalCoords > 0) {
            vector<float> floatParams(vsiteLocalCoordsParamVec.size());
            for (int i = 0; i < (int) vsiteLocalCoordsParamVec.size(); i++)
                floatParams[i] = (float) vsiteLocalCoordsParamVec[i];
            vsiteLocalCoordsParams->upload(floatParams);
        }
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
    if (posDelta != NULL)
        delete posDelta;
    if (settleAtoms != NULL)
        delete settleAtoms;
    if (settleParams != NULL)
        delete settleParams;
    if (shakeAtoms != NULL)
        delete shakeAtoms;
    if (shakeParams != NULL)
        delete shakeParams;
    if (random != NULL)
        delete random;
    if (randomSeed != NULL)
        delete randomSeed;
    if (stepSize != NULL)
        delete stepSize;
    if (ccmaAtoms != NULL)
        delete ccmaAtoms;
    if (ccmaDistance != NULL)
        delete ccmaDistance;
    if (ccmaReducedMass != NULL)
        delete ccmaReducedMass;
    if (ccmaAtomConstraints != NULL)
        delete ccmaAtomConstraints;
    if (ccmaNumAtomConstraints != NULL)
        delete ccmaNumAtomConstraints;
    if (ccmaConstraintMatrixColumn != NULL)
        delete ccmaConstraintMatrixColumn;
    if (ccmaConstraintMatrixValue != NULL)
        delete ccmaConstraintMatrixValue;
    if (ccmaDelta1 != NULL)
        delete ccmaDelta1;
    if (ccmaDelta2 != NULL)
        delete ccmaDelta2;
    if (ccmaConverged != NULL)
        delete ccmaConverged;
    if (ccmaConvergedMemory != NULL)
        cuMemFreeHost(ccmaConvergedMemory);
    if (vsite2AvgAtoms != NULL)
        delete vsite2AvgAtoms;
    if (vsite2AvgWeights != NULL)
        delete vsite2AvgWeights;
    if (vsite3AvgAtoms != NULL)
        delete vsite3AvgAtoms;
    if (vsite3AvgWeights != NULL)
        delete vsite3AvgWeights;
    if (vsiteOutOfPlaneAtoms != NULL)
        delete vsiteOutOfPlaneAtoms;
    if (vsiteOutOfPlaneWeights != NULL)
        delete vsiteOutOfPlaneWeights;
    if (vsiteLocalCoordsAtoms != NULL)
        delete vsiteLocalCoordsAtoms;
    if (vsiteLocalCoordsParams != NULL)
        delete vsiteLocalCoordsParams;
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
    if (settleAtoms != NULL) {
        int numClusters = settleAtoms->getSize();
        void* args[] = {&numClusters, tolPointer, &context.getPosq().getDevicePointer(), &posCorrection,
                &posDelta->getDevicePointer(), &context.getVelm().getDevicePointer(),
                &settleAtoms->getDevicePointer(), &settleParams->getDevicePointer()};
        context.executeKernel(settleKernel, args, settleAtoms->getSize());
    }
    if (shakeAtoms != NULL) {
        int numClusters = shakeAtoms->getSize();
        void* args[] = {&numClusters, tolPointer, &context.getPosq().getDevicePointer(), &posCorrection,
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta->getDevicePointer(),
                &shakeAtoms->getDevicePointer(), &shakeParams->getDevicePointer()};
        context.executeKernel(shakeKernel, args, shakeAtoms->getSize());
    }
    if (ccmaAtoms != NULL) {
        void* directionsArgs[] = {&ccmaAtoms->getDevicePointer(), &ccmaDistance->getDevicePointer(), &context.getPosq().getDevicePointer(), &posCorrection, &ccmaConverged->getDevicePointer()};
        context.executeKernel(ccmaDirectionsKernel, directionsArgs, ccmaAtoms->getSize());
        int i;
        void* forceArgs[] = {&ccmaAtoms->getDevicePointer(), &ccmaDistance->getDevicePointer(),
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta->getDevicePointer(),
                &ccmaReducedMass->getDevicePointer(), &ccmaDelta1->getDevicePointer(), &ccmaConverged->getDevicePointer(),
                &ccmaConvergedDeviceMemory, tolPointer, &i};
        void* multiplyArgs[] = {&ccmaDelta1->getDevicePointer(), &ccmaDelta2->getDevicePointer(),
                &ccmaConstraintMatrixColumn->getDevicePointer(), &ccmaConstraintMatrixValue->getDevicePointer(), &ccmaConverged->getDevicePointer(), &i};
        void* updateArgs[] = {&ccmaNumAtomConstraints->getDevicePointer(), &ccmaAtomConstraints->getDevicePointer(), &ccmaDistance->getDevicePointer(),
                constrainVelocities ? &context.getVelm().getDevicePointer() : &posDelta->getDevicePointer(),
                &context.getVelm().getDevicePointer(), &ccmaDelta1->getDevicePointer(), &ccmaDelta2->getDevicePointer(),
                &ccmaConverged->getDevicePointer(), &i};
        const int checkInterval = 4;
        ccmaConvergedMemory[0] = 0;
        for (i = 0; i < 150; i++) {
            context.executeKernel(ccmaForceKernel, forceArgs, ccmaAtoms->getSize());
            if ((i+1)%checkInterval == 0)
                CHECK_RESULT2(cuEventRecord(ccmaEvent, 0), "Error recording event for CCMA");
            context.executeKernel(ccmaMultiplyKernel, multiplyArgs, ccmaAtoms->getSize());
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
        void* args[] = {&context.getPosq().getDevicePointer(), &posCorrection, &vsite2AvgAtoms->getDevicePointer(), &vsite2AvgWeights->getDevicePointer(),
                &vsite3AvgAtoms->getDevicePointer(), &vsite3AvgWeights->getDevicePointer(),
                &vsiteOutOfPlaneAtoms->getDevicePointer(), &vsiteOutOfPlaneWeights->getDevicePointer(),
                &vsiteLocalCoordsAtoms->getDevicePointer(), &vsiteLocalCoordsParams->getDevicePointer()};
        context.executeKernel(vsitePositionKernel, args, numVsites);
    }
}

void CudaIntegrationUtilities::distributeForcesFromVirtualSites() {
    if (numVsites > 0) {
        CUdeviceptr posCorrection = (context.getUseMixedPrecision() ? context.getPosqCorrection().getDevicePointer() : 0);
        void* args[] = {&context.getPosq().getDevicePointer(), &posCorrection, &context.getForce().getDevicePointer(),
                &vsite2AvgAtoms->getDevicePointer(), &vsite2AvgWeights->getDevicePointer(),
                &vsite3AvgAtoms->getDevicePointer(), &vsite3AvgWeights->getDevicePointer(),
                &vsiteOutOfPlaneAtoms->getDevicePointer(), &vsiteOutOfPlaneWeights->getDevicePointer(),
                &vsiteLocalCoordsAtoms->getDevicePointer(), &vsiteLocalCoordsParams->getDevicePointer()};
        context.executeKernel(vsiteForceKernel, args, numVsites);
    }
}

void CudaIntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random != NULL) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("CudaIntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random = CudaArray::create<float4>(context, 4*context.getPaddedNumAtoms(), "random");
    randomSeed = CudaArray::create<int4>(context, context.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
    randomPos = random->getSize();

    // Use a quick and dirty RNG to pick seeds for the real random number generator.

    vector<int4> seed(randomSeed->getSize());
    unsigned int r = randomNumberSeed;
    for (int i = 0; i < randomSeed->getSize(); i++) {
        seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
    }
    randomSeed->upload(seed);
}

int CudaIntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random->getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    if (numValues > random->getSize()) {
        delete random;
        random = CudaArray::create<float4>(context, numValues, "random");
    }
    int size = random->getSize();
    void* args[] = {&size, &random->getDevicePointer(), &randomSeed->getDevicePointer()};
    context.executeKernel(randomKernel, args, random->getSize());
    randomPos = numValues;
    return 0;
}

void CudaIntegrationUtilities::createCheckpoint(ostream& stream) {
    if(random == NULL) 
        return;
    stream.write((char*) &randomPos, sizeof(int));
    vector<float4> randomVec;
    random->download(randomVec);
    stream.write((char*) &randomVec[0], sizeof(float4)*random->getSize());
    vector<int4> randomSeedVec;
    randomSeed->download(randomSeedVec);
    stream.write((char*) &randomSeedVec[0], sizeof(int4)*randomSeed->getSize());
}

void CudaIntegrationUtilities::loadCheckpoint(istream& stream) {
    if(random == NULL) 
        return;
    stream.read((char*) &randomPos, sizeof(int));
    vector<float4> randomVec(random->getSize());
    stream.read((char*) &randomVec[0], sizeof(float4)*random->getSize());
    random->upload(randomVec);
    vector<int4> randomSeedVec(randomSeed->getSize());
    stream.read((char*) &randomSeedVec[0], sizeof(int4)*randomSeed->getSize());
    randomSeed->upload(randomSeedVec);
}

double CudaIntegrationUtilities::computeKineticEnergy(double timeShift) {
    int numParticles = context.getNumAtoms();
    if (timeShift != 0) {
        float timeShiftFloat = (float) timeShift;
        void* timeShiftPtr = (context.getUseDoublePrecision() ? (void*) &timeShift : (void*) &timeShiftFloat);

        // Copy the velocities into the posDelta array while we temporarily modify them.

        context.getVelm().copyTo(*posDelta);

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
        posDelta->copyTo(context.getVelm());
    return 0.5*energy;
}
