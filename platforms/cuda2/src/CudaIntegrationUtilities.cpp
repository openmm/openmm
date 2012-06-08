/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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
        ccmaConstraintMatrixValue(NULL), ccmaDelta1(NULL), ccmaDelta2(NULL), ccmaConverged(NULL),
        ccmaConvergedMemory(NULL), vsite2AvgAtoms(NULL), vsite2AvgWeights(NULL), vsite3AvgAtoms(NULL), vsite3AvgWeights(NULL),
        vsiteOutOfPlaneAtoms(NULL), vsiteOutOfPlaneWeights(NULL), hasInitializedPosConstraintKernels(false), hasInitializedVelConstraintKernels(false) {
    // Create workspace arrays.

    if (context.getUseDoublePrecision()) {
        posDelta = CudaArray::create<double4>(context.getPaddedNumAtoms(), "posDelta");
        vector<double4> deltas(posDelta->getSize(), make_double4(0.0, 0.0, 0.0, 0.0));
        posDelta->upload(deltas);
        stepSize = CudaArray::create<double2>(1, "stepSize");
        vector<double2> step(1, make_double2(0.0f, 0.0f));
        stepSize->upload(step);
    }
    else {
        posDelta = CudaArray::create<float4>(context.getPaddedNumAtoms(), "posDelta");
        vector<float4> deltas(posDelta->getSize(), make_float4(0.0, 0.0, 0.0, 0.0));
        posDelta->upload(deltas);
        stepSize = CudaArray::create<float2>(1, "stepSize");
        vector<float2> step(1, make_float2(0.0f, 0.0f));
        stepSize->upload(step);
    }

    // Create kernels for enforcing constraints.

    map<string, string> velocityDefines;
    velocityDefines["CONSTRAIN_VELOCITIES"] = "1";
//    CUmodule settleModule = context.createModule(CudaKernelSources::settle);
//    settlePosKernel = context.getKernel(settleModule, "applySettle");
//    settleVelKernel = context.getKernel(settleModule, "constrainVelocities");
//    CUmodule shakeModule = context.createModule(CudaKernelSources::shakeHydrogens);
//    shakePosKernel = context.getKernel(shakeModule, "applyShakeToHydrogens");
//    shakeModule = context.createModule(CudaKernelSources::shakeHydrogens, velocityDefines);
//    shakeVelKernel = context.getKernel(shakeModule, "applyShakeToHydrogens");

    // Record the set of constraints and how many constraints each atom is involved in.

    int numConstraints = system.getNumConstraints();
    vector<int> atom1(numConstraints);
    vector<int> atom2(numConstraints);
    vector<double> distance(numConstraints);
    vector<int> constraintCount(context.getNumAtoms(), 0);
    for (int i = 0; i < numConstraints; i++) {
        system.getConstraintParameters(i, atom1[i], atom2[i], distance[i]);
        constraintCount[atom1[i]]++;
        constraintCount[atom2[i]]++;
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
        settleAtoms = CudaArray::create<int4>(atoms.size(), "settleAtoms");
        settleParams = CudaArray::create<float2>(params.size(), "settleParams");
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
        shakeAtoms = CudaArray::create<int4>(atoms.size(), "shakeAtoms");
        shakeParams = CudaArray::create<float4>(params.size(), "shakeParams");
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

        ccmaAtoms = CudaArray::create<int2>(numCCMA, "CcmaAtoms");
        ccmaDistance = CudaArray::create<float4>(numCCMA, "CcmaDistance");
        ccmaAtomConstraints = CudaArray::create<int>(numAtoms*maxAtomConstraints, "CcmaAtomConstraints");
        ccmaNumAtomConstraints = CudaArray::create<int>(numAtoms, "CcmaAtomConstraintsIndex");
        ccmaDelta1 = CudaArray::create<float>(numCCMA, "CcmaDelta1");
        ccmaDelta2 = CudaArray::create<float>(numCCMA, "CcmaDelta2");
        ccmaConverged = CudaArray::create<int>(2, "CcmaConverged");
        CHECK_RESULT2(cuMemHostAlloc((void**) &ccmaConvergedMemory, 2*sizeof(int), 0), "Error allocating pinned memory");
        ccmaReducedMass = CudaArray::create<float>(numCCMA, "CcmaReducedMass");
        ccmaConstraintMatrixColumn = CudaArray::create<int>(numCCMA*maxRowElements, "ConstraintMatrixColumn");
        ccmaConstraintMatrixValue = CudaArray::create<float>(numCCMA*maxRowElements, "ConstraintMatrixValue");
        vector<int2> atomsVec(ccmaAtoms->getSize());
        vector<float4> distanceVec(ccmaDistance->getSize());
        vector<int> atomConstraintsVec(ccmaAtomConstraints->getSize());
        vector<int> numAtomConstraintsVec(ccmaNumAtomConstraints->getSize());
        vector<float> reducedMassVec(ccmaReducedMass->getSize());
        vector<int> constraintMatrixColumnVec(ccmaConstraintMatrixColumn->getSize());
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
        for (unsigned int i = 0; i < atomConstraints.size(); i++) {
            numAtomConstraintsVec[i] = atomConstraints[i].size();
            for (unsigned int j = 0; j < atomConstraints[i].size(); j++) {
                bool forward = (atom1[ccmaConstraints[atomConstraints[i][j]]] == i);
                atomConstraintsVec[i+j*numAtoms] = (forward ? inverseOrder[atomConstraints[i][j]]+1 : -inverseOrder[atomConstraints[i][j]]-1);
            }
        }
        ccmaAtoms->upload(atomsVec);
        ccmaDistance->upload(distanceVec);
        ccmaAtomConstraints->upload(atomConstraintsVec);
        ccmaNumAtomConstraints->upload(numAtomConstraintsVec);
        ccmaReducedMass->upload(reducedMassVec);
        ccmaConstraintMatrixColumn->upload(constraintMatrixColumnVec);
        ccmaConstraintMatrixValue->upload(constraintMatrixValueVec);

        // Create the CCMA kernels.

        map<string, string> defines;
        defines["NUM_CONSTRAINTS"] = context.intToString(numCCMA);
        defines["NUM_ATOMS"] = context.intToString(numAtoms);
//        CUmodule ccmaModule = context.createModule(CudaKernelSources::ccma, defines);
//        ccmaDirectionsKernel = context.getKernel(ccmaModule, "computeConstraintDirections");
//        ccmaPosForceKernel = context.getKernel(ccmaModule, "computeConstraintForce");
//        ccmaMultiplyKernel = context.getKernel(ccmaModule, "multiplyByConstraintMatrix");
//        ccmaPosUpdateKernel = context.getKernel(ccmaModule, "updateAtomPositions");
//        defines["CONSTRAIN_VELOCITIES"] = "1";
//        ccmaModule = context.createModule(CudaKernelSources::ccma, defines);
//        ccmaVelForceKernel = context.getKernel(ccmaModule, "computeConstraintForce");
//        ccmaVelUpdateKernel = context.getKernel(ccmaModule, "updateAtomPositions");
    }
    
    // Build the list of virtual sites.
    
    vector<int4> vsite2AvgAtomVec;
    vector<float2> vsite2AvgWeightVec;
    vector<int4> vsite3AvgAtomVec;
    vector<float4> vsite3AvgWeightVec;
    vector<int4> vsiteOutOfPlaneAtomVec;
    vector<float4> vsiteOutOfPlaneWeightVec;
    for (int i = 0; i < numAtoms; i++) {
        if (system.isVirtualSite(i)) {
            if (dynamic_cast<const TwoParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A two particle average.
                
                const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(system.getVirtualSite(i));
                vsite2AvgAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), 0));
                vsite2AvgWeightVec.push_back(make_float2((float) site.getWeight(0), (float) site.getWeight(1)));
            }
            else if (dynamic_cast<const ThreeParticleAverageSite*>(&system.getVirtualSite(i)) != NULL) {
                // A three particle average.
                
                const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(system.getVirtualSite(i));
                vsite3AvgAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsite3AvgWeightVec.push_back(make_float4((float) site.getWeight(0), (float) site.getWeight(1), (float) site.getWeight(2), 0.0f));
            }
            else if (dynamic_cast<const OutOfPlaneSite*>(&system.getVirtualSite(i)) != NULL) {
                // An out of plane site.
                
                const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(system.getVirtualSite(i));
                vsiteOutOfPlaneAtomVec.push_back(make_int4(i, site.getParticle(0), site.getParticle(1), site.getParticle(2)));
                vsiteOutOfPlaneWeightVec.push_back(make_float4((float) site.getWeight12(), (float) site.getWeight13(), (float) site.getWeightCross(), 0.0f));
            }
        }
    }
    int num2Avg = vsite2AvgAtomVec.size();
    int num3Avg = vsite3AvgAtomVec.size();
    int numOutOfPlane = vsiteOutOfPlaneAtomVec.size();
    vsite2AvgAtoms = CudaArray::create<int4>(max(1, num2Avg), "vsite2AvgAtoms");
    vsite2AvgWeights = CudaArray::create<float2>(max(1, num2Avg), "vsite2AvgWeights");
    vsite3AvgAtoms = CudaArray::create<int4>(max(1, num3Avg), "vsite3AvgAtoms");
    vsite3AvgWeights = CudaArray::create<float4>(max(1, num3Avg), "vsite3AvgWeights");
    vsiteOutOfPlaneAtoms = CudaArray::create<int4>(max(1, numOutOfPlane), "vsiteOutOfPlaneAtoms");
    vsiteOutOfPlaneWeights = CudaArray::create<float4>(max(1, numOutOfPlane), "vsiteOutOfPlaneWeights");
    if (num2Avg > 0) {
        vsite2AvgAtoms->upload(vsite2AvgAtomVec);
        vsite2AvgWeights->upload(vsite2AvgWeightVec);
    }
    if (num3Avg > 0) {
        vsite3AvgAtoms->upload(vsite3AvgAtomVec);
        vsite3AvgWeights->upload(vsite3AvgWeightVec);
    }
    if (numOutOfPlane > 0) {
        vsiteOutOfPlaneAtoms->upload(vsiteOutOfPlaneAtomVec);
        vsiteOutOfPlaneWeights->upload(vsiteOutOfPlaneWeightVec);
    }
    
    // Create the kernels for virtual sites.

    map<string, string> defines;
    defines["NUM_2_AVERAGE"] = context.intToString(num2Avg);
    defines["NUM_3_AVERAGE"] = context.intToString(num3Avg);
    defines["NUM_OUT_OF_PLANE"] = context.intToString(numOutOfPlane);
//    CUmodule vsiteModule = context.createModule(CudaKernelSources::virtualSites, defines);
//    vsitePositionKernel = context.getKernel(vsiteModule, "computeVirtualSites");
//    vsitePositionKernel.setArg<cl::Buffer>(0, context.getPosq().getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(1, vsite2AvgAtoms->getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(2, vsite2AvgWeights->getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(3, vsite3AvgAtoms->getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(4, vsite3AvgWeights->getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(5, vsiteOutOfPlaneAtoms->getDeviceBuffer());
//    vsitePositionKernel.setArg<cl::Buffer>(6, vsiteOutOfPlaneWeights->getDeviceBuffer());
//    vsiteForceKernel = context.getKernel(vsiteModule, "distributeForces");
//    vsiteForceKernel.setArg<cl::Buffer>(0, context.getPosq().getDeviceBuffer());
//    // Skip argument 1: the force array hasn't been created yet.
//    vsiteForceKernel.setArg<cl::Buffer>(2, vsite2AvgAtoms->getDeviceBuffer());
//    vsiteForceKernel.setArg<cl::Buffer>(3, vsite2AvgWeights->getDeviceBuffer());
//    vsiteForceKernel.setArg<cl::Buffer>(4, vsite3AvgAtoms->getDeviceBuffer());
//    vsiteForceKernel.setArg<cl::Buffer>(5, vsite3AvgWeights->getDeviceBuffer());
//    vsiteForceKernel.setArg<cl::Buffer>(6, vsiteOutOfPlaneAtoms->getDeviceBuffer());
//    vsiteForceKernel.setArg<cl::Buffer>(7, vsiteOutOfPlaneWeights->getDeviceBuffer());
    numVsites = num2Avg+num3Avg+numOutOfPlane;
}

CudaIntegrationUtilities::~CudaIntegrationUtilities() {
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
}

void CudaIntegrationUtilities::applyConstraints(double tol) {
    applyConstraints(false, tol);
}

void CudaIntegrationUtilities::applyVelocityConstraints(double tol) {
    applyConstraints(true, tol);
}

void CudaIntegrationUtilities::applyConstraints(bool constrainVelocities, double tol) {
//    bool hasInitialized;
//    CUfunction settleKernel, shakeKernel, ccmaForceKernel, ccmaUpdateKernel;
//    if (constrainVelocities) {
//        hasInitialized = hasInitializedVelConstraintKernels;
//        settleKernel = settleVelKernel;
//        shakeKernel = shakeVelKernel;
//        ccmaForceKernel = ccmaVelForceKernel;
//        ccmaUpdateKernel = ccmaVelUpdateKernel;
//        hasInitializedVelConstraintKernels = true;
//    }
//    else {
//        hasInitialized = hasInitializedPosConstraintKernels;
//        settleKernel = settlePosKernel;
//        shakeKernel = shakePosKernel;
//        ccmaForceKernel = ccmaPosForceKernel;
//        ccmaUpdateKernel = ccmaPosUpdateKernel;
//        hasInitializedPosConstraintKernels = true;
//    }
//    if (settleAtoms != NULL) {
//        if (!hasInitialized) {
//            settleKernel.setArg<int>(0, settleAtoms->getSize());
//            settleKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
//            settleKernel.setArg<cl::Buffer>(3, posDelta->getDeviceBuffer());
//            settleKernel.setArg<cl::Buffer>(4, context.getVelm().getDeviceBuffer());
//            settleKernel.setArg<cl::Buffer>(5, settleAtoms->getDeviceBuffer());
//            settleKernel.setArg<cl::Buffer>(6, settleParams->getDeviceBuffer());
//        }
//        settleKernel.setArg<float>(1, (float) tol);
//        context.executeKernel(settleKernel, settleAtoms->getSize());
//    }
//    if (shakeAtoms != NULL) {
//        if (!hasInitialized) {
//            shakeKernel.setArg<int>(0, shakeAtoms->getSize());
//            shakeKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
//            shakeKernel.setArg<cl::Buffer>(3, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta->getDeviceBuffer());
//            shakeKernel.setArg<cl::Buffer>(4, shakeAtoms->getDeviceBuffer());
//            shakeKernel.setArg<cl::Buffer>(5, shakeParams->getDeviceBuffer());
//        }
//        shakeKernel.setArg<float>(1, (float) tol);
//        context.executeKernel(shakeKernel, shakeAtoms->getSize());
//    }
//    if (ccmaAtoms != NULL) {
//        if (!hasInitialized) {
//            ccmaDirectionsKernel.setArg<cl::Buffer>(0, ccmaAtoms->getDeviceBuffer());
//            ccmaDirectionsKernel.setArg<cl::Buffer>(1, ccmaDistance->getDeviceBuffer());
//            ccmaDirectionsKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(0, ccmaAtoms->getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(1, ccmaDistance->getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(2, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta->getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(3, ccmaReducedMass->getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(4, ccmaDelta1->getDeviceBuffer());
//            ccmaForceKernel.setArg<cl::Buffer>(5, ccmaConverged->getDeviceBuffer());
//            ccmaMultiplyKernel.setArg<cl::Buffer>(0, ccmaDelta1->getDeviceBuffer());
//            ccmaMultiplyKernel.setArg<cl::Buffer>(1, ccmaDelta2->getDeviceBuffer());
//            ccmaMultiplyKernel.setArg<cl::Buffer>(2, ccmaConstraintMatrixColumn->getDeviceBuffer());
//            ccmaMultiplyKernel.setArg<cl::Buffer>(3, ccmaConstraintMatrixValue->getDeviceBuffer());
//            ccmaMultiplyKernel.setArg<cl::Buffer>(4, ccmaConverged->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(0, ccmaNumAtomConstraints->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(1, ccmaAtomConstraints->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(2, ccmaDistance->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(3, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(4, context.getVelm().getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(5, ccmaDelta1->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(6, ccmaDelta2->getDeviceBuffer());
//            ccmaUpdateKernel.setArg<cl::Buffer>(7, ccmaConverged->getDeviceBuffer());
//        }
//        ccmaForceKernel.setArg<float>(6, (float) tol);
//        context.executeKernel(ccmaDirectionsKernel, ccmaAtoms->getSize());
//        const int checkInterval = 4;
//        cl::Event event;
//        for (int i = 0; i < 150; i++) {
//            ccmaForceKernel.setArg<int>(7, i);
//            if (i == 0) {
//                ccmaConvergedMemory[0] = 1;
//                ccmaConvergedMemory[1] = 0;
//                context.getQueue().enqueueWriteBuffer(ccmaConverged->getDeviceBuffer(), CL_FALSE, 0, 2*sizeof(int), ccmaConvergedMemory);
//            }
//            context.executeKernel(ccmaForceKernel, ccmaAtoms->getSize());
//            if ((i+1)%checkInterval == 0)
//                context.getQueue().enqueueReadBuffer(ccmaConverged->getDeviceBuffer(), CL_FALSE, 0, 2*sizeof(int), ccmaConvergedMemory, NULL, &event);
//            ccmaMultiplyKernel.setArg<int>(5, i);
//            context.executeKernel(ccmaMultiplyKernel, ccmaAtoms->getSize());
//            ccmaUpdateKernel.setArg<int>(8, i);
//            context.executeKernel(ccmaUpdateKernel, context.getNumAtoms());
//            if ((i+1)%checkInterval == 0) {
//                event.wait();
//                if (ccmaConvergedMemory[i%2])
//                    break;
//            }
//        }
//    }
}

void CudaIntegrationUtilities::computeVirtualSites() {
//    if (numVsites > 0)
//        context.executeKernel(vsitePositionKernel, numVsites);
}

void CudaIntegrationUtilities::distributeForcesFromVirtualSites() {
//    if (numVsites > 0) {
//        vsiteForceKernel.setArg<cl::Buffer>(1, context.getForce().getDeviceBuffer());
//        context.executeKernel(vsiteForceKernel, numVsites);
//    }
}

void CudaIntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random != NULL) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("CudaIntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random = CudaArray::create<float4>(32*context.getPaddedNumAtoms(), "random");
    randomSeed = CudaArray::create<int4>(context.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
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

    // Create the kernel.

    CUmodule randomModule = context.createModule(CudaKernelSources::random);
    randomKernel = context.getKernel(randomModule, "generateRandomNumbers");
}

int CudaIntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random->getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    if (numValues > random->getSize()) {
        delete random;
        random = CudaArray::create<float4>(numValues, "random");
    }
    int size = random->getSize();
    void* args[] = {&size, &random->getDevicePointer(), &randomSeed->getDevicePointer()};
    context.executeKernel(randomKernel, args, random->getSize());
    randomPos = numValues;
    return 0;
}

void CudaIntegrationUtilities::createCheckpoint(ostream& stream) {
    stream.write((char*) &randomPos, sizeof(int));
    vector<float4> randomVec;
    random->download(randomVec);
    stream.write((char*) &randomVec[0], sizeof(float4)*random->getSize());
    vector<int4> randomSeedVec;
    randomSeed->download(randomSeedVec);
    stream.write((char*) &randomSeedVec[0], sizeof(int4)*randomSeed->getSize());
}

void CudaIntegrationUtilities::loadCheckpoint(istream& stream) {
    stream.read((char*) &randomPos, sizeof(int));
    vector<float4> randomVec(random->getSize());
    stream.read((char*) &randomVec[0], sizeof(float4)*random->getSize());
    random->upload(randomVec);
    vector<int4> randomSeedVec(randomSeed->getSize());
    stream.read((char*) &randomSeedVec[0], sizeof(int4)*randomSeed->getSize());
    randomSeed->upload(randomSeedVec);
}
