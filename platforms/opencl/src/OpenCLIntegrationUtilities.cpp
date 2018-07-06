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

#include "OpenCLIntegrationUtilities.h"
#include "OpenCLKernelSources.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/VirtualSite.h"
#include "quern.h"
#include "OpenCLExpressionUtilities.h"
#include "ReferenceCCMAAlgorithm.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>

using namespace OpenMM;
using namespace std;

struct OpenCLIntegrationUtilities::ShakeCluster {
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

struct OpenCLIntegrationUtilities::ConstraintOrderer : public binary_function<int, int, bool> {
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

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseMixedPrecision())
        kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
    else
        kernel.setArg<void*>(index, NULL);
}

OpenCLIntegrationUtilities::OpenCLIntegrationUtilities(OpenCLContext& context, const System& system) : context(context),
        randomPos(0), hasInitializedPosConstraintKernels(false), hasInitializedVelConstraintKernels(false), hasOverlappingVsites(false) {
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
    
    // Create the time shift kernel for calculating kinetic energy.
    
    map<string, string> timeShiftDefines;
    timeShiftDefines["NUM_ATOMS"] = context.intToString(system.getNumParticles());
    cl::Program utilitiesProgram = context.createProgram(OpenCLKernelSources::integrationUtilities, timeShiftDefines);
    timeShiftKernel = cl::Kernel(utilitiesProgram, "timeShiftVelocities");

    // Create kernels for enforcing constraints.

    map<string, string> velocityDefines;
    velocityDefines["CONSTRAIN_VELOCITIES"] = "1";
    cl::Program settleProgram = context.createProgram(OpenCLKernelSources::settle);
    settlePosKernel = cl::Kernel(settleProgram, "applySettle");
    settleVelKernel = cl::Kernel(settleProgram, "constrainVelocities");
    cl::Program shakeProgram = context.createProgram(OpenCLKernelSources::shakeHydrogens);
    shakePosKernel = cl::Kernel(shakeProgram, "applyShakeToHydrogens");
    shakeProgram = context.createProgram(OpenCLKernelSources::shakeHydrogens, velocityDefines);
    shakeVelKernel = cl::Kernel(shakeProgram, "applyShakeToHydrogens");

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
            params.push_back(mm_float4((cl_float) cluster.centralInvMass, (cl_float) (0.5/(cluster.centralInvMass+cluster.peripheralInvMass)), (cl_float) (cluster.distance*cluster.distance), (cl_float) cluster.peripheralInvMass));
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
            refMasses[i] = (double) system.getParticleMass(i);

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
        ccmaAtomConstraints.initialize<cl_int>(context, numAtoms*maxAtomConstraints, "CcmaAtomConstraints");
        ccmaNumAtomConstraints.initialize<cl_int>(context, numAtoms, "CcmaAtomConstraintsIndex");
        ccmaConstraintMatrixColumn.initialize<cl_int>(context, numCCMA*maxRowElements, "ConstraintMatrixColumn");
        ccmaConverged.initialize<cl_int>(context, 2, "CcmaConverged");
        ccmaConvergedHostBuffer.initialize<cl_int>(context, 1, "CcmaConvergedHostBuffer", CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR);
        // Different communication mechanisms give optimal performance on AMD and on NVIDIA.
        string vendor = context.getDevice().getInfo<CL_DEVICE_VENDOR>();
        ccmaUseDirectBuffer = (vendor.size() >= 28 && vendor.substr(0, 28) == "Advanced Micro Devices, Inc.");
        vector<mm_int2> atomsVec(ccmaAtoms.getSize());
        vector<cl_int> atomConstraintsVec(ccmaAtomConstraints.getSize());
        vector<cl_int> numAtomConstraintsVec(ccmaNumAtomConstraints.getSize());
        vector<cl_int> constraintMatrixColumnVec(ccmaConstraintMatrixColumn.getSize());
        int elementSize = (context.getUseDoublePrecision() || context.getUseMixedPrecision() ? sizeof(cl_double) : sizeof(cl_float));
        ccmaDistance.initialize(context, numCCMA, 4*elementSize, "CcmaDistance");
        ccmaDelta1.initialize(context, numCCMA, elementSize, "CcmaDelta1");
        ccmaDelta2.initialize(context, numCCMA, elementSize, "CcmaDelta2");
        ccmaReducedMass.initialize(context, numCCMA, elementSize, "CcmaReducedMass");
        ccmaConstraintMatrixValue.initialize(context, numCCMA*maxRowElements, elementSize, "ConstraintMatrixValue");
        vector<mm_double4> distanceVec(ccmaDistance.getSize());
        vector<cl_double> reducedMassVec(ccmaReducedMass.getSize());
        vector<cl_double> constraintMatrixValueVec(ccmaConstraintMatrixValue.getSize());
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
        for (unsigned int i = 0; i < atomConstraints.size(); i++) {
            numAtomConstraintsVec[i] = atomConstraints[i].size();
            for (unsigned int j = 0; j < atomConstraints[i].size(); j++) {
                bool forward = (atom1[ccmaConstraints[atomConstraints[i][j]]] == i);
                atomConstraintsVec[i+j*numAtoms] = (forward ? inverseOrder[atomConstraints[i][j]]+1 : -inverseOrder[atomConstraints[i][j]]-1);
            }
        }
        ccmaDistance.upload(distanceVec, true, true);
        ccmaReducedMass.upload(reducedMassVec, true, true);
        ccmaConstraintMatrixValue.upload(constraintMatrixValueVec, true, true);
        ccmaAtoms.upload(atomsVec);
        ccmaAtomConstraints.upload(atomConstraintsVec);
        ccmaNumAtomConstraints.upload(numAtomConstraintsVec);
        ccmaConstraintMatrixColumn.upload(constraintMatrixColumnVec);

        // Create the CCMA kernels.

        map<string, string> defines;
        defines["NUM_CONSTRAINTS"] = context.intToString(numCCMA);
        defines["NUM_ATOMS"] = context.intToString(numAtoms);
        cl::Program ccmaProgram = context.createProgram(OpenCLKernelSources::ccma, defines);
        ccmaDirectionsKernel = cl::Kernel(ccmaProgram, "computeConstraintDirections");
        ccmaPosForceKernel = cl::Kernel(ccmaProgram, "computeConstraintForce");
        ccmaMultiplyKernel = cl::Kernel(ccmaProgram, "multiplyByConstraintMatrix");
        ccmaPosUpdateKernel = cl::Kernel(ccmaProgram, "updateAtomPositions");
        defines["CONSTRAIN_VELOCITIES"] = "1";
        ccmaProgram = context.createProgram(OpenCLKernelSources::ccma, defines);
        ccmaVelForceKernel = cl::Kernel(ccmaProgram, "computeConstraintForce");
        ccmaVelUpdateKernel = cl::Kernel(ccmaProgram, "updateAtomPositions");
    }
    
    // Build the list of virtual sites.
    
    vector<mm_int4> vsite2AvgAtomVec;
    vector<mm_double2> vsite2AvgWeightVec;
    vector<mm_int4> vsite3AvgAtomVec;
    vector<mm_double4> vsite3AvgWeightVec;
    vector<mm_int4> vsiteOutOfPlaneAtomVec;
    vector<mm_double4> vsiteOutOfPlaneWeightVec;
    vector<cl_int> vsiteLocalCoordsIndexVec;
    vector<cl_int> vsiteLocalCoordsAtomVec;
    vector<cl_int> vsiteLocalCoordsStartVec;
    vector<cl_double> vsiteLocalCoordsWeightVec;
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
    vsiteLocalCoordsIndex.initialize<cl_int>(context, max(1, (int) vsiteLocalCoordsIndexVec.size()), "vsiteLocalCoordsIndex");
    vsiteLocalCoordsAtoms.initialize<cl_int>(context, max(1, (int) vsiteLocalCoordsAtomVec.size()), "vsiteLocalCoordsAtoms");
    vsiteLocalCoordsStartIndex.initialize<cl_int>(context, max(1, (int) vsiteLocalCoordsStartVec.size()), "vsiteLocalCoordsStartIndex");
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
    int elementSize = (context.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
    vsite2AvgWeights.initialize(context, max(1, num2Avg), 2*elementSize, "vsite2AvgWeights");
    vsite3AvgWeights.initialize(context, max(1, num3Avg), 4*elementSize, "vsite3AvgWeights");
    vsiteOutOfPlaneWeights.initialize(context, max(1, numOutOfPlane), 4*elementSize, "vsiteOutOfPlaneWeights");
    vsiteLocalCoordsWeights.initialize(context, max(1, (int) vsiteLocalCoordsWeightVec.size()), elementSize, "vsiteLocalCoordsWeights");
    vsiteLocalCoordsPos.initialize(context, max(1, (int) vsiteLocalCoordsPosVec.size()), 4*elementSize, "vsiteLocalCoordsPos");
    if (num2Avg > 0)
        vsite2AvgWeights.upload(vsite2AvgWeightVec, true, true);
    if (num3Avg > 0)
        vsite3AvgWeights.upload(vsite3AvgWeightVec, true, true);
    if (numOutOfPlane > 0)
        vsiteOutOfPlaneWeights.upload(vsiteOutOfPlaneWeightVec, true, true);
    if (numLocalCoords > 0) {
        vsiteLocalCoordsWeights.upload(vsiteLocalCoordsWeightVec, true, true);
        vsiteLocalCoordsPos.upload(vsiteLocalCoordsPosVec, true, true);
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
    if (hasOverlappingVsites && context.getUseDoublePrecision() && !context.getSupports64BitGlobalAtomics())
        throw OpenMMException("This device does not support 64 bit atomics.  Cannot use double precision when multiple virtual sites depend on the same atom.");
    
    // Create the kernels for virtual sites.

    map<string, string> defines;
    defines["NUM_2_AVERAGE"] = context.intToString(num2Avg);
    defines["NUM_3_AVERAGE"] = context.intToString(num3Avg);
    defines["NUM_OUT_OF_PLANE"] = context.intToString(numOutOfPlane);
    defines["NUM_LOCAL_COORDS"] = context.intToString(numLocalCoords);
    defines["NUM_ATOMS"] = context.intToString(numAtoms);
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    if (hasOverlappingVsites)
        defines["HAS_OVERLAPPING_VSITES"] = "1";
    cl::Program vsiteProgram = context.createProgram(OpenCLKernelSources::virtualSites, defines);
    vsitePositionKernel = cl::Kernel(vsiteProgram, "computeVirtualSites");
    int index = 0;
    vsitePositionKernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
    if (context.getUseMixedPrecision())
        vsitePositionKernel.setArg<cl::Buffer>(index++, context.getPosqCorrection().getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsite2AvgAtoms.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsite2AvgWeights.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsite3AvgAtoms.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsite3AvgWeights.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteOutOfPlaneAtoms.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteOutOfPlaneWeights.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsIndex.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsAtoms.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsWeights.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsPos.getDeviceBuffer());
    vsitePositionKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsStartIndex.getDeviceBuffer());
    vsiteForceKernel = cl::Kernel(vsiteProgram, "distributeForces");
    index = 0;
    vsiteForceKernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
    index++; // Skip argument 1: the force array hasn't been created yet.
    if (context.getSupports64BitGlobalAtomics())
        index++; // Skip argument 2: the force array hasn't been created yet.
    if (context.getUseMixedPrecision())
        vsiteForceKernel.setArg<cl::Buffer>(index++, context.getPosqCorrection().getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsite2AvgAtoms.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsite2AvgWeights.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsite3AvgAtoms.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsite3AvgWeights.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteOutOfPlaneAtoms.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteOutOfPlaneWeights.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsIndex.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsAtoms.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsWeights.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsPos.getDeviceBuffer());
    vsiteForceKernel.setArg<cl::Buffer>(index++, vsiteLocalCoordsStartIndex.getDeviceBuffer());
    if (hasOverlappingVsites && context.getSupports64BitGlobalAtomics())
        vsiteAddForcesKernel = cl::Kernel(vsiteProgram, "addDistributedForces");
}

void OpenCLIntegrationUtilities::setNextStepSize(double size) {
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

double OpenCLIntegrationUtilities::getLastStepSize() {
    if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
        stepSize.download(&lastStepSize);
    else {
        mm_float2 lastStepSizeFloat;
        stepSize.download(&lastStepSizeFloat);
        lastStepSize = mm_double2(lastStepSizeFloat.x, lastStepSizeFloat.y);
    }
    return lastStepSize.y;
}

void OpenCLIntegrationUtilities::applyConstraints(double tol) {
    applyConstraints(false, tol);
}

void OpenCLIntegrationUtilities::applyVelocityConstraints(double tol) {
    applyConstraints(true, tol);
}

void OpenCLIntegrationUtilities::applyConstraints(bool constrainVelocities, double tol) {
    bool hasInitialized;
    cl::Kernel settleKernel, shakeKernel, ccmaForceKernel, ccmaUpdateKernel;
    if (constrainVelocities) {
        hasInitialized = hasInitializedVelConstraintKernels;
        settleKernel = settleVelKernel;
        shakeKernel = shakeVelKernel;
        ccmaForceKernel = ccmaVelForceKernel;
        ccmaUpdateKernel = ccmaVelUpdateKernel;
        hasInitializedVelConstraintKernels = true;
    }
    else {
        hasInitialized = hasInitializedPosConstraintKernels;
        settleKernel = settlePosKernel;
        shakeKernel = shakePosKernel;
        ccmaForceKernel = ccmaPosForceKernel;
        ccmaUpdateKernel = ccmaPosUpdateKernel;
        hasInitializedPosConstraintKernels = true;
    }
    if (settleAtoms.isInitialized()) {
        if (!hasInitialized) {
            settleKernel.setArg<cl_int>(0, settleAtoms.getSize());
            settleKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
            if (context.getUseMixedPrecision())
                settleKernel.setArg<cl::Buffer>(3, context.getPosqCorrection().getDeviceBuffer());
            else
                settleKernel.setArg<void*>(3, NULL);
            settleKernel.setArg<cl::Buffer>(4, posDelta.getDeviceBuffer());
            settleKernel.setArg<cl::Buffer>(5, context.getVelm().getDeviceBuffer());
            settleKernel.setArg<cl::Buffer>(6, settleAtoms.getDeviceBuffer());
            settleKernel.setArg<cl::Buffer>(7, settleParams.getDeviceBuffer());
        }
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            settleKernel.setArg<cl_double>(1, (cl_double) tol);
        else
            settleKernel.setArg<cl_float>(1, (cl_float) tol);
        context.executeKernel(settleKernel, settleAtoms.getSize());
    }
    if (shakeAtoms.isInitialized()) {
        if (!hasInitialized) {
            shakeKernel.setArg<cl_int>(0, shakeAtoms.getSize());
            shakeKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
            if (context.getUseMixedPrecision())
                shakeKernel.setArg<cl::Buffer>(3, context.getPosqCorrection().getDeviceBuffer());
            else
                shakeKernel.setArg<void*>(3, NULL);
            shakeKernel.setArg<cl::Buffer>(4, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta.getDeviceBuffer());
            shakeKernel.setArg<cl::Buffer>(5, shakeAtoms.getDeviceBuffer());
            shakeKernel.setArg<cl::Buffer>(6, shakeParams.getDeviceBuffer());
        }
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            shakeKernel.setArg<cl_double>(1, (cl_double) tol);
        else
            shakeKernel.setArg<cl_float>(1, (cl_float) tol);
        context.executeKernel(shakeKernel, shakeAtoms.getSize());
    }
    if (ccmaAtoms.isInitialized()) {
        if (!hasInitialized) {
            ccmaDirectionsKernel.setArg<cl::Buffer>(0, ccmaAtoms.getDeviceBuffer());
            ccmaDirectionsKernel.setArg<cl::Buffer>(1, ccmaDistance.getDeviceBuffer());
            ccmaDirectionsKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
            if (context.getUseMixedPrecision())
                ccmaDirectionsKernel.setArg<cl::Buffer>(3, context.getPosqCorrection().getDeviceBuffer());
            else
                ccmaDirectionsKernel.setArg<void*>(3, NULL);
            ccmaDirectionsKernel.setArg<cl::Buffer>(4, ccmaConverged.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(0, ccmaAtoms.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(1, ccmaDistance.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(2, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(3, ccmaReducedMass.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(4, ccmaDelta1.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(5, ccmaConverged.getDeviceBuffer());
            ccmaForceKernel.setArg<cl::Buffer>(6, ccmaConvergedHostBuffer.getDeviceBuffer());
            ccmaMultiplyKernel.setArg<cl::Buffer>(0, ccmaDelta1.getDeviceBuffer());
            ccmaMultiplyKernel.setArg<cl::Buffer>(1, ccmaDelta2.getDeviceBuffer());
            ccmaMultiplyKernel.setArg<cl::Buffer>(2, ccmaConstraintMatrixColumn.getDeviceBuffer());
            ccmaMultiplyKernel.setArg<cl::Buffer>(3, ccmaConstraintMatrixValue.getDeviceBuffer());
            ccmaMultiplyKernel.setArg<cl::Buffer>(4, ccmaConverged.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(0, ccmaNumAtomConstraints.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(1, ccmaAtomConstraints.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(2, ccmaDistance.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(3, constrainVelocities ? context.getVelm().getDeviceBuffer() : posDelta.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(4, context.getVelm().getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(5, ccmaDelta1.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(6, ccmaDelta2.getDeviceBuffer());
            ccmaUpdateKernel.setArg<cl::Buffer>(7, ccmaConverged.getDeviceBuffer());
        }
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            ccmaForceKernel.setArg<cl_double>(7, (cl_double) tol);
        else
            ccmaForceKernel.setArg<cl_float>(7, (cl_float) tol);
        context.executeKernel(ccmaDirectionsKernel, ccmaAtoms.getSize());
        const int checkInterval = 4;
        int* converged = (int*) context.getPinnedBuffer();
        int* ccmaConvergedHostMemory = (int*) context.getQueue().enqueueMapBuffer(ccmaConvergedHostBuffer.getDeviceBuffer(), CL_TRUE, CL_MAP_WRITE, 0, sizeof(cl_int));
        ccmaConvergedHostMemory[0] = 0;
        context.getQueue().enqueueUnmapMemObject(ccmaConvergedHostBuffer.getDeviceBuffer(), ccmaConvergedHostMemory);
        for (int i = 0; i < 150; i++) {
            ccmaForceKernel.setArg<cl_int>(8, i);
            context.executeKernel(ccmaForceKernel, ccmaAtoms.getSize());
            cl::Event event;
            if ((i+1)%checkInterval == 0 && !ccmaUseDirectBuffer)
                context.getQueue().enqueueReadBuffer(ccmaConverged.getDeviceBuffer(), CL_FALSE, 0, 2*sizeof(cl_int), converged, NULL, &event);
            ccmaMultiplyKernel.setArg<cl_int>(5, i);
            context.executeKernel(ccmaMultiplyKernel, ccmaAtoms.getSize());
            ccmaUpdateKernel.setArg<cl_int>(8, i);
            context.executeKernel(ccmaUpdateKernel, context.getNumAtoms());
            if ((i+1)%checkInterval == 0) {
                if (ccmaUseDirectBuffer) {
                    ccmaConvergedHostMemory = (int*) context.getQueue().enqueueMapBuffer(ccmaConvergedHostBuffer.getDeviceBuffer(), CL_FALSE, CL_MAP_READ, 0, sizeof(cl_int), NULL, &event);
                    context.getQueue().flush();
                    while (event.getInfo<CL_EVENT_COMMAND_EXECUTION_STATUS>() != CL_COMPLETE)
                        ;
                    converged[i%2] = ccmaConvergedHostMemory[0];
                    context.getQueue().enqueueUnmapMemObject(ccmaConvergedHostBuffer.getDeviceBuffer(), ccmaConvergedHostMemory);
                }
                else
                    event.wait();
                if (converged[i%2])
                    break;
            }
        }
    }
}

void OpenCLIntegrationUtilities::computeVirtualSites() {
    if (numVsites > 0)
        context.executeKernel(vsitePositionKernel, numVsites);
}

void OpenCLIntegrationUtilities::distributeForcesFromVirtualSites() {
    if (numVsites > 0) {
        // Set arguments that didn't exist yet in the constructor.
        
        vsiteForceKernel.setArg<cl::Buffer>(1, context.getForce().getDeviceBuffer());
        if (context.getSupports64BitGlobalAtomics()) {
            vsiteForceKernel.setArg<cl::Buffer>(2, context.getLongForceBuffer().getDeviceBuffer());
            if (hasOverlappingVsites) {
                // We'll be using 64 bit atomics for the force redistribution, so clear the buffer.
                
                context.clearBuffer(context.getLongForceBuffer());
            }
        }
        context.executeKernel(vsiteForceKernel, numVsites);
        if (context.getSupports64BitGlobalAtomics() && hasOverlappingVsites) {
            // Add the redistributed forces from the virtual sites to the main force array.
            
            vsiteAddForcesKernel.setArg<cl::Buffer>(0, context.getLongForceBuffer().getDeviceBuffer());
            vsiteAddForcesKernel.setArg<cl::Buffer>(1, context.getForce().getDeviceBuffer());
            context.executeKernel(vsiteAddForcesKernel, context.getNumAtoms());
        }
    }
}

void OpenCLIntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random.isInitialized()) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("OpenCLIntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random.initialize<mm_float4>(context, 4*context.getPaddedNumAtoms(), "random");
    randomSeed.initialize<mm_int4>(context, context.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize, "randomSeed");
    randomPos = random.getSize();

    // Use a quick and dirty RNG to pick seeds for the real random number generator.

    vector<mm_int4> seed(randomSeed.getSize());
    unsigned int r = randomNumberSeed;
    // A seed of 0 means use a unique one
    if (r == 0) r = (unsigned int) osrngseed();
    for (int i = 0; i < randomSeed.getSize(); i++) {
        seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
    }
    randomSeed.upload(seed);

    // Create the kernel.

    cl::Program randomProgram = context.createProgram(OpenCLKernelSources::random);
    randomKernel = cl::Kernel(randomProgram, "generateRandomNumbers");
}

int OpenCLIntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random.getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    if (numValues > random.getSize()) {
        random.resize(numValues);
    }
    randomKernel.setArg<cl_int>(0, random.getSize());
    randomKernel.setArg<cl::Buffer>(1, random.getDeviceBuffer());
    randomKernel.setArg<cl::Buffer>(2, randomSeed.getDeviceBuffer());
    context.executeKernel(randomKernel, random.getSize());
    randomPos = numValues;
    return 0;
}

void OpenCLIntegrationUtilities::createCheckpoint(ostream& stream) {
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

void OpenCLIntegrationUtilities::loadCheckpoint(istream& stream) {
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

double OpenCLIntegrationUtilities::computeKineticEnergy(double timeShift) {
    int numParticles = context.getNumAtoms();
    if (timeShift != 0) {
        // Copy the velocities into the posDelta array while we temporarily modify them.

        context.getVelm().copyTo(posDelta);

        // Apply the time shift.

        timeShiftKernel.setArg<cl::Buffer>(0, context.getVelm().getDeviceBuffer());
        timeShiftKernel.setArg<cl::Buffer>(1, context.getForce().getDeviceBuffer());
        if (context.getUseDoublePrecision())
            timeShiftKernel.setArg<cl_double>(2, timeShift);
        else
            timeShiftKernel.setArg<cl_float>(2, (cl_float) timeShift);
        context.executeKernel(timeShiftKernel, numParticles);
        applyConstraints(true, 1e-4);
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
