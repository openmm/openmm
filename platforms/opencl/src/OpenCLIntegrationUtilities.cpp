/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include "OpenCLArray.h"
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
        if (size == 3 || (size > 0 && dist != distance) || (size > 0 && invMass != peripheralInvMass))
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

OpenCLIntegrationUtilities::OpenCLIntegrationUtilities(OpenCLContext& context, const System& system) : context(context),
        posDelta(NULL), settleAtoms(NULL), settleParams(NULL), shakeAtoms(NULL), shakeParams(NULL),
        random(NULL), randomSeed(NULL), randomPos(NULL), stepSize(NULL) {
    // Create workspace arrays.

    posDelta = new OpenCLArray<mm_float4>(context, context.getPaddedNumAtoms(), "posDelta");
    stepSize = new OpenCLArray<mm_float2>(context, 1, "stepSize", true);
    stepSize->set(0, (mm_float2) {0.0f, 0.0f});
    stepSize->upload();

    // Create kernels for enforcing constraints.

    cl::Program settleProgram = context.createProgram(context.loadSourceFromFile("settle.cl"));
    settleKernel = cl::Kernel(settleProgram, "applySettle");
    cl::Program shakeProgram = context.createProgram(context.loadSourceFromFile("shakeHydrogens.cl"));
    shakeKernel = cl::Kernel(shakeProgram, "applyShakeToHydrogens");

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

    vector<map<int, float> > settleConstraints(system.getNumParticles());
    for (int i = 0; i < (int)atom1.size(); i++) {
        if (constraintCount[atom1[i]] == 2 && constraintCount[atom2[i]] == 2) {
            settleConstraints[atom1[i]][atom2[i]] = distance[i];
            settleConstraints[atom2[i]][atom1[i]] = distance[i];
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

    vector<bool> isShakeAtom(system.getNumParticles(), false);
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
                atoms.push_back((mm_int4) {atom1, atom2, atom3, 0});
                params.push_back((mm_float2) {dist12, dist23});
            }
            else if (dist12 == dist23) {
                // atom2 is the central atom
                atoms.push_back((mm_int4) {atom2, atom1, atom3, 0});
                params.push_back((mm_float2) {dist12, dist13});
            }
            else if (dist13 == dist23) {
                // atom3 is the central atom
                atoms.push_back((mm_int4) {atom3, atom1, atom2, 0});
                params.push_back((mm_float2) {dist13, dist12});
            }
            else
                throw OpenMMException("Two of the three distances constrained with SETTLE must be the same.");
            isShakeAtom[atom1] = true;
            isShakeAtom[atom2] = true;
            isShakeAtom[atom3] = true;
        }
        settleAtoms = new OpenCLArray<mm_int4>(context, atoms.size(), "settleAtoms");
        settleParams = new OpenCLArray<mm_float2>(context, params.size(), "settleParams");
        settleAtoms->upload(atoms);
        settleParams->upload(params);
    }

    // Find clusters consisting of a central atom with up to three peripheral atoms.

    map<int, ShakeCluster> clusters;
    vector<bool> invalidForShake(system.getNumParticles(), false);
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
            cluster.valid = !invalidForShake[cluster.centralID];
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
            atoms.push_back((mm_int4) {cluster.centralID, cluster.peripheralID[0], (cluster.size > 1 ? cluster.peripheralID[1] : -1), (cluster.size > 2 ? cluster.peripheralID[2] : -1)});
            params.push_back((mm_float4) {cluster.centralInvMass, 0.5f/(cluster.centralInvMass+cluster.peripheralInvMass), cluster.distance*cluster.distance, cluster.peripheralInvMass});
            isShakeAtom[cluster.centralID] = true;
            isShakeAtom[cluster.peripheralID[0]] = true;
            if (cluster.size > 1)
                isShakeAtom[cluster.peripheralID[1]] = true;
            if (cluster.size > 2)
                isShakeAtom[cluster.peripheralID[2]] = true;
            ++index;
        }
        shakeAtoms = new OpenCLArray<mm_int4>(context, atoms.size(), "shakeAtoms");
        shakeParams = new OpenCLArray<mm_float4>(context, params.size(), "shakeParams");
        shakeAtoms->upload(atoms);
        shakeParams->upload(params);
    }
}

OpenCLIntegrationUtilities::~OpenCLIntegrationUtilities() {
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
}

void OpenCLIntegrationUtilities::applyConstraints(double tol) {
    if (settleAtoms != NULL) {
        settleKernel.setArg<cl_int>(0, settleAtoms->getSize());
        settleKernel.setArg<cl_float>(1, tol);
        settleKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
        settleKernel.setArg<cl::Buffer>(3, posDelta->getDeviceBuffer());
        settleKernel.setArg<cl::Buffer>(4, posDelta->getDeviceBuffer());
        settleKernel.setArg<cl::Buffer>(5, context.getVelm().getDeviceBuffer());
        settleKernel.setArg<cl::Buffer>(6, settleAtoms->getDeviceBuffer());
        settleKernel.setArg<cl::Buffer>(7, settleParams->getDeviceBuffer());
        context.executeKernel(settleKernel, settleAtoms->getSize());
    }
    if (shakeAtoms != NULL) {
        shakeKernel.setArg<cl_int>(0, shakeAtoms->getSize());
        shakeKernel.setArg<cl_float>(1, tol);
        shakeKernel.setArg<cl::Buffer>(2, context.getPosq().getDeviceBuffer());
        shakeKernel.setArg<cl::Buffer>(3, posDelta->getDeviceBuffer());
        shakeKernel.setArg<cl::Buffer>(4, posDelta->getDeviceBuffer());
        shakeKernel.setArg<cl::Buffer>(5, shakeAtoms->getDeviceBuffer());
        shakeKernel.setArg<cl::Buffer>(6, shakeParams->getDeviceBuffer());
        context.executeKernel(shakeKernel, shakeAtoms->getSize());
    }
}

void OpenCLIntegrationUtilities::initRandomNumberGenerator(unsigned int randomNumberSeed) {
    if (random != NULL) {
        if (randomNumberSeed != lastSeed)
           throw OpenMMException("OpenCLIntegrationUtilities::initRandomNumberGenerator(): Requested two different values for the random number seed");
        return;
    }

    // Create the random number arrays.

    lastSeed = randomNumberSeed;
    random = new OpenCLArray<mm_float4>(context, 32*context.getPaddedNumAtoms(), "random");
    randomSeed = new OpenCLArray<mm_int4>(context, context.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize, "randomSeed");
    randomPos = random->getSize();

    // Use a quick and dirty RNG to pick seeds for the real random number generator.

    vector<mm_int4> seed(randomSeed->getSize());
    unsigned int r = randomNumberSeed;
    for (int i = 0; i < randomSeed->getSize(); i++) {
        seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
    }
    randomSeed->upload(seed);

    // Create the kernel.

    cl::Program randomProgram = context.createProgram(context.loadSourceFromFile("random.cl"));
    randomKernel = cl::Kernel(randomProgram, "generateRandomNumbers");
}

int OpenCLIntegrationUtilities::prepareRandomNumbers(int numValues) {
    if (randomPos+numValues <= random->getSize()) {
        int oldPos = randomPos;
        randomPos += numValues;
        return oldPos;
    }
    randomKernel.setArg<cl_int>(0, random->getSize());
    randomKernel.setArg<cl::Buffer>(1, random->getDeviceBuffer());
    randomKernel.setArg<cl::Buffer>(2, randomSeed->getDeviceBuffer());
    context.executeKernel(randomKernel, random->getSize());
    randomPos = numValues;
    return 0;
}
