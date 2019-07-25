/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2015 Stanford University and the Authors.      *
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

#include "CpuSETTLE.h"
#include "openmm/internal/gmx_atomic.h"

using namespace OpenMM;
using namespace std;

class CpuSETTLE::ApplyToPositionsTask : public ThreadPool::Task {
public:
    ApplyToPositionsTask(vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& atomCoordinatesP, vector<RealOpenMM>& inverseMasses,
            RealOpenMM tolerance, vector<ReferenceSETTLEAlgorithm*>& threadSettle) : atomCoordinates(atomCoordinates), atomCoordinatesP(atomCoordinatesP),
            inverseMasses(inverseMasses), tolerance(tolerance), threadSettle(threadSettle) {
        gmx_atomic_set(&atomicCounter, 0);
    }
    void execute(ThreadPool& threads, int threadIndex) {
        while (true) {
            int index = gmx_atomic_fetch_add(&atomicCounter, 1);
            if (index >= threadSettle.size())
                break;
            threadSettle[index]->apply(atomCoordinates, atomCoordinatesP, inverseMasses, tolerance);
        }
    }
    vector<OpenMM::RealVec>& atomCoordinates;
    vector<OpenMM::RealVec>& atomCoordinatesP;
    vector<RealOpenMM>& inverseMasses;
    RealOpenMM tolerance;
    vector<ReferenceSETTLEAlgorithm*>& threadSettle;
    gmx_atomic_t atomicCounter;
};

class CpuSETTLE::ApplyToVelocitiesTask : public ThreadPool::Task {
public:
    ApplyToVelocitiesTask(vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& velocities, vector<RealOpenMM>& inverseMasses,
            RealOpenMM tolerance, vector<ReferenceSETTLEAlgorithm*>& threadSettle) : atomCoordinates(atomCoordinates), velocities(velocities),
            inverseMasses(inverseMasses), tolerance(tolerance), threadSettle(threadSettle) {
        gmx_atomic_set(&atomicCounter, 0);
    }
    void execute(ThreadPool& threads, int threadIndex) {
        while (true) {
            int index = gmx_atomic_fetch_add(&atomicCounter, 1);
            if (index >= threadSettle.size())
                break;
            threadSettle[index]->applyToVelocities(atomCoordinates, velocities, inverseMasses, tolerance);
        }
    }
    vector<OpenMM::RealVec>& atomCoordinates;
    vector<OpenMM::RealVec>& velocities;
    vector<RealOpenMM>& inverseMasses;
    RealOpenMM tolerance;
    vector<ReferenceSETTLEAlgorithm*>& threadSettle;
    gmx_atomic_t atomicCounter;
};

CpuSETTLE::CpuSETTLE(const System& system, const ReferenceSETTLEAlgorithm& settle, ThreadPool& threads) : threads(threads) {
    int numBlocks = 10*threads.getNumThreads();
    int numClusters = settle.getNumClusters();
    vector<RealOpenMM> mass(system.getNumParticles());
    for (int i = 0; i < system.getNumParticles(); i++)
        mass[i] = system.getParticleMass(i);
    for (int i = 0; i < numBlocks; i++) {
        int start = i*numClusters/numBlocks;
        int end = (i+1)*numClusters/numBlocks;
        if (start != end) {
            int numThreadClusters = end-start;
            vector<int> atom1(numThreadClusters), atom2(numThreadClusters), atom3(numThreadClusters);
            vector<RealOpenMM> distance1(numThreadClusters), distance2(numThreadClusters);
            for (int j = 0; j < numThreadClusters; j++)
                settle.getClusterParameters(start+j, atom1[j], atom2[j], atom3[j], distance1[j], distance2[j]);
            threadSettle.push_back(new ReferenceSETTLEAlgorithm(atom1, atom2, atom3, distance1, distance2, mass));
        }
    }
}

CpuSETTLE::~CpuSETTLE() {
    for (int i = 0; i < (int) threadSettle.size(); i++)
        delete threadSettle[i];
}

void CpuSETTLE::apply(vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& atomCoordinatesP, vector<RealOpenMM>& inverseMasses, RealOpenMM tolerance) {
    ApplyToPositionsTask task(atomCoordinates, atomCoordinatesP, inverseMasses, tolerance, threadSettle);
    threads.execute(task);
    threads.waitForThreads();
}

void CpuSETTLE::applyToVelocities(vector<OpenMM::RealVec>& atomCoordinates, vector<OpenMM::RealVec>& velocities, vector<RealOpenMM>& inverseMasses, RealOpenMM tolerance) {
    ApplyToVelocitiesTask task(atomCoordinates, velocities, inverseMasses, tolerance, threadSettle);
    threads.execute(task);
    threads.waitForThreads();
}
