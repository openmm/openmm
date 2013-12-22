
/* Portions copyright (c) 2006-2013 Stanford University and Simbios.
 * Authors: Peter Eastman
 * Contributors: 
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuLangevinDynamics.h"

using namespace OpenMM;
using namespace std;

class CpuLangevinDynamics::Update1Task : public ThreadPool::Task {
public:
    Update1Task(CpuLangevinDynamics& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadUpdate1(threadIndex);
    }
    CpuLangevinDynamics& owner;
};

class CpuLangevinDynamics::Update2Task : public ThreadPool::Task {
public:
    Update2Task(CpuLangevinDynamics& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadUpdate2(threadIndex);
    }
    CpuLangevinDynamics& owner;
};

CpuLangevinDynamics::CpuLangevinDynamics(int numberOfAtoms, RealOpenMM deltaT, RealOpenMM tau, RealOpenMM temperature, ThreadPool& threads, CpuRandom& random) : 
           ReferenceStochasticDynamics(numberOfAtoms, deltaT, tau, temperature), threads(threads), random(random) {
}

CpuLangevinDynamics::~CpuLangevinDynamics() {
}

void CpuLangevinDynamics::updatePart1(int numberOfAtoms, vector<RealVec>& atomCoordinates, vector<RealVec>& velocities,
                                      vector<RealVec>& forces, vector<RealOpenMM>& inverseMasses, vector<RealVec>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->forces = &forces[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    Update1Task task(*this);
    threads.execute(task);
    threads.waitForThreads();
}

void CpuLangevinDynamics::updatePart2(int numberOfAtoms, vector<RealVec>& atomCoordinates, vector<RealVec>& velocities,
                                      vector<RealVec>& forces, vector<RealOpenMM>& inverseMasses, vector<RealVec>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->forces = &forces[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    Update2Task task(*this);
    threads.execute(task);
    threads.waitForThreads();
}

void CpuLangevinDynamics::threadUpdate1(int threadIndex) {
    const RealOpenMM tau = getTau();
    const RealOpenMM vscale = EXP(-getDeltaT()/tau);
    const RealOpenMM fscale = (1-vscale)*tau;
    const RealOpenMM kT = BOLTZ*getTemperature();
    const RealOpenMM noisescale = SQRT(2*kT/tau)*SQRT(0.5*(1-vscale*vscale)*tau);
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            RealOpenMM sqrtInvMass = SQRT(inverseMasses[i]);
            RealVec noise(random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex));
            velocities[i]  = velocities[i]*vscale + forces[i]*(fscale*inverseMasses[i]) + noise*(noisescale*sqrtInvMass);
        }
   }
}

void CpuLangevinDynamics::threadUpdate2(int threadIndex) {
    const RealOpenMM dt = getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            RealOpenMM sqrtInvMass = SQRT(inverseMasses[i]);
            xPrime[i] = atomCoordinates[i]+velocities[i]*dt;
        }
   }
}
