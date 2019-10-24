
/* Portions copyright (c) 2006-2019 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "CpuBAOABDynamics.h"

using namespace OpenMM;
using namespace std;

CpuBAOABDynamics::CpuBAOABDynamics(int numberOfAtoms, double deltaT, double friction, double temperature, ThreadPool& threads, CpuRandom& random) : 
           ReferenceBAOABDynamics(numberOfAtoms, deltaT, friction, temperature), threads(threads), random(random) {
}

CpuBAOABDynamics::~CpuBAOABDynamics() {
}

void CpuBAOABDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
                                   vector<Vec3>& forces, vector<double>& inverseMasses, vector<Vec3>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->forces = &forces[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate1(threadIndex); });
    threads.waitForThreads();
}

void CpuBAOABDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
                                   vector<double>& inverseMasses, vector<Vec3>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate2(threadIndex); });
    threads.waitForThreads();
}

void CpuBAOABDynamics::updatePart3(ContextImpl& context, int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
                                   vector<double>& inverseMasses, vector<Vec3>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate3(threadIndex); });
    threads.waitForThreads();
    context.calcForcesAndEnergy(true, false);
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate4(threadIndex); });
    threads.waitForThreads();
}

void CpuBAOABDynamics::threadUpdate1(int threadIndex) {
    const double halfdt = 0.5*getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (halfdt*inverseMasses[i])*forces[i];
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void CpuBAOABDynamics::threadUpdate2(int threadIndex) {
    const double halfdt = 0.5*getDeltaT();
    const double kT = BOLTZ*getTemperature();
    const double friction = getFriction();
    const double vscale = exp(-getDeltaT()*friction);
    const double noisescale = sqrt(1-vscale*vscale);
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/halfdt;
            Vec3 noise(random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex));
            velocities[i] = vscale*velocities[i] + noisescale*sqrt(kT*inverseMasses[i])*noise;
            atomCoordinates[i] = xPrime[i];
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void CpuBAOABDynamics::threadUpdate3(int threadIndex) {
    const double halfdt = 0.5*getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; ++i)
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/halfdt;
            atomCoordinates[i] = xPrime[i];
        }
}

void CpuBAOABDynamics::threadUpdate4(int threadIndex) {
    const double halfdt = 0.5*getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; ++i)
        if (inverseMasses[i] != 0.0)
            velocities[i] += (halfdt*inverseMasses[i])*forces[i];
}
