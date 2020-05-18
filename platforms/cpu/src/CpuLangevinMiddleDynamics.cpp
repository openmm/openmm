/* Portions copyright (c) 2006-2020 Stanford University and Simbios.
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
#include "CpuLangevinMiddleDynamics.h"

using namespace OpenMM;
using namespace std;

CpuLangevinMiddleDynamics::CpuLangevinMiddleDynamics(int numberOfAtoms, double deltaT, double friction, double temperature, ThreadPool& threads, CpuRandom& random) : 
           ReferenceLangevinMiddleDynamics(numberOfAtoms, deltaT, friction, temperature), threads(threads), random(random) {
}

CpuLangevinMiddleDynamics::~CpuLangevinMiddleDynamics() {
}

void CpuLangevinMiddleDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& velocities, vector<Vec3>& forces, vector<double>& inverseMasses) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->velocities = &velocities[0];
    this->forces = &forces[0];
    this->inverseMasses = &inverseMasses[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate1(threadIndex); });
    threads.waitForThreads();
}

void CpuLangevinMiddleDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
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

void CpuLangevinMiddleDynamics::updatePart3(ContextImpl& context, int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
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
}

void CpuLangevinMiddleDynamics::threadUpdate1(int threadIndex) {
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++)
        if (inverseMasses[i] != 0.0)
            velocities[i] += (getDeltaT()*inverseMasses[i])*forces[i];
}

void CpuLangevinMiddleDynamics::threadUpdate2(int threadIndex) {
    const double halfdt = 0.5*getDeltaT();
    const double kT = BOLTZ*getTemperature();
    const double friction = getFriction();
    const double vscale = exp(-getDeltaT()*friction);
    const double noisescale = sqrt(1-vscale*vscale);
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = atomCoordinates[i] + velocities[i]*halfdt;
            Vec3 noise(random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex));
            velocities[i] = vscale*velocities[i] + noisescale*sqrt(kT*inverseMasses[i])*noise;
            xPrime[i] = xPrime[i] + velocities[i]*halfdt;
            oldx[i] = xPrime[i];
        }
    }
}

void CpuLangevinMiddleDynamics::threadUpdate3(int threadIndex) {
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; ++i)
        if (inverseMasses[i] != 0.0) {
            velocities[i] += (xPrime[i]-oldx[i])/getDeltaT();
            atomCoordinates[i] = xPrime[i];
        }
}
