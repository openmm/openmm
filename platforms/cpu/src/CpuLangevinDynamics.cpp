
/* Portions copyright (c) 2006-2017 Stanford University and Simbios.
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
#include "CpuLangevinDynamics.h"

using namespace OpenMM;
using namespace std;

CpuLangevinDynamics::CpuLangevinDynamics(int numberOfAtoms, double deltaT, double friction, double temperature, ThreadPool& threads, CpuRandom& random) : 
           ReferenceStochasticDynamics(numberOfAtoms, deltaT, friction, temperature), threads(threads), random(random) {
}

CpuLangevinDynamics::~CpuLangevinDynamics() {
}

void CpuLangevinDynamics::updatePart1(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
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

void CpuLangevinDynamics::updatePart2(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
                                      vector<Vec3>& forces, vector<double>& inverseMasses, vector<Vec3>& xPrime) {
    // Record the parameters for the threads.
    
    this->numberOfAtoms = numberOfAtoms;
    this->atomCoordinates = &atomCoordinates[0];
    this->velocities = &velocities[0];
    this->forces = &forces[0];
    this->inverseMasses = &inverseMasses[0];
    this->xPrime = &xPrime[0];
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadUpdate2(threadIndex); });
    threads.waitForThreads();
}

void CpuLangevinDynamics::updatePart3(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<Vec3>& velocities,
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

void CpuLangevinDynamics::threadUpdate1(int threadIndex) {
    double dt = getDeltaT();
    double friction = getFriction();
    const double vscale = exp(-dt*friction);
    const double fscale = (friction == 0 ? dt : (1-vscale)/friction);
    const double kT = BOLTZ*getTemperature();
    const double noisescale = sqrt(kT*(1-vscale*vscale));
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            double sqrtInvMass = sqrt(inverseMasses[i]);
            Vec3 noise(random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex), random.getGaussianRandom(threadIndex));
            velocities[i]  = velocities[i]*vscale + forces[i]*(fscale*inverseMasses[i]) + noise*(noisescale*sqrtInvMass);
        }
   }
}

void CpuLangevinDynamics::threadUpdate2(int threadIndex) {
    const double dt = getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

    for (int i = start; i < end; i++) {
        if (inverseMasses[i] != 0.0) {
            xPrime[i] = atomCoordinates[i]+velocities[i]*dt;
        }
   }
}

void CpuLangevinDynamics::threadUpdate3(int threadIndex) {
   const double invStepSize = 1.0/getDeltaT();
    int start = threadIndex*numberOfAtoms/threads.getNumThreads();
    int end = (threadIndex+1)*numberOfAtoms/threads.getNumThreads();

   for (int i = start; i < end; ++i)
       if (inverseMasses[i] != 0.0) {
            velocities[i] = (xPrime[i]-atomCoordinates[i])*invStepSize;
            atomCoordinates[i] = xPrime[i];
       }
}

