/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <sys/time.h>
#include <iostream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

void testPerformance() {
    const int xsize = 20;
    const int ysize = 21;
    const int zsize = 21;
    const int numParticles = xsize*ysize*zsize;
    const double spacing = 0.3;
    CudaPlatform platform;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(xsize*spacing, 0, 0), Vec3(0, ysize*spacing, 0), Vec3(0, 0, zsize*spacing));
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    LangevinIntegrator integrator(1.0, 0.1, 0.001);
    NonbondedForce* nonbonded = new NonbondedForce();
    vector<Vec3> positions;
    vector<Vec3> velocities;
    double charge = 0.1;
    for (int i = 0; i < xsize; ++i)
        for (int j = 0; j < ysize; ++j)
            for (int k = 0; k < zsize; ++k) {
                nonbonded->addParticle(charge, 0.2, 0.1);
                charge = -charge;
                positions.push_back(Vec3(i*spacing, j*spacing, k*spacing));
                velocities.push_back(Vec3(0, 0, 0));
            }
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(3*spacing);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    timeval startTime;
    gettimeofday(&startTime, NULL);
    integrator.step(5000);
    State state = context.getState(State::Positions | State::Velocities | State::Forces | State::Energy);
    timeval endTime;
    gettimeofday(&endTime, NULL);
    double dt = endTime.tv_sec-startTime.tv_sec+1e-6*(endTime.tv_usec-startTime.tv_usec);
    std::cout << "Elapsed time: " << dt << std::endl;
    std::cout << "Final energy: " << state.getPotentialEnergy()+state.getKineticEnergy() << std::endl;
}

int main() {
    try {
        testPerformance();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


