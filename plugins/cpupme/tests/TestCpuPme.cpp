/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/**
 * This tests the CPU implementation of PME.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "../src/CpuPmeKernels.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

class IO : public CalcPmeReciprocalForceKernel::IO {
public:
    vector<float> posq;
    float* force;
    float* getPosq() {
        return &posq[0];
    }
    void setForce(float* force) {
        this->force = force;
    }
};

void testPME() {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 5.0;
    const double cutoff = 1.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::PME);
    force->setCutoffDistance(cutoff);
    force->setReciprocalSpaceForceGroup(1);
    force->setEwaldErrorTolerance(1e-4);
    
    // Compute the reciprocal space forces with the reference platform.
    
    Platform& platform = Platform::getPlatformByName("Reference");
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State refState = context.getState(State::Forces | State::Energy, false, 1<<1);
    
    // Now compute them with the optimized kernel.
    
    double alpha;
    int gridx, gridy, gridz;
    NonbondedForceImpl::calcPMEParameters(system, *force, alpha, gridx, gridy, gridz);
    CpuCalcPmeReciprocalForceKernel pme(CalcPmeReciprocalForceKernel::Name(), platform);
    IO io;
    double sumSquaredCharges = 0;
    for (int i = 0; i < numParticles; i++) {
        io.posq.push_back(positions[i][0]);
        io.posq.push_back(positions[i][1]);
        io.posq.push_back(positions[i][2]);
        double charge, sigma, epsilon;
        force->getParticleParameters(i, charge, sigma, epsilon);
        io.posq.push_back(charge);
        sumSquaredCharges += charge*charge;
    }
    double ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI);
    pme.initialize(gridx, gridy, gridz, numParticles, alpha);
    pme.beginComputation(io, Vec3(boxWidth, boxWidth, boxWidth), true);
    double energy = pme.finishComputation(io);
    
    // See if they match.
    
    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), energy+ewaldSelfEnergy, 1e-3);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(refState.getForces()[i], Vec3(io.force[4*i], io.force[4*i+1], io.force[4*i+2]), 1e-3);
}

int main(int argc, char* argv[]) {
    try {
        if (!CpuCalcPmeReciprocalForceKernel::isProcessorSupported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 0;
        }
        testPME();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
