/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the Ewald summation method cuda implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "kernels/gputypes.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testEwaldPME() {

//      Use amorphous NaCl system for the tests

    const int numParticles = 216;
    const double cutoff = 0.8;
    const double boxSize = 1.86206;
    const double tol = 0.001;

    CudaPlatform cuda;
    ReferencePlatform reference;
    System system;
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setEwaldErrorTolerance(tol);

    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    system.setPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(nonbonded);

    vector<Vec3> positions(numParticles);
    #include "nacl_amorph.dat"

//    (1)  Check whether the Reference and Cuda platforms agree when using Ewald Method 
 
    Context cudaContext(system, integrator, cuda);
    Context referenceContext(system, integrator, reference);
    cudaContext.setPositions(positions);
    referenceContext.setPositions(positions);
    State cudaState = cudaContext.getState(State::Forces | State::Energy);
    State referenceState = referenceContext.getState(State::Forces | State::Energy);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(cudaState.getForces()[i], referenceState.getForces()[i], tol);
    }
    ASSERT_EQUAL_TOL(cudaState.getPotentialEnergy(), referenceState.getPotentialEnergy(), tol);

//    (2) Check whether Ewald method in Cuda is self-consistent

    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = cudaState.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }

    norm = std::sqrt(norm);
    const double delta = 1e-3;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = cudaState.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    Context cudaContext2(system, integrator, cuda);
    cudaContext2.setPositions(positions);
    
    State cudaState2 = cudaContext2.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (cudaState2.getPotentialEnergy()-cudaState.getPotentialEnergy())/delta, tol)

//    (3)  Check whether the Reference and Cuda platforms agree when using PME

    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    cudaContext.reinitialize();
    referenceContext.reinitialize();
    cudaContext.setPositions(positions);
    referenceContext.setPositions(positions);
    cudaState = cudaContext.getState(State::Forces | State::Energy);
    referenceState = referenceContext.getState(State::Forces | State::Energy);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(cudaState.getForces()[i], referenceState.getForces()[i], tol);
    }
    ASSERT_EQUAL_TOL(cudaState.getPotentialEnergy(), referenceState.getPotentialEnergy(), tol);

//    (4) Check whether PME method in Cuda is self-consistent

    norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = cudaState.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }

    norm = std::sqrt(norm);
    step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = cudaState.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    Context cudaContext3(system, integrator, cuda);
    cudaContext3.setPositions(positions);
     
    State cudaState3 = cudaContext3.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (cudaState3.getPotentialEnergy()-cudaState.getPotentialEnergy())/delta, tol)

}

void testEwald2Ions() {
    CudaPlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->addParticle(1.0, 1, 0);
    nonbonded->addParticle(-1.0, 1, 0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setEwaldErrorTolerance(TOL);
    system.setPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_VEC(Vec3(-123.711,  64.1877, -302.716), forces[0], 10*TOL);
    ASSERT_EQUAL_VEC(Vec3( 123.711, -64.1877,  302.716), forces[1], 10*TOL);
    ASSERT_EQUAL_TOL(-217.276, state.getPotentialEnergy(), 0.01/*10*TOL*/);
}


int main() {
    try {
     testEwaldPME();
//     testEwald2Ions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
