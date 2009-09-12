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
 * This tests the Eewald summation method reference implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "openmm/HarmonicBondForce.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testLargeSystem() {
    ReferencePlatform platform;
    System system;
    for (int i = 0; i < 500; i++)
        system.addParticle(22.99);
    for (int i = 0; i < 500; i++)
        system.addParticle(35.45);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0; i < 500; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
//        nonbonded->addParticle(1.0, 0.33284,0.0115897);
    for (int i = 0; i < 500; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
//        nonbonded->addParticle(-1.0, 0.440104,0.4184);
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    const double cutoff = 1.0;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(2.82, 0, 0), Vec3(0, 2.82, 0), Vec3(0, 0, 2.82));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1000);
    #include "nacl_crystal.dat"
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
//    cout << "force 0: " << forces[0] << endl;
//    cout << "force 1: " << forces[1] << endl;
    cout << "PotentialEnergy: " << state.getPotentialEnergy() << endl;
    ASSERT_EQUAL_TOL(-430355.0, state.getPotentialEnergy(), 100*TOL);
//    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], 10*TOL);
//    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], 10*TOL);
}

void testEwald() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->addParticle(1.0, 1, 0);
    nonbonded->addParticle(-1.0, 1, 0);
// Sodium Chloride
//    nonbonded->addParticle(1.0, 0.33284,0.0115897);
//    nonbonded->addParticle(-1.0, 0.440104,0.4184);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
//    cout << "force 0: " << forces[0] << endl;
//    cout << "force 1: " << forces[1] << endl;
//    cout << "PotentialEnergy: " << state.getPotentialEnergy() << endl;
    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], 10*TOL);
    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], 10*TOL);
}

void testPME() {
    ReferencePlatform platform;
    System system;
    for (int i = 0 ; i < 42 ; i++)
    {
       system.addParticle(1.0);
    }
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0 ; i < 14 ; i++)
    {
      nonbonded->addParticle(-0.82, 1, 0);
      nonbonded->addParticle(0.41, 1, 0);
      nonbonded->addParticle(0.41, 1, 0);
    }
    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    const double cutoff = 0.8;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(1.86206, 0, 0), Vec3(0, 1.86206, 0), Vec3(0, 0, 1.86206));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(42);
    #include "water.dat"
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
   // for (int i = 0 ; i < 42 ; i++)
   //   cout << "f [" << i << " : ]" << forces[i] << endl;
   // cout << "PotentialEnergy: " << state.getPotentialEnergy() << endl;
//    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], 10*TOL);
//    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], 10*TOL);
}


int main() {
    try {
     testLargeSystem();
     testEwald();
     testPME();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
