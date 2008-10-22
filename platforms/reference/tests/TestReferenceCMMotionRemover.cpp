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
 * This tests the reference implementation of AndersenThermostat.
 */

#include "../../../tests/AssertionUtilities.h"
#include "CMMotionRemover.h"
#include "OpenMMContext.h"
#include "ReferencePlatform.h"
#include "HarmonicBondForce.h"
#include "NonbondedForce.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

Vec3 calcCM(const vector<Vec3>& values, System& system) {
    Vec3 cm;
    for (int j = 0; j < system.getNumAtoms(); ++j) {
        cm[0] += values[j][0]*system.getAtomMass(j);
        cm[1] += values[j][1]*system.getAtomMass(j);
        cm[2] += values[j][2]*system.getAtomMass(j);
    }
    return cm;
}

void testMotionRemoval() {
    const int numAtoms = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    ReferencePlatform platform;
    System system(numAtoms, 0);
    VerletIntegrator integrator(0.01);
    HarmonicBondForce* bonds = new HarmonicBondForce(1);
    bonds->setBondParameters(0, 2, 3, 2.0, 0.5);
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce(numAtoms, 0);
    for (int i = 0; i < numAtoms; ++i) {
        system.setAtomMass(i, i+1);
        nonbonded->setAtomParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(nonbonded);
    CMMotionRemover* remover = new CMMotionRemover();
    system.addForce(remover);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numAtoms);
    vector<Vec3> velocities(numAtoms);
    init_gen_rand(0);
    for (int i = 0; i < numAtoms; ++i) {
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
        velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
    }
    context.setPositions(positions);
    context.setVelocities(velocities);
    
    // Now run it for a while and see if the center of mass remains fixed.
    
    Vec3 cmPos = calcCM(context.getState(State::Positions).getPositions(), system);
    for (int i = 0; i < 1000; ++i) {
        integrator.step(1);
        State state = context.getState(State::Positions | State::Velocities);
        Vec3 pos = calcCM(state.getPositions(), system);
        ASSERT_EQUAL_VEC(cmPos, pos, 1e-2);
        Vec3 vel = calcCM(state.getVelocities(), system);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), vel, 1e-2);
    }
}

int main() {
    try {
        testMotionRemoval();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
