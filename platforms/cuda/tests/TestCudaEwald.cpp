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
 * This tests all the different force terms in the reference implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "OpenMMContext.h"
#include "CudaPlatform.h"
#include "ReferencePlatform.h"
#include "HarmonicBondForce.h"
#include "NonbondedForce.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include "VerletIntegrator.h"
#include "internal/OpenMMContextImpl.h"
#include "kernels/gputypes.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testEwald() {
    CudaPlatform platform;
    System system(2, 0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce(2, 0);
    nonbonded->setParticleParameters(0, 1.0, 1, 0);
    nonbonded->setParticleParameters(1, -1.0, 1, 0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    system.addForce(nonbonded);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    cout << "    force 0: " << forces[0] << endl;
    cout << "    force 1: " << forces[1] << endl;
    cout << "    energyPoten: " << state.getPotentialEnergy() << endl;

    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], TOL);

    //ASSERT_EQUAL_TOL(2*138.935485*(1.0)*(1.0+krf*1.0-crf), state.getPotentialEnergy(), TOL);
}

void testPeriodic() {
    CudaPlatform platform;
    System system(3, 0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    HarmonicBondForce* bonds = new HarmonicBondForce(1);
    bonds->setBondParameters(0, 0, 1, 1, 0);
    system.addForce(bonds);
    NonbondedForce* nonbonded = new NonbondedForce(3, 0);
    nonbonded->setParticleParameters(0, 1.0, 1, 0);
    nonbonded->setParticleParameters(1, 1.0, 1, 0);
    nonbonded->setParticleParameters(2, 1.0, 1, 0);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    system.addForce(nonbonded);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    positions[2] = Vec3(3, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    cout << ""  << endl;
    cout << "    force 0: " << forces[0] << endl;
    cout << "    force 1: " << forces[1] << endl;
    cout << "    force 2: " << forces[2] << endl;
    cout << "energyPoten: " << state.getPotentialEnergy() << endl;
    cout << ""  << endl;
    cout << "    no cutoff force: 138.935485"  << endl;
    cout << ""  << endl;
    const double eps = 78.3;
    const double krf = (1.0/(cutoff*cutoff*cutoff))*(eps-1.0)/(2.0*eps+1.0);
    const double crf = (1.0/cutoff)*(3.0*eps)/(2.0*eps+1.0);
    const double force = 138.935485*(1.0)*(1.0-2.0*krf*1.0);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(2*138.935485*(1.0)*(1.0+krf*1.0-crf), state.getPotentialEnergy(), TOL);
}



int main() {
    try {
    cout << ""  << endl;
    cout << "Executing Periodic..." << endl;
        testPeriodic();
    cout << ""  << endl;
    cout << "Executing Ewald..." << endl;
        testEwald();
    cout << ""  << endl;
    cout << "Done" << endl;
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
