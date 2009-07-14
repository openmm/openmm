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

void testEwald() {
    ReferencePlatform platform;
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
    nonbonded->setPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    cout << "force 0: " << forces[0] << endl;
    cout << "force 1: " << forces[1] << endl;
    cout << "PotentialEnergy: " << state.getPotentialEnergy() << endl;
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
    nonbonded->setPeriodicBoxVectors(Vec3(1.86206, 0, 0), Vec3(0, 1.86206, 0), Vec3(0, 0, 1.86206));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(42);
positions[0] = Vec3( 0.23,0.628,0.113);
positions[1] = Vec3(0.137,0.626, 0.15);
positions[2] = Vec3(0.231,0.589,0.021);
positions[3] = Vec3(-0.307,-0.351,0.703);
positions[4] = Vec3(-0.364,-0.367,0.784);
positions[5] = Vec3(-0.366,-0.341,0.623);
positions[6] = Vec3(-0.569,-0.634,-0.439);
positions[7] = Vec3(-0.532,-0.707,-0.497);
positions[8] = Vec3(-0.517,-0.629,-0.354);
positions[9] = Vec3(-0.871, 0.41,-0.62);
positions[10] = Vec3(-0.948,0.444,-0.566);
positions[11] = Vec3(-0.905,0.359,-0.699);
positions[12] = Vec3(0.249,-0.077,-0.621);
positions[13] = Vec3(0.306,-0.142,-0.571);
positions[14] = Vec3(0.233,-0.11,-0.714);
positions[15] = Vec3(0.561,0.222,-0.715);
positions[16] = Vec3(0.599,0.138,-0.678);
positions[17] = Vec3(0.473,0.241,-0.671);
positions[18] = Vec3(-0.515,-0.803,-0.628);
positions[19] = Vec3(-0.491,-0.866,-0.702);
positions[20] = Vec3(-0.605,-0.763,-0.646);
positions[21] = Vec3(-0.021,0.175,-0.899);
positions[22] = Vec3(0.018, 0.09,-0.935);
positions[23] = Vec3(-0.119,0.177,-0.918);
positions[24] = Vec3(-0.422,0.856,-0.464);
positions[25] = Vec3(-0.479,0.908,-0.527);
positions[26] = Vec3(-0.326,0.868,-0.488);
positions[27] = Vec3(-0.369,-0.095,-0.903);
positions[28] = Vec3(-0.336,-0.031,-0.972);
positions[29] = Vec3(-0.303,-0.101,-0.828);
positions[30] = Vec3(0.594,0.745,0.652);
positions[31] = Vec3(0.644, 0.83,0.633);
positions[32] = Vec3(0.506,0.747,0.604);
positions[33] = Vec3(-0.157,-0.375,-0.758);
positions[34] = Vec3(-0.25, -0.4,-0.785);
positions[35] = Vec3(-0.131,-0.425,-0.676);
positions[36] = Vec3(0.618,-0.295,-0.578);
positions[37] = Vec3(0.613,-0.213,-0.521);
positions[38] = Vec3(0.707,-0.298,-0.623);
positions[39] = Vec3(0.039,-0.785,  0.3);
positions[40] = Vec3(0.138,-0.796,0.291);
positions[41] = Vec3(-0.001,-0.871,0.332);
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
