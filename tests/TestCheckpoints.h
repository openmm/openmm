/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2015 Stanford University and the Authors.      *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/AndersenThermostat.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void compareStates(State& s1, State& s2) {
    ASSERT_EQUAL_TOL(s1.getTime(), s2.getTime(), TOL);
    int numParticles = s1.getPositions().size();
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(s1.getPositions()[i], s2.getPositions()[i], TOL);
        ASSERT_EQUAL_VEC(s1.getVelocities()[i], s2.getVelocities()[i], TOL);
        Vec3 a1, b1, c1, a2, b2, c2;
        s1.getPeriodicBoxVectors(a1, b1, c1);
        s2.getPeriodicBoxVectors(a2, b2, c2);
        ASSERT_EQUAL_VEC(a1, a2, TOL);
        ASSERT_EQUAL_VEC(b1, b2, TOL);
        ASSERT_EQUAL_VEC(c1, c2, TOL);
        for (map<string, double>::const_iterator iter = s1.getParameters().begin(); iter != s1.getParameters().end(); ++iter)
            ASSERT_EQUAL(iter->second, (*s2.getParameters().find(iter->first)).second);
    }
}

void testSetState() {
    const int numParticles = 10;
    const double boxSize = 3.0;
    const double temperature = 200.0;
    System system;
    system.addForce(new AndersenThermostat(0.0, 100.0));
    NonbondedForce* nonbonded = new NonbondedForce();
    system.addForce(nonbonded);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(i%2 == 0 ? 0.1 : -0.1, 0.2, 0.1);
        positions[i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
    }
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    context.setParameter(AndersenThermostat::Temperature(), temperature);
    
    // Run for a little while.
    
    integrator.step(100);
    
    // Record the current state.
    
    State s1 = context.getState(State::Positions | State::Velocities | State::Parameters);
    
    // Continue the simulation for a few more steps and record a partial state.
    
    integrator.step(10);
    State s2 = context.getState(State::Positions);
    
    // Restore the original state and see if everything gets restored correctly.
    
    context.setPeriodicBoxVectors(Vec3(2*boxSize, 0, 0), Vec3(0, 2*boxSize, 0), Vec3(0, 0, 2*boxSize));
    context.setParameter(AndersenThermostat::Temperature(), temperature+10);
    context.setState(s1);
    State s3 = context.getState(State::Positions | State::Velocities | State::Parameters);
    compareStates(s1, s3);
    
    // Set the partial state and see if the correct things were set.
    
    context.setState(s2);
    State s4 = context.getState(State::Positions | State::Velocities | State::Parameters);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(s2.getPositions()[i], s4.getPositions()[i], TOL);
        ASSERT_EQUAL_VEC(s1.getVelocities()[i], s4.getVelocities()[i], TOL);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSetState();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
