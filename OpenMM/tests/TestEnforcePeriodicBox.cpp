/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
 * Authors: Robert McGibbon                                                   *
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
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testTruncatedOctahedron() {
    const int numMolecules = 50;
    const int numParticles = numMolecules*2;
    const float cutoff = 2.0;
    Vec3 a(6.7929, 0, 0);
    Vec3 b(-2.264163559406279, 6.404455775962287, 0);
    Vec3 c(-2.264163559406279, -3.2019384603140684, 5.54658849047036);
        
    System system;
    system.setDefaultPeriodicBoxVectors(a, b, c);
    NonbondedForce* force = new NonbondedForce();
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    
    force->setCutoffDistance(cutoff);
    force->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        force->addParticle(-1, 0.2, 0.2);
        force->addParticle(1, 0.2, 0.2);
        positions[2*i] = a*(5*genrand_real2(sfmt)-2) + b*(5*genrand_real2(sfmt)-2) + c*(5*genrand_real2(sfmt)-2);
        positions[2*i+1] = positions[2*i] + Vec3(1.0, 0.0, 0.0);
        system.addConstraint(2*i, 2*i+1, 1.0);
    }
    system.addForce(force);
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));
    context.setPositions(positions);
    State initialState = context.getState(State::Positions | State::Energy, true);
    for (int i = 0; i < numMolecules; i++) {
        Vec3 center = (initialState.getPositions()[2*i]+initialState.getPositions()[2*i+1])*0.5;
        ASSERT(center[0] >= 0.0);
        ASSERT(center[1] >= 0.0);
        ASSERT(center[2] >= 0.0);
        ASSERT(center[0] <= a[0]);
        ASSERT(center[1] <= b[1]);
        ASSERT(center[2] <= c[2]);
    }
    double initialEnergy = initialState.getPotentialEnergy();

    context.setState(initialState);
    State finalState = context.getState(State::Positions | State::Energy, true);
    double finalEnergy = finalState.getPotentialEnergy();

    ASSERT_EQUAL_TOL(initialEnergy, finalEnergy, 1e-4);
}

int main(int argc, char* argv[]) {
    try {
        testTruncatedOctahedron();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
