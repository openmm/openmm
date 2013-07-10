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

/**
 * This tests the reference implementation of CustomExternalForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testForce() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomExternalForce* forceField = new CustomExternalForce("scale*(x+yscale*(y-y0)^2)");
    forceField->addPerParticleParameter("y0");
    forceField->addPerParticleParameter("yscale");
    forceField->addGlobalParameter("scale", 0.5);
    vector<double> parameters(2);
    parameters[0] = 0.5;
    parameters[1] = 2.0;
    forceField->addParticle(0, parameters);
    parameters[0] = 1.5;
    parameters[1] = 3.0;
    forceField->addParticle(2, parameters);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 1);
    positions[2] = Vec3(1, 0, 1);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(-0.5, -0.5*2.0*2.0*1.5, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5, 0.5*3.0*2.0*1.5, 0), forces[2], TOL);
        ASSERT_EQUAL_TOL(0.5*(1.0 + 2.0*1.5*1.5 + 3.0*1.5*1.5), state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the parameters and make sure it's still correct.
    
    parameters[0] = 1.4;
    parameters[1] = 3.5;
    forceField->setParticleParameters(1, 2, parameters);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(-0.5, -0.5*2.0*2.0*1.5, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
        ASSERT_EQUAL_VEC(Vec3(-0.5, 0.5*3.5*2.0*1.4, 0), forces[2], TOL);
        ASSERT_EQUAL_TOL(0.5*(1.0 + 2.0*1.5*1.5 + 3.5*1.4*1.4), state.getPotentialEnergy(), TOL);
    }
}

int main() {
    try {
        testForce();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


