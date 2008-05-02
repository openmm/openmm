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
#include "AndersenThermostat.h"
#include "OpenMMContext.h"
#include "ReferencePlatform.h"
#include "StandardMMForceField.h"
#include "System.h"
#include "VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testTemperature() {
    const int numAtoms = 8;
    const double temp = 100.0;
    const double collisionFreq = 10.0;
    ReferencePlatform platform;
    System system(numAtoms, 0);
    VerletIntegrator integrator(0.01);
    StandardMMForceField* forceField = new StandardMMForceField(numAtoms, 0, 0, 0, 0);
    for (int i = 0; i < numAtoms; ++i) {
        system.setAtomMass(i, 2.0);
        forceField->setAtomParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    AndersenThermostat* thermstat = new AndersenThermostat(temp, collisionFreq);
    system.addForce(thermstat);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numAtoms);
    for (int i = 0; i < numAtoms; ++i)
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
    context.setPositions(positions);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double ke = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(1);
    }
    ke /= 1000;
    double expected = 0.5*numAtoms*3*BOLTZ*temp;
    ASSERT_EQUAL_TOL(expected, ke, 3*expected/std::sqrt(1000.0));
}

int main() {
    try {
        testTemperature();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
