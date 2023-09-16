/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2023 Stanford University and the Authors.           *
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
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

class TestForce : public Force {
protected:
    ForceImpl* createImpl() const;
};

class TestForceImpl : public CustomCPPForceImpl {
public:
    const TestForce& owner;
    TestForceImpl(const TestForce& owner) : CustomCPPForceImpl(owner), owner(owner) {
    }
    double computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
        for (int i = 0; i < forces.size(); i++)
            forces[i] = Vec3(i, i+1, i+2);
        return 10.5;
    }
    const TestForce& getOwner() const {
        return owner;
    }
};

ForceImpl* TestForce::createImpl() const {
    return new TestForceImpl(*this);
}

void testForce() {
    int numParticles = 3;
    System system;
    NonbondedForce* nb = new NonbondedForce();
    TestForce* test = new TestForce();
    nb->setForceGroup(0);
    test->setForceGroup(1);
    system.addForce(test);
    system.addForce(nb);
    vector<Vec3> positions;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nb->addParticle(0.0, 0.5, 1.0);
        positions.push_back(Vec3(0.6*i, 0, 0));
    }
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    State state2 = context.getState(State::Forces | State::Energy, false, 1<<0);
    State state3 = context.getState(State::Forces | State::Energy, false, 1<<1);
    ASSERT_EQUAL_TOL(10.5, state3.getPotentialEnergy(), 1e-6);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy()+state3.getPotentialEnergy(), 1e-6);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(Vec3(i, i+1, i+2), state3.getForces()[i], 1e-6);
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i]+state3.getForces()[i], 1e-6);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testForce();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
