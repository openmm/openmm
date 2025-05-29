/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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
#include "openmm/CustomVolumeForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testVolume() {
    System system;
    Vec3 a(2, 0, 0);
    Vec3 b(0.1, 2, 0);
    Vec3 c(0.1, 0.1, 2);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    system.addParticle(1.0);
    CustomVolumeForce* force = new CustomVolumeForce("2*v");
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatform("Reference"));
    context.setPositions({Vec3()});
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < 10; i++) {
        a[0] = 2.0 + 0.1*genrand_real2(sfmt);
        b[1] = 2.0 + 0.1*genrand_real2(sfmt);
        c[2] = 2.0 + 0.1*genrand_real2(sfmt);
        b[0] = 0.1*genrand_real2(sfmt);
        c[0] = 0.1*genrand_real2(sfmt);
        c[1] = 0.1*genrand_real2(sfmt);
        context.setPeriodicBoxVectors(a, b, c);
        double energy = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL(2*a[0]*b[1]*c[2], energy, 1e-6);
    }
}

void testBoxVectors() {
    System system;
    Vec3 a(2, 0, 0);
    Vec3 b(0.1, 2, 0);
    Vec3 c(0.1, 0.1, 2);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    system.addParticle(1.0);
    CustomVolumeForce* force = new CustomVolumeForce("1 + ax + 2*bx + 3*by + 4*cx + 5*cy + 6*cz");
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatform("Reference"));
    context.setPositions({Vec3()});
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < 10; i++) {
        a[0] = 2.0 + 0.1*genrand_real2(sfmt);
        b[1] = 2.0 + 0.1*genrand_real2(sfmt);
        c[2] = 2.0 + 0.1*genrand_real2(sfmt);
        b[0] = 0.1*genrand_real2(sfmt);
        c[0] = 0.1*genrand_real2(sfmt);
        c[1] = 0.1*genrand_real2(sfmt);
        context.setPeriodicBoxVectors(a, b, c);
        double energy = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL(1+a[0]+2*b[0]+3*b[1]+4*c[0]+5*c[1]+6*c[2], energy, 1e-6);
    }
}

void testGlobalParameters() {
    System system;
    Vec3 a(2, 0, 0);
    Vec3 b(0, 2, 0);
    Vec3 c(0, 0, 3);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    system.addParticle(1.0);
    CustomVolumeForce* force = new CustomVolumeForce("p1*ax + p2*cz");
    force->addGlobalParameter("p1", 1.0);
    force->addGlobalParameter("p2", 1.0);
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatform("Reference"));
    context.setPositions({Vec3()});
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < 10; i++) {
        double p1 = genrand_real2(sfmt);
        double p2 = genrand_real2(sfmt);
        context.setParameter("p1", p1);
        context.setParameter("p2", p2);
        double energy = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL(2*p1+3*p2, energy, 1e-6);
    }
}

int main(int argc, char* argv[]) {
    try {
        testVolume();
        testBoxVectors();
        testGlobalParameters();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
