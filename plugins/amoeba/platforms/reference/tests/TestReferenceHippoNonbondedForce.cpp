/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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
 * This tests the Reference implementation of HippoNonbondedForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerAmoebaReferenceKernelFactories();

void buildWaterSystem(System& system, int numWaters) {
    HippoNonbondedForce* hippo = new HippoNonbondedForce();
    system.addForce(hippo);
    for (int mol = 0; mol < numWaters; mol++) {
        system.addParticle(15.995);
        system.addParticle(1.008);
        system.addParticle(1.008);
        double bohr = 0.52917720859;
        double ds = 0.1*bohr;
        double qs = 0.01*bohr*bohr/3.0;
        double c6s = sqrt(4.184)*0.001;
        double ps = sqrt(4.184*0.1);
        hippo->addParticle(-0.38280, {0.0, 0.0, ds*0.05477}, {qs*0.69866, 0.0, 0.0, 0.0, qs*-0.60471, 0.0, 0.0, 0.0, qs*-0.09395}, 6.0,
                    10*4.7075, 4.184*1326.0, 10*40.0, c6s*18.7737, ps*2.7104, -2.4233, 10*4.3097,
                    0.001*0.795, HippoNonbondedForce::Bisector, 3*mol+1, 3*mol+2, -1);
        hippo->addParticle(0.19140, {0.0, 0.0, ds*-0.20097}, {qs*0.03881, 0.0, 0.0, 0.0, qs*0.02214, 0.0, 0.0, 0.0, qs*-0.06095}, 1.0,
                    10*4.7909, 0.0, 10*3.5582, c6s*4.5670, ps*2.0037, -0.8086, 10*4.6450,
                    0.001*0.341, HippoNonbondedForce::ZThenX, 3*mol, 3*mol+2, -1);
        hippo->addParticle(0.19140, {0.0, 0.0, ds*-0.20097}, {qs*0.03881, 0.0, 0.0, 0.0, qs*0.02214, 0.0, 0.0, 0.0, qs*-0.06095}, 1.0,
                    10*4.7909, 0.0, 10*3.5582, c6s*4.5670, ps*2.0037, -0.8086, 10*4.6450,
                    0.001*0.341, HippoNonbondedForce::ZThenX, 3*mol, 3*mol+1, -1);
        hippo->addException(3*mol, 3*mol+1, 0.0, 0.0, 0.2, 0.0, 0.0);
        hippo->addException(3*mol, 3*mol+2, 0.0, 0.0, 0.2, 0.0, 0.0);
        hippo->addException(3*mol+1, 3*mol+2, 0.0, 0.0, 1.0, 0.0, 0.0);
    }
}

void checkForceEnergyConsistency(Context& context) {
    // Find the direction and magnitude of the force.

    State state = context.getState(State::Positions | State::Forces | State::Energy);
    int numParticles = context.getSystem().getNumParticles();
    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }
    norm = std::sqrt(norm);

    // Take a small step in the direction of the force and see whether the potential energy changes by the expected amount.

    const double delta = 1e-3;
    double step = 0.5*delta/norm;
    vector<Vec3> positions2(numParticles), positions3(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = state.getPositions()[i];
        Vec3 f = state.getForces()[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(state3.getPotentialEnergy()-state2.getPotentialEnergy(), norm*delta, 1e-4)
}

void testWaterDimer() {
    System system;
    buildWaterSystem(system, 2);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));
    vector<Vec3> positions = {
        0.1*Vec3(1.505434, 0.0, -0.065656),
        0.1*Vec3(0.553912, 0.0, 0.057710),
        0.1*Vec3(1.907155, 0.0, 0.801980),
        0.1*Vec3(-1.436029, 0.0, 0.060505),
        0.1*Vec3(-1.781197, 0.772272, -0.388976),
        0.1*Vec3(-1.781197, -0.772272, -0.388976)
    };
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    
    // Compare the forces and energy to reference values computed with Tinker.

    ASSERT_EQUAL_TOL(-18.393623712669680, state.getPotentialEnergy(), 1e-5);
    vector<Vec3> expectedForces = {
        Vec3(-162.94090034728887, 0.0, 35.06615691195519),
        Vec3(127.50063696213348, 0.0, -46.51857483822334),
        Vec3(39.59601328153432, 0.0, 11.805509637931072),
        Vec3(-73.52341534248339, 0.0, -92.08855312751808),
        Vec3(34.68383272305204, -26.35219958830841, 45.867730707927564),
        Vec3(34.68383272305204, 26.35219958830841, 45.867730707927564),
    };
    for (int i = 0; i < system.getNumParticles(); i++)
        ASSERT_EQUAL_VEC(expectedForces[i], state.getForces()[i], 1e-5);
    checkForceEnergyConsistency(context);
}

int main(int numberOfArguments, char* argv[]) {

    try {
        registerAmoebaReferenceKernelFactories();
        testWaterDimer();
    }
    catch (const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
