/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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
 * This tests the CUDA implementation of CustomTorsionForce.
 */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/CustomTorsionForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testTorsions() {
    CudaPlatform platform;

    // Create a system using a CustomTorsionForce.

    System customSystem;
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    CustomTorsionForce* custom = new CustomTorsionForce("k*(1+cos(n*theta-theta0))");
    custom->addPerTorsionParameter("theta0");
    custom->addPerTorsionParameter("n");
    custom->addGlobalParameter("k", 0.5);
    vector<double> parameters(2);
    parameters[0] = 1.5;
    parameters[1] = 1;
    custom->addTorsion(0, 1, 2, 3, parameters);
    parameters[0] = 2.0;
    parameters[1] = 2;
    custom->addTorsion(1, 2, 3, 4, parameters);
    customSystem.addForce(custom);

    // Create an identical system using a PeriodicTorsionForce.

    System harmonicSystem;
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    PeriodicTorsionForce* periodic = new PeriodicTorsionForce();
    periodic->addTorsion(0, 1, 2, 3, 1, 1.5, 0.5);
    periodic->addTorsion(1, 2, 3, 4, 2, 2.0, 0.5);
    harmonicSystem.addForce(periodic);

    // Set the atoms in various positions, and verify that both systems give identical forces and energy.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<Vec3> positions(5);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        double energy1, energy2;
        vector<Vec3> forces1, forces2;
        {
            VerletIntegrator integrator(0.01);
            Context c(customSystem, integrator, platform);
            c.setPositions(positions);
            State s = c.getState(State::Forces | State::Energy);
            energy1 = s.getPotentialEnergy();
            forces1 = s.getForces();
        }
        {
            VerletIntegrator integrator(0.01);
            Context c(harmonicSystem, integrator, platform);
            c.setPositions(positions);
            State s = c.getState(State::Forces | State::Energy);
            energy2 = s.getPotentialEnergy();
            forces2 = s.getForces();
        }
        for (int i = 0; i < customSystem.getNumParticles(); i++)
            ASSERT_EQUAL_VEC(forces2[i], forces1[i], TOL);
        ASSERT_EQUAL_TOL(energy2, energy1, TOL);
    }
}

void testRange() {
    CudaPlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomTorsionForce* custom = new CustomTorsionForce("theta");
    custom->addTorsion(0, 1, 2, 3, vector<double>());
    system.addForce(custom);

    // Set the atoms in various positions, and verify that the angle is always in the expected range.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(4);
    VerletIntegrator integrator(0.01);
    double minAngle = 1000;
    double maxAngle = -1000;
    Context context(system, integrator, platform);
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        context.setPositions(positions);
        double angle = context.getState(State::Energy).getPotentialEnergy();
        if (angle < minAngle)
            minAngle = angle;
        if (angle > maxAngle)
            maxAngle = angle;
    }
    ASSERT(minAngle >= -M_PI);
    ASSERT(maxAngle <= M_PI);
}

int main() {
    try {
        testTorsions();
        testRange();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

