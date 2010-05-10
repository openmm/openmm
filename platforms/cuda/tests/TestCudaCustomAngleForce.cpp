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
 * This tests the CUDA implementation of CustomAngleForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/CustomAngleForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testAngles() {
    CudaPlatform platform;

    // Create a system using a CustomAngleForce.

    System customSystem;
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    CustomAngleForce* custom = new CustomAngleForce("scale*k*(theta-theta0)^2");
    custom->addPerAngleParameter("theta0");
    custom->addPerAngleParameter("k");
    custom->addGlobalParameter("scale", 0.5);
    vector<double> parameters(2);
    parameters[0] = 1.5;
    parameters[1] = 0.8;
    custom->addAngle(0, 1, 2, parameters);
    parameters[0] = 2.0;
    parameters[1] = 0.5;
    custom->addAngle(1, 2, 3, parameters);
    customSystem.addForce(custom);

    // Create an identical system using a HarmonicAngleForce.

    System harmonicSystem;
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    HarmonicAngleForce* harmonic = new HarmonicAngleForce();
    harmonic->addAngle(0, 1, 2, 1.5, 0.8);
    harmonic->addAngle(1, 2, 3, 2.0, 0.5);
    harmonicSystem.addForce(harmonic);

    // Set the atoms in various positions, and verify that both systems give identical forces and energy.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<Vec3> positions(4);
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        double energy1, energy2;
        vector<Vec3> forces1, forces2;
        {
            Context c(customSystem, integrator1, platform);
            c.setPositions(positions);
            State s = c.getState(State::Forces | State::Energy);
            energy1 = s.getPotentialEnergy();
            forces1 = s.getForces();
        }
        {
            Context c(harmonicSystem, integrator1, platform);
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

int main() {
    try {
        testAngles();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}



