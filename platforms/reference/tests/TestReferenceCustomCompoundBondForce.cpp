/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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
 * This tests the reference implementation of CustomCompoundBondForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBond() {
    ReferencePlatform platform;

    // Create a system using a CustomCompoundBondForce.

    System customSystem;
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    CustomCompoundBondForce* custom = new CustomCompoundBondForce(4, "0.5*kb*((distance(p1,p2)-b0)^2+(distance(p2,p3)-b0)^2)+0.5*ka*(angle(p2,p3,p4)-a0)^2+kt*(1+cos(dihedral(p1,p2,p3,p4)-t0))");
    custom->addPerBondParameter("kb");
    custom->addPerBondParameter("ka");
    custom->addPerBondParameter("kt");
    custom->addPerBondParameter("b0");
    custom->addPerBondParameter("a0");
    custom->addPerBondParameter("t0");
    vector<int> particles(4);
    particles[0] = 0;
    particles[1] = 1;
    particles[2] = 3;
    particles[3] = 2;
    vector<double> parameters(6);
    parameters[0] = 1.5;
    parameters[1] = 0.8;
    parameters[2] = 0.6;
    parameters[3] = 1.1;
    parameters[4] = 2.9;
    parameters[5] = 1.3;
    custom->addBond(particles, parameters);
    customSystem.addForce(custom);

    // Create an identical system using standard forces.

    System standardSystem;
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.1, 1.5);
    bonds->addBond(1, 3, 1.1, 1.5);
    standardSystem.addForce(bonds);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    angles->addAngle(1, 3, 2, 2.9, 0.8);
    standardSystem.addForce(angles);
    PeriodicTorsionForce* torsions = new PeriodicTorsionForce();
    torsions->addTorsion(0, 1, 3, 2, 1, 1.3, 0.6);
    standardSystem.addForce(torsions);

    // Set the atoms in various positions, and verify that both systems give identical forces and energy.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<Vec3> positions(4);
    for (int i = 0; i < 10; i++) {
        VerletIntegrator integrator1(0.01);
        VerletIntegrator integrator2(0.01);
        Context c1(customSystem, integrator1, platform);
        Context c2(standardSystem, integrator2, platform);
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        c1.setPositions(positions);
        c2.setPositions(positions);
        State s1 = c1.getState(State::Forces | State::Energy);
        State s2 = c2.getState(State::Forces | State::Energy);
        for (int i = 0; i < customSystem.getNumParticles(); i++)
            ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], TOL);
        ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), s2.getPotentialEnergy(), TOL);
    }
}

void testPositionDependence() {
    ReferencePlatform platform;
    System customSystem;
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    CustomCompoundBondForce* custom = new CustomCompoundBondForce(2, "scale1*distance(p1,p2)+scale2*x1+2*y2");
    custom->addGlobalParameter("scale1", 0.3);
    custom->addGlobalParameter("scale2", 0.2);
    vector<int> particles(2);
    particles[0] = 0;
    particles[1] = 1;
    vector<double> parameters;
    custom->addBond(particles, parameters);
    customSystem.addForce(custom);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0.5, 1, 0);
    positions[1] = Vec3(1.5, 1, 0);
    VerletIntegrator integrator(0.01);
    Context context(customSystem, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(0.3*1.0+0.2*0.5+2*1, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(Vec3(0.3-0.2, 0, 0), state.getForces()[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(-0.3, -2, 0), state.getForces()[1], 1e-5);
}

int main() {
    try {
        testBond();
        testPositionDependence();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


