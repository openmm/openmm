/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
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
    ASSERT(!custom->usesPeriodicBoundaryConditions());
    ASSERT(!customSystem.usesPeriodicBoundaryConditions());

    // Create an identical system using a HarmonicAngleForce.

    System harmonicSystem;
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
    harmonicSystem.addParticle(1.0);
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
    Context c1(customSystem, integrator1, platform);
    Context c2(harmonicSystem, integrator2, platform);
    for (int i = 0; i < 10; i++) {
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
    
    // Try changing the angle parameters and make sure it's still correct.
    
    parameters[0] = 1.6;
    parameters[1] = 0.9;
    custom->setAngleParameters(0, 0, 1, 2, parameters);
    parameters[0] = 2.1;
    parameters[1] = 0.6;
    custom->setAngleParameters(1, 1, 2, 3, parameters);
    custom->updateParametersInContext(c1);
    harmonic->setAngleParameters(0, 0, 1, 2, 1.6, 0.9);
    harmonic->setAngleParameters(1, 1, 2, 3, 2.1, 0.6);
    harmonic->updateParametersInContext(c2);
    {
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

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomAngleForce* force = new CustomAngleForce("theta+none");
    force->addAngle(0, 1, 2);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    bool threwException = false;
    try {
        Context(system, integrator, platform);
    }
    catch (const exception& e) {
        threwException = true;
    }
    ASSERT(threwException);
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 1.5, 0), Vec3(0, 0, 3));
    VerletIntegrator integrator(0.01);
    CustomAngleForce* angles = new CustomAngleForce("0.5*k*(theta-theta0)^2");
    angles->addPerAngleParameter("theta0");
    angles->addPerAngleParameter("k");
    vector<double> parameters(2);
    parameters[0] = M_PI/3;
    parameters[1] = 1.1;
    angles->addAngle(0, 1, 2, parameters);
    angles->setUsesPeriodicBoundaryConditions(true);
    system.addForce(angles);
    angles->setUsesPeriodicBoundaryConditions(true);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double torque = 1.1*M_PI/6;
    ASSERT_EQUAL_VEC(Vec3(2*torque, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -torque, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(0.5*1.1*(M_PI/6)*(M_PI/6), state.getPotentialEnergy(), TOL);
}

void testEnergyParameterDerivatives() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomAngleForce* angles = new CustomAngleForce("k*(theta-theta0)^2");
    angles->addGlobalParameter("theta0", 0.0);
    angles->addGlobalParameter("k", 0.0);
    angles->addEnergyParameterDerivative("theta0");
    angles->addEnergyParameterDerivative("k");
    vector<double> parameters;
    angles->addAngle(0, 1, 2, parameters);
    system.addForce(angles);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 1, 0);
    context.setPositions(positions);
    double theta = M_PI/4;
    for (int i = 0; i < 10; i++) {
        double theta0 = 0.1*i;
        double k = 10-i;
        context.setParameter("theta0", theta0);
        context.setParameter("k", k);
        State state = context.getState(State::ParameterDerivatives);
        map<string, double> derivs = state.getEnergyParameterDerivatives();
        double dEdtheta0 = -2*k*(theta-theta0);
        double dEdk = (theta-theta0)*(theta-theta0);
        ASSERT_EQUAL_TOL(dEdtheta0, derivs["theta0"], 1e-5);
        ASSERT_EQUAL_TOL(dEdk, derivs["k"], 1e-5);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testAngles();
        testIllegalVariable();
        testPeriodic();
        testEnergyParameterDerivatives();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


