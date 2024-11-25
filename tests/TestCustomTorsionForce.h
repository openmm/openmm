/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
    ASSERT(!custom->usesPeriodicBoundaryConditions());
    ASSERT(!customSystem.usesPeriodicBoundaryConditions());

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
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context c1(customSystem, integrator1, platform);
    Context c2(harmonicSystem, integrator2, platform);
    for (int i = 0; i < 50; i++) {
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
    
    // Try changing the torsion parameters and make sure it's still correct.
    
    parameters[0] = 1.6;
    parameters[1] = 2;
    custom->setTorsionParameters(0, 0, 1, 2, 3, parameters);
    parameters[0] = 2.1;
    parameters[1] = 3;
    custom->setTorsionParameters(1, 1, 2, 3, 4, parameters);
    custom->updateParametersInContext(c1);
    periodic->setTorsionParameters(0, 0, 1, 2, 3, 2, 1.6, 0.5);
    periodic->setTorsionParameters(1, 1, 2, 3, 4, 3, 2.1, 0.5);
    periodic->updateParametersInContext(c2);
    {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        c1.setPositions(positions);
        c2.setPositions(positions);
        State s1 = c1.getState(State::Forces | State::Energy);
        State s2 = c2.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = s1.getForces();
        for (int i = 0; i < customSystem.getNumParticles(); i++)
            ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], TOL);
        ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), s2.getPotentialEnergy(), TOL);
    }
}

void testRange() {
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

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomTorsionForce* force = new CustomTorsionForce("theta+none");
    force->addTorsion(0, 1, 2, 3);
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
    system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    VerletIntegrator integrator(0.01);
    CustomTorsionForce* torsions = new CustomTorsionForce("k*(1+cos(n*theta-theta0))");
    torsions->addPerTorsionParameter("theta0");
    torsions->addPerTorsionParameter("n");
    torsions->addGlobalParameter("k", 1.1);
    vector<double> parameters(2);
    parameters[0] = M_PI/3;
    parameters[1] = 2;
    torsions->addTorsion(0, 1, 2, 3, parameters);
    torsions->setUsesPeriodicBoundaryConditions(true);
    system.addForce(torsions);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 0, 2);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double torque = -2*1.1*std::sin(2*M_PI/3);
    ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -torque, 0), forces[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0], forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1], forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]), Vec3(0, 0, 0), TOL);
    ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*M_PI/3)), state.getPotentialEnergy(), TOL);
}

void testEnergyParameterDerivatives() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomTorsionForce* torsions = new CustomTorsionForce("k*(theta-theta0)^2");
    torsions->addGlobalParameter("theta0", 0.0);
    torsions->addGlobalParameter("k", 0.0);
    torsions->addEnergyParameterDerivative("theta0");
    torsions->addEnergyParameterDerivative("k");
    vector<double> parameters;
    torsions->addTorsion(0, 1, 2, 3, parameters);
    system.addForce(torsions);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    positions[3] = Vec3(1, 1, 1);
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

void testParallelComputation() {
    System system;
    const int numParticles = 200;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    CustomTorsionForce* force = new CustomTorsionForce("k*sin(theta-theta0)");
    force->addPerTorsionParameter("k");
    force->addPerTorsionParameter("theta0");
    for (int i = 3; i < numParticles; i++)
        force->addTorsion(i-3, i-2, i-1, i, {1.0, 1.1});
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(i, i%2, i%3);
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    VerletIntegrator integrator2(0.01);
    string deviceIndex = platform.getPropertyValue(context1, "DeviceIndex");
    map<string, string> props;
    props["DeviceIndex"] = deviceIndex+","+deviceIndex;
    Context context2(system, integrator2, platform, props);
    context2.setPositions(positions);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-5);

    // Try updating some parameters and see if they still match.

    vector<double> params;
    for (int i = 95; i < 102; i++) {
        int p1, p2, p3, p4;
        force->getTorsionParameters(i, p1, p2, p3, p4, params);
        force->setTorsionParameters(i, p1, p2, p3, p4, {2.0, 1.2});
    }
    force->updateParametersInContext(context1);
    force->updateParametersInContext(context2);
    State state3 = context1.getState(State::Energy);
    State state4 = context2.getState(State::Energy);
    ASSERT_EQUAL_TOL(state3.getPotentialEnergy(), state4.getPotentialEnergy(), 1e-5);
    ASSERT(fabs(state1.getPotentialEnergy()-state3.getPotentialEnergy()) > 0.1);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testTorsions();
        testRange();
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



