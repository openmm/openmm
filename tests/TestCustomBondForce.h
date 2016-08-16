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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/CustomBondForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testBonds() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomBondForce* forceField = new CustomBondForce("scale*k*(r-r0)^2");
    forceField->addPerBondParameter("r0");
    forceField->addPerBondParameter("k");
    forceField->addGlobalParameter("scale", 0.5);
    vector<double> parameters(2);
    parameters[0] = 1.5;
    parameters[1] = 0.8;
    forceField->addBond(0, 1, parameters);
    parameters[0] = 1.2;
    parameters[1] = 0.7;
    forceField->addBond(1, 2, parameters);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, -0.8*0.5, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0.7*0.2, 0, 0), forces[2], TOL);
        ASSERT_EQUAL_VEC(Vec3(-forces[0][0]-forces[2][0], -forces[0][1]-forces[2][1], -forces[0][2]-forces[2][2]), forces[1], TOL);
        ASSERT_EQUAL_TOL(0.5*0.8*0.5*0.5 + 0.5*0.7*0.2*0.2, state.getPotentialEnergy(), TOL);
    }
    
    // Try changing the bond parameters and make sure it's still correct.
    
    parameters[0] = 1.6;
    parameters[1] = 0.9;
    forceField->setBondParameters(0, 0, 1, parameters);
    parameters[0] = 1.3;
    parameters[1] = 0.8;
    forceField->setBondParameters(1, 1, 2, parameters);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    {
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, -0.9*0.4, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0.8*0.3, 0, 0), forces[2], TOL);
        ASSERT_EQUAL_VEC(Vec3(-forces[0][0]-forces[2][0], -forces[0][1]-forces[2][1], -forces[0][2]-forces[2][2]), forces[1], TOL);
        ASSERT_EQUAL_TOL(0.5*0.9*0.4*0.4 + 0.5*0.8*0.3*0.3, state.getPotentialEnergy(), TOL);
    }
}

void testManyParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomBondForce* forceField = new CustomBondForce("(a+b+c+d+e+f+g+h+i)*r");
    forceField->addPerBondParameter("a");
    forceField->addPerBondParameter("b");
    forceField->addPerBondParameter("c");
    forceField->addPerBondParameter("d");
    forceField->addPerBondParameter("e");
    forceField->addPerBondParameter("f");
    forceField->addPerBondParameter("g");
    forceField->addPerBondParameter("h");
    forceField->addPerBondParameter("i");
    vector<double> parameters(forceField->getNumPerBondParameters());
    for (int i = 0; i < parameters.size(); i++)
        parameters[i] = i;
    forceField->addBond(0, 1, parameters);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2.5, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double f = 1+2+3+4+5+6+7+8;
    ASSERT_EQUAL_VEC(Vec3(0, f, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -f, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(f*2.5, state.getPotentialEnergy(), TOL);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomBondForce* force = new CustomBondForce("r+none");
    force->addBond(0, 1);
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
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    VerletIntegrator integrator(0.01);
    CustomBondForce* forceField = new CustomBondForce("scale*k*(r-r0)^2");
    forceField->addPerBondParameter("r0");
    forceField->addPerBondParameter("k");
    forceField->addGlobalParameter("scale", 0.5);
    vector<double> parameters(2);
    parameters[0] = 1.9;
    parameters[1] = 0.8;
    forceField->addBond(0, 1, parameters);
    forceField->setUsesPeriodicBoundaryConditions(true);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, -0.8*0.9, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0.8*0.9, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(0.5*0.8*0.9*0.9, state.getPotentialEnergy(), TOL);
}

void testEnergyParameterDerivatives() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomBondForce* bonds = new CustomBondForce("k*(r-r0)^2");
    bonds->addGlobalParameter("r0", 0.0);
    bonds->addGlobalParameter("k", 0.0);
    bonds->addEnergyParameterDerivative("k");
    bonds->addEnergyParameterDerivative("r0");
    vector<double> parameters;
    bonds->addBond(0, 1, parameters);
    bonds->addBond(1, 2, parameters);
    system.addForce(bonds);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    for (int i = 0; i < 10; i++) {
        double r0 = 0.1*i;
        double k = 10-i;
        context.setParameter("r0", r0);
        context.setParameter("k", k);
        State state = context.getState(State::ParameterDerivatives);
        map<string, double> derivs = state.getEnergyParameterDerivatives();
        double dEdr0 = -2*k*((2-r0)+(1-r0));
        double dEdk = (2-r0)*(2-r0) + (1-r0)*(1-r0);
        ASSERT_EQUAL_TOL(dEdr0, derivs["r0"], 1e-5);
        ASSERT_EQUAL_TOL(dEdk, derivs["k"], 1e-5);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testBonds();
        testManyParameters();
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

