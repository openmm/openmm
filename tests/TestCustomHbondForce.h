/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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
#include "openmm/CustomHbondForce.h"
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

void testHbond() {
    // Create a system using a CustomHbondForce.

    System customSystem;
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    customSystem.addParticle(1.0);
    CustomHbondForce* custom = new CustomHbondForce("0.5*kr*(distance(d1,a1)-r0)^2 + 0.5*ktheta*(angle(a1,d1,d2)-theta0)^2 + 0.5*kpsi*(angle(d1,a1,a2)-psi0)^2 + kchi*(1+cos(n*dihedral(a3,a2,a1,d1)-chi0))");
    custom->addPerDonorParameter("r0");
    custom->addPerDonorParameter("theta0");
    custom->addPerDonorParameter("psi0");
    custom->addPerAcceptorParameter("chi0");
    custom->addPerAcceptorParameter("n");
    custom->addGlobalParameter("kr", 0.4);
    custom->addGlobalParameter("ktheta", 0.5);
    custom->addGlobalParameter("kpsi", 0.6);
    custom->addGlobalParameter("kchi", 0.7);
    vector<double> parameters(3);
    parameters[0] = 1.5;
    parameters[1] = 1.7;
    parameters[2] = 1.9;
    custom->addDonor(1, 0, -1, parameters);
    parameters.resize(2);
    parameters[0] = 2.1;
    parameters[1] = 2;
    custom->addAcceptor(2, 3, 4, parameters);
    custom->setCutoffDistance(10.0);
    customSystem.addForce(custom);
    ASSERT(!custom->usesPeriodicBoundaryConditions());
    ASSERT(!customSystem.usesPeriodicBoundaryConditions());

    // Create an identical system using HarmonicBondForce, HarmonicAngleForce, and PeriodicTorsionForce.

    System standardSystem;
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    standardSystem.addParticle(1.0);
    HarmonicBondForce* bond = new HarmonicBondForce();
    bond->addBond(1, 2, 1.5, 0.4);
    standardSystem.addForce(bond);
    HarmonicAngleForce* angle = new HarmonicAngleForce();
    angle->addAngle(0, 1, 2, 1.7, 0.5);
    angle->addAngle(1, 2, 3, 1.9, 0.6);
    standardSystem.addForce(angle);
    PeriodicTorsionForce* torsion = new PeriodicTorsionForce();
    torsion->addTorsion(1, 2, 3, 4, 2, 2.1, 0.7);
    standardSystem.addForce(torsion);

    // Set the atoms in various positions, and verify that both systems give identical forces and energy.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<Vec3> positions(5);
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context c1(customSystem, integrator1, platform);
    Context c2(standardSystem, integrator2, platform);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(2.0*genrand_real2(sfmt), 2.0*genrand_real2(sfmt), 2.0*genrand_real2(sfmt));
        c1.setPositions(positions);
        c2.setPositions(positions);
        State s1 = c1.getState(State::Forces | State::Energy);
        State s2 = c2.getState(State::Forces | State::Energy);
        for (int i = 0; i < customSystem.getNumParticles(); i++)
            ASSERT_EQUAL_VEC(s2.getForces()[i], s1.getForces()[i], TOL);
        ASSERT_EQUAL_TOL(s2.getPotentialEnergy(), s1.getPotentialEnergy(), TOL);
    }
    
    // Try changing the parameters and make sure it's still correct.
    
    parameters.resize(3);
    parameters[0] = 1.4;
    parameters[1] = 1.7;
    parameters[2] = 1.9;
    custom->setDonorParameters(0, 1, 0, -1, parameters);
    parameters.resize(2);
    parameters[0] = 2.2;
    parameters[1] = 2;
    custom->setAcceptorParameters(0, 2, 3, 4, parameters);
    bond->setBondParameters(0, 1, 2, 1.4, 0.4);
    torsion->setTorsionParameters(0, 1, 2, 3, 4, 2, 2.2, 0.7);
    custom->updateParametersInContext(c1);
    bond->updateParametersInContext(c2);
    torsion->updateParametersInContext(c2);
    State s1 = c1.getState(State::Forces | State::Energy);
    State s2 = c2.getState(State::Forces | State::Energy);
    for (int i = 0; i < customSystem.getNumParticles(); i++)
        ASSERT_EQUAL_VEC(s2.getForces()[i], s1.getForces()[i], TOL);
    ASSERT_EQUAL_TOL(s2.getPotentialEnergy(), s1.getPotentialEnergy(), TOL);
}

void testExclusions() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomHbondForce* custom = new CustomHbondForce("(distance(d1,a1)-1)^2");
    custom->addDonor(0, 1, -1, vector<double>());
    custom->addDonor(1, 0, -1, vector<double>());
    custom->addAcceptor(2, 0, -1, vector<double>());
    custom->addExclusion(1, 0);
    system.addForce(custom);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(2, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(-2, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(1.0, state.getPotentialEnergy(), TOL);
}

void testCutoff() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomHbondForce* custom = new CustomHbondForce("(distance(d1,a1)-1)^2");
    custom->addDonor(0, 1, -1, vector<double>());
    custom->addDonor(1, 0, -1, vector<double>());
    custom->addAcceptor(2, 0, -1, vector<double>());
    custom->setNonbondedMethod(CustomHbondForce::CutoffNonPeriodic);
    custom->setCutoffDistance(2.5);
    system.addForce(custom);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 3, 0);
    positions[2] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(2, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(-2, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(1.0, state.getPotentialEnergy(), TOL);
}

void testCustomFunctions() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomHbondForce* custom = new CustomHbondForce("foo(distance(d1,a1))");
    custom->addDonor(1, 0, -1, vector<double>());
    custom->addDonor(2, 0, -1, vector<double>());
    custom->addAcceptor(0, 1, -1, vector<double>());
    vector<double> function(2);
    function[0] = 0;
    function[1] = 1;
    custom->addTabulatedFunction("foo", new Continuous1DFunction(function, 0, 10));
    system.addForce(custom);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0.1, 0.1, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -0.1, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(-0.1, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(0.1*2+0.1*2, state.getPotentialEnergy(), TOL);
}

void test2DFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomHbondForce* custom = new CustomHbondForce("tab(dtype, atype)");
    custom->addPerDonorParameter("dtype");
    custom->addPerAcceptorParameter("atype");
    custom->addDonor(1, 0, -1, {0.0});
    custom->addDonor(2, 0, -1, {2.0});
    custom->addAcceptor(0, 1, -1, {1.0});
    vector<double> function = {0.0, 1.0, 2.0, 5.0, 6.0, 7.0};
    custom->addTabulatedFunction("tab", new Discrete2DFunction(3, 2, function));
    system.addForce(custom);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(12.0, state.getPotentialEnergy(), TOL);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomHbondForce* force = new CustomHbondForce("1+none");
    force->addDonor(0, -1, -1);
    force->addAcceptor(1, -1, -1);
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

void testParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomHbondForce* custom = new CustomHbondForce("(2*d+a)*distance(d1,a1)");
    custom->addPerDonorParameter("d");
    custom->addPerAcceptorParameter("a");
    custom->addDonor(1, 0, -1, vector<double>({1.5}));
    custom->addDonor(2, 0, -1, vector<double>({1.8}));
    custom->addAcceptor(0, 1, -1, vector<double>({2.1}));
    system.addForce(custom);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3((2*1.8+2.1), (2*1.5+2.1), 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -(2*1.5+2.1), 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(-(2*1.8+2.1), 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(2*(2*1.8+2.1)+2*(2*1.5+2.1), state.getPotentialEnergy(), TOL);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHbond();
        testExclusions();
        testCutoff();
        testCustomFunctions();
        test2DFunction();
        testIllegalVariable();
        testParameters();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

