/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2016 Stanford University and the Authors.      *
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
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/System.h"
#include "openmm/TabulatedFunction.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testHarmonicBond() {
    System system;
    system.addParticle(1.0);
    system.addParticle(2.0);
    system.addParticle(3.0);
    system.addParticle(4.0);
    system.addParticle(5.0);
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "k*distance(g1,g2)^2");
    force->addPerBondParameter("k");
    vector<int> particles1;
    particles1.push_back(0);
    particles1.push_back(1);
    vector<int> particles2;
    particles2.push_back(2);
    particles2.push_back(3);
    particles2.push_back(4);
    force->addGroup(particles1);
    force->addGroup(particles2);
    vector<int> groups;
    groups.push_back(0);
    groups.push_back(1);
    vector<double> parameters;
    parameters.push_back(1.0);
    force->addBond(groups, parameters);
    system.addForce(force);
    ASSERT(!system.usesPeriodicBoundaryConditions());

    // The center of mass of group 0 is (1.5, 0, 0).

    vector<Vec3> positions(5);
    positions[0] = Vec3(2.5, 0, 0);
    positions[1] = Vec3(1, 0, 0);

    // The center of mass of group 1 is (-1, 0, 0).

    positions[2] = Vec3(-6, 0, 0);
    positions[3] = Vec3(-1, 0, 0);
    positions[4] = Vec3(2, 0, 0);

    // Check the forces and energy.

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(2.5*2.5, state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(Vec3(-2*2.5*(1.0/3.0), 0, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-2*2.5*(2.0/3.0), 0, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*2.5*(3.0/12.0), 0, 0), state.getForces()[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*2.5*(4.0/12.0), 0, 0), state.getForces()[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*2.5*(5.0/12.0), 0, 0), state.getForces()[4], TOL);

    // Update the per-bond parameter and see if the results change.

    parameters[0] = 2.0;
    force->setBondParameters(0, groups, parameters);
    force->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(2*2.5*2.5, state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(Vec3(-4*2.5*(1.0/3.0), 0, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-4*2.5*(2.0/3.0), 0, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(4*2.5*(3.0/12.0), 0, 0), state.getForces()[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(4*2.5*(4.0/12.0), 0, 0), state.getForces()[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(4*2.5*(5.0/12.0), 0, 0), state.getForces()[4], TOL);
    
    // All the particles should be treated as a single molecule.
    
    vector<std::vector<int> > molecules = context.getMolecules();
    ASSERT_EQUAL(1, molecules.size());
    ASSERT_EQUAL(5, molecules[0].size());
}

void testComplexFunction() {
    int numParticles = 5;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(2.0);
    vector<double> table(20);
    for (int i = 0; i < 20; i++)
        table[i] = sin(0.11*i);

    // When every group contains only one particle, a CustomCentroidBondForce is identical to a
    // CustomCompoundBondForce.  Use that to test a complicated energy function with lots of terms.

    CustomCompoundBondForce* compound = new CustomCompoundBondForce(4, "x1+y2+z4+fn(distance(p1,p2))*angle(p3,p2,p4)+scale*dihedral(p2,p1,p4,p3)");
    CustomCentroidBondForce* centroid = new CustomCentroidBondForce(4, "x1+y2+z4+fn(distance(g1,g2))*angle(g3,g2,g4)+scale*dihedral(g2,g1,g4,g3)");
    compound->addGlobalParameter("scale", 0.5);
    centroid->addGlobalParameter("scale", 0.5);
    compound->addTabulatedFunction("fn", new Continuous1DFunction(table, -1, 10));
    centroid->addTabulatedFunction("fn", new Continuous1DFunction(table, -1, 10));

    // Add two bonds to the CustomCompoundBondForce.

    vector<int> particles(4);
    vector<double> parameters;
    particles[0] = 0;
    particles[1] = 1;
    particles[2] = 2;
    particles[3] = 3;
    compound->addBond(particles, parameters);
    particles[0] = 2;
    particles[1] = 4;
    particles[2] = 3;
    particles[3] = 1;
    compound->addBond(particles, parameters);

    // Add identical bonds to the CustomCentroidBondForce.  As a stronger test, make sure that
    // group number is different from particle number.

    vector<int> groupMembers(1);
    groupMembers[0] = 3;
    centroid->addGroup(groupMembers);
    groupMembers[0] = 0;
    centroid->addGroup(groupMembers);
    groupMembers[0] = 1;
    centroid->addGroup(groupMembers);
    groupMembers[0] = 2;
    centroid->addGroup(groupMembers);
    groupMembers[0] = 4;
    centroid->addGroup(groupMembers);
    vector<int> groups(4);
    groups[0] = 1;
    groups[1] = 2;
    groups[2] = 3;
    groups[3] = 0;
    centroid->addBond(groups, parameters);
    groups[0] = 3;
    groups[1] = 4;
    groups[2] = 0;
    groups[3] = 2;
    centroid->addBond(groups, parameters);

    // Add both forces as different force groups, and create a context.

    centroid->setForceGroup(1);
    system.addForce(compound);
    system.addForce(centroid);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);

    // Evaluate the force and energy for various positions and see if they match.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < numParticles; j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        context.setPositions(positions);
        State state1 = context.getState(State::Forces | State::Energy, false, 1<<0);
        State state2 = context.getState(State::Forces | State::Energy, false, 1<<1);
        ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), TOL);
        for (int i = 0; i < numParticles; i++)
            ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], TOL);
    }
}

void testCustomWeights() {
    System system;
    system.addParticle(1.0);
    system.addParticle(2.0);
    system.addParticle(3.0);
    system.addParticle(4.0);
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "distance(g1,g2)^2");
    vector<int> particles(2);
    vector<double> weights(2);
    particles[0] = 0;
    particles[1] = 1;
    weights[0] = 0.5;
    weights[1] = 1.5;
    force->addGroup(particles, weights);
    particles[0] = 2;
    particles[1] = 3;
    weights[0] = 2.0;
    weights[1] = 1.0;
    force->addGroup(particles, weights);
    vector<int> groups;
    groups.push_back(0);
    groups.push_back(1);
    vector<double> parameters;
    force->addBond(groups, parameters);
    system.addForce(force);

    // The center of mass of group 0 is (0, 1, 0).

    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 4, 0);
    positions[1] = Vec3(0, 0, 0);

    // The center of mass of group 1 is (0, 10, 0).

    positions[2] = Vec3(0, 9, 0);
    positions[3] = Vec3(0, 12, 0);

    // Check the forces and energy.

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(9.0*9.0, state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 2*9*(0.5/2.0), 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 2*9*(1.5/2.0), 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -2*9*(2.0/3.0), 0), state.getForces()[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -2*9*(1.0/3.0), 0), state.getForces()[3], TOL);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "1+none");
    vector<int> particles;
    particles.push_back(0);
    force->addGroup(particles);
    force->addGroup(particles);
    vector<int> groups;
    groups.push_back(0);
    groups.push_back(1);
    force->addBond(groups);
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
    system.addParticle(2.0);
    system.addParticle(3.0);
    system.addParticle(4.0);
    system.addParticle(5.0);
    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "k*distance(g1,g2)^2");
    force->addPerBondParameter("k");
    vector<int> particles1;
    particles1.push_back(0);
    particles1.push_back(1);
    vector<int> particles2;
    particles2.push_back(2);
    particles2.push_back(3);
    particles2.push_back(4);
    force->addGroup(particles1);
    force->addGroup(particles2);
    vector<int> groups;
    groups.push_back(0);
    groups.push_back(1);
    vector<double> parameters;
    parameters.push_back(1.0);
    force->addBond(groups, parameters);
    force->setUsesPeriodicBoundaryConditions(true);
    system.addForce(force);

    // The center of mass of group 0 is (1.5, 0, 0).

    vector<Vec3> positions(5);
    positions[0] = Vec3(2.5, 0, 0);
    positions[1] = Vec3(1, 0, 0);

    // The center of mass of group 1 is (-1, 0, 0).

    positions[2] = Vec3(-6, 0, 0);
    positions[3] = Vec3(-1, 0, 0);
    positions[4] = Vec3(2, 0, 0);

    // Check the forces and energy.

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(0.5*0.5, state.getPotentialEnergy(), TOL);
    ASSERT_EQUAL_VEC(Vec3(-2*0.5*(1.0/3.0), 0, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-2*0.5*(2.0/3.0), 0, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*0.5*(3.0/12.0), 0, 0), state.getForces()[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*0.5*(4.0/12.0), 0, 0), state.getForces()[3], TOL);
    ASSERT_EQUAL_VEC(Vec3(2*0.5*(5.0/12.0), 0, 0), state.getForces()[4], TOL);
}

void testEnergyParameterDerivatives() {
    System system;
    system.addParticle(1.0);
    system.addParticle(2.0);
    system.addParticle(3.0);
    system.addParticle(4.0);
    system.addParticle(5.0);
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "k*(distance(g1,g2)-r0)^2");
    force->addGlobalParameter("r0", 0.0);
    force->addGlobalParameter("k", 0.0);
    force->addEnergyParameterDerivative("r0");
    force->addEnergyParameterDerivative("k");
    vector<int> particles1;
    particles1.push_back(0);
    particles1.push_back(1);
    vector<int> particles2;
    particles2.push_back(2);
    particles2.push_back(3);
    particles2.push_back(4);
    force->addGroup(particles1);
    force->addGroup(particles2);
    vector<int> groups;
    groups.push_back(0);
    groups.push_back(1);
    vector<double> parameters;
    force->addBond(groups, parameters);
    system.addForce(force);

    // The center of mass of group 0 is (1.5, 0, 0).

    vector<Vec3> positions(5);
    positions[0] = Vec3(2.5, 0, 0);
    positions[1] = Vec3(1, 0, 0);

    // The center of mass of group 1 is (-1, 0, 0).

    positions[2] = Vec3(-6, 0, 0);
    positions[3] = Vec3(-1, 0, 0);
    positions[4] = Vec3(2, 0, 0);
    
    // Check the derivatives.

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    for (int i = 0; i < 10; i++) {
        double r0 = 0.1*i;
        double k = 10-i;
        context.setParameter("r0", r0);
        context.setParameter("k", k);
        State state = context.getState(State::ParameterDerivatives);
        map<string, double> derivs = state.getEnergyParameterDerivatives();
        double dEdr0 = -2*k*(2.5-r0);
        double dEdk = (2.5-r0)*(2.5-r0);
        ASSERT_EQUAL_TOL(dEdr0, derivs["r0"], 1e-5);
        ASSERT_EQUAL_TOL(dEdk, derivs["k"], 1e-5);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHarmonicBond();
        testComplexFunction();
        testCustomWeights();
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
