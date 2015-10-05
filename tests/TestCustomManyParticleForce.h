/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2015 Stanford University and the Authors.      *
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
#include "openmm/CustomCompoundBondForce.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/System.h"
#include "openmm/TabulatedFunction.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

Vec3 computeDelta(const Vec3& pos1, const Vec3& pos2, bool periodic, const Vec3* periodicBoxVectors) {
    Vec3 diff = pos1-pos2;
    if (periodic) {
        diff -= periodicBoxVectors[2]*floor(diff[2]/periodicBoxVectors[2][2]+0.5);
        diff -= periodicBoxVectors[1]*floor(diff[1]/periodicBoxVectors[1][1]+0.5);
        diff -= periodicBoxVectors[0]*floor(diff[0]/periodicBoxVectors[0][0]+0.5);
    }
    return diff;
}

void validateAxilrodTeller(CustomManyParticleForce* force, const vector<Vec3>& positions, const vector<const int*>& expectedSets, double boxSize, bool triclinic) {
    // Create a System and Context.
    
    int numParticles = force->getNumParticles();
    CustomManyParticleForce::NonbondedMethod nonbondedMethod = force->getNonbondedMethod();
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    Vec3 boxVectors[3];
    if (triclinic) {
        boxVectors[0] = Vec3(boxSize, 0, 0);
        boxVectors[1] = Vec3(0.2*boxSize, boxSize, 0);
        boxVectors[2] = Vec3(-0.3*boxSize, -0.1*boxSize, boxSize);
    }
    else {
        boxVectors[0] = Vec3(boxSize, 0, 0);
        boxVectors[1] = Vec3(0, boxSize, 0);
        boxVectors[2] = Vec3(0, 0, boxSize);
    }
    system.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    system.addForce(force);
    if (force->getNonbondedMethod() == CustomManyParticleForce::CutoffPeriodic) {
        ASSERT(force->usesPeriodicBoundaryConditions());
        ASSERT(system.usesPeriodicBoundaryConditions());
    }
    else {
        ASSERT(!force->usesPeriodicBoundaryConditions());
        ASSERT(!system.usesPeriodicBoundaryConditions());
    }
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    double c = context.getParameter("C");
    
    // See if the energy matches the expected value.
    
    double expectedEnergy = 0;
    bool periodic = (nonbondedMethod == CustomManyParticleForce::CutoffPeriodic);
    for (int i = 0; i < (int) expectedSets.size(); i++) {
        int p1 = expectedSets[i][0];
        int p2 = expectedSets[i][1];
        int p3 = expectedSets[i][2];
        Vec3 d12 = computeDelta(positions[p2], positions[p1], periodic, boxVectors);
        Vec3 d13 = computeDelta(positions[p3], positions[p1], periodic, boxVectors);
        Vec3 d23 = computeDelta(positions[p3], positions[p2], periodic, boxVectors);
        double r12 = sqrt(d12.dot(d12));
        double r13 = sqrt(d13.dot(d13));
        double r23 = sqrt(d23.dot(d23));
        double ctheta1 = d12.dot(d13)/(r12*r13);
        double ctheta2 = -d12.dot(d23)/(r12*r23);
        double ctheta3 = d13.dot(d23)/(r13*r23);
        double rprod = r12*r13*r23;
        expectedEnergy += c*(1+3*ctheta1*ctheta2*ctheta3)/(rprod*rprod*rprod);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state1.getPotentialEnergy(), 1e-5);

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    const vector<Vec3>& forces = state1.getForces();
    double norm = 0.0;
    for (int i = 0; i < (int) forces.size(); ++i)
        norm += forces[i].dot(forces[i]);
    norm = std::sqrt(norm);
    const double stepSize = 1e-3;
    double step = 0.5*stepSize/norm;
    vector<Vec3> positions2(numParticles), positions3(numParticles);
    for (int i = 0; i < (int) positions.size(); ++i) {
        Vec3 p = positions[i];
        Vec3 f = forces[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-4);
}

void validateStillingerWeber(CustomManyParticleForce* force, const vector<Vec3>& positions, const vector<const int*>& expectedSets, double boxSize) {
    // Create a System and Context.
    
    int numParticles = force->getNumParticles();
    CustomManyParticleForce::NonbondedMethod nonbondedMethod = force->getNonbondedMethod();
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    double L = context.getParameter("L");
    double eps = context.getParameter("eps");
    double a = context.getParameter("a");
    double gamma = context.getParameter("gamma");
    double sigma = context.getParameter("sigma");
    
    // See if the energy matches the expected value.
    
    double expectedEnergy = 0;
    for (int i = 0; i < (int) expectedSets.size(); i++) {
        int p1 = expectedSets[i][0];
        int p2 = expectedSets[i][1];
        int p3 = expectedSets[i][2];
        Vec3 d12 = positions[p2]-positions[p1];
        Vec3 d13 = positions[p3]-positions[p1];
        Vec3 d23 = positions[p3]-positions[p2];
        if (nonbondedMethod == CustomManyParticleForce::CutoffPeriodic) {
            for (int j = 0; j < 3; j++) {
                d12[j] -= floor(d12[j]/boxSize+0.5f)*boxSize;
                d13[j] -= floor(d13[j]/boxSize+0.5f)*boxSize;
                d23[j] -= floor(d23[j]/boxSize+0.5f)*boxSize;
            }
        }
        double r12 = sqrt(d12.dot(d12));
        double r13 = sqrt(d13.dot(d13));
        double r23 = sqrt(d23.dot(d23));
        double ctheta1 = d12.dot(d13)/(r12*r13);
        double ctheta2 = -d12.dot(d23)/(r12*r23);
        double ctheta3 = d13.dot(d23)/(r13*r23);
        expectedEnergy += L*eps*(ctheta1+1.0/3.0)*(ctheta1+1.0/3.0)*exp(sigma*gamma/(r12-a*sigma))*exp(sigma*gamma/(r13-a*sigma));
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state1.getPotentialEnergy(), 1e-5);

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    const vector<Vec3>& forces = state1.getForces();
    double norm = 0.0;
    for (int i = 0; i < (int) forces.size(); ++i)
        norm += forces[i].dot(forces[i]);
    norm = std::sqrt(norm);
    const double stepSize = 1e-3;
    double step = 0.5*stepSize/norm;
    vector<Vec3> positions2(numParticles), positions3(numParticles);
    for (int i = 0; i < (int) positions.size(); ++i) {
        Vec3 p = positions[i];
        Vec3 f = forces[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-4);
}

void testNoCutoff() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    int sets[4][3] = {{0,1,2}, {1,2,3}, {2,3,0}, {3,0,1}};
    vector<const int*> expectedSets(&sets[0], &sets[4]);
    validateAxilrodTeller(force, positions, expectedSets, 2.0, false);
}

void testCutoff() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffNonPeriodic);
    force->setCutoffDistance(1.55);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    positions.push_back(Vec3(0.2, 0.5, -0.1));
    int sets[7][3] = {{0,1,2}, {0,1,3}, {0,1,4}, {0,2,4}, {0,3,4}, {1,2,4}, {1,3,4}};
    vector<const int*> expectedSets(&sets[0], &sets[7]);
    validateAxilrodTeller(force, positions, expectedSets, 2.0, false);
}

void testPeriodic() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffPeriodic);
    force->setCutoffDistance(1.05);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    positions.push_back(Vec3(0.2, 0.5, -0.1));
    double boxSize = 2.1;
    int sets[5][3] = {{0,1,3}, {0,1,4}, {0,2,4}, {0,3,4}, {1,3,4}};
    vector<const int*> expectedSets(&sets[0], &sets[5]);
    validateAxilrodTeller(force, positions, expectedSets, boxSize, false);
}

void testTriclinic() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffPeriodic);
    force->setCutoffDistance(1.05);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    positions.push_back(Vec3(0.2, 0.5, -0.1));
    double boxSize = 2.1;
    int sets[4][3] = {{0,1,3}, {0,1,4}, {0,3,4}, {1,3,4}};
    vector<const int*> expectedSets(&sets[0], &sets[4]);
    validateAxilrodTeller(force, positions, expectedSets, boxSize, true);
}

void testExclusions() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    positions.push_back(Vec3(0.2, 0.5, -0.1));
    force->addExclusion(0, 2);
    force->addExclusion(0, 3);
    int sets[5][3] = {{0,1,4}, {1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}};
    vector<const int*> expectedSets(&sets[0], &sets[5]);
    validateAxilrodTeller(force, positions, expectedSets, 2.0, false);
}

void testAllTerms() {
    int numParticles = 4;
    
    // Create a system with a CustomManyParticleForce.
    
    System system1;
    CustomManyParticleForce* force1 = new CustomManyParticleForce(4,
        "distance(p1,p2)+angle(p1,p4,p3)+dihedral(p1,p3,p2,p4)+x1+y4+z3");
    system1.addForce(force1);
    vector<double> params;
    for (int i = 0; i < numParticles; i++) {
        system1.addParticle(1.0);
        force1->addParticle(params, i);
    }
    set<int> filter;
    filter.insert(0);
    force1->setTypeFilter(0, filter);
    filter.clear();
    filter.insert(1);
    force1->setTypeFilter(1, filter);
    filter.clear();
    filter.insert(3);
    force1->setTypeFilter(2, filter);
    filter.clear();
    filter.insert(2);
    force1->setTypeFilter(3, filter);
    
    // Create a system that use a CustomCompoundBondForce to compute exactly the same interactions.
    
    System system2;
    CustomCompoundBondForce* force2 = new CustomCompoundBondForce(4,
        "distance(p1,p2)+angle(p1,p3,p4)+dihedral(p1,p4,p2,p3)+x1+y3+z4");
    system2.addForce(force2);
    vector<int> particles;
    particles.push_back(0);
    particles.push_back(1);
    particles.push_back(2);
    particles.push_back(3);
    force2->addBond(particles, params);
    for (int i = 0; i < numParticles; i++)
        system2.addParticle(1.0);
    
    // Create contexts for both of them.

    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++)
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
    VerletIntegrator integrator1(0.001);
    VerletIntegrator integrator2(0.001);
    Context context1(system1, integrator1, platform);
    Context context2(system2, integrator2, platform);
    context1.setPositions(positions);
    context2.setPositions(positions);
    
    // See if they produce identical forces and energies.
    
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state1.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state2.getForces()[i], state1.getForces()[i], 1e-4);
}

void testParameters() {
    // Create a system.
    
    int numParticles = 5;
    System system;
    CustomManyParticleForce* force = new CustomManyParticleForce(3, "C*scale1*scale2*scale3*(distance(p1,p2)+distance(p2,p3)+distance(p1,p3))");
    force->addGlobalParameter("C", 2.0);
    force->addPerParticleParameter("scale");
    vector<double> params(1);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        params[0] = i+1;
        force->addParticle(params);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
        system.addParticle(1.0);
    }
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the energy is correct.

    State state = context.getState(State::Energy);
    double expectedEnergy = 0;
    for (int i = 0; i < numParticles; i++)
        for (int j = i+1; j < numParticles; j++)
            for (int k = j+1; k < numParticles; k++) {
                Vec3 d12 = positions[j]-positions[i];
                Vec3 d13 = positions[k]-positions[i];
                Vec3 d23 = positions[k]-positions[j];
                double r12 = sqrt(d12.dot(d12));
                double r13 = sqrt(d13.dot(d13));
                double r23 = sqrt(d23.dot(d23));
                expectedEnergy += 2.0*(i+1)*(j+1)*(k+1)*(r12+r13+r23);
            }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);
    
    // Modify the parameters.
    
    context.setParameter("C", 3.5);
    for (int i = 0; i < numParticles; i++) {
        params[0] = 0.5*i-0.1;
        force->setParticleParameters(i, params, 0);
    }
    force->updateParametersInContext(context);
    
    // See if the energy is still correct.
    
    state = context.getState(State::Energy);
    expectedEnergy = 0;
    for (int i = 0; i < numParticles; i++)
        for (int j = i+1; j < numParticles; j++)
            for (int k = j+1; k < numParticles; k++) {
                Vec3 d12 = positions[j]-positions[i];
                Vec3 d13 = positions[k]-positions[i];
                Vec3 d23 = positions[k]-positions[j];
                double r12 = sqrt(d12.dot(d12));
                double r13 = sqrt(d13.dot(d13));
                double r23 = sqrt(d23.dot(d23));
                expectedEnergy += 3.5*(0.5*i-0.1)*(0.5*j-0.1)*(0.5*k-0.1)*(r12+r13+r23);
            }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);
}

void testTabulatedFunctions() {
    int numParticles = 5;
    
    // Create two tabulated functions.
    
    vector<double> values;
    values.push_back(0.0);
    values.push_back(50.0);
    Continuous1DFunction* f1 = new Continuous1DFunction(values, 0, 100);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<double> c(numParticles);
    for (int i = 0; i < numParticles; i++)
        c[i] = genrand_real2(sfmt);
    values.resize(numParticles*numParticles*numParticles);
    for (int i = 0; i < numParticles; i++)
        for (int j = 0; j < numParticles; j++)
            for (int k = 0; k < numParticles; k++)
                values[i+numParticles*j+numParticles*numParticles*k] = c[i]+c[j]+c[k];
    Discrete3DFunction* f2 = new Discrete3DFunction(numParticles, numParticles, numParticles, values);
    
    // Create a system.
    
    System system;
    CustomManyParticleForce* force = new CustomManyParticleForce(3, "f1(distance(p1,p2)+distance(p2,p3)+distance(p1,p3))*f2(atom1, atom2, atom3)");
    force->addPerParticleParameter("atom");
    force->addTabulatedFunction("f1", f1);
    force->addTabulatedFunction("f2", f2);
    vector<double> params(1);
    vector<Vec3> positions;
    for (int i = 0; i < numParticles; i++) {
        params[0] = i;
        force->addParticle(params);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
        system.addParticle(1.0);
    }
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the energy is correct.

    State state = context.getState(State::Energy);
    double expectedEnergy = 0;
    for (int i = 0; i < numParticles; i++)
        for (int j = i+1; j < numParticles; j++)
            for (int k = j+1; k < numParticles; k++) {
                Vec3 d12 = positions[j]-positions[i];
                Vec3 d13 = positions[k]-positions[i];
                Vec3 d23 = positions[k]-positions[j];
                double r12 = sqrt(d12.dot(d12));
                double r13 = sqrt(d13.dot(d13));
                double r23 = sqrt(d23.dot(d23));
                expectedEnergy += 0.5*(r12+r13+r23)*(c[i]+c[j]+c[k]);
            }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);
}

void testTypeFilters() {
    // Create a system.
    
    System system;
    for (int i = 0; i < 5; i++)
        system.addParticle(1.0);
    CustomManyParticleForce* force = new CustomManyParticleForce(3, "c1*(distance(p1,p2)+distance(p1,p3))");
    force->addPerParticleParameter("c");
    double c[] = {1.0, 2.0, 1.3, 1.5, -2.1};
    int type[] = {0, 1, 0, 1, 5};
    vector<double> params(1);
    for (int i = 0; i < 5; i++) {
        params[0] = c[i];
        force->addParticle(params, type[i]);
    }
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(1, 0, 0));
    positions.push_back(Vec3(0, 1.1, 0.3));
    positions.push_back(Vec3(0.4, 0, -0.8));
    positions.push_back(Vec3(0.2, 0.5, -0.1));
    set<int> f1, f2;
    f1.insert(0);
    f2.insert(1);
    f2.insert(5);
    force->setTypeFilter(0, f1);
    force->setTypeFilter(1, f2);
    force->setTypeFilter(2, f2);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the energy is correct.

    State state = context.getState(State::Energy);
    double expectedEnergy = 0;
    int sets[6][3] = {{0,1,3}, {0,1,4}, {0,3,4}, {2,1,3}, {2,1,4}, {2,3,4}};
    for (int i = 0; i < 6; i++) {
        int p1 = sets[i][0];
        int p2 = sets[i][1];
        int p3 = sets[i][2];
            Vec3 d12 = positions[p2]-positions[p1];
            Vec3 d13 = positions[p3]-positions[p1];
            double r12 = sqrt(d12.dot(d12));
            double r13 = sqrt(d13.dot(d13));
            expectedEnergy += c[p1]*(r12+r13);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);
}

void testLargeSystem() {
    int gridSize = 8;
    int numParticles = gridSize*gridSize*gridSize;
    double boxSize = 3.0;
    double spacing = boxSize/gridSize;
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "C*(1+3*cos(theta1)*cos(theta2)*cos(theta3))/(r12*r13*r23)^3;"
        "theta1=angle(p1,p2,p3); theta2=angle(p2,p3,p1); theta3=angle(p3,p1,p2);"
        "r12=distance(p1,p2); r13=distance(p1,p3); r23=distance(p2,p3)");
    force->addGlobalParameter("C", 1.5);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffPeriodic);
    force->setCutoffDistance(0.6);
    vector<double> params;
    vector<Vec3> positions;
    System system;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                force->addParticle(params);
                positions.push_back(Vec3((i+0.4*genrand_real2(sfmt))*spacing, (j+0.4*genrand_real2(sfmt))*spacing, (k+0.4*genrand_real2(sfmt))*spacing));
                system.addParticle(1.0);
            }
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(force);
    VerletIntegrator integrator1(0.001);
    VerletIntegrator integrator2(0.001);
    Context context1(system, integrator1, Platform::getPlatformByName("Reference"));
    Context context2(system, integrator2, platform);
    context1.setPositions(positions);
    context2.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
}

void testCentralParticleModeNoCutoff() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "L*eps*(cos(theta1)+1/3)^2*exp(sigma*gamma/(r12-a*sigma))*exp(sigma*gamma/(r13-a*sigma));"
        "r12 = distance(p1,p2); r13 = distance(p1,p3); theta1 = angle(p3,p1,p2)");
    force->setPermutationMode(CustomManyParticleForce::UniqueCentralParticle);
    force->addGlobalParameter("L", 23.13);
    force->addGlobalParameter("eps", 25.894776);
    force->addGlobalParameter("a", 1.8);
    force->addGlobalParameter("sigma", 0.23925);
    force->addGlobalParameter("gamma", 1.2);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(0.1, 0, 0));
    positions.push_back(Vec3(0, 0.11, 0.03));
    positions.push_back(Vec3(0.04, 0, -0.08));
    int sets[12][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,0,2}, {1,0,3}, {1, 2, 3}, {2,0,1}, {2,0,3}, {2, 1, 3}, {3,0,1}, {3,0,2}, {3,1,2}};
    vector<const int*> expectedSets(&sets[0], &sets[12]);
    validateStillingerWeber(force, positions, expectedSets, 2.0);
}

void testCentralParticleModeCutoff() {
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "L*eps*(cos(theta1)+1/3)^2*exp(sigma*gamma/(r12-a*sigma))*exp(sigma*gamma/(r13-a*sigma));"
        "r12 = distance(p1,p2); r13 = distance(p1,p3); theta1 = angle(p3,p1,p2)");
    force->setPermutationMode(CustomManyParticleForce::UniqueCentralParticle);
    force->addGlobalParameter("L", 23.13);
    force->addGlobalParameter("eps", 25.894776);
    force->addGlobalParameter("a", 1.8);
    force->addGlobalParameter("sigma", 0.23925);
    force->addGlobalParameter("gamma", 1.2);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffNonPeriodic);
    force->setCutoffDistance(0.155);
    vector<double> params;
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    force->addParticle(params);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(0.1, 0, 0));
    positions.push_back(Vec3(0, 0.11, 0.03));
    positions.push_back(Vec3(0.04, 0, -0.08));
    int sets[8][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,0,2}, {1,0,3}, {1, 2, 3}, {2,0,1}, {3,0,1}};
    vector<const int*> expectedSets(&sets[0], &sets[8]);
    validateStillingerWeber(force, positions, expectedSets, 2.0);
}

void testCentralParticleModeLargeSystem() {
    int gridSize = 8;
    int numParticles = gridSize*gridSize*gridSize;
    double boxSize = 2.0;
    double spacing = boxSize/gridSize;
    CustomManyParticleForce* force = new CustomManyParticleForce(3,
        "L*eps*(cos(theta1)+1/3)^2*exp(sigma*gamma/(r12-a*sigma))*exp(sigma*gamma/(r13-a*sigma));"
        "r12 = distance(p1,p2); r13 = distance(p1,p3); theta1 = angle(p3,p1,p2)");
    force->setPermutationMode(CustomManyParticleForce::UniqueCentralParticle);
    force->addGlobalParameter("L", 23.13);
    force->addGlobalParameter("eps", 25.894776);
    force->addGlobalParameter("a", 1.8);
    force->addGlobalParameter("sigma", 0.23925);
    force->addGlobalParameter("gamma", 1.2);
    force->setNonbondedMethod(CustomManyParticleForce::CutoffPeriodic);
    force->setCutoffDistance(1.8*0.23925);
    vector<double> params;
    vector<Vec3> positions;
    System system;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                force->addParticle(params);
                positions.push_back(Vec3((i+0.4*genrand_real2(sfmt))*spacing, (j+0.4*genrand_real2(sfmt))*spacing, (k+0.4*genrand_real2(sfmt))*spacing));
                system.addParticle(1.0);
            }
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(force);
    VerletIntegrator integrator1(0.001);
    VerletIntegrator integrator2(0.001);
    Context context1(system, integrator1, Platform::getPlatformByName("Reference"));
    Context context2(system, integrator2, platform);
    context1.setPositions(positions);
    context2.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    CustomManyParticleForce* force = new CustomManyParticleForce(2, "x1+y2+none");
    force->addParticle();
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

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testNoCutoff();
        testCutoff();
        testPeriodic();
        testTriclinic();
        testExclusions();
        testAllTerms();
        testParameters();
        testTabulatedFunctions();
        testTypeFilters();
        testLargeSystem();
        testCentralParticleModeNoCutoff();
        testCentralParticleModeCutoff();
        testCentralParticleModeLargeSystem();
        testIllegalVariable();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
