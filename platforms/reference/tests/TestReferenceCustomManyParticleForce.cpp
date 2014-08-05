/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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
 * This tests the reference implementation of CustomManyParticleForce.
 */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void validateAxilrodTeller(CustomManyParticleForce* force, const vector<Vec3>& positions, const vector<const int*>& expectedSets, double boxSize) {
    // Create a System and Context.
    
    int numParticles = force->getNumParticles();
    CustomManyParticleForce::NonbondedMethod nonbondedMethod = force->getNonbondedMethod();
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    ReferencePlatform platform;
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    double c = context.getParameter("C");
    
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
    validateAxilrodTeller(force, positions, expectedSets, 2.0);
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
    validateAxilrodTeller(force, positions, expectedSets, 2.0);
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
    validateAxilrodTeller(force, positions, expectedSets, boxSize);
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
    ReferencePlatform platform;
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // See if the energy is correct.

    State state1 = context.getState(State::Energy);
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
    ASSERT_EQUAL_TOL(expectedEnergy, state1.getPotentialEnergy(), 1e-5);
}

int main() {
    try {
        testNoCutoff();
        testCutoff();
        testPeriodic();
        testParameters();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
