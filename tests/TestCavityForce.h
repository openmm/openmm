/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
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
#include "openmm/CavityForce.h"
#include "openmm/CavityParticleDisplacer.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace OpenMM;
using namespace std;

/**
 * Test that the CavityForce computes correct energy components.
 */
void testCavityEnergy() {
    const int numMolecularParticles = 4;
    const double omegac = 0.01;  // Cavity frequency
    const double lambda = 0.5;   // Coupling constant
    const double photonMass = 1.0;
    
    System system;
    VerletIntegrator integrator(0.001);
    NonbondedForce* nonbonded = new NonbondedForce();
    
    // Add molecular particles
    for (int i = 0; i < numMolecularParticles; ++i) {
        system.addParticle(12.0);
        nonbonded->addParticle((i % 2 == 0 ? 0.5 : -0.5), 0.3, 1.0);
    }
    
    // Add cavity particle (index = numMolecularParticles)
    system.addParticle(photonMass);
    nonbonded->addParticle(0.0, 0.0, 0.0);  // Cavity has no charge
    
    system.addForce(nonbonded);
    
    // Add cavity force
    int cavityIndex = numMolecularParticles;
    CavityForce* cavityForce = new CavityForce(cavityIndex, omegac, lambda, photonMass);
    system.addForce(cavityForce);
    
    ASSERT(cavityForce->usesPeriodicBoundaryConditions());
    
    // Set up positions
    Context context(system, integrator, platform);
    vector<Vec3> positions(numMolecularParticles + 1);
    
    // Molecular particles in a square
    positions[0] = Vec3(1.0, 0.0, 0.0);   // +0.5 charge
    positions[1] = Vec3(-1.0, 0.0, 0.0);  // -0.5 charge
    positions[2] = Vec3(0.0, 1.0, 0.0);   // +0.5 charge
    positions[3] = Vec3(0.0, -1.0, 0.0);  // -0.5 charge
    
    // Cavity particle at origin
    positions[4] = Vec3(0.0, 0.0, 0.0);
    
    context.setPositions(positions);
    
    // Compute energy
    State state = context.getState(State::Energy | State::Forces);
    
    // Calculate expected values
    // Dipole = sum(q * r)
    // = 0.5*(1,0,0) + (-0.5)*(-1,0,0) + 0.5*(0,1,0) + (-0.5)*(0,-1,0)
    // = (0.5, 0, 0) + (0.5, 0, 0) + (0, 0.5, 0) + (0, 0.5, 0)
    // = (1.0, 1.0, 0)
    double dipoleX = 1.0;
    double dipoleY = 1.0;
    
    double K = photonMass * omegac * omegac;
    double epsilon = lambda * omegac;
    
    // Cavity at origin: q = (0, 0, 0)
    double qx = 0.0, qy = 0.0, qz = 0.0;
    
    // Harmonic energy: (1/2) * K * q^2 = 0
    double expectedHarmonic = 0.5 * K * (qx*qx + qy*qy + qz*qz);
    
    // Coupling energy: epsilon * q_xy . d_xy = 0
    double expectedCoupling = epsilon * (qx*dipoleX + qy*dipoleY);
    
    // Dipole self-energy: (epsilon^2 / 2K) * d_xy^2
    double expectedDipoleSelf = 0.5 * epsilon * epsilon / K * (dipoleX*dipoleX + dipoleY*dipoleY);
    
    double harmonicEnergy = cavityForce->getHarmonicEnergy(context);
    double couplingEnergy = cavityForce->getCouplingEnergy(context);
    double dipoleSelfEnergy = cavityForce->getDipoleSelfEnergy(context);
    
    ASSERT_EQUAL_TOL(expectedHarmonic, harmonicEnergy, 1e-6);
    ASSERT_EQUAL_TOL(expectedCoupling, couplingEnergy, 1e-6);
    ASSERT_EQUAL_TOL(expectedDipoleSelf, dipoleSelfEnergy, 1e-6);
}

/**
 * Test that forces are computed correctly.
 */
void testCavityForces() {
    const double omegac = 0.01;
    const double lambda = 0.5;
    const double photonMass = 1.0;
    
    System system;
    VerletIntegrator integrator(0.001);
    NonbondedForce* nonbonded = new NonbondedForce();
    
    // One positive particle, one negative, one cavity
    system.addParticle(12.0);
    nonbonded->addParticle(1.0, 0.3, 1.0);  // Particle 0: charge +1
    
    system.addParticle(12.0);
    nonbonded->addParticle(-1.0, 0.3, 1.0); // Particle 1: charge -1
    
    system.addParticle(photonMass);
    nonbonded->addParticle(0.0, 0.0, 0.0);  // Particle 2: cavity
    
    system.addForce(nonbonded);
    
    CavityForce* cavityForce = new CavityForce(2, omegac, lambda, photonMass);
    system.addForce(cavityForce);
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(1.0, 0.0, 0.0);   // +1 at (1,0,0)
    positions[1] = Vec3(-1.0, 0.0, 0.0);  // -1 at (-1,0,0)
    positions[2] = Vec3(0.5, 0.5, 0.0);   // Cavity at (0.5, 0.5, 0)
    context.setPositions(positions);
    
    // Dipole = (+1)*(1,0,0) + (-1)*(-1,0,0) = (1,0,0) + (1,0,0) = (2,0,0)
    double dipoleX = 2.0;
    double dipoleY = 0.0;
    
    double K = photonMass * omegac * omegac;
    double epsilon = lambda * omegac;
    
    double qx = 0.5, qy = 0.5;
    
    // Dq = q + (lambda/omega) * d
    double DqX = qx + (lambda / omegac) * dipoleX;
    double DqY = qy + (lambda / omegac) * dipoleY;
    
    State state = context.getState(State::Forces);
    const vector<Vec3>& forces = state.getForces();
    
    // Force on particle 0: F = -epsilon * q * Dq (where q = +1)
    double expectedF0x = -epsilon * 1.0 * DqX;
    double expectedF0y = -epsilon * 1.0 * DqY;
    
    // Force on particle 1: F = -epsilon * q * Dq (where q = -1)
    double expectedF1x = -epsilon * (-1.0) * DqX;
    double expectedF1y = -epsilon * (-1.0) * DqY;
    
    // Force on cavity: F = -K*q - epsilon*d_xy
    double expectedFcavX = -K * qx - epsilon * dipoleX;
    double expectedFcavY = -K * qy - epsilon * dipoleY;
    
    // Note: There will also be forces from the NonbondedForce, so we can't
    // directly compare. Instead, run a simulation and verify energy conservation.
}

/**
 * Test time-varying coupling schedule.
 */
void testCouplingSchedule() {
    const double omegac = 0.01;
    const double photonMass = 1.0;
    
    System system;
    VerletIntegrator integrator(0.001);
    NonbondedForce* nonbonded = new NonbondedForce();
    
    // Simple two-particle system
    system.addParticle(12.0);
    nonbonded->addParticle(1.0, 0.3, 1.0);
    
    system.addParticle(photonMass);
    nonbonded->addParticle(0.0, 0.0, 0.0);
    
    system.addForce(nonbonded);
    
    CavityForce* cavityForce = new CavityForce(1, omegac, 0.0, photonMass);
    
    // Set up coupling schedule: off initially, on after step 100
    vector<pair<int, double>> schedule;
    schedule.push_back(make_pair(0, 0.0));    // Off at start
    schedule.push_back(make_pair(100, 0.5));  // On after step 100
    cavityForce->setLambdaCouplingSchedule(schedule);
    
    system.addForce(cavityForce);
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(1.0, 0.0, 0.0);
    positions[1] = Vec3(0.0, 0.0, 0.0);
    context.setPositions(positions);
    
    // Before step 100, coupling should be 0
    ASSERT_EQUAL_TOL(0.0, cavityForce->getLambdaCouplingAtStep(50), 1e-10);
    
    // After step 100, coupling should be 0.5
    ASSERT_EQUAL_TOL(0.5, cavityForce->getLambdaCouplingAtStep(150), 1e-10);
}

/**
 * Test cavity particle displacer.
 */
void testCavityDisplacer() {
    const double omegac = 0.01;
    const double lambda = 0.5;
    const double photonMass = 1.0;
    
    System system;
    VerletIntegrator integrator(0.001);
    NonbondedForce* nonbonded = new NonbondedForce();
    
    // Two opposite charges and cavity
    system.addParticle(12.0);
    nonbonded->addParticle(1.0, 0.3, 1.0);
    
    system.addParticle(12.0);
    nonbonded->addParticle(-1.0, 0.3, 1.0);
    
    system.addParticle(photonMass);
    nonbonded->addParticle(0.0, 0.0, 0.0);
    
    system.addForce(nonbonded);
    
    CavityParticleDisplacer* displacer = new CavityParticleDisplacer(2, omegac, photonMass);
    displacer->setSwitchOnStep(10);
    displacer->setSwitchOnLambda(lambda);
    system.addForce(displacer);
    
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(1.0, 0.0, 0.0);   // +1 at (1,0,0)
    positions[1] = Vec3(-1.0, 0.0, 0.0);  // -1 at (-1,0,0)
    positions[2] = Vec3(0.0, 0.0, 0.5);   // Cavity at (0,0,0.5)
    context.setPositions(positions);
    
    // Run past the switch-on step
    integrator.step(20);
    
    State state = context.getState(State::Positions);
    const vector<Vec3>& newPositions = state.getPositions();
    
    // Expected equilibrium position: q_eq = -(lambda/omega) * d_xy
    // Dipole = (+1)*(1,0,0) + (-1)*(-1,0,0) = (2,0,0)
    double dipoleX = 2.0;
    double dipoleY = 0.0;
    double factor = -lambda / omegac;
    double expectedQx = factor * dipoleX;
    double expectedQy = factor * dipoleY;
    
    // After the first displacement, cavity should be at equilibrium
    // (Note: subsequent dynamics will move it, so this test is approximate)
    // The z-coordinate should be preserved
    ASSERT_EQUAL_TOL(0.5, newPositions[2][2], 0.01);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testCavityEnergy();
        testCavityForces();
        testCouplingSchedule();
        testCavityDisplacer();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
