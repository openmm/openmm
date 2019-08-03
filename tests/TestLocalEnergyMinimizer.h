
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2018 Stanford University and the Authors.      *
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
#include "openmm/HarmonicBondForce.h"
#include "openmm/LocalEnergyMinimizer.h"
#include "openmm/NonbondedForce.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VirtualSite.h"
#include "sfmt/SFMT.h"
#include <algorithm>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testHarmonicBonds() {
    const int numParticles = 10;
    System system;
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);

    // Create a chain of particles connected by harmonic bonds.

    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i] = Vec3(i, 0, 0);
        if (i > 0)
            bonds->addBond(i-1, i, 1+0.1*i, 1);
    }

    // Minimize it and check that all bonds are at their equilibrium distances.

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    LocalEnergyMinimizer::minimize(context, 1e-5);
    State state = context.getState(State::Positions);
    for (int i = 1; i < numParticles; i++) {
        Vec3 delta = state.getPositions()[i]-state.getPositions()[i-1];
        ASSERT_EQUAL_TOL(1+0.1*i, sqrt(delta.dot(delta)), 1e-4);
    }
}

void testLargeSystem() {
    const int numMolecules = 25;
    const int numParticles = numMolecules*2;
    const double cutoff = 2.0;
    const double boxSize = 4.0;
    const double tolerance = 15;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(nonbonded);

    // Create a cloud of molecules.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        nonbonded->addParticle(-1.0, 0.2, 0.2);
        nonbonded->addParticle(1.0, 0.2, 0.2);
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        system.addConstraint(2*i, 2*i+1, 1.0);
    }

    // Minimize it and verify that the energy has decreased.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State initialState = context.getState(State::Forces | State::Energy);
    LocalEnergyMinimizer::minimize(context, tolerance);
    State finalState = context.getState(State::Forces | State::Energy | State::Positions);
    ASSERT(finalState.getPotentialEnergy() < initialState.getPotentialEnergy());

    // Compute the force magnitude, subtracting off any component parallel to a constraint, and
    // check that it satisfies the requested tolerance.

    double forceNorm = 0.0;
    for (int i = 0; i < numParticles; i += 2) {
        Vec3 dir = finalState.getPositions()[i+1]-finalState.getPositions()[i];
        double distance = sqrt(dir.dot(dir));
        dir *= 1.0/distance;
        Vec3 f = finalState.getForces()[i];
        f -= dir*dir.dot(f);
        forceNorm += f.dot(f);
        f = finalState.getForces()[i+1];
        f -= dir*dir.dot(f);
        forceNorm += f.dot(f);
    }
    forceNorm = sqrt(forceNorm/(5*numMolecules));
    ASSERT(forceNorm < 2*tolerance);
}

void testVirtualSites() {
    const int numMolecules = 25;
    const int numParticles = numMolecules*3;
    const double cutoff = 2.0;
    const double boxSize = 4.0;
    const double tolerance = 10;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setCutoffDistance(cutoff);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    system.addForce(nonbonded);

    // Create a cloud of molecules.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(0.0);
        nonbonded->addParticle(-1.0, 0.2, 0.2);
        nonbonded->addParticle(0.5, 0.2, 0.2);
        nonbonded->addParticle(0.5, 0.2, 0.2);
        positions[3*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[3*i+1] = Vec3(positions[3*i][0]+1.0, positions[3*i][1], positions[3*i][2]);
        positions[3*i+2] = Vec3();
        system.addConstraint(3*i, 3*i+1, 1.0);
        system.setVirtualSite(3*i+2, new TwoParticleAverageSite(3*i, 3*i+1, 0.5, 0.5));
    }

    // Minimize it and verify that the energy has decreased.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.applyConstraints(1e-5);
    State initialState = context.getState(State::Forces | State::Energy);
    LocalEnergyMinimizer::minimize(context, tolerance);
    State finalState = context.getState(State::Forces | State::Energy | State::Positions);
    ASSERT(finalState.getPotentialEnergy() < initialState.getPotentialEnergy());

    // Compute the force magnitude, subtracting off any component parallel to a constraint, and
    // check that it satisfies the requested tolerance.

    double forceNorm = 0.0;
    for (int i = 0; i < numParticles; i += 3) {
        Vec3 dir = finalState.getPositions()[i+1]-finalState.getPositions()[i];
        double distance = sqrt(dir.dot(dir));
        dir *= 1.0/distance;
        Vec3 f = finalState.getForces()[i];
        f -= dir*dir.dot(f);
        forceNorm += f.dot(f);
        f = finalState.getForces()[i+1];
        f -= dir*dir.dot(f);
        forceNorm += f.dot(f);
        
        // Check the virtual site location.
        
        ASSERT_EQUAL_VEC((finalState.getPositions()[i+1]+finalState.getPositions()[i])*0.5, finalState.getPositions()[i+2], 1e-5);
    }
    forceNorm = sqrt(forceNorm/(5*numMolecules));
    ASSERT(forceNorm < 2*tolerance);
}

void testLargeForces() {
    // Create a set of particles that are almost on top of each other so the initial
    // forces are huge.
    
    const int numParticles = 10;
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    system.addForce(nonbonded);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(1.0, 0.2, 1.0);
    }
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*1e-10;

    // Minimize it and verify that it didn't blow up.                                                                               

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    LocalEnergyMinimizer::minimize(context, 1.0);
    State state = context.getState(State::Positions);
    double maxdist = 0.0;
    for (int i = 0; i < numParticles; i++) {
        Vec3 r = state.getPositions()[i];
        maxdist = max(maxdist, sqrt(r.dot(r)));
    }
    ASSERT(maxdist > 0.1);
    ASSERT(maxdist < 10.0);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testHarmonicBonds();
        testLargeSystem();
        testVirtualSites();
        testLargeForces();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
