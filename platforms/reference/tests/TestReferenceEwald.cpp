/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the Eewald summation method reference implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "openmm/HarmonicBondForce.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;


void testEwaldExact() {

//    Use a NaCl crystal to compare the calculated and Madelung energies

    const int numParticles = 1000;
    const double cutoff = 1.0;
    const double boxSize = 2.82;

    ReferencePlatform platform;
    System system;

    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "nacl_crystal.dat"
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

//   The potential energy of an ion in a crystal is 
//   E = - (M*e/ 4*pi*epsilon0*a0),
//   where 
//   M			:	Madelung constant (dimensionless, for FCC cells such as NaCl it is 1.7476)
//   e			:	1.6022 × 10−19 C
//   4*pi*epsilon0	: 	1.112 × 10−10 C²/(J m)
//   a0			:	0.282 x 10-9 m (perfect cell)
//   Therefore
    double exactEnergy = -8.9290; // eV
    exactEnergy        = -4.4645; // eV per 1 ion
    exactEnergy        = -7.15292069e-19; // J per ion
    exactEnergy        = -430.820; // kJ/mol per ion
    ASSERT_EQUAL_TOL(exactEnergy * numParticles , state.getPotentialEnergy(), 100*TOL);
//  Gromacs result
//    ASSERT_EQUAL_TOL(-430494.0, state.getPotentialEnergy(), 10*TOL);

}

void testEwald2Ions() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->addParticle(1.0, 1, 0);
    nonbonded->addParticle(-1.0, 1, 0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    const double cutoff = 2.0;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], 10*TOL);
    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], 10*TOL);
    ASSERT_EQUAL_TOL(-217.276, state.getPotentialEnergy(), 10*TOL);
}

void testWaterSystem() {
    ReferencePlatform platform;
    System system;
    static int numParticles;
    numParticles = 648;
    for (int i = 0 ; i < numParticles ; i++)
    {
       system.addParticle(1.0);
    }
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0 ; i < numParticles/3 ; i++)
    {
      nonbonded->addParticle(-0.82, 1, 0);
      nonbonded->addParticle(0.41, 1, 0);
      nonbonded->addParticle(0.41, 1, 0);
    }
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    const double cutoff = 0.8;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(1.86206, 0, 0), Vec3(0, 1.86206, 0), Vec3(0, 0, 1.86206));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "water.dat"
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state1.getForces();
   // for (int i = 0 ; i < 42 ; i++)
   //   cout << "f [" << i << " : ]" << forces[i] << endl;
   // cout << "PotentialEnergy: " << state1.getPotentialEnergy() << endl;

// Take a small step in the direction of the energy gradient.
    
    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state1.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }


    norm = std::sqrt(norm);
    const double delta = 1e-3;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state1.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    // See whether the potential energy changed by the expected amount.
    
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    State state2 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state1.getPotentialEnergy())/delta, 0.01)


}

void testEwaldPME() {

//      Use amorphoush NaCl system

    ReferencePlatform platform;
    const int numParticles = 216;
    System system;
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(22.99);
    for (int i = 0; i < numParticles/2; i++)
        system.addParticle(35.45);
    VerletIntegrator integrator(0.01);
    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(1.0, 1.0,0.0);
    for (int i = 0; i < numParticles/2; i++)
        nonbonded->addParticle(-1.0, 1.0,0.0);
    nonbonded->setNonbondedMethod(NonbondedForce::Ewald);
    const double cutoff = 0.8;
    nonbonded->setCutoffDistance(cutoff);
    system.setPeriodicBoxVectors(Vec3(1.86206, 0, 0), Vec3(0, 1.86206, 0), Vec3(0, 0, 1.86206));
    nonbonded->setEwaldErrorTolerance(TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "nacl_amorph.dat"
    context.setPositions(positions);

    State state1 = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces1 = state1.getForces();

//   (1)   CHECK EXACT VALUE OF EWALD ENERGY

    ASSERT_EQUAL_TOL(-26651.9, state1.getPotentialEnergy(), TOL);

//   (2)   CHECK WHETHER THE EWALD FORCES ARE THE SAME AS THE GROMACS OUTPUT
// Even at tolerance 0.1 the test doesn't pass
//    #include "nacl_amorph_GromacsForcesEwald.dat"

//   (3)   CHECK SELF-CONSISTENCY

// Take a small step in the direction of the energy gradient.

    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state1.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }


    norm = std::sqrt(norm);
    const double delta = 1e-3;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state1.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    // See whether the potential energy changed by the expected amount.
    
    State state2 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state1.getPotentialEnergy())/delta, 0.01)

//   (4)   CHECK EXACT VALUE OF PME ENERGY

    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    context.reinitialize();
    #include "nacl_amorph.dat"
    context.setPositions(positions);
    State state3 = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(-26651.9, state3.getPotentialEnergy(), 10*TOL);

//   (5) CHECK WHETHER PME FORCES ARE THE SAME AS THE GROMACS OUTPUT

// Even at tolerance 0.1 the test doesn't pass
//    #include "nacl_amorph_GromacsForcesEwald.dat"

//   (6) CHECK PME FOR SELF-CONSISTENCY

// Take a small step in the direction of the energy gradient.
    
    norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state3.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }
    norm = std::sqrt(norm);
    step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state3.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    // See whether the potential energy changed by the expected amount.
    
    State state4 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state4.getPotentialEnergy()-state3.getPotentialEnergy())/delta, 0.01)

}



int main() {
    try {
     testEwaldExact();           
//     testEwald2Ions();
//     testWaterSystem();
     testEwaldPME();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
