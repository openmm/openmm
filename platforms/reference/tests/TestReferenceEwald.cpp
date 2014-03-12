/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include "openmm/HarmonicBondForce.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double EWALD_TOL = 1e-5;
const double PME_TOL = 5e-5;

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
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setEwaldErrorTolerance(EWALD_TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "nacl_crystal.dat"
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

//   The potential energy of an ion in a crystal is 
//   E = - (M*e^2/ 4*pi*epsilon0*a0),
//   where 
//   M            :    Madelung constant (dimensionless, for FCC cells such as NaCl it is 1.7476)
//   e            :    1.6022 × 10−19 C
//   4*pi*epsilon0:     1.112 × 10−10 C²/(J m)
//   a0           :    0.282 x 10-9 m (perfect cell)
// 
//   E is then the energy per pair of ions, so for our case
//   E has to be divided by 2 (per ion), multiplied by N(avogadro), multiplied by number of particles, and divided by 1000 for kJ
    double exactEnergy        = - (1.7476 * 1.6022e-19 * 1.6022e-19  * 6.02214e+23 * numParticles) / (1.112e-10 * 0.282e-9 * 2 * 1000);
    //cout << "exact\t\t: " << exactEnergy << endl;
    //cout << "calc\t\t: " << state.getPotentialEnergy() << endl;
    ASSERT_EQUAL_TOL(exactEnergy, state.getPotentialEnergy(), 100*EWALD_TOL);

}

void testEwaldPME() {

    double tol = 1e-5;
    const double boxSize = 3.00646;
    const double cutoff = 1.2;
    const int numParticles = 894;

//      Use amorphous NaCl system
//      The particles are simple charges, no VdW interactions

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
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setEwaldErrorTolerance(EWALD_TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "nacl_amorph.dat"
    context.setPositions(positions);

    State state1 = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces1 = state1.getForces();

//   (1)   CHECK EXACT VALUE OF EWALD ENERGY (Against Gromacs output)

    tol = 1e-4;
    ASSERT_EQUAL_TOL(-3.82047e+05, state1.getPotentialEnergy(), tol);

//   (2)   CHECK WHETHER THE EWALD FORCES ARE THE SAME AS THE GROMACS OUTPUT
//         these are forces for alpha: 2.82756, kmax(x/y/z) = 11
  
    tol = 1e-2;
//    #include "nacl_amorph_GromacsForcesEwald.dat"

//   (3)   CHECK SELF-CONSISTENCY

    // Take a small step in the direction of the energy gradient.

    double norm = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state1.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    }


    norm = std::sqrt(norm);
    const double delta = 1e-2;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state1.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    // See whether the potential energy changed by the expected amount.
    
    State state2 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state1.getPotentialEnergy())/delta, fabs(EWALD_TOL*state2.getPotentialEnergy()/(state2.getPotentialEnergy()-state1.getPotentialEnergy())))

//   (4)   CHECK EXACT VALUE OF PME ENERGY 

    nonbonded->setNonbondedMethod(NonbondedForce::PME);
    nonbonded->setEwaldErrorTolerance(PME_TOL);
    context.reinitialize();
    #include "nacl_amorph.dat"
    context.setPositions(positions);
    State state3 = context.getState(State::Forces | State::Energy);

//  Gromacs PME energy for the same mesh
    tol = 1e-4;
    ASSERT_EQUAL_TOL(-3.82047e+05, state3.getPotentialEnergy(), tol);

//   (5) CHECK WHETHER PME FORCES ARE THE SAME AS THE GROMACS OUTPUT USING EWALD

    tol = 1e-1;
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
    ASSERT_EQUAL_TOL(norm, (state4.getPotentialEnergy()-state3.getPotentialEnergy())/delta, fabs(PME_TOL*state4.getPotentialEnergy()/(state4.getPotentialEnergy()-state3.getPotentialEnergy())))
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
    system.setDefaultPeriodicBoxVectors(Vec3(6, 0, 0), Vec3(0, 6, 0), Vec3(0, 0, 6));
    nonbonded->setEwaldErrorTolerance(EWALD_TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(3.048000,2.764000,3.156000);
    positions[1] = Vec3(2.809000,2.888000,2.571000);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(-123.711, 64.1877, -302.716), forces[0], 10*EWALD_TOL);
    ASSERT_EQUAL_VEC(Vec3(123.711, -64.1877, 302.716), forces[1], 10*EWALD_TOL);
    ASSERT_EQUAL_TOL(-217.276, state.getPotentialEnergy(), 10*EWALD_TOL);
}

void testWaterSystem() {
    ReferencePlatform platform;
    System system;
    static int numParticles = 648;
    const double boxSize = 1.86206;

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
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setEwaldErrorTolerance(EWALD_TOL);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    #include "water.dat"
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state1.getForces();

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

void testErrorTolerance(NonbondedForce::NonbondedMethod method) {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 5.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(method);
    ReferencePlatform platform;

    // For various values of the cutoff and error tolerance, see if the actual error is reasonable.

    for (double cutoff = 1.0; cutoff < boxWidth/2; cutoff *= 1.2) {
        force->setCutoffDistance(cutoff);
        vector<Vec3> refForces;
        double norm = 0.0;
        for (double tol = 5e-5; tol < 1e-3; tol *= 2.0) {
            force->setEwaldErrorTolerance(tol);
            VerletIntegrator integrator(0.01);
            Context context(system, integrator, platform);
            context.setPositions(positions);
            State state = context.getState(State::Forces);
            if (refForces.size() == 0) {
                refForces = state.getForces();
                for (int i = 0; i < numParticles; i++)
                    norm += refForces[i].dot(refForces[i]);
                norm = sqrt(norm);
            }
            else {
                double diff = 0.0;
                for (int i = 0; i < numParticles; i++) {
                    Vec3 delta = refForces[i]-state.getForces()[i];
                    diff += delta.dot(delta);
                }
                diff = sqrt(diff)/norm;
                ASSERT(diff < 2*tol);
            }
        }
    }
}

void testPMEParameters() {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 4.7;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::PME);
    ReferencePlatform platform;
    
    // Compute the energy with an error tolerance of 1e-3.

    force->setEwaldErrorTolerance(1e-3);
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    double energy1 = context1.getState(State::Energy).getPotentialEnergy();
    
    // Try again with an error tolerance of 1e-4.

    force->setEwaldErrorTolerance(1e-4);
    VerletIntegrator integrator2(0.01);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    double energy2 = context2.getState(State::Energy).getPotentialEnergy();
    
    // Now explicitly set the parameters.  These should match the values that were
    // used for tolerance 1e-3.

    force->setPMEParameters(2.49291157051793, 32, 32, 32);
    VerletIntegrator integrator3(0.01);
    Context context3(system, integrator3, platform);
    context3.setPositions(positions);
    double energy3 = context3.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(energy1, energy3, 1e-6);
    ASSERT(fabs((energy1-energy2)/energy1) > 1e-5);
}

int main() {
    try {
     testEwaldExact();
     testEwaldPME();
//     testEwald2Ions();
//     testWaterSystem();
     testErrorTolerance(NonbondedForce::Ewald);
     testErrorTolerance(NonbondedForce::PME);
     testPMEParameters();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
