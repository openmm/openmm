/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2023 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Emilio Gallicchio                                  *
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
#include "openmm/ATMForce.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include "sfmt/SFMT.h"
#include <algorithm>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void test2Particles() {
    // A pair of particles tethered by an harmonic bond.
    // Displace the second one to test energy and forces at different lambda values

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    double lambda1 = 0.5;
    double lambda2 = 0.5;
    double alpha = 0.0;
    double u0 = 0.0;
    double w0 = 0.0;
    double umax = 1.e6;
    double ubcore= 0.5*1.e6;
    double acore = 1./16.;
    double direction = 1.0;

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);

    ATMForce* atm = new ATMForce(lambda1, lambda2, alpha, u0, w0, umax, ubcore, acore, direction);
    vector<double> nodispl({0., 0., 0.});
    vector<double>   displ({1., 0., 0.});
    atm->addParticle(nodispl);
    atm->addParticle(  displ);
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    for (double lmbd : {0.0, 0.5, 1.0}) {
        context.setParameter(ATMForce::Lambda1(), lmbd);
        context.setParameter(ATMForce::Lambda2(), lmbd);
        State state = context.getState(State::Energy | State::Forces);
        double epot = state.getPotentialEnergy();
        double epert = atm->getPerturbationEnergy(context);

        ASSERT_EQUAL_TOL(lmbd, context.getParameter(atm->Lambda1()), 1e-6);
        ASSERT_EQUAL_TOL(lmbd, context.getParameter(atm->Lambda2()), 1e-6);
        ASSERT_EQUAL_TOL(lmbd*0.5*displ[0]*displ[0], epot, 1e-6);
        ASSERT_EQUAL_TOL(0.5*displ[0]*displ[0], epert, 1e-6);
        ASSERT_EQUAL_VEC(Vec3(-lmbd*displ[0], 0.0, 0.0), state.getForces()[1], 1e-6);
    }
    
}

void test2Particles2Displacements() {
    // A pair of particles tethered by an harmonic bond. 
    // Displace the second one to test energy and forces at different lambda values
    // In this version the second particle is displaced in both the initial and final states
    // by different amounts.
  
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    double lambda1 = 0.5;
    double lambda2 = 0.5;
    double alpha = 0.0;
    double u0 = 0.0;
    double w0 = 0.0;
    double umax = 1.e6;
    double ubcore= 0.5*1.e6;
    double acore = 1./16.;
    double direction = 1.0;

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);
    
    ATMForce* atm = new ATMForce(lambda1, lambda2, alpha, u0, w0, umax, ubcore, acore, direction);
    //first particle is not displaced at either state
    vector<double> nodispl({0., 0., 0.});
    atm->addParticle( nodispl );
    //second particle is displaced at both states but by the same amount (1,0,0)
    vector<double> displ({1., 0., 0.,   1., 0., 0.});
    atm->addParticle( displ );
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state;
    double epot, epert;
    
    // U = U0 + lambda*epert; epert = U1 - U0

    // When the second particle is displaced by the same amount,
    // the perturbation energy should be zero since the second particle
    // is at the same position in the target and initial states,
    // and potential energy should be U0, the energy of the bond with the
    // second particle displaced
    state = context.getState(State::Energy | State::Forces);
    epot = state.getPotentialEnergy();
    epert = atm->getPerturbationEnergy(context);
    double d = displ[0];
    ASSERT_EQUAL_TOL(0.5*d*d, epot, 1e-6);
    ASSERT_EQUAL_TOL(0.0, epert, 1e-6);

    //Displace the second particle further in the target state
    displ[0] += 1.;
    atm->setParticleParameters(1, displ );
    atm->updateParametersInContext(context);
    state = context.getState(State::Energy | State::Forces);
    epot = state.getPotentialEnergy();
    epert = atm->getPerturbationEnergy(context);
    double d1 = displ[0];
    double d0 = displ[3];
    ASSERT_EQUAL_TOL(0.5*d1*d1 - 0.5*d0*d0, epert, 1e-6);
    ASSERT_EQUAL_TOL(0.5*d0*d0 + lambda2*epert, epot, 1e-6);
}



void testLargeSystem() {
    // Create a system with lots of particles, each displaced differently.
    
    int numParticles = 1000;
    System system;
    CustomExternalForce* external = new CustomExternalForce("x^2 + 2*y^2 + 3*z^2");
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
    atm->addForce(external);
    system.addForce(atm);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions, displacements;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions.push_back(3*Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
        Vec3 d(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5);
        displacements.push_back(d);
        external->addParticle(i);
	vector<double> displ({d[0], d[1], d[2]});
        atm->addParticle( displ );
    }

    // Also add a nonbonded force to trigger atom reordering on the GPU.

    CustomNonbondedForce* nb = new CustomNonbondedForce("a*r^2");
    nb->addGlobalParameter("a", 0.0);
    for (int i = 0; i < numParticles; i++)
        nb->addParticle();
    system.addForce(nb);
    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    // Evaluate the forces to see if the particles are at the correct positions.

    for (double lambda : {0.0, 1.0}) {
        context.setParameter(ATMForce::Lambda1(), lambda);
        context.setParameter(ATMForce::Lambda2(), lambda);
        State state = context.getState(State::Forces);
        for (int i = 0; i < numParticles; i++) {
            Vec3 expectedPos = positions[i] + lambda*displacements[i];
            Vec3 expectedForce(-2*expectedPos[0], -4*expectedPos[1], -6*expectedPos[2]);
            ASSERT_EQUAL_VEC(expectedForce, state.getForces()[i], 1e-6);
        }
    }
}

void testNonbonded() {
    // Create a system with lots of particles, each displaced differently.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    
    NonbondedForce* nbforce = new NonbondedForce();
    nbforce->addParticle( 1.0, 1.0, 1.0);
    nbforce->addParticle(-1.0, 1.0, 1.0);

    system.addForce(nbforce);

    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
    vector<double> nodispl({0., 0., 0.});
    vector<double>  displ({1., 0., 0.});
    atm->addParticle( nodispl );
    atm->addParticle(   displ );
    //atm->addForce(nbforce);
    atm->addForce(XmlSerializer::clone<Force>(*nbforce));
    system.removeForce(0);
    system.addForce(atm);

    vector<Vec3> positions;
    positions.push_back(Vec3(0., 0., 0.));
    positions.push_back(Vec3(1., 0., 0.));

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    
    double lambda = 0.5;
    context.setParameter(ATMForce::Lambda1(), lambda);
    context.setParameter(ATMForce::Lambda2(), lambda);
    State state = context.getState( State::Energy );
    double epot = state.getPotentialEnergy();
    double epert = atm->getPerturbationEnergy(context);
    ASSERT_EQUAL_TOL(-104.2320, epot,  1e-3);
    ASSERT_EQUAL_TOL(  69.4062, epert, 1e-3);
    //std::cout << "Nonbonded: epot = " << epot << std::endl;
    //std::cout << "Nonbonded: epert = " << epert << std::endl;
    
}


void testMolecules() {
    // Verify that ATMForce correctly propagates information about molecules
    // from the forces it contains.
    
    System system;
    for (int i = 0; i < 5; i++)
        system.addParticle(1.0);
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
    vector<double> nodispl({0., 0., 0.});
    for (int i = 0; i < 5; i++)
      atm->addParticle( nodispl );
    system.addForce(atm);
    HarmonicBondForce* bonds1 = new HarmonicBondForce();
    bonds1->addBond(0, 1, 1.0, 1.0);
    bonds1->addBond(2, 3, 1.0, 1.0);
    atm->addForce(bonds1);
    HarmonicBondForce* bonds2 = new HarmonicBondForce();
    bonds2->addBond(1, 2, 1.0, 1.0);
    atm->addForce(bonds2);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    vector<vector<int> > molecules = context.getMolecules();
    ASSERT_EQUAL(2, molecules.size());
    for (auto& mol : molecules) {
        if (mol.size() == 1) {
            ASSERT_EQUAL(4, mol[0]);
        }
        else {
            ASSERT_EQUAL(4, mol.size());
            for (int i = 0; i < 4; i++)
                ASSERT(find(mol.begin(), mol.end(), i) != mol.end());
        }
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        test2Particles();
	test2Particles2Displacements();
        testLargeSystem();
	testNonbonded();
        testMolecules();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
