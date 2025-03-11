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
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/LangevinMiddleIntegrator.h"
#include "openmm/LocalEnergyMinimizer.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include "sfmt/SFMT.h"
#include <algorithm>
#include <random>
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

    double lmbd = 0.5;
    double umax =  0.;
    double ubcore= 0.;
    double acore = 0.;
    double direction = 1.0;

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);

    ATMForce* atm = new ATMForce(lmbd, lmbd, 0., 0, 0, umax, ubcore, acore, direction);
    Vec3 nodispl = Vec3(0., 0., 0.);
    Vec3   displ = Vec3(1., 0., 0.);
    atm->addParticle( nodispl );
    atm->addParticle(   displ );
    atm->addForce(bond);
    atm->addEnergyParameterDerivative(ATMForce::Lambda1());
    atm->addEnergyParameterDerivative(ATMForce::Lambda2());
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    for (double lm : {0.0, 0.5, 1.0}) {
        context.setParameter(ATMForce::Lambda1(), lm);
        context.setParameter(ATMForce::Lambda2(), lm);
        State state = context.getState(State::Energy | State::Forces | State::ParameterDerivatives);
        double epot = state.getPotentialEnergy();
        double u0, u1, energy;
        atm->getPerturbationEnergy(context, u1, u0, energy);
        double epert = u1 - u0;
        ASSERT_EQUAL_TOL(lm, context.getParameter(atm->Lambda1()), 1e-6);
        ASSERT_EQUAL_TOL(lm, context.getParameter(atm->Lambda2()), 1e-6);
        ASSERT_EQUAL_TOL(energy, epot, 1e-6);
        ASSERT_EQUAL_TOL(lm*0.5*displ[0]*displ[0], epot, 1e-6);
        ASSERT_EQUAL_TOL(0.0, u0, 1e-6);
        ASSERT_EQUAL_TOL(0.5*displ[0]*displ[0], u1, 1e-6);
        ASSERT_EQUAL_TOL(0.5*displ[0]*displ[0], epert, 1e-6);
        ASSERT_EQUAL_VEC(Vec3(-lm*displ[0], 0.0, 0.0), state.getForces()[1], 1e-6);
        ASSERT_EQUAL_TOL(0.0, state.getEnergyParameterDerivatives().at(ATMForce::Lambda1()), 1e-6);
        ASSERT_EQUAL_TOL(0.5*displ[0]*displ[0], state.getEnergyParameterDerivatives().at(ATMForce::Lambda2()), 1e-6);
    }
}

void test2Particles2Displacement0() {
    // A pair of particles tethered by an harmonic bond. 
    // Displace the second one to test energy and forces at different lambda values
    // In this version the second particle is displaced in both the initial and final states
    // by different amounts.
  
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    double lmbd = 0.5;
    double umax = 0.;
    double ubcore= 0.;
    double acore = 0.;
    double direction = 1.0;

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);
    
    ATMForce* atm = new ATMForce(lmbd, lmbd, 0., 0., 0., umax, ubcore, acore, direction);
    //first particle is not displaced at either state
    Vec3 nodispl = Vec3(0., 0., 0.);
    atm->addParticle( nodispl );
    //second particle is displaced at both states but by the same amount (1,0,0)
    Vec3 displ0 = Vec3(1., 0., 0.);
    atm->addParticle( displ0, displ0 );
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state;
    double epot, epert;
    double u0, u1, energy;
    
    // U = U0 + lambda*epert; epert = U1 - U0

    // When the second particle is displaced by the same amount at each state,
    // the perturbation energy should be zero since the second particle
    // is at the same position in the target and initial states,
    // and the potential energy should be U0, the energy of the bond with the
    // second particle displaced
    state = context.getState(State::Energy | State::Forces);
    epot = state.getPotentialEnergy();
    atm->getPerturbationEnergy(context, u1, u0, energy);
    epert = u1 - u0;
    ASSERT_EQUAL_TOL(0.5*displ0[0]*displ0[0], epot, 1e-6);
    ASSERT_EQUAL_TOL(0.0, epert, 1e-6);

    //Displace the second particle further in the target state
    Vec3 displ1 = Vec3(2., 0., 0.);
    atm->setParticleParameters(1, displ1, displ0 );
    atm->updateParametersInContext(context);
    state = context.getState(State::Energy | State::Forces);
    epot = state.getPotentialEnergy();
    atm->getPerturbationEnergy(context, u1, u0, energy);
    epert = u1 - u0;
    ASSERT_EQUAL_TOL(0.5*displ1[0]*displ1[0] - 0.5*displ0[0]*displ0[0], epert, 1e-6);
    ASSERT_EQUAL_TOL(0.5*displ0[0]*displ0[0] + lmbd*epert, epot, 1e-6);
}


double softCoreFunc(double u, double umax, double ub, double a, double& df) {
    double usc = u;
    df = 1.;

    if(u > ub) {
        double gu = (u-ub)/(a*(umax-ub)); //this is y/alpha
        double zeta = 1. + 2.*gu*(gu + 1.) ;
        double zetap = pow( zeta , a);
        double s = 4.*(2.*gu + 1.)/zeta;
        df = s*zetap/pow(1.+zetap,2);
        usc = (umax-ub)*(zetap - 1.)/(zetap + 1.) + ub;
    }
    return usc;
}


void test2ParticlesSoftCore() {
    // Similar to test2Particles() but employing a soft-core function

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    double lmbd = 0.5;
    double umax =  10.;
    double ubcore= 3.;
    double acore = 0.125;
    double direction = 1.0;

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);

    ATMForce* atm = new ATMForce(lmbd, lmbd, 0., 0, 0, umax, ubcore, acore, direction);
    Vec3 nodispl = Vec3(0., 0., 0.);
    Vec3   displ = Vec3(5., 0., 0.);
    atm->addParticle( nodispl );
    atm->addParticle(   displ );
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Energy | State::Forces);
    double epot = state.getPotentialEnergy();
    double u0, u1, energy;
    atm->getPerturbationEnergy(context, u1, u0, energy);
    double epert = u1 - u0;
    double ee = 0.5*displ[0]*displ[0];
    double df;
    ASSERT_EQUAL_TOL(energy, epot, 1e-6);
    ASSERT_EQUAL_TOL(0.0, u0, 1e-6);
    ASSERT_EQUAL_TOL(ee,  u1, 1e-6);
    ASSERT_EQUAL_TOL(ee,  epert, 1e-6);
    ASSERT_EQUAL_TOL(u0 + lmbd*softCoreFunc(epert, umax, ubcore, acore, df), epot, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(-lmbd*df*displ[0], 0.0, 0.0), state.getForces()[1], 1e-6);
}


void testNonbonded() {
    // Tests a system with a nonbonded Force

    System system;
    double u0, u1, energy;
    double lambda = 0.5;
    double width = 4.0;

    system.setDefaultPeriodicBoxVectors(Vec3(width, 0, 0), Vec3(0, width, 0), Vec3(0, 0, width));
    NonbondedForce* nbforce = new NonbondedForce();
    nbforce->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nbforce->setCutoffDistance(0.7);
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    double spacing = width/6.0;
    double offset = spacing/5.0;
    vector<Vec3> positions;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++) {
	        positions.push_back(Vec3(spacing*i+offset, spacing*j+offset, spacing*k+offset));
		system.addParticle(10.0);
		nbforce->addParticle(0, 0.3, 1.0);
		atm->addParticle(Vec3());
            }
    auto rng = std::default_random_engine {};
    std::shuffle(std::begin(positions), std::end(positions), rng);
    atm->setParticleParameters(0, Vec3(0.5, 0, 0), Vec3(0.0, 0, 0));

    //in this scenario the non-bonded force is added to the System, a copy is added to ATMForce and
    //the System's copy is disabled by giving it a force group that is not evaluated.
    //This used to be needed to ensure atoms would be reordered.  It isn't anymore, but
    //this test is left in to make sure it still works.
    system.addForce(nbforce);
    atm->addForce(XmlSerializer::clone<Force>(*nbforce));
    nbforce->setForceGroup(1);
    system.addForce(atm);
    LangevinMiddleIntegrator integrator1(300, 1.0, 0.004);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    context1.setParameter(ATMForce::Lambda1(), lambda);
    context1.setParameter(ATMForce::Lambda2(), lambda);
    State state1 = context1.getState( State::Energy, false, 0 );
    double epot1 = state1.getPotentialEnergy();
    atm->getPerturbationEnergy(context1, u1, u0, energy);
    double epert1 = u1 - u0;

    //in this second scenario the non-bonded force is remove from the System
    system.removeForce(0);
    LangevinMiddleIntegrator integrator2(300, 1.0, 0.004);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    context2.setParameter(ATMForce::Lambda1(), lambda);
    context2.setParameter(ATMForce::Lambda2(), lambda);
    State state2 = context2.getState( State::Energy );
    double epot2 = state2.getPotentialEnergy();
    atm->getPerturbationEnergy(context2, u1, u0, energy);
    double epert2 = u1 - u0;
    ASSERT_EQUAL_TOL(epert1, epert2,  1e-3);
}


void testNonbondedwithEndpointClash() {
    // Similar to testNonbonded() but at the initial alchemical state
    // and with an invalid potential at the final state due to a clash
    // between two particles.

    System system;
    double u0, u1, energy;
    double lambda = 0.0; //U(lambda) = u0; it does not depend on u1
    double width = 4.0;

    system.setDefaultPeriodicBoxVectors(Vec3(width, 0, 0), Vec3(0, width, 0), Vec3(0, 0, width));
    NonbondedForce* nbforce = new NonbondedForce();
    nbforce->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nbforce->setCutoffDistance(0.7);
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    double spacing = width/6.0;
    double offset = spacing/5.0;
    vector<Vec3> positions;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++) {
	        positions.push_back(Vec3(spacing*i+offset, spacing*j+offset, spacing*k+offset));
		system.addParticle(10.0);
		nbforce->addParticle(0, 0.3, 1.0);
		atm->addParticle(Vec3(0,0,0));
            }
    //places first particle almost on top of another particle in displaced system
    atm->setParticleParameters(0, Vec3(spacing+1.e-4, 0, 0), Vec3(0.0, 0, 0));
    atm->addForce(nbforce);
    system.addForce(atm);
    LangevinMiddleIntegrator integrator(300, 1.0, 0.004);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setParameter(ATMForce::Lambda1(), lambda);
    context.setParameter(ATMForce::Lambda2(), lambda);

    State state = context.getState(State::Energy, false);
    double epot = state.getPotentialEnergy();
    ASSERT(!isnan(epot) && !isinf(epot));

    atm->getPerturbationEnergy(context, u1, u0, energy);
    double epert = u1 - u0;
    ASSERT(!isnan(energy) && !isinf(energy));
    ASSERT(!isnan(epert) && !isinf(epert));

    integrator.step(10);
    state = context.getState(State::Energy | State::Positions, false);
    vector<Vec3> positions2 = state.getPositions();
    ASSERT(fabs(positions[0][0] - positions2[0][0]) < width);
}

void testParticlesCustomExpressionLinear() {
    // Similar to test2Particles() but employing a custom alchemical energy expression

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);

    double lmbd = 0.5;
    ATMForce* atm = new ATMForce("u0 + Lambda*(u1 - u0)");
    atm->addGlobalParameter("Lambda", lmbd);
    Vec3 nodispl = Vec3(0., 0., 0.);
    Vec3   displ = Vec3(5., 0., 0.);
    atm->addParticle( nodispl );
    atm->addParticle(   displ );
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Energy | State::Forces);
    double epot = state.getPotentialEnergy();
    double u0, u1, energy;
    atm->getPerturbationEnergy(context, u1, u0, energy);
    double epert = u1 - u0;
    double ee = 0.5*displ[0]*displ[0];
    ASSERT_EQUAL_TOL(energy, epot, 1e-6);
    ASSERT_EQUAL_TOL(0.0, u0, 1e-6);
    ASSERT_EQUAL_TOL(ee,  u1, 1e-6);
    ASSERT_EQUAL_TOL(ee,  epert, 1e-6);
    ASSERT_EQUAL_TOL(u0 + lmbd*epert, epot, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(-lmbd*displ[0], 0.0, 0.0), state.getForces()[1], 1e-6);
}

void testParticlesCustomExpressionSoftplus() {
    // Similar to test2Particles() but employing a custom alchemical energy expression

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);

    Vec3 nodispl = Vec3(0., 0., 0.);
    Vec3   displ = Vec3(2., 0., 0.);

    CustomBondForce* bond = new CustomBondForce("0.5*r^2");
    bond->addBond(0, 1);

    ATMForce* atm = new ATMForce("u0 + ((Lambda2-Lambda1)/Alpha)*log(1.  + exp(-Alpha*((u1-u0) - Uh))) + Lambda2*(u1-u0) + W0");
    double lambda1 = 0.2;
    double lambda2 = 0.5;
    double alpha = 0.1;
    double uh = 0;
    double w0 = 0;

    atm->addGlobalParameter("Lambda1", lambda1);
    atm->addGlobalParameter("Lambda2", lambda2);
    atm->addGlobalParameter("Alpha", alpha);
    atm->addGlobalParameter("Uh", uh);
    atm->addGlobalParameter("W0", w0);

    atm->addParticle( nodispl );
    atm->addParticle(   displ );
    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Energy | State::Forces);
    double epot = state.getPotentialEnergy();
    double u0, u1, energy;
    atm->getPerturbationEnergy(context, u1, u0, energy);
    double epert = u1 - u0;

    double ebias = 0.0;
    double ee = 1.0 + exp(-alpha*(epert  - uh));
    if(alpha > 0){
      ebias = ((lambda2 - lambda1)/alpha) * log(ee);
    }
    ebias += lambda2 * epert  + w0;
    double bfp = (lambda2 - lambda1)/ee + lambda1;

    double ep = 0.5*displ[0]*displ[0];
    ASSERT_EQUAL_TOL(energy, epot, 1e-6);
    ASSERT_EQUAL_TOL(0.0, u0, 1e-6);
    ASSERT_EQUAL_TOL(ep,  u1, 1e-6);
    ASSERT_EQUAL_TOL(ep,  epert, 1e-6);
    ASSERT_EQUAL_TOL(u0 + ebias, epot, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(-bfp*displ[0], 0.0, 0.0), state.getForces()[1], 1e-6);
}

void testLargeSystem() {
    // Create a system with lots of particles, each displaced differently.
    
    int numParticles = 1000;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    CustomExternalForce* external = new CustomExternalForce("x^2 + 2*y^2 + 3*z^2");
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.1, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
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
        atm->addParticle(d);
    }

    // Also add nonbonded forces to trigger atom reordering on the GPU.

    CustomNonbondedForce* nb = new CustomNonbondedForce("a*r^2");
    nb->addGlobalParameter("a", 0.0);
    for (int i = 0; i < numParticles; i++)
        nb->addParticle();
    nb->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    system.addForce(nb);
    CustomNonbondedForce* nb1 = new CustomNonbondedForce("0");
    nb1->addPerParticleParameter("b");
    for (int i = 0; i < numParticles; i++)
        nb1->addParticle({(double) (i%3)});
    nb1->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    atm->addForce(nb1);
    
    // Evaluate the forces to see if the particles are at the correct positions.

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
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

void testMolecules() {
    // Verify that ATMForce correctly propagates information about molecules
    // from the forces it contains.
    
    System system;
    for (int i = 0; i < 5; i++)
        system.addParticle(1.0);
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.1, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
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

void testSimulation() {
    // Create a box of Lennard-Jones spheres, including an ATMForce that displaces
    // one particle to two different locations.
    
    int numParticles = 27;
    double width = 2.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(width, 0, 0), Vec3(0, width, 0), Vec3(0, 0, width));
    ATMForce* atm = new ATMForce("(u0+u1)/2");
    system.addForce(atm);
    NonbondedForce* nb = new NonbondedForce();
    nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nb->setCutoffDistance(1.0);
    atm->addForce(nb);
    vector<Vec3> positions;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) {
                system.addParticle(10.0);
                positions.push_back(Vec3(0.6*i, 0.6*j, 0.6*k));
                nb->addParticle(0, 0.3, 1.0);
                atm->addParticle(Vec3());
            }
    atm->setParticleParameters(0, Vec3(0.3, 0, 0), Vec3(-0.3, 0, 0));

    // Simulate it and make sure that the other particles avoid the displaced positions.

    LangevinMiddleIntegrator integrator(300, 1.0, 0.004);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(300);
    for (int i = 0; i < 100; i++) {
        integrator.step(10);
        vector<Vec3> pos = context.getState(State::Positions).getPositions();
        for (int j = 1; j < numParticles; j++) {
            for (double displacement : {-0.3, 0.3}) {
                Vec3 d = pos[0]-pos[j];
                d[0] += displacement;
                for (int k = 0; k < 3; k++)
                    d[k] -= round(d[k]/width)*width;
                assert(sqrt(d.dot(d)) > 0.2);
            }
        }
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        test2Particles();
        test2Particles2Displacement0();
        test2ParticlesSoftCore();
        testNonbonded();
	testNonbondedwithEndpointClash();
        testParticlesCustomExpressionLinear();
        testParticlesCustomExpressionSoftplus();
        testLargeSystem();
        testMolecules();
        testSimulation();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

