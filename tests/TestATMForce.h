/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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
#include <string>

using namespace OpenMM;
using namespace std;

void testAPI(){
    ATMForce* atm = new ATMForce(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.6, 0.8, -1.0);
    atm->addParticle(Vec3(1, 2, 3), Vec3(4, 5, 6)); //old interface
    atm->addParticle(new ATMForce::FixedDisplacement(Vec3(7, 8, 9), Vec3(10, 11, 12)));
    atm->addParticle(new ATMForce::ParticleOffsetDisplacement(1, 0));
    atm->addParticle();

    Vec3 d1, d0;
    atm->getParticleParameters(0, d1, d0);
    ASSERT_EQUAL_VEC(Vec3(1, 2, 3), d1, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(4, 5, 6), d0, 1e-6);

    const ATMForce::FixedDisplacement* fd = (dynamic_cast<const ATMForce::FixedDisplacement*>(&(atm->getParticleTransformation(1))));
    d1 = fd->getFixedDisplacement1();
    d0 = fd->getFixedDisplacement0();
    ASSERT_EQUAL_VEC(Vec3(7, 8, 9), d1, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(10, 11, 12), d0, 1e-6);

    const ATMForce::ParticleOffsetDisplacement* vt = (dynamic_cast<const ATMForce::ParticleOffsetDisplacement*>(&(atm->getParticleTransformation(2))));
    int j1 = vt->getDestinationParticle1();
    int i1 = vt->getOriginParticle1();
    int j0 = vt->getDestinationParticle0();
    int i0 = vt->getOriginParticle0();
    ASSERT_EQUAL( 1, j1);
    ASSERT_EQUAL( 0, i1);
    ASSERT_EQUAL(-1, j0);
    ASSERT_EQUAL(-1, i0);

    atm->getParticleParameters(3, d1, d0);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), d1, 1e-6);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), d0, 1e-6);
}

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
    Vec3 displ = Vec3(1., 0., 0.);
    atm->addParticle();
    atm->addParticle(new ATMForce::FixedDisplacement(displ));
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

void test3ParticlesSwap() {
    // A pair of particles tethered by harmonic bonds to a central particle.
    // Swap the pair and test energy and forces at different lambda values

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    double lmbd = 0.5;
    double umax =  0.;
    double ubcore= 0.;
    double acore = 0.;
    double direction = 1.0;

    Vec3 origin = Vec3(0., 0., 0.);
    Vec3   r1 = Vec3(1., 0., 0.);
    double r1sq = r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
    Vec3   r2 = Vec3(-2., 0., 0.);
    double r2sq = r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2];

    vector<Vec3> positions(3);
    positions[0] = origin;
    positions[1] = r1;
    positions[2] = r2;

    CustomBondForce* bond = new CustomBondForce("0.5*kf*r^2");
    double kf1 = 0.31;
    double kf2 = 0.17;
    bond->addPerBondParameter("kf");
    std::vector<double> kf1v = {kf1};
    bond->addBond(0, 1, kf1v);
    std::vector<double> kf2v = {kf2};
    bond->addBond(0, 2, kf2v);

    ATMForce* atm = new ATMForce(lmbd, lmbd, 0., 0, 0, umax, ubcore, acore, direction);
    //swap particles 1 and 2
    atm->addParticle( ); //particle 0 is not displaced
    atm->addParticle(new ATMForce::ParticleOffsetDisplacement(2,  1) );
    atm->addParticle(new ATMForce::ParticleOffsetDisplacement(1,  2) );

    atm->addForce(bond);
    system.addForce(atm);

    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    for (double lm : {0.0, 0.5, 1.0}) {
        context.setParameter(ATMForce::Lambda1(), lm);
        context.setParameter(ATMForce::Lambda2(), lm);
        State state = context.getState(State::Energy | State::Forces );
        double epot = state.getPotentialEnergy();
        double u0, u1, energy;
        atm->getPerturbationEnergy(context, u1, u0, energy);
        double epert = u1 - u0;
        ASSERT_EQUAL_TOL(energy, epot, 1e-6);
        ASSERT_EQUAL_TOL(0.5*kf1*r1sq + 0.5*kf2*r2sq, u0, 1e-6);
        ASSERT_EQUAL_TOL(0.5*kf1*r2sq + 0.5*kf2*r1sq, u1, 1e-6);
        ASSERT_EQUAL_TOL(0.5*kf1*(r2sq-r1sq) + 0.5*kf2*(r1sq-r2sq), epert, 1e-6);
        ASSERT_EQUAL_TOL(u0 + lm*epert, epot, 1e-6);
        ASSERT_EQUAL_VEC(- ( ((1.-lm)*kf1+lm*kf2)*r1 ), state.getForces()[1], 1e-6);
        ASSERT_EQUAL_VEC(- ( ((1.-lm)*kf2+lm*kf1)*r2 ), state.getForces()[2], 1e-6);
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
    atm->addParticle();
    //second particle is displaced at both states but by the same amount (1,0,0)
    Vec3 displ0 = Vec3(1., 0., 0.);
    atm->addParticle(new ATMForce::FixedDisplacement(displ0, displ0));
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
    Vec3   displ = Vec3(5., 0., 0.);
    atm->addParticle();
    atm->addParticle(new ATMForce::FixedDisplacement(displ));
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
                atm->addParticle();
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

    //in this second scenario the non-bonded force is removed from the System
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
                atm->addParticle();
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
    Vec3   displ = Vec3(5., 0., 0.);
    atm->addParticle();
    atm->addParticle(new ATMForce::FixedDisplacement(displ));
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

    atm->addParticle();
    atm->addParticle(new ATMForce::FixedDisplacement(displ));
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
        atm->addParticle(new ATMForce::FixedDisplacement(d));
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

void testLargeSystemSwap() {
    // Create a system with lots of particles in an external field
    // that depends on atom indexes. Swap their positions, check
    // energies and forces.

    int numParticles = 1000;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    CustomExternalForce* external = new CustomExternalForce("qf*(x^2 + 2*y^2 + 3*z^2)");
    external->addPerParticleParameter("qf");
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    atm->addForce(external);
    system.addForce(atm);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    Vec3 nodispl = Vec3(0,0,0);
    vector<Vec3> positions;
    vector<int> target_particle(numParticles);
    for (int i = 0; i < numParticles; i++) {
        target_particle[i] = i;
    }
    auto rng = default_random_engine {};
    shuffle(begin(target_particle), end(target_particle), rng);
    vector<int> target_particle_inv(numParticles);
    for (int i = 0; i < numParticles; i++) {
        target_particle_inv[target_particle[i]] = i;
    }
    vector<double> qf(numParticles);
    for (int i = 0; i < numParticles; i++)
      qf[i] = (double)i/(double)numParticles;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions.push_back(3*Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
        external->addParticle(i, {qf[i]});
        atm->addParticle(new ATMForce::ParticleOffsetDisplacement(target_particle[i], i));
    }

    double energy0 = 0.;
    for (int i = 0; i < numParticles; i++) {
        Vec3 pos = positions[i];
        energy0 += qf[i]*(pos[0]*pos[0]+2*pos[1]*pos[1]+3*pos[2]*pos[2]);
    }
    double energy1 = 0.;
    for (int i = 0; i < numParticles; i++) {
        Vec3 pos = positions[target_particle[i]];
        energy1 += qf[i]*(pos[0]*pos[0]+2*pos[1]*pos[1]+3*pos[2]*pos[2]);
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

    // Evaluate energies and forces at lambda 0 and 1
    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    for (double lambda : {0.0, 1.0}) {
        context.setParameter(ATMForce::Lambda1(), lambda);
        context.setParameter(ATMForce::Lambda2(), lambda);
        State state = context.getState(State::Energy | State::Forces);
        double u1, u0, energy;
        double epot = state.getPotentialEnergy();
        atm->getPerturbationEnergy(context, u1, u0, energy);
        ASSERT_EQUAL_TOL(u0, energy0, 1e-6);
        ASSERT_EQUAL_TOL(u1, energy1, 1e-6);
        ASSERT_EQUAL_TOL(u0+lambda*(u1-u0), epot, 1e-6);
        for (int i = 0; i < numParticles; i++) {
            int l;
            if (lambda > 0){
                l = target_particle_inv[i];
            }else{
                l = i;
            }
            Vec3 pos = positions[i];
            Vec3 expectedForce(-2*pos[0], -4*pos[1], -6*pos[2]);
            ASSERT_EQUAL_VEC(qf[l]*expectedForce, state.getForces()[i], 1e-6);
        }
    }
}

void testChangingBoxVectors() {
    // Create a periodic system with incorrect default box vectors.

    int numParticles = 500;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    NonbondedForce* force = new NonbondedForce();
    force->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    ATMForce* atm = new ATMForce(0.0, 0.0, 0.1, 0.0, 0.0, 1e6, 5e5, 1.0/16, 1.0);
    atm->addForce(force);
    system.addForce(atm);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions.push_back(3*Vec3(genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5, genrand_real2(sfmt)-0.5));
        force->addParticle(0.0, 0.1, 1.0);
        atm->addParticle();
        for (int j = 0; j < i; j++) {
            Vec3 delta = positions[i]-positions[j];
            for (int k = 0; k < 3; k++)
                delta[k] -= round(delta[k]/2.0)*2.0;
            if (sqrt(delta.dot(delta)) < 0.1)
                force->addException(i, j, 0.0, 0.1, 0.0);
        }
    }

    // Set the correct box vectors in the context and check that energy is calculated correctly.

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 2, 0), Vec3(0, 0, 2));
    double energy1 = context.getState(State::Energy).getPotentialEnergy();
    double energy2 = context.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(energy1, energy2, 1e-6);
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
                atm->addParticle();
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
        testAPI();
        test2Particles();
        test3ParticlesSwap();
        test2Particles2Displacement0();
        test2ParticlesSoftCore();
        testNonbonded();
        testNonbondedwithEndpointClash();
        testParticlesCustomExpressionLinear();
        testParticlesCustomExpressionSoftplus();
        testLargeSystem();
        testLargeSystemSwap();
        testChangingBoxVectors();
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

