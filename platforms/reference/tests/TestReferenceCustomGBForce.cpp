
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2013 Stanford University and the Authors.      *
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
 * This tests all the different force terms in the reference implementation of CustomGBForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomGBForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/GBVIForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testOBC(GBSAOBCForce::NonbondedMethod obcMethod, CustomGBForce::NonbondedMethod customMethod) {
    const int numMolecules = 70;
    const int numParticles = numMolecules*2;
    const double boxSize = 10.0;
    const double cutoff = 2.0;
    ReferencePlatform platform;

    // Create two systems: one with a GBSAOBCForce, and one using a CustomGBForce to implement the same interaction.

    System standardSystem;
    System customSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0.0, 0.0), Vec3(0.0, boxSize, 0.0), Vec3(0.0, 0.0, boxSize));
    customSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0.0, 0.0), Vec3(0.0, boxSize, 0.0), Vec3(0.0, 0.0, boxSize));
    GBSAOBCForce* obc = new GBSAOBCForce();
    CustomGBForce* custom = new CustomGBForce();
    obc->setCutoffDistance(cutoff);
    custom->setCutoffDistance(cutoff);
    custom->addPerParticleParameter("q");
    custom->addPerParticleParameter("radius");
    custom->addPerParticleParameter("scale");
    custom->addGlobalParameter("solventDielectric", obc->getSolventDielectric());
    custom->addGlobalParameter("soluteDielectric", obc->getSoluteDielectric());
    custom->addComputedValue("I", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                                  "U=r+sr2;"
                                  "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                                  "L=max(or1, D);"
                                  "D=abs(r-sr2);"
                                  "sr2 = scale2*or2;"
                                  "or1 = radius1-0.009; or2 = radius2-0.009", CustomGBForce::ParticlePairNoExclusions);
    custom->addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                                  "psi=I*or; or=radius-0.009", CustomGBForce::SingleParticle);
    custom->addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B", CustomGBForce::SingleParticle);
    string invCutoffString = "";
    if (obcMethod != GBSAOBCForce::NoCutoff) {
        stringstream s;
        s<<(1.0/cutoff);
        invCutoffString = s.str();
    }
    custom->addEnergyTerm("138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2*("+invCutoffString+"-1/f);"
                          "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePairNoExclusions);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<double> params(3);
    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {
            obc->addParticle(1.0, 0.2, 0.5);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.5;
            custom->addParticle(params);
            obc->addParticle(-1.0, 0.1, 0.5);
            params[0] = -1.0;
            params[1] = 0.1;
            custom->addParticle(params);
        }
        else {
            obc->addParticle(1.0, 0.2, 0.8);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.8;
            custom->addParticle(params);
            obc->addParticle(-1.0, 0.1, 0.8);
            params[0] = -1.0;
            params[1] = 0.1;
            custom->addParticle(params);
        }
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
    }
    obc->setNonbondedMethod(obcMethod);
    custom->setNonbondedMethod(customMethod);
    standardSystem.addForce(obc);
    customSystem.addForce(custom);
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context context1(standardSystem, integrator1, platform);
    context1.setPositions(positions);
    context1.setVelocities(velocities);
    State state1 = context1.getState(State::Forces | State::Energy);
    Context context2(customSystem, integrator2, platform);
    context2.setPositions(positions);
    context2.setVelocities(velocities);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
    }
    
    // Try changing the particle parameters and make sure it's still correct.
    
    for (int i = 0; i < numMolecules/2; i++) {
        obc->setParticleParameters(2*i, 1.1, 0.3, 0.6);
        params[0] = 1.1;
        params[1] = 0.3;
        params[2] = 0.6;
        custom->setParticleParameters(2*i, params);
        obc->setParticleParameters(2*i+1, -1.1, 0.2, 0.4);
        params[0] = -1.1;
        params[1] = 0.2;
        params[2] = 0.4;
        custom->setParticleParameters(2*i+1, params);
    }
    obc->updateParametersInContext(context1);
    custom->updateParametersInContext(context2);
    state1 = context1.getState(State::Forces | State::Energy);
    state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
    }
}

void testMembrane() {
    const int numMolecules = 70;
    const int numParticles = numMolecules*2;
    const double boxSize = 10.0;
    ReferencePlatform platform;

    // Create a system with an implicit membrane.

    System system;
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
    }
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0.0, 0.0), Vec3(0.0, boxSize, 0.0), Vec3(0.0, 0.0, boxSize));
    CustomGBForce* custom = new CustomGBForce();
    custom->setCutoffDistance(2.0);
    custom->addPerParticleParameter("q");
    custom->addPerParticleParameter("radius");
    custom->addPerParticleParameter("scale");
    custom->addGlobalParameter("thickness", 3);
    custom->addGlobalParameter("solventDielectric", 78.3);
    custom->addGlobalParameter("soluteDielectric", 1);
    custom->addComputedValue("Imol", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                             "U=r+sr2;"
                             "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                             "L=max(or1, D);"
                             "D=abs(r-sr2);"
                             "sr2 = scale2*or2;"
                             "or1 = radius1-0.009; or2 = radius2-0.009", CustomGBForce::ParticlePairNoExclusions);
    custom->addComputedValue("Imem", "(1/radius+2*log(2)/thickness)/(1+exp(7.2*(abs(z)+radius-0.5*thickness)))", CustomGBForce::SingleParticle);
    custom->addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                             "psi=max(Imol,Imem)*or; or=radius-0.009", CustomGBForce::SingleParticle);
    custom->addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B", CustomGBForce::SingleParticle);
    custom->addEnergyTerm("-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                          "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePairNoExclusions);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<double> params(3);
    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.5;
            custom->addParticle(params);
            params[0] = -1.0;
            params[1] = 0.1;
            custom->addParticle(params);
        }
        else {
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.8;
            custom->addParticle(params);
            params[0] = -1.0;
            params[1] = 0.1;
            custom->addParticle(params);
        }
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
    }
    system.addForce(custom);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocities(velocities);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    double norm = 0.0;
    for (int i = 0; i < (int) forces.size(); ++i)
        norm += forces[i].dot(forces[i]);
    norm = std::sqrt(norm);
    const double stepSize = 1e-2;
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
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-3);
}

void testTabulatedFunction() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomGBForce* force = new CustomGBForce();
    force->addComputedValue("a", "0", CustomGBForce::ParticlePair);
    force->addEnergyTerm("fn(r)+1", CustomGBForce::ParticlePair);
    force->addParticle(vector<double>());
    force->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < 21; i++)
        table.push_back(std::sin(0.25*i));
    force->addTabulatedFunction("fn", new Continuous1DFunction(table, 1.0, 6.0));
    system.addForce(force);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 1; i < 30; i++) {
        double x = (7.0/30.0)*i;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double force = (x < 1.0 || x > 6.0 ? 0.0 : -std::cos(x-1.0));
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : std::sin(x-1.0))+1.0;
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], 0.1);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], 0.1);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.02);
    }
}

void testMultipleChainRules() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomGBForce* force = new CustomGBForce();
    force->addComputedValue("a", "2*r", CustomGBForce::ParticlePair);
    force->addComputedValue("b", "a+1", CustomGBForce::SingleParticle);
    force->addComputedValue("c", "2*b+a", CustomGBForce::SingleParticle);
    force->addEnergyTerm("0.1*a+1*b+10*c", CustomGBForce::SingleParticle); // 0.1*(2*r) + 2*r+1 + 10*(3*a+2) = 0.2*r + 2*r+1 + 40*r+20+20*r = 62.2*r+21
    force->addParticle(vector<double>());
    force->addParticle(vector<double>());
    system.addForce(force);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 1; i < 5; i++) {
        positions[1] = Vec3(i, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(124.4, 0, 0), forces[0], 1e-4);
        ASSERT_EQUAL_VEC(Vec3(-124.4, 0, 0), forces[1], 1e-4);
        ASSERT_EQUAL_TOL(2*(62.2*i+21), state.getPotentialEnergy(), 0.02);
    }
}

void testPositionDependence() {
    ReferencePlatform platform;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomGBForce* force = new CustomGBForce();
    force->addComputedValue("a", "r", CustomGBForce::ParticlePair);
    force->addComputedValue("b", "a+x*y", CustomGBForce::SingleParticle);
    force->addEnergyTerm("b*z", CustomGBForce::SingleParticle);
    force->addEnergyTerm("b1+b2", CustomGBForce::ParticlePair); // = 2*r+x1*y1+x2*y2
    force->addParticle(vector<double>());
    force->addParticle(vector<double>());
    system.addForce(force);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    vector<Vec3> forces(2);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < 5; i++) {
        positions[0] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        positions[1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        Vec3 delta = positions[0]-positions[1];
        double r = sqrt(delta.dot(delta));
        double energy = 2*r+positions[0][0]*positions[0][1]+positions[1][0]*positions[1][1];
        for (int j = 0; j < 2; j++)
            energy += positions[j][2]*(r+positions[j][0]*positions[j][1]);
        Vec3 force1(-(1+positions[0][2])*delta[0]/r-(1+positions[0][2])*positions[0][1]-(1+positions[1][2])*delta[0]/r,
                    -(1+positions[0][2])*delta[1]/r-(1+positions[0][2])*positions[0][0]-(1+positions[1][2])*delta[1]/r,
                    -(1+positions[0][2])*delta[2]/r-(r+positions[0][0]*positions[0][1])-(1+positions[1][2])*delta[2]/r);
        Vec3 force2((1+positions[0][2])*delta[0]/r+(1+positions[1][2])*delta[0]/r-(1+positions[1][2])*positions[1][1],
                    (1+positions[0][2])*delta[1]/r+(1+positions[1][2])*delta[1]/r-(1+positions[1][2])*positions[1][0],
                    (1+positions[0][2])*delta[2]/r+(1+positions[1][2])*delta[2]/r-(r+positions[1][0]*positions[1][1]));
        ASSERT_EQUAL_VEC(force1, forces[0], 1e-4);
        ASSERT_EQUAL_VEC(force2, forces[1], 1e-4);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.02);

        // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

        double norm = 0.0;
        for (int i = 0; i < (int) forces.size(); ++i)
            norm += forces[i].dot(forces[i]);
        norm = std::sqrt(norm);
        const double stepSize = 1e-3;
        double step = 0.5*stepSize/norm;
        vector<Vec3> positions2(2), positions3(2);
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
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-3);
    }
}

void testExclusions() {
    ReferencePlatform platform;
    for (int i = 3; i < 4; i++) {
        System system;
        system.addParticle(1.0);
        system.addParticle(1.0);
        VerletIntegrator integrator(0.01);
        CustomGBForce* force = new CustomGBForce();
        force->addComputedValue("a", "r", i < 2 ? CustomGBForce::ParticlePair : CustomGBForce::ParticlePairNoExclusions);
        force->addEnergyTerm("a", CustomGBForce::SingleParticle);
        force->addEnergyTerm("(1+a1+a2)*r", i%2 == 0 ? CustomGBForce::ParticlePair : CustomGBForce::ParticlePairNoExclusions);
        force->addParticle(vector<double>());
        force->addParticle(vector<double>());
        force->addExclusion(0, 1);
        system.addForce(force);
        Context context(system, integrator, platform);
        vector<Vec3> positions(2);
        positions[0] = Vec3(0, 0, 0);
        positions[1] = Vec3(1, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double f, energy;
        switch (i)
        {
            case 0: // e = 0
                f = 0;
                energy = 0;
                break;
            case 1: // e = r
                f = 1;
                energy = 1;
                break;
            case 2: // e = 2r
                f = 2;
                energy = 2;
                break;
            case 3: // e = 3r + 2r^2
                f = 7;
                energy = 5;
                break;
            default:
                ASSERT(false);
        }
        ASSERT_EQUAL_VEC(Vec3(f, 0, 0), forces[0], 1e-4);
        ASSERT_EQUAL_VEC(Vec3(-f, 0, 0), forces[1], 1e-4);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-4);

        // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

        double norm = 0.0;
        for (int i = 0; i < (int) forces.size(); ++i)
            norm += forces[i].dot(forces[i]);
        norm = std::sqrt(norm);
        const double stepSize = 1e-3;
        double step = stepSize/norm;
        for (int i = 0; i < (int) positions.size(); ++i) {
            Vec3 p = positions[i];
            Vec3 f = forces[i];
            positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        }
        context.setPositions(positions);
        State state2 = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/stepSize, 1e-3*abs(state.getPotentialEnergy()));
    }
}

// create custom GB/VI force

static CustomGBForce* createCustomGBVI( double solventDielectric, double soluteDielectric ) {

    CustomGBForce* customGbviForce  = new CustomGBForce();

    customGbviForce->setCutoffDistance(2.0);

    customGbviForce->addPerParticleParameter("q");
    customGbviForce->addPerParticleParameter("radius");
    customGbviForce->addPerParticleParameter("scaleFactor"); // derived in GBVIForce implmentation, but parameter here
    customGbviForce->addPerParticleParameter("gamma");

    customGbviForce->addGlobalParameter("solventDielectric", solventDielectric);
    customGbviForce->addGlobalParameter("soluteDielectric", soluteDielectric);

    customGbviForce->addComputedValue("V", "                uL - lL + factor3/(radius1*radius1*radius1);"
                                      "uL                   = 1.5*x2uI*(0.25*rI-0.33333*xuI+0.125*(r2-S2)*rI*x2uI);"
                                      "lL                   = 1.5*x2lI*(0.25*rI-0.33333*xlI+0.125*(r2-S2)*rI*x2lI);"
                                      "x2lI                 = 1.0/(xl*xl);"
                                      "xlI                  = 1.0/(xl);"
                                      "xuI                  = 1.0/(xu);"
                                      "x2uI                 = 1.0/(xu*xu);"
                                      "xu                   = (r+scaleFactor2);"
                                      "rI                   = 1.0/(r);"
                                      "r2                   = (r*r);"
                                      "xl                   = factor1*lMax + factor2*xuu + factor3*(r-scaleFactor2);"
                                      "xuu                  = (r+scaleFactor2);"
                                      "S2                   = (scaleFactor2*scaleFactor2);"
                                      "factor1              = step(r-absRadiusScaleDiff);"
                                      "absRadiusScaleDiff   = abs(radiusScaleDiff);"
                                      "radiusScaleDiff      = (radius1-scaleFactor2);"
                                      "factor2              = step(radius1-scaleFactor2-r);"
                                      "factor3              = step(scaleFactor2-radius1-r);"
                                      "lMax                 = max(radius1,r-scaleFactor2);"
                                      , CustomGBForce::ParticlePairNoExclusions);

    customGbviForce->addComputedValue("B", "(1.0/(radius*radius*radius)-V)^(-0.33333333)", CustomGBForce::SingleParticle);

    // nonpolar term + polar self energy

    customGbviForce->addEnergyTerm("(-138.935485*0.5*((1.0/soluteDielectric)-(1.0/solventDielectric))*q^2/B)-((1.0/soluteDielectric)-(1.0/solventDielectric))*((gamma*(radius/B)^3))", CustomGBForce::SingleParticle);

    // polar pair energy

    customGbviForce->addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                                   "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePairNoExclusions);

    return customGbviForce;
}

// ethance GB/VI test case

static void buildEthane( GBVIForce* gbviForce, std::vector<Vec3>& positions ) {

    const int numParticles = 8;

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;
    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
       C_gamma  = -0.2863;
       H_radius =  0.125;
       H_gamma  =  0.2437;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    for( int i = 0; i < numParticles; i++ ){
       gbviForce->addParticle( H_charge, H_radius, H_gamma);
    }
    gbviForce->setParticleParameters( 1, C_charge, C_radius, C_gamma);
    gbviForce->setParticleParameters( 4, C_charge, C_radius, C_gamma);
 
    gbviForce->addBond( 0, 1, C_HBondDistance );
    gbviForce->addBond( 2, 1, C_HBondDistance );
    gbviForce->addBond( 3, 1, C_HBondDistance );
    gbviForce->addBond( 1, 4, C_CBondDistance );
    gbviForce->addBond( 5, 4, C_HBondDistance );
    gbviForce->addBond( 6, 4, C_HBondDistance );
    gbviForce->addBond( 7, 4, C_HBondDistance );
    
    std::vector<pair<int, int> > bondExceptions;
    std::vector<double> bondDistances;
    
    bondExceptions.push_back(pair<int, int>(0, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(2, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(3, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(1, 4)); 
    bondDistances.push_back( C_CBondDistance );
    
    bondExceptions.push_back(pair<int, int>(5, 4)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(6, 4)); 
    bondDistances.push_back( C_HBondDistance );
 
    bondExceptions.push_back(pair<int, int>(7, 4));
    bondDistances.push_back( C_HBondDistance );
 
    positions.resize(numParticles);
    positions[0] = Vec3(0.5480,    1.7661,    0.0000);
    positions[1] = Vec3(0.7286,    0.8978,    0.6468);
    positions[2] = Vec3(0.4974,    0.0000,    0.0588);
    positions[3] = Vec3(0.0000,    0.9459,    1.4666);
    positions[4] = Vec3(2.1421,    0.8746,    1.1615);
    positions[5] = Vec3(2.3239,    0.0050,    1.8065);
    positions[6] = Vec3(2.8705,    0.8295,    0.3416);
    positions[7] = Vec3(2.3722,    1.7711,    1.7518);

}

// dimer GB/VI test case

static void buildDimer( GBVIForce* gbviForce, std::vector<Vec3>& positions ) {

    const int numParticles = 2;

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;
    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;

    H_charge    = 0.0;
    C_charge    = 0.0;
    if( AM1_BCC ){
       C_radius =  0.180;
       C_gamma  = -0.2863;
       H_radius =  0.125;
       H_gamma  =  0.2437;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    for( int i = 0; i < numParticles; i++ ){
       gbviForce->addParticle( H_charge, H_radius, H_gamma);
    }
    gbviForce->setParticleParameters( 1, C_charge, C_radius, C_gamma);
 
    gbviForce->addBond( 0, 1, C_HBondDistance );
    std::vector<pair<int, int> > bondExceptions;
    std::vector<double> bondDistances;
    
    bondExceptions.push_back(pair<int, int>(0, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    positions.resize(numParticles);
    positions[0] = Vec3(0.0,       0.0,       0.0);
    positions[1] = Vec3(0.15,       0.0,       0.0);
}

// monomer GB/VI test case

static void buildMonomer( GBVIForce* gbviForce, std::vector<Vec3>& positions ) {

    const int numParticles = 1;

    double H_radius, H_gamma, H_charge;

    H_charge = 1.0;
    H_radius = 0.125;
    H_gamma  = 0.2437;

    for( int i = 0; i < numParticles; i++ ){
       gbviForce->addParticle( H_charge, H_radius, H_gamma);
    }
    positions.resize(numParticles);
    positions[0] = Vec3(0.0,    0.0,    0.0);
}

// taken from gbviForceImpl class
// computes the scaled radii based on covalent info and atomic radii

static void findScaledRadii( GBVIForce& gbviForce, std::vector<double> & scaledRadii) {

    int     numberOfParticles = gbviForce.getNumParticles();
    int numberOfBonds         = gbviForce.getNumBonds();
    
    // load 1-2 atom pairs along w/ bond distance using HarmonicBondForce & constraints
    // numberOfBonds < 1, indicating they were not set by the user
    
    if( numberOfBonds < 1 && numberOfParticles > 1 ){
        (void) fprintf( stderr, "Warning: no covalent bonds set for GB/VI force!\n" );
    }
    
    std::vector< std::vector<int> > bondIndices;
    bondIndices.resize( numberOfBonds );
    
    std::vector<double> bondLengths;
    bondLengths.resize( numberOfBonds );

    scaledRadii.resize(numberOfParticles);
    for (int i = 0; i < numberOfParticles; i++) {
        double charge, radius, gamma;
        gbviForce.getParticleParameters(i, charge, radius, gamma);
        scaledRadii[i] = radius;
    }

    for (int i = 0; i < numberOfBonds; i++) {
        int particle1, particle2;
        double bondLength;
        gbviForce.getBondParameters(i, particle1, particle2, bondLength);
        if (particle1 < 0 || particle1 >= gbviForce.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= gbviForce.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (bondLength < 0 ) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: negative bondlength: ";
            msg << bondLength;
            throw OpenMMException(msg.str());
        }
        bondIndices[i].push_back( particle1 );
        bondIndices[i].push_back( particle2 );
        bondLengths[i] = bondLength;
    }


    // load 1-2 indicies for each atom 

    std::vector<std::vector<int> > bonded12(numberOfParticles);

    for (int i = 0; i < (int) bondIndices.size(); ++i) {
        bonded12[bondIndices[i][0]].push_back(i);
        bonded12[bondIndices[i][1]].push_back(i);
    }

    int errors = 0;

    // compute scaled radii (Eq. 5 of Labute paper [JCC 29 p. 1693-1698 2008])

    for (int j = 0; j < (int) bonded12.size(); ++j){

        double charge;
        double gamma;
        double radiusJ;
        double scaledRadiusJ;
     
        gbviForce.getParticleParameters(j, charge, radiusJ, gamma); 

        if(  bonded12[j].size() == 0 ){
            if( numberOfParticles > 1 ){
                (void) fprintf( stderr, "Warning GBVIForceImpl::findScaledRadii atom %d has no covalent bonds; using atomic radius=%.3f.\n", j, radiusJ );
            }
            scaledRadiusJ = radiusJ;
//             errors++;
        } else {

            double rJ2    = radiusJ*radiusJ;
    
            // loop over bonded neighbors of atom j, applying Eq. 5 in Labute

            scaledRadiusJ = 0.0;
            for (int i = 0; i < (int) bonded12[j].size(); ++i){
    
               int index            = bonded12[j][i];
               int bondedAtomIndex  = (j == bondIndices[index][0]) ? bondIndices[index][1] : bondIndices[index][0];
              
               double radiusI;
               gbviForce.getParticleParameters(bondedAtomIndex, charge, radiusI, gamma); 
               double rI2           = radiusI*radiusI;
    
               double a_ij          = (radiusI - bondLengths[index]);
                      a_ij         *= a_ij;
                      a_ij          = (rJ2 - a_ij)/(2.0*bondLengths[index]);
    
               double a_ji          = radiusJ - bondLengths[index];
                      a_ji         *= a_ji;
                      a_ji          = (rI2 - a_ji)/(2.0*bondLengths[index]);
    
               scaledRadiusJ       += a_ij*a_ij*(3.0*radiusI - a_ij) + a_ji*a_ji*( 3.0*radiusJ - a_ji );
            }
    
            scaledRadiusJ  = (radiusJ*radiusJ*radiusJ) - 0.125*scaledRadiusJ; 
            if( scaledRadiusJ > 0.0 ){
                scaledRadiusJ  = 0.95*pow( scaledRadiusJ, (1.0/3.0) );
            } else {
                scaledRadiusJ  = 0.0;
            }
        }
        //(void) fprintf( stderr, "scaledRadii %d %12.4f\n", j, scaledRadiusJ );
        scaledRadii[j] = scaledRadiusJ;

    }

    // abort if errors

    if( errors ){
        throw OpenMMException("GBVIForceImpl::findScaledRadii errors -- aborting");
    }

#if GBVIDebug
    (void) fprintf( stderr, "                  R              q          gamma   scaled radii no. bnds\n" );
    double totalQ = 0.0;
    for( int i = 0; i < (int) scaledRadii.size(); i++ ){

        double charge;
        double gamma;
        double radiusI;
     
        gbviForce.getParticleParameters(i, charge, radiusI, gamma); 
        totalQ += charge;
        (void) fprintf( stderr, "%4d %14.5e %14.5e %14.5e %14.5e %d\n", i, radiusI, charge, gamma, scaledRadii[i], (int) bonded12[i].size() );
    }
    (void) fprintf( stderr, "Total charge=%e\n", totalQ );
    (void) fflush( stderr );
#endif

#undef GBVIDebug

}

// load parameters from gbviForce to customGbviForce
// findScaledRadii() is called to calculate the scaled radii (S)
// S is derived quantity in GBVIForce, not a parameter is the case here

static void loadGbviParameters( GBVIForce* gbviForce, CustomGBForce* customGbviForce ) {

    int numParticles = gbviForce->getNumParticles();

    // charge, radius, scale factor, gamma

    vector<double> params(4);
    std::vector<double> scaledRadii;
    findScaledRadii( *gbviForce, scaledRadii);

    for( int ii = 0; ii < numParticles; ii++) {
        double charge, radius, gamma;
        gbviForce->getParticleParameters( ii, charge, radius, gamma );
        params[0] = charge;
        params[1] = radius;
        params[2] = scaledRadii[ii];
        params[3] = gamma;
        customGbviForce->addParticle(params);
    }

}

void testGBVI(GBVIForce::NonbondedMethod gbviMethod, CustomGBForce::NonbondedMethod customGbviMethod, std::string molecule) {

    const int numMolecules = 1;
    const double boxSize   = 10.0;
    ReferencePlatform platform;

    GBVIForce*        gbvi = new GBVIForce();
    std::vector<Vec3> positions;

    // select molecule

    if( molecule == "Monomer" ){
        buildMonomer( gbvi, positions );
    } else if( molecule == "Dimer" ){
        buildDimer( gbvi, positions );
    } else {
        buildEthane( gbvi, positions );
    }

    int numParticles = gbvi->getNumParticles();
    System standardSystem;
    System customGbviSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customGbviSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0.0, 0.0), Vec3(0.0, boxSize, 0.0), Vec3(0.0, 0.0, boxSize));
    customGbviSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0.0, 0.0), Vec3(0.0, boxSize, 0.0), Vec3(0.0, 0.0, boxSize));
    gbvi->setCutoffDistance(2.0);

    // create customGbviForce GBVI force

    CustomGBForce* customGbviForce  = createCustomGBVI( gbvi->getSolventDielectric(), gbvi->getSoluteDielectric() );
    customGbviForce->setCutoffDistance(2.0);

    // load parameters from gbvi to customGbviForce

    loadGbviParameters( gbvi, customGbviForce );

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<Vec3> velocities(numParticles);
    for (int ii = 0; ii < numParticles; ii++) {
        velocities[ii] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
    }
    gbvi->setNonbondedMethod(gbviMethod);
    customGbviForce->setNonbondedMethod(customGbviMethod);

    standardSystem.addForce(gbvi);
    customGbviSystem.addForce(customGbviForce);

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);

    Context context1(standardSystem, integrator1, platform);
    context1.setPositions(positions);
    context1.setVelocities(velocities);
    State state1 = context1.getState(State::Forces | State::Energy);

    Context context2(customGbviSystem, integrator2, platform);
    context2.setPositions(positions);
    context2.setVelocities(velocities);
    State state2 = context2.getState(State::Forces | State::Energy);

    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);

    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
    }
}

int main() {

    try {
        testOBC(GBSAOBCForce::NoCutoff, CustomGBForce::NoCutoff);
        testOBC(GBSAOBCForce::CutoffNonPeriodic, CustomGBForce::CutoffNonPeriodic);
        testOBC(GBSAOBCForce::CutoffPeriodic, CustomGBForce::CutoffPeriodic);
        testMembrane();
        testTabulatedFunction();
        testMultipleChainRules();
        testPositionDependence();
        testExclusions();

        // GBVI tests

        testGBVI(GBVIForce::NoCutoff, CustomGBForce::NoCutoff, "Monomer");
        testGBVI(GBVIForce::NoCutoff, CustomGBForce::NoCutoff, "Dimer");
        testGBVI(GBVIForce::NoCutoff, CustomGBForce::NoCutoff, "Ethane");

    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

