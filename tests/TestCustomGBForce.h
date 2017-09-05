
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "openmm/CustomGBForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testOBC(GBSAOBCForce::NonbondedMethod obcMethod, CustomGBForce::NonbondedMethod customMethod) {
    const int numMolecules = 70;
    const int numParticles = numMolecules*2;
    const double boxSize = 10.0;
    const double cutoff = 2.0;

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
        custom->addExclusion(2*i, 2*i+1);
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
    }
    obc->setNonbondedMethod(obcMethod);
    custom->setNonbondedMethod(customMethod);
    standardSystem.addForce(obc);
    customSystem.addForce(custom);
    if (customMethod == CustomGBForce::CutoffPeriodic) {
        ASSERT(custom->usesPeriodicBoundaryConditions());
        ASSERT(obc->usesPeriodicBoundaryConditions());
        ASSERT(standardSystem.usesPeriodicBoundaryConditions());
        ASSERT(customSystem.usesPeriodicBoundaryConditions());
    }
    else {
        ASSERT(!custom->usesPeriodicBoundaryConditions());
        ASSERT(!obc->usesPeriodicBoundaryConditions());
        ASSERT(!standardSystem.usesPeriodicBoundaryConditions());
        ASSERT(!customSystem.usesPeriodicBoundaryConditions());
    }
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
        // Check interpolated values.
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
    for (int i = 1; i < 20; i++) {
        // These values are exactly on table points, so they should match more precisely.
        double x = 0.25*i+1.0;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Energy);
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : std::sin(x-1.0))+1.0;
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-4);
    }
}

void testMultipleChainRules() {
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
    for (int i = 0; i < 4; i++) {
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
        if (norm > 0) {
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
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    CustomGBForce* force = new CustomGBForce();
    force->addComputedValue("a", "r+none", CustomGBForce::ParticlePair);
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

void testEnergyParameterDerivatives() {
    // Create a box of particles.
    
    const int numParticles = 40;
    const int numParameters = 4;
    const double boxSize = 2.0;
    const double delta = 1e-3;
    const string paramNames[] = {"A", "B", "C", "D"};
    const double paramValues[] = {0.8, 2.1, 3.2, 1.3};
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    CustomGBForce* force = new CustomGBForce();
    system.addForce(force);
    force->addComputedValue("a", "0.5*(r-A)^2", CustomGBForce::ParticlePair);
    force->addComputedValue("b", "a+B", CustomGBForce::SingleParticle);
    force->addEnergyTerm("C*(a1+b1+a2+b2+r)^0.8", CustomGBForce::ParticlePair);
    force->addEnergyTerm("(D-B)*b", CustomGBForce::SingleParticle);
    for (int i = 0; i < numParameters; i++)
        force->addGlobalParameter(paramNames[i], paramValues[i]);
    for (int i = numParameters-1; i >= 0; i--)
        force->addEnergyParameterDerivative(paramNames[i]);
    force->setNonbondedMethod(CustomGBForce::CutoffPeriodic);
    force->setCutoffDistance(1.0);
    vector<Vec3> positions;
    vector<double> parameters;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(parameters);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*boxSize);
    }
    
    // Compute the energy derivative and compare it to a finite difference approximation.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    map<string, double> derivs = context.getState(State::ParameterDerivatives).getEnergyParameterDerivatives();
    for (int i = 0; i < numParameters; i++) {
        context.setParameter(paramNames[i], paramValues[i]+delta);
        double energy1 = context.getState(State::Energy).getPotentialEnergy();
        context.setParameter(paramNames[i], paramValues[i]-delta);
        double energy2 = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL((energy1-energy2)/(2*delta), derivs[paramNames[i]], 5e-3);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testOBC(GBSAOBCForce::NoCutoff, CustomGBForce::NoCutoff);
        testOBC(GBSAOBCForce::CutoffNonPeriodic, CustomGBForce::CutoffNonPeriodic);
        testOBC(GBSAOBCForce::CutoffPeriodic, CustomGBForce::CutoffPeriodic);
        testMembrane();
        testTabulatedFunction();
        testMultipleChainRules();
        testPositionDependence();
        testExclusions();
        testIllegalVariable();
        testEnergyParameterDerivatives();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

