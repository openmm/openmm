/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2014 Stanford University and the Authors.      *
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
 * This tests the CUDA implementation of virtual sites.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VirtualSite.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

CudaPlatform platform;

/**
 * Check that massless particles are handled correctly.
 */
void testMasslessParticle() {
    System system;
    system.addParticle(0.0);
    system.addParticle(1.0);
    CustomBondForce* bonds = new CustomBondForce("-1/r");
    system.addForce(bonds);
    vector<double> params;
    bonds->addBond(0, 1, params);
    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    vector<Vec3> velocities(2);
    velocities[0] = Vec3(0, 0, 0);
    velocities[1] = Vec3(0, 1, 0);
    context.setVelocities(velocities);
    
    // The second particle should move in a circular orbit around the first one.
    // Compare it to the analytical solution.
    
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Positions | State::Velocities | State::Forces);
        double time = state.getTime();
        ASSERT_EQUAL_VEC(Vec3(), state.getPositions()[0], 0.0);
        ASSERT_EQUAL_VEC(Vec3(), state.getVelocities()[0], 0.0);
        ASSERT_EQUAL_VEC(Vec3(cos(time), sin(time), 0), state.getPositions()[1], 0.01);
        ASSERT_EQUAL_VEC(Vec3(-sin(time), cos(time), 0), state.getVelocities()[1], 0.01);
        integrator.step(1);
    }
}

/**
 * Test a TwoParticleAverageSite virtual site.
 */
void testTwoParticleAverage() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(2, new TwoParticleAverageSite(0, 1, 0.8, 0.2));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        ASSERT_EQUAL_VEC(pos[0]*0.8+pos[1]*0.2, pos[2], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.3*0.8, 0, 0), state.getForces()[0], 1e-4);
        ASSERT_EQUAL_VEC(Vec3(0.2+0.3*0.2, 0, 0), state.getForces()[1], 1e-4);
        integrator.step(1);
    }
}

/**
 * Test a ThreeParticleAverageSite virtual site.
 */
void testThreeParticleAverage() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new ThreeParticleAverageSite(0, 1, 2, 0.2, 0.3, 0.5));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    params[0] = 0.4;
    forceField->addParticle(3, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        ASSERT_EQUAL_VEC(pos[0]*0.2+pos[1]*0.3+pos[2]*0.5, pos[3], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.4*0.2, 0, 0), state.getForces()[0], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.2+0.4*0.3, 0, 0), state.getForces()[1], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.3+0.4*0.5, 0, 0), state.getForces()[2], 1e-5);
        integrator.step(1);
    }
}

/**
 * Test an OutOfPlaneSite virtual site.
 */
void testOutOfPlane() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new OutOfPlaneSite(0, 1, 2, 0.3, 0.4, 0.5));
    CustomExternalForce* forceField = new CustomExternalForce("-a*x");
    system.addForce(forceField);
    forceField->addPerParticleParameter("a");
    vector<double> params(1);
    params[0] = 0.1;
    forceField->addParticle(0, params);
    params[0] = 0.2;
    forceField->addParticle(1, params);
    params[0] = 0.3;
    forceField->addParticle(2, params);
    params[0] = 0.4;
    forceField->addParticle(3, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1, 0, 0);
    positions[2] = Vec3(0, 1, 0);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        Vec3 v12 = pos[1]-pos[0];
        Vec3 v13 = pos[2]-pos[0];
        Vec3 cross = v12.cross(v13);
        ASSERT_EQUAL_VEC(pos[0]+v12*0.3+v13*0.4+cross*0.5, pos[3], 1e-5);
        const vector<Vec3>& f = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0.1+0.2+0.3+0.4, 0, 0), f[0]+f[1]+f[2], 1e-5);
        Vec3 f2(0.4*0.3, 0.4*0.5*v13[2], -0.4*0.5*v13[1]);
        Vec3 f3(0.4*0.4, -0.4*0.5*v12[2], 0.4*0.5*v12[1]);
        ASSERT_EQUAL_VEC(Vec3(0.1+0.4, 0, 0)-f2-f3, f[0], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.2, 0, 0)+f2, f[1], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(0.3, 0, 0)+f3, f[2], 1e-5);
        integrator.step(1);
    }
}

/**
 * Test a LocalCoordinatesSite virtual site.
 */
void testLocalCoordinates() {
    const Vec3 originWeights(0.2, 0.3, 0.5);
    const Vec3 xWeights(-1.0, 0.5, 0.5);
    const Vec3 yWeights(0.0, -1.0, 1.0);
    const Vec3 localPosition(0.4, 0.3, 0.2);
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(3, new LocalCoordinatesSite(0, 1, 2, originWeights, xWeights, yWeights, localPosition));
    CustomExternalForce* forceField = new CustomExternalForce("2*x^2+3*y^2+4*z^2");
    system.addForce(forceField);
    vector<double> params;
    forceField->addParticle(0, params);
    forceField->addParticle(1, params);
    forceField->addParticle(2, params);
    forceField->addParticle(3, params);
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    vector<Vec3> positions(4), positions2(4), positions3(4);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < 100; i++) {
        // Set the particles at random positions.
        
        Vec3 xdir, ydir, zdir;
        do {
            for (int j = 0; j < 3; j++)
                positions[j] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
            xdir = positions[0]*xWeights[0] + positions[1]*xWeights[1] + positions[2]*xWeights[2];
            ydir = positions[0]*yWeights[0] + positions[1]*yWeights[1] + positions[2]*yWeights[2];
            zdir = xdir.cross(ydir);
            if (sqrt(xdir.dot(xdir)) > 0.1 && sqrt(ydir.dot(ydir)) > 0.1 && sqrt(zdir.dot(zdir)) > 0.1)
                break; // These positions give a reasonable coordinate system.
        } while (true);
        context.setPositions(positions);
        context.applyConstraints(0.0001);
        
        // See if the virtual site is positioned correctly.
        
        State state = context.getState(State::Positions | State::Forces);
        const vector<Vec3>& pos = state.getPositions();
        Vec3 origin = pos[0]*originWeights[0] + pos[1]*originWeights[1] + pos[2]*originWeights[2];
        xdir /= sqrt(xdir.dot(xdir));
        zdir /= sqrt(zdir.dot(zdir));
        ydir = zdir.cross(xdir);
        ASSERT_EQUAL_VEC(origin+xdir*localPosition[0]+ydir*localPosition[1]+zdir*localPosition[2], pos[3], 1e-5);

        // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

        double norm = 0.0;
        for (int i = 0; i < 3; ++i) {
            Vec3 f = state.getForces()[i];
            norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        }
        norm = std::sqrt(norm);
        const double delta = 1e-2;
        double step = 0.5*delta/norm;
        for (int i = 0; i < 3; ++i) {
            Vec3 p = positions[i];
            Vec3 f = state.getForces()[i];
            positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
            positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
        }
        context.setPositions(positions2);
        context.applyConstraints(0.0001);
        State state2 = context.getState(State::Energy);
        context.setPositions(positions3);
        context.applyConstraints(0.0001);
        State state3 = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/delta, 1e-3)
    }
}

/**
 * Make sure that energy, linear momentum, and angular momentum are all conserved
 * when using virtual sites.
 */
void testConservationLaws() {
    System system;
    NonbondedForce* forceField = new NonbondedForce();
    system.addForce(forceField);
    vector<Vec3> positions;
    
    // Create a linear molecule with a TwoParticleAverage virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(2, new TwoParticleAverageSite(0, 1, 0.4, 0.6));
    system.addConstraint(0, 1, 2.0);
    for (int i = 0; i < 3; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i, j, 0, 1, 0);
    }
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(2, 0, 0));
    positions.push_back(Vec3());
    
    // Create a planar molecule with a ThreeParticleAverage virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(6, new ThreeParticleAverageSite(3, 4, 5, 0.3, 0.5, 0.2));
    system.addConstraint(3, 4, 1.0);
    system.addConstraint(3, 5, 1.0);
    system.addConstraint(4, 5, sqrt(2.0));
    for (int i = 0; i < 4; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i+3, j+3, 0, 1, 0);
    }
    positions.push_back(Vec3(0, 0, 1));
    positions.push_back(Vec3(1, 0, 1));
    positions.push_back(Vec3(0, 1, 1));
    positions.push_back(Vec3());
    
    // Create a tetrahedral molecule with an OutOfPlane virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(10, new OutOfPlaneSite(7, 8, 9, 0.3, 0.5, 0.2));
    system.addConstraint(7, 8, 1.0);
    system.addConstraint(7, 9, 1.0);
    system.addConstraint(8, 9, sqrt(2.0));
    for (int i = 0; i < 4; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i+7, j+7, 0, 1, 0);
    }
    positions.push_back(Vec3(1, 0, -1));
    positions.push_back(Vec3(2, 0, -1));
    positions.push_back(Vec3(1, 1, -1));
    positions.push_back(Vec3());
    
    // Create a molecule with a LocalCoordinatesSite virtual site.
    
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.setVirtualSite(14, new LocalCoordinatesSite(11, 12, 13, Vec3(0.3, 0.3, 0.4), Vec3(1.0, -0.5, -0.5), Vec3(0, -1.0, 1.0), Vec3(0.2, 0.2, 1.0)));
    system.addConstraint(11, 12, 1.0);
    system.addConstraint(11, 13, 1.0);
    system.addConstraint(12, 13, sqrt(2.0));
    for (int i = 0; i < 4; i++) {
        forceField->addParticle(0, 1, 10);
        for (int j = 0; j < i; j++)
            forceField->addException(i+11, j+11, 0, 1, 0);
    }
    positions.push_back(Vec3(1, 2, 0));
    positions.push_back(Vec3(2, 2, 0));
    positions.push_back(Vec3(1, 3, 0));
    positions.push_back(Vec3());

    // Simulate it and check conservation laws.
    
    VerletIntegrator integrator(0.002);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    int numParticles = system.getNumParticles();
    double initialEnergy;
    Vec3 initialMomentum, initialAngularMomentum;
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions | State::Velocities | State::Forces | State::Energy);
        const vector<Vec3>& pos = state.getPositions();
        const vector<Vec3>& vel = state.getVelocities();
        const vector<Vec3>& f = state.getForces();
        double energy = state.getPotentialEnergy();
        for (int j = 0; j < numParticles; j++) {
            Vec3 v = vel[j] + f[j]*0.5*integrator.getStepSize();
            energy += 0.5*system.getParticleMass(j)*v.dot(v);
        }
        if (i == 0)
            initialEnergy = energy;
        else
            ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
        Vec3 momentum;
        for (int j = 0; j < numParticles; j++)
            momentum += vel[j]*system.getParticleMass(j);
        if (i == 0)
            initialMomentum = momentum;
        else
            ASSERT_EQUAL_VEC(initialMomentum, momentum, 0.02);
        Vec3 angularMomentum;
        for (int j = 0; j < numParticles; j++)
            angularMomentum += pos[j].cross(vel[j])*system.getParticleMass(j);
        if (i == 0)
            initialAngularMomentum = angularMomentum;
        else
            ASSERT_EQUAL_VEC(initialAngularMomentum, angularMomentum, 0.03);
        integrator.step(1);
    }
}

/**
 * Make sure that atom reordering respects virtual sites.
 */
void testReordering() {
    const double cutoff = 2.0;
    const double boxSize = 20.0;
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    system.addForce(nonbonded);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    nonbonded->setCutoffDistance(cutoff);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    
    // Create linear molecules with TwoParticleAverage virtual sites.
    
    for (int i = 0; i < 50; i++) {
        int start = system.getNumParticles();
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(0.0);
        system.setVirtualSite(start+2, new TwoParticleAverageSite(start, start+1, 0.4, 0.6));
        system.addConstraint(start, start+1, 2.0);
        for (int i = 0; i < 3; i++) {
            nonbonded->addParticle(0, 0.2, 1);
            for (int j = 0; j < i; j++)
                nonbonded->addException(start+i, start+j, 0, 1, 0);
        }
        Vec3 pos(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions.push_back(pos);
        positions.push_back(pos+Vec3(2, 0, 0));
        positions.push_back(Vec3());
    }
    
    // Create planar molecules with ThreeParticleAverage virtual sites.
    
    for (int i = 0; i < 50; i++) {
        int start = system.getNumParticles();
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(0.0);
        system.setVirtualSite(start+3, new ThreeParticleAverageSite(start, start+1, start+2, 0.3, 0.5, 0.2));
        system.addConstraint(start, start+1, 1.0);
        system.addConstraint(start, start+2, 1.0);
        system.addConstraint(start+1, start+2, sqrt(2.0));
        for (int i = 0; i < 4; i++) {
            nonbonded->addParticle(0, 0.2, 1);
            for (int j = 0; j < i; j++)
                nonbonded->addException(start+i, start+j, 0, 1, 0);
        }
        Vec3 pos(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions.push_back(pos);
        positions.push_back(pos+Vec3(1, 0, 0));
        positions.push_back(pos+Vec3(0, 1, 0));
        positions.push_back(Vec3());
    }
    
    // Create tetrahedral molecules with OutOfPlane virtual sites.
    
    for (int i = 0; i < 50; i++) {
        int start = system.getNumParticles();
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(1.0);
        system.addParticle(0.0);
        system.setVirtualSite(start+3, new OutOfPlaneSite(start, start+1, start+2, 0.3, 0.5, 0.2));
        system.addConstraint(start, start+1, 1.0);
        system.addConstraint(start, start+2, 1.0);
        system.addConstraint(start+1, start+2, sqrt(2.0));
        for (int i = 0; i < 4; i++) {
            nonbonded->addParticle(0, 0.2, 1);
            for (int j = 0; j < i; j++)
                nonbonded->addException(start+i, start+j, 0, 1, 0);
        }
        Vec3 pos(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions.push_back(pos);
        positions.push_back(pos+Vec3(1, 0, 0));
        positions.push_back(pos+Vec3(0, 1, 0));
        positions.push_back(Vec3());
    }

    // Simulate it and check conservation laws.
    
    LangevinIntegrator integrator(300.0, 0.1, 0.002);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.applyConstraints(0.0001);
    for (int i = 0; i < 1000; i++) {
        State state = context.getState(State::Positions);
        const vector<Vec3>& pos = state.getPositions();
        for (int j = 0; j < 150; j += 3)
            ASSERT_EQUAL_VEC(pos[j]*0.4+pos[j+1]*0.6, pos[j+2], 1e-5);
        for (int j = 150; j < 350; j += 4)
            ASSERT_EQUAL_VEC(pos[j]*0.3+pos[j+1]*0.5+pos[j+2]*0.2, pos[j+3], 1e-5);
        for (int j = 350; j < 550; j += 4) {
            Vec3 v12 = pos[j+1]-pos[j];
            Vec3 v13 = pos[j+2]-pos[j];
            Vec3 cross = v12.cross(v13);
            ASSERT_EQUAL_VEC(pos[j]+v12*0.3+v13*0.5+cross*0.2, pos[j+3], 1e-5);
        }
        integrator.step(1);
    }
}

/**
 * Test a System where multiple virtual sites are all calculated from the same particles.
 */
void testOverlappingSites() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    NonbondedForce* nonbonded = new NonbondedForce();
    system.addForce(nonbonded);
    nonbonded->addParticle(1.0, 0.0, 0.0);
    nonbonded->addParticle(-0.5, 0.0, 0.0);
    nonbonded->addParticle(-0.5, 0.0, 0.0);
    vector<Vec3> positions;
    positions.push_back(Vec3(0, 0, 0));
    positions.push_back(Vec3(10, 0, 0));
    positions.push_back(Vec3(0, 10, 0));
    for (int i = 0; i < 20; i++) {
        system.addParticle(0.0);
        double u = 0.1*((i+1)%4);
        double v = 0.05*i;
        system.setVirtualSite(3+i, new ThreeParticleAverageSite(0, 1, 2, u, v, 1-u-v));
        nonbonded->addParticle(i%2 == 0 ? -1.0 : 1.0, 0.0, 0.0);
        positions.push_back(Vec3());
    }
    VerletIntegrator i1(0.002);
    VerletIntegrator i2(0.002);
    Context c1(system, i1, Platform::getPlatformByName("Reference"));
    Context c2(system, i2, platform);
    c1.setPositions(positions);
    c2.setPositions(positions);
    c1.applyConstraints(0.0001);
    c2.applyConstraints(0.0001);
    State s1 = c1.getState(State::Positions | State::Forces);
    State s2 = c2.getState(State::Positions | State::Forces);
    for (int i = 0; i < system.getNumParticles(); i++)
        ASSERT_EQUAL_VEC(s1.getPositions()[i], s2.getPositions()[i], 1e-5);
    for (int i = 0; i < 3; i++)
        ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], 1e-5);
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        testMasslessParticle();
        testTwoParticleAverage();
        testThreeParticleAverage();
        testOutOfPlane();
        testLocalCoordinates();
        testConservationLaws();
        testReordering();
        testOverlappingSites();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
