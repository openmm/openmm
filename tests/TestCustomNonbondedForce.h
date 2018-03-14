
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2018 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSimpleExpression() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("-0.1*r^3");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    double force = 0.1*3*(2*2);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(-0.1*(2*2*2), state.getPotentialEnergy(), TOL);
}

void testParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("scale*a*(r*b)^3; a=a1*a2; b=c+b1+b2");
    forceField->addPerParticleParameter("a");
    forceField->addPerParticleParameter("b");
    forceField->addGlobalParameter("scale", 3.0);
    forceField->addGlobalParameter("c", -1.0);
    vector<double> params(2);
    params[0] = 1.5;
    params[1] = 2.0;
    forceField->addParticle(params);
    params[0] = 2.0;
    params[1] = 3.0;
    forceField->addParticle(params);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    context.setParameter("scale", 1.0);
    context.setParameter("c", 0.0);
    State state = context.getState(State::Forces | State::Energy);
    vector<Vec3> forces = state.getForces();
    double force = -3.0*3*5.0*(10*10);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(3.0*(10*10*10), state.getPotentialEnergy(), TOL);
    
    // Try changing the global parameters and make sure it's still correct.
    
    context.setParameter("scale", 1.5);
    context.setParameter("c", 1.0);
    state = context.getState(State::Forces | State::Energy);
    forces = state.getForces();
    force = -1.5*3.0*3*6.0*(12*12);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(1.5*3.0*(12*12*12), state.getPotentialEnergy(), TOL);
    
    // Try changing the per-particle parameters and make sure it's still correct.
    
    params[0] = 1.6;
    params[1] = 2.1;
    forceField->setParticleParameters(0, params);
    params[0] = 1.9;
    params[1] = 2.8;
    forceField->setParticleParameters(1, params);
    forceField->updateParametersInContext(context);
    state = context.getState(State::Forces | State::Energy);
    forces = state.getForces();
    force = -1.5*1.6*1.9*3*5.9*(11.8*11.8);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(1.5*1.6*1.9*(11.8*11.8*11.8), state.getPotentialEnergy(), TOL);
}

void testManyParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("(a1*a2+b1*b2+c1*c2+d1*d2+e1*e2)*r");
    forceField->addPerParticleParameter("a");
    forceField->addPerParticleParameter("b");
    forceField->addPerParticleParameter("c");
    forceField->addPerParticleParameter("d");
    forceField->addPerParticleParameter("e");
    vector<double> params(5);
    params[0] = 1.0;
    params[1] = 2.0;
    params[2] = 3.0;
    params[3] = 4.0;
    params[4] = 5.0;
    forceField->addParticle(params);
    params[0] = 1.1;
    params[1] = 1.2;
    params[2] = 1.3;
    params[3] = 1.4;
    params[4] = 1.5;
    forceField->addParticle(params);
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    vector<Vec3> forces = state.getForces();
    double force = 1*1.1 + 2*1.2 + 3*1.3 + 4*1.4 + 5*1.5;
    ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_TOL(2*force, state.getPotentialEnergy(), TOL);
}

void testExclusions() {
    System system;
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("a*r; a=a1+a2");
    nonbonded->addPerParticleParameter("a");
    vector<double> params(1);
    vector<Vec3> positions(4);
    for (int i = 0; i < 4; i++) {
        system.addParticle(1.0);
        params[0] = i+1;
        nonbonded->addParticle(params);
        positions[i] = Vec3(i, 0, 0);
    }
    nonbonded->addExclusion(0, 1);
    nonbonded->addExclusion(1, 2);
    nonbonded->addExclusion(2, 3);
    nonbonded->addExclusion(0, 2);
    nonbonded->addExclusion(1, 3);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(1+4, 0, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(-(1+4), 0, 0), forces[3], TOL);
    ASSERT_EQUAL_TOL((1+4)*3.0, state.getPotentialEnergy(), TOL);
}

void testCutoff() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("r");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    forceField->setCutoffDistance(2.5);
    system.addForce(forceField);
    ASSERT(!forceField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, 1, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -1, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(2.0+1.0, state.getPotentialEnergy(), TOL);
}

void testPeriodic() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("r");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    forceField->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    forceField->setCutoffDistance(2.0);
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    system.addForce(forceField);
    ASSERT(forceField->usesPeriodicBoundaryConditions());
    ASSERT(system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2.1, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, -2, 0), forces[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 2, 0), forces[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[2], TOL);
    ASSERT_EQUAL_TOL(1.9+1+0.9, state.getPotentialEnergy(), TOL);
}

void testTriclinic() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    Vec3 a(3.1, 0, 0);
    Vec3 b(0.4, 3.5, 0);
    Vec3 c(-0.1, -0.5, 4.0);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("r");
    nonbonded->addParticle(vector<double>());
    nonbonded->addParticle(vector<double>());
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    const double cutoff = 1.5;
    nonbonded->setCutoffDistance(cutoff);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int iteration = 0; iteration < 50; iteration++) {
        // Generate random positions for the two particles.

        positions[0] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        positions[1] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        context.setPositions(positions);

        // Loop over all possible periodic copies and find the nearest one.

        Vec3 delta;
        double distance2 = 100.0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                for (int k = -1; k < 2; k++) {
                    Vec3 d = positions[1]-positions[0]+a*i+b*j+c*k;
                    if (d.dot(d) < distance2) {
                        delta = d;
                        distance2 = d.dot(d);
                    }
                }
        double distance = sqrt(distance2);

        // See if the force and energy are correct.

        State state = context.getState(State::Forces | State::Energy);
        if (distance >= cutoff) {
            ASSERT_EQUAL(0.0, state.getPotentialEnergy());
            ASSERT_EQUAL_VEC(Vec3(0, 0, 0), state.getForces()[0], 0);
            ASSERT_EQUAL_VEC(Vec3(0, 0, 0), state.getForces()[1], 0);
        }
        else {
            const Vec3 force = delta/sqrt(delta.dot(delta));
            ASSERT_EQUAL_TOL(distance, state.getPotentialEnergy(), TOL);
            ASSERT_EQUAL_VEC(force, state.getForces()[0], TOL);
            ASSERT_EQUAL_VEC(-force, state.getForces()[1], TOL);
        }
    }
}

void testContinuous1DFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r)+1");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < 21; i++)
        table.push_back(sin(0.25*i));
    forceField->addTabulatedFunction("fn", new Continuous1DFunction(table, 1.0, 6.0));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 1; i < 30; i++) {
        double x = (7.0/30.0)*i;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double force = (x < 1.0 || x > 6.0 ? 0.0 : -cos(x-1.0));
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : sin(x-1.0))+1.0;
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], 0.1);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], 0.1);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.02);
    }
    for (int i = 1; i < 20; i++) {
        double x = 0.25*i+1.0;
        positions[1] = Vec3(x, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Energy);
        double energy = (x < 1.0 || x > 6.0 ? 0.0 : sin(x-1.0))+1.0;
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 1e-4);
    }
}

void testContinuous2DFunction() {
    const int xsize = 20;
    const int ysize = 21;
    const double xmin = 0.4;
    const double xmax = 1.5;
    const double ymin = 0.0;
    const double ymax = 2.1;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r,a)+1");
    forceField->addGlobalParameter("a", 0.0);
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table(xsize*ysize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            double x = xmin + i*(xmax-xmin)/xsize;
            double y = ymin + j*(ymax-ymin)/ysize;
            table[i+xsize*j] = sin(0.25*x)*cos(0.33*y);
        }
    }
    forceField->addTabulatedFunction("fn", new Continuous2DFunction(xsize, ysize, table, xmin, xmax, ymin, ymax));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (double x = xmin-0.15; x < xmax+0.2; x += 0.1) {
        for (double y = ymin-0.15; y < ymax+0.2; y += 0.1) {
            positions[1] = Vec3(x, 0, 0);
            context.setParameter("a", y);
            context.setPositions(positions);
            State state = context.getState(State::Forces | State::Energy);
            const vector<Vec3>& forces = state.getForces();
            double energy = 1;
            double force = 0;
            if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
                energy = sin(0.25*x)*cos(0.33*y)+1.0;
                force = -0.25*cos(0.25*x)*cos(0.33*y);
            }
            ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], 0.1);
            ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], 0.1);
            ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.02);
        }
    }
}

void testContinuous3DFunction() {
    const int xsize = 10;
    const int ysize = 11;
    const int zsize = 12;
    const double xmin = 0.6;
    const double xmax = 1.1;
    const double ymin = 0.0;
    const double ymax = 0.7;
    const double zmin = 0.2;
    const double zmax = 0.9;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r,a,b)+1");
    forceField->addGlobalParameter("a", 0.0);
    forceField->addGlobalParameter("b", 0.0);
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table(xsize*ysize*zsize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            for (int k = 0; k < zsize; k++) {
                double x = xmin + i*(xmax-xmin)/xsize;
                double y = ymin + j*(ymax-ymin)/ysize;
                double z = zmin + k*(zmax-zmin)/zsize;
                table[i+xsize*j+xsize*ysize*k] = sin(0.25*x)*cos(0.33*y)*(1+z);
            }
        }
    }
    forceField->addTabulatedFunction("fn", new Continuous3DFunction(xsize, ysize, zsize, table, xmin, xmax, ymin, ymax, zmin, zmax));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (double x = xmin-0.15; x < xmax+0.2; x += 0.1) {
        for (double y = ymin-0.15; y < ymax+0.2; y += 0.1) {
            for (double z = zmin-0.15; z < zmax+0.2; z += 0.1) {
                positions[1] = Vec3(x, 0, 0);
                context.setParameter("a", y);
                context.setParameter("b", z);
                context.setPositions(positions);
                State state = context.getState(State::Forces | State::Energy);
                const vector<Vec3>& forces = state.getForces();
                double energy = 1;
                double force = 0;
                if (x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin && z <= zmax) {
                    energy = sin(0.25*x)*cos(0.33*y)*(1.0+z)+1.0;
                    force = -0.25*cos(0.25*x)*cos(0.33*y)*(1.0+z);
                }
                ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], 0.1);
                ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], 0.1);
                ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.05);
            }
        }
    }
}

void testDiscrete1DFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r-1)+1");
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < 21; i++)
        table.push_back(sin(0.25*i));
    forceField->addTabulatedFunction("fn", new Discrete1DFunction(table));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 0; i < (int) table.size(); i++) {
        positions[1] = Vec3(i+1, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[0], 1e-6);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], 1e-6);
        ASSERT_EQUAL_TOL(table[i]+1.0, state.getPotentialEnergy(), 1e-6);
    }
}

void testDiscrete2DFunction() {
    const int xsize = 10;
    const int ysize = 5;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r-1,a)+1");
    forceField->addGlobalParameter("a", 0.0);
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++)
            table.push_back(sin(0.25*i)+cos(0.33*j));
    forceField->addTabulatedFunction("fn", new Discrete2DFunction(xsize, ysize, table));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 0; i < (int) table.size(); i++) {
        positions[1] = Vec3((i%xsize)+1, 0, 0);
        context.setPositions(positions);
        context.setParameter("a", i/xsize);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[0], 1e-6);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], 1e-6);
        ASSERT_EQUAL_TOL(table[i]+1.0, state.getPotentialEnergy(), 1e-6);
    }
}

void testDiscrete3DFunction() {
    const int xsize = 8;
    const int ysize = 5;
    const int zsize = 6;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r-1,a,b)+1");
    forceField->addGlobalParameter("a", 0.0);
    forceField->addGlobalParameter("b", 0.0);
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < xsize; i++)
        for (int j = 0; j < ysize; j++)
            for (int k = 0; k < zsize; k++)
                table.push_back(sin(0.25*i)+cos(0.33*j)+0.12345*k);
    forceField->addTabulatedFunction("fn", new Discrete3DFunction(xsize, ysize, zsize, table));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 0; i < (int) table.size(); i++) {
        positions[1] = Vec3((i%xsize)+1, 0, 0);
        context.setPositions(positions);
        context.setParameter("a", (i/xsize)%ysize);
        context.setParameter("b", i/(xsize*ysize));
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[0], 1e-6);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], 1e-6);
        ASSERT_EQUAL_TOL(table[i]+1.0, state.getPotentialEnergy(), 1e-6);
    }
}

void testCoulombLennardJones() {
    const int numMolecules = 300;
    const int numParticles = numMolecules*2;
    const double boxSize = 20.0;

    // Create two systems: one with a NonbondedForce, and one using a CustomNonbondedForce to implement the same interaction.

    System standardSystem;
    System customSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customSystem.addParticle(1.0);
    }
    NonbondedForce* standardNonbonded = new NonbondedForce();
    CustomNonbondedForce* customNonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r; q=q1*q2; sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    customNonbonded->addPerParticleParameter("q");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<double> params(3);
    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {
            standardNonbonded->addParticle(1.0, 0.2, 0.1);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.1;
            customNonbonded->addParticle(params);
            standardNonbonded->addParticle(-1.0, 0.1, 0.1);
            params[0] = -1.0;
            params[1] = 0.1;
            customNonbonded->addParticle(params);
        }
        else {
            standardNonbonded->addParticle(1.0, 0.2, 0.2);
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.2;
            customNonbonded->addParticle(params);
            standardNonbonded->addParticle(-1.0, 0.1, 0.2);
            params[0] = -1.0;
            params[1] = 0.1;
            customNonbonded->addParticle(params);
        }
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        standardNonbonded->addException(2*i, 2*i+1, 0.0, 1.0, 0.0);
        customNonbonded->addExclusion(2*i, 2*i+1);
    }
    standardNonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);
    customNonbonded->setNonbondedMethod(CustomNonbondedForce::NoCutoff);
    standardSystem.addForce(standardNonbonded);
    customSystem.addForce(customNonbonded);
    ASSERT(!customNonbonded->usesPeriodicBoundaryConditions());
    ASSERT(!customSystem.usesPeriodicBoundaryConditions());
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context context1(standardSystem, integrator1, platform);
    Context context2(customSystem, integrator2, platform);
    context1.setPositions(positions);
    context2.setPositions(positions);
    context1.setVelocities(velocities);
    context2.setVelocities(velocities);
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
    }
}

void testSwitchingFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("10/r^2");
    vector<double> params;
    nonbonded->addParticle(params);
    nonbonded->addParticle(params);
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded->setCutoffDistance(2.0);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setSwitchingDistance(1.5);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    
    // Compute the interaction at various distances.
    
    for (double r = 1.0; r < 2.5; r += 0.1) {
        positions[1] = Vec3(r, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        
        // See if the energy is correct.
        
        double expectedEnergy = 10/(r*r);
        double switchValue;
        if (r <= 1.5)
            switchValue = 1;
        else if (r >= 2.0)
            switchValue = 0;
        else {
            double t = (r-1.5)/0.5;
            switchValue = 1+t*t*t*(-10+t*(15-t*6));
        }
        ASSERT_EQUAL_TOL(switchValue*expectedEnergy, state.getPotentialEnergy(), TOL);
        
        // See if the force is the gradient of the energy.
        
        double delta = 1e-3;
        positions[1] = Vec3(r-delta, 0, 0);
        context.setPositions(positions);
        double e1 = context.getState(State::Energy).getPotentialEnergy();
        positions[1] = Vec3(r+delta, 0, 0);
        context.setPositions(positions);
        double e2 = context.getState(State::Energy).getPotentialEnergy();
        ASSERT_EQUAL_TOL((e2-e1)/(2*delta), state.getForces()[0][0], 1e-3);
    }
}

void testLongRangeCorrection() {
    // Create a box of particles.

    int gridSize = 5;
    int numParticles = gridSize*gridSize*gridSize;
    double boxSize = gridSize*0.7;
    double cutoff = boxSize/3;
    System standardSystem;
    System customSystem;
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    NonbondedForce* standardNonbonded = new NonbondedForce();
    CustomNonbondedForce* customNonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    int index = 0;
    vector<double> params1(2);
    params1[0] = 1.1;
    params1[1] = 0.5;
    vector<double> params2(2);
    params2[0] = 1;
    params2[1] = 1;
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                standardSystem.addParticle(1.0);
                customSystem.addParticle(1.0);
                if (index%2 == 0) {
                    standardNonbonded->addParticle(0, params1[0], params1[1]);
                    customNonbonded->addParticle(params1);
                }
                else {
                    standardNonbonded->addParticle(0, params2[0], params2[1]);
                    customNonbonded->addParticle(params2);
                }
                positions[index] = Vec3(i*boxSize/gridSize, j*boxSize/gridSize, k*boxSize/gridSize);
                index++;
            }
    standardNonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    customNonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    standardNonbonded->setCutoffDistance(cutoff);
    customNonbonded->setCutoffDistance(cutoff);
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    customSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    standardNonbonded->setUseDispersionCorrection(true);
    customNonbonded->setUseLongRangeCorrection(true);
    standardNonbonded->setUseSwitchingFunction(true);
    customNonbonded->setUseSwitchingFunction(true);
    standardNonbonded->setSwitchingDistance(0.8*cutoff);
    customNonbonded->setSwitchingDistance(0.8*cutoff);
    standardSystem.addForce(standardNonbonded);
    customSystem.addForce(customNonbonded);

    // Compute the correction for the standard force.

    Context context1(standardSystem, integrator1, platform);
    context1.setPositions(positions);
    double standardEnergy1 = context1.getState(State::Energy).getPotentialEnergy();
    standardNonbonded->setUseDispersionCorrection(false);
    context1.reinitialize();
    context1.setPositions(positions);
    double standardEnergy2 = context1.getState(State::Energy).getPotentialEnergy();

    // Compute the correction for the custom force.

    Context context2(customSystem, integrator2, platform);
    context2.setPositions(positions);
    double customEnergy1 = context2.getState(State::Energy).getPotentialEnergy();
    customNonbonded->setUseLongRangeCorrection(false);
    context2.reinitialize();
    context2.setPositions(positions);
    double customEnergy2 = context2.getState(State::Energy).getPotentialEnergy();

    // See if they agree.

    ASSERT_EQUAL_TOL(standardEnergy1-standardEnergy2, customEnergy1-customEnergy2, 1e-4);
}

void testInteractionGroups() {
    const int numParticles = 6;
    System system;
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("v1+v2");
    nonbonded->addPerParticleParameter("v");
    vector<double> params(1, 0.001);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(params);
        params[0] *= 10;
    }
    set<int> set1, set2, set3, set4;
    set1.insert(2);
    set2.insert(0);
    set2.insert(1);
    set2.insert(2);
    set2.insert(3);
    set2.insert(4);
    set2.insert(5);
    nonbonded->addInteractionGroup(set1, set2); // Particle 2 interacts with every other particle.
    set3.insert(0);
    set3.insert(1);
    set4.insert(4);
    set4.insert(5);
    nonbonded->addInteractionGroup(set3, set4); // Particles 0 and 1 interact with 4 and 5.
    nonbonded->addExclusion(1, 2); // Add an exclusion to make sure it gets skipped.
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    double expectedEnergy = 331.423; // Each digit is the number of interactions a particle particle is involved in.
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), TOL);
}

void testLargeInteractionGroup() {
    const int numMolecules = 300;
    const int numParticles = numMolecules*2;
    const double boxSize = 20.0;
    
    // Create a large system.
    
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r; q=q1*q2; sigma=0.5*(sigma1+sigma2); eps=sqrt(eps1*eps2)");
    nonbonded->addPerParticleParameter("q");
    nonbonded->addPerParticleParameter("sigma");
    nonbonded->addPerParticleParameter("eps");
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<double> params(3);
    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.1;
            nonbonded->addParticle(params);
            params[0] = -1.0;
            params[1] = 0.1;
            nonbonded->addParticle(params);
        }
        else {
            params[0] = 1.0;
            params[1] = 0.2;
            params[2] = 0.2;
            nonbonded->addParticle(params);
            params[0] = -1.0;
            params[1] = 0.1;
            nonbonded->addParticle(params);
        }
        positions[2*i] = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1] = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        nonbonded->addExclusion(2*i, 2*i+1);
    }
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    system.addForce(nonbonded);
    
    // Compute the forces.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces);
    
    // Modify the force so only one particle interacts with everything else.
    
    set<int> set1, set2;
    set1.insert(151);
    for (int i = 0; i < numParticles; i++)
        set2.insert(i);
    nonbonded->addInteractionGroup(set1, set2);
    context.reinitialize();
    context.setPositions(positions);
    State state2 = context.getState(State::Forces);
    
    // The force on that one particle should be the same.
    
    ASSERT_EQUAL_VEC(state1.getForces()[151], state2.getForces()[151], 1e-4);
    
    // Modify the interaction group so it includes all interactions.  This should now reproduce the original forces
    // on all atoms.

    for (int i = 0; i < numParticles; i++)
        set1.insert(i);
    nonbonded->setInteractionGroupParameters(0, set1, set2);
    context.reinitialize();
    context.setPositions(positions);
    State state3 = context.getState(State::Forces);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i], state3.getForces()[i], 1e-4);
    
    // Now make the interaction group empty.  The energy should then be zero.
    
    set1.clear();
    set2.clear();
    nonbonded->setInteractionGroupParameters(0, set1, set2);
    context.reinitialize();
    context.setPositions(positions);
    State state4 = context.getState(State::Energy);
    ASSERT_EQUAL(0.0, state4.getPotentialEnergy());
}

void testInteractionGroupLongRangeCorrection() {
    const int numParticles = 10;
    const double boxSize = 10.0;
    const double cutoff = 0.5;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("c1*c2*r^-4");
    nonbonded->addPerParticleParameter("c");
    vector<Vec3> positions(numParticles);
    vector<double> params(1);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        params[0] = (i%2 == 0 ? 1.1 : 2.0);
        nonbonded->addParticle(params);
        positions[i] = Vec3(0.5*i, 0, 0);
    }
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(cutoff);
    system.addForce(nonbonded);
    
    // Setup nonbonded groups.  They involve 1 interaction of type AA,
    // 2 of type BB, and 5 of type AB.
    
    set<int> set1, set2, set3, set4, set5;
    set1.insert(0);
    set1.insert(1);
    set1.insert(2);
    nonbonded->addInteractionGroup(set1, set1);
    set2.insert(3);
    set3.insert(4);
    set3.insert(6);
    set3.insert(8);
    nonbonded->addInteractionGroup(set2, set3);
    set4.insert(5);
    set5.insert(7);
    set5.insert(9);
    nonbonded->addInteractionGroup(set4, set5);
    
    // Compute energy with and without the correction.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    double energy1 = context.getState(State::Energy).getPotentialEnergy();
    nonbonded->setUseLongRangeCorrection(true);
    context.reinitialize();
    context.setPositions(positions);
    double energy2 = context.getState(State::Energy).getPotentialEnergy();
    
    // Check the result.
    
    double sum = (1.1*1.1 + 2*2.0*2.0 + 5*1.1*2.0)*2.0;
    int numPairs = (numParticles*(numParticles+1))/2;
    double expected = 2*M_PI*numParticles*numParticles*sum/(numPairs*boxSize*boxSize*boxSize);
    ASSERT_EQUAL_TOL(expected, energy2-energy1, 1e-4);
}

void testInteractionGroupTabulatedFunction() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* forceField = new CustomNonbondedForce("fn(r-1)+1");
    set<int> set1, set2;
    set1.insert(0);
    set2.insert(1);
    forceField->addInteractionGroup(set1, set2);
    forceField->addParticle(vector<double>());
    forceField->addParticle(vector<double>());
    vector<double> table;
    for (int i = 0; i < 21; i++)
        table.push_back(sin(0.25*i));
    forceField->addTabulatedFunction("fn", new Discrete1DFunction(table));
    system.addForce(forceField);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    for (int i = 0; i < (int) table.size(); i++) {
        positions[1] = Vec3(i+1, 0, 0);
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[0], 1e-6);
        ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forces[1], 1e-6);
        ASSERT_EQUAL_TOL(table[i]+1.0, state.getPotentialEnergy(), 1e-6);
    }
}

void testInteractionGroupWithCutoff() {
    const int numParticles = 1000;
    const double boxSize = 10.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    NonbondedForce* standard = new NonbondedForce();
    CustomNonbondedForce* custom = new CustomNonbondedForce("100/(r+0.1)");
    system.addForce(standard);
    system.addForce(custom);
    standard->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    custom->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    standard->setCutoffDistance(1.0);
    custom->setCutoffDistance(1.0);
    standard->setUseSwitchingFunction(true);
    custom->setUseSwitchingFunction(true);
    standard->setSwitchingDistance(0.9);
    custom->setSwitchingDistance(0.8);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(10.0);
        standard->addParticle(0.0, 0.2, 0.1);
        custom->addParticle();
        while (true) {
            positions[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*boxSize;
            bool tooClose = false;
            for (int j = 0; j < i; j++) {
                Vec3 delta = positions[i]-positions[j];
                if (delta.dot(delta) < 0.5*0.5)
                    tooClose = true;
            }
            if (!tooClose)
                break;
        }
    }
    set<int> set1, set2;
    for (int i = 0; i < 10; i++)
        set1.insert(2*i);
    for (int i = 0; i < numParticles; i++)
        set2.insert(i);
    custom->addInteractionGroup(set1, set2);
    custom->setForceGroup(1);
    
    // Try simulating it and see if energy is conserved (indicating that any optimizations
    // for combining the cutoff with the interaction group are behaving consistently).

    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    context.setVelocitiesToTemperature(100);
    ASSERT(context.getState(State::Energy, false, 1<<1).getPotentialEnergy() != 0.0);
    State initialState = context.getState(State::Energy);
    double initialEnergy = initialState.getPotentialEnergy()+initialState.getKineticEnergy();
    for (int i = 0; i < 100; i++) {
        integrator.step(10);
        State state = context.getState(State::Energy);
        double energy = state.getPotentialEnergy()+state.getKineticEnergy();
        ASSERT_EQUAL_TOL(initialEnergy, energy, 0.001);
    }
}

void testMultipleCutoffs() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    
    // Add multiple nonbonded forces that have different cutoffs.
    
    CustomNonbondedForce* nonbonded1 = new CustomNonbondedForce("2*r");
    nonbonded1->addParticle(vector<double>());
    nonbonded1->addParticle(vector<double>());
    nonbonded1->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded1->setCutoffDistance(2.5);
    system.addForce(nonbonded1);
    CustomNonbondedForce* nonbonded2 = new CustomNonbondedForce("3*r");
    nonbonded2->addParticle(vector<double>());
    nonbonded2->addParticle(vector<double>());
    nonbonded2->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded2->setCutoffDistance(2.9);
    nonbonded2->setForceGroup(1);
    system.addForce(nonbonded2);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);
    for (double r = 2.4; r < 3.2; r += 0.2) {
        positions[1][1] = r;
        context.setPositions(positions);
        double e1 = (r < 2.5 ? 2.0*r : 0.0);
        double e2 = (r < 2.9 ? 3.0*r : 0.0);
        double f1 = (r < 2.5 ? 2.0 : 0.0);
        double f2 = (r < 2.9 ? 3.0 : 0.0);
        
        // Check the first force.
        
        State state = context.getState(State::Forces | State::Energy, false, 1);
        ASSERT_EQUAL_VEC(Vec3(0, f1, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f1, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e1, state.getPotentialEnergy(), TOL);
        
        // Check the second force.
        
        state = context.getState(State::Forces | State::Energy, false, 2);
        ASSERT_EQUAL_VEC(Vec3(0, f2, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f2, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e2, state.getPotentialEnergy(), TOL);
        
        // Check the sum of both forces.

        state = context.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_VEC(Vec3(0, f1+f2, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f1-f2, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e1+e2, state.getPotentialEnergy(), TOL);
    }
}

void testMultipleSwitches() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    
    // Add multiple CustomNonbondedForces, one with a switching function and one without.
    
    CustomNonbondedForce* nonbonded1 = new CustomNonbondedForce("2*r");
    nonbonded1->addParticle(vector<double>());
    nonbonded1->addParticle(vector<double>());
    nonbonded1->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded1->setCutoffDistance(1.0);
    system.addForce(nonbonded1);
    CustomNonbondedForce* nonbonded2 = new CustomNonbondedForce("3*r");
    nonbonded2->addParticle(vector<double>());
    nonbonded2->addParticle(vector<double>());
    nonbonded2->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    nonbonded2->setCutoffDistance(1.0);
    nonbonded2->setSwitchingDistance(0.5);
    nonbonded2->setUseSwitchingFunction(true);
    nonbonded2->setForceGroup(1);
    system.addForce(nonbonded2);
    Context context(system, integrator, platform);
    double r = 0.8;
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, r, 0);
    context.setPositions(positions);
    double t = (r-0.5)/(1.0-0.5);
    double switchValue = 1+t*t*t*(-10+t*(15-t*6));
    double switchDeriv = t*t*(-30+t*(60-t*30))/(1.0-0.5);
    double e1 = 2.0*r;
    double e2 = 3.0*r*switchValue;
    double f1 = 2.0;
    double f2 = 3.0*switchValue + 3.0*switchDeriv*r;

    // Check the first force.

    State state = context.getState(State::Forces | State::Energy, false, 1);
    ASSERT_EQUAL_VEC(Vec3(0, f1, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -f1, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_TOL(e1, state.getPotentialEnergy(), TOL);

    // Check the second force.

    state = context.getState(State::Forces | State::Energy, false, 2);
    ASSERT_EQUAL_VEC(Vec3(0, f2, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -f2, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_TOL(e2, state.getPotentialEnergy(), TOL);

    // Check the sum of both forces.

    state = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_VEC(Vec3(0, f1+f2, 0), state.getForces()[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -f1-f2, 0), state.getForces()[1], TOL);
    ASSERT_EQUAL_TOL(e1+e2, state.getPotentialEnergy(), TOL);
}

void testIllegalVariable() {
    System system;
    system.addParticle(1.0);
    CustomNonbondedForce* force = new CustomNonbondedForce("r+none");
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
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("k*(r-r0)^2");
    nonbonded->addGlobalParameter("r0", 0.0);
    nonbonded->addGlobalParameter("k", 0.0);
    nonbonded->addEnergyParameterDerivative("k");
    nonbonded->addEnergyParameterDerivative("r0");
    vector<double> parameters;
    nonbonded->addParticle(parameters);
    nonbonded->addParticle(parameters);
    nonbonded->addParticle(parameters);
    nonbonded->addExclusion(0, 2);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    for (int i = 0; i < 10; i++) {
        double r0 = 0.1*i;
        double k = 10-i;
        context.setParameter("r0", r0);
        context.setParameter("k", k);
        State state = context.getState(State::ParameterDerivatives);
        map<string, double> derivs = state.getEnergyParameterDerivatives();
        double dEdr0 = -2*k*((2-r0)+(1-r0));
        double dEdk = (2-r0)*(2-r0) + (1-r0)*(1-r0);
        ASSERT_EQUAL_TOL(dEdr0, derivs["r0"], 1e-5);
        ASSERT_EQUAL_TOL(dEdk, derivs["k"], 1e-5);
    }
}

void testEnergyParameterDerivatives2() {
    // Create a box of particles.
    
    const int numParticles = 30;
    const double boxSize = 2.0;
    const double a = 1.0;
    const double delta = 1e-3;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("(r+a)^-4");
    system.addForce(nonbonded);
    nonbonded->addGlobalParameter("a", a);
    nonbonded->addEnergyParameterDerivative("a");
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    nonbonded->setSwitchingDistance(0.9);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setUseLongRangeCorrection(true);
    vector<Vec3> positions;
    vector<double> parameters;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        nonbonded->addParticle(parameters);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*boxSize);
    }
    
    // Compute the energy derivative and compare it to a finite difference approximation.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    map<string, double> derivs = context.getState(State::ParameterDerivatives).getEnergyParameterDerivatives();
    context.setParameter("a", a+delta);
    double energy1 = context.getState(State::Energy).getPotentialEnergy();
    context.setParameter("a", a-delta);
    double energy2 = context.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL((energy1-energy2)/(2*delta), derivs["a"], 1e-4);
}

void testEnergyParameterDerivativesWithGroups() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomNonbondedForce* nonbonded = new CustomNonbondedForce("k*(r-r0)^2");
    nonbonded->addGlobalParameter("r0", 0.0);
    nonbonded->addGlobalParameter("k", 0.0);
    nonbonded->addEnergyParameterDerivative("k");
    nonbonded->addEnergyParameterDerivative("r0");
    vector<double> parameters;
    nonbonded->addParticle(parameters);
    nonbonded->addParticle(parameters);
    nonbonded->addParticle(parameters);
    set<int> set1, set2;
    set1.insert(1);
    set2.insert(0);
    set2.insert(2);
    nonbonded->addInteractionGroup(set1, set2);
    system.addForce(nonbonded);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 2, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 0);
    context.setPositions(positions);
    for (int i = 0; i < 10; i++) {
        double r0 = 0.1*i;
        double k = 10-i;
        context.setParameter("r0", r0);
        context.setParameter("k", k);
        State state = context.getState(State::ParameterDerivatives);
        map<string, double> derivs = state.getEnergyParameterDerivatives();
        double dEdr0 = -2*k*((2-r0)+(1-r0));
        double dEdk = (2-r0)*(2-r0) + (1-r0)*(1-r0);
        ASSERT_EQUAL_TOL(dEdr0, derivs["r0"], 1e-5);
        ASSERT_EQUAL_TOL(dEdk, derivs["k"], 1e-5);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSimpleExpression();
        testParameters();
        testExclusions();
        testCutoff();
        testPeriodic();
        testTriclinic();
        testContinuous1DFunction();
        testContinuous2DFunction();
        testContinuous3DFunction();
        testDiscrete1DFunction();
        testDiscrete2DFunction();
        testDiscrete3DFunction();
        testCoulombLennardJones();
        testSwitchingFunction();
        testLongRangeCorrection();
        testInteractionGroups();
        testLargeInteractionGroup();
        testInteractionGroupLongRangeCorrection();
        testInteractionGroupTabulatedFunction();
        testInteractionGroupWithCutoff();
        testMultipleCutoffs();
        testMultipleSwitches();
        testIllegalVariable();
        testEnergyParameterDerivatives();
        testEnergyParameterDerivatives2();
        testEnergyParameterDerivativesWithGroups();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
