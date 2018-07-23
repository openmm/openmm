/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2017 Stanford University and the Authors.           *
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
#include "openmm/Context.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomCVForce.h"
#include "openmm/CustomExternalForce.h"
#include "openmm/CustomNonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

void testCVs() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    
    // Create a CustomCVForce with two collective variables.
    
    CustomCVForce* cv = new CustomCVForce("v1+3*v2");
    system.addForce(cv);
    CustomBondForce* v1 = new CustomBondForce("2*r");
    v1->addBond(0, 1);
    cv->addCollectiveVariable("v1", v1);
    CustomBondForce* v2 = new CustomBondForce("r");
    v2->addBond(0, 2);
    cv->addCollectiveVariable("v2", v2);
    ASSERT(!cv->usesPeriodicBoundaryConditions());
    
    // Create a context for evaluating it.
    
    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(1.5, 0, 0);
    positions[2] = Vec3(0, 0, 2.5);
    context.setPositions(positions);
    
    // Verify the values of the collective variables.
    
    vector<double> values;
    cv->getCollectiveVariableValues(context, values);
    ASSERT_EQUAL_TOL(2*1.5, values[0], 1e-6);
    ASSERT_EQUAL_TOL(2.5, values[1], 1e-6);
    
    // Verify the energy and forces.
    
    State state = context.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(2*1.5+3*2.5, state.getPotentialEnergy(), 1e-6);
    ASSERT_EQUAL_VEC(Vec3(2.0, 0.0, 3.0), state.getForces()[0], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(-2.0, 0.0, 0.0), state.getForces()[1], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(0.0, 0.0, -3.0), state.getForces()[2], 1e-6);
}

void testEnergyParameterDerivatives() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    
    // Create a CustomCVForce with one collective variable and two global parameters.
    // The CV in turn depends on a global parameter.
    
    CustomCVForce* cv = new CustomCVForce("v1*g1+g2");
    system.addForce(cv);
    CustomBondForce* v1 = new CustomBondForce("r*g3");
    v1->addGlobalParameter("g3", 2.0);
    v1->addEnergyParameterDerivative("g3");
    v1->addBond(0, 1);
    cv->addCollectiveVariable("v1", v1);
    cv->addGlobalParameter("g1", 1.5);
    cv->addGlobalParameter("g2", -1.0);
    cv->addEnergyParameterDerivative("g2");
    cv->addEnergyParameterDerivative("g1");
    
    // Create a context for evaluating it.
    
    VerletIntegrator integrator(1.0);
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2.0, 0, 0);
    context.setPositions(positions);

    // Verify the energy and parameter derivatives.

    State state = context.getState(State::Energy | State::ParameterDerivatives);
    ASSERT_EQUAL_TOL(4.0*1.5-1.0, state.getPotentialEnergy(), 1e-6);
    map<string, double> derivs = state.getEnergyParameterDerivatives();
    ASSERT_EQUAL_TOL(4.0, derivs["g1"], 1e-6);
    ASSERT_EQUAL_TOL(1.0, derivs["g2"], 1e-6);
    ASSERT_EQUAL_TOL(2.0*1.5, derivs["g3"], 1e-6);
}

void testTabulatedFunction() {
    const int xsize = 10;
    const int ysize = 11;
    const double xmin = 0.4;
    const double xmax = 1.1;
    const double ymin = 0.0;
    const double ymax = 0.95;
    System system;
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomCVForce* cv = new CustomCVForce("fn(x,y)+1");
    CustomExternalForce* v1 = new CustomExternalForce("x");
    v1->addParticle(0);
    cv->addCollectiveVariable("x", v1);
    CustomExternalForce* v2 = new CustomExternalForce("y");
    v2->addParticle(0);
    cv->addCollectiveVariable("y", v2);
    vector<double> table(xsize*ysize);
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            double x = xmin + i*(xmax-xmin)/xsize;
            double y = ymin + j*(ymax-ymin)/ysize;
            table[i+xsize*j] = sin(0.25*x)*cos(0.33*y);
        }
    }
    cv->addTabulatedFunction("fn", new Continuous2DFunction(xsize, ysize, table, xmin, xmax, ymin, ymax));
    system.addForce(cv);
    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    double scale = 1.0;
    for (int i = 0; i < 2; i++) {
        for (double x = xmin-0.15; x < xmax+0.2; x += 0.1) {
            for (double y = ymin-0.15; y < ymax+0.2; y += 0.1) {
                positions[0] = Vec3(x, y, 1.5);
                context.setPositions(positions);
                State state = context.getState(State::Forces | State::Energy);
                const vector<Vec3>& forces = state.getForces();
                double energy = 1;
                Vec3 force(0, 0, 0);
                if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
                    energy = scale*sin(0.25*x)*cos(0.33*y)+1;
                    force[0] = -scale*0.25*cos(0.25*x)*cos(0.33*y);
                    force[1] = scale*0.3*sin(0.25*x)*sin(0.33*y);
                }
                ASSERT_EQUAL_VEC(force, forces[0], 0.1);
                ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), 0.05);
            }
        }
        
        // Now update the tabulated function, call updateParametersInContext(),
        // and see if it's still correct.
        
        for (int j = 0; j < table.size(); j++)
            table[j] *= 2;
        dynamic_cast<Continuous2DFunction&>(cv->getTabulatedFunction(0)).setFunctionParameters(xsize, ysize, table, xmin, xmax, ymin, ymax);
        cv->updateParametersInContext(context);
        scale *= 2.0;
    }
}

void testReordering() {
    // Create a larger system with a nonbonded force, since that will trigger atom
    // reordering on the GPU.
    
    const int numParticles = 100;
    System system;
    CustomCVForce* cv = new CustomCVForce("2*v2");
    CustomNonbondedForce* v1 = new CustomNonbondedForce("r");
    v1->addPerParticleParameter("a");
    CustomBondForce* v2 = new CustomBondForce("r+1");
    v2->addBond(5, 10);
    cv->addCollectiveVariable("v1", v1);
    cv->addCollectiveVariable("v2", v2);
    cv->setForceGroup(2);
    system.addForce(cv);
    CustomNonbondedForce* nb = new CustomNonbondedForce("r^2");
    nb->addPerParticleParameter("a");
    system.addForce(nb);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions;
    vector<double> params(1);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        params[0] = i%2;
        v1->addParticle(params);
        params[0] = (i < numParticles/2 ? 2.0 : 3.0);
        nb->addParticle(params);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*10);
    }
    
    // Make sure it works correctly.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces, false, 1<<2);
    Vec3 delta = positions[5]-positions[10];
    double r = sqrt(delta.dot(delta));
    ASSERT_EQUAL_TOL(2*(r+1), state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(-delta*2/r, state.getForces()[5], 1e-5);
    ASSERT_EQUAL_VEC(delta*2/r, state.getForces()[10], 1e-5);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testCVs();
        testEnergyParameterDerivatives();
        testTabulatedFunction();
        testReordering();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
