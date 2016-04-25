/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the CUDA implementation of CudaAmoebaStretchBendForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

extern "C" void registerAmoebaCudaKernelFactories();

const double TOL = 1e-4;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308
const double DegreesToRadians = PI_M/180.0;

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
static void crossProductVector3(double* vectorX, double* vectorY, double* vectorZ) {

    vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
    vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
    vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

    return;
}

static double dotVector3(double* vectorX, double* vectorY) {
    return vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1] + vectorX[2]*vectorY[2];
}


static void computeAmoebaStretchBendForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaStretchBendForce& amoebaStretchBendForce,
                                          std::vector<Vec3>& forces, double* energy) {

    int particle1, particle2, particle3;
    double abBondLength, cbBondLength, angleStretchBend, kStretchBend, k2StretchBend;

    amoebaStretchBendForce.getStretchBendParameters(bondIndex, particle1, particle2, particle3, abBondLength, cbBondLength, angleStretchBend, kStretchBend, k2StretchBend);
    angleStretchBend *= RADIAN;
    enum { A, B, C, LastAtomIndex };
    enum { AB, CB, CBxAB, ABxP, CBxP, LastDeltaIndex };
 
    // ---------------------------------------------------------------------------------------
 
    // get deltaR between various combinations of the 3 atoms
    // and various intermediate terms
 
    double deltaR[LastDeltaIndex][3];
    double rAB2 = 0.0;
    double rCB2 = 0.0;
    for (int ii = 0; ii < 3; ii++) {
         deltaR[AB][ii]  = positions[particle1][ii] - positions[particle2][ii];
         rAB2           += deltaR[AB][ii]*deltaR[AB][ii];

         deltaR[CB][ii]  = positions[particle3][ii] - positions[particle2][ii];
         rCB2           += deltaR[CB][ii]*deltaR[CB][ii];
    }
    double rAB   = sqrt(rAB2);
    double rCB   = sqrt(rCB2);

    crossProductVector3(deltaR[CB], deltaR[AB], deltaR[CBxAB]);
    double  rP   = dotVector3(deltaR[CBxAB], deltaR[CBxAB]);
            rP   = sqrt(rP);
 
    if (rP <= 0.0) {
       return;
    }
    double dot    = dotVector3(deltaR[CB], deltaR[AB]);
    double cosine = dot/(rAB*rCB);
 
    double angle;
    if (cosine >= 1.0) {
       angle = 0.0;
    }
    else if (cosine <= -1.0) {
       angle = PI_M;
    }
    else {
       angle = RADIAN*acos(cosine);
    }
 
    double termA = -RADIAN/(rAB2*rP);
    double termC =  RADIAN/(rCB2*rP);
 
    // P = CBxAB
 
    crossProductVector3(deltaR[AB], deltaR[CBxAB], deltaR[ABxP]);
    crossProductVector3(deltaR[CB], deltaR[CBxAB], deltaR[CBxP]);
    for (int ii = 0; ii < 3; ii++) {
       deltaR[ABxP][ii] *= termA;
       deltaR[CBxP][ii] *= termC;
    }
 
    double dr1   = rAB - abBondLength;
    double dr2   = rCB - cbBondLength;
 
    termA        = 1.0/rAB;
    termC        = 1.0/rCB;
 
    double drkk = dr1 * kStretchBend + dr2 * k2StretchBend;
 
    // ---------------------------------------------------------------------------------------
 
    // forces
 
    // calculate forces for atoms a, b, c
    // the force for b is then -(a + c)
 
    double subForce[LastAtomIndex][3];
    double dt = angle - angleStretchBend;
    for (int jj = 0; jj < 3; jj++) {
        subForce[A][jj] = kStretchBend*dt*termA*deltaR[AB][jj] + drkk*deltaR[ABxP][jj];
        subForce[C][jj] = k2StretchBend*dt*termC*deltaR[CB][jj] + drkk*deltaR[CBxP][jj];
        subForce[B][jj] = -(subForce[A][jj] + subForce[C][jj]);
    }
 
    // ---------------------------------------------------------------------------------------
 
    // accumulate forces and energy
 
    forces[particle1][0]       -= subForce[0][0];
    forces[particle1][1]       -= subForce[0][1];
    forces[particle1][2]       -= subForce[0][2];

    forces[particle2][0]       -= subForce[1][0];
    forces[particle2][1]       -= subForce[1][1];
    forces[particle2][2]       -= subForce[1][2];

    forces[particle3][0]       -= subForce[2][0];
    forces[particle3][1]       -= subForce[2][1];
    forces[particle3][2]       -= subForce[2][2];

    *energy                    += dt*drkk;
}
 
static void computeAmoebaStretchBendForces(Context& context, AmoebaStretchBendForce& amoebaStretchBendForce,
                                          std::vector<Vec3>& expectedForces, double* expectedEnergy) {

    // get positions and zero forces

    State state                 = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize(positions.size());
    
    for (unsigned int ii = 0; ii < expectedForces.size(); ii++) {
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for (int ii = 0; ii < amoebaStretchBendForce.getNumStretchBends(); ii++) {
        computeAmoebaStretchBendForce(ii, positions, amoebaStretchBendForce, expectedForces, expectedEnergy);
    }
}

void compareWithExpectedForceAndEnergy(Context& context, AmoebaStretchBendForce& amoebaStretchBendForce,
                                        double tolerance, const std::string& idString) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaStretchBendForces(context, amoebaStretchBendForce, expectedForces, &expectedEnergy);
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
}

void testOneStretchBend() {

    System system;
    int numberOfParticles = 3;
    for (int ii = 0; ii < numberOfParticles; ii++) {
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaStretchBendForce* amoebaStretchBendForce = new AmoebaStretchBendForce();

    double abLength         = 0.144800000E+01;
    double cbLength         = 0.101500000E+01;
    double angleStretchBend = 0.108500000E+03*DegreesToRadians;
    //double kStretchBend     = 0.750491578E-01;
    double kStretchBend     = 1.0;

    amoebaStretchBendForce->addStretchBend(0, 1, 2, abLength, cbLength, angleStretchBend, kStretchBend, kStretchBend);

    system.addForce(amoebaStretchBendForce);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(0.262660000E+02,  0.254130000E+02,  0.284200000E+01);
    positions[1] = Vec3(0.273400000E+02,  0.244300000E+02,  0.261400000E+01);
    positions[2] = Vec3(0.269573220E+02,  0.236108860E+02,  0.216376800E+01);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy(context, *amoebaStretchBendForce, TOL, "testOneStretchBend");
    
    // Try changing the stretch-bend parameters and make sure it's still correct.
    
    amoebaStretchBendForce->setStretchBendParameters(0, 0, 1, 2, 1.1*abLength, 1.2*cbLength, 1.3*angleStretchBend, 1.4*kStretchBend, 1.4*kStretchBend);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy(context, *amoebaStretchBendForce, TOL, "testOneStretchBend");
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaStretchBendForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy(context, *amoebaStretchBendForce, TOL, "testOneStretchBend");
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    int numberOfParticles = 3;
    for (int ii = 0; ii < numberOfParticles; ii++)
        system.addParticle(1.0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    AmoebaStretchBendForce* amoebaStretchBendForce = new AmoebaStretchBendForce();
    double abLength         = 0.144800000E+01;
    double cbLength         = 0.101500000E+01;
    double angleStretchBend = 0.108500000E+03*DegreesToRadians;
    double kStretchBend     = 1.0;
    amoebaStretchBendForce->addStretchBend(0, 1, 2, abLength, cbLength, angleStretchBend, kStretchBend, kStretchBend);
    amoebaStretchBendForce->setUsesPeriodicBoundaryConditions(true);
    system.addForce(amoebaStretchBendForce);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    std::vector<Vec3> positions(numberOfParticles);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 1);
    context.setPositions(positions);
    State s1 = context.getState(State::Forces | State::Energy);
    
    // Move one atom to a position that should give identical results.

    positions[2] = Vec3(0, 0, -2);
    context.setPositions(positions);
    State s2 = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), s2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], 1e-5);
}

int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaAmoebaStretchBendForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));
        testOneStretchBend();
        testPeriodic();
    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}



