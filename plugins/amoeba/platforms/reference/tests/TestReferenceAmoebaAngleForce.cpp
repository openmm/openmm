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

/**
 * This tests the Cuda implementation of CudaAmoebaAngleForce.
 */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/CustomAngleForce.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerAmoebaReferenceKernelFactories();

const double TOL = 1e-5;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308
#define RADIAN_TO_DEGREE  57.29577951308
#define DEGREE_TO_RADIAN   0.01745329252
#define RADIAN_INVERSE     0.01745329252

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

static void getPrefactorsGivenAngleCosine(double cosine, double idealAngle, double quadraticK, double cubicK,
                                          double quarticK, double penticK, double sexticK,
                                          double* dEdR, double* energyTerm) {

    double angle;
    if (cosine >= 1.0) {
        angle = 0.0f;
    }
    else if (cosine <= -1.0) {
        angle = RADIAN*PI_M;
    }
    else {
        angle = RADIAN*acos(cosine);
    }

    double deltaIdeal         = angle - idealAngle;
    double deltaIdeal2        = deltaIdeal*deltaIdeal;
    double deltaIdeal3        = deltaIdeal*deltaIdeal2;
    double deltaIdeal4        = deltaIdeal2*deltaIdeal2;
 
    // deltaIdeal = r - r_0
 
    *dEdR        = (2.0                        +
                    3.0*cubicK*  deltaIdeal    +
                    4.0*quarticK*deltaIdeal2   +
                    5.0*penticK* deltaIdeal3   +
                    6.0*sexticK* deltaIdeal4   );
 
    *dEdR       *= RADIAN*quadraticK*deltaIdeal;
 

    *energyTerm  = 1.0f + cubicK* deltaIdeal    +
                          quarticK*deltaIdeal2   +
                          penticK* deltaIdeal3   +
                          sexticK* deltaIdeal4;
    *energyTerm *= quadraticK*deltaIdeal2;

    return;
}

static void computeAmoebaAngleForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaAngleForce& AmoebaAngleForce,
                                             std::vector<Vec3>& forces, double* energy) {

    int particle1, particle2, particle3;
    double idealAngle;
    double quadraticK;
    AmoebaAngleForce.getAngleParameters(bondIndex, particle1, particle2, particle3, idealAngle, quadraticK);

    double cubicK         = AmoebaAngleForce.getAmoebaGlobalAngleCubic();
    double quarticK       = AmoebaAngleForce.getAmoebaGlobalAngleQuartic();
    double penticK        = AmoebaAngleForce.getAmoebaGlobalAnglePentic();
    double sexticK        = AmoebaAngleForce.getAmoebaGlobalAngleSextic();

    double deltaR[2][3];
    double r2_0 = 0.0;
    double r2_1 = 0.0;
    for (int ii = 0; ii < 3; ii++) {

           deltaR[0][ii]    = positions[particle1][ii] - positions[particle2][ii];
           r2_0            += deltaR[0][ii]*deltaR[0][ii];

           deltaR[1][ii]    = positions[particle3][ii] - positions[particle2][ii];
           r2_1            += deltaR[1][ii]*deltaR[1][ii];

    }

    double pVector[3];
    crossProductVector3(deltaR[0], deltaR[1], pVector);
    double rp      = sqrt(pVector[0]*pVector[0] + pVector[1]*pVector[1] + pVector[2]*pVector[2]);
    if (rp < 1.0e-06) {
       rp = 1.0e-06;
    }
    double dot    = deltaR[0][0]*deltaR[1][0] + deltaR[0][1]*deltaR[1][1] + deltaR[0][2]*deltaR[1][2];
    double cosine = dot/sqrt(r2_0*r2_1);

    double dEdR;
    double energyTerm;
    getPrefactorsGivenAngleCosine(cosine, idealAngle, quadraticK, cubicK,
                                  quarticK, penticK, sexticK, &dEdR, &energyTerm);

    double termA  = -dEdR/(r2_0*rp);
    double termC  =  dEdR/(r2_1*rp);

    double deltaCrossP[3][3];
    crossProductVector3(deltaR[0], pVector, deltaCrossP[0]);
    crossProductVector3(deltaR[1], pVector, deltaCrossP[2]);
    for (int ii = 0; ii < 3; ii++) {
        deltaCrossP[0][ii] *= termA;
        deltaCrossP[2][ii] *= termC;
        deltaCrossP[1][ii]  = -1.0*(deltaCrossP[0][ii] + deltaCrossP[2][ii]);
    }

    forces[particle1][0]       += deltaCrossP[0][0];
    forces[particle1][1]       += deltaCrossP[0][1];
    forces[particle1][2]       += deltaCrossP[0][2];

    forces[particle2][0]       += deltaCrossP[1][0];
    forces[particle2][1]       += deltaCrossP[1][1];
    forces[particle2][2]       += deltaCrossP[1][2];

    forces[particle3][0]       += deltaCrossP[2][0];
    forces[particle3][1]       += deltaCrossP[2][1];
    forces[particle3][2]       += deltaCrossP[2][2];

    *energy                    += energyTerm;
}

static void computeAmoebaAngleForces(Context& context, AmoebaAngleForce& AmoebaAngleForce,
                                             std::vector<Vec3>& expectedForces, double* expectedEnergy) {

    // get positions and zero forces

    State state = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize(positions.size());
    
    for (unsigned int ii = 0; ii < expectedForces.size(); ii++) {
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for (int ii = 0; ii < AmoebaAngleForce.getNumAngles(); ii++) {
        computeAmoebaAngleForce(ii, positions, AmoebaAngleForce, expectedForces, expectedEnergy);
    }
}

void compareWithExpectedForceAndEnergy(Context& context, AmoebaAngleForce& AmoebaAngleForce,
                                       double tolerance, const std::string& idString) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaAngleForces(context, AmoebaAngleForce, expectedForces, &expectedEnergy);
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
}

void testOneAngle() {

    System system;
    int numberOfParticles = 3;
    for (int ii = 0; ii < numberOfParticles; ii++) {
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaAngleForce* amoebaAngleForce = new AmoebaAngleForce();

    double angle      = 100.0;
    double quadraticK = 1.0;
    double cubicK     = 1.0e-01;
    double quarticK   = 1.0e-02;
    double penticK    = 1.0e-03;
    double sexticK    = 1.0e-04;
    amoebaAngleForce->addAngle(0, 1, 2, angle, quadraticK);

    amoebaAngleForce->setAmoebaGlobalAngleCubic(cubicK);
    amoebaAngleForce->setAmoebaGlobalAngleQuartic(quarticK);
    amoebaAngleForce->setAmoebaGlobalAnglePentic(penticK);
    amoebaAngleForce->setAmoebaGlobalAngleSextic(sexticK);

    system.addForce(amoebaAngleForce);
    ASSERT(!amoebaAngleForce->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, Platform::getPlatformByName("Reference"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy(context, *amoebaAngleForce, TOL, "testOneAngle");
    
    // Try changing the angle parameters and make sure it's still correct.
    
    amoebaAngleForce->setAngleParameters(0, 0, 1, 2, 1.1*angle, 1.4*quadraticK);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy(context, *amoebaAngleForce, TOL, "testOneAngle");
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaAngleForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy(context, *amoebaAngleForce, TOL, "testOneAngle");
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions, then compare to an identical custom force.
    
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    int numParticles = 3;
    for (int ii = 0; ii < numParticles; ii++)
        system.addParticle(1.0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    AmoebaAngleForce* amoebaAngleForce = new AmoebaAngleForce();
    double angle      = 100.0;
    double quadraticK = 1.0;
    double cubicK     = 1.0e-01;
    double quarticK   = 1.0e-02;
    double penticK    = 1.0e-03;
    double sexticK    = 1.0e-04;
    amoebaAngleForce->addAngle(0, 1, 2, angle, quadraticK);
    amoebaAngleForce->setAmoebaGlobalAngleCubic(cubicK);
    amoebaAngleForce->setAmoebaGlobalAngleQuartic(quarticK);
    amoebaAngleForce->setAmoebaGlobalAnglePentic(penticK);
    amoebaAngleForce->setAmoebaGlobalAngleSextic(sexticK);
    amoebaAngleForce->setUsesPeriodicBoundaryConditions(true);
    system.addForce(amoebaAngleForce);
    CustomAngleForce* customForce = new CustomAngleForce("k2*delta^2 + k3*delta^3 + k4*delta^4 + k5*delta^5 + k6*delta^6; delta=theta-theta0");
    customForce->addGlobalParameter("theta0", angle*M_PI/180);
    customForce->addGlobalParameter("k2", quadraticK*pow(180/M_PI, 2.0));
    customForce->addGlobalParameter("k3", cubicK*pow(180/M_PI, 3.0));
    customForce->addGlobalParameter("k4", quarticK*pow(180/M_PI, 4.0));
    customForce->addGlobalParameter("k5", penticK*pow(180/M_PI, 5.0));
    customForce->addGlobalParameter("k6", sexticK*pow(180/M_PI, 6.0));
    customForce->addAngle(0, 1, 2);
    customForce->setUsesPeriodicBoundaryConditions(true);
    customForce->setForceGroup(1);
    system.addForce(customForce);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));

    std::vector<Vec3> positions(numParticles);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 2);

    context.setPositions(positions);
    State s1 = context.getState(State::Forces | State::Energy, true, 1);
    State s2 = context.getState(State::Forces | State::Energy, true, 2);
    ASSERT_EQUAL_TOL(s2.getPotentialEnergy(), s1.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(s2.getForces()[i], s1.getForces()[i], 1e-5);
}

int main(int numberOfArguments, char* argv[]) {

    try {
        std::cout << "TestReferenceAmoebaAngleForce running test..." << std::endl;
        registerAmoebaReferenceKernelFactories();
        testOneAngle();
        testPeriodic();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    //std::cout << "PASS - Test succeeded." << std::endl;
    std::cout << "Done" << std::endl;
    return 0;
}
