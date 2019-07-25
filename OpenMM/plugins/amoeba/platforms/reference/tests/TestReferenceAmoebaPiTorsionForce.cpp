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
 * This tests the Reference implementation of ReferenceAmoebaPiTorsionForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
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


static void computeAmoebaPiTorsionForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaPiTorsionForce& amoebaPiTorsionForce,
                                        std::vector<Vec3>& forces, double* energy) {

    int particle1, particle2, particle3, particle4, particle5, particle6;
    double kTorsion;

    amoebaPiTorsionForce.getPiTorsionParameters(bondIndex, particle1, particle2, particle3, particle4,  particle5, particle6, kTorsion);

    enum { AD, BD, EC, FC, P, Q, CP, DC, QD, T, U, TU, DP, QC, dT, dU, dP, dQ, dC1, dC2, dD1, dD2, LastDeltaIndex };
    double deltaR[LastDeltaIndex][3];

    enum { A, B, C, D, E, F, LastAtomIndex };
    double d[LastAtomIndex][3];
 
    for (int ii = 0; ii < 3; ii++) {
         deltaR[AD][ii] = positions[particle1][ii] - positions[particle4][ii];
         deltaR[BD][ii] = positions[particle2][ii] - positions[particle4][ii];
         deltaR[EC][ii] = positions[particle5][ii] - positions[particle3][ii];
         deltaR[FC][ii] = positions[particle6][ii] - positions[particle3][ii];
    }
 
    crossProductVector3(deltaR[AD], deltaR[BD], deltaR[P]);
    crossProductVector3(deltaR[EC], deltaR[FC], deltaR[Q]);
    for (int ii = 0; ii < 3; ii++) {
        deltaR[CP][ii]  = -deltaR[P][ii];
        deltaR[DC][ii]  =  positions[particle4][ii] - positions[particle3][ii];
        deltaR[QD][ii]  =  deltaR[Q][ii];
  
        deltaR[P][ii]  += positions[particle3][ii];
        deltaR[Q][ii]  += positions[particle4][ii];
    }
    crossProductVector3(deltaR[CP], deltaR[DC], deltaR[T]);
    crossProductVector3(deltaR[DC], deltaR[QD], deltaR[U]);
    crossProductVector3(deltaR[T],  deltaR[U],  deltaR[TU]);
 
    double rT2  = dotVector3(deltaR[T], deltaR[T]);
    double rU2  = dotVector3(deltaR[U], deltaR[U]);
    double rTrU = sqrt(rT2*rU2);
    if (rTrU <= 0.0) {
        return;
    }
 
    double rDC     = dotVector3(deltaR[DC], deltaR[DC]);
         rDC       = sqrt(rDC);
  
    double cosine  = dotVector3(deltaR[T], deltaR[U]);
         cosine   /= rTrU;
 
    double sine   = dotVector3(deltaR[DC], deltaR[TU]);
         sine    /= (rDC*rTrU);
 
    double cosine2 = cosine*cosine - sine*sine;
    double sine2   = 2.0*cosine*sine;
 
    double phi2    = 1.0 - cosine2;
    double dphi2   = 2.0*sine2;
 
    double dedphi = kTorsion*dphi2; 
 
    for (int ii = 0; ii < 3; ii++) {
        deltaR[DP][ii] = positions[particle4][ii] -    deltaR[P][ii];
        deltaR[QC][ii] = deltaR[Q][ii]    - positions[particle3][ii];
    }
 
    double factorT =  dedphi/(rDC*rT2);
    double factorU = -dedphi/(rDC*rU2);
 
    crossProductVector3(deltaR[T], deltaR[DC], deltaR[dT]);
    crossProductVector3(deltaR[U], deltaR[DC], deltaR[dU]);
    for (int ii = 0; ii < 3; ii++) {
        deltaR[dT][ii] *= factorT;
        deltaR[dU][ii] *= factorU;
    }
 
    crossProductVector3(deltaR[dT], deltaR[DC], deltaR[dP] );
    crossProductVector3(deltaR[dU], deltaR[DC], deltaR[dQ] );
 
    crossProductVector3(deltaR[DP], deltaR[dT], deltaR[dC1]);
    crossProductVector3(deltaR[dU], deltaR[QD], deltaR[dC2]);
 
    crossProductVector3(deltaR[dT], deltaR[CP], deltaR[dD1]);
    crossProductVector3(deltaR[QC], deltaR[dU], deltaR[dD2]);
 
    crossProductVector3(deltaR[BD], deltaR[dP], d[A]       );
    crossProductVector3(deltaR[dP], deltaR[AD], d[B]       );
 
    crossProductVector3(deltaR[FC], deltaR[dQ], d[E]       );
    crossProductVector3(deltaR[dQ], deltaR[EC], d[F]       );
 
    for (int ii = 0; ii < 3; ii++) {
        d[C][ii] = deltaR[dC1][ii] + deltaR[dC2][ii] + deltaR[dP][ii] - d[E][ii] - d[F][ii];
        d[D][ii] = deltaR[dD1][ii] + deltaR[dD2][ii] + deltaR[dQ][ii] - d[A][ii] - d[B][ii];
    }
 
    // ---------------------------------------------------------------------------------------
 
    // accumulate forces and energy
 
    forces[particle1][0]       -= d[0][0];
    forces[particle1][1]       -= d[0][1];
    forces[particle1][2]       -= d[0][2];

    forces[particle2][0]       -= d[1][0];
    forces[particle2][1]       -= d[1][1];
    forces[particle2][2]       -= d[1][2];

    forces[particle3][0]       -= d[2][0];
    forces[particle3][1]       -= d[2][1];
    forces[particle3][2]       -= d[2][2];

    forces[particle4][0]       -= d[3][0];
    forces[particle4][1]       -= d[3][1];
    forces[particle4][2]       -= d[3][2];

    forces[particle5][0]       -= d[4][0];
    forces[particle5][1]       -= d[4][1];
    forces[particle5][2]       -= d[4][2];

    forces[particle6][0]       -= d[5][0];
    forces[particle6][1]       -= d[5][1];
    forces[particle6][2]       -= d[5][2];

    *energy                    += kTorsion*phi2;
 
    return;
}
 
static void computeAmoebaPiTorsionForces(Context& context, AmoebaPiTorsionForce& amoebaPiTorsionForce,
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
    for (int ii = 0; ii < amoebaPiTorsionForce.getNumPiTorsions(); ii++) {
        computeAmoebaPiTorsionForce(ii, positions, amoebaPiTorsionForce, expectedForces, expectedEnergy);
    }
}

void compareWithExpectedForceAndEnergy(Context& context, AmoebaPiTorsionForce& amoebaPiTorsionForce,
                                       double tolerance, const std::string& idString) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaPiTorsionForces(context, amoebaPiTorsionForce, expectedForces, &expectedEnergy);
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), tolerance);
}

void testOnePiTorsion() {

    System system;
    int numberOfParticles = 6;
    for (int ii = 0; ii < numberOfParticles; ii++) {
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaPiTorsionForce* amoebaPiTorsionForce = new AmoebaPiTorsionForce();

    double kTorsion = 6.85;
    amoebaPiTorsionForce->addPiTorsion(0, 1, 2, 3, 4, 5, kTorsion);

    system.addForce(amoebaPiTorsionForce);
    ASSERT(!amoebaPiTorsionForce->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Context context(system, integrator, Platform::getPlatformByName("Reference"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0] = Vec3(0.262660000E+02,  0.254130000E+02,  0.284200000E+01);
    positions[1] = Vec3(0.278860000E+02,  0.264630000E+02,  0.426300000E+01);
    positions[2] = Vec3(0.269130000E+02,  0.266390000E+02,  0.353100000E+01);

    positions[3] = Vec3(0.245568230E+02,  0.250215290E+02,  0.796852800E+01);
    positions[4] = Vec3(0.261000000E+02,  0.292530000E+02,  0.520200000E+01);
    positions[5] = Vec3(0.254124630E+02,  0.234691880E+02,  0.773335400E+01);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy(context, *amoebaPiTorsionForce, TOL, "testOnePiTorsion");
    
    // Try changing the torsion parameters and make sure it's still correct.
    
    amoebaPiTorsionForce->setPiTorsionParameters(0, 0, 1, 2, 3, 4, 5, 1.2*kTorsion);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy(context, *amoebaPiTorsionForce, TOL, "testOnePiTorsion");
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaPiTorsionForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy(context, *amoebaPiTorsionForce, TOL, "testOnePiTorsion");
}

void testPeriodic() {
    // Create a force that uses periodic boundary conditions.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3));
    int numberOfParticles = 6;
    for (int ii = 0; ii < numberOfParticles; ii++)
        system.addParticle(1.0);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    AmoebaPiTorsionForce* amoebaPiTorsionForce = new AmoebaPiTorsionForce();
    double kTorsion = 6.85;
    amoebaPiTorsionForce->addPiTorsion(0, 1, 2, 3, 4, 5, kTorsion);
    amoebaPiTorsionForce->setUsesPeriodicBoundaryConditions(true);
    system.addForce(amoebaPiTorsionForce);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));
    std::vector<Vec3> positions(numberOfParticles);
    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(0, 0, 0.5);
    positions[3] = Vec3(0.4, 0.4, 0.4);
    positions[4] = Vec3(1, 0, 1);
    positions[5] = Vec3(1, 1, 0);
    context.setPositions(positions);
    State s1 = context.getState(State::Forces | State::Energy);
    
    // Move one atom to a position that should give identical results.

    positions[0] = Vec3(0, -2, 0);
    context.setPositions(positions);
    State s2 = context.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), s2.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], 1e-5);
}

int main(int numberOfArguments, char* argv[]) {

    try {
        std::cout << "TestReferenceAmoebaPiTorsionForce running test..." << std::endl;
        registerAmoebaReferenceKernelFactories();
        testOnePiTorsion();
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
