/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
 * This tests the Reference implementation of AmoebaBondForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <math.h>

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerAmoebaReferenceKernelFactories();

const double TOL = 1e-5;

static void computeAmoebaBondForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaBondForce& AmoebaBondForce,
                                           std::vector<Vec3>& forces, double* energy ) {

    int particle1, particle2;
    double bondLength;
    double quadraticK;
    double cubicK    = AmoebaBondForce.getAmoebaGlobalBondCubic();
    double quarticK  = AmoebaBondForce.getAmoebaGlobalBondQuartic();
    AmoebaBondForce.getBondParameters(bondIndex, particle1, particle2,  bondLength,  quadraticK );

    double deltaR[3];
    double r2 = 0.0;
    for( int ii = 0; ii < 3; ii++ ){
           deltaR[ii]    = positions[particle2][ii] - positions[particle1][ii];
           r2           += deltaR[ii]*deltaR[ii];
    }
    double r                   = sqrt( r2 );

    double bondDelta           = (r - bondLength);
    double bondDelta2          = bondDelta*bondDelta;
    double dEdR                = 1.0 + 1.5*cubicK*bondDelta + 2.0*quarticK*bondDelta2;

           dEdR               *= (r > 0.0) ? (2.0*quadraticK*bondDelta)/r : 0.0;

   forces[particle1][0]       += dEdR*deltaR[0];
   forces[particle1][1]       += dEdR*deltaR[1];
   forces[particle1][2]       += dEdR*deltaR[2];

   forces[particle2][0]       -= dEdR*deltaR[0];
   forces[particle2][1]       -= dEdR*deltaR[1];
   forces[particle2][2]       -= dEdR*deltaR[2];

   *energy                    += (1.0f + cubicK*bondDelta + quarticK*bondDelta2)*quadraticK*bondDelta2;

}

static void computeAmoebaBondForces( Context& context, AmoebaBondForce& AmoebaBondForce,
                                             std::vector<Vec3>& expectedForces, double* expectedEnergy, FILE* log ) {

    // get positions and zero forces

    State state = context.getState(State::Positions);
    std::vector<Vec3> positions = state.getPositions();
    expectedForces.resize( positions.size() );
    
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii][0] = expectedForces[ii][1] = expectedForces[ii][2] = 0.0;
    }

    // calculates forces/energy

    *expectedEnergy = 0.0;
    for( int ii = 0; ii < AmoebaBondForce.getNumBonds(); ii++ ){
        computeAmoebaBondForce(ii, positions, AmoebaBondForce, expectedForces, expectedEnergy );
    }
#ifdef AMOEBA_DEBUG
    if( log ){
        (void) fprintf( log, "computeAmoebaBondForces: expected energy=%15.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
#endif
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, AmoebaBondForce& AmoebaBondForce, double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaBondForces( context, AmoebaBondForce, expectedForces, &expectedEnergy, NULL );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();
#ifdef AMOEBA_DEBUG
    if( log ){
        (void) fprintf( log, "computeAmoebaBondForces: expected energy=%15.7e %15.7e\n", expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }
#endif

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, state.getPotentialEnergy(), tolerance );
}

void testOneBond( FILE* log ) {

    System system;

    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaBondForce* amoebaBondForce = new AmoebaBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaBondForce->setAmoebaGlobalBondCubic( cubicK );
    amoebaBondForce->setAmoebaGlobalBondQuartic( quarticicK );
    amoebaBondForce->addBond(0, 1, bondLength, quadraticK);

    system.addForce(amoebaBondForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));
    std::vector<Vec3> positions(2);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaBondForce, TOL, "testOneBond", log );
}

void testTwoBond( FILE* log ) {

    System system;

    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaBondForce* amoebaBondForce = new AmoebaBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaBondForce->setAmoebaGlobalBondCubic( cubicK );
    amoebaBondForce->setAmoebaGlobalBondQuartic( quarticicK );
    amoebaBondForce->addBond(0, 1, bondLength, quadraticK);
    amoebaBondForce->addBond(1, 2, bondLength, quadraticK);

    system.addForce(amoebaBondForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));
    std::vector<Vec3> positions(3);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaBondForce, TOL, "testTwoBond", log );
    
    // Try changing the bond parameters and make sure it's still correct.
    
    amoebaBondForce->setBondParameters(0, 0, 1, 1.1*bondLength, 1.4*quadraticK);
    amoebaBondForce->setBondParameters(1, 1, 2, 1.2*bondLength, 0.9*quadraticK);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareWithExpectedForceAndEnergy( context, *amoebaBondForce, TOL, "testTwoBond", log );
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaBondForce->updateParametersInContext(context);
    compareWithExpectedForceAndEnergy( context, *amoebaBondForce, TOL, "testTwoBond", log );
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceAmoebaBondForce running test..." << std::endl;
        registerAmoebaReferenceKernelFactories();
        FILE* log = NULL;
        //FILE* log = stderr;

        //testOneBond( log );
        testTwoBond( log );
#ifdef AMOEBA_DEBUG
        if( log && log != stderr )
            (void) fclose( log );
#endif

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
