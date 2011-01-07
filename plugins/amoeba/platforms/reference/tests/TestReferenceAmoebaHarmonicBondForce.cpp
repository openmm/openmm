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
 * This tests the Reference implementation of HarmonicBondForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace OpenMM;

const double TOL = 1e-5;

static void computeAmoebaHarmonicBondForce(int bondIndex,  std::vector<Vec3>& positions, AmoebaHarmonicBondForce& amoebaHarmonicBondForce,
                                           std::vector<Vec3>& forces, double* energy ) {

    int particle1, particle2;
    double bondLength;
    double quadraticK;
    double cubicK    = amoebaHarmonicBondForce.getAmoebaGlobalHarmonicBondCubic();
    double quarticK  = amoebaHarmonicBondForce.getAmoebaGlobalHarmonicBondQuartic();
    amoebaHarmonicBondForce.getBondParameters(bondIndex, particle1, particle2,  bondLength,  quadraticK );

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

static void computeAmoebaHarmonicBondForces( Context& context, AmoebaHarmonicBondForce& amoebaHarmonicBondForce,
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
    for( int ii = 0; ii < amoebaHarmonicBondForce.getNumBonds(); ii++ ){
        computeAmoebaHarmonicBondForce(ii, positions, amoebaHarmonicBondForce, expectedForces, expectedEnergy );
    }

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicBondForces: expected energy=%15.7e\n", *expectedEnergy );
        for( unsigned int ii = 0; ii < positions.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]\n", ii, expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2] );
        }
        (void) fflush( log );
    }
    return;

}

void compareWithExpectedForceAndEnergy( Context& context, AmoebaHarmonicBondForce& amoebaHarmonicBondForce, double tolerance, const std::string& idString, FILE* log) {

    std::vector<Vec3> expectedForces;
    double expectedEnergy;
    computeAmoebaHarmonicBondForces( context, amoebaHarmonicBondForce, expectedForces, &expectedEnergy, NULL );
   
    State state                      = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces   = state.getForces();

    if( log ){
        (void) fprintf( log, "computeAmoebaHarmonicBondForces: expected energy=%15.7e %15.7e\n", expectedEnergy, state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n", ii,
                            expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2], forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }

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

    AmoebaHarmonicBondForce* amoebaHarmonicBondForce = new AmoebaHarmonicBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondCubic( cubicK );
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondQuartic( quarticicK );
    amoebaHarmonicBondForce->addBond(0, 1, bondLength, quadraticK);

    system.addForce(amoebaHarmonicBondForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));
    std::vector<Vec3> positions(2);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaHarmonicBondForce, TOL, "testOneBond", log );
}

void testTwoBond( FILE* log ) {

    System system;

    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    AmoebaHarmonicBondForce* amoebaHarmonicBondForce = new AmoebaHarmonicBondForce();

    double bondLength = 1.5;
    double quadraticK = 1.0;
    double cubicK     = 2.0;
    double quarticicK = 3.0;
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondCubic( cubicK );
    amoebaHarmonicBondForce->setAmoebaGlobalHarmonicBondQuartic( quarticicK );
    amoebaHarmonicBondForce->addBond(0, 1, bondLength, quadraticK);
    amoebaHarmonicBondForce->addBond(1, 2, bondLength, quadraticK);

    system.addForce(amoebaHarmonicBondForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));
    std::vector<Vec3> positions(3);

    positions[0] = Vec3(0, 1, 0);
    positions[1] = Vec3(0, 0, 0);
    positions[2] = Vec3(1, 0, 1);

    context.setPositions(positions);
    compareWithExpectedForceAndEnergy( context, *amoebaHarmonicBondForce, TOL, "testTwoBond", log );
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceAmoebaHarmonicBondForce running test..." << std::endl;
        Platform::loadPluginsFromDirectory( Platform::getDefaultPluginsDirectory() );
        FILE* log = NULL;
        //FILE* log = stderr;

        //testOneBond( log );
        testTwoBond( log );

        if( log && log != stderr )
            (void) fclose( log );

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "PASS - Test succeeded." << std::endl;
    return 0;
}
