/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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
 * This tests the Ewald summation method cuda implementation of NonbondedForce.
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testPMEWater() {
    Platform& platform = Platform::getPlatformByName("Cuda");
    System system;
    system.addParticle(16);
    system.addParticle(1);
    system.addParticle(1);
    VerletIntegrator integrator(0.01);
    AmoebaMultipoleForce* mp = new AmoebaMultipoleForce();
    mp->setNonbondedMethod(AmoebaMultipoleForce::PME);
    vector<double> dipole(3, 0.0);
    dipole[2] = 7.556121361e-2;
    vector<double> quadrupole(9, 0.0);
    quadrupole[0] = 3.540307211e-2;
    quadrupole[4] = -3.902570771e-2;
    quadrupole[8] = 3.622635596e-3;
    double damp = 9.707801995e-01*sqrt(0.1);
    double polarity = 0.837*0.001;
    mp->addParticle(-0.51966, dipole, quadrupole, 1, 1, 2, -1, 0.39, damp, polarity);
    dipole[0] = -2.042094848e-2;
    dipole[2] = -3.078753000e-2;
    quadrupole[0] = -3.428482490e-3;
    quadrupole[2] = -1.894859639e-4;
    quadrupole[4] = -1.002408752e-2;
    quadrupole[6] = -1.894859639e-4;
    quadrupole[8] = 1.345257001e-2;
    damp          = 8.897068742e-01*sqrt(0.1);
    polarity      = 0.496*0.001;
    mp->addParticle(0.25983, dipole, quadrupole, 0, 0, 2, -1, 0.39, damp, polarity);
    mp->addParticle(0.25983, dipole, quadrupole, 0, 0, 1, -1, 0.39, damp, polarity);
    mp->setCutoffDistance(1.0);

    std::vector<int> intVector;
    intVector.push_back( 0 );
    intVector.push_back( 1 );
    intVector.push_back( 2 );
    mp->setCovalentMap( 0, AmoebaMultipoleForce::PolarizationCovalent11, intVector );
    mp->setCovalentMap( 1, AmoebaMultipoleForce::PolarizationCovalent11, intVector );
    mp->setCovalentMap( 2, AmoebaMultipoleForce::PolarizationCovalent11, intVector );

    intVector.resize(0); intVector.push_back( 1 ); intVector.push_back( 2 );
    mp->setCovalentMap( 0, AmoebaMultipoleForce::Covalent12, intVector );

    intVector.resize(0); intVector.push_back( 0 ); intVector.push_back( 2 );
    mp->setCovalentMap( 1, AmoebaMultipoleForce::Covalent12, intVector );

    intVector.resize(0); intVector.push_back( 0 ); intVector.push_back( 1 );
    mp->setCovalentMap( 2, AmoebaMultipoleForce::Covalent12, intVector );

    mp->setEwaldErrorTolerance(TOL);
    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0, 2, 0), Vec3(0, 0, 2));
    system.addForce(mp);
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    const double angle = 109.47*M_PI/180;
    const double dOH = 0.1;
    positions[0] = Vec3();
    positions[1] = Vec3(dOH, 0, 0);
    positions[2] = Vec3(dOH*std::cos(angle), dOH*std::sin(angle), 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces = state.getForces();

#ifdef AMOEBA_DEBUG
    (void) fprintf( stderr, "PME forces\n" );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( stderr, "%6u [%14.7e %14.7e %14.7e]\n", ii, 
                            forces[ii][0], forces[ii][1], forces[ii][2] );
    }
    (void) fflush( stderr );
#endif

}

int main() {
    try {
        Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
        testPMEWater();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

