/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
 * This tests the reference implementation of GBVIForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/GBVIForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testSingleParticle() {
    const int log          = 0;
    ReferencePlatform platform;
    System system;
    system.addParticle(2.0);
    LangevinIntegrator integrator(0, 0.1, 0.01);

    GBVIForce* forceField = new GBVIForce();

    double charge         = -1.0;
    double radius         =  0.15;
    double gamma          =  1.0;
    forceField->addParticle(charge, radius, gamma);
    system.addForce(forceField);

    Context context(system, integrator, platform);
    vector<Vec3> positions(1);
    positions[0] = Vec3(0, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Energy);

    double bornRadius     = radius; 
    double eps0           = EPSILON0;
    double tau            = (1.0/forceField->getSoluteDielectric()-1.0/forceField->getSolventDielectric());

    double bornEnergy     = (-charge*charge/(8*PI_M*eps0))*tau/bornRadius;
    double nonpolarEnergy = -gamma*tau*std::pow( radius/bornRadius, 3.0);

    double expectedE      = (bornEnergy+nonpolarEnergy); 
    double obtainedE      = state.getPotentialEnergy(); 
    double diff           = fabs( (obtainedE - expectedE)/expectedE );
    if( log ){
        (void) fprintf( stderr, "testSingleParticle expected=%14.6e obtained=%14.6e diff=%14.6e breakdown:[%14.6e %14.6e]\n",
                        expectedE, obtainedE, diff, bornEnergy, nonpolarEnergy );
    }

    ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);
}

void testEnergyEthane( int applyBornRadiiScaling ) {

    ReferencePlatform platform;
    const int numParticles = 8;
    const int log          = 0;
    System system;
    LangevinIntegrator integrator(0, 0.1, 0.01);

    // harmonic bond

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, C_HBondDistance, 0.0);
    bonds->addBond(2, 1, C_HBondDistance, 0.0);
    bonds->addBond(3, 1, C_HBondDistance, 0.0);

    bonds->addBond(1, 4, C_CBondDistance, 0.0);

    bonds->addBond(5, 4, C_HBondDistance, 0.0);
    bonds->addBond(6, 4, C_HBondDistance, 0.0);
    bonds->addBond(7, 4, C_HBondDistance, 0.0);

    system.addForce(bonds);

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
       C_gamma  = -0.2863;
       H_radius =  0.125;
       H_gamma  =  0.2437;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    if( log ){
       (void) fprintf( stderr, "Applying GB/VI\n" );
    }
    GBVIForce* forceField = new GBVIForce();
    if( applyBornRadiiScaling ){
        forceField->setBornRadiusScalingMethod( GBVIForce::QuinticSpline );
    } else {
        forceField->setBornRadiusScalingMethod( GBVIForce::NoScaling );
    }
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
       forceField->addParticle( H_charge, H_radius, H_gamma);
       nonbonded->addParticle(  H_charge, H_radius, 0.0);
    }
 
    forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma);
    forceField->setParticleParameters( 4, C_charge, C_radius, C_gamma);
 
    nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);
    nonbonded->setParticleParameters(  4, C_charge, C_radius, 0.0);
 
    forceField->addBond( 0, 1, C_HBondDistance );
    forceField->addBond( 2, 1, C_HBondDistance );
    forceField->addBond( 3, 1, C_HBondDistance );
    forceField->addBond( 1, 4, C_CBondDistance );
    forceField->addBond( 5, 4, C_HBondDistance );
    forceField->addBond( 6, 4, C_HBondDistance );
    forceField->addBond( 7, 4, C_HBondDistance );
    
    std::vector<pair<int, int> > bondExceptions;
    std::vector<double> bondDistances;
    
    bondExceptions.push_back(pair<int, int>(0, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(2, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(3, 1)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(1, 4)); 
    bondDistances.push_back( C_CBondDistance );
    
    bondExceptions.push_back(pair<int, int>(5, 4)); 
    bondDistances.push_back( C_HBondDistance );
    
    bondExceptions.push_back(pair<int, int>(6, 4)); 
    bondDistances.push_back( C_HBondDistance );
 
    bondExceptions.push_back(pair<int, int>(7, 4));
    bondDistances.push_back( C_HBondDistance );
 
    nonbonded->createExceptionsFromBonds(bondExceptions, 0.0, 0.0);
 
    system.addForce(forceField);
    system.addForce(nonbonded);

    Context context(system, integrator, platform);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0.5480,    1.7661,    0.0000);
    positions[1] = Vec3(0.7286,    0.8978,    0.6468);
    positions[2] = Vec3(0.4974,    0.0000,    0.0588);
    positions[3] = Vec3(0.0000,    0.9459,    1.4666);
    positions[4] = Vec3(2.1421,    0.8746,    1.1615);
    positions[5] = Vec3(2.3239,    0.0050,    1.8065);
    positions[6] = Vec3(2.8705,    0.8295,    0.3416);
    positions[7] = Vec3(2.3722,    1.7711,    1.7518);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    if( log ){
        (void) fprintf( stderr, "Energy %.4e\n", state.getPotentialEnergy() );
    }
    
    // Take a small step in the direction of the energy gradient.
    
    double norm        = 0.0;
    double forceSum[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f  = state.getForces()[i];
        if( log ){
            (void) fprintf( stderr, "F %d [%14.6e %14.6e %14.6e]\n", i, f[0], f[1], f[2] );
        }
        norm        += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        forceSum[0] += f[0];
        forceSum[1] += f[1];
        forceSum[2] += f[2];
    }
    norm               = std::sqrt(norm);
    if( log ){
        (void) fprintf( stderr, "Fsum [%14.6e %14.6e %14.6e] norm=%14.6e\n", forceSum[0], forceSum[1], forceSum[2], norm );
    }

    const double delta = 1e-4;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    State state2 = context.getState(State::Energy);

    if( log ){
        double deltaE = fabs( state.getPotentialEnergy() - state2.getPotentialEnergy() )/delta;
        double diff   = (deltaE - norm)/norm;
        (void) fprintf( stderr, "Energies %.8e %.8e deltaE=%14.7e %14.7e diff=%14.7e\n", state.getPotentialEnergy(), state2.getPotentialEnergy(), deltaE, norm, diff );
    }

    // See whether the potential energy changed by the expected amount.
    
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)
}

int main() {
    try {
        testSingleParticle();
        testEnergyEthane( 0 );
        testEnergyEthane( 1 );
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
