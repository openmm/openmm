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

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "ReferencePlatform.h"
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/NonbondedSoftcoreForce.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/sfmt/SFMT.h"
#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"
#include "ReferenceFreeEnergyKernelFactory.h"
#include "CudaFreeEnergyKernelFactory.h"

#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

#define PRINT_ON 1

int compareForcesOfTwoStates( int numParticles, State& state1, State& state2, double relativeTolerance, double absoluteTolerance ) {

    int error = 0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f1       = state1.getForces()[i];
        Vec3 f2       = state2.getForces()[i];
        double diff   = (f1[0] - f2[0])*(f1[0] - f2[0]) +
                        (f1[1] - f2[1])*(f1[1] - f2[1]) +
                        (f1[2] - f2[2])*(f1[2] - f2[2]); 
        double denom1 = fabs( f1[0] ) + fabs( f1[1] ) +fabs( f1[2] );
        double denom2 = fabs( f2[0] ) + fabs( f2[1] ) +fabs( f2[2] );
        int        ok = 1;
        if( (denom1 > 0.0 || denom2 > 0.0) && (sqrt( diff )/(denom1+denom2)) > relativeTolerance ){
           error++;
           ok = 0;
        }
#if PRINT_ON == 1
        (void) fprintf( stderr, "F %d [%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e] %s\n", i, 
                        f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], (ok ? "":"XXXXXX") );
#endif
    }

    return error;
}

void testSingleParticle() {
    CudaPlatform platform;
    System system;
    system.addParticle(2.0);
    LangevinIntegrator integrator(0, 0.1, 0.01);

    GBVISoftcoreForce* forceField = new GBVISoftcoreForce;

    double charge         = 1.0;
    double radius         = 0.15;
    double gamma          = 1.0;
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
    double nonpolarEnergy = -0.1*CAL2JOULE*gamma*tau*std::pow( radius/bornRadius, 3.0);

    double expectedE      = (bornEnergy+nonpolarEnergy); 
    double obtainedE      = state.getPotentialEnergy(); 
    double diff           = fabs( obtainedE - expectedE );
#if PRINT_ON == 1
    (void) fprintf( stderr, "testSingleParticle expected=%14.6e obtained=%14.6e diff=%14.6e breakdown:[%14.6e %14.6e]\n",
                    expectedE, obtainedE, diff, bornEnergy, nonpolarEnergy );
#endif
    ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);
}

void testCutoffAndPeriodic() {

    CudaPlatform cuda;

    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);

    LangevinIntegrator integrator(0, 0.1, 0.01);

    GBVISoftcoreForce* gbsa    = new GBVISoftcoreForce();
    NonbondedForce* nonbonded  = new NonbondedForce();

    gbsa->addParticle(-1, 0.15, 1.0);
    nonbonded->addParticle(-1, 1, 0);

    gbsa->addParticle(1, 0.15, 1.0);
    nonbonded->addParticle(1, 1, 0);

    const double cutoffDistance = 3.0;
    const double boxSize = 10.0;

    nonbonded->setCutoffDistance(cutoffDistance);
    system.setPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    system.addForce(gbsa);
    system.addForce(nonbonded);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(2, 0, 0);

    // Calculate the forces for both cutoff and periodic with two different atom positions.

    nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    Context context(system, integrator, cuda);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    context.reinitialize();
    context.setPositions(positions);
    State state2 = context.getState(State::Forces);
    positions[1][0]+= boxSize;
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    context.reinitialize();
    context.setPositions(positions);
    State state3 = context.getState(State::Forces);
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    context.reinitialize();
    context.setPositions(positions);
    State state4 = context.getState(State::Forces);

    // All forces should be identical, exception state3 which should be zero.

#if PRINT_ON == 1
    (void) fprintf( stderr, "testCutoffAndPeriodic\n" );
#endif
    ASSERT_EQUAL_VEC(state1.getForces()[0], state2.getForces()[0], 0.01);
    ASSERT_EQUAL_VEC(state1.getForces()[1], state2.getForces()[1], 0.01);
    ASSERT_EQUAL_VEC(state1.getForces()[0], state4.getForces()[0], 0.01);
    ASSERT_EQUAL_VEC(state1.getForces()[1], state4.getForces()[1], 0.01);
    ASSERT_EQUAL_VEC(state3.getForces()[0], Vec3(0, 0, 0), 0.01);
    ASSERT_EQUAL_VEC(state3.getForces()[1], Vec3(0, 0, 0), 0.01);
}

void testEnergyEthane() {

    std::string methodName = "testEnergyEthane";

#if 0
    CudaPlatform platform;
    CudaFreeEnergyKernelFactory* factory  = new CudaFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
#else

    ReferencePlatform platform;

    ReferenceFreeEnergyKernelFactory* referenceFactoryT  = new ReferenceFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactoryT);
#endif

    ReferencePlatform referencePlatform;

    ReferenceFreeEnergyKernelFactory* referenceFactory  = new ReferenceFreeEnergyKernelFactory();

    referencePlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactory);

    System system;
    const int numParticles = 8;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator(0, 0.1, 0.01);

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
//C_radius =  0.360;
       C_gamma  = -0.2863;
       C_gamma  =  1.0;
       H_radius =  0.125;
//H_radius =  0.25;
       H_gamma  =  0.2437;
       H_gamma  =  1.0;
//H_charge = C_charge = 0.0;
//H_gamma = C_gamma = 0.0;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    // for ethane all Coulomb forces are excluded since all atoms 3 or
    // fewer bonds away from all other atoms -- is this true for H's on
    // difference carbons? -- should be computed in 14 ixn 
  
    int VI = 1;
    if( VI ){

       //double bornRadiusScaleFactorsEven = 0.5;
       double bornRadiusScaleFactorsEven = 1.0;
       //double bornRadiusScaleFactorsOdd  = 0.5;
       double bornRadiusScaleFactorsOdd  = 1.0;
#if PRINT_ON == 1
       (void) fprintf( stderr, "%s: Applying GB/VI\n", methodName.c_str() );
       (void) fprintf( stderr, "C[%14.7e %14.7e %14.7e] H[%14.7e %14.7e %14.7e] scale[%.1f %.1f]\n",
                       C_charge, C_radius, C_gamma, H_charge, H_radius, H_gamma,
                       bornRadiusScaleFactorsEven, bornRadiusScaleFactorsOdd);
#endif

       GBVISoftcoreForce* forceField             = new GBVISoftcoreForce();
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_gamma, (i%2) ? bornRadiusScaleFactorsOdd : bornRadiusScaleFactorsEven);
          nonbonded->addParticle(  H_charge, H_radius, 0.0);
       }

       forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsOdd);
       nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);

       forceField->setParticleParameters( 4, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsEven);
       nonbonded->setParticleParameters(  4, C_charge, C_radius, 0.0);

       forceField->addBond( 0, 1, C_HBondDistance );
       forceField->addBond( 2, 1, C_HBondDistance );
       forceField->addBond( 3, 1, C_HBondDistance );
       forceField->addBond( 1, 4, C_CBondDistance );
       forceField->addBond( 5, 4, C_HBondDistance );
       forceField->addBond( 6, 4, C_HBondDistance );
       forceField->addBond( 7, 4, C_HBondDistance );
   
       std::vector<pair<int, int> > bonds;
       std::vector<double> bondDistances;
   
       bonds.push_back(pair<int, int>(0, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(2, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(3, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(1, 4));
       bondDistances.push_back( C_CBondDistance );
   
       bonds.push_back(pair<int, int>(5, 4));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(6, 4));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(7, 4));
       bondDistances.push_back( C_HBondDistance );
   
       nonbonded->createExceptionsFromBonds(bonds, 0.0, 0.0);

       system.addForce(forceField);

    } else {

#if PRINT_ON == 1
       (void) fprintf( stderr, "testEnergyEthane: Applying GBSA OBC\n" );
#endif
       GBSAOBCForce* forceField = new GBSAOBCForce();
       double H_scale           =  0.85;
       double C_scale           =  0.72;
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_scale );
          nonbonded->addParticle(  H_charge, 1, 0);
       }

       forceField->setParticleParameters(1, C_charge, C_radius, C_scale);
       forceField->setParticleParameters(4, C_charge, C_radius, C_scale);

       nonbonded->setParticleParameters( 1, C_charge, C_radius, 0.0);
       nonbonded->setParticleParameters( 4, C_charge, C_radius, 0.0);

       system.addForce(forceField);
    }

    system.addForce(nonbonded);

    Context referenceContext(system, integrator, referencePlatform);
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
    referenceContext.setPositions(positions);

    State state           = context.getState(State::Forces | State::Energy);
    State referenceState  = referenceContext.getState(State::Forces | State::Energy);

#if PRINT_ON == 1
    (void) fprintf( stderr, "cudaE=%14.7e refE=%14.7e\n", state.getPotentialEnergy(), referenceState.getPotentialEnergy() );
#endif
    
    // Take a small step in the direction of the energy gradient.
    
    if( compareForcesOfTwoStates( numParticles, state, referenceState, 0.001, 0.001 ) ){
       ASSERT_EQUAL_TOL(0.0, 1.0, 0.01)
    }

    double norm        = 0.0;
    double forceSum[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f       = state.getForces()[i];
        norm        += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        forceSum[0] += f[0];
        forceSum[1] += f[1];
        forceSum[2] += f[2];
    }
    norm               = std::sqrt(norm);

#if PRINT_ON == 1
    (void) fprintf( stderr, "Fsum [%14.7e %14.7e %14.7e] norm=%14.7e\n", forceSum[0], forceSum[1], forceSum[2], norm );
#endif

    const double delta = 1e-3;
    double step        = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    State state2 = context.getState(State::Energy);

    double diff  = (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta;
    double off   = fabs( diff - norm );

#if PRINT_ON == 1
    (void) fprintf( stderr, "X Energies %.8e %.8e norms[%14.7e %14.7e] deltaNorms=%14.7e delta=%.2e\n",
                    state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, off, delta );
#endif

    // See whether the potential energy changed by the expected amount.
    
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)
}

void testEnergyEthaneSwitchingFunction() {

    std::string methodName = "testEnergyEthaneSwitchingFunction";

#if 0
    CudaPlatform platform;
    CudaFreeEnergyKernelFactory* factory  = new CudaFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);

#else

    ReferencePlatform platform;

    ReferenceFreeEnergyKernelFactory* referenceFactoryT  = new ReferenceFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactoryT);
#endif

    ReferencePlatform referencePlatform;
    ReferenceFreeEnergyKernelFactory* referenceFactory  = new ReferenceFreeEnergyKernelFactory();

    referencePlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactory);

    System system;
    const int numParticles = 9;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator(0, 0.1, 0.01);

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;

    NonbondedSoftcoreForce* nonbonded = new NonbondedSoftcoreForce();
    nonbonded->setNonbondedMethod(NonbondedSoftcoreForce::NoCutoff);

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
//C_radius =  0.360;
       C_gamma  = -0.2863;
       C_gamma  =  1.0;
       H_radius =  0.125;
//H_radius =  0.25;
       H_gamma  =  0.2437;
       H_gamma  =  1.0;
//H_charge = C_charge = 0.0;
//H_gamma = C_gamma = 0.0;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    // for ethane all Coulomb forces are excluded since all atoms 3 or
    // fewer bonds away from all other atoms -- is this true for H's on
    // difference carbons? -- should be computed in 14 ixn 
  
    int VI = 1;
    if( VI ){

       //double bornRadiusScaleFactorsEven = 0.5;
       double bornRadiusScaleFactorsEven = 1.0;
       //double bornRadiusScaleFactorsOdd  = 0.5;
       double bornRadiusScaleFactorsOdd  = 1.0;
#if PRINT_ON == 1
       (void) fprintf( stderr, "%s: Applying GB/VI\n", methodName.c_str() );
       (void) fprintf( stderr, "C[%14.7e %14.7e %14.7e] H[%14.7e %14.7e %14.7e] scale[%.1f %.1f]\n",
                       C_charge, C_radius, C_gamma, H_charge, H_radius, H_gamma,
                       bornRadiusScaleFactorsEven, bornRadiusScaleFactorsOdd);
#endif

       GBVISoftcoreForce* forceField             = new GBVISoftcoreForce();
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_gamma, (i%2) ? bornRadiusScaleFactorsOdd : bornRadiusScaleFactorsEven);
          nonbonded->addParticle(  H_charge, H_radius, 0.0);
       }

       forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsOdd);
       nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);

       forceField->setParticleParameters( 4, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsEven);
       nonbonded->setParticleParameters(  4, C_charge, C_radius, 0.0);

       forceField->setParticleParameters( 8, C_charge, (C_radius+0.5), C_gamma, bornRadiusScaleFactorsEven);
       nonbonded->setParticleParameters(  8, C_charge, C_radius, 0.0);

       forceField->setBornRadiusScalingMethod( GBVISoftcoreForce::NoScaling );
//       forceField->setBornRadiusScalingMethod( GBVISoftcoreForce::QuinticSpline );

       forceField->addBond( 0, 1, C_HBondDistance );
       forceField->addBond( 2, 1, C_HBondDistance );
       forceField->addBond( 3, 1, C_HBondDistance );
       forceField->addBond( 1, 4, C_CBondDistance );
       forceField->addBond( 5, 4, C_HBondDistance );
       forceField->addBond( 6, 4, C_HBondDistance );
       forceField->addBond( 7, 4, C_HBondDistance );
   
       std::vector<pair<int, int> > bonds;
       std::vector<double> bondDistances;
   
       bonds.push_back(pair<int, int>(0, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(2, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(3, 1));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(1, 4));
       bondDistances.push_back( C_CBondDistance );
   
       bonds.push_back(pair<int, int>(5, 4));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(6, 4));
       bondDistances.push_back( C_HBondDistance );
   
       bonds.push_back(pair<int, int>(7, 4));
       bondDistances.push_back( C_HBondDistance );
   
       nonbonded->createExceptionsFromBonds(bonds, 0.0, 0.0);

       system.addForce(forceField);

    } else {

#if PRINT_ON == 1
       (void) fprintf( stderr, "testEnergyEthane: Applying GBSA OBC\n" );
#endif
       GBSAOBCForce* forceField = new GBSAOBCForce();
       double H_scale           =  0.85;
       double C_scale           =  0.72;
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_scale );
          nonbonded->addParticle(  H_charge, 1, 0);
       }

       forceField->setParticleParameters(1, C_charge, C_radius, C_scale);
       forceField->setParticleParameters(4, C_charge, C_radius, C_scale);

       nonbonded->setParticleParameters( 1, C_charge, C_radius, 0.0);
       nonbonded->setParticleParameters( 4, C_charge, C_radius, 0.0);

       system.addForce(forceField);
    }

    system.addForce(nonbonded);

    Context referenceContext(system, integrator, referencePlatform);
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

    positions[8] = Vec3(2.1421,    0.8746,    2.1615);

    vector<Vec3> originalPositions(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){
       originalPositions[ii][0] = positions[ii][0];
       originalPositions[ii][1] = positions[ii][1];
       originalPositions[ii][2] = positions[ii][2];
    }

    int tries                = 7;
    double positionIncrement = 0.15;
    for( int ii = 0; ii < tries; ii++ ){

       context.setPositions(positions);
       referenceContext.setPositions(positions);
   
       State state           = context.getState(State::Forces | State::Energy);
       State referenceState  = referenceContext.getState(State::Forces | State::Energy);
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "cudaE=%14.7e refE=%14.7e\n", state.getPotentialEnergy(), referenceState.getPotentialEnergy() );
#endif
       
       // Take a small step in the direction of the energy gradient.
       
       if( compareForcesOfTwoStates( numParticles, state, referenceState, 0.001, 0.001 ) ){
          ASSERT_EQUAL_TOL(0.0, 1.0, 0.01)
       }
   
       double norm        = 0.0;
       double forceSum[3] = { 0.0, 0.0, 0.0 };
       for (int i = 0; i < numParticles; ++i) {
           Vec3 f       = state.getForces()[i];
           norm        += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
           forceSum[0] += f[0];
           forceSum[1] += f[1];
           forceSum[2] += f[2];
       }
       norm               = std::sqrt(norm);
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "Fsum [%14.7e %14.7e %14.7e] norm=%14.7e\n", forceSum[0], forceSum[1], forceSum[2], norm );
#endif
   
       const double delta = 1e-03;
       double step        = delta/norm;
       for (int i = 0; i < numParticles; ++i) {
           Vec3 p = positions[i];
           Vec3 f = state.getForces()[i];
           positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
       }
       context.setPositions(positions);
       
       State state2 = context.getState(State::Energy);
   
       double diff  = (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta;
       double off   = fabs( diff - norm )/norm;
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "%2d Energies %.8e %.8e norms[%13.7e %13.7e] deltaNorms=%13.7e delta=%.2e\n",
                       ii, state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, off, delta );
#endif
   
       // See whether the potential energy changed by the expected amount.
       
       ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 1e-3*abs(state.getPotentialEnergy()) );

       if( ii < (tries-1) ){
           for( int jj = 0; jj < numParticles; jj++ ){
              positions[jj][0]  = originalPositions[jj][0];
              positions[jj][1]  = originalPositions[jj][1];
              positions[jj][2]  = originalPositions[jj][2];
           }
       
           positions[8][2] -=  static_cast<double>(ii+1)*0.1;
           positions[8][2] -=  0.001;
           (void) fprintf( stderr, "r48=%14.6e r28=%14.6e r24=%14.6e\n", positions[8][2]-positions[4][2], positions[8][2], positions[4][2] );
       }
#if 0
       int carbonIndex    = 1;
       int hydrogenIndex  = 0;
       while( hydrogenIndex < 8 ){
          Vec3 carbonDelta;
          for( int kk = 0; kk < 3; kk++ ){
             positions[hydrogenIndex][kk] += positionIncrement*(positions[carbonIndex][kk] - positions[hydrogenIndex][kk] );
          }
          double dist = 0.0;
          for( int kk = 0; kk < 3; kk++ ){
             dist += (positions[carbonIndex][kk] - positions[hydrogenIndex][kk] )*(positions[carbonIndex][kk] - positions[hydrogenIndex][kk]);
          }
           (void) fprintf( stderr, "H=%d C=%d r=%14.6e\n", hydrogenIndex, carbonIndex, dist );
          hydrogenIndex++;
          if( hydrogenIndex == carbonIndex ){
             hydrogenIndex++;
          }
          if( carbonIndex == 1 && hydrogenIndex == 4 ){
             carbonIndex    = 4;
             hydrogenIndex  = 5;
          }
       }
#endif

   }
}

void testTwoParticleEnergyEthaneSwitchingFunction() {

    std::string methodName = "testTwoParticleEnergyEthaneSwitchingFunction";

#if 0
    CudaPlatform platform;
    CudaFreeEnergyKernelFactory* factory  = new CudaFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
#else

    ReferencePlatform platform;

    ReferenceFreeEnergyKernelFactory* referenceFactoryT  = new ReferenceFreeEnergyKernelFactory();
    platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactoryT);
    platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactoryT);
#endif

    ReferencePlatform referencePlatform;
    ReferenceFreeEnergyKernelFactory* referenceFactory  = new ReferenceFreeEnergyKernelFactory();

    referencePlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), referenceFactory);
    referencePlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), referenceFactory);

    System system;
    const int numParticles = 3;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator(0, 0.1, 0.01);

    double C_HBondDistance   = 0.1097;
    double C_CBondDistance   = 0.1504;

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;

    int AM1_BCC = 1;
    H_charge    = -0.053;
    C_charge    = -3.0*H_charge;
    if( AM1_BCC ){
       C_radius =  0.180;
//C_radius =  0.360;
       C_gamma  = -0.2863;
       C_gamma  =  1.0;
       H_radius =  0.125;
//H_radius =  0.25;
       H_gamma  =  0.2437;
       H_gamma  =  1.0;
//H_charge = C_charge = 0.0;
//H_gamma = C_gamma = 0.0;
    } else {
       C_radius =  0.215;
       C_gamma  = -1.1087;
       H_radius =  0.150;
       H_gamma  =  0.1237;
    }

    // for ethane all Coulomb forces are excluded since all atoms 3 or
    // fewer bonds away from all other atoms -- is this true for H's on
    // difference carbons? -- should be computed in 14 ixn 
  
    int VI = 1;
    if( VI ){

       //double bornRadiusScaleFactorsEven = 0.5;
       double bornRadiusScaleFactorsEven = 1.0;
       //double bornRadiusScaleFactorsOdd  = 0.5;
       double bornRadiusScaleFactorsOdd  = 1.0;
#if PRINT_ON == 1
       (void) fprintf( stderr, "%s: Applying GB/VI\n", methodName.c_str() );
       (void) fprintf( stderr, "C[%14.7e %14.7e %14.7e] H[%14.7e %14.7e %14.7e] scale[%.1f %.1f]\n",
                       C_charge, C_radius, C_gamma, H_charge, H_radius, H_gamma,
                       bornRadiusScaleFactorsEven, bornRadiusScaleFactorsOdd);
#endif

       GBVISoftcoreForce* forceField             = new GBVISoftcoreForce();
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_gamma, (i%2) ? bornRadiusScaleFactorsOdd : bornRadiusScaleFactorsEven);
          nonbonded->addParticle(  H_charge, H_radius, 0.0);
       }

       forceField->setParticleParameters( 0, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsOdd);
       nonbonded->setParticleParameters(  0, C_charge, C_radius, 0.0);

       forceField->setParticleParameters( 1, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsOdd);
       nonbonded->setParticleParameters(  1, C_charge, C_radius, 0.0);

       forceField->setParticleParameters( 2, C_charge, C_radius, C_gamma, bornRadiusScaleFactorsOdd);
       nonbonded->setParticleParameters(  2, C_charge, C_radius, 0.0);

       system.addForce(forceField);

    } else {

#if PRINT_ON == 1
       (void) fprintf( stderr, "testEnergyEthane: Applying GBSA OBC\n" );
#endif
       GBSAOBCForce* forceField = new GBSAOBCForce();
       double H_scale           =  0.85;
       double C_scale           =  0.72;
       for( int i = 0; i < numParticles; i++ ){
          forceField->addParticle( H_charge, H_radius, H_scale );
          nonbonded->addParticle(  H_charge, 1, 0);
       }

       forceField->setParticleParameters(1, C_charge, C_radius, C_scale);
       forceField->setParticleParameters(4, C_charge, C_radius, C_scale);

       nonbonded->setParticleParameters( 1, C_charge, C_radius, 0.0);
       nonbonded->setParticleParameters( 4, C_charge, C_radius, 0.0);

       system.addForce(forceField);
    }

    system.addForce(nonbonded);

    Context referenceContext(system, integrator, referencePlatform);
    Context context(system, integrator, platform);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3( 0.0000,    0.0000,    0.0000);
    positions[1] = Vec3( 1.0000,    0.0000,    0.0000);
    positions[2] = Vec3(-1.0000,    0.0000,    0.0000);

    vector<Vec3> originalPositions(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){
       originalPositions[ii][0] = positions[ii][0];
       originalPositions[ii][1] = positions[ii][1];
       originalPositions[ii][2] = positions[ii][2];
    }

    int tries                = 11;
    double positionIncrement = 0.15;
    for( int ii = 0; ii < tries; ii++ ){

       context.setPositions(positions);
       referenceContext.setPositions(positions);
   
       State state           = context.getState(State::Forces | State::Energy);
       State referenceState  = referenceContext.getState(State::Forces | State::Energy);
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "cudaE=%14.7e refE=%14.7e\n", state.getPotentialEnergy(), referenceState.getPotentialEnergy() );
#endif
       
       // Take a small step in the direction of the energy gradient.
       
       if( compareForcesOfTwoStates( numParticles, state, referenceState, 0.001, 0.001 ) ){
          ASSERT_EQUAL_TOL(0.0, 1.0, 0.01)
       }
   
       double norm        = 0.0;
       double forceSum[3] = { 0.0, 0.0, 0.0 };
       for (int i = 0; i < numParticles; ++i) {
           Vec3 f       = state.getForces()[i];
           norm        += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
           forceSum[0] += f[0];
           forceSum[1] += f[1];
           forceSum[2] += f[2];
       }
       norm               = std::sqrt(norm);
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "Fsum [%14.7e %14.7e %14.7e] norm=%14.7e\n", forceSum[0], forceSum[1], forceSum[2], norm );
#endif
   
       const double delta = 1e-3;
       double step        = delta/norm;
       for (int i = 0; i < numParticles; ++i) {
           Vec3 p = positions[i];
           Vec3 f = state.getForces()[i];
           positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
       }
       context.setPositions(positions);
       
       State state2 = context.getState(State::Energy);
   
       double diff  = (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta;
       double off   = fabs( diff - norm );
   
#if PRINT_ON == 1
       (void) fprintf( stderr, "X Energies %.8e %.8e norms[%14.7e %14.7e] deltaNorms=%14.7e delta=%.2e\n",
                       state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, off, delta );
#endif
   
       // See whether the potential energy changed by the expected amount.
       
//       ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)

       for( int jj = 0; jj < numParticles; jj++ ){
          positions[jj][0]  = originalPositions[jj][0];
          positions[jj][1]  = originalPositions[jj][1];
          positions[jj][2]  = originalPositions[jj][2];
       }
   
       positions[1][0] -=  static_cast<double>(ii+1)*0.1;
       positions[2][0] +=  static_cast<double>(ii+1)*0.1;
       positions[1][0] -=  0.001;
       positions[2][0] +=  0.001;
       (void) fprintf( stderr, "r12=%14.6e\n", positions[1][0]);
#if 0
       int carbonIndex    = 1;
       int hydrogenIndex  = 0;
       while( hydrogenIndex < 8 ){
          Vec3 carbonDelta;
          for( int kk = 0; kk < 3; kk++ ){
             positions[hydrogenIndex][kk] += positionIncrement*(positions[carbonIndex][kk] - positions[hydrogenIndex][kk] );
          }
          double dist = 0.0;
          for( int kk = 0; kk < 3; kk++ ){
             dist += (positions[carbonIndex][kk] - positions[hydrogenIndex][kk] )*(positions[carbonIndex][kk] - positions[hydrogenIndex][kk]);
          }
           (void) fprintf( stderr, "H=%d C=%d r=%14.6e\n", hydrogenIndex, carbonIndex, dist );
          hydrogenIndex++;
          if( hydrogenIndex == carbonIndex ){
             hydrogenIndex++;
          }
          if( carbonIndex == 1 && hydrogenIndex == 4 ){
             carbonIndex    = 4;
             hydrogenIndex  = 5;
          }
       }
#endif

   }
}

void testEnergyTwoParticle() {

    CudaPlatform platform;
    const int numParticles = 2;
    System system;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator(0, 0.1, 0.01);

    //void HarmonicBondForce::getBondParameters(int index, int& particle1, int& particle2, double& length, double& k)
    double C_HBondDistance   = 3.0;

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;
/*
    H_charge    = -1.0;
    C_charge    =  1.0;

    H_gamma     =  1.0;
    C_gamma     =  1.0;

    H_radius    =  1.0;
    C_radius    =  1.0;
*/ 
    H_charge    = -0.5;
    C_charge    =  0.5;

    H_gamma     =  0.5;
    C_gamma     =  0.5;

    H_radius    =  0.15;
    C_radius    =  0.15;
 
    int VI = 1;
    if( VI ){
       (void) fprintf( stderr, "Applying GB/VI\n" );
       GBVISoftcoreForce* forceField = new GBVISoftcoreForce();
       forceField->addParticle( H_charge, H_radius, H_gamma);
       forceField->addParticle( C_charge, C_radius, C_gamma);
       system.addForce(forceField);
    } else {
       (void) fprintf( stderr, "Applying GBSA OBC\n" );
       GBSAOBCForce* forceField = new GBSAOBCForce();
       forceField->addParticle( H_charge, H_radius, 0.8);
       forceField->addParticle( C_charge, C_radius, 0.8);
       system.addForce(forceField);
    }

    NonbondedForce* nonbonded = new NonbondedForce();
    for (int i = 0; i < numParticles; ++i) {
        double charge = i%2 == 0 ? -1 : 1;
        nonbonded->addParticle( charge, 1, 0);
    }
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);
    system.addForce(nonbonded);

    Context context(system, integrator, platform);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(         0.0000,    0.0000,    0.0000);
    positions[1] = Vec3(C_HBondDistance,    0.0000,    0.0000);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    
    // Take a small step in the direction of the energy gradient.
    
    double norm        = 0.0;
    double forceSum[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f  = state.getForces()[i];
#if PRINT_ON == 1
        (void) fprintf( stderr, "F %d [%14.6e %14.6e %14.6e]\n", i, f[0], f[1], f[2] );
#endif
        norm   += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        forceSum[0] += f[0];
        forceSum[1] += f[1];
        forceSum[2] += f[2];
    }
    norm               = std::sqrt(norm);

#if PRINT_ON == 1
    (void) fprintf( stderr, "Fsum [%14.6e %14.6e %14.6e] norm=%14.6e\n", forceSum[0], forceSum[1], forceSum[2], norm );
#endif

    const double delta = 1e-4;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    State state2 = context.getState(State::Energy);

    double diff = fabs( norm - (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta );
#if PRINT_ON == 1
    (void) fprintf( stderr, "Energies %14.6e %14.6e diff=%14.6e [%14.6e %14.6e]\n",
                    state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta );
#endif

    // See whether the potential energy changed by the expected amount.
    
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)
}


void testEnergyManyParticles( int numParticles ) {

    CudaPlatform platform;
    System system;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator(0, 0.1, 0.01);

    //void HarmonicBondForce::getBondParameters(int index, int& particle1, int& particle2, double& length, double& k)
/* 
    double C_HBondDistance   = 3.0;
    HarmonicBondForce* bonds = new HarmonicBondForce(numParticles-1);
    for( int ii = 1; ii < numParticles; ii++ ){
       bonds->setBondParameters(ii-1, ii-1, ii, C_HBondDistance, 0.0); 
    }
    system.addForce(bonds);
*/

    double C_radius, C_gamma, C_charge, H_radius, H_gamma, H_charge;
/*
    H_charge    = -1.0;
    C_charge    =  1.0;

    H_gamma     =  1.0;
    C_gamma     =  1.0;

    H_radius    =  1.0;
    C_radius    =  1.0;
*/ 
    H_charge    = -0.5;
    C_charge    =  0.5;

    H_gamma     =  0.5;
    C_gamma     =  0.5;

    H_radius    =  0.15;
    C_radius    =  0.15;
 
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::NoCutoff);

    int VI = 0;
    if( VI ){

#if PRINT_ON == 1
       (void) fprintf( stderr, "testEnergyManyParticles: Applying GB/VI\n" );
#endif
       GBVISoftcoreForce* forceField = new GBVISoftcoreForce();
       for( int ii = 0; ii < numParticles; ii++ ){
          forceField->addParticle( H_charge, H_radius, H_gamma);
          nonbonded->addParticle( H_charge, H_radius, 0.0);
       }
       system.addForce(forceField);
    } else {
#if PRINT_ON == 1
       (void) fprintf( stderr, "testEnergyManyParticles: Applying GBSA OBC\n" );
#endif
       GBSAOBCForce* forceField = new GBSAOBCForce();
       for( int ii = 0; ii < numParticles; ii++ ){
          forceField->addParticle(H_charge, H_radius, 0.8);
          nonbonded->addParticle( H_charge, H_radius, 0.0);
       }
       system.addForce(forceField);
    }
    system.addForce(nonbonded);

    Context context(system, integrator, platform);
    vector<Vec3> positions(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){
       positions[ii] = Vec3(  (double) ii,    0.0000,    0.0000);
    }
    context.setPositions(positions);

    //State state = context.getState(State::Forces | State::Energy);
    
   /* 
    // Take a small step in the direction of the energy gradient.
    
    double norm        = 0.0;
    double forceSum[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f  = state.getForces()[i];
        (void) fprintf( stderr, "F %d [%14.6e %14.6e %14.6e]\n", i, f[0], f[1], f[2] );
        norm   += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        forceSum[0] += f[0];
        forceSum[1] += f[1];
        forceSum[2] += f[2];
    }
    norm               = std::sqrt(norm);
    (void) fprintf( stderr, "Fsum [%14.6e %14.6e %14.6e] norm=%14.6e\n", forceSum[0], forceSum[1], forceSum[2], norm );

    const double delta = 1e-4;
    double step = delta/norm;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
    }
    context.setPositions(positions);
    
    State state2 = context.getState(State::Energy);

    double diff = fabs( norm - (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta );
    (void) fprintf( stderr, "Energies %14.6e %14.6e diff=%14.6e [%14.6e %14.6e]\n",
                    state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta );

    // See whether the potential energy changed by the expected amount.
    
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 0.01)
*/
}


void testForce(int numParticles, NonbondedForce::NonbondedMethod method) {
    CudaPlatform cuda;
    ReferencePlatform reference;
    System system;
    LangevinIntegrator integrator(0, 0.1, 0.01);

    GBVISoftcoreForce* gbsa    = new GBVISoftcoreForce();
    NonbondedForce* nonbonded  = new NonbondedForce();

    double radius              = 0.15;
    double gamma               = 0.0;

    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        double charge = i%2 == 0 ? -1 : 1;
        gbsa->addParticle(charge, radius, gamma);
        nonbonded->addParticle(charge, 1, 0);
    }

    nonbonded->setNonbondedMethod(method);
    nonbonded->setCutoffDistance(3.0);

    int grid = (int) floor(0.5+pow(numParticles, 1.0/3.0));
    if (method == NonbondedForce::CutoffPeriodic) {
        double boxSize = (grid+1)*2.0;
        system.setPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    }
    system.addForce(gbsa);
    system.addForce(nonbonded);

    Context context(system, integrator, cuda);
    Context refContext(system, integrator, reference);
    
    // Set random (but uniformly distributed) positions for all the particles.
    
    vector<Vec3> positions(numParticles);
    init_gen_rand(0);
    for (int i = 0; i < grid; i++)
        for (int j = 0; j < grid; j++)
            for (int k = 0; k < grid; k++)
                //positions[i*grid*grid+j*grid+k] = Vec3(i*2.0, j*2.0, k*2.0);
                positions[i*grid*grid+j*grid+k] = Vec3(i*0.5, j*0.5, k*0.5);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = positions[i] + Vec3(0.5*genrand_real2(), 0.5*genrand_real2(), 0.5*genrand_real2());
    context.setPositions(positions);
    refContext.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    State refState = refContext.getState(State::Forces | State::Energy);

    // Make sure the Cuda and Reference platforms agree.

    double norm = 0.0;
    double diff = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        Vec3 f = state.getForces()[i];
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
        Vec3 delta = f-refState.getForces()[i];
        Vec3 g = refState.getForces()[i];
#if PRINT_ON == 1
fprintf( stderr, "FFF %d fcud[%14.6e %14.6e %14.6e] [%14.6e %14.6e %14.6e]\n", i, f[0], f[1], f[2], g[0], g[1], g[2] );
        diff += delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
#endif
    }
    norm = std::sqrt(norm);
    diff = std::sqrt(diff);

#if PRINT_ON == 1
    (void) fprintf( stderr, "F norm%14.6e diff w/ ref=%14.6e\n", norm, diff );
#endif
    ASSERT_EQUAL_TOL(0.0, diff, 0.001*norm); 

    // Take a small step in the direction of the energy gradient.  (This doesn't work with cutoffs, since the energy
    // changes discontinuously at the cutoff distance.)

    if (method == NonbondedForce::NoCutoff)
    {
        const double delta = 1e-2;
        double step = delta/norm;
        for (int i = 0; i < numParticles; ++i) {
            Vec3 p = positions[i];
            Vec3 f = state.getForces()[i];
            positions[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        }
        context.setPositions(positions);

        // See whether the potential energy changed by the expected amount.

        State state2 = context.getState(State::Energy);
        ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state.getPotentialEnergy())/delta, 1e-3*abs(state.getPotentialEnergy()));
    }
}

int main() {
    try {
//        testEnergyEthane();
        testEnergyEthaneSwitchingFunction();
//        testTwoParticleEnergyEthaneSwitchingFunction();
//        testSingleParticle();
//        testCutoffAndPeriodic();
//        testEnergyTwoParticle();
//       for (int i = 2; i < 8; i++) {
//            testForce(i*i*i, NonbondedForce::NoCutoff);
//        }

    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
