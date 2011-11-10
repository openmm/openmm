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
#include "ReferencePlatform.h"
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/CustomGBForce.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/NonbondedForce.h"
#include "openmm/NonbondedSoftcoreForce.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../../../../../libraries/sfmt/include/sfmt/SFMT.h"
#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"
#include "ReferenceFreeEnergyKernelFactory.h"

#define USE_SOFTCORE

#include "TestReferenceSoftcoreForce.h"

#define OBC_FLAG  1
#define GBVI_FLAG 2
#define IMPLICIT_SOLVENT GBVI_FLAG

#include <iomanip>

using namespace OpenMM;
using namespace std;

#define PRINT_ON 0

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

    System system;
    system.addParticle(2.0);
    LangevinIntegrator integrator(0, 0.1, 0.01);

    GBVISoftcoreForce* forceField = new GBVISoftcoreForce();

    double charge         = 1.0;
    double radius         = 0.15;
    double gamma          = 1.0;
    forceField->addParticle(charge, radius, gamma);
    system.addForce(forceField);

    NonbondedSoftcoreForce* nonbonded = new NonbondedSoftcoreForce();
    nonbonded->setNonbondedMethod(NonbondedSoftcoreForce::NoCutoff);
    nonbonded->addParticle( charge, 1.0, 0.0);
    system.addForce(nonbonded);

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
    double diff           = fabs( obtainedE - expectedE );
#if PRINT_ON == 1
    (void) fprintf( stderr, "testSingleParticle expected=%14.6e obtained=%14.6e diff=%14.6e breakdown:[%14.6e %14.6e]\n",
                    expectedE, obtainedE, diff, bornEnergy, nonpolarEnergy );
#endif
    ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);
}

void testEnergyEthaneSwitchingFunction( int useSwitchingFunction ) {

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
    const int numParticles = 8;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }
    LangevinIntegrator integrator1(0, 0.1, 0.01);
    LangevinIntegrator integrator2(0, 0.1, 0.01);

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
  
    double bornRadiusScaleFactorsEven = 0.5;
    //double bornRadiusScaleFactorsEven = 1.0;
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

//       forceField->setParticleParameters( 8, C_charge, (C_radius+0.5), C_gamma, bornRadiusScaleFactorsEven);
//       nonbonded->setParticleParameters(  8, C_charge, C_radius, 0.0);

    if( useSwitchingFunction ){
       forceField->setBornRadiusScalingMethod( GBVISoftcoreForce::QuinticSpline );
    } else {
       forceField->setBornRadiusScalingMethod( GBVISoftcoreForce::NoScaling );
    }

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

    system.addForce(nonbonded);

    Context referenceContext(system, integrator1, referencePlatform);
    Context context(system, integrator2, platform);
    
    vector<Vec3> positions(numParticles);
    positions[0] = Vec3(0.5480,    1.7661,    0.0000);
    positions[1] = Vec3(0.7286,    0.8978,    0.6468);
    positions[2] = Vec3(0.4974,    0.0000,    0.0588);
    positions[3] = Vec3(0.0000,    0.9459,    1.4666);
    positions[4] = Vec3(2.1421,    0.8746,    1.1615);
    positions[5] = Vec3(2.3239,    0.0050,    1.8065);
    positions[6] = Vec3(2.8705,    0.8295,    0.3416);
    positions[7] = Vec3(2.3722,    1.7711,    1.7518);

    //positions[8] = Vec3(2.1421,    0.8746,    2.1615);

    vector<Vec3> originalPositions(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){
       originalPositions[ii][0] = positions[ii][0];
       originalPositions[ii][1] = positions[ii][1];
       originalPositions[ii][2] = positions[ii][2];
    }

    int tries                = 1;
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
       
//           positions[8][2] -=  static_cast<double>(ii+1)*0.1;
//           positions[8][2] -=  0.001;
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

// computes the scaled radii based on covalent info and atomic radii

#ifdef USE_SOFTCORE
static void findScaledRadii( const GBVISoftcoreForce& gbviForce, std::vector<double> & scaledRadii,  FILE* log ) {
#else
static void findScaledRadii( const GBVIForce& gbviForce, std::vector<double> & scaledRadii,  FILE* log) {
#endif

    int     numberOfParticles = gbviForce.getNumParticles();
    int numberOfBonds         = gbviForce.getNumBonds();
    FILE* errorLog            = log ? log : stderr;
    
    // load 1-2 atom pairs along w/ bond distance using HarmonicBondForce & constraints
    // numberOfBonds < 1, indicating they were not set by the user
    
    if( numberOfBonds < 1 && numberOfParticles > 1 ){
        (void) fprintf( errorLog, "Warning: no covalent bonds set for GB/VI force!\n" );
    }
    
    std::vector< std::vector<int> > bondIndices( numberOfBonds );
    std::vector<double> bondLengths( numberOfBonds );
    std::vector<double> radii( numberOfParticles);

    scaledRadii.resize(numberOfParticles);

    for (int i = 0; i < numberOfParticles; i++) {
        double charge, radius, gamma;
#ifdef USE_SOFTCORE
        double lambda;
        gbviForce.getParticleParameters(i, charge, radius, gamma, lambda);
#else
        gbviForce.getParticleParameters(i, charge, radius, gamma);
#endif
        radii[i]       = radius;
        scaledRadii[i] = radius;
    }

    for (int i = 0; i < numberOfBonds; i++) {
        int particle1, particle2;
        double bondLength;
        gbviForce.getBondParameters(i, particle1, particle2, bondLength);
        if (particle1 < 0 || particle1 >= gbviForce.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= gbviForce.getNumParticles()) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: Illegal particle index: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (bondLength < 0 ) {
            std::stringstream msg;
            msg << "GBVISoftcoreForce: negative bondlength: ";
            msg << bondLength;
            throw OpenMMException(msg.str());
        }
        bondIndices[i].push_back( particle1 );
        bondIndices[i].push_back( particle2 );
        bondLengths[i] = bondLength;
    }


    // load 1-2 indicies for each atom 

    std::vector<std::vector<int> > bonded12(numberOfParticles);

    for (int i = 0; i < (int) bondIndices.size(); ++i) {
        bonded12[bondIndices[i][0]].push_back(i);
        bonded12[bondIndices[i][1]].push_back(i);
    }

    int errors = 0;

    // compute scaled radii (Eq. 5 of Labute paper [JCC 29 p. 1693-1698 2008])

    for (int j = 0; j < (int) bonded12.size(); ++j){

        double radiusJ = radii[j];
        double scaledRadiusJ;
        if(  bonded12[j].size() == 0 ){
            if( numberOfParticles > 1 ){
                (void) fprintf( errorLog, "Warning GBVIForceImpl::findScaledRadii atom %d has no covalent bonds; using atomic radius=%.3f.\n", j, radiusJ );
            }
            scaledRadiusJ = radiusJ;
//             errors++;
        } else {

            double rJ2    = radiusJ*radiusJ;
    
            // loop over bonded neighbors of atom j, applying Eq. 5 in Labute

            scaledRadiusJ = 0.0;
            for (int i = 0; i < (int) bonded12[j].size(); ++i){
    
               int index            = bonded12[j][i];
               int bondedAtomIndex  = (j == bondIndices[index][0]) ? bondIndices[index][1] : bondIndices[index][0];
              
               double radiusI       = radii[bondedAtomIndex];
               double rI2           = radiusI*radiusI;
    
               double a_ij          = (radiusI - bondLengths[index]);
                      a_ij         *= a_ij;
                      a_ij          = (rJ2 - a_ij)/(2.0*bondLengths[index]);
    
               double a_ji          = radiusJ - bondLengths[index];
                      a_ji         *= a_ji;
                      a_ji          = (rI2 - a_ji)/(2.0*bondLengths[index]);
    
               scaledRadiusJ       += a_ij*a_ij*(3.0*radiusI - a_ij) + a_ji*a_ji*( 3.0*radiusJ - a_ji );
            }
    
            scaledRadiusJ  = (radiusJ*radiusJ*radiusJ) - 0.125*scaledRadiusJ; 
            if( scaledRadiusJ > 0.0 ){
                scaledRadiusJ  = 0.95*pow( scaledRadiusJ, (1.0/3.0) );
            } else {
                scaledRadiusJ  = 0.0;
            }
        }
        if( log ){
            //(void) fprintf( stderr, "scaledRadii %d %12.4f\n", j, scaledRadiusJ );
        }
        scaledRadii[j] = scaledRadiusJ;

    }

    // abort if errors

    if( errors ){
        throw OpenMMException("findScaledRadii errors -- aborting");
    }

    if( log ){
        (void) fprintf( log, "                  R              q          gamma   scaled radii no. bnds\n" );
        double totalQ = 0.0;
        for( int i = 0; i < (int) scaledRadii.size(); i++ ){
    
            double charge;
            double gamma;
            double radiusI;
         
            gbviForce.getParticleParameters(i, charge, radiusI, gamma); 
            totalQ += charge;
            (void) fprintf( log, "%4d %14.5e %14.5e %14.5e %14.5e %d\n", i, radiusI, charge, gamma, scaledRadii[i], (int) bonded12[i].size() );
        }
        (void) fprintf( log, "Total charge=%e\n", totalQ );
        (void) fflush( log );
    }

    return;

}

// load parameters from gbviForce to customGbviForce
// findScaledRadii() is called to calculate the scaled radii (S)
// S is derived quantity in GBVIForce, not a parameter is the case here

#ifdef USE_SOFTCORE
static void loadGbviParameters( const GBVISoftcoreForce& gbviSoftcoreForce, CustomGBForce& customGbviForce, FILE* log ) {
    vector<double> params(5);
    double lambda;
#else
static void loadGbviParameters( const GBVIForce& gbviSoftcoreForce, CustomGBForce& customGbviForce,  FILE* log ) {
    vector<double> params(4);
#endif

    int numParticles = gbviSoftcoreForce.getNumParticles();

    // get scaled radii

    std::vector<double> scaledRadii( numParticles );
    findScaledRadii( gbviSoftcoreForce, scaledRadii, log);

    for( int ii = 0; ii < numParticles; ii++) {
        double charge, radius, gamma;
#ifdef USE_SOFTCORE
        gbviSoftcoreForce.getParticleParameters(ii, charge, radius, gamma, lambda);
        params[4] = lambda;
#else
        gbviSoftcoreForce.getParticleParameters(ii, charge, radius, gamma);
#endif
        params[0] = charge;
        params[1] = radius;
        params[2] = scaledRadii[ii];
        params[3] = gamma;
        customGbviForce.addParticle(params);
    }

}

// create custom GB/VI force

#ifdef USE_SOFTCORE
static void createCustomGBVI( CustomGBForce& customGbviForce, const GBVISoftcoreForce& gbviSoftcoreForce, FILE* log ){
#else
static void createCustomGBVI( CustomGBForce& customGbviForce, const GBVIForce& gbviSoftcoreForce, FILE* log ){
#endif

    customGbviForce.setCutoffDistance( gbviSoftcoreForce.getCutoffDistance() );

    customGbviForce.addPerParticleParameter("q");
    customGbviForce.addPerParticleParameter("radius");
    customGbviForce.addPerParticleParameter("scaleFactor"); // derived in GBVIForceImpl implmentation, but parameter here
    customGbviForce.addPerParticleParameter("gamma");

    customGbviForce.addGlobalParameter("solventDielectric", gbviSoftcoreForce.getSolventDielectric() );
    customGbviForce.addGlobalParameter("soluteDielectric",  gbviSoftcoreForce.getSoluteDielectric() );

#ifdef USE_SOFTCORE
    customGbviForce.addPerParticleParameter("lambda");
    customGbviForce.addComputedValue("V", "                lambda2*(uL - lL + factor3/(radius1*radius1*radius1));"
                                      "uL                   = 1.5*x2uI*(0.25*rI-0.33333*xuI+0.125*(r2-S2)*rI*x2uI);"
                                      "lL                   = 1.5*x2lI*(0.25*rI-0.33333*xlI+0.125*(r2-S2)*rI*x2lI);"
                                      "x2lI                 = 1.0/(xl*xl);"
                                      "xlI                  = 1.0/(xl);"
                                      "xuI                  = 1.0/(xu);"
                                      "x2uI                 = 1.0/(xu*xu);"
                                      "xu                   = (r+scaleFactor2);"
                                      "rI                   = 1.0/(r);"
                                      "r2                   = (r*r);"
                                      "xl                   = factor1*lMax + factor2*xuu + factor3*(r-scaleFactor2);"
                                      "xuu                  = (r+scaleFactor2);"
                                      "S2                   = (scaleFactor2*scaleFactor2);"
                                      "factor1              = step(r-absRadiusScaleDiff);"
                                      "absRadiusScaleDiff   = abs(radiusScaleDiff);"
                                      "radiusScaleDiff      = (radius1-scaleFactor2);"
                                      "factor2              = step(radius1-scaleFactor2-r);"
                                      "factor3              = step(scaleFactor2-radius1-r);"
                                      "lMax                 = max(radius1,r-scaleFactor2);"
                                      , CustomGBForce::ParticlePairNoExclusions);

    customGbviForce.addComputedValue("B", "(1.0/(radius*radius*radius)-V)^(-0.33333333)", CustomGBForce::SingleParticle);

    // nonpolar term + polar self energy

    customGbviForce.addEnergyTerm("(-138.935485*0.5*((1.0/soluteDielectric)-(1.0/solventDielectric))*q^2/B)-((1.0/soluteDielectric)-(1.0/solventDielectric))*((gamma*(radius/B)^3))", CustomGBForce::SingleParticle);

    // polar pair energy

    customGbviForce.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                                   "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePairNoExclusions);

#else

    customGbviForce.addComputedValue("V", "                uL - lL + factor3/(radius1*radius1*radius1);"
                                      "uL                   = 1.5*x2uI*(0.25*rI-0.33333*xuI+0.125*(r2-S2)*rI*x2uI);"
                                      "lL                   = 1.5*x2lI*(0.25*rI-0.33333*xlI+0.125*(r2-S2)*rI*x2lI);"
                                      "x2lI                 = 1.0/(xl*xl);"
                                      "xlI                  = 1.0/(xl);"
                                      "xuI                  = 1.0/(xu);"
                                      "x2uI                 = 1.0/(xu*xu);"
                                      "xu                   = (r+scaleFactor2);"
                                      "rI                   = 1.0/(r);"
                                      "r2                   = (r*r);"
                                      "xl                   = factor1*lMax + factor2*xuu + factor3*(r-scaleFactor2);"
                                      "xuu                  = (r+scaleFactor2);"
                                      "S2                   = (scaleFactor2*scaleFactor2);"
                                      "factor1              = step(r-absRadiusScaleDiff);"
                                      "absRadiusScaleDiff   = abs(radiusScaleDiff);"
                                      "radiusScaleDiff      = (radius1-scaleFactor2);"
                                      "factor2              = step(radius1-scaleFactor2-r);"
                                      "factor3              = step(scaleFactor2-radius1-r);"
                                      "lMax                 = max(radius1,r-scaleFactor2);"
                                      , CustomGBForce::ParticlePairNoExclusions);

    customGbviForce.addComputedValue("B", "(1.0/(radius*radius*radius)-V)^(-0.33333333)", CustomGBForce::SingleParticle);

    // nonpolar term + polar self energy

    customGbviForce.addEnergyTerm("(-138.935485*0.5*((1.0/soluteDielectric)-(1.0/solventDielectric))*q^2/B)-((1.0/soluteDielectric)-(1.0/solventDielectric))*((gamma*(radius/B)^3))", CustomGBForce::SingleParticle);

    // polar pair energy

    customGbviForce.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                                   "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", CustomGBForce::ParticlePairNoExclusions);

#endif
    // load energies

    loadGbviParameters( gbviSoftcoreForce, customGbviForce, log );

    return;
}

void testGBVISoftcore( MapStringToDouble inputArgumentMap, FILE* log ){

    double lambda1                       = 1.0;
    double lambda2                       = 1.0;
    int nonbondedMethod                  = 0;
    int numMolecules                     = 1;
    int numParticlesPerMolecule          = 2;
    int useQuinticSpline                 = 1;
    int applyAssert                      = 1;
    int positionPlacementMethod          = 0;
    int serialize                        = 0;
    double boxSize                       = 10.0;
    double relativeTolerance             = 1.0e-04;

    setDoubleFromMapStringToDouble( inputArgumentMap, "lambda1",                      lambda1 );
    setDoubleFromMapStringToDouble( inputArgumentMap, "lambda2",                      lambda2 );
    setDoubleFromMapStringToDouble( inputArgumentMap, "boxSize",                      boxSize );
    double cutoffDistance                = boxSize*0.4;;
    setDoubleFromMapStringToDouble( inputArgumentMap, "cutoffDistance",               cutoffDistance);
    setDoubleFromMapStringToDouble( inputArgumentMap, "relativeTolerance",            relativeTolerance );

    setIntFromMapStringToDouble(    inputArgumentMap, "positionPlacementMethod",      positionPlacementMethod ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "nonbondedMethod",              nonbondedMethod );
    setIntFromMapStringToDouble(    inputArgumentMap, "numMolecules",                 numMolecules );
    setIntFromMapStringToDouble(    inputArgumentMap, "numParticlesPerMolecule",      numParticlesPerMolecule );
    setIntFromMapStringToDouble(    inputArgumentMap, "serialize",                    serialize );

    int numParticles                     = numMolecules*numParticlesPerMolecule;
    int includeGbvi                      = 1;
    double reactionFieldDielectric       = 80.0;

    if( log ){
        double particleDensity = static_cast<double>(numParticles)/(boxSize*boxSize*boxSize);
        double particleCube    = pow( particleDensity, (-1.0/3.0) );

        (void) fprintf( log, "\n--------------------------------------------------------------------------------------\n" );
        (void) fprintf( log, "Input arguments\n" );
        (void) fflush( log );
        (void) fprintf( log, "    includeGbvi                 %d\n", includeGbvi );
        (void) fprintf( log, "    nonbondedMethod             %d\n", nonbondedMethod );
        (void) fprintf( log, "    numParticles                %d\n", numParticles );
        (void) fprintf( log, "    numMolecules                %d\n", numMolecules );
        (void) fprintf( log, "    numParticlesPerMolecule     %d\n", numParticlesPerMolecule );
        (void) fprintf( log, "    useQuinticSpline            %d\n", useQuinticSpline );
        (void) fprintf( log, "    positionPlacementMethod     %d\n", positionPlacementMethod);

#ifdef USE_SOFTCORE
        (void) fprintf( log, "    lambda1                     %8.3f\n", lambda1 );
        (void) fprintf( log, "    lambda2                     %8.3f\n", lambda2 );
#endif
        (void) fprintf( log, "    boxSize                     %8.3f\n", boxSize );
        (void) fprintf( log, "    cutoffDistance              %8.3f\n", cutoffDistance );
        (void) fprintf( log, "    reactionFieldDielectric     %8.3f\n", reactionFieldDielectric );
        (void) fprintf( log, "    relativeTolerance           %8.1e\n", relativeTolerance );
        (void) fprintf( log, "    particleDensity             %8.2e\n", particleDensity );
        (void) fprintf( log, "    particleCube                %8.2e\n", particleCube );
    }

    // set positions

    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    PositionGenerator positionGenerator( numMolecules, numParticlesPerMolecule, boxSize );
    if( log ){
        positionGenerator.setLog( log );
    }
    if( positionPlacementMethod == 1 ){
        positionGenerator.setPositions( PositionGenerator::SimpleGrid, sfmt, positions );
    } else {
        positionGenerator.setBondDistance( 0.3 );
        positionGenerator.setPositions( PositionGenerator::Random, sfmt, positions );
    }

    // create GBSAGBVISoftcoreForce and populate w/ parameters

#ifdef USE_SOFTCORE
    GBVISoftcoreForce* gbviSoftcoreForce = new GBVISoftcoreForce();
#else
    GBVIForce* gbviSoftcoreForce         = new GBVIForce();
#endif
    gbviSoftcoreForce->setSolventDielectric( 78.3 );
    //gbviSoftcoreForce.setSolventDielectric( 1.0e+10 );
    //gbviSoftcoreForce.setSolventDielectric( 1.0 );
    gbviSoftcoreForce->setSoluteDielectric( 1.0 );
    gbviSoftcoreForce->setCutoffDistance( cutoffDistance );

    const int numberOfParameters             = 4;

    const int ChargeIndex                    = 0;
    const int SigmaIndex                     = 1;
    const int GammaIndex                     = 2;
    const int LambdaIndex                    = 3;

    std::vector<double> parameterLowerBound( numberOfParameters, 0.0 );

    double fixedCharge                       = 1.0;
    parameterLowerBound[ChargeIndex]         = fixedCharge;  // charge
    parameterLowerBound[SigmaIndex]          = 0.1;          // sigma
    parameterLowerBound[GammaIndex]          = 0.1;          // gamma
    parameterLowerBound[LambdaIndex]         = lambda1;      // lambda

    std::vector<double> parameterUpperBound( parameterLowerBound );
    parameterUpperBound[ChargeIndex]         = fixedCharge;  // charge
    parameterUpperBound[SigmaIndex]          = 0.3;          // sigma
    parameterUpperBound[GammaIndex]          = 40.0;         // gamma

    std::vector<double> parameters( numberOfParameters );
    double charge                            = fixedCharge;

    for( int ii = 0; ii < numMolecules; ii++) {

        charge       *= -1.0;

        double lambda =  ii < (numMolecules/2) ? lambda1 : lambda2;
        randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );

#ifdef USE_SOFTCORE
        gbviSoftcoreForce->addParticle( charge,  parameters[SigmaIndex],  parameters[GammaIndex],  lambda );
#else
        gbviSoftcoreForce->addParticle( charge,  parameters[SigmaIndex],  parameters[GammaIndex] );
#endif

        int baseParticleIndex = ii*numParticlesPerMolecule;
        for( int jj = 1; jj < numParticlesPerMolecule; jj++) {

            // alternate charges

            charge *= -1.0;

            randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );

#ifdef USE_SOFTCORE
            gbviSoftcoreForce->addParticle( charge,  parameters[SigmaIndex],  parameters[GammaIndex],  lambda );
#else
            gbviSoftcoreForce->addParticle( charge,  parameters[SigmaIndex],  parameters[GammaIndex] );
#endif

#if IMPLICIT_SOLVENT == GBVI_FLAG
            double bondDistance  = positionGenerator.getDistance( baseParticleIndex, baseParticleIndex+jj, positions );
            gbviSoftcoreForce->addBond( baseParticleIndex, baseParticleIndex+jj,  bondDistance );
#endif
        }

        // alternate charge if numParticlesPerMolecule is odd

        if( (numParticlesPerMolecule % 2) ){
            charge *= -1.0;
        }
    }

    // Create system

    System standardSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    standardSystem.addForce(gbviSoftcoreForce);

    // copy system and forces

    System systemCopy;
    copySystem( standardSystem, systemCopy );
    CustomGBForce* customGbviForce = new CustomGBForce();
    createCustomGBVI( *customGbviForce, *gbviSoftcoreForce, log );
    systemCopy.addForce(customGbviForce);

    // perform comparison

    std::stringstream idString;
    idString << "Nb " << nonbondedMethod << " l2 " << std::fixed << setprecision(2) << lambda2;
    runSystemComparisonTest( standardSystem, systemCopy, "Reference", "Reference", positions, inputArgumentMap, idString.str(), log );

    // serialize

    std::stringstream baseFileName;
    if( serialize ){
        baseFileName  << "_N"     << positions.size();
        baseFileName  << "_Nb"    << nonbondedMethod;
        serializeSystemAndPositions( standardSystem, positions, baseFileName.str(), log);
    }

}

int main() {
    try {
        //FILE* log = stderr;
        FILE* log = NULL;
        //testSingleParticle();
        //testEnergyEthaneSwitchingFunction( 0 );
        //testEnergyEthaneSwitchingFunction( 1 );

        MapStringToDouble inputArgumentMap;
        inputArgumentMap["lambda2"]                         = 0.5;
        inputArgumentMap["nonbondedMethod"]                 = 0;
        inputArgumentMap["numMolecules"]                    = 10;
        inputArgumentMap["boxSize"]                         = 5.0;
        inputArgumentMap["positionPlacementMethod"]         = 0;
        inputArgumentMap["cutoffDistance"]                  = 0.3*inputArgumentMap["boxSize"];
        //inputArgumentMap["cutoffDistance"]                  = 1.0;
        inputArgumentMap["relativeTolerance"]               = 5.0e-04;
        inputArgumentMap["serialize"]                       = 1;
        //inputArgumentMap["numParticlesPerMolecule"]         = 2;

        testGBVISoftcore( inputArgumentMap, log );
    } catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
