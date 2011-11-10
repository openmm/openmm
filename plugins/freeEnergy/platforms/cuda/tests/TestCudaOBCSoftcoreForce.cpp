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
 * This tests the Cuda implementations of GBSAOBCSoftcoreForce
 */

#include "openmm/System.h"
#include "../../../tests/AssertionUtilities.h"

#include "sfmt/SFMT.h"
#include "openmm/Context.h"

#include "openmm/NonbondedSoftcoreForce.h"
#include "openmm/GBSAOBCSoftcoreForce.h"
 #include "../src/SimTKUtilities/SimTKOpenMMRealType.h"

#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

typedef std::vector<double> DoubleVector;
typedef DoubleVector::iterator DoubleVectorI;
typedef DoubleVector::const_iterator DoubleVectorCI;
typedef std::vector<DoubleVector> DoubleVectorVector;

typedef std::pair<int, double> IntDoublePair;
typedef std::vector<IntDoublePair> IntDoublePairVector;
typedef IntDoublePairVector::iterator IntDoublePairVectorI;
typedef IntDoublePairVector::const_iterator IntDoublePairVectorCI;

extern "C" void registerFreeEnergyCudaKernelFactories();

using namespace OpenMM;
using namespace std;

void testSingleParticle( FILE* log ) {

    System system;
    system.addParticle(2.0);
    VerletIntegrator integrator(0.01);

    GBSAOBCSoftcoreForce* forceField = new GBSAOBCSoftcoreForce();

    double charge         = 0.5;
    double radius         = 0.15;
    double scaleFactor    = 1.0;
    forceField->addParticle(charge, radius, scaleFactor);
    system.addForce(forceField);

    NonbondedSoftcoreForce* nonbonded = new NonbondedSoftcoreForce();
    nonbonded->setNonbondedMethod(NonbondedSoftcoreForce::NoCutoff);
    nonbonded->addParticle( charge, 1.0, 0.0);
    system.addForce(nonbonded);

    Context context(system, integrator, Platform::getPlatformByName( "Cuda") );
    vector<Vec3> positions(1);
    positions[0] = Vec3(0, 0, 0);
    context.setPositions(positions);
    State state = context.getState(State::Energy);

    double bornRadius     = radius-0.009; 
    double eps0           = EPSILON0;
    double tau            = (1.0/forceField->getSoluteDielectric()-1.0/forceField->getSolventDielectric());

    double bornEnergy     = (-0.5*0.5/(8*PI_M*eps0))*(1.0/forceField->getSoluteDielectric()-1.0/forceField->getSolventDielectric())/bornRadius;
    double extendedRadius = bornRadius+0.14; // probe radius
    double nonpolarEnergy = CAL2JOULE*PI_M*0.0216*(10*extendedRadius)*(10*extendedRadius)*std::pow(0.15/bornRadius, 6.0); // Where did this formula come from?  Just copied it from CpuImplicitSolvent.cpp

    double expectedE      = (bornEnergy+nonpolarEnergy); 
    double obtainedE      = state.getPotentialEnergy(); 
    double diff           = fabs( obtainedE - expectedE );
    if( log ){
        (void) fprintf( log, "testSingleParticle expected=%14.6e obtained=%14.6e diff=%14.6e breakdown:[%14.6e %14.6e]\n",
                        expectedE, obtainedE, diff, bornEnergy, nonpolarEnergy );
    }
    ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);
}

/** 
 * Predicate for sorting <int,double> pair
 *
 * @param d1 first  IntDoublePair to compare
 * @param d2 second IntDoublePair to compare
 *
 */

bool TestIntDoublePair( const IntDoublePair& d1, const IntDoublePair& d2 ){
   return d1.second < d2.second;
}

/**---------------------------------------------------------------------------------------
 *
 * Get relative difference between two forces
 * 
 * @param  f1                   force1
 * @param  f2                   force2
 * @param  forceNorm1           output norm of force1
 * @param  forceNorm2           output norm of force2
 * @param  relativeDiff         output relative difference between force norms
 * @param  log                  if set, output forces
 *
 *
   --------------------------------------------------------------------------------------- */

static void getForceRelativeDifference( const Vec3& f1, const Vec3& f2, double& forceNorm1, double& forceNorm2,
                                        double& relativeDiff, FILE* log ) {

    double diff     = (f1[0] - f2[0])*(f1[0] - f2[0]) +
                      (f1[1] - f2[1])*(f1[1] - f2[1]) +
                      (f1[2] - f2[2])*(f1[2] - f2[2]);

    forceNorm1      = sqrt( f1[0]*f1[0] + f1[1]*f1[1] + f1[2]*f1[2] );
    forceNorm2      = sqrt( f2[0]*f2[0] + f2[1]*f2[1] + f2[2]*f2[2] );

    if( forceNorm1 > 0.0 || forceNorm2 > 0.0 ){
        relativeDiff = 2.0*sqrt( diff )/(forceNorm1+forceNorm2);
    } else {
        relativeDiff = 0.0;
    }

    return;
}

/**---------------------------------------------------------------------------------------
 *
 * Compare forces from two states
 * 
 * @param  state1               state1
 * @param  state2               state2
 * @param  relativeTolerance    relative tolerance
 * @param  log                  if set, output forces
 *
 * @return number of entries with relative difference > tolerance 
 *
   --------------------------------------------------------------------------------------- */

int compareForcesOfTwoStates( State& state1, State& state2, double relativeTolerance,
                              DoubleVector& stats, FILE* log ) {

    int error                             = 0;
    vector<Vec3> force1                   = state1.getForces();
    vector<Vec3> force2                   = state2.getForces();
    double averageRelativeDifference      = 0.0;
    double count                          = 0.0;

    DoubleVector medians1( force1.size() );
    DoubleVector medians2( force1.size() );

    IntDoublePairVector relativeDifferences;

    for( unsigned int ii = 0; ii < force1.size(); ii++ ){

        double forceNorm1;
        double forceNorm2;
        double relativeDiff;
        getForceRelativeDifference( force1[ii], force2[ii], forceNorm1, forceNorm2, relativeDiff, log );

        medians1[ii]               = forceNorm1;
        medians2[ii]               = forceNorm2;
 
        relativeDifferences.push_back( IntDoublePair(ii, relativeDiff ) );
        averageRelativeDifference += relativeDiff;
        count                     += 1.0;

        if( relativeDiff > relativeTolerance ){
           error++;
        }

        if( log ){
            (void) fprintf( log, "F %6u %15.7e [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] %15.7e %15.7e %s\n", static_cast<unsigned int>(ii), 
                            relativeDiff, force1[ii][0], force1[ii][1], force1[ii][2], force2[ii][0], force2[ii][1], force2[ii][2],
                            forceNorm1, forceNorm2, (relativeDiff < relativeTolerance ? "":"XXXXXX") );
        }
    }

    // sort relative differences

    std::sort( relativeDifferences.begin(), relativeDifferences.end(), TestIntDoublePair );

    if( log ){
        (void) fprintf( log, "\nEntries w/ largest relative differences.\n" );
        for( unsigned int ii = relativeDifferences.size()-1; ii >= relativeDifferences.size()-10 && ii >= 0; ii-- ){
            double forceNorm1;
            double forceNorm2;
            double relativeDiff;
            int index = relativeDifferences[ii].first;
            getForceRelativeDifference( force1[index], force2[index], forceNorm1, forceNorm2, relativeDiff, log );
            (void) fprintf( log, "Fs %6u %15.7e [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] %15.7e %15.7e %s\n",
                            static_cast<unsigned int>(index), relativeDiff, 
                            force1[index][0], force1[index][1], force1[index][2],
                            force2[index][0], force2[index][1], force2[index][2], 
                            forceNorm1, forceNorm2, (relativeDiff < relativeTolerance ? "":"XXXXXX") );
        }
    }

    if( count > 0.0 ){
        averageRelativeDifference /= count;
    }

    std::sort( medians1.begin(), medians1.end() );
    std::sort( medians2.begin(), medians2.end() );
    double median1 = medians1[medians1.size()/2];
    double median2 = medians2[medians2.size()/2];

    stats.resize( 4 );
    stats[0]           = averageRelativeDifference;
    IntDoublePair pair = relativeDifferences[relativeDifferences.size()-1];
    stats[1]           = pair.second;
    stats[2]           = static_cast<double>(pair.first);
    stats[3]           = median1 < median2 ? median1 : median2;
    
    return error;
}

void testEnergyEthaneSwitchingFunction( int useSwitchingFunction, FILE* log ) {

    std::string methodName = "testEnergyEthaneSwitchingFunction";

    System system;
    const int numParticles = 8;
    for( int i = 0; i < numParticles; i++ ){
       system.addParticle(1.0);
    }

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
    if( log ){
        (void) fprintf( log, "%s: Applying GB/VI\n", methodName.c_str() );
        (void) fprintf( log, "C[%14.7e %14.7e %14.7e] H[%14.7e %14.7e %14.7e] scale[%.1f %.1f]\n",
                    C_charge, C_radius, C_gamma, H_charge, H_radius, H_gamma,
                    bornRadiusScaleFactorsEven, bornRadiusScaleFactorsOdd);
    }

    GBSAOBCSoftcoreForce* forceField             = new GBSAOBCSoftcoreForce();
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

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context referenceContext(system, integrator1,  Platform::getPlatformByName( "Reference") );
    Context context(system, integrator2,  Platform::getPlatformByName( "Cuda") );
    
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
   
       
       if( log ){
           (void) fprintf( log, "cudaE=%14.7e refE=%14.7e\n", state.getPotentialEnergy(), referenceState.getPotentialEnergy() );
       }
       
       // Take a small step in the direction of the energy gradient.
       
       DoubleVector stats;
       if( compareForcesOfTwoStates( state, referenceState, 0.001, stats, log ) ){
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
   
       if( log ){
           (void) fprintf( log, "Fsum [%14.7e %14.7e %14.7e] norm=%14.7e\n", forceSum[0], forceSum[1], forceSum[2], norm );
       }
   
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
   
       if( log ){
           (void) fprintf( log, "%2d Energies %.8e %.8e norms[%13.7e %13.7e] deltaNorms=%13.7e delta=%.2e\n",
                           ii, state.getPotentialEnergy(), state2.getPotentialEnergy(), diff, norm, off, delta );
       }
   
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
           if( log ){
               (void) fprintf( log, "r48=%14.6e r28=%14.6e r24=%14.6e\n", positions[8][2]-positions[4][2], positions[8][2], positions[4][2] );
           }
       }
   }
}

int main() {

    try {

        registerFreeEnergyCudaKernelFactories( );

        //FILE* log = stderr;
        FILE* log = NULL;

        testSingleParticle( log );

        testEnergyEthaneSwitchingFunction( 0, log );
        testEnergyEthaneSwitchingFunction( 1, log );

    } catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
