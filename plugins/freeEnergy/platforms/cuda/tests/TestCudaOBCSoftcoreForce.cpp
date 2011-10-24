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

#include "TestCudaSoftcoreForce.h"

//#define USE_SOFTCORE
//#define IMPLICIT_SOLVENT GBVI
//#define IMPLICIT_SOLVENT OBC

#define OBC_FLAG  1
#define GBVI_FLAG 2

#include "openmm/GBVIForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/NonbondedForce.h"

#ifdef USE_SOFTCORE
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/GBSAOBCSoftcoreForce.h"
#include "openmm/NonbondedSoftcoreForce.h"
#endif

#include <iomanip>

void testSingleParticle( FILE* log ) {

    System system;
    system.addParticle(2.0);
    VerletIntegrator integrator(0.01);

    GBVISoftcoreForce* forceField = new GBVISoftcoreForce;

    double charge         = 1.0;
    double radius         = 0.15;
    double gamma          = 1.0;
    forceField->addParticle(charge, radius, gamma);
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

    double bornRadius     = radius; 
    double eps0           = EPSILON0;
    double tau            = (1.0/forceField->getSoluteDielectric()-1.0/forceField->getSolventDielectric());

    double bornEnergy     = (-charge*charge/(8*PI_M*eps0))*tau/bornRadius;
    double nonpolarEnergy = -gamma*tau*std::pow( radius/bornRadius, 3.0);

    double expectedE      = (bornEnergy+nonpolarEnergy); 
    double obtainedE      = state.getPotentialEnergy(); 
    double diff           = fabs( obtainedE - expectedE );
    if( log ){
        (void) fprintf( log, "testSingleParticle expected=%14.6e obtained=%14.6e diff=%14.6e breakdown:[%14.6e %14.6e]\n",
                        expectedE, obtainedE, diff, bornEnergy, nonpolarEnergy );
    }
    ASSERT_EQUAL_TOL((bornEnergy+nonpolarEnergy), state.getPotentialEnergy(), 0.01);
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
           (void) fprintf( log, "H=%d C=%d r=%14.6e\n", hydrogenIndex, carbonIndex, dist );
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

static GBVISoftcoreForce* copyGbviSoftcoreForce( const GBVISoftcoreForce& gbviSoftcoreForce ){

    GBVISoftcoreForce* copyGbviSoftcoreForce = new GBVISoftcoreForce(gbviSoftcoreForce);
/*
    GBVISoftcoreForce* copyGbviSoftcoreForce = new GBVISoftcoreForce();

    copyGbviSoftcoreForce->setNonbondedMethod( gbviSoftcoreForce.getNonbondedMethod() );

    copyGbviSoftcoreForce->setCutoffDistance( gbviSoftcoreForce.getCutoffDistance() );

    copyGbviSoftcoreForce->setSolventDielectric( gbviSoftcoreForce.getSolventDielectric() );
    copyGbviSoftcoreForce->setSoluteDielectric( gbviSoftcoreForce.getSoluteDielectric() );

    copyGbviSoftcoreForce->setBornRadiusScalingMethod( gbviSoftcoreForce.getBornRadiusScalingMethod() );
    copyGbviSoftcoreForce->setQuinticLowerLimitFactor( gbviSoftcoreForce.getQuinticLowerLimitFactor() );
    copyGbviSoftcoreForce->setQuinticUpperBornRadiusLimit( gbviSoftcoreForce.getQuinticUpperBornRadiusLimit() );

    // particle parameters

    for( unsigned int ii = 0; ii < gbviSoftcoreForce.getNumParticles(); ii++ ){

        double charge;
        double sigma;
        double gamma;
        double softcoreLJLambda;
        gbviSoftcoreForce.getParticleParameters(ii, charge, sigma, gamma, softcoreLJLambda);
        copyGbviSoftcoreForce->addParticle( charge, sigma, gamma, softcoreLJLambda);
    }

    // bonds

    for( unsigned int ii = 0; ii < gbviSoftcoreForce.getNumBonds(); ii++ ){
        int particle1, particle2;
        double distance;
        gbviSoftcoreForce.getBondParameters( ii, particle1, particle2, distance);
        copyGbviSoftcoreForce->addBond( particle1, particle2, distance );
    }
*/
    return copyGbviSoftcoreForce;
}

static GBVIForce* copyGbviForce( const GBVIForce& gbviForce ){
    return new GBVIForce(gbviForce);
}

static GBSAOBCSoftcoreForce* copyGBSAOBCSoftcoreForce( const GBSAOBCSoftcoreForce& gbviSoftcoreForce ){
    return new GBSAOBCSoftcoreForce(gbviSoftcoreForce);
}

static GBSAOBCForce* copyGbsaObcForce( const GBSAOBCForce& gbviForce ){
    return new GBSAOBCForce(gbviForce);
}

void testGbviSoftcore( MapStringToDouble& inputArgumentMap, FILE* log ){

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
   
    if( nonbondedMethod == 2 && cutoffDistance > boxSize*0.5 ){
        cutoffDistance = boxSize*0.5;
    }

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

    // Create two systems: one with GbviSoftcoreForce NonbondedSoftcoreForce forces, and one using a CustomNonbondedForce, CustomGBVI force to implement the same interaction.

    System standardSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));

#ifdef USE_SOFTCORE
    NonbondedSoftcoreForce* nonbondedSoftcoreForce   = new NonbondedSoftcoreForce();
    if( nonbondedMethod == NoCutoff ){
        nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::CutoffNonPeriodic );
        } else {
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::CutoffPeriodic );
        }
    }
#else
    NonbondedForce* nonbondedSoftcoreForce = new NonbondedForce();
    if( nonbondedMethod == NoCutoff ){
        nonbondedSoftcoreForce->setNonbondedMethod( NonbondedForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedForce::CutoffNonPeriodic );
        } else {
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedForce::CutoffPeriodic );
        }
    }
#endif
    nonbondedSoftcoreForce->setCutoffDistance( cutoffDistance );
    nonbondedSoftcoreForce->setReactionFieldDielectric( reactionFieldDielectric );

#ifdef USE_SOFTCORE

#if IMPLICIT_SOLVENT == GBVI_FLAG
    GBVISoftcoreForce* gbviSoftcoreForce             = new GBVISoftcoreForce();
    if( nonbondedMethod == NoCutoff ){
        gbviSoftcoreForce->setNonbondedMethod( GBVISoftcoreForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            gbviSoftcoreForce->setNonbondedMethod( GBVISoftcoreForce::CutoffNonPeriodic );
        } else {
            gbviSoftcoreForce->setNonbondedMethod( GBVISoftcoreForce::CutoffPeriodic );
        }
    }
#else
    GBSAOBCSoftcoreForce* gbviSoftcoreForce          = new GBSAOBCSoftcoreForce();
    if( nonbondedMethod == NoCutoff ){
        gbviSoftcoreForce->setNonbondedMethod( GBSAOBCSoftcoreForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            gbviSoftcoreForce->setNonbondedMethod( GBSAOBCSoftcoreForce::CutoffNonPeriodic );
        } else {
            gbviSoftcoreForce->setNonbondedMethod( GBSAOBCSoftcoreForce::CutoffPeriodic );
        }
    }
#endif

#else

#if IMPLICIT_SOLVENT == GBVI_FLAG
    GBVIForce* gbviSoftcoreForce           = new GBVIForce();
    if( nonbondedMethod == NoCutoff ){
        gbviSoftcoreForce->setNonbondedMethod( GBVIForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            gbviSoftcoreForce->setNonbondedMethod( GBVIForce::CutoffNonPeriodic );
        } else {
            gbviSoftcoreForce->setNonbondedMethod( GBVIForce::CutoffPeriodic );
        }
    }

#else

    GBSAOBCForce* gbviSoftcoreForce           = new GBSAOBCForce();
    if( nonbondedMethod == NoCutoff ){
        gbviSoftcoreForce->setNonbondedMethod( GBSAOBCForce::NoCutoff );
    } else {
        if( nonbondedMethod == CutoffNonPeriodic ){
            gbviSoftcoreForce->setNonbondedMethod( GBSAOBCForce::CutoffNonPeriodic );
        } else {
            gbviSoftcoreForce->setNonbondedMethod( GBSAOBCForce::CutoffPeriodic );
        }
    }

#endif

#endif

#if IMPLICIT_SOLVENT == GBVI_FLAG
#ifdef USE_SOFTCORE
    if( useQuinticSpline ){
        gbviSoftcoreForce->setBornRadiusScalingMethod( GBVISoftcoreForce::QuinticSpline );
    } else {
        gbviSoftcoreForce->setBornRadiusScalingMethod( GBVISoftcoreForce::NoScaling );
    }
#else
    if( useQuinticSpline ){
        gbviSoftcoreForce->setBornRadiusScalingMethod( GBVIForce::QuinticSpline );
    } else {
        gbviSoftcoreForce->setBornRadiusScalingMethod( GBVIForce::NoScaling );
    }
#endif
#endif

    gbviSoftcoreForce->setSolventDielectric( 78.3 );
    //gbviSoftcoreForce->setSolventDielectric( 1.0e+10 );
    //gbviSoftcoreForce->setSolventDielectric( 1.0 );
    gbviSoftcoreForce->setSoluteDielectric( 1.0 );
    gbviSoftcoreForce->setCutoffDistance( nonbondedSoftcoreForce->getCutoffDistance( ) );

    std::vector<Vec3> positions(numParticles);

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

    // show info on particle positions

    if( log ){
        Vec3 box[2];
        positionGenerator.getEnclosingBox( positions, box );
        (void) fprintf( log, "Enclosing Box (in A): [%15.7e %15.7e] [%15.7e %15.7e] [%15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n",
                        box[0][0], box[1][0], box[0][1], box[1][1], box[0][2], box[1][2],
                        (box[1][0] - box[0][0]), (box[1][1] - box[0][1]), (box[1][2] - box[0][2]) );

        int showIndex                        = 5;
        int periodicBoundaryConditions       = (nonbondedMethod == 2) ? 1 : 0;

        IntVector positionIndexVector;
        positionIndexVector.push_back( 0 );
        positionIndexVector.push_back( static_cast<int>(positions.size())-1 );
        //positionIndexVector.push_back( 542 );

        for( unsigned int ii = 0; ii < positionIndexVector.size(); ii++ ){
            if( positionIndexVector[ii] < positions.size() ){
                int positionIndex = positionIndexVector[ii];
                IntDoublePairVector sortVector;
                positionGenerator.getSortedDistances( periodicBoundaryConditions, positionIndex, positions, sortVector );
                (void) fprintf( log, "Min/max distance from %6d:\n    ", positionIndex );
                for( unsigned int jj = 0; jj < sortVector.size() && jj < showIndex; jj++ ){
                    IntDoublePair pair = sortVector[jj];
                    (void) fprintf( log, "[%6d %15.7e] ", pair.first, pair.second);
                }
                (void) fprintf( log, "\n    " );
                for( unsigned int jj = (sortVector.size() - showIndex); jj < sortVector.size() && jj >= 0; jj++ ){
                    IntDoublePair pair = sortVector[jj];
                    (void) fprintf( log, "[%6d %15.7e] ", pair.first, pair.second);
                }
                (void) fprintf( log, "\n" );
            }
        }
        IntIntPairVector pairs;
        pairs.push_back( IntIntPair( 732, 0 ) );
        pairs.push_back( IntIntPair( 732, 1 ) );
        pairs.push_back( IntIntPair( 732, 2 ) );
        pairs.push_back( IntIntPair( 732, 3 ) );
        pairs.push_back( IntIntPair( 732, 4 ) );
        for( IntIntPairVectorCI ii = pairs.begin(); ii != pairs.end(); ii++ ){
            if( ii->first < positions.size() && ii->second < positions.size() ){
                 double d = positionGenerator.getDistance( ii->first, ii->second, positions );
                 (void) fprintf( log, "Distance %6d %6d  %15.7e d2=%15.7e\n", ii->first, ii->second,  d, d*d );
            }
        }
    }    

    const int numberOfParameters             = 5;

    const int ChargeIndex                    = 0;
    const int SigmaIndex                     = 1;
    const int EpsIndex                       = 2;
    const int GammaIndex                     = 3;
    const int LambdaIndex                    = 4;

    std::vector<double> parameterLowerBound( numberOfParameters, 0.0 );

    double fixedCharge                       = 1.0;
    parameterLowerBound[ChargeIndex]         = fixedCharge;  // charge
    parameterLowerBound[SigmaIndex]          = 0.1;          // sigma
    parameterLowerBound[EpsIndex]            = 0.5;          // eps
    parameterLowerBound[GammaIndex]          = 0.1;          // gamma
    parameterLowerBound[LambdaIndex]         = lambda1;      // lambda

    std::vector<double> parameterUpperBound( parameterLowerBound );
    parameterUpperBound[ChargeIndex]         = fixedCharge;  // charge
    parameterUpperBound[SigmaIndex]          = 0.3;          // sigma
    parameterUpperBound[EpsIndex]            = 40.0;         // eps
    parameterUpperBound[GammaIndex]          = 40.0;         // gamma

#if IMPLICIT_SOLVENT == OBC_FLAG
    parameterLowerBound[GammaIndex]          = 0.1;          // overlap factor
    parameterUpperBound[GammaIndex]          = 1.5;        
#endif

    std::vector<double> parameters( numberOfParameters );
    double charge = fixedCharge;

    for( int ii = 0; ii < numMolecules; ii++) {

        charge       *= -1.0;

        double lambda =  ii < (numMolecules/2) ? lambda1 : lambda2;
        randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );

#ifdef USE_SOFTCORE
        nonbondedSoftcoreForce->addParticle(   charge,  parameters[SigmaIndex],  parameters[EpsIndex],    lambda );
        gbviSoftcoreForce->addParticle(        charge,  parameters[SigmaIndex],  parameters[GammaIndex],  lambda );
#else
        nonbondedSoftcoreForce->addParticle(   charge,  parameters[SigmaIndex],  parameters[EpsIndex] );
        gbviSoftcoreForce->addParticle(        charge,  parameters[SigmaIndex],  parameters[GammaIndex] );
#endif

        int baseParticleIndex                    = ii*numParticlesPerMolecule;
        for( int jj = 1; jj < numParticlesPerMolecule; jj++) {

            // alternate charges

            charge *= -1.0;

            randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );

#ifdef USE_SOFTCORE
            nonbondedSoftcoreForce->addParticle(   charge,  parameters[SigmaIndex],  parameters[EpsIndex],    lambda );
            gbviSoftcoreForce->addParticle(        charge,  parameters[SigmaIndex],  parameters[GammaIndex],  lambda );
#else
            nonbondedSoftcoreForce->addParticle(   charge,  parameters[SigmaIndex],  parameters[EpsIndex] );
            gbviSoftcoreForce->addParticle(        charge,  parameters[SigmaIndex],  parameters[GammaIndex] );
#endif

            nonbondedSoftcoreForce->addException( baseParticleIndex, baseParticleIndex+jj, 0.0f, 1.0, 0.0f );

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

    standardSystem.addForce(nonbondedSoftcoreForce);
    if( includeGbvi ){
        standardSystem.addForce(gbviSoftcoreForce);
    }

    // copy system and forces

    System* systemCopy = copySystem( standardSystem );

#ifdef USE_SOFTCORE
    NonbondedSoftcoreForce* nonbondedSoftcoreForceCopy;
    nonbondedSoftcoreForceCopy = copyNonbondedSoftcoreForce( *nonbondedSoftcoreForce );
#else
    NonbondedForce* nonbondedSoftcoreForceCopy;
    nonbondedSoftcoreForceCopy = copyNonbondedForce( *nonbondedSoftcoreForce );
#endif
    systemCopy->addForce( nonbondedSoftcoreForceCopy );
    std::stringstream baseFileName;

    if( includeGbvi ){
#ifdef USE_SOFTCORE

#if IMPLICIT_SOLVENT == GBVI_FLAG
        GBVISoftcoreForce* gBVISoftcoreForceCopy  = copyGbviSoftcoreForce( *gbviSoftcoreForce );
        baseFileName  << "GBVISoftcore";
#endif
#if IMPLICIT_SOLVENT == OBC_FLAG
        baseFileName  << "GBSAObcSoftcore";
        GBSAOBCSoftcoreForce* gBVISoftcoreForceCopy       = copyGBSAOBCSoftcoreForce( *gbviSoftcoreForce );
#endif
        baseFileName  << "_lbda" << std::fixed << setprecision(2) << lambda2;

#else

#if IMPLICIT_SOLVENT == GBVI_FLAG
        GBVIForce* gBVISoftcoreForceCopy          = copyGbviForce( *gbviSoftcoreForce );
        baseFileName  << "Gbvi";
#endif
#if IMPLICIT_SOLVENT == OBC_FLAG
        GBSAOBCForce* gBVISoftcoreForceCopy       = copyGbsaObcForce( *gbviSoftcoreForce );
        baseFileName  << "GBSAOBC";
#endif

#endif
        systemCopy->addForce( gBVISoftcoreForceCopy );
    }

    // perform comparison

    std::stringstream idString;
    idString << "Nb " << nonbondedMethod << " l2 " << std::fixed << setprecision(2) << lambda2;
    runSystemComparisonTest( standardSystem, *systemCopy, "Cuda", "Reference", positions, inputArgumentMap, idString.str(), log );

    // serialize

    baseFileName  << "_N"     << positions.size();
    baseFileName  << "_Nb"    << nonbondedMethod;
    serializeSystemAndPositions( standardSystem, positions, baseFileName.str(), log);

    delete systemCopy;

}

int main() {

    try {

        registerFreeEnergyCudaKernelFactories( );

        VectorOfMapStringToDouble vectorOfMapStringToDouble;
        MapStringToDouble inputArgumentMap;
        MapStringToDoubleVector generativeArgumentMaps;
        //FILE* log = stderr;
        FILE* log = NULL;
/*
        testSingleParticle( log );

        testEnergyEthaneSwitchingFunction( 0, log );
        testEnergyEthaneSwitchingFunction( 1, log );
*/

        inputArgumentMap["lambda2"]                         = 1.0;
        inputArgumentMap["nonbondedMethod"]                 = 0;
        inputArgumentMap["numMolecules"]                    = 10;
        inputArgumentMap["boxSize"]                         = 5.0;
        inputArgumentMap["positionPlacementMethod"]         = 0;
        inputArgumentMap["cutoffDistance"]                  = 0.3*inputArgumentMap["boxSize"];
        //inputArgumentMap["cutoffDistance"]                  = 1.0;
        inputArgumentMap["relativeTolerance"]               = 5.0e-04;
        inputArgumentMap["serialize"]                       = 1;
        //inputArgumentMap["numParticlesPerMolecule"]         = 2;

#ifdef USE_SOFTCORE
        DoubleVector lamda2;
        lamda2.push_back( 1.0 );
        lamda2.push_back( 0.5 );
        lamda2.push_back( 0.0 );
        if( lamda2.size() > 0 ){
            generativeArgumentMaps["lambda2"] = lamda2;
            inputArgumentMap["lambda2"]       = lamda2[0];
        }   
#endif

        DoubleVector numberOfMolecules;
        numberOfMolecules.push_back( 10 );
        numberOfMolecules.push_back( 100 );
        numberOfMolecules.push_back( 1000 );
        //numberOfMolecules.push_back( 2000 );
        //numberOfMolecules.push_back( 4000 );
        //numberOfMolecules.push_back( 8000 );
        if( numberOfMolecules.size() > 0 ){
            generativeArgumentMaps["numMolecules"] = numberOfMolecules;
            inputArgumentMap["numMolecules"]       = numberOfMolecules[0];
        }   

        DoubleVector nonbondedMethod;
        nonbondedMethod.push_back( 0 );
        nonbondedMethod.push_back( 1 );
        nonbondedMethod.push_back( 2 );
        if( nonbondedMethod.size() > 0 ){
            generativeArgumentMaps["nonbondedMethod"] = nonbondedMethod;
            inputArgumentMap["nonbondedMethod"]       = nonbondedMethod[0];
        }

        vectorOfMapStringToDouble.push_back( inputArgumentMap );
        generateInputArgumentMapsFromStringVectors( generativeArgumentMaps, vectorOfMapStringToDouble ); 

        // big box/many particle tests

        //bool bigBox = true;
        bool bigBox = false;
        if( bigBox ){
            MapStringToDouble inputArgumentMapBig;
            VectorOfMapStringToDouble vectorOfMapStringToDoubleBig;
            inputArgumentMapBig["lambda2"]                         = 1.0;
            inputArgumentMapBig["nonbondedMethod"]                 = 1;
            inputArgumentMapBig["numMolecules"]                    = 10;
            inputArgumentMapBig["boxSize"]                         = 20.0;
            inputArgumentMapBig["relativeTolerance"]               = 6.0e-04;
            vectorOfMapStringToDoubleBig.push_back( inputArgumentMapBig );
            //MapStringToDoubleVector generativeArgumentMapsBig;

            numberOfMolecules.resize( 0 );
            numberOfMolecules.push_back( 4000 );
            generativeArgumentMaps["numMolecules"] = numberOfMolecules;

            nonbondedMethod.resize( 0 );
            nonbondedMethod.push_back( 1 );
            nonbondedMethod.push_back( 2 );
            generativeArgumentMaps["nonbondedMethod"] = nonbondedMethod;
            generateInputArgumentMapsFromStringVectors( generativeArgumentMaps, vectorOfMapStringToDoubleBig ); 
            vectorOfMapStringToDouble.resize( 0 );
            vectorOfMapStringToDouble.insert( vectorOfMapStringToDouble.end(), vectorOfMapStringToDoubleBig.begin(), vectorOfMapStringToDoubleBig.end() );
        }

        if( log ){
            MapStringToInt exclude;
            exclude["lambda1"]                 = 1;
            exclude["numParticlesPerMolecule"] = 1;
            std::stringstream outputStream;
            std::sort( vectorOfMapStringToDouble.begin(), vectorOfMapStringToDouble.end(), TestMapSortPredicate);
            StringVector printOrder;
            printOrder.push_back( "numMolecules" );
            printOrder.push_back( "nonbondedMethod" );
            printOrder.push_back( "lambda2" );
            printOrder.push_back( "boxSize" );
            for( unsigned int kk = 0; kk < vectorOfMapStringToDouble.size(); kk++ ){
                streamArgumentMapOneLine( vectorOfMapStringToDouble[kk], exclude, printOrder, kk, outputStream );
            }
            (void) fprintf( log, "Initial argument maps: %u\n%s", static_cast<unsigned int>(vectorOfMapStringToDouble.size()), outputStream.str().c_str() );
        }

        // run tests

        for( unsigned int kk = 0; kk < vectorOfMapStringToDouble.size(); kk++ ){
            testGbviSoftcore( vectorOfMapStringToDouble[kk], log );
            sleep(2);
        }

    } catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
