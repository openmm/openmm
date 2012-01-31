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
 * This tests all the different force terms in the reference implementation of CustomGBForce.
 */

#include "openmm/internal/AssertionUtilities.h"

#include "sfmt/SFMT.h"
#include "openmm/Context.h"
#include "openmm/CustomBondForce.h"
#include "openmm/CustomNonbondedForce.h"

#include "openmm/NonbondedForce.h"
#include "openmm/NonbondedSoftcoreForce.h"

#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>
#include <algorithm>

extern "C" void registerFreeEnergyCudaKernelFactories();

using namespace OpenMM;
using namespace std;

const double TOL = 1e-4;

static const int NoCutoff          = 0;
static const int CutoffNonPeriodic = 1;
static const int CutoffPeriodic    = 2;

void testNonbondedSoftcore( double lambda1, double lambda2, int nonbondedMethod, FILE* log  ){

    const int numMolecules               = 70;
    const int numParticles               = numMolecules*2;
    const double boxSize                 = 10.0;
    const double reactionFieldDielectric = 80.0;
    const double cutoffDistance          = 0.4*boxSize;

    // Create two systems: one with a NonbondedSoftcoreForce, and one using a CustomNonbondedForce to implement the same interaction.

    System standardSystem;
    System customSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
        customSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    customSystem.setDefaultPeriodicBoxVectors(  Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));

    NonbondedSoftcoreForce* nonbondedSoftcoreForce   = new NonbondedSoftcoreForce();
    CustomNonbondedForce* customNonbonded;
    CustomBondForce* customBond;
    if( nonbondedMethod == NoCutoff ){

        nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::NoCutoff );

        customNonbonded          = new CustomNonbondedForce("lambda*4*eps*(dem^2-dem)+138.935456*q/r;"
                                                            "q=q1*q2;"
                                                            "dem=1.0/(soft+rsig);"
                                                            "rsig=(r/sigma)^6;"
                                                            "rsig=(r/sigma)^6;"
                                                            "soft=0.5*(1.0-lambda);"
                                                            "sigma=0.5*(sigma1+sigma2);"
                                                            "eps=sqrt(eps1*eps2);"
                                                            "lambda=min(lambda1,lambda2)");

        customNonbonded->setNonbondedMethod( CustomNonbondedForce::NoCutoff );

        customBond               = new CustomBondForce("lambda*4*eps*(dem^2-dem)+138.935456*q/r;"
                                                       "dem=1.0/(soft+rsig);"
                                                       "rsig=(r/sigma)^6;"
                                                       "soft=0.5*(1.0-lambda)");

    } else {

        nonbondedSoftcoreForce->setCutoffDistance( cutoffDistance );
        nonbondedSoftcoreForce->setReactionFieldDielectric( reactionFieldDielectric );
        if( nonbondedMethod == CutoffNonPeriodic ){
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::CutoffNonPeriodic );
        } else {
            nonbondedSoftcoreForce->setNonbondedMethod( NonbondedSoftcoreForce::CutoffPeriodic );
        }

        customNonbonded          = new CustomNonbondedForce("lambda*4*eps*(dem^2-dem)+138.935456*q*(1.0/r+(krf*r*r)-crf);"
                                                            "q=q1*q2;"
                                                            "dem=1.0/(soft+rsig);"
                                                            "rsig=(r/sigma)^6;"
                                                            "rsig=(r/sigma)^6;"
                                                            "soft=0.5*(1.0-lambda);"
                                                            "sigma=0.5*(sigma1+sigma2);"
                                                            "eps=sqrt(eps1*eps2);"
                                                            "lambda=min(lambda1,lambda2)");

        customBond               = new CustomBondForce("withinCutoff*(lambda*4*eps*(dem^2-dem)+138.935456*q*(1.0/r+(krf*r*r)-crf));"
                                                       "withinCutoff=step(cutoff-r);"
                                                       "dem=1.0/(soft+rsig);"
                                                       "rsig=(r/sigma)^6;"
                                                       "soft=0.5*(1.0-lambda)");
 
        customNonbonded->setCutoffDistance( cutoffDistance );
        if( nonbondedMethod == CutoffNonPeriodic ){
            customNonbonded->setNonbondedMethod( CustomNonbondedForce::CutoffNonPeriodic );
        } else {
            customNonbonded->setNonbondedMethod( CustomNonbondedForce::CutoffPeriodic );
        }

        double eps2               = (reactionFieldDielectric - 1.0)/(2.0*reactionFieldDielectric+1.0);
        double kValue             = eps2/(cutoffDistance*cutoffDistance*cutoffDistance);
        customNonbonded->addGlobalParameter("krf", kValue );

        customBond->addGlobalParameter("krf", kValue );

        double cValue             = (1.0/cutoffDistance)*(3.0*reactionFieldDielectric)/(2.0*reactionFieldDielectric + 1.0); 
        customNonbonded->addGlobalParameter("crf", cValue );
        customBond->addGlobalParameter("crf", cValue );
        customBond->addGlobalParameter("cutoff", cutoffDistance );
    }

    customNonbonded->addPerParticleParameter("q");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    customNonbonded->addPerParticleParameter("lambda");

    customBond->addPerBondParameter("q");
    customBond->addPerBondParameter("sigma");
    customBond->addPerBondParameter("eps");
    customBond->addPerBondParameter("lambda");

    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<double> params(4);

    // periodic boundary conditions not possible w/ CustomBond?

    int includeExceptions = nonbondedMethod != NoCutoff ? 0 : 1;

    for (int i = 0; i < numMolecules; i++) {
        if (i < numMolecules/2) {

            double charge = 1.0;

            double sigma1 = 0.2;
            double sigma2 = 0.2;

            double eps1   = 0.5;
            double eps2   = 0.5;
//eps1 = eps2 = 0.0;
            nonbondedSoftcoreForce->addParticle( charge, sigma1, eps1, lambda1);
            nonbondedSoftcoreForce->addParticle(-charge, sigma2, eps2, lambda1);

            params[0] = charge;
            params[1] = sigma1;
            params[2] = eps1;
            params[3] = lambda1;

            customNonbonded->addParticle(params);

            params[0] = -charge;
            params[1] = sigma2;
            params[2] = eps2;
            customNonbonded->addParticle(params);

            if( includeExceptions && i && ((i%4) == 0) ){
                vector<double> bondParams(4);
                nonbondedSoftcoreForce->addException(i-4, i, charge*charge, sigma1, eps1, false, lambda1);
                customNonbonded->addExclusion( i-4,i);
                bondParams[0] = charge*charge;
                bondParams[1] = sigma1;
                bondParams[2] = eps1;
                bondParams[3] = lambda1;
                customBond->addBond(i-4,i, bondParams );
            }

        } else {

            double charge = 1.0;

            double sigma1 = 0.2;
            double sigma2 = 0.1;

            double eps1   = 0.8;
            double eps2   = 0.8;

//eps1 = eps2 = 0.0;
            nonbondedSoftcoreForce->addParticle( charge, sigma1, eps1, lambda2);
            nonbondedSoftcoreForce->addParticle(-charge, sigma2, eps2, lambda2);

            params[0] = charge;
            params[1] = sigma1;
            params[2] = eps1;
            params[3] = lambda2;
            customNonbonded->addParticle(params);

            params[0] = -charge;
            params[1] = sigma2;
            params[2] = eps2;
            customNonbonded->addParticle(params);

            if( includeExceptions && i && ((i%4) == 0) ){

                vector<double> bondParams(4);
                nonbondedSoftcoreForce->addException(i-4, i, charge*charge, sigma1, eps1, false, lambda2);
                customNonbonded->addExclusion( i-4,i);
                bondParams[0] = charge*charge;
                bondParams[1] = sigma1;
                bondParams[2] = eps1;
                bondParams[3] = lambda2;
                customBond->addBond(i-4,i, bondParams );
            }
        }

        positions[2*i]    = Vec3(boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt), boxSize*genrand_real2(sfmt));
        positions[2*i+1]  = Vec3(positions[2*i][0]+1.0, positions[2*i][1], positions[2*i][2]);
        velocities[2*i]   = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));
        velocities[2*i+1] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));

    }

    standardSystem.addForce(nonbondedSoftcoreForce);
    customSystem.addForce(customNonbonded);
    customSystem.addForce(customBond);

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);

    Context context1(standardSystem, integrator1, Platform::getPlatformByName( "Cuda"));
    context1.setPositions(positions);
    context1.setVelocities(velocities);

    State state1 = context1.getState(State::Forces | State::Energy);

    Context context2(customSystem, integrator2, Platform::getPlatformByName( "Reference"));
    context2.setPositions(positions);
    context2.setVelocities(velocities);
    State state2 = context2.getState(State::Forces | State::Energy);

    double diff = 0.0;
    static int previousNonbondedMethod = -1;
    if( state1.getPotentialEnergy() - state2.getPotentialEnergy() != 0.0 ){
        diff = fabs( state1.getPotentialEnergy() - state2.getPotentialEnergy() )/( fabs( state1.getPotentialEnergy() ) + fabs( state2.getPotentialEnergy() ) );
    }
    if( previousNonbondedMethod != nonbondedMethod ){
        if( log ){
            (void) fprintf( log, "\n" );
        }
        previousNonbondedMethod = nonbondedMethod;
    }

    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), 1e-4);
    double maxDiff   = -1.0;
    int maxDiffIndex = -1;
    for (int i = 0; i < numParticles; i++) {
        Vec3 f1     = state1.getForces()[i];
        Vec3 f2     = state2.getForces()[i];

        double f1N  = sqrt( (f1[0]*f1[0]) + (f1[1]*f1[1]) + (f1[2]*f1[2]) );
        double f2N  = sqrt( (f2[0]*f2[0]) + (f2[1]*f2[1]) + (f2[2]*f2[2]) );

        double diff = (f1[0]-f2[0])*(f1[0]-f2[0]) +
                      (f1[1]-f2[1])*(f1[1]-f2[1]) +
                      (f1[2]-f2[2])*(f1[2]-f2[2]);
        if( f1N > 0.0 || f1N > 0.0 ){
            diff        = 2.0*sqrt( diff )/(f1N + f2N);
        }
        if( diff > maxDiff ){
            maxDiff      = diff;
            maxDiffIndex = i;
        }
/*
        (void) fprintf( log, "%4d %15.7e [%15.7e %15.7e %15.7e]  [%15.7e %15.7e %15.7e] \n", i, diff,
                        f1[0], f1[1],f1[2], f2[0], f2[1],f2[2] );
*/
//        ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], 1e-4);
         ASSERT( diff < 1.0e-3);
    }
    
    if( log ){
        (void) fprintf( log, "%d %d %10.1f %10.1f  %15.7e %15.7e %15.7e maxFDiff=%15.7e %d\n",
                        nonbondedMethod, includeExceptions, lambda1, lambda2, diff,
                        state1.getPotentialEnergy(), state2.getPotentialEnergy(), maxDiff, maxDiffIndex); fflush( log );
    }

}

int main() {

    try {

        registerFreeEnergyCudaKernelFactories( );

        // test various combinations of lambdas and boundary conditions/cutoffs
        FILE* log = NULL;
        testNonbondedSoftcore( 1.0, 1.0 , NoCutoff, log );
        testNonbondedSoftcore( 1.0, 0.0 , NoCutoff, log );
        testNonbondedSoftcore( 1.0, 0.5 , NoCutoff, log );
        testNonbondedSoftcore( 0.0, 0.0 , NoCutoff, log );

        testNonbondedSoftcore( 1.0, 1.0 , CutoffNonPeriodic, log );
        testNonbondedSoftcore( 1.0, 0.0 , CutoffNonPeriodic, log );
        testNonbondedSoftcore( 1.0, 0.5 , CutoffNonPeriodic, log );
        testNonbondedSoftcore( 0.0, 0.0 , CutoffNonPeriodic, log );


        testNonbondedSoftcore( 1.0, 1.0 , CutoffPeriodic, log );
        testNonbondedSoftcore( 1.0, 0.0 , CutoffPeriodic, log );
        testNonbondedSoftcore( 1.0, 0.5 , CutoffPeriodic, log );
        testNonbondedSoftcore( 0.0, 0.0 , CutoffPeriodic, log );

    } catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

