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

#include "BrookIntegrateLangevinStepKernel.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookIntegrateLangevinStepKernel constructor
 * 
 * @param name                  name of the stream to create
 * @param platform              platform
 * @param openMMBrookInterface  OpenMMBrookInterface reference
 * @param system                System reference  
 *
 */

BrookIntegrateLangevinStepKernel::BrookIntegrateLangevinStepKernel( std::string name, const Platform& platform,
                                  OpenMMBrookInterface& openMMBrookInterface, System& system ) :
                                  IntegrateLangevinStepKernel( name, platform ), _openMMBrookInterface( openMMBrookInterface ), _system( system ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateLangevinStepKernel::BrookIntegrateLangevinStepKernel";

// ---------------------------------------------------------------------------------------

   _brookLangevinDynamics          = NULL;
   _brookShakeAlgorithm            = NULL;
   _brookRandomNumberGenerator     = NULL;
   _log                            = NULL;

   const BrookPlatform& brookPlatform        = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }

}

/** 
 * BrookIntegrateVerletStepKernel destructor
 * 
 */
  
BrookIntegrateLangevinStepKernel::~BrookIntegrateLangevinStepKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateLangevinStepKernel::~BrookIntegrateLangevinStepKernel";

// ---------------------------------------------------------------------------------------
   
   delete _brookLangevinDynamics;
   delete _brookShakeAlgorithm;
   delete _brookRandomNumberGenerator;

}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookIntegrateLangevinStepKernel::getLog( void ) const {
   return _log;
}

/** 
 * Set log file reference
 * 
 * @param  log file reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookIntegrateLangevinStepKernel::setLog( FILE* log ){
   _log = log;
   return DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up all parameters related to integrator.
 * 
 * @param system                System reference  
 * @param integrator            LangevinIntegrator reference
 *
 */

void BrookIntegrateLangevinStepKernel::initialize( const System& system, const LangevinIntegrator& integrator ){

// ---------------------------------------------------------------------------------------

   int printOn                               = 0;
   static const std::string methodName       = "BrookIntegrateLangevinStepKernel::initialize";
   FILE* log                                 = NULL;

// ---------------------------------------------------------------------------------------
   
//setLog( stderr );
   printOn               = (printOn && getLog()) ? printOn : 0;

   if( printOn ){
      log = getLog();
      (void) fprintf( log, "%s\n", methodName.c_str() );
      (void) fflush( log );
   }

   int numberOfParticles = system.getNumParticles();

   // masses

   std::vector<double> masses;
   masses.resize( numberOfParticles );

   if( printOn ){
      (void) fprintf( log, "%s %d\n", methodName.c_str(), numberOfParticles );
      (void) fflush( log );
   }

   for( int ii = 0; ii < numberOfParticles; ii++ ){
      masses[ii] = static_cast<double>(system.getParticleMass(ii));
   }

   // constraints

   int numberOfConstraints = system.getNumConstraints();

   if( printOn ){
      (void) fprintf( log, "%s const=%d\n", methodName.c_str(), numberOfConstraints );
      (void) fflush( log );
   }

   std::vector<std::vector<int> > constraintIndicesVector;
   constraintIndicesVector.resize( numberOfConstraints );
   std::vector<double> constraintLengths;

   for( int ii = 0; ii < numberOfConstraints; ii++ ){

      int particle1, particle2;
      double distance;

      system.getConstraintParameters( ii, particle1, particle2, distance );

      constraintIndicesVector[ii].push_back( particle1 );
      constraintIndicesVector[ii].push_back( particle2 );
      constraintLengths.push_back( distance );

//(void) fprintf( log, "%s shake setup const=%d ", methodName.c_str(), ii ); fflush( log );
//(void) fprintf( log, "[ %d %d %f]\n", particle1, particle2, distance ); fflush( log );
   }

   _brookLangevinDynamics        = new BrookLangevinDynamics( );
   _brookLangevinDynamics->setup( masses, getPlatform() );

   _brookShakeAlgorithm          = new BrookShakeAlgorithm( );
   _brookShakeAlgorithm->setup( masses, constraintIndicesVector, constraintLengths, getPlatform() );

   // tolerance

   BrookOpenMMFloat tolerance = static_cast<BrookOpenMMFloat>( integrator.getConstraintTolerance() );
   _brookShakeAlgorithm->setShakeTolerance( tolerance );
   _brookShakeAlgorithm->setMaxIterations( 40 );
   if( log ){
      _brookShakeAlgorithm->setLog( log );
   }

   // random number generator

   _brookRandomNumberGenerator   = new BrookRandomNumberGenerator( );
   _brookRandomNumberGenerator->setup( (int) masses.size(), getPlatform() );

   if( printOn ){
      (void) fprintf( log, "%s done setup:\nBrookShakeAlgorithm:\n%s\nBrookRandomNumberGenerator:\n%s\n\n", methodName.c_str(),
                      _brookShakeAlgorithm->getContentsString().c_str(), 
                      _brookRandomNumberGenerator->getContentsString().c_str() );
      (void) fflush( log );
   }

}

/** 
 * Execute kernel
 * 
 * @param context            OpenMMContextImpl reference
 * @param integrator         LangevinIntegrator reference
 *
 */

void BrookIntegrateLangevinStepKernel::execute( OpenMMContextImpl& context, const LangevinIntegrator& integrator ){

// ---------------------------------------------------------------------------------------

   double epsilon                           = 1.0e-06;
   static const std::string methodName      = "BrookIntegrateLangevinStepKernel::execute";

// ---------------------------------------------------------------------------------------
   
   // first time through initialize _brookLangevinDynamics

   // for each subsequent call, check if parameters need to be updated due to a change
   // in T, gamma, or the step size

   // take step

   double differences[3];
   differences[0] = integrator.getTemperature() - (double) _brookLangevinDynamics->getTemperature();
   differences[1] = integrator.getFriction()    - (double) _brookLangevinDynamics->getFriction();
   differences[2] = integrator.getStepSize()    - (double) _brookLangevinDynamics->getStepSize();
   if( fabs( differences[0] ) > epsilon || fabs( differences[1] ) > epsilon || fabs( differences[2] ) > epsilon ){
      _brookLangevinDynamics->updateParameters( integrator.getTemperature(), integrator.getFriction(), integrator.getStepSize() );
   }

   _brookLangevinDynamics->update( *(_openMMBrookInterface.getParticlePositions()), *(_openMMBrookInterface.getParticleVelocities()), 
                                   *(_openMMBrookInterface.getParticleForces()), *_brookShakeAlgorithm, *_brookRandomNumberGenerator );

}
