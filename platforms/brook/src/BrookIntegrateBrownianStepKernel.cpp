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

#include "BrookIntegrateBrownianStepKernel.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookIntegrateBrownianStepKernel constructor
 * 
 * @param name        name of the stream to create
 * @param platform    platform
 *
 */

BrookIntegrateBrownianStepKernel::BrookIntegrateBrownianStepKernel( std::string name, const Platform& platform ) :
                                  IntegrateBrownianStepKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateBrownianStepKernel::BrookIntegrateBrownianStepKernel";

// ---------------------------------------------------------------------------------------

   _brookBrownianDynamics          = NULL;
   _brookShakeAlgorithm            = NULL;
   _brookRandomNumberGenerator     = NULL;

}

/** 
 * BrookIntegrateVerletStepKernel destructor
 * 
 */
  
BrookIntegrateBrownianStepKernel::~BrookIntegrateBrownianStepKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateBrownianStepKernel::~BrookIntegrateBrownianStepKernel";

// ---------------------------------------------------------------------------------------
   
   delete _brookBrownianDynamics;
   delete _brookShakeAlgorithm;
   delete _brookRandomNumberGenerator;

}

/** 
 * Initialize the kernel, setting up all parameters related to integrator.
 * 
 * @param masses             the mass of each particle
 * @param constraintIndices  each element contains the indices of two particles whose distance should be constrained
 * @param constraintLengths  the required distance between each pair of constrained particles
 *
 */

void BrookIntegrateBrownianStepKernel::initialize( const vector<double>& masses,
                                                   const vector<vector<int> >& constraintIndices,
                                                   const vector<double>& constraintLengths ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateBrownianStepKernel::initialize";

// ---------------------------------------------------------------------------------------
   
   _brookBrownianDynamics        = new BrookBrownianDynamics( );
   _brookBrownianDynamics->setup( masses, getPlatform() );

   _brookShakeAlgorithm          = new BrookShakeAlgorithm( );
   _brookShakeAlgorithm->setup( masses, constraintIndices, constraintLengths, getPlatform() );

   // assert( (_brookShakeAlgorithm->getNumberOfConstraints() > 0) );

   _brookRandomNumberGenerator   = new BrookRandomNumberGenerator( );
   _brookRandomNumberGenerator->setup( (int) masses.size(), getPlatform() );
//   _brookRandomNumberGenerator->setVerbosity( 1 );
}

/** 
 * Execute kernel
 * 
 * @param positions          particle coordinates
 * @param velocities         particle velocities
 * @param forces             particle forces
 * @param temperature        heat bath temperature
 * @param friction           friction coefficient coupling the system to the heat bath
 * @param stepSize           integration step size
 *
 */

void BrookIntegrateBrownianStepKernel::execute( Stream& positions, Stream& velocities,
                                                const Stream& forces, double temperature,
                                                double friction, double stepSize ){

// ---------------------------------------------------------------------------------------

   double epsilon                           = 1.0e-04;
   static const std::string methodName      = "BrookIntegrateBrownianStepKernel::execute";

// ---------------------------------------------------------------------------------------
   
   // first time through initialize _brookBrownianDynamics

   // for each subsequent call, check if parameters need to be updated due to a change
   // in T, gamma, or the step size

   // take step

   double differences[3];
   differences[0] = temperature - (double) _brookBrownianDynamics->getTemperature();
   differences[1] = friction    - (double) _brookBrownianDynamics->getFriction();
   differences[2] = stepSize    - (double) _brookBrownianDynamics->getStepSize();
   if( fabs( differences[0] ) > epsilon || fabs( differences[1] ) > epsilon || fabs( differences[2] ) > epsilon ){
//printf( "%s calling updateParameters\n", methodName.c_str() );
      _brookBrownianDynamics->updateParameters( temperature, friction, stepSize );
   } else {
//printf( "%s NOT calling updateParameters\n", methodName.c_str() );
   }

   _brookBrownianDynamics->update( positions, velocities, forces, *_brookShakeAlgorithm, *_brookRandomNumberGenerator );

}
