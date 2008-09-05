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

BrookIntegrateLangevinStepKernel::BrookIntegrateLangevinStepKernel( std::string name, const Platform& platform ) :
                                  IntegrateLangevinStepKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateLangevinStepKernel::BrookIntegrateLangevinStepKernel";

// ---------------------------------------------------------------------------------------

   _brookStochasticDynamics        = NULL;
   _brookShakeAlgorithm            = NULL;
   _brookRandomNumberGenerator     = NULL;

}

BrookIntegrateLangevinStepKernel::~BrookIntegrateLangevinStepKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateLangevinStepKernel::~BrookIntegrateLangevinStepKernel";

// ---------------------------------------------------------------------------------------
   
   delete _brookStochasticDynamics;
   delete _brookShakeAlgorithm;
   delete _brookRandomNumberGenerator;

}

void BrookIntegrateLangevinStepKernel::initialize( const vector<double>& masses,
                                                   const vector<vector<int> >& constraintIndices,
                                                   const vector<double>& constraintLengths ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateLangevinStepKernel::initialize";

// ---------------------------------------------------------------------------------------
   
   _brookStochasticDynamics      = new BrookStochasticDynamics( );
   _brookStochasticDynamics->setup( masses, getPlatform() );

   _brookShakeAlgorithm          = new BrookShakeAlgorithm( );
   _brookShakeAlgorithm->setup( masses, constraintIndices, constraintLengths, getPlatform() );

   _brookRandomNumberGenerator   = new BrookRandomNumberGenerator( );
   _brookRandomNumberGenerator->setup( (int) masses.size(), getPlatform() );
}

/** 
 * Execute kernel
 * 
 * @param positions          atom coordinates
 * @param velocities         atom velocities
 * @param forces             atom forces
 * @param temperature        heat bath temperature
 * @param friction           friction coefficient coupling the system to the heat bath
 * @param stepSize           integration step size
 *
 */

void BrookIntegrateLangevinStepKernel::execute( Stream& positions, Stream& velocities,
                                                const Stream& forces, double temperature,
                                                double friction, double stepSize ){

// ---------------------------------------------------------------------------------------

   double epsilon = 1.0e-04;
   static const std::string methodName      = "BrookIntegrateLangevinStepKernel::execute";

// ---------------------------------------------------------------------------------------
   
   // first time through initialize _brookStochasticDynamics

   // for each subsequent call, check if parameters need to be updated due to a change
   // in T, gamma, or the step size

   // take step

   double differences[3];
   differences[0] = temperature - (double) _brookStochasticDynamics->getTemperature();
   differences[1] = friction    - (double) _brookStochasticDynamics->getFriction();
   differences[2] = stepSize    - (double) _brookStochasticDynamics->getStepSize();
   if( fabs( differences[0] ) < epsilon || fabs( differences[1] ) < epsilon || fabs( differences[2] ) < epsilon ){
      _brookStochasticDynamics->updateParameters( temperature, friction, stepSize );
   } 

   assert( _brookShakeAlgorithm );
/*
   if( _brookShakeAlgorithm == NULL ){
      std::stringstream message;
      message << methodName << " _brookShakeAlgorithm is not set -- case not handled.";
      throw OpenMMException( message.str() );
      return NULL;
   }
*/

   _brookStochasticDynamics->update( positions, velocities, forces, *_brookShakeAlgorithm, *_brookRandomNumberGenerator );

}
