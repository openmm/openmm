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

#include "BrookIntegrateVerletStepKernel.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookIntegrateVerletStepKernel constructor
 * 
 * @param name        name of the stream to create
 * @param platform    platform
 *
 */

BrookIntegrateVerletStepKernel::BrookIntegrateVerletStepKernel( std::string name, const Platform& platform ) :
                                IntegrateVerletStepKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::BrookIntegrateVerletStepKernel";

// ---------------------------------------------------------------------------------------
   
   _brookVerletDynamics   = NULL;
   _brookShakeAlgorithm   = NULL;
}

/** 
 * BrookIntegrateVerletStepKernel destructor
 * 
 */
  
BrookIntegrateVerletStepKernel::~BrookIntegrateVerletStepKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::~BrookIntegrateVerletStepKernel";

// ---------------------------------------------------------------------------------------

   delete _brookVerletDynamics;
   delete _brookShakeAlgorithm;
   
}

/** 
 * Initialize the kernel, setting up all parameters related to integrator.
 * 
 * @param masses             the mass of each atom
 * @param constraintIndices  each element contains the indices of two atoms whose distance should be constrained
 * @param constraintLengths  the required distance between each pair of constrained atoms
 *
 */

void BrookIntegrateVerletStepKernel::initialize( const vector<double>& masses,
                                                 const vector<vector<int> >& constraintIndices,
                                                 const vector<double>& constraintLengths ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::initialize";

// ---------------------------------------------------------------------------------------
   
   _brookVerletDynamics          = new BrookVerletDynamics( );
   _brookVerletDynamics->setup( masses, getPlatform() );

   _brookShakeAlgorithm          = new BrookShakeAlgorithm( );
   _brookShakeAlgorithm->setup( masses, constraintIndices, constraintLengths, getPlatform() );

}

/** 
 * Execute kernel
 * 
 * @param positions          atom coordinates
 * @param velocities         atom velocities
 * @param forces             atom forces
 * @param stepSize           integration step size
 *
 */

void BrookIntegrateVerletStepKernel::execute( Stream& positions, Stream& velocities,
                                              const Stream& forces, double stepSize ){

// ---------------------------------------------------------------------------------------

   double epsilon                           = 1.0e-04;
   static const std::string methodName      = "BrookIntegrateVerletStepKernel::execute";

// ---------------------------------------------------------------------------------------

   // first time through initialize _brookVerletDynamics

   // for each subsequent call, check if parameters need to be updated due to a change
   // in the step size

   // take step

   double difference = stepSize - (double) _brookVerletDynamics->getStepSize();
   if( fabs( difference ) > epsilon ){
//printf( "%s calling updateParameters\n", methodName.c_str() );
      _brookVerletDynamics->updateParameters( stepSize );
   }
   _brookVerletDynamics->update( positions, velocities, forces, *_brookShakeAlgorithm );

}

