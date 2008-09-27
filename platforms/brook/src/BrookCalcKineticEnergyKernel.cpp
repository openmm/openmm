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

#include <sstream>
#include "OpenMMException.h"
#include "BrookCalcKineticEnergyKernel.h"
#include "BrookStreamImpl.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookCalcKineticEnergyKernel constructor
 * 
 * @param name        name of the stream to create
 * @param platform    platform
 *
 */

BrookCalcKineticEnergyKernel::BrookCalcKineticEnergyKernel( std::string name, const Platform& platform ) :
                              CalcKineticEnergyKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::BrookCalcKineticEnergyKernel";

// ---------------------------------------------------------------------------------------

   _masses   = NULL;

}

/** 
 * BrookCalcKineticEnergyKernel destructor
 * 
 */
  
BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel";

// ---------------------------------------------------------------------------------------

   delete[] _masses;
}

/** 
 * Initialize the kernel
 * 
 * @param masses   mass of each atom
 *
 */

void BrookCalcKineticEnergyKernel::initialize( const vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::initialize";

// ---------------------------------------------------------------------------------------

   // masses

   if( _masses ){
      delete[] _masses;
   }

   _masses = new BrookOpenMMFloat[masses.size()];

   for( unsigned int ii = 0; ii < masses.size(); ii++ ){
      _masses[ii]  = static_cast<BrookOpenMMFloat> (masses[ii]);
   }

   return;
}

/** 
 * Execute kernel
 * 
 * @param velocities  stream of atom velocities
 *
 * @return kinetic energy of the system
 *
 */

double BrookCalcKineticEnergyKernel::execute( const Stream& velocities ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcKineticEnergyKernel::execute";

// ---------------------------------------------------------------------------------------

   const BrookStreamImpl& velocityStreamC = dynamic_cast<const BrookStreamImpl&> (velocities.getImpl());
         BrookStreamImpl& velocityStream  = const_cast<BrookStreamImpl&> (velocityStreamC);
   void* dataV                            = velocityStream.getData( );
   float* velocity                        = (float*) dataV;

   double energy                          = 0.0;
   int index                              = 0;

   if( _masses == NULL ){
      std::stringstream message;
      message << methodName << " masses not set.";
      throw OpenMMException( message.str() );
   }    

/*
printf( "   BrookCalcKineticEnergyKernel Masses=%12.5e %12.5e", _masses[0], _masses[1] );
printf( " [%12.5e %12.5e %12.5e]", velocity[index], velocity[index+1], velocity[index+2] );
index += 3;
printf( " [%12.5e %12.5e %12.5e]\n", velocity[index], velocity[index+1], velocity[index+2] );
index = 0;
*/

   for ( int ii = 0; ii < velocityStream.getSize(); ii++, index += 3 ){
      energy += _masses[ii]*(velocity[index]*velocity[index] + velocity[index + 1]*velocity[index + 1] + velocity[index + 2]*velocity[index + 2]);
   }

   return 0.5*energy;
}
