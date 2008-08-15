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

#include "BrookRemoveCMMotionKernel.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookRemoveCMMotionKernel constructor
 * 
 * @param name        name of the stream to create
 * @param platform    platform
 *
 */

BrookRemoveCMMotionKernel::BrookRemoveCMMotionKernel( std::string name, const Platform& platform ) :
                              RemoveCMMotionKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::BrookRemoveCMMotionKernel";

// ---------------------------------------------------------------------------------------

}

/** 
 * BrookRemoveCMMotionKernel destructor
 * 
 */
  
BrookRemoveCMMotionKernel::~BrookRemoveCMMotionKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::~BrookRemoveCMMotionKernel";

// ---------------------------------------------------------------------------------------

}

/** 
 * Initialize the kernel
 * 
 * @param masses   array of atom masses
 *
 */

void BrookRemoveCMMotionKernel::initialize( const vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::initialize";

// ---------------------------------------------------------------------------------------

/*
    this->masses.resize(masses.size());
    for (size_t i = 0; i < masses.size(); ++i)
        this->masses[i] = masses[i];
*/

   return;

}

/** 
 * Execute kernel
 * 
 * @param velocities  array of atom velocities
 *
 */

void BrookRemoveCMMotionKernel::execute( Stream& velocities ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookRemoveCMMotionKernel::execute";

// ---------------------------------------------------------------------------------------

/*
    RealOpenMM** velData = ((BrookFloatStreamImpl&) velocities.getImpl()).getData();
    
    // Calculate the center of mass momentum.
    
    RealOpenMM momentum[] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < masses.size(); ++i) {
        momentum[0] += static_cast<RealOpenMM>( masses[i]*velData[i][0] );
        momentum[1] += static_cast<RealOpenMM>( masses[i]*velData[i][1] );
        momentum[2] += static_cast<RealOpenMM>( masses[i]*velData[i][2] );
    }
    
    // Adjust the atom velocities.
    
    momentum[0] /= static_cast<RealOpenMM>( masses.size() );
    momentum[1] /= static_cast<RealOpenMM>( masses.size() );
    momentum[2] /= static_cast<RealOpenMM>( masses.size() );
    for (size_t i = 0; i < masses.size(); ++i) {
        velData[i][0] -= static_cast<RealOpenMM>( momentum[0]/masses[i] );
        velData[i][1] -= static_cast<RealOpenMM>( momentum[1]/masses[i] );
        velData[i][2] -= static_cast<RealOpenMM>( momentum[2]/masses[i] );
    }
*/
   return;
}
