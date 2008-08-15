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

#include "BrookCalcKineticEnergyKernel.h"
#include "BrookStreamInternal.h"

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

}

/** 
 * BrookCalcKineticEnergyKernel destructor
 * 
 */
  
BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::~BrookCalcKineticEnergyKernel";

// ---------------------------------------------------------------------------------------

/*
    if (dynamics)
        delete dynamics;
    if (shake)
        delete shake;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (shakeParameters)
        disposeRealArray(shakeParameters, numConstraints);
*/

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

    // this->masses = masses;

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

   // static const std::string methodName      = "BrookCalcKineticEnergyKernel::execute";

// ---------------------------------------------------------------------------------------

/*
    RealOpenMM** velData = const_cast<RealOpenMM**>(((BrookFloatStreamImpl&) velocities.getImpl()).getData()); // Brook code needs to be made const correct
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        energy += masses[i]*(velData[i][0]*velData[i][0]+velData[i][1]*velData[i][1]+velData[i][2]*velData[i][2]);
    return 0.5*energy;
*/

   return 0.0;
}
