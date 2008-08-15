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

#include "BrookIntegrateKernals.h"
#include "BrookStreamInternal.h"

using namespace OpenMM;
using namespace std;

BrookIntegrateVerletStepKernel::BrookIntegrateVerletStepKernel( std::string name, const Platform& platform ) :
                                IntegrateVerletStepKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::BrookIntegrateVerletStepKernel";

// ---------------------------------------------------------------------------------------
   
}

BrookIntegrateVerletStepKernel::~BrookIntegrateVerletStepKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::~BrookIntegrateVerletStepKernel";

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

void BrookIntegrateVerletStepKernel::initialize( const vector<double>& masses,
                                                 const vector<vector<int> >& constraintIndices,
                                                 const vector<double>& constraintLengths ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::initialize";

// ---------------------------------------------------------------------------------------
   
/*
    this->masses = new RealOpenMM[masses.size()];
    for (size_t i = 0; i < masses.size(); ++i)
        this->masses[i] = static_cast<RealOpenMM>( masses[i] );
    numConstraints = constraintIndices.size();
    this->constraintIndices = allocateIntArray(numConstraints, 2);
    for (int i = 0; i < numConstraints; ++i) {
        this->constraintIndices[i][0] = constraintIndices[i][0];
        this->constraintIndices[i][1] = constraintIndices[i][1];
    }
    shakeParameters = allocateRealArray(constraintLengths.size(), 1);
    for (size_t i = 0; i < constraintLengths.size(); ++i)
        shakeParameters[i][0] = static_cast<RealOpenMM>( constraintLengths[i] );
*/
}

void BrookIntegrateVerletStepKernel::execute( Stream& positions, Stream& velocities,
                                              const Stream& forces, double stepSize ){
// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntegrateVerletStepKernel::execute";

// ---------------------------------------------------------------------------------------
   
/*
    RealOpenMM** posData = ((BrookFloatStreamImpl&) positions.getImpl()).getData();
    RealOpenMM** velData = ((BrookFloatStreamImpl&) velocities.getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((BrookFloatStreamImpl&) forces.getImpl()).getData()); // Brook code needs to be made const correct
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete shake;
        }
        dynamics = new BrookVerletDynamics(positions.getSize(), static_cast<RealOpenMM>(stepSize) );
        shake = new BrookShakeAlgorithm(numConstraints, constraintIndices, shakeParameters);
        dynamics->setBrookShakeAlgorithm(shake);
        prevStepSize = stepSize;
    }
    dynamics->update(positions.getSize(), posData, velData, forceData, masses);
*/
}

