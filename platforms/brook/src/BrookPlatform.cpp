/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "BrookPlatform.h"
#include "BrookKernelFactory.h"
//#include "BrookKernels.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;

BrookPlatform* registerBrookPlatform( void ){
   BrookPlatform* platform = new BrookPlatform();
   Platform::registerPlatform(platform);
   return platform;
}

BrookPlatform* staticPlatform = registerBrookPlatform( );

BrookPlatform::BrookPlatform( ){
   _defaultAtomStreamWidth  = DefaultAtomStreamWidth;
   _initializeFactory();
}

BrookPlatform::BrookPlatform( int defaultAtomStreamWidth ){
   _defaultAtomStreamWidth  = defaultAtomStreamWidth;
   _initializeFactory();
}

BrookPlatform::~BrookPlatform( ){
}

void BrookPlatform::_initializeFactory( void ){
   //BrookKernelFactory* factory = new BrookKernelFactory();
   /*
   registerKernelFactory( CalcStandardMMForceFieldKernel::Name(), factory);
   registerKernelFactory( CalcGBSAOBCForceFieldKernel::Name(),    factory);
   registerKernelFactory( IntegrateVerletStepKernel::Name(),      factory);
   registerKernelFactory( IntegrateLangevinStepKernel::Name(),    factory);
   registerKernelFactory( IntegrateBrownianStepKernel::Name(),    factory);
   registerKernelFactory( ApplyAndersenThermostatKernel::Name(),  factory);
   registerKernelFactory( CalcKineticEnergyKernel::Name(),        factory);
   */
}

bool BrookPlatform::supportsDoublePrecision( void ) const {
    return (sizeof(RealOpenMM) >= sizeof(double));
}

const StreamFactory& BrookPlatform::getDefaultStreamFactory( void ) const {
    return defaultStreamFactory;
}

int BrookPlatform::getStreamSize( int size, int streamWidth, int* outputHeight ) const {

   if( streamWidth < 1 ){
      return -1;
   }

   int height = size/streamWidth;
   if( streamWidth*height < size ){
      height++;
   }
   if( outputHeight ){
      *outputHeight = height;
   }
   return height*streamWidth;
}
