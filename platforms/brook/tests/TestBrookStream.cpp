/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the Brook stream.
 */

#include "../../../tests/AssertionUtilities.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include <iostream>
#include <vector>
#include <fstream>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testWriteRead( void ){

   static const int ArraySz = 3000;
   
   BrookPlatform platform;

   // create and initialize arrays

   float* array      = new float[ArraySz];
   float* saveArray  = new float[ArraySz];
   for( int ii = 0; ii < ArraySz; ii++ ){
      array[ii] = (float) ii;
   }
   memset( saveArray, 0, sizeof( float )*ArraySz );

   // get factory & create stream

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());
   StreamImpl* testStream                       = brookStreamFactory.createStreamImpl( "TestStream", ArraySz, Stream::Float, platform );

   // load & retreive data

   testStream->loadFromArray( array );
   testStream->saveToArray( saveArray );

   // test for equality

   for( int ii = 0; ii < ArraySz; ii++ ){
      ASSERT_EQUAL( array[ii], saveArray[ii] );
   }

   delete[] saveArray;
   delete[] array;
}

int main( ){

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
       testWriteRead();
    }  catch( const exception& e ){
      cout << "exception: " << e.what() << endl;
      return 1;
   }   
   cout << "Done" << endl;

   return 0;
}
