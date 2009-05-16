/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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
