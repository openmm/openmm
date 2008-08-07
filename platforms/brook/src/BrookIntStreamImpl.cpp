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

#include "BrookIntStreamImpl.h"
#include "OpenMMException.h"
#include <sstream>

using namespace OpenMM;

/** 
 * BrookIntStreamImpl constructor
 * 
 * @param name                      stream name
 * @param size                      stream size
 * @param platform                  platform
 *
 */

BrookIntStreamImpl::BrookIntStreamImpl( std::string name, int size, Stream::DataType type, const Platform& platform ) : StreamImpl( name, size, type, platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookIntStreamImpl::BrookIntStreamImpl";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   switch( type ){

      case Stream::Integer:
          width = 1;
          break;

      case Stream::Integer2:
          width = 2;
          break;

      case Stream::Integer3:
          width = 3;
          break;

      case Stream::Integer4:
          width = 4;
          break;

      default:
         std::stringstream message;
         message << methodName << " type=" << type << " not recognized.";
         throw OpenMMException( message.str() );
   }

   data = new int*[size];

   for( int ii = 0; ii < size; ii++ ){
       data[ii] = new int[width];
   }
}

/** 
 * BrookIntStreamImpl destructor
 * 
 */

BrookIntStreamImpl::~BrookIntStreamImpl() {
   delete[] data;
}

/** 
 * Load data from array into stream
 * 
 * @param array                     array to load (length=size*width)
 *
 */

void BrookIntStreamImpl::loadFromArray( const void* array ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamImpl::loadFromArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int* arrayData = (int*) array;
   int index      = 0;
   for( int ii = 0; ii < getSize(); ii++ ){
      for( int jj = 0; jj < width; jj++ ){
         data[ii][jj] = arrayData[index++];
      }
   }
}

/** 
 * Save data from stream to array 
 * 
 * @param array                     array to save data to (length=size*width)
 *
 */

void BrookIntStreamImpl::saveToArray( void* array ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamImpl::saveToArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int* arrayData = (int*) array;
   int index      = 0;
   for( int ii = 0; ii < getSize(); ii++ ){
      for( int jj = 0; jj < width; jj++ ){
         arrayData[index++] = data[ii][jj];
      }
   }
}

/** 
 * Set all stream entries to input value
 * 
 * @param value                     value to load into stream
 *
 */

void BrookIntStreamImpl::fillWithValue( void* value ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamImpl::fillWithValue";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int valueData = *((int*) value);
   for( int ii = 0; ii < getSize(); ii++ ){
      for (int jj = 0; jj < width; jj++ ){
         data[ii][jj] = valueData;
      }
   }
}

/** 
 * Get data
 * 
 * @return data ptr
 *
 */

const int* const * BrookIntStreamImpl::getData( void ) const {
   return data;
}

/** 
 * Get data
 * 
 * @return data ptr
 *
 */

int** BrookIntStreamImpl::getData( void ){
   return data;
}

