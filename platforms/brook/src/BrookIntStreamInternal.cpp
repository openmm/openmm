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

#include "BrookIntStreamInternal.h"
#include "OpenMMException.h"
#include <sstream>

using namespace OpenMM;

/** 
 * BrookIntStreamInternal constructor
 * 
 * @param name                      stream name
 * @param size                      stream size
 * @param platform                  platform
 *
 */

BrookIntStreamInternal::BrookIntStreamInternal( std::string name, int size, int width, BrookStreamInternal::DataType type,
                                                int dangleValue ) : BrookStreamInternal( name, size, width, type ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookIntStreamInternal::BrookIntStreamInternal";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _dangleValue = dangleValue;

   switch( type ){

      case BrookStreamInternal::Integer:
          width = 1;
          break;

      case BrookStreamInternal::Integer2:
          width = 2;
          break;

      case BrookStreamInternal::Integer3:
          width = 3;
          break;

      case BrookStreamInternal::Integer4:
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
 * BrookIntStreamInternal destructor
 * 
 */

BrookIntStreamInternal::~BrookIntStreamInternal() {
   delete[] data;
}

/** 
 * Load data from array into stream
 * 
 * @param array                     array to load (length=size*width)
 *
 */

void BrookIntStreamInternal::loadFromArray( const void* array ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamInternal::loadFromArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return loadFromArray( array, getBaseDataType() );
}

/** 
 * Load data from array into stream
 * 
 * @param array                     array to load (length=size*width)
 *
 */

void BrookIntStreamInternal::loadFromArray( const void* array, BrookStreamInternal::DataType baseType ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookIntStreamInternal::loadFromArray(1)";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   if( baseType != BrookStreamInternal::Integer ){
      std::stringstream message;
      message << methodName << " stream=" << getName() << " base type=" << getTypeString( baseType ) << " not handled -- add code.";
      throw OpenMMException( message.str() );
   }

   int* arrayData = (int*) array;
   int index      = 0;
   for( int ii = 0; ii < getSize(); ii++ ){
      for( int jj = 0; jj < getWidth(); jj++ ){
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

void BrookIntStreamInternal::saveToArray( void* array ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamInternal::saveToArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int* arrayData = (int*) array;
   int index      = 0;
   for( int ii = 0; ii < getSize(); ii++ ){
      for( int jj = 0; jj < getWidth(); jj++ ){
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

void BrookIntStreamInternal::fillWithValue( void* value ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamInternal::fillWithValue";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int valueData = *((int*) value);
   for( int ii = 0; ii < getSize(); ii++ ){
      for (int jj = 0; jj < getWidth(); jj++ ){
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

const int* const * BrookIntStreamInternal::getData( void ) const {
   return data;
}

/** 
 * Get data
 * 
 * @return data ptr
 *
 */

int** BrookIntStreamInternal::getData( void ){
   return data;
}

