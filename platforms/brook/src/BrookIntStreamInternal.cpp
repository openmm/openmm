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
 * @param size                      array size
 * @param streamWidth               stream width
 * @param type                      stream type (Integer, Integer2, ...)
 * @param dangleValue               fill value for tail of stream beyond array size
 *
 */

BrookIntStreamInternal::BrookIntStreamInternal( std::string name, int size, int streamWidth,
                                                BrookStreamInternal::DataType type,
                                                int dangleValue ) : 
                        BrookStreamInternal( name, size, streamWidth, type ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookIntStreamInternal::BrookIntStreamInternal";

// ---------------------------------------------------------------------------------------

   _dangleValue = dangleValue;

   switch( type ){

      case BrookStreamInternal::Integer:

          _width = 1;
          break;

      case BrookStreamInternal::Integer2:

          _width = 2;
          break;

      case BrookStreamInternal::Integer3:

          _width = 3;
          break;

      case BrookStreamInternal::Integer4:

          _width = 4;
          break;

      default:

         std::stringstream message;
         message << methodName << " type=" << type << " not recognized.";
         throw OpenMMException( message.str() );
   }

   _data = new int[size*_width];

}

/** 
 * BrookIntStreamInternal destructor
 * 
 */

BrookIntStreamInternal::~BrookIntStreamInternal() {
   delete[] _data;
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

// ---------------------------------------------------------------------------------------

   if( baseType != BrookStreamInternal::Integer ){
      std::stringstream message;
      message << methodName << " stream=" << getName() << " base type=" << getTypeString( baseType ) << " not handled -- add code.";
      throw OpenMMException( message.str() );
   }

   int* arrayData = (int*) array;
   int totalSize  = getSize()*getWidth();

   for( int ii = 0; ii < totalSize; ii++ ){
      _data[ii] = arrayData[ii];
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

// ---------------------------------------------------------------------------------------

   int* arrayData = (int*) array;
   int totalSize  = getSize()*getWidth();
   for( int ii = 0; ii < totalSize; ii++ ){
      arrayData[ii] = _data[ii];
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

// ---------------------------------------------------------------------------------------

   int valueData = *((int*) value);
   int totalSize  = getSize()*getWidth();

   for( int ii = 0; ii < totalSize; ii++ ){
         _data[ii] = valueData;
   }
}

/** 
 * Get data
 * 
 * @return data ptr
 *
 */

void* BrookIntStreamInternal::getData( void ){
   return _data;
}

