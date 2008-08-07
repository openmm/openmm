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

#include <sstream>
#include "BrookFloatStreamImpl.h"
#include "OpenMMException.h"

using namespace OpenMM;

/** 
 * BrookFloatStreamImpl constructor
 * 
 * @param name                      stream name
 * @param size                      stream size
 * @param type                      stream type (float, float2, ...)
 * @param platform                  platform
 * @param inputStreamWidth          stream width
 * @param inputDefaultDangleValue   default dangle value
 *
 */

BrookFloatStreamImpl::BrookFloatStreamImpl( std::string name, int size, Stream::DataType type,
                                            const Platform& platform, int inputStreamWidth,
                                            double inputDefaultDangleValue = 0.0 ) : StreamImpl( name, size, type, platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamImpl::BrookFloatStreamImpl";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // set base type (currently only FLOAT supported)

   switch( type ){

      case Stream::Float:
      case Stream::Float2:
      case Stream::Float3:
      case Stream::Float4:
   
          _baseType = Stream::Float;
          break;
   
      case Stream::Double:
      case Stream::Double2:
      case Stream::Double3:
      case Stream::Double4:

          _baseType = Stream::Float;
          break;

      default:
         std::stringstream message;
         message << methodName << " stream=" << name << " input type=" << type << " not recognized.";
         throw OpenMMException( message.str() );
         break;
   }

   // set _width (FLOAT, FLOAT2, ... )

   switch( type ){

      case Stream::Float:
      case Stream::Double:

          _width = 1;
          break;

      case Stream::Float2:
      case Stream::Double2:

          _width = 2;
          break;

      case Stream::Float3:
      case Stream::Double3:

          _width = 3;
          break;

      case Stream::Float4:
      case Stream::Double4:

          _width = 4;
          break;
   }

   _defaultDangleValue = (BrookOpenMMFloat) inputDefaultDangleValue;

   // set stream height based on specified stream _width

   if( inputStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " stream=" << name << " input stream width=" << type << " is less than 1.";
      throw OpenMMException( message.str() );
   }

   _streamWidth    = inputStreamWidth;
   _streamHeight   = size/_streamWidth + ((size % _streamWidth) ? 1 : 0);
   int streamSize  = getStreamSize();

   // create Brook stream handle

   switch( _width ){

      case 1:
         _aStream      = brook::stream::create<float>(  _streamHeight, _streamWidth );
         break;
   
      case 2:
         _aStream      = brook::stream::create<float2>( _streamHeight, _streamWidth );
         break;
   
      case 3:
         _aStream      = brook::stream::create<float3>( _streamHeight, _streamWidth );
         break;
   
      case 4:
         _aStream      = brook::stream::create<float4>( _streamHeight, _streamWidth );
         break;
   }

   // allocate memory for data buffer

   _data = new float*[streamSize];
   for( int ii = 0; ii < streamSize; ii++ ){
       _data[ii] = new float[_width];
   }

   // check if we are using float or double

   if( sizeof( float ) != sizeof( RealOpenMM ) ){
   
      _realOpenMMData = new RealOpenMM*[streamSize];
      for( int ii = 0; ii < streamSize; ii++ ){
         _realOpenMMData[ii] = new RealOpenMM[_width];
      }
   } else {
      _realOpenMMData = NULL;
   }
}

/** 
 * BrookFloatStreamImpl destructor
 * 
 */

BrookFloatStreamImpl::~BrookFloatStreamImpl( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookFloatStreamImpl::~BrookFloatStreamImpl";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   //delete _aStream;

   int streamSize = getStreamSize();

   for( int ii = 0; ii < streamSize; ii++ ){
      delete[] _data[ii];
   }
   delete[] _data;

   if( _realOpenMMData ){
      for( int ii = 0; ii < streamSize; ii++ ){
         delete[] _realOpenMMData[ii];
      }
      delete[] _realOpenMMData;
   }

}

/** 
 * Load _data from input array
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 *
 */

void BrookFloatStreamImpl::loadFromArray( const void* array ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamImpl::loadFromArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   if( _baseType == Stream::Float ){
      loadFromArray( array, Stream::Float );
   } else {
      loadFromArray( array, Stream::Double );
   }

}

/** 
 * Get width
 * 
 * @return width
 */

int BrookFloatStreamImpl::getWidth( void ) const {
   return _width;
}

/** 
 * Get stream width
 * 
 * @return stream width
 */

int BrookFloatStreamImpl::getStreamWidth( void ) const {
   return _streamWidth;
}

/** 
 * Get stream height
 * 
 * @return stream height
 */

int BrookFloatStreamImpl::getStreamHeight( void ) const {
   return _streamHeight;
}

/** 
 * Get Brook stream
 * 
 * @return Brook stream 
 */

brook::stream& BrookFloatStreamImpl::getBrookStream( void ){
   return _aStream;
}

/** 
 * Get stream size
 * 
 * @return stream size
 */

int BrookFloatStreamImpl::getStreamSize( void ) const {
   return _streamWidth*_streamHeight;
}

/** 
 * Load _data from input array
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 *
 * @param inputType  type of input array (Stream::Float, Stream::Double, Stream::Integer)
 *
 * @throw if input type not recognized, an exception is thrown
 */

void BrookFloatStreamImpl::loadFromArray( const void* array, Stream::DataType inputType ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamImpl::loadFromArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int index = 0;
   if( inputType == Stream::Float ){

   float* arrayData = (float*) array;
      for( int ii = 0; ii < getSize(); ii++ ){
         for( int jj = 0; jj < getWidth(); jj++ ){
            _data[ii][jj] = (BrookOpenMMFloat) arrayData[index++];
         }
      }

   } else if( inputType == Stream::Double ){

      double* arrayData = (double*) array;
      for( int ii = 0; ii < getSize(); ii++ ){
         for( int jj = 0; jj < getWidth(); jj++ ){
            _data[ii][jj] = (BrookOpenMMFloat) arrayData[index++];
         }
      }

   } else if( inputType == Stream::Integer ){

      int* arrayData = (int*) array;
 
      for( int ii = 0; ii < getSize(); ii++ ){
         for( int jj = 0; jj < getWidth(); jj++ ){
            _data[ii][jj] = (BrookOpenMMFloat) arrayData[index++];
         }
      }

   } else {
      
      std::stringstream message;
      message << methodName << " stream=" << getName() << " input type=" << inputType << " not recognized.";
      throw OpenMMException( message.str() );

   }

   // set dangling values

   _loadDanglingValues();

   // write to GPU

   _aStream.read( _data );
}

void BrookFloatStreamImpl::saveToArray( void* array ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookFloatStreamImpl::saveToArray";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

  // get _data from GPU

  _aStream.write( _data );

  // load into array

  int index = 0;
  if( _baseType == Stream::Float ){

     float* arrayData = (float*) array;
     for( int ii = 0; ii < getSize(); ii++ ){
        for( int jj = 0; jj < getWidth(); jj++ ){
           arrayData[index++] = (float) _data[ii][jj];
        }
     }

   } else {

      double* arrayData = (double*) array;
      for( int ii = 0; ii < getSize(); ii++ ){
         for( int jj = 0; jj < getWidth(); jj++ ){
            arrayData[index++] = _data[ii][jj];
         }
      }

   }
}

/** 
 * Fill stream w/ specified value
 * 
 * @param value                     value to fill stream w/
 *
 */

void BrookFloatStreamImpl::fillWithValue( void* value ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamImpl::fillWithValue";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   BrookOpenMMFloat valueData;
   if( _baseType == Stream::Float) {
      valueData = (BrookOpenMMFloat) *((float*) value);
   } else {
      valueData = (BrookOpenMMFloat) *((double*) value);
   }
   for( int ii = 0; ii < getSize(); ii++ ){
      for (int jj = 0; jj < getWidth(); jj++ ){
         _data[ii][jj] = valueData;
      }
   }

   _loadDanglingValues();
   _aStream.read( _data );
}

const RealOpenMM* const * BrookFloatStreamImpl::getData() const {
   return NULL;
}

// problem w/ const here -- _data is modified (cast away?)

/*
const RealOpenMM* const * BrookFloatStreamImpl::getData() const {

   // retrieve _data from GPU

   _aStream.write( _data );

   // check if RealOpenMM is float; if not, then 
   // copy into realOpenMMData[][] array

   if( realOpenMMData ){
      for( int i = 0; i < getSize(); i++ ){
         for( int j = 0; j < _width; j++ ){
            realOpenMMData[i][j] = (RealOpenMM) _data[i][j];
         }
      }
      return realOpenMMData;
   } else {
      return _data;
   }
}
*/

/** 
 * Get data
 * 
 * @return data array
 *
 */

RealOpenMM** BrookFloatStreamImpl::getData( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamImpl::getData";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _aStream.write( _data );
   if( _realOpenMMData ){
   for( int ii = 0; ii < getSize(); ii++ ){
      for (int jj = 0; jj < getWidth(); jj++ ){
            _realOpenMMData[ii][jj] = (RealOpenMM) _data[ii][jj];
         }
      }
      return _realOpenMMData;
   } else {
      return _data;
   }
}


/** 
 * Load dangling value into stream
 * 
 * @param danglingValue   dangling value to load
 *
 */

void BrookFloatStreamImpl::_loadDanglingValues( float danglingValue ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamImpl::_loadDanglingValues";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int streamSize = getStreamSize();
   for( int ii = getSize(); ii < streamSize; ii++ ){
      for (int jj = 0; jj < getWidth(); jj++ ){
          _data[ii][jj] = danglingValue;
       }
    }
}

/** 
 * Load default dangling value into stream
 * 
 */

void BrookFloatStreamImpl::_loadDanglingValues( void ){
   _loadDanglingValues( _defaultDangleValue );
}
