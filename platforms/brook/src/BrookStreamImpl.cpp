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

#include <sstream>
#include "OpenMMException.h"
#include "BrookStreamImpl.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookStreamImpl constructor
 *
 * @param name        stream name
 * @param size        stream size
 * @param streamWidth stream width
 * @param type        StreamImpl data type (float, float2, ...)
 * @param platform    platform reference
 * 
 */

BrookStreamImpl::BrookStreamImpl( const std::string& name, int size, int streamWidth, Stream::DataType type, const Platform& platform ) :
                                  StreamImpl( name, size, type, platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamImpl::BrookStreamImpl";

// ---------------------------------------------------------------------------------------
   
   int isFloat;
   BrookStreamInternal::DataType internalType = getTypeMap( type, &isFloat );

   if( isFloat == 1 || isFloat == 2 ){
      double dangleValue = 0.0;
      _brookStreamInternal = new BrookFloatStreamInternal(  name, size, streamWidth, internalType, dangleValue );
   } else if( isFloat == 0 ){
      int dangleValue = 0;
      _brookStreamInternal = new BrookIntStreamInternal(    name, size, streamWidth, internalType, dangleValue );
   } else {
      std::stringstream message;
      message << methodName << " type=" << type << " for stream=" << name << " is invalid.";
      throw OpenMMException( message.str() );
   }   

}

/** 
 * BrookStreamImpl destructor
 *
 */
BrookStreamImpl::~BrookStreamImpl( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamImpl::~BrookStreamImpl";

// ---------------------------------------------------------------------------------------
   
}

/** 
 * BrookStreamImpl constructor
 *
 * @param type        StreamImpl data type (float, float2, ...)
 * @param isFloat     on output = 1 if float
 *                                2 if double
 *                                0 if integer
 *
 * @return BrookStreamInternal::DataType mapping to Stream::DataType type
 *         if no match, return BrookStreamInternal::Unknown
 * 
 */

BrookStreamInternal::DataType BrookStreamImpl::getTypeMap( Stream::DataType type, int* isFloat ) const { 

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookStreamImpl::getTypeMap";

// ---------------------------------------------------------------------------------------
   
   *isFloat = 0;
   BrookStreamInternal::DataType internalType = BrookStreamInternal::Unknown;

   switch ( type ){

      case Stream::Float:

         internalType = BrookStreamInternal::Float;
         *isFloat      = 1;
         break;

      case Stream::Float2:

         internalType = BrookStreamInternal::Float2;
         *isFloat      = 1;
         break;

      case Stream::Float3:

         internalType = BrookStreamInternal::Float3;
         *isFloat      = 1;
         break;

      case Stream::Float4:

         internalType = BrookStreamInternal::Float4;
         *isFloat      = 1;
         break;

      case Stream::Double:

         internalType = BrookStreamInternal::Double;
         *isFloat      = 2;
         break;

      case Stream::Double2:

         internalType = BrookStreamInternal::Double2;
         *isFloat      = 2;
         break;

      case Stream::Double3:

         internalType = BrookStreamInternal::Double3;
         *isFloat      = 2;
         break;

      case Stream::Double4:

         internalType = BrookStreamInternal::Double4;
         *isFloat      = 2;
         break;

      case Stream::Integer:

         internalType = BrookStreamInternal::Integer;
         *isFloat      = 0;
         break;

      case Stream::Integer2:

         internalType = BrookStreamInternal::Integer2;
         *isFloat      = 0;
         break;

      case Stream::Integer3:

         internalType = BrookStreamInternal::Integer3;
         *isFloat      = 0;
         break;

      case Stream::Integer4:

         internalType = BrookStreamInternal::Integer4;
         *isFloat      = 0;
         break;
   }

   return internalType;
}

/** 
 * Get width
 * 
 * @return width
 */

int BrookStreamImpl::getWidth( void ) const {
   return _brookStreamInternal->getWidth();
}

/** 
 * Get stream width
 * 
 * @return stream width
 */

int BrookStreamImpl::getStreamWidth( void ) const {
   return _brookStreamInternal->getStreamWidth();
}

/** 
 * Get stream height
 * 
 * @return stream height
 */

int BrookStreamImpl::getStreamHeight( void ) const {
   return _brookStreamInternal->getStreamHeight();
}

/** 
 * Get stream size
 * 
 * @return stream size
 */

int BrookStreamImpl::getStreamSize( void ) const {
   return _brookStreamInternal->getStreamSize();
}

/** 
 * Copy the contents of an array into this stream.
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 */
void BrookStreamImpl::loadFromArray( const void* array ){
   return _brookStreamInternal->loadFromArray( array );
}
  
/** 
 * Copy the contents of this stream into an array.
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct data type for this stream.  If the stream has a compound data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 */
void BrookStreamImpl::saveToArray( void* array ){
   return _brookStreamInternal->saveToArray( array );
}
  
/** 
 * Set every element of this stream to the same value.
 * 
 * @param a pointer to the value.  It is assumed to be of the correct data type for this stream.
 */
void BrookStreamImpl::fillWithValue( void* value ){
   return _brookStreamInternal->fillWithValue( value );
}

/** 
 * Get Brook stream
 * 
 * @return Brook stream reference
 */

brook::stream& BrookStreamImpl::getBrookStream( void ){
   return _brookStreamInternal->getBrookStream( );
}

