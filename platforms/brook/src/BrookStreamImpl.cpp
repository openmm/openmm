/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
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

#include "openmm/OpenMMException.h"
#include "BrookStreamImpl.h"

#include <sstream>

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
 * Set every element of this stream to the same value.
 * 
 * @return data array
 */
void* BrookStreamImpl::getData( void ){
   return _brookStreamInternal->getData( );
}

/** 
 * Set every element of this stream to the same value.
 *
 * @param readFromBoard if set, read data from board 
 *
 * @return data array
 */
void* BrookStreamImpl::getData( int readFromBoard ){
   return _brookStreamInternal->getData( readFromBoard );
}

/** 
 * Get Brook stream
 * 
 * @return Brook stream reference
 */

brook::stream& BrookStreamImpl::getBrookStream( void ){
   return _brookStreamInternal->getBrookStream( );
}

/** 
 * Get Brook stream impl 
 * 
 * @return Brook stream impl
 */

BrookStreamInternal* BrookStreamImpl::getBrookStreamInternal( void ) const {
   return _brookStreamInternal;
}
