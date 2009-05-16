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

#include "BrookIntStreamInternal.h"
#include "openmm/OpenMMException.h"
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
 * Get array of appropritate size for loading data
 *
 * @return data array -- user's responsibility to free
 */

void* BrookIntStreamInternal::getDataArray( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookIntStreamInternal::getDataArray";

// ---------------------------------------------------------------------------------------

   int totalSize                          = getStreamSize()*getWidth();
   return new int[totalSize];
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

/** 
 * Get data
 * 
 * @param readFromBoard if set, read values on board 
 *
 * @return data array
 *
 */

void* BrookIntStreamInternal::getData( int readFromBoard ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamInternal::getData";

// ---------------------------------------------------------------------------------------

   if( readFromBoard ){
      _aStream.write( _data );
   }

   return (void*) _data;
}

/* 
 * Print array contents of object to file
 *
 * @param log         file to print to
 *
 * @return DefaultReturnValue
 *
 * */

int BrookIntStreamInternal::_bodyPrintToFile( FILE* log, int maxPrint ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookIntStreamInternal::_bodyPrintToFile";

// ---------------------------------------------------------------------------------------

   void* dataArrayV = getDataArray( );
   saveToArray( dataArrayV );

   int streamSize   = getStreamSize();
   int width        = getWidth();
   int index        = 0;
   int* dataArray   = (int*) dataArrayV;
   for( int ii = 0; ii < streamSize; ii++ ){
      std::stringstream message;
      message.width( 10 );
      message << ii << " [ ";
      for( int jj = 0; jj < width; jj++ ){
         message << dataArray[index++] << " ";
      }   
      message << "]\n";
      (void) fprintf( log, "%s", message.str().c_str() );    
   }   

   delete[] dataArrayV;

   return DefaultReturnValue;

}

/* 
 * Get contents of object
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

const std::string BrookIntStreamInternal::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookIntStreamInternal::getContentsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   //static const char* Set                   = "Set";
   //static const char* NotSet                = "Not set";


// ---------------------------------------------------------------------------------------

   std::stringstream message;
   std::string tab   = "   ";

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   (void) LOCAL_SPRINTF( value, "%s", getName().c_str() );
   message << _getLine( tab, "Name:", value );

   (void) LOCAL_SPRINTF( value, "%d", getWidth() );
   message << _getLine( tab, "Width:", value );

   (void) LOCAL_SPRINTF( value, "%d", getStreamSize() );
   message << _getLine( tab, "Stream size:", value );

   (void) LOCAL_SPRINTF( value, "%d", getStreamWidth() );
   message << _getLine( tab, "Stream width:", value );

   (void) LOCAL_SPRINTF( value, "%d", getStreamHeight() );
   message << _getLine( tab, "Stream height:", value );

   return message.str();
}

/** 
 * BrookFloatStreamInternal constructor
 * 
 * @param stopIndex                 index to stop sum
 * @param sum                       array of size=getWidth()
 *
 * @return DefaultReturnValue
 *
 * @throw exception if stopIndex is too large
 */

int BrookIntStreamInternal::sumByDimension( int stopIndex, double* sum ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookIntStreamInternal::sumByDimension";

// ---------------------------------------------------------------------------------------

   if( stopIndex > getSize() ){
      std::stringstream message;
      message << methodName << " stream=" << getName() << " input topIndex" << stopIndex << " is too large: stream size=" << getSize();
      throw OpenMMException( message.str() );
   }   
   
   // get _data from GPU

   _aStream.write( _data );

   int width                                 = getWidth();
   int widthM1                               = getWidth() - 1;
   stopIndex                                *= width;

   for( int ii = 0; ii < width; ii++ ){
      sum[ii] = 0.0;
   }   

   int index = 0;
   for( int ii = 0; ii < stopIndex; ii++ ){
      sum[index] += (double) _data[ii];
      if( index == widthM1 ){
         index = 0;
      } else {
         index++;
      }   
   }   

   return DefaultReturnValue;
}
