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
#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"

using namespace OpenMM;

/** 
 * BrookFloatStreamInternal constructor
 * 
 * @param name                      stream name
 * @param size                      stream size
 * @param type                      stream type (float, float2, ...)
 * @param platform                  platform
 * @param inputStreamWidth          stream width
 * @param inputDefaultDangleValue   default dangle value
 *
 * @throw exception if stream type not recognized or stream width < 1
 *
 */

BrookFloatStreamInternal::BrookFloatStreamInternal( const std::string& name, int size, int streamWidth, BrookStreamInternal::DataType type,
                                                    double inputDefaultDangleValue ) : BrookStreamInternal( name, size, streamWidth, type ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamInternal::BrookFloatStreamInternal";

// ---------------------------------------------------------------------------------------

//fprintf( stderr,"%s %s\n", methodName.c_str(), getName().c_str() );
//fflush( stderr );

   // set base type (currently only FLOAT supported)

   switch( type ){

      case BrookStreamInternal::Float:
      case BrookStreamInternal::Float2:
      case BrookStreamInternal::Float3:
      case BrookStreamInternal::Float4:
   
          _baseType = BrookStreamInternal::Float;
          break;
   
      case BrookStreamInternal::Double:
      case BrookStreamInternal::Double2:
      case BrookStreamInternal::Double3:
      case BrookStreamInternal::Double4:

          _baseType = BrookStreamInternal::Double;
          break;

      default:
         std::stringstream message;
         message << methodName << " stream=" << name << " input type=" << type << " not recognized.";
         throw OpenMMException( message.str() );
         break;
   }

   // set _width (FLOAT, FLOAT2, ... )

   switch( type ){

      case BrookStreamInternal::Float:
      case BrookStreamInternal::Double:

          _width = 1;
          break;

      case BrookStreamInternal::Float2:
      case BrookStreamInternal::Double2:

          _width = 2;
          break;

      case BrookStreamInternal::Float3:
      case BrookStreamInternal::Double3:

          _width = 3;
          break;

      case BrookStreamInternal::Float4:
      case BrookStreamInternal::Double4:

          _width = 4;
          break;
   }

   _defaultDangleValue = (float) inputDefaultDangleValue;

   // set stream height based on specified stream _width

   if( streamWidth < 1 ){

      std::stringstream message;
      message << methodName << " stream=" << name << " input stream width=" << streamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
   }

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

   int streamSize = getStreamSize();

   // allocate memory for data buffer

   _data = new float[streamSize*_width];

//printf( "%s %s data=%u stream=%d [%d %d] width=%d\n", methodName.c_str(), getName().c_str(), (unsigned int) _data, streamSize, _streamHeight, _streamWidth, _width );
//fflush( stdout );

}

/** 
 * BrookFloatStreamInternal destructor
 * 
 */

BrookFloatStreamInternal::~BrookFloatStreamInternal( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamInternal::~BrookFloatStreamInternal";

// ---------------------------------------------------------------------------------------

//printf( "%s %s data=%u stream=%d [%d %d] width=%d\n", methodName.c_str(), getName().c_str(), (unsigned int) _data, getStreamSize(), _streamHeight, _streamWidth, _width );
//fflush( stdout );

   delete[] _data;

}

/** 
 * Get dangle value
 * 
 * @return  dangle value
 */

double BrookFloatStreamInternal::getDangleValue( void ) const {
   return _defaultDangleValue;
}

/** 
 * Load data from input array
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 *
 * @throw exception if baseType not recognized
 *
 */

void BrookFloatStreamInternal::loadFromArray( const void* array ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookFloatStreamInternal::loadFromArray";

// ---------------------------------------------------------------------------------------

   return loadFromArray( array, getBaseDataType() );

}

/** 
 * Load data from input array
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 *
 * @throw exception if baseType not float, double, or integer
 *
 */

void BrookFloatStreamInternal::loadFromArray( const void* array, BrookStreamInternal::DataType baseType ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamInternal::loadFromArray";

// ---------------------------------------------------------------------------------------

   int totalSize                           = getSize()*getWidth();
   //int totalSize                           = getSize();

   if( baseType == BrookStreamInternal::Float ){

      memcpy( _data, array, sizeof( float )*totalSize );
/*
      float* arrayData = (float*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         _data[ii] = (BrookOpenMMFloat) arrayData[ii];
      }
*/

   } else if( baseType == BrookStreamInternal::Double ){

      double* arrayData = (double*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         _data[ii] = (BrookOpenMMFloat) arrayData[ii];
      }

   } else if( baseType == BrookStreamInternal::Integer ){

      int* arrayData = (int*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         _data[ii] = (BrookOpenMMFloat) arrayData[ii];
      }

   } else {
      
      std::stringstream message;
      message << methodName << " stream=" << getName() << " base type=" << getTypeString( baseType ) << " not recognized.";
      throw OpenMMException( message.str() );

   }

   // set dangling values

   _loadDanglingValues();

   // write to GPU

   _aStream.read( _data );
}

/** 
 * Save data to input array
 * 
 * @param  array a pointer to the start of the array.  The array is assumed to have the same length as this stream,
 * and to contain elements of the correct _data type for this stream.  If the stream has a compound _data type, all
 * the values should be packed into a single array: all the values for the first element, followed by all the values
 * for the next element, etc.
 *
 * @throw exception if baseType not float, double, or integer
 *
 */

void BrookFloatStreamInternal::saveToArray( void* array ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamInternal::saveToArray";

// ---------------------------------------------------------------------------------------

   // get _data from GPU

   _aStream.write( _data );

   // load into array

   int totalSize                             = getSize()*getWidth();
   BrookStreamInternal::DataType  baseType   = getBaseDataType();

   if( baseType == BrookStreamInternal::Float ){

//printf( "%s Basetype is float\n", methodName.c_str() );
//fflush( stdout );
      memcpy( array, _data, sizeof( float )*totalSize );
/*
      float* arrayData = (float*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         arrayData[ii] = (float) _data[ii];
      }
*/

   } else if( baseType == BrookStreamInternal::Double ){

//printf( "%s Basetype is double\n", methodName.c_str() );
//fflush( stdout );

      double* arrayData = (double*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         arrayData[ii] = (double) _data[ii];
      }

   } else if( baseType == BrookStreamInternal::Integer ){

//printf( "%s Basetype is int\n", methodName.c_str() );
//fflush( stdout );

      int* arrayData = (int*) array;
      for( int ii = 0; ii < totalSize; ii++ ){
         arrayData[ii] = (int) _data[ii];
      }
   } else {
      
      std::stringstream message;
      message << methodName << " stream=" << getName() << " base type=" << getTypeString( baseType ) << " not recognized.";
      throw OpenMMException( message.str() );

   }
}

/** 
 * Fill stream w/ specified value
 * 
 * @param value                     value to fill stream w/
 *
 */

void BrookFloatStreamInternal::fillWithValue( void* value ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamInternal::fillWithValue";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   BrookOpenMMFloat valueData;
   if( _baseType == BrookStreamInternal::Float) {
      valueData = (BrookOpenMMFloat) *((float*) value);
   } else {
      valueData = (BrookOpenMMFloat) *((double*) value);
   }
   //memset( _data, valueData, sizeof( float )*getSize()*getWidth() );

   int totalSize = getSize()*getWidth();
   for( int ii = 0; ii < totalSize; ii++ ){
      _data[ii] = valueData;
   }

   _loadDanglingValues();
   _aStream.read( _data );

}

/** 
 * Get data
 * 
 * @param readFromBoard if set, read values on board 
 *
 * @return data array
 *
 */

void* BrookFloatStreamInternal::getData( int readFromBoard ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamInternal::getData";

// ---------------------------------------------------------------------------------------

   if( readFromBoard ){
      _aStream.write( _data );
   }

   return (void*) _data;
}

/** 
 * Get data
 * 
 * @return data array
 *
 */

void* BrookFloatStreamInternal::getData( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamInternal::getData";

// ---------------------------------------------------------------------------------------

   return getData( 0 );
}

/** 
 * Load dangling value into stream
 * 
 * @param danglingValue   dangling value to load
 *
 */

void BrookFloatStreamInternal::_loadDanglingValues( float danglingValue ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "BrookFloatStreamInternal::_loadDanglingValues";

// ---------------------------------------------------------------------------------------

   int width      = getWidth();

   int arraySize  = getSize()*width;
   int streamSize = getStreamSize()*width;

//printf( "%s array=%d stream=%d width=%d %s\n", methodName.c_str(), arraySize, streamSize, width, getName().c_str() );
//fflush( stdout );

   if( arraySize < streamSize ){
      for( int ii = arraySize; ii < streamSize; ii++ ){
          _data[ii] = danglingValue;
       }
    }
}

/** 
 * Load default dangling value into stream
 * 
 */

void BrookFloatStreamInternal::_loadDanglingValues( void ){
   _loadDanglingValues( _defaultDangleValue );
}

/* 
 * Get contents of object
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

const std::string BrookFloatStreamInternal::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookFloatStreamInternal::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%3e", getDangleValue() );
   message << _getLine( tab, "Dangle value:", value ); 

   return message.str();
}

/* 
 * Print array contents of object to file
 *
 * @param log         file to print to
 *
 * @return DefaultReturnValue
 *
 * */

int BrookFloatStreamInternal::_bodyPrintToFile( FILE* log, int maxPrint ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamInternal::_bodyPrintToFile";
   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];


// ---------------------------------------------------------------------------------------

   assert( log );

   void* dataArrayV = getData( 1 );

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   int streamSize         = getStreamSize();
   int width              = getWidth();
   int index              = 0;
   maxPrint               = maxPrint < 0 ? streamSize : maxPrint;
   maxPrint              *= width;
   const float* dataArray = (float*) dataArrayV;
   for( int ii = 0; ii < streamSize && ii < maxPrint; ii++ ){

      std::stringstream message;
      (void) LOCAL_SPRINTF( value, "%6d ", ii );
      message << value << " [ ";

      for( unsigned int jj = 0; jj < width; jj++ ){
         (void) LOCAL_SPRINTF( value, "%16.7e ", dataArray[index++] );
         message << value;
      }   
      message << "]\n";
      if( index == (_size+1)*width ){
         (void) fprintf( log, "\n" );
      }
      (void) fprintf( log, "%s", message.str().c_str() );    
   }

   return DefaultReturnValue;

}

/* 
 * Get stats
 *
 * @return statistics vector
 *
 * */

int  BrookFloatStreamInternal::getStatistics( std::vector<std::vector<double> >& statistics, int maxScan ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::getStatistics";
   static const int MinIndex                = 4;
   static const int MinIndexIndex           = 5;
   static const int MaxIndex                = 6;
   static const int MaxIndexIndex           = 7;
   static const int CountIndex              = 8;
   static const int vectorSize              = CountIndex + 1;
   static const double bigValue             = 1.0e+10;

// ---------------------------------------------------------------------------------------

   void* dataArrayV = getData( 1 );

   statistics.resize( vectorSize );
   int streamSize         = getStreamSize();
   int width              = getWidth();
   int index              = 0;
   const float* dataArray = (float*) dataArrayV;
   
   for( int ii = 0; ii < vectorSize; ii++ ){
      for( int jj = 0; jj < width; jj++ ){
         if( ii == MinIndex ){
            statistics[ii].push_back( bigValue );
         } else if( ii == MaxIndex ){
            statistics[ii].push_back( -bigValue );
         } else {
            statistics[ii].push_back( 0.0 );
         }
      }
   }

   for( int ii = 0; ii < streamSize && ii < maxScan; ii++ ){
      for( int jj = 0; jj < width; jj++ ){

         double value                  =  (double) dataArray[index++];

         statistics[0][jj]            += value;
         statistics[1][jj]            += value*value;

         statistics[2][jj]            += fabs( value );

         statistics[CountIndex][jj]   += 1.0;
         if( value < statistics[MinIndex][jj] ){
            statistics[MinIndex][jj]       = value;
            statistics[MinIndexIndex][jj]  = ii;
         }

         if( value > statistics[MaxIndex][jj] ){
            statistics[MaxIndex][jj]       = value;
            statistics[MaxIndexIndex][jj]  = ii;
         }
      }
   }


   for( int jj = 0; jj < width; jj++ ){
      if( statistics[CountIndex][jj] > 0.0 ){
         statistics[3][jj]   = statistics[0][jj];
         statistics[0][jj]  /= statistics[CountIndex][jj];
         statistics[2][jj]  /= statistics[CountIndex][jj];
         statistics[1][jj]   = statistics[1][jj] - statistics[0][jj]*statistics[0][jj]*statistics[CountIndex][jj];
         if( statistics[CountIndex][jj] > 1.0 ){
            statistics[1][jj] = sqrt( statistics[1][jj]/( statistics[CountIndex][jj] - 1.0 ) );
         }
      }
   }

   return DefaultReturnValue;

}

/** 
 * Get array of appropritate size for loading data
 *
 * @return data array -- user's responsibility to free
 */

void* BrookFloatStreamInternal::getDataArray( void ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamInternal::getDataArray";

// ---------------------------------------------------------------------------------------

   int totalSize                          = getStreamSize()*getWidth();
   BrookStreamInternal::DataType baseType = getBaseDataType();
   
   if( baseType == Double || baseType == Double2 ||  baseType == Double3 ||  baseType == Double4 ){
      totalSize *= 2;
   }
   return new float[totalSize];
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

int BrookFloatStreamInternal::sumByDimension( int stopIndex, double* sum ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookFloatStreamInternal::sumByDimension";

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

