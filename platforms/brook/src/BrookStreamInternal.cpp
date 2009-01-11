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
#include "BrookStreamInternal.h"


using namespace OpenMM;
using namespace std;

BrookStreamInternal::BrookStreamInternal( const std::string& name, int size, int streamWidth, BrookStreamInternal::DataType type ) :
                                           _name(name), _size(size), _streamWidth(streamWidth), _type(type) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::BrookStreamInternal";

// ---------------------------------------------------------------------------------------
   
   _streamHeight   = _size/_streamWidth + ( (_size % _streamWidth) ? 1 : 0);
   _streamSize     = _streamWidth*_streamHeight;
   _baseType       = BrookStreamInternal::Unknown;

//   _aStream        = 0;
}

BrookStreamInternal::~BrookStreamInternal( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStreamInternal::~BrookStreamInternal";

// ---------------------------------------------------------------------------------------

   // _aStream is not ptr
   // delete _aStream;
   
}

/** 
 * Get name 
 * 
 * @return name
 */

const std::string& BrookStreamInternal::getName( void ) const {
   return _name;
}

/** 
 * Get size
 * 
 * @return size
 */

int BrookStreamInternal::getSize( void ) const {
   return _size;
}

/** 
 * Get data type
 * 
 * @return data type
 */

BrookStreamInternal::DataType BrookStreamInternal::getDataType( void ) const {
   return _type;
}
 
/** 
 * Get base data type
 * 
 * @return base data type ( float, double, int )
 */

BrookStreamInternal::DataType BrookStreamInternal::getBaseDataType( void ) const {
   return _baseType;
}
 
/** 
 * Get width
 * 
 * @return width
 */

int BrookStreamInternal::getWidth( void ) const {
   return _width;
}

/** 
 * Get stream width
 * 
 * @return stream width
 */

int BrookStreamInternal::getStreamWidth( void ) const {
   return _streamWidth;
}

/** 
 * Get stream height
 * 
 * @return stream height
 */

int BrookStreamInternal::getStreamHeight( void ) const {
   return _streamHeight;
}

/** 
 * Get Brook stream
 * 
 * @return Brook stream 
 */

brook::stream& BrookStreamInternal::getBrookStream( void ){
   return _aStream;
}

/** 
 * Get stream size
 * 
 * @return stream size
 */

int BrookStreamInternal::getStreamSize( void ) const {
   return _streamSize;
}

/** 
 * Get type string
 *
 * @param type        BrookStreamInternal data type (float, float2, ...)
 *
 * @return string matching type or "Unknown"
 * 
 */

std::string BrookStreamInternal::getTypeString( BrookStreamInternal::DataType type ) const { 

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookStreamImpl::getTypeString";

// ---------------------------------------------------------------------------------------
   
   std::string typeString;

   switch ( type ){

      case BrookStreamInternal::Float:

         typeString = "Float";
         break;

      case BrookStreamInternal::Float2:

         typeString = "Float2";
         break;

      case BrookStreamInternal::Float3:

         typeString = "Float3";
         break;

      case BrookStreamInternal::Float4:

         typeString = "Float4";
         break;

      case BrookStreamInternal::Double:

         typeString = "Double";
         break;

      case BrookStreamInternal::Double2:

         typeString = "Double2";
         break;

      case BrookStreamInternal::Double3:

         typeString = "Double3";
         break;

      case BrookStreamInternal::Double4:

         typeString = "Double4";
         break;

      case BrookStreamInternal::Integer:

         typeString = "Integer";
         break;

      case BrookStreamInternal::Integer2:

         typeString = "Integer2";
         break;

      case BrookStreamInternal::Integer3:

         typeString = "Integer3";
         break;

      case BrookStreamInternal::Integer4:

         typeString = "Integer4";
         break;

      default:

         typeString = "Unknown";
         break;

   }

   return typeString;
}

/* 
 * Get contents of object
 *
 * @param tab         tab
 * @param description description
 * @param value       value
 *
 * @return string containing contents
 *
 * */

std::string BrookStreamInternal::_getLine( const std::string& tab,
                                           const std::string& description,
                                           const std::string& value ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::_getLine";

   static const unsigned int MAX_LINE_CHARS = 256;
   char line[MAX_LINE_CHARS];

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   memset( line, ' ', MAX_LINE_CHARS );  
#ifdef WIN32
   (void) sprintf_s( line, MAX_LINE_CHARS, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#else
   (void) sprintf( line, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#endif
   message << std::string( line ) << std::endl;

   return message.str();

}

/* 
 * Print contents of object to file
 *
 * @param log         file to print to
 *
 * @return DefaultReturnValue
 *
 * */

int BrookStreamInternal::printToFile( FILE* log, int maxPrint ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::printToFile";

// ---------------------------------------------------------------------------------------

   if( log == NULL ){
      log = stderr;
   }
   std::string contents = getContentsString();
   (void) fprintf( log, "%s\n", contents.c_str() );

   _bodyPrintToFile( log, maxPrint );

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

const std::string BrookStreamInternal::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookStreamInternal::getContentsString";

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

/* 
 * Get stats
 *
 * @return statistics vector
 *
 * */

int BrookStreamInternal::getStatistics( std::vector<std::vector<double> >& statistics, int maxScan ){
   return 0;
}

/* 
 * Get stat string
 *
 * @param         tag    id tag
 * @param  statistics    stat vector
 * @return stat string
 *
 * */

std::string BrookStreamInternal::printStatistics( std::string tag, std::vector<std::vector<double> >& statistics ) const {
   
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::getContentsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   std::string tab                          = "   ";

   static int initialized                   = 0;
   static std::vector<std::string> header;

// ---------------------------------------------------------------------------------------

   // row headers

   if( !initialized ){
      initialized = 1;
      header.push_back( "Average" );
      header.push_back( "StdDev" );
      header.push_back( "|Average|" );
      header.push_back( "RawSum" );
      header.push_back( "Min" );
      header.push_back( "MinIndex" );
      header.push_back( "Max" );
      header.push_back( "MaxIndex" );
      header.push_back( "Count" );
   }

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   std::stringstream message;

   // loop over rows (Average, StdDev, ... )
   // building message string: tag + header + tab + values

   for( unsigned int ii = 0; ii < statistics.size(); ii++ ){

      std::vector<double> row = statistics[ii];

      std::stringstream head;
      std::stringstream valueString;
      head << tag << " ";
      if( ii < header.size() ){
         head << header[ii];
      }
      head << " ";

      valueString << "[ ";
      for( unsigned int jj = 0; jj < row.size(); jj++ ){
         (void) LOCAL_SPRINTF( value, "%16.7e ", row[jj] );
         valueString << value;
      }   
      valueString << "]";

      message << _getLine( tab, head.str(), valueString.str() );
   }

   return message.str();
}

/* 
 * Print streams to file
 *
 * @param fileName     file name
 * @param streams      streams to print
 *
 * @return DefaultReturnValue
 *
 * */

int BrookStreamInternal::printStreamsToFile( std::string fileName, std::vector<BrookStreamInternal*>& streams ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::printStreamsToFile";

// ---------------------------------------------------------------------------------------

   FILE* filePtr = fopen( fileName.c_str(), "w" );
   if( !filePtr ){
      (void) fprintf( stderr, "%s could not open file=<%s>\n", methodName.c_str(), fileName.c_str() );
      (void) fflush( stderr );
      return ErrorReturnValue;
   }

   // gather arrays, widths for eah stream, and set index for each stream
   // also set minimum of stream sizes

   int minIndex    = 10000000;
   float** arrays  = new float*[streams.size()];
   float** sums    = new float*[streams.size()];
   int* widths     = new int[streams.size()];
   int* indices    = new int[streams.size()];
   for( unsigned int ii = 0; ii < streams.size(); ii++ ){
      BrookStreamInternal* stream  = streams[ii];
      void* dataArrayV         = stream->getData( 1 );
      arrays[ii]               = (float*) dataArrayV;
      widths[ii]               =  stream->getWidth();
      indices[ii]              = 0;
      sums[ii]                 = new float[4];
      sums[ii][0] = sums[ii][1] = sums[ii][2] = sums[ii][3] = 0.0f;
      if( minIndex > stream->getSize() ){
         minIndex = stream->getSize();
      }
   }

   // sum columns 

   for( int ii = 0; ii < minIndex; ii++ ){
      for( unsigned int kk = 0; kk < streams.size(); kk++ ){
         for( int jj = 0; jj < widths[kk]; jj++ ){
            sums[kk][jj] += arrays[kk][indices[kk]++];
         }
      }
   }

   // reinitialize indices

   for( unsigned int kk = 0; kk < streams.size(); kk++ ){
      indices[kk] = 0;
   }

   // show column sums

   for( unsigned int kk = 0; kk < streams.size(); kk++ ){
      (void) fprintf( filePtr, "Sms " );
      for( int jj = 0; jj < widths[kk]; jj++ ){
         (void) fprintf( filePtr, "%15.5e ", sums[kk][jj] );
      }
   }
   (void) fprintf( filePtr, "\n" );

   for( int ii = 0; ii < minIndex; ii++ ){
      (void) fprintf( filePtr, "%6d ", ii );

      // streams

      for( unsigned int kk = 0; kk < streams.size(); kk++ ){

//         (void) fprintf( filePtr, "[ " );
 
         // ii elements of stream kk

         for( int jj = 0; jj < widths[kk]; jj++ ){
            (void) fprintf( filePtr, "%15.5e ", arrays[kk][indices[kk]++] );
         }

 //        (void) fprintf( filePtr, " ]", ii );
      }
      (void) fprintf( filePtr, "\n", ii );
   }


   // cleanup

   (void) fclose( filePtr );

   delete[] arrays;
   delete[] widths;
   delete[] indices;
   for( int ii = 0; ii < 4; ii++ ){
      delete[] sums[ii];
   }
   delete[] sums;

   return DefaultReturnValue;

}

typedef struct {
  unsigned int type;
  unsigned int dimensions;
  unsigned int dims[4];
} STREAM_HEADER;

#include <winsock.h>
#include <stdarg.h>
#include <limits>
#include <set>

int BrookStreamInternal::loadStreamGivenFileName( std::string& filename ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::loadStreamGivenFileName";

// ---------------------------------------------------------------------------------------

   FILE* log = stderr;
   
   /* open file, read header */

   FILE * filePtr = fopen( filename.c_str(), "rb" );
   if(  filePtr == NULL ){
      (void) fprintf( log, "%s Unable to open/read %s for stream creation\n", methodName.c_str(), filename.c_str() );
      return ErrorReturnValue;
   }

   STREAM_HEADER header;
   if( 1 != fread( &header, sizeof( header ), 1, filePtr ) ){
      (void) fprintf( log, "%s Unable to read %s header for stream %s\n", methodName.c_str(), filename.c_str() );
      (void) fclose( filePtr );
      return ErrorReturnValue;
   }
 
   /* fix endian */

//   header.type = ntohl( header.type );
//  header.dimensions = ntohl( header.dimensions );
/*
   for( int ii = 0; ii < 4; ii++ ){
   //   header.dims[ii] = ntohl( header.dims[ii] );
      (void) fprintf( log, "%s header %d %d for stream %s from %s\n", methodName.c_str(), ii, header.dims[ii], getName().c_str(), filename.c_str() );
   }

   (void) fflush( log );
*/
/*
   if( header.dimensions[0]
     (unsigned int) header.dimensions, (const ::brook::StreamType *) &type,
     (bool) false );
*/
 
   /* ok, load in the data */

/*
   if( getStreamSize()*getWidth() != header.dims[0]*header.dims[1] ){
      (void) fprintf( log, "%s dimension inconsistency for stream %s from %s dim:[%d %d %d %d] sz=%d != sz=%d containerW=%d\n", 
                      methodName.c_str(), getName().c_str(), filename.c_str(),
                      header.dims[0], header.dims[1], header.dims[2], header.dims[3], header.dims[0]*header.dims[1],
                      getStreamSize(), getWidth() );
      (void) fclose( filePtr );
      return ErrorReturnValue;
   }
*/

   // always float

   int bytesToRead   = getStreamSize()*getWidth()*sizeof( float );
   void* dataBuffer  = (void*) malloc( bytesToRead );

   if( dataBuffer == NULL ){
      (void) fprintf( log, "%s Memory Error for file=%s\n", methodName.c_str(), filename.c_str() );
      (void) fclose( filePtr );
      return ErrorReturnValue;
   }

   size_t bytesRead = fread( dataBuffer, bytesToRead, 1, filePtr );
   if( bytesRead != 1 ){
      (void) fprintf( log, "%s Unable to read %d bytes=%d stream from %s\n", methodName.c_str(), bytesRead, bytesToRead, filename.c_str() );
      free( dataBuffer );
      (void) fclose( filePtr );
      return ErrorReturnValue;
   }

   // float/integer case -- nothing more to do but load; if double, then convert float to double

   if( getBaseDataType() == Float || getBaseDataType() == Integer ){
      loadFromArray( dataBuffer );
   } else {
      double* loadBuffer  = (double*) malloc( bytesToRead*2 );
      float* readBuffer   = (float*) dataBuffer;
      for( int ii = 0; ii < getStreamSize()*getWidth(); ii++ ){
         loadBuffer[ii] = (double) readBuffer[ii];
      }
      loadFromArray( loadBuffer );
      free( loadBuffer );
   }
         
   (void) fclose( filePtr );
   free( dataBuffer );
 
   (void) fprintf( log, "%s read %d bytes for stream %s from %s dim:[%d %d %d %d] container=%d %d\n", methodName.c_str(), bytesToRead, getName().c_str(), filename.c_str(),
                   header.dims[0], header.dims[0], header.dims[0], header.dims[0], getStreamSize(), getWidth() );
   (void) fflush( log );

   return DefaultReturnValue;
}
