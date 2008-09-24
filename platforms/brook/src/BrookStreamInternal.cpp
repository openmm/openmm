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

int BrookStreamInternal::printToFile( FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStreamInternal::printToFile";

// ---------------------------------------------------------------------------------------

   if( log == NULL ){
      log = stderr;
   }
   std::string contents = getContentsString();
   (void) fprintf( log, "%s\n", contents.c_str() );

   _bodyPrintToFile( log );

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

