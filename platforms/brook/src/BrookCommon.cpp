/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs                                                   *
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

#include <math.h>
#include <sstream>
#include "BrookCommon.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include "OpenMMException.h"

using namespace OpenMM;
using namespace std;

/** 
 * Constructor
 * 
 */

BrookCommon::BrookCommon(  ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::BrookCommon";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms             = 0;
   _atomSizeModified          = 0;

   _atomStreamWidth           = -1;
   _atomStreamHeight          = -1; 
   _atomStreamSize            = -1;

   _log                       = NULL;

}   
 
/** 
 * Destructor
 * 
 */

BrookCommon::~BrookCommon( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCommon::~BrookCommon";
   //static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

}

/** 
 * Get number of atoms
 * 
 * @return  number of atoms
 *
 */

int BrookCommon::getNumberOfAtoms( void ) const {
   return _numberOfAtoms;
}

/** 
 * Get number of atoms
 * 
 * @param  numberOfAtoms number of atoms
 * @return  number of atoms
 *
 */

int BrookCommon::setNumberOfAtoms( int numberOfAtoms ){
   if( numberOfAtoms != _numberOfAtoms ){
      _atomSizeModified = numberOfAtoms;
   }
   _numberOfAtoms = numberOfAtoms;
   return _numberOfAtoms;
}

/** 
 * Get atom stream width
 * 
 * @param platform  platform
 *
 * @return  atom stream width
 *
 */

int BrookCommon::getAtomStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream width

   if( _atomStreamWidth < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamWidth;
}

/** 
 * Get atom stream width
 * 
 * @return  atom stream width
 *
 */

int BrookCommon::getAtomStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _atomStreamWidth;
}

/** 
 * Get atom stream height
 * 
 * @param platform platform
 *
 * @return  atom stream height 
 *
 */

int BrookCommon::getAtomStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream height

   if( _atomStreamHeight < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamHeight;
}

/** 
 * Get atom stream height
 * 
 * @return  atom stream height 
 *
 */

int BrookCommon::getAtomStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _atomStreamHeight;
}

/** 
 * Get atom stream size
 * 
 * @param platform  platform
 *
 * @return  atom stream size
 *
 */

int BrookCommon::getAtomStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::getAtomStreamSize";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream size

   if( _atomStreamSize < 0 ){
      _getAtomStreamDimensions( platform );
   }
   return _atomStreamSize;
}

/** 
 * Get atom stream size
 * 
 * @return  atom stream size
 *
 */

int BrookCommon::getAtomStreamSize( void ) const {
   return _atomStreamSize;
}

/** 
 * Get atom stream dimensions
 * 
 * @param platform                  platform
 *
 */

void BrookCommon::_getAtomStreamDimensions( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCommon::_getAtomStreamDimensions";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get atom stream size

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory() );
   BrookStreamInfo* brookStreamInfo             = brookStreamFactory.getBrookStreamInfo( BrookStreamFactory::AtomPositions );
   _atomStreamWidth                             = brookStreamInfo->getStreamWidth();
   _atomStreamSize                              = brookPlatform.getStreamSize( getNumberOfAtoms(), _atomStreamWidth, NULL );
   _atomStreamHeight                            = (int) ( ((float) _atomStreamSize)/( (float) _atomStreamWidth) + 0.001);

   return;
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCommon::getLog( void ) const {
   return _log;
}

/** 
 * Set log file reference
 * 
 * @param  log file reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookCommon::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
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

std::string BrookCommon::_getLine( const std::string& tab, const std::string& description, const std::string& value ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCommon::_getLine";

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
