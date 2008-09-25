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

#include <sstream>
#include "BrookShakeAlgorithm.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookShakeAlgorithm::BrookShakeAlgorithm( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::BrookShakeAlgorithm";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms                = -1;
   _numberOfConstraints          = -1;

   // mark stream dimension variables as unset

   _shakeAtomStreamWidth         = -1;
   _shakeAtomStreamHeight        = -1;
   _shakeAtomStreamSize          = -1;

   _shakeConstraintStreamSize    = -1;
   _shakeConstraintStreamWidth   = -1;
   _shakeConstraintStreamHeight  = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _shakeStreams[ii]   = NULL;
   }

   _inverseHydrogenMass    = one/( (BrookOpenMMFloat) 1.008);

   // setup inverse sqrt masses

   _inverseSqrtMasses = NULL;

}   
 
/** 
 * Destructor
 * 
 */

BrookShakeAlgorithm::~BrookShakeAlgorithm( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::~BrookShakeAlgorithm";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _shakeStreams[ii];
   }

   delete[] _inverseSqrtMasses;

}

/** 
 * Get number of constraints
 * 
 * @return   number of constraints
 *
 */

int BrookShakeAlgorithm::getNumberOfConstraints( void ) const {
   return _numberOfConstraints;
}

/** 
 * Get inverse of hydrogen mass
 * 
 * @return  inverse of hydrogen mass
 *
 */

BrookOpenMMFloat BrookShakeAlgorithm::getInverseHydrogenMass( void ) const {
   return _inverseHydrogenMass;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookShakeAlgorithm::getShakeAtomStreamSize( void ) const {
   return _shakeAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookShakeAlgorithm::getShakeAtomStreamWidth( void ) const {
   return _shakeAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookShakeAlgorithm::getShakeAtomStreamHeight( void ) const {
   return _shakeAtomStreamHeight;
}

/** 
 * Get Constraint stream size
 *
 * @return  Constraint stream size
 *
 */

int BrookShakeAlgorithm::getShakeConstraintStreamSize( void ) const {
   return _shakeConstraintStreamSize;
}

/** 
 * Get constraint stream width
 *
 * @return  constraint stream width
 *
 */

int BrookShakeAlgorithm::getShakeConstraintStreamWidth( void ) const {
   return _shakeConstraintStreamWidth;
}

/** 
 * Get constraint stream height
 *
 * @return constraint stream height
 */

int BrookShakeAlgorithm::getShakeConstraintStreamHeight( void ) const {
   return _shakeConstraintStreamHeight;
}

/** 
 * Get Shake atom indices stream 
 *
 * @return  Shake atom indices stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeAtomIndicesStream( void ) const {
   return _shakeStreams[ShakeAtomIndicesStream];
}

/** 
 * Get Shake atom parameter stream
 *
 * @return  Shake atom parameter sStream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeAtomParameterStream( void ) const {
   return _shakeStreams[ShakeAtomParameterStream];
}

/** 
 * Get Shake XCons0 stream
 *
 * @return  XCons0 stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeXCons0Stream( void ) const {
   return _shakeStreams[ShakeXCons0Stream];
}

/** 
 * Get Shake XCons1 stream
 *
 * @return  XCons1 stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeXCons1Stream( void ) const {
   return _shakeStreams[ShakeXCons1Stream];
}

/** 
 * Get Shake XCons2 stream
 *
 * @return  XCons2 stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeXCons2Stream( void ) const {
   return _shakeStreams[ShakeXCons2Stream];
}

/** 
 * Get Shake XCons3 stream
 *
 * @return  XCons3 stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeXCons3Stream( void ) const {
   return _shakeStreams[ShakeXCons3Stream];
}

/** 
 * Get Shake inverse map stream
 *
 * @return  Shake inverse map stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeInverseMapStream( void ) const {
   return _shakeStreams[ShakeInverseMapStream];
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfAtoms             number of atoms
 * @param numberOfConstraints       number of constraints
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookShakeAlgorithm::_initializeStreamSizes( int numberOfAtoms, int numberOfConstraints,
                                                 const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _shakeAtomStreamSize           = getAtomStreamSize( platform );
   _shakeAtomStreamWidth          = getAtomStreamWidth( platform );
   _shakeAtomStreamHeight         = getAtomStreamHeight( platform );

   // get constraint stream width & height, and then set stream size to the product

   BrookCommon::getStreamDimensions( numberOfConstraints, &_shakeConstraintStreamWidth, &_shakeConstraintStreamHeight );
   _shakeConstraintStreamSize     = _shakeConstraintStreamWidth*_shakeConstraintStreamHeight;

   return DefaultReturnValue;
}

/** 
 * Initialize streams
 * 
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookShakeAlgorithm::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::_initializeStreams";

   BrookOpenMMFloat dangleValue              = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int shakeAtomStreamSize                   = getShakeAtomStreamSize();
   int shakeAtomStreamWidth                  = getShakeAtomStreamWidth();

   int shakeConstraintStreamSize             = getShakeConstraintStreamSize();
   int shakeConstraintStreamWidth            = getShakeConstraintStreamWidth();

    _shakeStreams[ShakeAtomIndicesStream]    = new BrookFloatStreamInternal( BrookCommon::ShakeAtomIndicesStream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeAtomParameterStream]  = new BrookFloatStreamInternal( BrookCommon::ShakeAtomParameterStream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeXCons0Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons0Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeXCons1Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons1Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeXCons2Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons2Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeXCons3Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons3Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeInverseMapStream]     = new BrookFloatStreamInternal( BrookCommon::ShakeInverseMapStream,
                                                                             shakeAtomStreamSize, shakeAtomStreamWidth,
                                                                             BrookStreamInternal::Float2, dangleValue );

   return DefaultReturnValue;
}
 
/*  
 * Set Shake streams
 *
 * @param masses                masses
 * @param constraintIndices     constraint atom indices
 * @param constraintLengths     constraint lengths
 *
 * @return ErrorReturnValue if error
 *
 * @throw OpenMMException if constraintIndices.size() != constraintLengths.size()
 *
 */
    
int BrookShakeAlgorithm::_setShakeStreams( const std::vector<double>& masses, const std::vector< std::vector<int> >& constraintIndices,
                                           const std::vector<double>& constraintLengths ){
    
// ---------------------------------------------------------------------------------------

   BrookOpenMMFloat one                      = (BrookOpenMMFloat)  1.0;
   BrookOpenMMFloat half                     = (BrookOpenMMFloat)  0.5;

   static const std::string methodName       = "BrookShakeAlgorithm::_updateSdStreams";

// ---------------------------------------------------------------------------------------

   int shakeAtomStreamSize                   = getShakeAtomStreamSize();
   int shakeConstraintStreamSize             = getShakeConstraintStreamSize();

   // check that number of constraints for two input vectors is consistent

   if( constraintIndices.size() != constraintLengths.size() ){
      std::stringstream message;
      message << methodName << " constraintIndices size=" << constraintIndices.size() << " does not equal constraintLengths size=" << constraintLengths.size();
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }   

   // allocate arrays to be read down to board

   BrookOpenMMFloat* atomIndices      = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   BrookOpenMMFloat* shakeParameters  = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   for( int ii = 0; ii < 4*shakeConstraintStreamSize; ii++ ){
      atomIndices[ii] = (BrookOpenMMFloat) -1;
   }
   memset( shakeParameters, 0, 4*shakeConstraintStreamSize*sizeof( BrookOpenMMFloat ) ); 

   std::vector< std::vector<int> >::const_iterator atomIterator = constraintIndices.begin();
   std::vector<double>::const_iterator distanceIterator         = constraintLengths.begin();

   int constraintIndex                = 0;
   while( atomIterator != constraintIndices.end() ){

      std::vector<int> atomVector = *atomIterator;

      // check that array of indices is not too small or large

      if( atomVector.size() < 2 ){
         std::stringstream message;
         message << methodName << " atomIndices size=" << atomVector.size() << " is too small at constraintIndex=" << (constraintIndex/4);
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }   
   
      if( atomVector.size() > 4 ){
         std::stringstream message;
         message << methodName << " atomIndices size=" << atomVector.size() << " is too large at constraintIndex=" << (constraintIndex/4);
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }   
   
      int index                    = 0;
      int atomIndex1               = -1;
      int atomIndex2               = -1;
      for( std::vector<int>::const_iterator ii = atomVector.begin(); ii != atomVector.end(); ii++, index++ ){
         atomIndices[constraintIndex + index]   = (BrookOpenMMFloat) *ii;
         if( index == 0 ){
            atomIndex1 = *ii;
         } else if( index == 1 ){
            atomIndex2 = *ii;
         }
      }

      // insure heavy atom is first

      if( masses[atomIndex1] < masses[atomIndex2] ){

         BrookOpenMMFloat swap           = atomIndices[constraintIndex];
         atomIndices[constraintIndex]    = atomIndices[constraintIndex+1];
         atomIndices[constraintIndex+1]  = swap;

         int swapI                       = atomIndex1;
         atomIndex1                      = atomIndex2;
         atomIndex2                      = swapI;
      }

      // set parameters:

      //    (1) 1/(heavy atom mass)
      //    (2) 0.5/(heavy atom mass+hydrogen mass)
      //    (3) constraint distance**2

      shakeParameters[constraintIndex]    = one/( (BrookOpenMMFloat) masses[atomIndex1] );
      shakeParameters[constraintIndex+1]  = half/( (BrookOpenMMFloat) (masses[atomIndex1] + masses[atomIndex2]) );
      shakeParameters[constraintIndex+2]  = (BrookOpenMMFloat) ( (*distanceIterator)*(*distanceIterator) );

      atomIterator++;
      distanceIterator++;
      constraintIndex += 4;
   }

   // write entries to board

   _shakeStreams[ShakeAtomIndicesStream]->loadFromArray( atomIndices );
   _shakeStreams[ShakeAtomParameterStream]->loadFromArray( shakeParameters );

   delete[] shakeParameters;

   // initialize inverse map

   BrookOpenMMFloat* inverseMap = new BrookOpenMMFloat[2*shakeAtomStreamSize];
   for( int ii = 0; ii < shakeAtomStreamSize*2; ii++ ){
      inverseMap[ii] = -1;
   }
   
   // build inverse map
 
   for( int ii = 0; ii < shakeConstraintStreamSize; ii++ ){
      int ii4 = ii << 2;
      for( int jj = 0; jj < 4; jj++ ){
         if( atomIndices[ii4+jj] != -1 ){
            int atomIndex             = (int) (atomIndices[ii4+jj] + 0.001);
            inverseMap[atomIndex*2]   = (float) ii;
            inverseMap[atomIndex*2+1] = (float) jj;
         }
      }
   }
   
   _shakeStreams[ShakeInverseMapStream]->loadFromArray( inverseMap );

   delete[] atomIndices;
   delete[] inverseMap;

   return DefaultReturnValue;

}
 
/*  
 * Setup of Shake parameters
 *
 * @param masses                masses
 * @param constraintIndices     constraint atom indices
 * @param constraintLengths     constraint lengths
 * @param platform              Brook platform
 *
 * @return ErrorReturnValue if error
 *
 */
    
int BrookShakeAlgorithm::setup( const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
                                const std::vector<double>& constraintLengths, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::setup";

// ---------------------------------------------------------------------------------------

   int numberOfAtoms  = (int) masses.size();
   setNumberOfAtoms( numberOfAtoms );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfAtoms, (int) constraintIndices.size(), platform );
   _initializeStreams( platform );

   _setShakeStreams( masses, constraintIndices, constraintLengths );

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

std::string BrookShakeAlgorithm::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookShakeAlgorithm::getContentsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   static const char* Set                   = "Set";
   static const char* NotSet                = "Not set";

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   std::string tab   = "   ";

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfConstraints() );
   message << _getLine( tab, "Number of constraints:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamWidth() );
   message << _getLine( tab, "Constraint stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamHeight() );
   message << _getLine( tab, "Constraint stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamSize() );
   message << _getLine( tab, "Constraint stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", 1.0f/getInverseHydrogenMass() );
   message << _getLine( tab, "H-mass:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                          ? Set : NotSet) ); 

   message << _getLine( tab, "AtomIndices:",          (getShakeAtomIndicesStream()       ? Set : NotSet) ); 
   message << _getLine( tab, "AtomParameters:",       (getShakeAtomParameterStream()     ? Set : NotSet) ); 
   message << _getLine( tab, "XCons0:",               (getShakeXCons0Stream()            ? Set : NotSet) ); 
   message << _getLine( tab, "XCons1:",               (getShakeXCons1Stream()            ? Set : NotSet) ); 
   message << _getLine( tab, "XCons2:",               (getShakeXCons2Stream()            ? Set : NotSet) ); 
   message << _getLine( tab, "XCons3:",               (getShakeXCons3Stream()            ? Set : NotSet) ); 
   message << _getLine( tab, "InverseMap:",           (getShakeInverseMapStream()        ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _shakeStreams[ii] ){
         message << _shakeStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
