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

   _numberOfParticles            = -1;
   _numberOfConstraints          = -1;

   // mark stream dimension variables as unset

   _shakeParticleStreamWidth     = -1;
   _shakeParticleStreamHeight    = -1;
   _shakeParticleStreamSize      = -1;

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
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookShakeAlgorithm::getShakeParticleStreamSize( void ) const {
   return _shakeParticleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookShakeAlgorithm::getShakeParticleStreamWidth( void ) const {
   return _shakeParticleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookShakeAlgorithm::getShakeParticleStreamHeight( void ) const {
   return _shakeParticleStreamHeight;
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
 * Get Shake particle indices stream 
 *
 * @return  Shake particle indices stream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeParticleIndicesStream( void ) const {
   return _shakeStreams[ShakeParticleIndicesStream];
}

/** 
 * Get Shake particle parameter stream
 *
 * @return  Shake particle parameter sStream
 *
 */

BrookFloatStreamInternal* BrookShakeAlgorithm::getShakeParticleParameterStream( void ) const {
   return _shakeStreams[ShakeParticleParameterStream];
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
 * @param numberOfParticles         number of particles
 * @param numberOfConstraints       number of constraints
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookShakeAlgorithm::_initializeStreamSizes( int numberOfParticles, int numberOfConstraints,
                                                 const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookShakeAlgorithm::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _shakeParticleStreamSize           = getParticleStreamSize( platform );
   _shakeParticleStreamWidth          = getParticleStreamWidth( platform );
   _shakeParticleStreamHeight         = getParticleStreamHeight( platform );

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

   int shakeParticleStreamSize               = getShakeParticleStreamSize();
   int shakeParticleStreamWidth              = getShakeParticleStreamWidth();

   int shakeConstraintStreamSize             = getShakeConstraintStreamSize();
   int shakeConstraintStreamWidth            = getShakeConstraintStreamWidth();

    _shakeStreams[ShakeParticleIndicesStream]
                                             = new BrookFloatStreamInternal( BrookCommon::ShakeParticleIndicesStream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeParticleParameterStream]  
                                             = new BrookFloatStreamInternal( BrookCommon::ShakeParticleParameterStream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float4, dangleValue );

    _shakeStreams[ShakeXCons0Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons0Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float3, dangleValue );

    _shakeStreams[ShakeXCons1Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons1Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float3, dangleValue );

    _shakeStreams[ShakeXCons2Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons2Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float3, dangleValue );

    _shakeStreams[ShakeXCons3Stream]         = new BrookFloatStreamInternal( BrookCommon::ShakeXCons3Stream,
                                                                             shakeConstraintStreamSize, shakeConstraintStreamWidth,
                                                                             BrookStreamInternal::Float3, dangleValue );

    _shakeStreams[ShakeInverseMapStream]     = new BrookFloatStreamInternal( BrookCommon::ShakeInverseMapStream,
                                                                             shakeParticleStreamSize, shakeParticleStreamWidth,
                                                                             BrookStreamInternal::Float2, dangleValue );

   return DefaultReturnValue;
}
 
/*  
 * Set Shake streams
 *
 * @param masses                masses
 * @param constraintIndices     constraint particle indices
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

   int shakeParticleStreamSize               = getShakeParticleStreamSize();
   int shakeConstraintStreamSize             = getShakeConstraintStreamSize();

   // check that number of constraints for two input vectors is consistent

   if( constraintIndices.size() != constraintLengths.size() ){
      std::stringstream message;
      message << methodName << " constraintIndices size=" << constraintIndices.size() << " does not equal constraintLengths size=" << constraintLengths.size();
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }   

   // allocate arrays to be read down to board

   BrookOpenMMFloat* particleIndices  = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   BrookOpenMMFloat* shakeParameters  = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   for( int ii = 0; ii < 4*shakeConstraintStreamSize; ii++ ){
      particleIndices[ii] = (BrookOpenMMFloat) -1;
   }
   memset( shakeParameters, 0, 4*shakeConstraintStreamSize*sizeof( BrookOpenMMFloat ) ); 

   std::vector< std::vector<int> >::const_iterator particleIterator = constraintIndices.begin();
   std::vector<double>::const_iterator distanceIterator             = constraintLengths.begin();

   int constraintIndex                = 0;
   _numberOfConstraints               = 0;
   while( particleIterator != constraintIndices.end() ){

      std::vector<int> particleVector = *particleIterator;

      // check that array of indices is not too small or large

      if( particleVector.size() < 2 ){
         std::stringstream message;
         message << methodName << " particleIndices size=" << particleVector.size() << " is too small at constraintIndex=" << (constraintIndex/4);
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }   
   
      if( particleVector.size() > 4 ){
         std::stringstream message;
         message << methodName << " particleIndices size=" << particleVector.size() << " is too large at constraintIndex=" << (constraintIndex/4);
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }   
   
      int index                    = 0;
      int particleIndex1           = -1;
      int particleIndex2           = -1;
      for( std::vector<int>::const_iterator ii = particleVector.begin(); ii != particleVector.end(); ii++, index++ ){
         particleIndices[constraintIndex + index]   = (BrookOpenMMFloat) *ii;
         if( index == 0 ){
            particleIndex1 = *ii;
         } else if( index == 1 ){
            particleIndex2 = *ii;
         }
      }

      // skip negative indices

      if( particleIndex1 < 0 || particleIndex2 < 0 ){
         particleIterator++;
         distanceIterator++;
         continue;
      }

      if( particleIndex1 >= (int) masses.size() || particleIndex2 >= (int) masses.size() ){
         std::stringstream message;
         message << methodName << " constraint indices=[" << particleIndex1 << ", " << particleIndex2 << "] greater than mass array size=" << masses.size();
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }   

      // insure heavy particle is first

      if( masses[particleIndex1] < masses[particleIndex2] ){

         BrookOpenMMFloat swap               = particleIndices[constraintIndex];
         particleIndices[constraintIndex]    = particleIndices[constraintIndex+1];
         particleIndices[constraintIndex+1]  = swap;

         int swapI                           = particleIndex1;
         particleIndex1                      = particleIndex2;
         particleIndex2                      = swapI;
      }

      // set parameters:

      //    (1) 1/(heavy particle mass)
      //    (2) 0.5/(heavy particle mass+hydrogen mass)
      //    (3) constraint distance**2

      shakeParameters[constraintIndex]    = one/( (BrookOpenMMFloat) masses[particleIndex1] );
      shakeParameters[constraintIndex+1]  = half/( (BrookOpenMMFloat) (masses[particleIndex1] + masses[particleIndex2]) );
      shakeParameters[constraintIndex+2]  = (BrookOpenMMFloat) ( (*distanceIterator)*(*distanceIterator) );

      particleIterator++;
      distanceIterator++;
      constraintIndex += 4;
      _numberOfConstraints++;
   }

   // write entries to board

   _shakeStreams[ShakeParticleIndicesStream]->loadFromArray( particleIndices );
   _shakeStreams[ShakeParticleParameterStream]->loadFromArray( shakeParameters );

   delete[] shakeParameters;

   // initialize inverse map

   BrookOpenMMFloat* inverseMap = new BrookOpenMMFloat[2*shakeParticleStreamSize];
   for( int ii = 0; ii < shakeParticleStreamSize*2; ii++ ){
      inverseMap[ii] = -1;
   }
   
   // build inverse map
 
   for( int ii = 0; ii < shakeConstraintStreamSize; ii++ ){
      int ii4 = ii << 2;
      for( int jj = 0; jj < 4; jj++ ){
         if( particleIndices[ii4+jj] != -1 ){
            int particleIndex             = (int) (particleIndices[ii4+jj] + 0.001);
            inverseMap[particleIndex*2]   = (float) ii;
            inverseMap[particleIndex*2+1] = (float) jj;
         }
      }
   }
   
   _shakeStreams[ShakeInverseMapStream]->loadFromArray( inverseMap );

   delete[] particleIndices;
   delete[] inverseMap;

   return DefaultReturnValue;

}
 
/*  
 * Setup of Shake parameters
 *
 * @param masses                masses
 * @param constraintIndices     constraint particle indices
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

   int numberOfParticles  = (int) masses.size();
   setNumberOfParticles( numberOfParticles );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, (int) constraintIndices.size(), platform );
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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfParticles() );
   message << _getLine( tab, "Number of particles:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 

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

   message << _getLine( tab, "ParticleIndices:",      (getShakeParticleIndicesStream()   ? Set : NotSet) ); 
   message << _getLine( tab, "ParticleParameters:",   (getShakeParticleParameterStream() ? Set : NotSet) ); 
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
