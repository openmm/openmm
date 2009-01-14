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

#include "BrookShakeAlgorithm.h"
#include "BrookPlatform.h"
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

   BrookOpenMMFloat zero                    = static_cast<BrookOpenMMFloat>( 0.0 );
   BrookOpenMMFloat one                     = static_cast<BrookOpenMMFloat>( 1.0 );

// ---------------------------------------------------------------------------------------

   _numberOfParticles            = -1;
   _numberOfConstraints          = -1;
   _maxIterations                = 25;
   _shakeTolerance               = static_cast<BrookOpenMMFloat>(1.0e-04);

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
 * Get max iterations
 * 
 * @return  max iterations
 *
 */
         
int BrookShakeAlgorithm::getMaxIterations( void ) const {
   return _maxIterations;
}
         
/** 
 * Set max iterations
 * 
 * @param  max iterations
 *
 * @return DefaultReturnValue
 *
 */
         
int BrookShakeAlgorithm::setMaxIterations( int maxIterations ){
   _maxIterations = maxIterations;
   return DefaultReturnValue;
}
    
/** 
 * Get SHAKE tolerance
 * 
 * @return  SHAKE tolerance  
 *
 */
         
BrookOpenMMFloat BrookShakeAlgorithm::getShakeTolerance( void ) const {
   return _shakeTolerance;
}
         
/** 
 * Set SHAKE tolerance
 * 
 * @param  SHAKE tolerance  
 *
 * @return DefaultReturnValue
 *
 */
         
int BrookShakeAlgorithm::setShakeTolerance( BrookOpenMMFloat tolerance ){
   _shakeTolerance = tolerance;
   return DefaultReturnValue;
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

   BrookOpenMMFloat dangleValue              = static_cast<BrookOpenMMFloat>( 0.0 );

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
 * @param platform              platform reference
 *
 * @return ErrorReturnValue if error
 *
 * @throw OpenMMException if constraintIndices.size() != constraintLengths.size()
 *
 */
    
int BrookShakeAlgorithm::_setShakeStreams( const std::vector<double>& masses, 
                                           const std::vector< std::vector<int> >& constraintIndices,
                                           const std::vector<double>& constraintLengths,
                                           const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   BrookOpenMMFloat one                      = static_cast<BrookOpenMMFloat>( 1.0 );
   BrookOpenMMFloat half                     = static_cast<BrookOpenMMFloat>( 0.5 );

   static const std::string methodName       = "BrookShakeAlgorithm::_updateSdStreams";

// ---------------------------------------------------------------------------------------

   FILE* log = getLog();

   // check that number of constraints for two input vectors is consistent

   if( constraintIndices.size() != constraintLengths.size() ){
      std::stringstream message;
      message << methodName << " constraintIndices size=" << constraintIndices.size() << " does not equal constraintLengths size=" << constraintLengths.size();
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }   

   // get constraint count for each atom

   vector<int> constraintCount( masses.size(), 0);
   std::vector< std::vector<int> >::const_iterator particleIterator = constraintIndices.begin();
   while( particleIterator != constraintIndices.end() ){

      std::vector<int> particleVector = *particleIterator;
      int atomI                       = particleVector[0]; 
      int atomJ                       = particleVector[1]; 

      // check that atom indices are not out of range

      if( atomI < 0 || atomI > (int) masses.size() ){
         std::stringstream message;
         message << methodName << " constraint I-index=" << atomI << " is  out of range [0, " << masses.size() << "]" << std::endl;
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }

      if( atomJ < 0 || atomJ > (int) masses.size() ){
         std::stringstream message;
         message << methodName << " constraint J-index=" << atomJ << " is out of range [0, " << masses.size() << "]" << std::endl;
         throw OpenMMException( message.str() );
         return ErrorReturnValue;
      }

      constraintCount[atomI]++;
      constraintCount[atomJ]++;

      particleIterator++;
   }
   
   // Find clusters consisting of a central atom with up to three peripheral atoms.
   
   particleIterator                                     = constraintIndices.begin();
   std::vector<double>::const_iterator distanceIterator = constraintLengths.begin();
   _numberOfConstraints                                 = static_cast<int>(constraintIndices.size());
   while( particleIterator != constraintIndices.end() ){

      float distance                  = static_cast<float>( *distanceIterator );
      std::vector<int> particleVector = *particleIterator;
      int atomI                       = particleVector[0]; 
      int atomJ                       = particleVector[1]; 

      // Determine the central atom.
       
      bool firstIsCentral;
      if( constraintCount[atomI] > 1 ){
         firstIsCentral = true;
      } else if( constraintCount[atomJ] > 1 ){
         firstIsCentral = false;
      } else if( atomI < atomJ ){
         firstIsCentral = true;
      } else {
         firstIsCentral = false;
      }
      int centralID, peripheralID;
      float centralInvMass, peripheralInvMass;
      if( firstIsCentral ){
         centralID         = atomI;
         peripheralID      = atomJ;
         centralInvMass    = static_cast<float>( masses[atomI] );
         peripheralInvMass = static_cast<float>( masses[atomJ] );
      } else {
         centralID         = atomJ;
         peripheralID      = atomI;
         centralInvMass    = static_cast<float>( masses[atomJ] );
         peripheralInvMass = static_cast<float>( masses[atomI] );
      }
      if( constraintCount[peripheralID] != 1 ){
          throw OpenMMException("Only bonds to hydrogens may be constrained");
      }
       
      // Add it to the cluster
       
      if( _clusters.find(centralID) == _clusters.end() ){
         _clusters[centralID] = ShakeCluster( centralID, 1.0f/centralInvMass );
      }
      _clusters[centralID].addAtom( peripheralID, distance, 1.0f/peripheralInvMass );
 
      particleIterator++;
      distanceIterator++;
   }
    
   // allocate space for streams

   _initializeStreamSizes( (int) masses.size(), (int) _clusters.size(), platform );
   _initializeStreams( platform );

   int shakeParticleStreamSize               = getShakeParticleStreamSize();
   int shakeConstraintStreamSize             = getShakeConstraintStreamSize();

   // allocate arrays to be read down to board and initialize

   BrookOpenMMFloat* particleIndices  = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   BrookOpenMMFloat* shakeParameters  = new BrookOpenMMFloat[4*shakeConstraintStreamSize];
   BrookOpenMMFloat minusOne          = static_cast<BrookOpenMMFloat>(-1.0);
   for( int ii = 0; ii < 4*shakeConstraintStreamSize; ii++ ){
      particleIndices[ii] = minusOne;
   }
   memset( shakeParameters, 0, 4*shakeConstraintStreamSize*sizeof( BrookOpenMMFloat ) ); 

   // load indices & parameters

   int constraintIndex                   = 0;
   for( map<int, ShakeCluster>::const_iterator iter = _clusters.begin(); iter != _clusters.end(); ++iter ){

      const ShakeCluster& cluster        = iter->second;

      particleIndices[constraintIndex]   = static_cast<BrookOpenMMFloat>( cluster._centralID );
      for( int ii = 0; ii < cluster._size; ii++ ){
         particleIndices[constraintIndex+ii+1] = static_cast<BrookOpenMMFloat>( cluster._peripheralID[ii] );
      }

      shakeParameters[constraintIndex]    = static_cast<BrookOpenMMFloat>(  cluster._centralInvMass );
      shakeParameters[constraintIndex+1]  = half/( static_cast<BrookOpenMMFloat>( cluster._centralInvMass + cluster._peripheralInvMass ) );
      shakeParameters[constraintIndex+2]  = static_cast<BrookOpenMMFloat>( cluster._distance*cluster._distance );
      shakeParameters[constraintIndex+3]  = static_cast<BrookOpenMMFloat>(  cluster._peripheralInvMass );

      constraintIndex += 4;
   }

   // write entries to board

   _shakeStreams[ShakeParticleIndicesStream]->loadFromArray( particleIndices );
   _shakeStreams[ShakeParticleParameterStream]->loadFromArray( shakeParameters );

   delete[] shakeParameters;

   // initialize inverse map

   BrookOpenMMFloat* inverseMap = new BrookOpenMMFloat[2*shakeParticleStreamSize];
   for( int ii = 0; ii < shakeParticleStreamSize*2; ii++ ){
      inverseMap[ii] = static_cast<BrookOpenMMFloat>(-1.0);
   }
   
   // build inverse map
 
   for( int ii = 0; ii < shakeConstraintStreamSize; ii++ ){
      int ii4 = ii << 2;
      for( int jj = 0; jj < 4; jj++ ){
         if( particleIndices[ii4+jj] > -1.0f ){
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

   _setShakeStreams( masses, constraintIndices, constraintLengths, platform );

   return DefaultReturnValue;
}

/*  
 * Check constraints
 *
 * @param positions             atom positions
 * @param outputString          output message
 * @param tolerance             tolerance to compare (if < 0, then use algorithm tolerance
 *
 * @return number of errors
 *
 */
    
int BrookShakeAlgorithm::checkConstraints( BrookStreamInternal* positions, std::string& outputString, float tolerance ) const {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookShakeAlgorithm::checkConstraints";
   static const int maxErrorToPrint         = -10;

// ---------------------------------------------------------------------------------------

   void* posArrayV           = positions->getData( 1 );
   const float* posArray     = (float*) posArrayV; 

/*
fprintf( log, "ShakeCoords\n" );
for( int ii = 0; ii < 30; ii++ ){
if( !(ii % 3 ) )fprintf( log, "%d [", ii );
fprintf( log, "%16.7e ", posArray[ii] );
if( !(ii % 2 ) )fprintf( log, "]\n", posArray[ii] );
} */

   // loop over constraints checking if distances are satisfied
   // build up message string if errors are detected
   // also track max difference in case all constraints satisfied
   // to within specified tolerance

   std::stringstream message;
   int errors                 = 0;
   float maxDiff              = -1.0f;
   int maxDiffCentralIndex    = -1;
   int maxDiffPeripheralIndex = -1;
   tolerance                  = tolerance > 0.0f ? tolerance : _shakeTolerance;
   for( map<int, ShakeCluster>::const_iterator iter = _clusters.begin(); iter != _clusters.end(); ++iter ){

      const ShakeCluster& cluster        = iter->second;

      // multiply by 3 to get index into 1d array

      int centralAtomIndex = 3*cluster._centralID;
      for( int ii = 0; ii < cluster._size; ii++ ){

         int peripheralAtomIndex = 3*cluster._peripheralID[ii];
         if( peripheralAtomIndex > -1 && centralAtomIndex > -1 ){

            float distance = sqrtf( 
                               (posArray[centralAtomIndex]   - posArray[peripheralAtomIndex]  )*(posArray[centralAtomIndex]   - posArray[peripheralAtomIndex]  ) +
                               (posArray[centralAtomIndex+1] - posArray[peripheralAtomIndex+1])*(posArray[centralAtomIndex+1] - posArray[peripheralAtomIndex+1]) +
                               (posArray[centralAtomIndex+2] - posArray[peripheralAtomIndex+2])*(posArray[centralAtomIndex+2] - posArray[peripheralAtomIndex+2]) );

            float diff = fabsf( distance - cluster._distance );
            if( diff > tolerance && (errors++ < maxErrorToPrint || maxErrorToPrint < 0) ){
               char value[1024];
#ifdef WIN32
               sprintf_s( value, 1024, "Error: Atom [%6d %6d] d[%16.7e %16.7e] diff=%16.7e [%16.7e %16.7e %16.7e] [%16.7e %16.7e %16.7e]\n",
                          centralAtomIndex/3, peripheralAtomIndex/3,  distance, cluster._distance, diff,
                          posArray[centralAtomIndex], posArray[centralAtomIndex+1], posArray[centralAtomIndex+2],
                          posArray[peripheralAtomIndex], posArray[peripheralAtomIndex+1], posArray[peripheralAtomIndex+2] );
#else
               sprintf( value, "Error: Atom [%6d %6d] diff=%16.7e d[%16.7e %16.7e] [%16.7e %16.7e %16.7e] [%16.7e %16.7e %16.7e]\n",
                          centralAtomIndex/3, peripheralAtomIndex/3, diff, distance, cluster._distance,
                          posArray[centralAtomIndex][0], posArray[centralAtomIndex][1], posArray[centralAtomIndex][2],
                          posArray[peripheralAtomIndex][0], posArray[peripheralAtomIndex][1], posArray[peripheralAtomIndex][2] );
#endif
               message << value;
            }
            if( diff > maxDiff ){
               maxDiffCentralIndex     = centralAtomIndex/3;
               maxDiffPeripheralIndex  = peripheralAtomIndex/3;
               maxDiff                 = diff;
            }
         }
      }

   }

   // report findings

   std::stringstream outputMessage;
   char text[1024];
   if( errors ){

#ifdef WIN32
      (void) sprintf_s( text, 1024, "Shake errors=%d tol=%.3e mxDff=%.3e atoms[%d %d]", errors, tolerance, maxDiff, maxDiffCentralIndex, maxDiffPeripheralIndex );
#else
      (void) sprintf( text, "Shake errors=%d tol=%.3e mxDff=%.3e atoms[%d %d]", errors, tolerance, maxDiff, maxDiffCentralIndex, maxDiffPeripheralIndex );
#endif
      outputMessage << text;

      if( errors >= maxErrorToPrint ){

#ifdef WIN32
         (void) sprintf_s( text, 1024, " only printing first %d errors", maxErrorToPrint );
#else
         (void) sprintf( text, " only printing first %d errors", maxErrorToPrint );
#endif
         outputMessage << text;
      }

      outputMessage << message.str();

   } else {

#ifdef WIN32
      (void) sprintf_s( text, 1024, "Shake no errors: tol=%.3e mxDff=%.3e", tolerance, maxDiff );
#else
      (void) sprintf( text, "Shake no errors: tol=%.3e mxDff=%.3e", tolerance, maxDiff );
#endif
      
      outputMessage << text;
   }

   outputString = outputMessage.str();

   return errors;
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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfConstraints() );
   message << _getLine( tab, "Number of constraints:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfParticles() );
   message << _getLine( tab, "Number of particles:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getMaxIterations() );
   message << _getLine( tab, "Max iterations:", value ); 

   (void) LOCAL_SPRINTF( value, "%2e", getShakeTolerance() );
   message << _getLine( tab, "Tolerance:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamWidth() );
   message << _getLine( tab, "Constraint stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamHeight() );
   message << _getLine( tab, "Constraint stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getShakeConstraintStreamSize() );
   message << _getLine( tab, "Constraint stream size:", value ); 
   message << _getLine( tab, "Constraints:", " " ); 

   int index = 0;
   for( map<int, ShakeCluster>::const_iterator iter = _clusters.begin(); iter != _clusters.end(); ++iter ){

      char buffer[1024];
      const ShakeCluster& cluster        = iter->second;

#ifdef WIN32
      sprintf_s( buffer, 1024, "%5d %6d [%6d %6d %6d] d=%.3f m[%12.5f %12.5f]", index++, cluster._centralID,
                 cluster._peripheralID[0], cluster._peripheralID[1], cluster._peripheralID[2], 
                 cluster._distance, cluster._centralInvMass, cluster._peripheralInvMass );

#endif

      message << _getLine( tab, "   ", buffer ); 
   }
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
