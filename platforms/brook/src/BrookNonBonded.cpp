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
#include "BrookNonBonded.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include "openmm/OpenMMException.h"
//#include "kernels/invmap.h"
#include "kernels/kforce.h"

using namespace OpenMM;
using namespace std;

/** 
 * Constructor
 * 
 */

BrookNonBonded::BrookNonBonded( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::BrookNonBonded";

// ---------------------------------------------------------------------------------------

   _particleSizeCeiling       = -1;
   _outerUnroll               =  4;
   _innerUnroll               =  4;

   _partialForceStreamWidth   = 64;
   _partialForceStreamHeight  = -1;
   _partialForceStreamSize    = -1;

   _duplicationFactor         =  4;

   _exclusionStreamWidth      = -1;
   _exclusionStreamHeight     = -1;
   _exclusionStreamSize       = -1;

   _jStreamWidth              = 16;
   _jStreamHeight             = -1;
   _jStreamSize               = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _nonbondedStreams[ii] = NULL;
   }

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      _nonbondedForceStreams[ii] = NULL;
   }

}   
 
/** 
 * Destructor
 * 
 */

BrookNonBonded::~BrookNonBonded( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookNonBonded::~BrookNonBonded";
   //static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _nonbondedStreams[ii];
   }

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      delete _nonbondedForceStreams[ii];
   }

}

/** 
 * Get number of force streams
 * 
 * @return  number of force streams (fixed value)
 *
 */

int BrookNonBonded::getNumberOfForceStreams( void ) const {
   return NumberOfForceStreams;
}

/** 
 * Get inner loop unroll
 * 
 * @return   inner loop unroll (fixed value)
 *
 */

int BrookNonBonded::getInnerLoopUnroll( void ) const {
   return _innerUnroll;
}

/** 
 * Get outer loop unroll
 * 
 * @return   outer loop unroll (fixed value)
 *
 */

int BrookNonBonded::getOuterLoopUnroll( void ) const {
   return _outerUnroll;
}

/** 
 * Set outer loop unroll
 * 
 * @param  outer loop unroll (fixed value)
 *
 * @return updated outer loop unroll (fixed value)
 *
 */

int BrookNonBonded::setOuterLoopUnroll( int outerUnroll ){
   if( outerUnroll != _outerUnroll ){
      _particleSizeCeiling = -1;
   }
   _outerUnroll = _outerUnroll;
   return _outerUnroll;
}

/** 
 * Get particle ceiling parameter
 * 
 * @return particle ceiling parameter
 *
 */

int BrookNonBonded::getParticleSizeCeiling( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookNonBonded::getParticleSizeCeiling";

// ---------------------------------------------------------------------------------------

   if( _particleSizeCeiling < 0 ){
      BrookNonBonded* localThis = const_cast<BrookNonBonded* const>(this);
      localThis->_particleSizeCeiling = localThis->getNumberOfParticles() % localThis->getOuterLoopUnroll();
      if( localThis->_particleSizeCeiling ){
         localThis->_particleSizeCeiling = localThis->getOuterLoopUnroll() - localThis->_particleSizeCeiling;
      }   
      localThis->_particleSizeCeiling += localThis->getNumberOfParticles();
   }

   return _particleSizeCeiling;
}

/** 
 * Get duplication factor
 * 
 * @return   duplication factor
 *
 */

int BrookNonBonded::getDuplicationFactor( void ) const {
   return _duplicationFactor;
}

/** 
 * Set duplication factor
 * 
 * @param   duplication factor
 *
 * @return  DefaultReturnValue
 *
 */

int BrookNonBonded::setDuplicationFactor( int duplicationFactor ){
   _duplicationFactor = duplicationFactor;
   return DefaultReturnValue;
}

/** 
 * Get j-stream width
 * 
 * @return  j-stream width
 *
 */

int BrookNonBonded::getJStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getJStreamWidth";

// ---------------------------------------------------------------------------------------

   return _jStreamWidth;
}

/** 
 * Get j-stream height
 * 
 * @return  j-stream height 
 *
 */

int BrookNonBonded::getJStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getJStreamHeight";

// ---------------------------------------------------------------------------------------

   return _jStreamHeight;
}

/** 
 * Get j-stream size
 * 
 * @return  j-stream size
 *
 */

int BrookNonBonded::getJStreamSize( void ) const {
   return _jStreamSize;
}

/** 
 * Get partial force stream width
 * 
 * @param platform  platform
 *
 * @return  partial force stream width
 *
 */

int BrookNonBonded::getPartialForceStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamWidth";

// ---------------------------------------------------------------------------------------

   // get partial force stream width

   if( _partialForceStreamWidth < 0 ){
      //_getPartialForceStreamDimensions( platform );
   }
   return _partialForceStreamWidth;
}

/** 
 * Get partial force stream width
 * 
 * @return  partial force stream width
 *
 */

int BrookNonBonded::getPartialForceStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamWidth";

// ---------------------------------------------------------------------------------------

   return _partialForceStreamWidth;
}

/** 
 * Get partial force stream height
 * 
 * @param platform platform
 *
 * @return  partial force stream height 
 *
 */

int BrookNonBonded::getPartialForceStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamHeight";

// ---------------------------------------------------------------------------------------

   // get partial force stream height

   if( _partialForceStreamHeight < 0 ){
      //_getPartialForceStreamDimensions( platform );
   }
   return _partialForceStreamHeight;
}

/** 
 * Get partial force stream height
 * 
 * @return  partial force stream height 
 *
 */

int BrookNonBonded::getPartialForceStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamHeight";

// ---------------------------------------------------------------------------------------

   return _partialForceStreamHeight;
}

/** 
 * Get partial force stream size
 * 
 * @param platform  platform
 *
 * @return  partial force stream size
 *
 */

int BrookNonBonded::getPartialForceStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamSize";

// ---------------------------------------------------------------------------------------

   // get partial force stream size

   if( _partialForceStreamSize < 0 ){
      //_getPartialForceStreamDimensions( platform );
   }
   return _partialForceStreamSize;
}

/** 
 * Get partial force stream size
 * 
 * @return  partial force stream size
 *
 */

int BrookNonBonded::getPartialForceStreamSize( void ) const {
   return _partialForceStreamSize;
}

/** 
 * Get exclusion stream size
 *
 * @return  exclusion stream size
 *
 */

int BrookNonBonded::getExclusionStreamSize( void ) const {
   return _exclusionStreamSize;
}

/** 
 * Get exclusion stream width
 *
 * @return  exclusion stream width
 *
 */

int BrookNonBonded::getExclusionStreamWidth( void ) const {
   return _exclusionStreamWidth;
}

/** 
 * Get exclusion stream 
 *
 * @return  exclusion stream
 *
 */

BrookFloatStreamInternal* BrookNonBonded::getExclusionStream( void ) const {
   return _nonbondedStreams[ExclusionStream];
}

/** 
 * Get vdw stream 
 *
 * @return  vdw stream
 *
 */

BrookFloatStreamInternal* BrookNonBonded::getOuterVdwStream( void ) const {
   return _nonbondedStreams[OuterVdwStream];
}

/** 
 * Get charge stream 
 *
 * @return  charge stream
 *
 */

BrookFloatStreamInternal* BrookNonBonded::getChargeStream( void ) const {
   return _nonbondedStreams[ChargeStream];
}

/** 
 * Get sigma stream 
 *
 * @return  sigma stream
 *
 */

BrookFloatStreamInternal* BrookNonBonded::getInnerSigmaStream( void ) const {
   return _nonbondedStreams[InnerSigmaStream];
}

/** 
 * Get epsilon stream 
 *
 * @return  epsilon stream
 *
 */

BrookFloatStreamInternal* BrookNonBonded::getInnerEpsilonStream( void ) const {
   return _nonbondedStreams[InnerEpsilonStream];
}

/** 
 * Get force streams 
 *
 * @return  force streams
 *
 */

BrookFloatStreamInternal** BrookNonBonded::getForceStreams( void ){
   return _nonbondedForceStreams;
}

/** 
 * Return true if force[index] stream is set
 *
 * @param    index into force stream
 * @return   true if index is valid && force[index] stream is set; else false
 *
 */

int BrookNonBonded::isForceStreamSet( int index ) const {
   return (index >= 0 && index < getNumberOfForceStreams() && _nonbondedForceStreams[index]) ? 1 : 0;
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfParticles         number of particles
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeStreamSizes";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int particleStreamSize                  = getParticleStreamSize( platform );
   int particleStreamWidth                 = getParticleStreamWidth( platform );

   _initializeExclusionStreamSize( particleStreamSize, particleStreamWidth );
   _initializeJStreamSize( particleStreamSize, particleStreamWidth );
   _initializeOuterVdwStreamSize( particleStreamSize, particleStreamWidth );
   _initializePartialForceStreamSize( particleStreamSize, particleStreamWidth );

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

int BrookNonBonded::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeStreams";
   static const double dangleValue          = 0.0;
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());

   // exclusion

   _nonbondedStreams[ExclusionStream]           = new BrookFloatStreamInternal( BrookCommon::NonBondedExclusionStream,
                                                                               getExclusionStreamSize(), getExclusionStreamWidth(),
                                                                               BrookStreamInternal::Float, dangleValue );

   // outer vdw

   _nonbondedStreams[OuterVdwStream]            = new BrookFloatStreamInternal( BrookCommon::OuterVdwStream, getParticleStreamSize(), 
                                                                                getParticleStreamWidth(), BrookStreamInternal::Float2, dangleValue );

   // inner sigma & epsilon

   _nonbondedStreams[InnerSigmaStream]          = new BrookFloatStreamInternal( BrookCommon::InnerSigmaStream, getJStreamSize(), 
                                                                                getJStreamWidth(), BrookStreamInternal::Float4, dangleValue );

   _nonbondedStreams[InnerEpsilonStream]        = new BrookFloatStreamInternal( BrookCommon::InnerEpsilonStream, getJStreamSize(),
                                                                                getJStreamWidth(), BrookStreamInternal::Float4, dangleValue );

   // charge stream

   _nonbondedStreams[ChargeStream]              = new BrookFloatStreamInternal( BrookCommon::NonBondedChargeStream, getParticleStreamSize(),
                                                                                getParticleStreamWidth(), BrookStreamInternal::Float, dangleValue );

   
   // partial force streams

   std::string partialForceStream = BrookCommon::PartialForceStream;
   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      std::stringstream name;
      name << partialForceStream << ii;
      _nonbondedForceStreams[ii] = new BrookFloatStreamInternal( name.str(), getPartialForceStreamSize(),
                                                                 getPartialForceStreamWidth(), BrookStreamInternal::Float3, dangleValue );
   }

   return DefaultReturnValue;
}

/** 
 * Set exclusion (4x4)
 * 
 * @param i                         particle i index
 * @param j                         particle j index
 * @param exclusionStreamWidth      exclusion stream width
 * @param exclusion                 array of packed exclusions
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_setExclusion( int i, int j, int exclusionStreamWidth, BrookOpenMMFloat* exclusion ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_setExclusion";

// ---------------------------------------------------------------------------------------

   int index    = ( (j/4)* exclusionStreamWidth + i ); // index is the x coordinate!
   int joffset  = j % 4; 

   BrookOpenMMFloat flag;
   switch( joffset ) {

      case 0: 
         flag = (BrookOpenMMFloat) 2.0; 
         break;

      case 1:
         flag = (BrookOpenMMFloat) 3.0; 
         break;

      case 2: 
         flag = (BrookOpenMMFloat) 5.0; 
         break;

      case 3 :  
         flag = (BrookOpenMMFloat) 7.0; 
         break;

   }

   exclusion[index] /= flag;

   return DefaultReturnValue;
}

/** 
 * Initialize exclusions
 * 
 * @param exclusions            vector of sets containing exclusions (1 set entry for every particle)
 * @param platform              Brook platform
 * @param log                   optional Log file reference
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_initializeExclusions( const std::vector<std::set<int> >& exclusionsVector, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeExclusions";
   static const int debug                   = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s exclusions vector size=%d\n", methodName.c_str(), exclusionsVector.size() );
   }

   int exclusionStreamSize = getExclusionStreamSize();
   if( exclusionStreamSize < 1 ){
      std::stringstream message;
      message << methodName << " exclusionStreamSize=" << exclusionStreamSize << " is not set.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   // allocate and initialize exclusions array

   BrookOpenMMFloat* exclusions = new BrookOpenMMFloat[exclusionStreamSize];
   for( unsigned int ii = 0; ii < (unsigned int) exclusionStreamSize; ii++ ){
      exclusions[ii] = (BrookOpenMMFloat) 210.0f;
   }

   // pack in values

   int  exclusionStreamWidth = getExclusionStreamWidth();
   int  numberOfParticles    = getNumberOfParticles();
   int  particleStreamSize       = getParticleStreamSize();
   for( unsigned int ii = 0; ii < exclusionsVector.size(); ii++ ){

      set<int> exclusionIndices  = exclusionsVector[ii];

//(void) fprintf( getLog(), "%s particles=%d excludeSz=%d\n", methodName.c_str(), ii, exclusionIndices.size() );

      int hitII = 0;
      for( set<int>::const_iterator jj = exclusionIndices.begin(); jj != exclusionIndices.end(); jj++ ){
//(void) fprintf( getLog(), "%s particles=%d exclude=%d\n", methodName.c_str(), ii, *jj );
         _setExclusion( ii, *jj, exclusionStreamWidth, exclusions );
         if( *jj == ii )hitII = 1;
      }
      if( !hitII ){
         _setExclusion( ii, ii, exclusionStreamWidth, exclusions );
      }

      // explicitly exclude junk particles from interacting w/ particle ii
      for( int jj = numberOfParticles; jj < particleStreamSize; jj++ ){
         _setExclusion( ii, jj, exclusionStreamWidth, exclusions );
      }
   }

   // explicitly exclude junk particles from interacting w/ all particles

   for( int ii = numberOfParticles; ii < particleStreamSize; ii++ ){
      for( int jj = 0; jj < particleStreamSize; jj++ ){
         _setExclusion( ii, jj, exclusionStreamWidth, exclusions );
      }
   }

   // diagnostics

   if( 1 && debug && getLog() ){
      FILE* log = getLog();
      (void) fprintf( log, "%s particles=%d excl=%d exStrSz=%d w=%d atmStrSz=%d\n", methodName.c_str(), numberOfParticles, exclusionsVector.size(), 
                      exclusionStreamSize, exclusionStreamWidth, particleStreamSize );
/*
      for( int ii = 0; ii < exclusionStreamSize; ii++ ){
int index  = ii/exclusionStreamWidth;
int offset = ii - index*exclusionStreamWidth;
         (void) fprintf( log, "   %6d %6d [%6d %6d] %10.0f\n", ii, offset, index*4, index*4+3, exclusions[ii] );
      }
*/

      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( log, "%6d ", ii );
         int count = 0;
         for( int jj = 0; jj <= numberOfParticles/4; jj++ ){
            int index = jj*exclusionStreamWidth + ii;
            (void) fprintf( log, " [%4d %4d] %5.1f", jj*4, jj*4+3, exclusions[index] );
            int excludeValue = (int) (exclusions[index] + 0.01);
            if( (excludeValue % 2) )count++;
            if( (jj*4+1) < numberOfParticles && (excludeValue % 3) )count++;
            if( (jj*4+2) < numberOfParticles && (excludeValue % 5) )count++;
            if( (jj*4+3) < numberOfParticles && (excludeValue % 7) )count++;
         }
         (void) fprintf( log, " TtlExcl=%d\n", count );
      }
         
      (void) fflush( log );
   }

   // load stream

   _nonbondedStreams[ExclusionStream]->loadFromArray( exclusions ); 
   delete[] exclusions;

   return DefaultReturnValue;
}

/** 
 * Set sigma & epsilon given c6 & c12 (geometric rule)
 * 
 * @param c6                        vdw c6
 * @param c12                       vdw c12
 * @param sigma                     massaged sigma
 * @param epsilon                   massaged epsilon
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_setSigmaEpsilon( double c6, double c12, double* sigma , double* epsilon ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::_setSigmaEpsilon";

// ---------------------------------------------------------------------------------------

   if( fabs( c12 ) < 1.0e-04 ){
      *sigma    = 1.0;
      *epsilon  = 0.0;
   } else {
      //*sigma    = 0.5*pow( c12/c6, 1.0/6.0 );
      //*epsilon  = sqrt( c6*c6/c12 );
      *sigma    = 0.5*c6;
      *epsilon  = 2.0*sqrt( c12 );
   }

   return DefaultReturnValue;
}

/** 
 * Initialize vdw & charge
 * 
 * @param exclusions                vector of sets containing exclusions (1 set entry for every particle)
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_initializeVdwAndCharge( const vector<vector<double> >& nonbondedParameters, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeVdwAndCharge";
   static const int debug                   = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nonbonded vector size=%d\n", methodName.c_str(), nonbondedParameters.size() );
   }

   int particleStreamSize        = getParticleStreamSize();

   // allocate and initialize vdw & charge array

   unsigned int vdwParametersSize  = 2*particleStreamSize;
   BrookOpenMMFloat* vdwParameters = new BrookOpenMMFloat[vdwParametersSize];
   memset( vdwParameters, 0, sizeof( BrookOpenMMFloat )*vdwParametersSize );
   
   BrookOpenMMFloat* charges       = new BrookOpenMMFloat[particleStreamSize];
   memset( charges, 0, sizeof( BrookOpenMMFloat )*particleStreamSize );
   
// ---------------------------------------------------------------------------------------

   // pack in values

   int  numberOfParticles    = getNumberOfParticles();
   int trackingIndex         = 0;
   double sigma, epsilon;
   for( unsigned int ii = 0; ii < nonbondedParameters.size(); ii++ ){

      vector<double> nonbondedParameterVector 
                                     = nonbondedParameters[ii];

      double c6                      = nonbondedParameterVector[1];
      double c12                     = nonbondedParameterVector[2];
      double charge                  = nonbondedParameterVector[0];

/*
      double charge                  = nonbondedParameterVector[0];
      double sigma                   = nonbondedParameterVector[1];
      double epsilon                 = nonbondedParameterVector[2];
*/

      // geometric combination rules

      _setSigmaEpsilon( c6, c12, &sigma, &epsilon );

      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) sigma;
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) epsilon;
      charges[ii]                    = (BrookOpenMMFloat) charge;
         
   }

   // set outlier particles

   for( unsigned int ii = nonbondedParameters.size(); ii < (unsigned int) particleStreamSize; ii++ ){
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) 1.0;
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) 0.0;
   }

// ---------------------------------------------------------------------------------------

   // load stream

   _nonbondedStreams[OuterVdwStream]->loadFromArray( vdwParameters ); 
   _nonbondedStreams[ChargeStream]->loadFromArray( charges ); 

   // diagnostics

   if( 1 && debug && getLog() ){
      FILE* log = getLog();
      int trackingIndex         = 0;
      (void) fprintf( log, "%s particles=%d strSz=%d\n   vdw[Sig Eps] Q\n", methodName.c_str(), numberOfParticles, particleStreamSize );
      for( int ii = 0; ii < particleStreamSize; ii++, trackingIndex += 2 ){
         (void) fprintf( log, "   %d %10.3f %10.3f %10.3f\n", ii, vdwParameters[trackingIndex], vdwParameters[trackingIndex+1], charges[ii] );
      }
      (void) fflush( log );
   }

   delete[] vdwParameters;
   delete[] charges;

   return DefaultReturnValue;
}


/** 
 * Initialize vdw & charge
 * 
 * @param exclusions                vector of sets containing exclusions (1 set entry for every particle)
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::_initializeJStreamVdw( const vector<vector<double> >& nonbondedParameters, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeJStreamVdw";
   static const int debug                   = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nonbonded vector size=%d\n", methodName.c_str(), nonbondedParameters.size() );
   }

   int jStreamSize        = getJStreamSize();

   // allocate and initialize vdw & charge array

   unsigned int jParametersSize        = 4*jStreamSize;

   BrookOpenMMFloat* sigmaParameters   = new BrookOpenMMFloat[jParametersSize];
   //memset( sigmaParameters, 0, sizeof( BrookOpenMMFloat )*jParametersSize );
   
   BrookOpenMMFloat* epsilonParameters = new BrookOpenMMFloat[jParametersSize];
   //memset( epsilonParameters, 0, sizeof( BrookOpenMMFloat )*jParametersSize );
   
// ---------------------------------------------------------------------------------------

   // pack in values

   int  numberOfParticles    = getNumberOfParticles();
   int trackingIndex         = 0;
   double sigma, epsilon;
   for( unsigned int ii = 0; ii < nonbondedParameters.size(); ii++ ){

      vector<double> nonbondedParameterVector 
                                     = nonbondedParameters[ii];

      double c6                      = nonbondedParameterVector[1];
      double c12                     = nonbondedParameterVector[2];

      //double sigma                   = nonbondedParameterVector[1];
      //double epsilon                 = nonbondedParameterVector[2];

      // geometric combination rules

      _setSigmaEpsilon( c6, c12, &sigma, &epsilon );

      sigmaParameters[ii]   = (BrookOpenMMFloat) sigma;
      epsilonParameters[ii] = (BrookOpenMMFloat) epsilon;
         
   }

   // set outlier particles

   for( unsigned int ii = nonbondedParameters.size(); ii < (unsigned int) jParametersSize; ii++ ){
      sigmaParameters[ii]   = (BrookOpenMMFloat) 1.0;
      epsilonParameters[ii] = (BrookOpenMMFloat) 0.0;
   }

// ---------------------------------------------------------------------------------------

   // load streams

   _nonbondedStreams[InnerSigmaStream]->loadFromArray( sigmaParameters ); 
   _nonbondedStreams[InnerEpsilonStream]->loadFromArray( epsilonParameters ); 

   // diagnostics

   if( 1 && debug && getLog() ){
      FILE* log = getLog();
      int trackingIndex         = 0;
      (void) fprintf( log, "%s particles=%d strSz=%d\n   Innervdw[Sig Eps]\n", methodName.c_str(), numberOfParticles, jStreamSize );
      for( unsigned int ii = 0; ii < jParametersSize; ii++ ){
         (void) fprintf( log, "   %d %10.3f %10.3f\n", ii, sigmaParameters[ii], epsilonParameters[ii] );
      }
      (void) fflush( log );
   }

   delete[] sigmaParameters;
   delete[] epsilonParameters;

   return DefaultReturnValue;
}

/* 
 * Setup of nonbonded ixns
 *
 * @param numberOfParticles     number of particles
 * @param nonbondedParameters   vector of nonbonded parameters [particleI][0=c6]
 *                                                             [particleI][1=c12]
 *                                                             [particleI][2=charge]
 * @param platform              Brook platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * */

int BrookNonBonded::setup( int numberOfParticles, const std::vector<std::vector<double> >& nonbondedParameters,
                           const std::vector<std::set<int> >& exclusions, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookNonBonded::setup";

// ---------------------------------------------------------------------------------------

   setNumberOfParticles( numberOfParticles );

   _initializeStreamSizes( numberOfParticles, platform );
   _initializeStreams( platform );

   _initializeExclusions( exclusions, platform );
   _initializeVdwAndCharge( nonbondedParameters, platform );
   _initializeJStreamVdw( nonbondedParameters, platform );

   setIsActive( 1 );

   return DefaultReturnValue;
}

/* 
 * Setup of stream dimensions for exclusion stream
 *
 * @param particleStreamSize        particle stream size
 * @param particleStreamWidth       particle stream width
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * */

int BrookNonBonded::_initializeExclusionStreamSize( int particleStreamSize, int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeExclusionStreamSize";
   //static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _exclusionStreamWidth     = particleStreamSize;
   int innerUnroll           = getInnerLoopUnroll();
   if( innerUnroll < 1 ){
      std::stringstream message;
      message << methodName << " innerUnrolls=" << innerUnroll << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _exclusionStreamHeight    = _exclusionStreamWidth/innerUnroll; 
   _exclusionStreamSize      = _exclusionStreamWidth*_exclusionStreamHeight;

   return DefaultReturnValue;
}

/* 
 * Setup of j-stream dimensions
 *
 * @param particleStreamSize        particle stream size
 * @param particleStreamWidth       particle stream width
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * @throw  OpenMMException  if jStreamWidth < 1 || innerUnroll < 1
 *
 * */

int BrookNonBonded::_initializeJStreamSize( int particleStreamSize, int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializeJStreamSize";

// ---------------------------------------------------------------------------------------

   // validate stream width & inner unroll

   int jStreamWidth                             = getJStreamWidth();
   if( jStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " jStreamWidth=" << jStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   int innerUnroll           = getInnerLoopUnroll();
   if( innerUnroll < 1 ){
      std::stringstream message;
      message << methodName << " innerUnroll=" << innerUnroll << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   // set dimesnions

   _jStreamSize                                 = particleStreamSize/innerUnroll;
   _jStreamHeight                               = _jStreamSize/jStreamWidth; 
   _jStreamHeight                              += ( (_jStreamSize % _jStreamWidth) ? 1 : 0 ); 
   _jStreamSize                                 = jStreamWidth*_jStreamHeight;

   return DefaultReturnValue;
}

/* 
 * Setup of outer vdw stream size
 *
 * @param particleStreamSize        particle stream size
 * @param particleStreamWidth       particle stream width
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * @throw  OpenMMException  if jStreamWidth < 1 || innerUnroll < 1
 *
 * */

int BrookNonBonded::_initializeOuterVdwStreamSize( int particleStreamSize, int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::initializeOuterVdwStreamSize";

// ---------------------------------------------------------------------------------------

   // validate stream width & inner unroll
/*
   int jStreamWidth                             = getJStreamWidth();
   if( jStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " jStreamWidth=" << jStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   int innerUnroll           = getInnerLoopUnroll();
   if( innerUnroll < 1 ){
      std::stringstream message;
      message << methodName << " innerUnroll=" << innerUnroll << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   // set dimesnions

   _jStreamSize                                 = particleStreamSize/innerUnroll;
   _jStreamHeight                               = _jStreamSize/jStreamWidth; 
   _jStreamHeight                              += ( (_jStreamSize % _jStreamWidth) ? 1 : 0 ); 
   _jStreamSize                                 = jStreamWidth*_jStreamHeight;
*/

   return DefaultReturnValue;
}

/* 
 * Setup of stream dimensions for partial force streams
 *
 * @param particleStreamSize        particle stream size
 * @param particleStreamWidth       particle stream width
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * */

int BrookNonBonded::_initializePartialForceStreamSize( int particleStreamSize, int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::_initializePartialForceStreamSize";

// ---------------------------------------------------------------------------------------

   int innerUnroll           = getInnerLoopUnroll();
   if( innerUnroll < 1 ){
      std::stringstream message;
      message << methodName << " innerUnrolls=" << innerUnroll << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   if( _partialForceStreamWidth < 1 ){
      std::stringstream message;
      message << methodName << " partial force stream width=" << _partialForceStreamWidth << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   _partialForceStreamSize    = particleStreamSize*getDuplicationFactor()/innerUnroll;
   _partialForceStreamHeight  = _partialForceStreamSize/_partialForceStreamWidth;
   _partialForceStreamHeight += ( (_partialForceStreamSize % _partialForceStreamWidth) ? 1 : 0);

   return DefaultReturnValue;
}

/* 
 * Get contents of object
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 **/

std::string BrookNonBonded::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfForceStreams() );
   message << _getLine( tab, "Number of force streams:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getDuplicationFactor() );
   message << _getLine( tab, "Duplication factor:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getInnerLoopUnroll () )
   message << _getLine( tab, "Inner loop unroll:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getOuterLoopUnroll() )
   message << _getLine( tab, "Outer loop unroll:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleSizeCeiling() );
   message << _getLine( tab, "Particle ceiling:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getJStreamWidth() );
   message << _getLine( tab, "J-stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getJStreamHeight() );
   message << _getLine( tab, "J-stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getJStreamSize() );
   message << _getLine( tab, "J-stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamWidth() );
   message << _getLine( tab, "Partial force stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamHeight() );
   message << _getLine( tab, "Partial force stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamSize() );
   message << _getLine( tab, "Partial force stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getExclusionStreamWidth() );
   message << _getLine( tab, "Exclusion stream width:", value ); 

   //(void) LOCAL_SPRINTF( value, "%d", getExclusionStreamHeight() );
   //message << _getLine( tab, "Exclusion stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getExclusionStreamSize() );
   message << _getLine( tab, "Exclusion stream size:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 
   message << _getLine( tab, "ExclusionStream:",      (getExclusionStream()    ? Set : NotSet) ); 
   message << _getLine( tab, "VdwStream:",            (getOuterVdwStream()     ? Set : NotSet) ); 
   message << _getLine( tab, "ChargeStream:",         (getChargeStream()       ? Set : NotSet) ); 
   message << _getLine( tab, "SigmaStream:",          (getInnerSigmaStream()   ? Set : NotSet) ); 
   message << _getLine( tab, "EpsilonStream:",        (getInnerEpsilonStream() ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _nonbondedStreams[ii] ){
         message << _nonbondedStreams[ii]->getContentsString( );
      }
   }

   // force streams

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      char description[256];
      (void) LOCAL_SPRINTF( description, "PartialForceStream %d", ii );
      message << _getLine( tab, description,  (isForceStreamSet(ii) ? Set : NotSet) ); 
   }
 
   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      message << std::endl;
      if( _nonbondedForceStreams[ii] ){
         message << _nonbondedForceStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}

/** 
 * Compute forces
 * 
 * @param context OpenMMContextImpl context
 *
 */

void BrookNonBonded::computeForces( BrookStreamImpl& positionStream, BrookStreamImpl& forceStream ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::computeForces";

   static       int printOn                 = 0;
   static const int MaxErrorMessages        = 2;
   static       int ErrorMessages           = 0;

   FILE* log;

// ---------------------------------------------------------------------------------------

    if( printOn && getLog() ){
       log = getLog();
    } else {
       printOn = 0;
    }

   // nonbonded forces

   float epsfac                                      = 138.935485f;
   BrookFloatStreamInternal**  nonbondedForceStreams = getForceStreams();

   knbforce_CDLJ4(
             (float) getNumberOfParticles(),
             (float) getParticleSizeCeiling(),
             (float) getDuplicationFactor(),
             (float) getParticleStreamHeight( ),
             (float) getParticleStreamWidth( ),
             (float) getJStreamWidth( ),
             (float) getPartialForceStreamWidth( ),
             epsfac,
             positionStream.getBrookStream(),
             getChargeStream()->getBrookStream(),
             getOuterVdwStream()->getBrookStream(),
             getInnerSigmaStream()->getBrookStream(),
             getInnerEpsilonStream()->getBrookStream(),
             getExclusionStream()->getBrookStream(),
             nonbondedForceStreams[0]->getBrookStream(),
             nonbondedForceStreams[1]->getBrookStream(),
             nonbondedForceStreams[2]->getBrookStream(),
             nonbondedForceStreams[3]->getBrookStream()
           );

/*
float zerof = 0.0f;
nonbondedForceStreams[0]->fillWithValue( &zerof );
nonbondedForceStreams[1]->fillWithValue( &zerof );
nonbondedForceStreams[2]->fillWithValue( &zerof );
nonbondedForceStreams[3]->fillWithValue( &zerof );
*/

   // diagnostics

   if( printOn ){
   //static int step = 0;
   //if( step++ < 1 ){
      FILE* log = getLog();
      //FILE* log = stderr;
      (void) fprintf( log, "%s\n", methodName.c_str() ); (void) fflush( log );
      (void) fprintf( log, "Post knbforce_CDLJ4: particles=%6d ceiling=%3d dupFac=%3d\n", getNumberOfParticles(),  
                                                                                          getParticleSizeCeiling(),
                                                                                          getDuplicationFactor()  );

      (void) fprintf( log, "                        hght=%6d   width=%3d   jWid=%3d\n", getParticleStreamHeight( ),
                                                                                        getParticleStreamWidth( ),
                                                                                        getJStreamWidth( ) );
      (void) fprintf( log, "                        pFrc=%6d     eps=%12.5e\n",         getPartialForceStreamWidth( ), epsfac );

      (void) fprintf( log, "%s Final NB & bonded forces\n", methodName.c_str() );
      (void) fprintf( log, "Forces here should be zero (pregather)\n" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
      brookStreamInternalF->printToFile( log );
 
      (void) fprintf( log, "Coords\n" );
      brookStreamInternalF   = positionStream.getBrookStreamInternal();
      brookStreamInternalF->printToFile( log );
/*
void* dataArrayV = brookStreamInternalF->getData( 1 );
float* ddd = (float*) dataArrayV;
int idx1 = 0;
int idx2 = 3;
      (void) fprintf( log, "%d %d %.5e\n", idx1/3, idx2/3, computeDistance( idx1, idx2, ddd ) );
*/
 
      (void) fprintf( log, "\nOuterVdwStreamd\n" );
      getOuterVdwStream()->printToFile( log );

      (void) fprintf( log, "\nInnerSigmaStream\n" );
      getInnerSigmaStream()->printToFile( log );

      (void) fprintf( log, "\nInnerEpsilonStream\n" );
      getInnerEpsilonStream()->printToFile( log );

//      (void) fprintf( log, "\nExclusionStream\n" );
//      getExclusionStream()->printToFile( log );

      (void) fprintf( log, "\nChargeStream\n" );
      getChargeStream()->printToFile( log );

      for( int ii = 0; ii < 4; ii++ ){
         (void) fprintf( log, "\nForce stream %d\n", ii );
         nonbondedForceStreams[ii]->printToFile( log );
      }
      (void) fflush( log );
   }

// ---------------------------------------------------------------------------------------

   // gather forces

   kMergeFloat3_4_nobranch( (float) getDuplicationFactor(),
                            (float) getParticleStreamWidth(),
                            (float) getPartialForceStreamWidth(),
                            (float) getNumberOfParticles(),
                            (float) getParticleSizeCeiling(),
                            (float) getOuterLoopUnroll(),
                            nonbondedForceStreams[0]->getBrookStream(),
                            nonbondedForceStreams[1]->getBrookStream(),
                            nonbondedForceStreams[2]->getBrookStream(),
                            nonbondedForceStreams[3]->getBrookStream(),
                            forceStream.getBrookStream() );

   // diagnostics

   if( printOn ){

      (void) fprintf( log, "\n%s NB forces\n", methodName.c_str() );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
      brookStreamInternalF->printToFile( log );
      (void) fprintf( log, "\n%s Done\n", methodName.c_str() ); (void) fflush( log );

   }

   return;
}

/*  
 * Utility to compute distances between two points
 *
 * @param idx1       index of point1 into array 
 * @param idx2       index of point2 into array 
 * @param coords     array of points p1={coord[0], coord[1], coord[2] },
 *                                   p2={coord[3], coord[4], coord[5] },
 *                                   p3={coord[6], coord[7], coord[8] }
 *
 * @return distance
 *
 **/
    

float BrookNonBonded::_computeDistance( int idx1, int idx2, float* coords ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookNonBonded::_computeDistance";

   float d1 = sqrtf( 
                      (coords[idx1]-coords[idx2])*(coords[idx1]-coords[idx2]) +
                      (coords[idx1+1]-coords[idx2+1])*(coords[idx1+1]-coords[idx2+1]) +
                      (coords[idx1+2]-coords[idx2+2])*(coords[idx1+2]-coords[idx2+2]) );
   return d1;

// ---------------------------------------------------------------------------------------
}
