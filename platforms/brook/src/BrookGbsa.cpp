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
#include "BrookGbsa.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"
#include "gpu/kgbsa.h"
#include "gpu/kforce.h"

using namespace OpenMM;
using namespace std;

/** 
 * Constructor
 * 
 */

BrookGbsa::BrookGbsa(  ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::BrookGbsa";

// ---------------------------------------------------------------------------------------

   _particleSizeCeiling       = -1;
   _outerUnroll               =  4;
   _innerUnroll               =  4;

   _partialForceStreamWidth   = 64;
   _partialForceStreamHeight  = -1;
   _partialForceStreamSize    = -1;

   _gbsaParticleStreamWidth   = -1;
   _gbsaParticleStreamHeight  = -1;
   _gbsaParticleStreamSize    = -1;

   _duplicationFactor         =  4;

   _includeAce                = 1;
   _solventDielectric         = 78.3;
   _soluteDielectric          = 1.0;
   _dielectricOffset          = 0.009;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _gbsaStreams[ii] = NULL;
   }

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      _gbsaForceStreams[ii] = NULL;
   }

   _bornRadiiInitialized      = 0;
   _cpuObc                    = NULL;
   _charges                   = NULL;

}   
 
/** 
 * Destructor
 * 
 */

BrookGbsa::~BrookGbsa( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookGbsa::~BrookGbsa";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _gbsaStreams[ii];
   }

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      delete _gbsaForceStreams[ii];
   }

   delete _cpuObc;

   delete[] _charges;

}

/** 
 * Get number of force streams
 * 
 * @return  number of force streams (fixed value)
 *
 */

int BrookGbsa::getNumberOfForceStreams( void ) const {
   return NumberOfForceStreams;
}


/** 
 * Include ACE approximation in calculation of force
 * 
 * @return true if ACE approximation is to be included in calculation of force
 *
 */

int BrookGbsa::includeAce( void ) const {
   return _includeAce;
}

/** 
 * Get inner loop unroll
 * 
 * @return   inner loop unroll (fixed value)
 *
 */

int BrookGbsa::getInnerLoopUnroll( void ) const {
   return _innerUnroll;
}

/** 
 * Get outer loop unroll
 * 
 * @return   outer loop unroll (fixed value)
 *
 */

int BrookGbsa::getOuterLoopUnroll( void ) const {
   return _outerUnroll;
}

/** 
 * Get solute dielectric
 * 
 * @return   solute dielectric
 *
 */

float BrookGbsa::getSoluteDielectric( void ) const {
   return (float) _soluteDielectric;
}

/** 
 * Get solvent dielectric
 * 
 * @return   solvent dielectric
 *
 */

float BrookGbsa::getSolventDielectric( void ) const {
   return (float) _solventDielectric;
}

/** 
 * Get OBC dielectric offset
 * 
 * @return  OBC dielectric offset
 *
 */

float BrookGbsa::getDielectricOffset( void ) const {
   return (float) _dielectricOffset;
}

/** 
 * Set outer loop unroll
 * 
 * @param  outer loop unroll (fixed value)
 *
 * @return updated outer loop unroll (fixed value)
 *
 */

int BrookGbsa::setOuterLoopUnroll( int outerUnroll ){
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

int BrookGbsa::getParticleSizeCeiling( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookGbsa::getParticleSizeCeiling";

// ---------------------------------------------------------------------------------------

   if( _particleSizeCeiling < 0 ){
      BrookGbsa* localThis = const_cast<BrookGbsa* const>(this);
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

int BrookGbsa::getDuplicationFactor( void ) const {
   return _duplicationFactor;
}

/** 
 * Get partial force stream width
 * 
 * @return  partial force stream width
 *
 */

int BrookGbsa::getPartialForceStreamWidth( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _partialForceStreamWidth;
}

/** 
 * Get partial force stream height
 * 
 * @return  partial force stream height 
 *
 */

int BrookGbsa::getPartialForceStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamHeight";

// ---------------------------------------------------------------------------------------

   return _partialForceStreamHeight;
}

/** 
 * Get partial force stream size
 * 
 * @return  partial force stream size
 *
 */

int BrookGbsa::getPartialForceStreamSize( void ) const {
   return _partialForceStreamSize;
}

/** 
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookGbsa::getGbsaParticleStreamSize( void ) const {
   return _gbsaParticleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookGbsa::getGbsaParticleStreamWidth( void ) const {
   return _gbsaParticleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookGbsa::getGbsaParticleStreamHeight( void ) const {
   return _gbsaParticleStreamHeight;
}

/** 
 * Get Obc particle radii stream 
 *
 * @return  Obc particle radii stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcParticleRadii( void ) const {
   return _gbsaStreams[ObcParticleRadiiStream];
}

/** 
 * Get Obc scaled particle radii stream 
 *
 * @return  Obc scaled particle radii stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcScaledParticleRadii( void ) const {
   return _gbsaStreams[ObcScaledParticleRadiiStream];
}

/** 
 * Get Obc particle radii w/ dielectric offset
 *
 * @return  Obc particle radii w/ dielectric offset
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcParticleRadiiWithDielectricOffset( void ) const {
   return _gbsaStreams[ObcParticleRadiiWithDielectricOffsetStream];
}

/** 
 * Get Obc Born radii
 *
 * @return  Obc Born radii
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcBornRadii( void ) const {
   return _gbsaStreams[ObcBornRadiiStream];
}

/** 
 * Get Obc Born radii2
 *
 * @return  Obc Born radii2
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcBornRadii2( void ) const {
   return _gbsaStreams[ObcBornRadii2Stream];
}

/** 
 * Get Obc intermediate force stream
 *
 * @return  Obcintermediate force stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcIntermediateForce( void ) const {
   return _gbsaStreams[ObcIntermediateForceStream];
}

/** 
 * Get Obc chain stream 
 *
 * @return  Obc chain stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcChain( void ) const {
   return _gbsaStreams[ObcChainStream];
}

/** 
 * Get force streams 
 *
 * @return  force streams
 *
 */

BrookFloatStreamInternal** BrookGbsa::getForceStreams( void ){
   return _gbsaForceStreams;
}

/** 
 * Return true if force[index] stream is set
 *
 * @param    index into force stream
 * @return   true if index is valid && force[index] stream is set; else false
 *
 */

int BrookGbsa::isForceStreamSet( int index ) const {
   return (index >= 0 && index < getNumberOfForceStreams() && _gbsaForceStreams[index]) ? 1 : 0;
}

/** 
 * Return true if Born radii have been initialized
 *
 * @return   true if Born radii have been initialized
 *
 */

int BrookGbsa::haveBornRadiiBeenInitialized( void ) const {
   return  _bornRadiiInitialized;
}

/** 
 * Calculate Born radii
 *
 * @return  calculate Born radii
 *
 */

int BrookGbsa::calculateBornRadii( const Stream& positions ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName                 = "BrookGbsa::calculateBornRadii";
   static const int PrintOn                            = 0;

// ---------------------------------------------------------------------------------------

   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (positions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   const RealOpenMM* coordinates                       = (RealOpenMM*) positionStream.getData( );

   // load coordinates into RealOpenMM 2d array

   int numberOfParticles                               = getNumberOfParticles();

   RealOpenMM** particleCoordinates                    = new RealOpenMM*[numberOfParticles];
   RealOpenMM* particleCoordinatesBlk                  = new RealOpenMM[3*numberOfParticles];

   // Born radii array size needs to match stream size since it will
   // be written down to board

   int streamSize                                      = getGbsaParticleStreamSize();
   RealOpenMM* bornRadii                               = new RealOpenMM[streamSize];
   memset( bornRadii, 0, sizeof( RealOpenMM )*streamSize );

   RealOpenMM* obcChain                                = new RealOpenMM[streamSize];
   memset( obcChain, 0, sizeof( RealOpenMM )*streamSize );

   int index                                           = 0;
   RealOpenMM* particleCoordinatesBlkPtr               = particleCoordinatesBlk;
   for( int ii = 0; ii < numberOfParticles; ii++ ){

      particleCoordinates[ii]    = particleCoordinatesBlkPtr;
      particleCoordinatesBlkPtr += 3;

      particleCoordinates[ii][0] = coordinates[index++];
      particleCoordinates[ii][1] = coordinates[index++];
      particleCoordinates[ii][2] = coordinates[index++];

   }

   // calculate Born radii

   _cpuObc->computeBornRadii( particleCoordinates, bornRadii, obcChain );

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\n%s: atms=%d\n", methodName.c_str(), numberOfParticles );
      for( int ii = 0; ii < numberOfParticles; ii++ ){
         (void) fprintf( getLog(), "%d coord=[%.5e %.5e %.5e]  bR=%.5e obcChain=%.6e\n", ii,
                         particleCoordinates[ii][0], particleCoordinates[ii][1], particleCoordinates[ii][2], bornRadii[ii], obcChain[ii] );
      }
   }
  
   // write radii to board and set flag to indicate radii calculated once

   _gbsaStreams[ObcBornRadiiStream]->loadFromArray( bornRadii );
   _gbsaStreams[ObcChainStream]->loadFromArray( obcChain );
   _bornRadiiInitialized   = 1;

   // free memory

   delete[] particleCoordinatesBlk;
   delete[] particleCoordinates;
   delete[] bornRadii;
   delete[] obcChain;

   return DefaultReturnValue;
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

int BrookGbsa::initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _gbsaParticleStreamSize     = getParticleStreamSize( platform );
   _gbsaParticleStreamWidth    = getParticleStreamWidth( platform );
   _gbsaParticleStreamHeight   = getParticleStreamHeight( platform );

   int innerUnroll         = getInnerLoopUnroll();
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

   _partialForceStreamSize    = _gbsaParticleStreamSize*getDuplicationFactor()/innerUnroll;
   _partialForceStreamHeight  = _partialForceStreamSize/_partialForceStreamWidth;
   _partialForceStreamHeight += ( (_partialForceStreamSize % _partialForceStreamWidth) ? 1 : 0);
   _partialForceStreamSize    = _partialForceStreamHeight*_partialForceStreamWidth;

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

int BrookGbsa::initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::initializeStreams";
   static const double dangleValue          = 0.0;

// ---------------------------------------------------------------------------------------

   int gbsaParticleStreamSize   = getGbsaParticleStreamSize();
   int gbsaParticleStreamWidth  = getGbsaParticleStreamWidth();

   // particle radii & charge

    _gbsaStreams[ObcParticleRadiiStream]                        = new BrookFloatStreamInternal( BrookCommon::ObcParticleRadiiStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float2, dangleValue );

   // scaled particle radii

    _gbsaStreams[ObcScaledParticleRadiiStream]                  = new BrookFloatStreamInternal( BrookCommon::ObcScaledParticleRadiiStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float2, dangleValue );

   // particle radii w/ DielectricOffset

    _gbsaStreams[ObcParticleRadiiWithDielectricOffsetStream]    = new BrookFloatStreamInternal( BrookCommon::ObcParticleRadiiWithDielectricOffsetStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // Born radii

    _gbsaStreams[ObcBornRadiiStream]                            = new BrookFloatStreamInternal( BrookCommon::ObcBornRadiiStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // Born2 radii

    _gbsaStreams[ObcBornRadii2Stream]                           = new BrookFloatStreamInternal( BrookCommon::ObcBornRadii2Stream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // IntermediateForce

    _gbsaStreams[ObcIntermediateForceStream]                    = new BrookFloatStreamInternal( BrookCommon::ObcIntermediateForceStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float4, dangleValue );

   // Obc chain

    _gbsaStreams[ObcChainStream]                                = new BrookFloatStreamInternal( BrookCommon::ObcChainStream,
                                                                                                gbsaParticleStreamSize, gbsaParticleStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // partial force streams

   std::string partialForceStream = BrookCommon::PartialForceStream;
   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      std::stringstream name;
      name << partialForceStream << ii;
      _gbsaForceStreams[ii] = new BrookFloatStreamInternal( name.str(), getPartialForceStreamSize(),
                                                            getPartialForceStreamWidth(), BrookStreamInternal::Float4, dangleValue );
   }

   return DefaultReturnValue;
}

/*  
 * Setup of Gbsa parameters
 *
 * @param particleParameters        vector of OBC parameters [particleI][0=charge]
 *                                                           [particleI][1=radius]
 *                                                           [particleI][2=scaling factor]
 * @param solventDielectric     solvent dielectric
 * @param soluteDielectric      solute dielectric
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookGbsa::setup( const std::vector<std::vector<double> >& vectorOfParticleParameters, 
                      double solventDielectric, double soluteDielectric, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const int particleParametersSize  = 3; 
   static const int maxErrors               = 20; 
   static const std::string methodName      = "BrookGbsa::setup";

// ---------------------------------------------------------------------------------------

   int numberOfParticles  = (int) vectorOfParticleParameters.size();
   setNumberOfParticles( numberOfParticles );

   _solventDielectric = solventDielectric;
   _soluteDielectric  = soluteDielectric;

   // initialize stream sizes and then Brook streams

   initializeStreamSizes( numberOfParticles, platform );
   initializeStreams( platform );

   int particleStreamSize                   = getGbsaParticleStreamSize();
   BrookOpenMMFloat* radiiAndCharge         = new BrookOpenMMFloat[particleStreamSize*2];
   BrookOpenMMFloat* scaledRadiiAndOffset   = new BrookOpenMMFloat[particleStreamSize*2];
   memset( radiiAndCharge, 0, particleStreamSize*2*sizeof( BrookOpenMMFloat ) );
   memset( scaledRadiiAndOffset, 0, particleStreamSize*2*sizeof( BrookOpenMMFloat ) );

   _charges                                 = new RealOpenMM[particleStreamSize];

   // used by CpuObc to calculate initial Born radii

   vector<RealOpenMM> particleRadii(numberOfParticles);
   vector<RealOpenMM> scaleFactors(numberOfParticles);

   float dielectricOffset                   = getDielectricOffset();

   // loop over particle parameters
   // track any errors and then throw exception
   //    check parameter vector is right size
   // set parameter entries or board and arrays used by CpuObc

   int vectorIndex  = 0;
   int errors       = 0;
   std::stringstream message;

   typedef std::vector< std::vector<double> > VectorOfDoubleVectors;
   typedef VectorOfDoubleVectors::const_iterator VectorOfDoubleVectorsCI;

   for( VectorOfDoubleVectorsCI ii = vectorOfParticleParameters.begin(); ii != vectorOfParticleParameters.end(); ii++ ){

      std::vector<double> particleParameters = *ii;

      if( particleParameters.size() != particleParametersSize && errors < maxErrors ){
         message << methodName << " parameter size=" << particleParameters.size() << " for parameter vector index=" << vectorIndex << " is less than expected.\n";
         errors++;
      } else {

         double charge                            = particleParameters[0];     
         double radius                            = particleParameters[1];     
         double scalingFactor                     = particleParameters[2];     

         int streamIndex                          = 2*vectorIndex;

         particleRadii[vectorIndex]               = static_cast<RealOpenMM> (radius);
         scaleFactors[vectorIndex]                = static_cast<RealOpenMM> (scalingFactor);
         _charges[vectorIndex]                    = static_cast<RealOpenMM> (charge);

         radiiAndCharge[streamIndex]              = static_cast<BrookOpenMMFloat> (radius);
         radiiAndCharge[streamIndex+1]            = static_cast<BrookOpenMMFloat> (charge);

         scaledRadiiAndOffset[streamIndex+1]      = static_cast<BrookOpenMMFloat> (radius - dielectricOffset);
         scaledRadiiAndOffset[streamIndex]        = static_cast<BrookOpenMMFloat> (scaledRadiiAndOffset[streamIndex+1]*scalingFactor);

//         scaledRadiiAndOffset[streamIndex]        = static_cast<BrookOpenMMFloat> (radius - dielectricOffset);
//         scaledRadiiAndOffset[streamIndex+1]      = static_cast<BrookOpenMMFloat> (scaledRadiiAndOffset[streamIndex]*scalingFactor);

      }

      vectorIndex++;
   }

   // throw exception if errors detected

   if( errors ){
      throw OpenMMException( message.str() );
   }

   // load streams

   _gbsaStreams[ObcParticleRadiiStream]->loadFromArray( radiiAndCharge );
   _gbsaStreams[ObcScaledParticleRadiiStream]->loadFromArray( scaledRadiiAndOffset );

   delete[] radiiAndCharge;
   delete[] scaledRadiiAndOffset;

   // setup for Born radii calculation

   ObcParameters* obcParameters  = new ObcParameters( numberOfParticles, ObcParameters::ObcTypeII );
   obcParameters->setAtomicRadii( particleRadii);

   obcParameters->setScaledRadiusFactors( scaleFactors );
   obcParameters->setSolventDielectric( static_cast<RealOpenMM>(solventDielectric) );
   obcParameters->setSoluteDielectric(  static_cast<RealOpenMM>(soluteDielectric)  );

   _cpuObc  = new CpuObc(obcParameters);
   _cpuObc->setIncludeAceApproximation( true );

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

int BrookGbsa::initializePartialForceStreamSize( int particleStreamSize, int particleStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::initializePartialForceStreamSize";
   //static const int debug                   = 1;

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
 * Setup of j-stream dimensions
 * 
 * Get contents of object
 *
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

std::string BrookGbsa::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamWidth() );
   message << _getLine( tab, "Partial force stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamHeight() );
   message << _getLine( tab, "Partial force stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getPartialForceStreamSize() );
   message << _getLine( tab, "Partial force stream size:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _gbsaStreams[ii] ){
         message << _gbsaStreams[ii]->getContentsString( );
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
      if( _gbsaForceStreams[ii] ){
         message << _gbsaForceStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}

/** 
 * Compute forces
 * 
 */

void BrookGbsa::computeForces( BrookStreamImpl& positionStream, BrookStreamImpl& forceStream ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName   = "BrookGbsa::executeForces";
   static const int PrintOn              = 1; 
   float mergeNonObcForces               = 1.0f;
   float kcalMolTokJNM                   = -0.4184f;

// ---------------------------------------------------------------------------------------

   float includeAceTerm                                = (float) (includeAce());
   BrookFloatStreamInternal**  gbsaForceStreams        = getForceStreams();

   // calculate Born radii 

(void) fprintf( getLog(), "\nPost kCalculateBornRadii: obcParticleRadiiWithDielectricOffset & obcScaledParticleRadii not set correctly!!!!!\n" );

   kCalculateBornRadii(   (float) getNumberOfParticles(),
                          (float) getParticleSizeCeiling(),
                          (float) getDuplicationFactor(),
                          (float) getParticleStreamWidth( ),
                          (float) getPartialForceStreamWidth( ),
                          positionStream.getBrookStream(),
                          getObcScaledParticleRadii()->getBrookStream(),
                          gbsaForceStreams[0]->getBrookStream() );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( 0 && PrintOn && getLog() ){

      (void) fprintf( getLog(), "\n%s Post kCalculateBornRadii: atms=%d ceil=%d dup=%d particleStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      methodName.c_str(), getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ) );

      BrookStreamInternal* brookStreamInternalF  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalF->printToFile( getLog() );

      (void) fprintf( getLog(), "\nRadii\n" );
      getObcParticleRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nObcScaledParticleRadii\n" );
      getObcScaledParticleRadii()->printToFile( getLog() );

   }

   // ---------------------------------------------------------------------------------------

   kPostCalculateBornRadii_nobranch(
                          (float) getDuplicationFactor(),
                          (float) getParticleStreamWidth( ),
                          (float) getPartialForceStreamWidth( ),
                          (float) getNumberOfParticles(),
                          (float) getParticleSizeCeiling(),
                          (float) getInnerLoopUnroll(),
                          kcalMolTokJNM,
                          (float) mergeNonObcForces,
                          gbsaForceStreams[0]->getBrookStream(),
                          getObcParticleRadii()->getBrookStream(),
                          getObcBornRadii()->getBrookStream(),
                          getObcChain()->getBrookStream() );   

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( 0 && PrintOn && getLog() ){

      (void) fprintf( getLog(), "\n%s Post kPostCalculateBornRadii_nobranch: atms=%d ceil=%d dup=%d particleStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      methodName.c_str(), getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ) );

      BrookStreamInternal* brookStreamInternalF  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalF->printToFile( getLog() );

      (void) fprintf( getLog(), "\nInput\n" );
      gbsaForceStreams[0]->printToFile( getLog() );

      (void) fprintf( getLog(), "\nObcParticleRadii\n" );
      getObcParticleRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      getObcBornRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nObcChain\n" );
      getObcChain()->printToFile( getLog() );

   }

// ---------------------------------------------------------------------------------------

   // seecond major loop

   kObcLoop1( (float) getNumberOfParticles(),
              (float) getParticleSizeCeiling(),
              (float) getDuplicationFactor(),
              (float) getParticleStreamWidth( ),
              (float) getPartialForceStreamWidth( ),
              getSoluteDielectric(),
              getSolventDielectric(),
              includeAceTerm,

              positionStream.getBrookStream(),

              getObcBornRadii()->getBrookStream(),
              getObcParticleRadii()->getBrookStream(),

              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream()
            );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( 1 && PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kObcLoop1: atms=%d ceil=%d dup=%d particleStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ),
                      getSoluteDielectric(),
                      getSolventDielectric(), includeAceTerm );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      getObcBornRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nParticleR\n" );
      getObcParticleRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStreams output\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

   }

// ---------------------------------------------------------------------------------------

   // gather for first loop

   kPostObcLoop1_nobranch(
              (float) getDuplicationFactor(),
              (float) getParticleStreamWidth( ),
              (float) getPartialForceStreamWidth( ),
              (float) getNumberOfParticles(),
              (float) getParticleSizeCeiling(),
              (float) getInnerLoopUnroll(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream(),
              getObcChain()->getBrookStream(),
              getObcBornRadii()->getBrookStream(),
              getObcIntermediateForce()->getBrookStream(),
              getObcBornRadii2()->getBrookStream() );
 
// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog()){

      (void) fprintf( getLog(), "\nPost kPostObcLoop1_nobranch: dup=%d aStrW=%d pStrW=%d no.atms=%3d ceil=%3d Unroll=%1d\n",
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ),
                      getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getInnerLoopUnroll() );

      (void) fprintf( getLog(), "\nForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

      (void) fprintf( getLog(), "\nObcChain\n" );
      getObcChain()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      getObcBornRadii()->printToFile( getLog() );

      // output

      (void) fprintf( getLog(), "\nObcIntermediateForce output\n" );
      getObcIntermediateForce()->printToFile( getLog() );

      // output

      (void) fprintf( getLog(), "\nObcBornRadii2 output\n" );
      getObcBornRadii2()->printToFile( getLog() );

   }

// ---------------------------------------------------------------------------------------

   // second major loop

(void) fprintf( getLog(), "\nkObcLoop2 messed up see /home/friedrim/src/openmmWork/trunk/OpenMM/platforms/brook/src/gpu/kObcBaseD2.br\n" );

   kObcLoop2( (float) getNumberOfParticles(),
              (float) getParticleSizeCeiling(),
              (float) getDuplicationFactor(),
              (float) getParticleStreamWidth( ),
              (float) getPartialForceStreamWidth( ),
              positionStream.getBrookStream(),
              getObcScaledParticleRadii()->getBrookStream(),
              getObcBornRadii2()->getBrookStream(),
              getObcBornRadii2()->getBrookStream(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream()
            );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kObcLoop2: no.atms=%5d ceil=%3d dup=%3d strW=%3d pStrW=%3d\n",
                      getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ) );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nObcScaledParticleRadii\n" );
      getObcScaledParticleRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\ngetObcBornRadii2\n" );
      getObcBornRadii2()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

   }

// ---------------------------------------------------------------------------------------

   // gather for second loop

    kPostObcLoop2_nobranch(
              (float) getDuplicationFactor(),
              (float) getParticleStreamWidth( ),
              (float) getPartialForceStreamWidth( ),
              (float) getNumberOfParticles(),
              (float) getParticleSizeCeiling(),
              (float) getInnerLoopUnroll(),
              kcalMolTokJNM,
              mergeNonObcForces,
              getObcIntermediateForce()->getBrookStream(),
              forceStream.getBrookStream(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream(),
              getObcParticleRadii()->getBrookStream(),
              getObcBornRadii()->getBrookStream(),
              getObcChain()->getBrookStream(),
              forceStream.getBrookStream()
           );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kPostObcLoop2_nobranch: atms=%d ceil=%d dup=%d particleStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      getNumberOfParticles(),
                      getParticleSizeCeiling(),
                      getDuplicationFactor(),
                      getParticleStreamWidth( ),
                      getPartialForceStreamWidth( ),
                      getSoluteDielectric(),
                      getSolventDielectric(), includeAceTerm );

      (void) fprintf( getLog(), "\nPartialForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

      BrookStreamInternal* brookStreamInternalF  = forceStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nForceStream\n" );
      brookStreamInternalF->printToFile( getLog() );

      (void) fprintf( getLog(), "\nChain\n" );
      getObcChain()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      getObcBornRadii()->printToFile( getLog() );

   }

   // ---------------------------------------------------------------------------------------
}
