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

   _atomSizeCeiling           = -1;
   _outerUnroll               =  4;
   _innerUnroll               =  4;

   _partialForceStreamWidth   = 64;
   _partialForceStreamHeight  = -1;
   _partialForceStreamSize    = -1;

   _gbsaAtomStreamWidth       = -1;
   _gbsaAtomStreamHeight      = -1;
   _gbsaAtomStreamSize        = -1;

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
      _atomSizeCeiling = -1;
   }
   _outerUnroll = _outerUnroll;
   return _outerUnroll;
}

/** 
 * Get atom ceiling parameter
 * 
 * @return atom ceiling parameter
 *
 */

int BrookGbsa::getAtomSizeCeiling( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookGbsa::getAtomSizeCeiling";

// ---------------------------------------------------------------------------------------

   if( _atomSizeCeiling < 0 ){
      BrookGbsa* localThis = const_cast<BrookGbsa* const>(this);
      localThis->_atomSizeCeiling = localThis->getNumberOfAtoms() % localThis->getOuterLoopUnroll();
      if( localThis->_atomSizeCeiling ){
         localThis->_atomSizeCeiling = localThis->getOuterLoopUnroll() - localThis->_atomSizeCeiling;
      }   
      localThis->_atomSizeCeiling += localThis->getNumberOfAtoms();
   }

   return _atomSizeCeiling;
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
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookGbsa::getGbsaAtomStreamSize( void ) const {
   return _gbsaAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookGbsa::getGbsaAtomStreamWidth( void ) const {
   return _gbsaAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookGbsa::getGbsaAtomStreamHeight( void ) const {
   return _gbsaAtomStreamHeight;
}

/** 
 * Get Obc atomic radii stream 
 *
 * @return  Obc atomic radii stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcAtomicRadii( void ) const {
   return _gbsaStreams[ObcAtomicRadiiStream];
}

/** 
 * Get Obc scaled atomic radii stream 
 *
 * @return  Obc scaled atomic radii stream
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcScaledAtomicRadii( void ) const {
   return _gbsaStreams[ObcScaledAtomicRadiiStream];
}

/** 
 * Get Obc atomic radii w/ dielectric offset
 *
 * @return  Obc atomic radii w/ dielectric offset
 *
 */

BrookFloatStreamInternal* BrookGbsa::getObcAtomicRadiiWithDielectricOffset( void ) const {
   return _gbsaStreams[ObcAtomicRadiiWithDielectricOffsetStream];
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

   int numberOfAtoms                                   = getNumberOfAtoms();

   RealOpenMM** atomCoordinates                        = new RealOpenMM*[numberOfAtoms];
   RealOpenMM* atomCoordinatesBlk                      = new RealOpenMM[3*numberOfAtoms];

   // Born radii array size needs to match stream size since it will
   // be written down to board

   int streamSize                                      = getGbsaAtomStreamSize();
   RealOpenMM* bornRadii                               = new RealOpenMM[streamSize];
   memset( bornRadii, 0, sizeof( RealOpenMM )*streamSize );

   RealOpenMM* obcChain                                = new RealOpenMM[streamSize];
   memset( obcChain, 0, sizeof( RealOpenMM )*streamSize );

   int index                                           = 0;
   RealOpenMM* atomCoordinatesBlkPtr                   = atomCoordinatesBlk;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      atomCoordinates[ii]    = atomCoordinatesBlkPtr;
      atomCoordinatesBlkPtr += 3;

      atomCoordinates[ii][0] = coordinates[index++];
      atomCoordinates[ii][1] = coordinates[index++];
      atomCoordinates[ii][2] = coordinates[index++];

   }

   // calculate Born radii

   _cpuObc->computeBornRadii( atomCoordinates, bornRadii, obcChain );

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\n%s: atms=%d\n", methodName.c_str(), numberOfAtoms );
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         (void) fprintf( getLog(), "%d coord=[%.5e %.5e %.5e]  bR=%.5e obcChain=%.6e\n", ii,
                         atomCoordinates[ii][0], atomCoordinates[ii][1], atomCoordinates[ii][2], bornRadii[ii], obcChain[ii] );
      }
   }
  
   // write radii to board and set flag to indicate radii calculated once

   _gbsaStreams[ObcBornRadiiStream]->loadFromArray( bornRadii );
   _gbsaStreams[ObcChainStream]->loadFromArray( obcChain );
   _bornRadiiInitialized   = 1;

   // free memory

   delete[] atomCoordinatesBlk;
   delete[] atomCoordinates;
   delete[] bornRadii;
   delete[] obcChain;

   return DefaultReturnValue;
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfAtoms             number of atoms
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookGbsa::initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _gbsaAtomStreamSize     = getAtomStreamSize( platform );
   _gbsaAtomStreamWidth    = getAtomStreamWidth( platform );
   _gbsaAtomStreamHeight   = getAtomStreamHeight( platform );

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

   _partialForceStreamSize    = _gbsaAtomStreamSize*getDuplicationFactor()/innerUnroll;
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

   int gbsaAtomStreamSize   = getGbsaAtomStreamSize();
   int gbsaAtomStreamWidth  = getGbsaAtomStreamWidth();

   // atomic radii & charge

    _gbsaStreams[ObcAtomicRadiiStream]                          = new BrookFloatStreamInternal( BrookCommon::ObcAtomicRadiiStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float2, dangleValue );

   // scaled atomic radii

    _gbsaStreams[ObcScaledAtomicRadiiStream]                    = new BrookFloatStreamInternal( BrookCommon::ObcScaledAtomicRadiiStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float2, dangleValue );

   // atomic radii w/ DielectricOffset

    _gbsaStreams[ObcAtomicRadiiWithDielectricOffsetStream]      = new BrookFloatStreamInternal( BrookCommon::ObcAtomicRadiiWithDielectricOffsetStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // Born radii

    _gbsaStreams[ObcBornRadiiStream]                            = new BrookFloatStreamInternal( BrookCommon::ObcBornRadiiStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // Born2 radii

    _gbsaStreams[ObcBornRadii2Stream]                           = new BrookFloatStreamInternal( BrookCommon::ObcBornRadii2Stream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float, dangleValue );

   // IntermediateForce

    _gbsaStreams[ObcIntermediateForceStream]                    = new BrookFloatStreamInternal( BrookCommon::ObcIntermediateForceStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
                                                                                                BrookStreamInternal::Float4, dangleValue );

   // Obc chain

    _gbsaStreams[ObcChainStream]                                = new BrookFloatStreamInternal( BrookCommon::ObcChainStream,
                                                                                                gbsaAtomStreamSize, gbsaAtomStreamWidth,
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
 * @param atomParameters        vector of OBC parameters [atomI][0=charge]
 *                                                       [atomI][1=radius]
 *                                                       [atomI][2=scaling factor]
 * @param solventDielectric     solvent dielectric
 * @param soluteDielectric      solute dielectric
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookGbsa::setup( const std::vector<std::vector<double> >& vectorOfAtomParameters, 
                      double solventDielectric, double soluteDielectric, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const int atomParametersSize      = 3; 
   static const int maxErrors               = 20; 
   static const std::string methodName      = "BrookGbsa::setup";

// ---------------------------------------------------------------------------------------

   int numberOfAtoms  = (int) vectorOfAtomParameters.size();
   setNumberOfAtoms( numberOfAtoms );

   _solventDielectric = solventDielectric;
   _soluteDielectric  = soluteDielectric;

   // initialize stream sizes and then Brook streams

   initializeStreamSizes( numberOfAtoms, platform );
   initializeStreams( platform );

   int atomStreamSize                       = getGbsaAtomStreamSize();
   BrookOpenMMFloat* radiiAndCharge         = new BrookOpenMMFloat[atomStreamSize*2];
   BrookOpenMMFloat* scaledRadiiAndOffset   = new BrookOpenMMFloat[atomStreamSize*2];
   memset( radiiAndCharge, 0, atomStreamSize*2*sizeof( BrookOpenMMFloat ) );
   memset( scaledRadiiAndOffset, 0, atomStreamSize*2*sizeof( BrookOpenMMFloat ) );

   _charges                                 = new RealOpenMM[atomStreamSize];

   // used by CpuObc to calculate initial Born radii

   vector<RealOpenMM> atomicRadii(numberOfAtoms);
   vector<RealOpenMM> scaleFactors(numberOfAtoms);

   float dielectricOffset                  = getDielectricOffset();

   // loop over atom parameters
   // track any errors and then throw exception
   //    check parameter vector is right size
   // set parameter entries or board and arrays used by CpuObc

   int vectorIndex  = 0;
   int errors       = 0;
   std::stringstream message;

   typedef std::vector< std::vector<double> > VectorOfDoubleVectors;
   typedef VectorOfDoubleVectors::const_iterator VectorOfDoubleVectorsCI;

   for( VectorOfDoubleVectorsCI ii = vectorOfAtomParameters.begin(); ii != vectorOfAtomParameters.end(); ii++ ){

      std::vector<double> atomParameters = *ii;

      if( atomParameters.size() != atomParametersSize && errors < maxErrors ){
         message << methodName << " parameter size=" << atomParameters.size() << " for parameter vector index=" << vectorIndex << " is less than expected.\n";
         errors++;
      } else {

         double charge                            = atomParameters[0];     
         double radius                            = atomParameters[1];     
         double scalingFactor                     = atomParameters[2];     

         int streamIndex                          = 2*vectorIndex;

         atomicRadii[vectorIndex]                 = static_cast<RealOpenMM> (radius);
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

   _gbsaStreams[ObcAtomicRadiiStream]->loadFromArray( radiiAndCharge );
   _gbsaStreams[ObcScaledAtomicRadiiStream]->loadFromArray( scaledRadiiAndOffset );

   delete[] radiiAndCharge;
   delete[] scaledRadiiAndOffset;

   // setup for Born radii calculation

   ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   obcParameters->setAtomicRadii( atomicRadii);

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
 * @param atomStreamSize        atom stream size
 * @param atomStreamWidth       atom stream width
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * */

int BrookGbsa::initializePartialForceStreamSize( int atomStreamSize, int atomStreamWidth ){

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

   _partialForceStreamSize    = atomStreamSize*getDuplicationFactor()/innerUnroll;
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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfForceStreams() );
   message << _getLine( tab, "Number of force streams:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getDuplicationFactor() );
   message << _getLine( tab, "Duplication factor:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getInnerLoopUnroll () )
   message << _getLine( tab, "Inner loop unroll:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getOuterLoopUnroll() )
   message << _getLine( tab, "Outer loop unroll:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomSizeCeiling() );
   message << _getLine( tab, "Atom ceiling:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 

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

/*  
 * Calculate OBC energy
 *
 * @param atomPositions        atom positions

 * @return energy
 *
 * @throw OpenMMException if _cpuObc or charges are not set
 *
 * */
    
double BrookGbsa::getEnergy( const Stream& atomPositions ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::getEnergy";

// ---------------------------------------------------------------------------------------

   // validate initialization

   if( _cpuObc == NULL ){
      std::stringstream message;
      message << methodName << " _cpuObc not set.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }   

   if( _charges == NULL ){
      std::stringstream message;
      message << methodName << " charges not set.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }   

   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (atomPositions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   BrookOpenMMFloat* positionsF                        = (BrookOpenMMFloat*) positionStream.getData();

   RealOpenMM** positions                              = copy1DArrayTo2DArray( positionStream.getSize(), 3, positionsF );
   RealOpenMM** forces                                 = allocateRealArray( positionStream.getSize(), 3 ); 

   _cpuObc->computeImplicitSolventForces( positions, _charges, forces, 1 ); 

   disposeRealArray( forces );
   disposeRealArray( positions );

   return _cpuObc->getEnergy();
}
