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

   _gbsaStreamWidth           = -1;
   _gbsaStreamHeight          = -1;
   _gbsaStreamSize            = -1;

   _duplicationFactor         =  4;

   _solventDielectric         = 78.3;
   _soluteDielectric          = 1.0;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _gbsaStreams[ii] = NULL;
   }

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      _gbsaForceStreams[ii] = NULL;
   }

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

double BrookGbsa::getSoluteDielectric( void ) const {
   return _soluteDielectric;
}

/** 
 * Get solvent dielectric
 * 
 * @return   solvent dielectric
 *
 */

double BrookGbsa::getSolventDielectric( void ) const {
   return _solventDielectric;
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
 * @param platform  platform
 *
 * @return  partial force stream width
 *
 */

int BrookGbsa::getPartialForceStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamWidth";

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
 * @param platform platform
 *
 * @return  partial force stream height 
 *
 */

int BrookGbsa::getPartialForceStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamHeight";

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

int BrookGbsa::getPartialForceStreamHeight( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamHeight";

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

int BrookGbsa::getPartialForceStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookGbsa::getPartialForceStreamSize";

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

int BrookGbsa::getPartialForceStreamSize( void ) const {
   return _partialForceStreamSize;
}

/** 
 * Get Gbsa stream size
 *
 * @return  Gbsa stream size
 *
 */

int BrookGbsa::getGbsaStreamSize( void ) const {
   return _gbsaStreamSize;
}

/** 
 * Get gbsa stream width
 *
 * @return  gbsa stream width
 *
 */

int BrookGbsa::getGbsaStreamWidth( void ) const {
   return _gbsaStreamWidth;
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
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int atomStreamSize                  = getAtomStreamSize( platform );
   int atomStreamWidth                 = getAtomStreamWidth( platform );

   initializeExclusionStreamSize( atomStreamSize, atomStreamWidth );
   initializeJStreamSize( atomStreamSize, atomStreamWidth );
   initializeOuterVdwStreamSize( atomStreamSize, atomStreamWidth );
   initializePartialForceStreamSize( atomStreamSize, atomStreamWidth );

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
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());

   // exclusion

    _gbsaStreams[ExclusionStream]          = new BrookFloatStreamInternal( BrookStreamFactory::NonBondedExclusionStream,
                                                                                getExclusionStreamSize(), getExclusionStreamWidth(),
                                                                                BrookStreamInternal::Float, dangleValue );

   // outer vdw

   _gbsaStreams[OuterVdwStream]            = new BrookFloatStreamInternal( BrookStreamFactory::OuterVdwStream, getAtomStreamSize(), 
                                                                                getAtomStreamWidth(), BrookStreamInternal::Float2, dangleValue );

   // inner sigma & epsilon

   _gbsaStreams[InnerSigmaStream]         = new BrookFloatStreamInternal( BrookStreamFactory::InnerSigmaStream, getJStreamSize(), 
                                                                               getJStreamWidth(), BrookStreamInternal::Float4, dangleValue );

   _gbsaStreams[InnerEpsilonStream]       = new BrookFloatStreamInternal( BrookStreamFactory::InnerEpsilonStream, getJStreamSize(),
                                                                               getJStreamWidth(), BrookStreamInternal::Float4, dangleValue );

   // charge stream

   _gbsaStreams[ChargeStream]              = new BrookFloatStreamInternal( BrookStreamFactory::NonBondedChargeStream, getAtomStreamSize(),
                                                                                getAtomStreamWidth(), BrookStreamInternal::Float, dangleValue );

   
   // partial force stream
   std::string partialForceStream = BrookStreamFactory::PartialForceStream;
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
 * @param atomParameters        vector of OBC parameters [atomI][0=index]
 *                                                       [atomI][1=charge]
 *                                                       [atomI][2=radius]
 *                                                       [atomI][2=scaling factor]
 * @param solventDielectric     solvent dielectric
 * @param soluteDielectric      solute dielectric
 *
 * @return nonzero value if error
 *
 * */
    
int BrookGbsa::setup( const std::vector<std::vector<double> >& atomParameters, 
                      double solventDielectric, double soluteDielectric );
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookGbsa::setup";

// ---------------------------------------------------------------------------------------

   setNumberOfAtoms( (int) atomParameters.size() );

   solventDielectric = _solventDielectric;
   soluteDielectric  = _soluteDielectric;

   initializeStreamSizes( getNumberOfAtoms(), platform );
   initializeStreams( platform );

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

   (void) LOCAL_SPRINTF( value, "%d", getExclusionStreamWidth() );
   message << _getLine( tab, "Exclusion stream width:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 
/*
   message << _getLine( tab, "ExclusionStream:",      (getExclusionStream()    ? Set : NotSet) ); 
   message << _getLine( tab, "VdwStream:",            (getOuterVdwStream()     ? Set : NotSet) ); 
   message << _getLine( tab, "ChargeStream:",         (getChargeStream()       ? Set : NotSet) ); 
   message << _getLine( tab, "SigmaStream:",          (getInnerSigmaStream()   ? Set : NotSet) ); 
   message << _getLine( tab, "EpsilonStream:",        (getInnerEpsilonStream() ? Set : NotSet) ); 
*/
 
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
