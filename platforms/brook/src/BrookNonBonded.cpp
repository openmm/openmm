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
#include "BrookNonBonded.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include "OpenMMException.h"
#include "gpu/invmap.h"
#include "gpu/kforce.h"

using namespace OpenMM;
using namespace std;

/** 
 * Constructor
 * 
 */

BrookNonBonded::BrookNonBonded(  ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::BrookNonBonded";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _atomSizeCeiling           = -1;
   _outerUnroll               = 4;
   _innerUnroll               = 4;

   _partialForceStreamWidth   = 64;
   _partialForceStreamHeight  = -1;
   _partialForceStreamSize    = -1;

   _duplicationFactor         = 4;

   _exclusionStreamWidth      = -1;
   _exclusionStreamHeight     = -1;
   _exclusionStreamSize       = -1;

   _jStreamWidth              = -1;
   _jStreamHeight             = -1;
   _jStreamSize               = -1;

   _exclusionStream           = NULL;
   _sigmaStream               = NULL;
   _epsilonStream             = NULL;

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

   delete _exclusionStream;
   delete _sigmaStream;
   delete _epsilonStream;

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

int BrookNonBonded::getAtomSizeCeiling( void ) const {

   if( _atomSizeCeiling < 0 ){
      BrookNonBonded* localThis = const_cast<BrookNonBonded* const>(this);
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

int BrookNonBonded::getDuplicationFactor( void ) const {
   return _duplicationFactor;
}

/** 
 * Get j-stream width
 * 
 * @param platform  platform
 *
 * @return  j-stream width
 *
 */

int BrookNonBonded::getJStreamWidth( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getJStreamWidth";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get j-stream width

   if( _jStreamWidth < 0 ){
      //_getJStreamDimensions( platform );
   }
   return _jStreamWidth;
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
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _jStreamWidth;
}

/** 
 * Get j-stream height
 * 
 * @param platform platform
 *
 * @return  j-stream height 
 *
 */

int BrookNonBonded::getJStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getJStreamHeight";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get j-stream height

   if( _jStreamHeight < 0 ){
      //_getJStreamDimensions( platform );
   }
   return _jStreamHeight;
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
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   return _jStreamHeight;
}

/** 
 * Get j-stream size
 * 
 * @param platform  platform
 *
 * @return  j-stream size
 *
 */

int BrookNonBonded::getJStreamSize( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getJStreamSize";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // get j-stream size

   if( _jStreamSize < 0 ){
      //_getJStreamDimensions( platform );
   }
   return _jStreamSize;
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
   // static const int debug                   = 1;

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

int BrookNonBonded::getPartialForceStreamHeight( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::getPartialForceStreamHeight";
   // static const int debug                   = 1;

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
   // static const int debug                   = 1;

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
   // static const int debug                   = 1;

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

BrookFloatStreamImpl* BrookNonBonded::getExclusionStream( void ) const {
   return _exclusionStream;
}

/** 
 * Get vdw stream 
 *
 * @return  vdw stream
 *
 */

BrookFloatStreamImpl* BrookNonBonded::getVdwStream( void ) const {
   return _vdwStream;
}

/** 
 * Get charge stream 
 *
 * @return  charge stream
 *
 */

BrookFloatStreamImpl* BrookNonBonded::getChargeStream( void ) const {
   return _chargeStream;
}

/** 
 * Get sigma stream 
 *
 * @return  sigma stream
 *
 */

BrookFloatStreamImpl* BrookNonBonded::getSigmaStream( void ) const {
   return _sigmaStream;
}

/** 
 * Get epsilon stream 
 *
 * @return  epsilon stream
 *
 */

BrookFloatStreamImpl* BrookNonBonded::getEpsilonStream( void ) const {
   return _epsilonStream;
}

/** 
 * Get force streams 
 *
 * @return  force streams
 *
 */

BrookFloatStreamImpl** BrookNonBonded::getForceStreams( void ){
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
 * Initialize stream dimensions and streams
 * 
 * @param numberOfAtoms             number of atoms
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::initializeStreams( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::initializeStreams";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int atomStreamSize        = getAtomStreamSize( platform );
   _exclusionStreamWidth     = atomStreamSize;

   int innerUnroll           = getInnerLoopUnroll();
   if( innerUnroll < 1 ){
      std::stringstream message;
      message << methodName << " innerUnrolls=" << innerUnroll << " is less than 1.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory() );

   _exclusionStreamHeight    = _exclusionStreamWidth/innerUnroll; 
   _exclusionStreamSize      = _exclusionStreamWidth*_exclusionStreamHeight;
   StreamImpl* exStream      = brookStreamFactory.createStreamImpl( BrookStreamFactory::NonBondedExclusionStream, _exclusionStreamSize, Stream::Float, platform );
   _exclusionStream          = dynamic_cast<BrookFloatStreamImpl*> ( exStream );


   StreamImpl* vStream       = brookStreamFactory.createStreamImpl( BrookStreamFactory::NonBondedVdwStream, atomStreamSize, Stream::Float3, platform );
   _vdwStream                = dynamic_cast<BrookFloatStreamImpl*> ( vStream );

   StreamImpl* cStream       = brookStreamFactory.createStreamImpl( BrookStreamFactory::NonBondedVdwStream, atomStreamSize, Stream::Float3, platform );
   _chargeStream             = dynamic_cast<BrookFloatStreamImpl*> ( cStream );

   if( getLog() ){
      int atomStreamWidth       = getAtomStreamWidth( platform );
      int atomStreamHeight      = getAtomStreamHeight( platform );
      (void) fprintf( getLog(), "%s initializeVdwStream: stream dimensions: [%d %d] size=%d\n", methodName.c_str(), atomStreamWidth, atomStreamHeight, atomStreamSize );
   }

   if( getLog() ){
      (void) fprintf( getLog(), "%s initializeExclusionStream: stream dimensions: [%d %d] size=%d\n", methodName.c_str(), _exclusionStreamWidth, _exclusionStreamHeight, _exclusionStreamSize );
   }

   return DefaultReturnValue;
}

/** 
 * Set exclusion (4x4)
 * 
 * @param i                         atom i index
 * @param j                         atom j index
 * @param exclusionStreamWidth      exclusion stream width
 * @param exclusion                 array of packed exclusions
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::setExclusion( int i, int j, int exclusionStreamWidth, BrookOpenMMFloat* exclusion ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::setExclusion";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   int index    = ( (j/4)* exclusionStreamWidth + i ); //iindex is the x coordinate!
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
 * @param exclusions            vector of sets containing exclusions (1 set entry for every atom)
 * @param platform              Brook platform
 * @param log                   optional Log file reference
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::initializeExclusions( const std::vector<std::set<int> >& exclusionsVector, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::initializeExclusions";
   static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "In initializeExclusions exclusions vector size=%d\n", exclusionsVector.size() );
   }

   int exclusionStreamSize = getExclusionStreamSize();
   if( exclusionStreamSize < 1 ){
      std::stringstream message;
      message << "BrookNonBonded::initializeExclusions exclusionStreamSize=" << exclusionStreamSize << " is not set.";
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   // allocate and initialize exclusions array

   BrookOpenMMFloat* exclusions = new BrookOpenMMFloat[exclusionStreamSize];
   for( unsigned int ii = 0; ii < (unsigned int) exclusionStreamSize; ii++ ){
      exclusions[ii] = (BrookOpenMMFloat) 210.0;
   }

   // pack in values

   int  exclusionStreamWidth = getExclusionStreamWidth();
   int  numberOfAtoms        = getNumberOfAtoms();
   int  atomStreamSize       = getAtomStreamSize();
   for( unsigned int ii = 0; ii < exclusionsVector.size(); ii++ ){

      set<int> exclusionIndices  = exclusionsVector[ii];
      for( set<int>::const_iterator jj = exclusionIndices.begin(); jj != exclusionIndices.end(); jj++ ){
         setExclusion( ii, *jj, exclusionStreamWidth, exclusions );
      }

      // explicitly exclude junk atoms from interacting w/ atom ii

      for( int jj = numberOfAtoms; jj < atomStreamSize; jj++ ){
         setExclusion( ii, jj, exclusionStreamWidth, exclusions );
      }
   }

   // explicitly exclude junk atoms from interacting w/ all atoms

   for( int ii = numberOfAtoms; ii < atomStreamSize; ii++ ){
      for( int jj = 0; jj < atomStreamSize; jj++ ){
         setExclusion( ii, jj, exclusionStreamWidth, exclusions );
      }
   }

   // load stream

   _exclusionStream->loadFromArray( exclusions ); 
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

int BrookNonBonded::setSigmaEpsilon( double c6, double c12, double* sigma , double* epsilon ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookNonBonded::setSigmaEpsilon";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   if( fabs( c12 ) < 1.0e-04 ){
      *sigma    = 1.0;
      *epsilon  = 0.0;
   } else {
      *sigma    = 0.5*pow( c12/c6, 1.0/6.0 );
      *epsilon  = sqrt( c6*c6/c12 );
   }

   return DefaultReturnValue;
}

/** 
 * Initialize vdw & charge
 * 
 * @param exclusions                vector of sets containing exclusions (1 set entry for every atom)
 * @param platform                  platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookNonBonded::initializeVdwAndCharge( const vector<vector<double> >& nonbondedParameters, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::initializeVdwAndCharge";
   static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "In initializeVdwAndCharge nonbonded vector size=%d\n", nonbondedParameters.size() );
   }

   int atomStreamSize        = getAtomStreamSize();

   // allocate and initialize vdw & charge array

   unsigned int vdwParametersSize  = 2*atomStreamSize;
   BrookOpenMMFloat* vdwParameters = new BrookOpenMMFloat[vdwParametersSize];
   memset( vdwParameters, 0, sizeof( BrookOpenMMFloat )*vdwParametersSize );
   
   BrookOpenMMFloat* charges       = new BrookOpenMMFloat[atomStreamSize];
   memset( charges, 0, sizeof( BrookOpenMMFloat )*atomStreamSize );
   
   // pack in values

   int  numberOfAtoms        = getNumberOfAtoms();
   double sigma, epsilon;
   int trackingIndex         = 0;
   for( unsigned int ii = 0; ii < nonbondedParameters.size(); ii++ ){

      vector<double> nonbondedParameterVector = nonbondedParameters[ii];

      double c6                      = nonbondedParameterVector[0];
      double c12                     = nonbondedParameterVector[1];
      double charge                  = nonbondedParameterVector[2];

      // geometric combination rules

      setSigmaEpsilon( c6, c12, &sigma, &epsilon );

      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) sigma;
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) epsilon;
      charges[ii]                    = (BrookOpenMMFloat) charge;
         
   }

   // set outlier atoms

   for( unsigned int ii = nonbondedParameters.size(); ii < (unsigned int) atomStreamSize; ii++ ){
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) 1.0;
      vdwParameters[trackingIndex++] = (BrookOpenMMFloat) 0.0;
   }

   // load stream

   _vdwStream->loadFromArray( vdwParameters ); 
   _chargeStream->loadFromArray( charges ); 

   delete[] vdwParameters;
   delete[] charges;

   return DefaultReturnValue;
}

/* 
 * Setup of nonbonded ixns
 *
 * @param numberOfAtoms         number of atoms
 * @param nonbondedParameters   vector of nonbonded parameters [atomI][0=c6]
 *                                                             [atomI][1=c12]
 *                                                             [atomI][2=charge]
 * @param platform              Brook platform
 *
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 * */

int BrookNonBonded::setup( int numberOfAtoms, const std::vector<std::vector<double> >& nonbondedParameters,
                           const std::vector<std::set<int> >& exclusions, const BrookPlatform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::setup";
   static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   setNumberOfAtoms( numberOfAtoms );

   initializeStreams( numberOfAtoms, platform );

   initializeExclusions( exclusions, platform );
   initializeVdwAndCharge( nonbondedParameters, platform );
   //initializeJStreamVdw( nonbondedParameters, platform );

   return DefaultReturnValue;
}

/* 
 * Get contents of object
 *
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

std::string BrookNonBonded::getContents( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookNonBonded::getContents";

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

   message << _getLine( tab, "Log:",                  (getLog()             ? Set : NotSet) ); 
   message << _getLine( tab, "ExclusionStream:",      (getExclusionStream() ? Set : NotSet) ); 
   message << _getLine( tab, "VdwStream:",            (getVdwStream()       ? Set : NotSet) ); 
   message << _getLine( tab, "ChargeStream:",         (getChargeStream()    ? Set : NotSet) ); 
   message << _getLine( tab, "SigmaStream:",          (getSigmaStream()     ? Set : NotSet) ); 
   message << _getLine( tab, "EpsilonStream:",        (getEpsilonStream()   ? Set : NotSet) ); 
 
   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      char description[256];
      (void) LOCAL_SPRINTF( description, "PartialForceStream %d", ii );
      message << _getLine( tab, description,  (isForceStreamSet(ii) ? Set : NotSet) ); 
   }
 
#undef LOCAL_SPRINTF

   return message.str();
}
