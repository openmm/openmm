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
#include "BrookStochasticDynamics.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"

// use random number generator

#include "SimTKOpenMMUtilities.h"

using namespace OpenMM;
using namespace std;

/** 
 * Constructor
 * 
 * @param masses atomic masses
 */

BrookStochasticDynamics::BrookStochasticDynamics( const std::vector<double>& masses, uint32_t randomNumberSeed ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::BrookStochasticDynamics";
   BrookOpenMMFloat zero                    = (BrookOpenMMFloat) 0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat) 1.0;

// ---------------------------------------------------------------------------------------

   _sdAtomStreamWidth         = -1;
   _sdAtomStreamHeight        = -1;
   _sdAtomStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _sdStreams[ii]   = NULL;
   }

   // setup inverse masses

   _inverseMasses = new std::vector<BrookOpenMMFloat>;
   for( std::vector<double>::const_interator ii = masses.begin(); ii != masses.end(); ii++ ){
      if( *ii != 0.0 ){
         _inverseMasses->push_back( SQRT( one/(*ii) ) );
      } else {
         _inverseMasses->push_back( zero );
      }
   }

   // randomNumberSeed 

   if( randomNumberSeed ){
      SimTKOpenMMUtilities::setRandomNumberSeed( randomNumberSeed );
   }
}   
 
/** 
 * Destructor
 * 
 */

BrookStochasticDynamics::~BrookStochasticDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookStochasticDynamics::~BrookStochasticDynamics";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _sdStreams[ii];
   }

   delete _inverseMasses;

}

/** 
 * Update fixed parameters
 * 
 * @return  DefaultReturn
 *
 */

int BrookStochasticDynamics::_updateFixedParameters( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::updateFixedParameters";

   static const BrookOpenMMFloat zero       =  0.0;
   static const BrookOpenMMFloat one        =  1.0;
   static const BrookOpenMMFloat two        =  2.0;
   static const BrookOpenMMFloat three      =  3.0;
   static const BrookOpenMMFloat four       =  4.0;
   static const BrookOpenMMFloat half       =  0.5;

   // ---------------------------------------------------------------------------------------

   setStepSize( stepSize );
   setFriction( friction );
   setTemperature( temperature );

   _fixedParameters[GDT]      = getDeltaT()/getTau();
   _fixedParameters[EPH]      = EXP(  half*_fixedParameters[GDT] );
   _fixedParameters[EMH]      = EXP( -half*_fixedParameters[GDT] );
   _fixedParameters[EM]       = EXP(      -_fixedParameters[GDT] );
   _fixedParameters[EP]       = EXP(       _fixedParameters[GDT] );

   if( _fixedParameters[GDT] >= (BrookOpenMMFloat) 0.1 ){

      BrookOpenMMFloat term1  = _fixedParameters[EPH] - one;
                 term1       *= term1;
      _fixedParameters[B]     = _fixedParameters[GDT]*(_fixedParameters[EP] - one) - four*term1;

      _fixedParameters[C]     = _fixedParameters[GDT] - three + four*_fixedParameters[EMH] - _fixedParameters[EM];
      _fixedParameters[D]     = two - _fixedParameters[EPH] - _fixedParameters[EMH];

    } else {

      // this has not been debugged

      BrookOpenMMFloat term1        = half*_fixedParameters[GDT];
      BrookOpenMMFloat term2        = term1*term1;
      BrookOpenMMFloat term4        = term2*term2;

      BrookOpenMMFloat third        = (RealOpenMM) ( 1.0/3.0 );
      BrookOpenMMFloat o7_9         = (RealOpenMM) ( 7.0/9.0 );
      BrookOpenMMFloat o1_12        = (RealOpenMM) ( 1.0/12.0 );
      BrookOpenMMFloat o17_90       = (RealOpenMM) ( 17.0/90.0 );
      BrookOpenMMFloat o7_30        = (RealOpenMM) ( 7.0/30.0 );
      BrookOpenMMFloat o31_1260     = (RealOpenMM) ( 31.0/1260.0 );
      BrookOpenMMFloat o_360        = (RealOpenMM) ( 1.0/360.0 );

      _fixedParameters[B]     = term4*( third  + term1*( third + term1*( o17_90 + term1*o7_9 )));
      _fixedParameters[C]     = term2*term1*( two*third + term1*( -half + term1*( o7_30 + term1*(-o1_12 + term1*o31_1260 ))));
      _fixedParameters[D]     = term2*( -one + term2*(-o1_12 - term2*o_360));
   }    

   RealOpenMM kT        = ((RealOpenMM) BOLTZ)*getTemperature();

   _fixedParameters[V]  = SQRT( kT*( one - _fixedParameters[EM]) );
   _fixedParameters[X]  = getTau()*SQRT( kT*_fixedParameters[C] );
   _fixedParameters[Yv] = SQRT( kT*_fixedParameters[B]/_fixedParameters[C] );
   _fixedParameters[Yx] = getTau()*SQRT( kT*_fixedParameters[B]/(one - _fixedParameters[EM]) );

   return DefaultReturn;

};


/** 
 * Update fixed parameters
 * 
 * @return  DefaultReturn
 *
 */

int BrookStochasticDynamics::_updateSdStreams( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::updateFixedParameters";

   static const BrookOpenMMFloat zero       =  0.0;
   static const BrookOpenMMFloat one        =  1.0;
   static const BrookOpenMMFloat two        =  2.0;
   static const BrookOpenMMFloat three      =  3.0;
   static const BrookOpenMMFloat four       =  4.0;
   static const BrookOpenMMFloat half       =  0.5;

   // ---------------------------------------------------------------------------------------

   
   return DefaultReturn;

};

/** 
 * Update parameters
 * 
 * @param  temperature     temperature
 * @param  friction        friction
 *
 * @return   solute dielectric
 *
 */

int BrookStochasticDynamics::updateParameters( double temperature, double friction, double stepSize ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookStochasticDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   setStepSize( stepSize );
   setFriction( friction );
   setTemperature( temperature );

   _updateFixedParameters( );
   _updateSdStreams( );

   return DefaultReturn;

};

/**---------------------------------------------------------------------------------------

   Get tau

   @return tau

   --------------------------------------------------------------------------------------- */

RealOpenMM BrookStochasticDynamics::getTau( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookStochasticDynamics::getTau";

   // ---------------------------------------------------------------------------------------

   return _tau;
}

/**---------------------------------------------------------------------------------------

   Get array of fixed parameters indexed by 'FixedParameters' enums

   @return array

   --------------------------------------------------------------------------------------- */
   
const RealOpenMM* BrookStochasticDynamics::getFixedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookStochasticDynamics::getFixedParameters";

   // ---------------------------------------------------------------------------------------

   return _fixedParameters;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamSize( void ) const {
   return _sdAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamWidth( void ) const {
   return _sdAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookStochasticDynamics::getStochasticDynamicsAtomStreamHeight( void ) const {
   return _sdAtomStreamHeight;
}

/** 
 * Get SDPC1 stream 
 *
 * @return  SDPC1 stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSDPC1( void ) const {
   return _sdStreams[SDPC1Stream];
}

/** 
 * Get SDPC2 stream 
 *
 * @return  SDPC2 stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSDPC2( void ) const {
   return _sdStreams[SDPC2Stream];
}

/** 
 * Get SD2X stream 
 *
 * @return  SD2X stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSD2X( void ) const {
   return _sdStreams[SD2XStream];
}

/** 
 * Get SD1V stream 
 *
 * @return  SD1V stream
 *
 */

BrookFloatStreamInternal* BrookStochasticDynamics::getSD1V( void ) const {
   return _sdStreams[SD1VStream];
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

int BrookStochasticDynamics::initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _sdAtomStreamSize     = getAtomStreamSize( platform );
   _sdAtomStreamWidth    = getAtomStreamWidth( platform );
   _sdAtomStreamHeight   = getAtomStreamHeight( platform );

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

int BrookStochasticDynamics::initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::initializeStreams";
   static const double dangleValue          = 0.0;

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize   = getStochasticDynamicsAtomStreamSize();
   int sdAtomStreamWidth  = getStochasticDynamicsAtomStreamWidth();

    _sdStreams[SDPC1Stream]                                   = new BrookFloatStreamInternal( BrookCommon::SDPC1Stream,
                                                                                              sdAtomStreamSize, sdAtomStreamWidth,
                                                                                              BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SDPC2Stream]                                   = new BrookFloatStreamInternal( BrookCommon::SDPC2Stream,
                                                                                              sdAtomStreamSize, sdAtomStreamWidth,
                                                                                              BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SD2XStream]                                    = new BrookFloatStreamInternal( BrookCommon::SD2XStream,
                                                                                              sdAtomStreamSize, sdAtomStreamWidth,
                                                                                              BrookStreamInternal::Float3, dangleValue );

    _sdStreams[ObcBornRadiiStream]                            = new BrookFloatStreamInternal( BrookCommon::SD1VStream,
                                                                                              sdAtomStreamSize, sdAtomStreamWidth,
                                                                                              BrookStreamInternal::Float3, dangleValue );

   return DefaultReturnValue;
}

/** 
 * Update sd streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookStochasticDynamics::_updateSdStreams( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::updateSdStreams";
   static const double dangleValue          = 0.0;

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize          = getStochasticDynamicsAtomStreamSize();

   BrookOpenMMFloat sdpc[2];
   for( int ii = 0; ii < 2; ii++ ){
      sdpc[ii] = new BrookOpenMMFloat[2*sdAtomStreamSize];
      memset( sdpc[ii], 0, 2*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 
   }

   const RealOpenMM* fixedParameters = getFixedParameters( );
   int numberOfAtoms                 = getNumberOfAtoms();
   int index                         = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){

      sdpc[0][index]      = _inverseMasses[ii]*( static_cast<BrookOpenMMFloat> (fixedParameters[Yv]) );
      sdpc[0][index+1]    = _inverseMasses[ii]*( static_cast<BrookOpenMMFloat> (fixedParameters[V])  );

      sdpc[1][index]      = _inverseMasses[ii]*( static_cast<BrookOpenMMFloat> (fixedParameters[Yx]) );
      sdpc[1][index+1]    = _inverseMasses[ii]*( static_cast<BrookOpenMMFloat> (fixedParameters[X])  );

      index              += 2;
   }

   _sdStreams[SDPC1Stream]->loadFromArray( sdpc[0] );
   _sdStreams[SDPC2Stream]->loadFromArray( sdpc[1] );

   for( int ii = 0; ii < 2; ii++ ){
      delete[] sdpc[ii];
   }

   // initialize SD2X

   sd2x = new BrookOpenMMFloat[3*sdAtomStreamSize];
   SimTKOpenMMUtilities::setRandomNumberSeed( (uint32_t) getRandomNumberSeed() );

   memset( sd2x, 0, 3*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      sd2x[index]        = _inverseMasses[ii]*fixedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+1]      = _inverseMasses[ii]*fixedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+2]      = _inverseMasses[ii]*fixedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
   }
   
   _sdStreams[SD2XStream]->loadFromArray( sd2x );

   delete[] sd2x;

   return DefaultReturnValue;

}

/*  
 * Setup of StochasticDynamics parameters
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
    
int BrookStochasticDynamics::setup( const std::vector<std::vector<double> >& vectorOfAtomParameters, 
                      double solventDielectric, double soluteDielectric, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const int atomParametersSize      = 4; 
   static const int maxErrors               = 20; 
   static const std::string methodName      = "BrookStochasticDynamics::setup";

// ---------------------------------------------------------------------------------------

   int numberOfAtoms  = (int) vectorOfAtomParameters.size();
   setNumberOfAtoms( numberOfAtoms );

   _solventDielectric = solventDielectric;
   _soluteDielectric  = soluteDielectric;

   // initialize stream sizes and then Brook streams

   initializeStreamSizes( numberOfAtoms, platform );
   initializeStreams( platform );

   BrookOpenMMFloat* radiiAndCharge         = new BrookOpenMMFloat[getNumberOfAtoms()*2];
   BrookOpenMMFloat* scaledRadiiAndOffset   = new BrookOpenMMFloat[getNumberOfAtoms()*2];

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

         radiiAndCharge[streamIndex]              = static_cast<BrookOpenMMFloat> (radius);
         radiiAndCharge[streamIndex+1]            = static_cast<BrookOpenMMFloat> (charge);

         scaledRadiiAndOffset[streamIndex]        = static_cast<BrookOpenMMFloat> (radius*scalingFactor);
         scaledRadiiAndOffset[streamIndex+1]      = static_cast<BrookOpenMMFloat> (radius - dielectricOffset);

      }

      vectorIndex++;
   }

   // throw exception if errors detected

   if( errors ){
      throw OpenMMException( message.str() );
   }

   // load streams

   _sdStreams[ObcAtomicRadiiStream]->loadFromArray( radiiAndCharge );
   _sdStreams[ObcScaledAtomicRadiiStream]->loadFromArray( scaledRadiiAndOffset );

   delete[] radiiAndCharge;
   delete[] scaledRadiiAndOffset;

   // setup for Born radii

   ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   obcParameters->setAtomicRadii( atomicRadii, SimTKOpenMMCommon::MdUnits);

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

int BrookStochasticDynamics::initializePartialForceStreamSize( int atomStreamSize, int atomStreamWidth ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::initializePartialForceStreamSize";
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

std::string BrookStochasticDynamics::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookStochasticDynamics::getContentsString";

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
/*
   message << _getLine( tab, "ExclusionStream:",      (getExclusionStream()    ? Set : NotSet) ); 
   message << _getLine( tab, "VdwStream:",            (getOuterVdwStream()     ? Set : NotSet) ); 
   message << _getLine( tab, "ChargeStream:",         (getChargeStream()       ? Set : NotSet) ); 
   message << _getLine( tab, "SigmaStream:",          (getInnerSigmaStream()   ? Set : NotSet) ); 
   message << _getLine( tab, "EpsilonStream:",        (getInnerEpsilonStream() ? Set : NotSet) ); 
*/
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _sdStreams[ii] ){
         message << _sdStreams[ii]->getContentsString( );
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
      if( _sdForceStreams[ii] ){
         message << _sdForceStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
