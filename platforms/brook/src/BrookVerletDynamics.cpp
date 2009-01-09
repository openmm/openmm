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
#include "BrookVerletDynamics.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"
#include "gpu/kshakeh.h"
#include "gpu/kupdatemd.h"
#include "gpu/kcommon.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookVerletDynamics::BrookVerletDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVerletDynamics::BrookVerletDynamics";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;
   BrookOpenMMFloat oneMinus                = (BrookOpenMMFloat) -1.0;

// ---------------------------------------------------------------------------------------

   _numberOfParticles                 = -1;
   _internalStepCount                 = 0;

   // mark stream dimension variables as unset

   _verletParticleStreamWidth         = -1;
   _verletParticleStreamHeight        = -1;
   _verletParticleStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _verletStreams[ii]   = NULL;
   }

   _stepSize                          = oneMinus;

   // setup inverse sqrt masses

   _inverseMasses                     = NULL;

}   
 
/** 
 * Destructor
 * 
 */

BrookVerletDynamics::~BrookVerletDynamics( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVerletDynamics::~BrookVerletDynamics";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _verletStreams[ii];
   }

   delete[] _inverseMasses;

}

/** 
 * Get stepSize
 * 
 * @return  stepSize
 *
 */

BrookOpenMMFloat BrookVerletDynamics::getStepSize( void ) const {
   return _stepSize;
}

/** 
 * Set stepSize
 * 
 * @param   stepSize
 *
 * @return      DefaultReturnValue
 *
 */

int BrookVerletDynamics::_setStepSize( BrookOpenMMFloat stepSize ){
   _stepSize = stepSize;
   return DefaultReturnValue;
}

/** 
 * Update parameters -- only way parameters can be set
 * 
 * @param  step size       step size
 *
 * @return   DefaultReturnValue
 *
 */

int BrookVerletDynamics::updateParameters( double stepSize ){

   // ---------------------------------------------------------------------------------------

   static int showUpdate               = 1;
   static int maxShowUpdate            = 3;
   static const std::string methodName = "BrookVerletDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize( (BrookOpenMMFloat) stepSize );
   _updateVerletStreams( );

   // show update

   if( showUpdate && getLog() && (showUpdate++ < maxShowUpdate) ){
      std::string contents = getContentsString( );
      (void) fprintf( getLog(), "%s contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( getLog() );

   }

   return DefaultReturnValue;

}

/** 
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookVerletDynamics::getVerletDynamicsParticleStreamSize( void ) const {
   return _verletParticleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookVerletDynamics::getVerletDynamicsParticleStreamWidth( void ) const {
   return _verletParticleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookVerletDynamics::getVerletDynamicsParticleStreamHeight( void ) const {
   return _verletParticleStreamHeight;
}

/** 
 * Get VPrime stream 
 *
 * @return  Vprime stream
 *
 */

BrookFloatStreamInternal* BrookVerletDynamics::getVPrimeStream( void ) const {
   return _verletStreams[VPrimeStream];
}

/** 
 * Get XPrime stream 
 *
 * @return  Xprime stream
 *
 */

BrookFloatStreamInternal* BrookVerletDynamics::getXPrimeStream( void ) const {
   return _verletStreams[XPrimeStream];
}

/** 
 * Get InverseMass stream 
 *
 * @return  inverse mass stream
 *
 */

BrookFloatStreamInternal* BrookVerletDynamics::getInverseMassStream( void ) const {
   return _verletStreams[InverseMassStream];
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfParticles             number of particles
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookVerletDynamics::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVerletDynamics::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _verletParticleStreamSize     = getParticleStreamSize( platform );
   _verletParticleStreamWidth    = getParticleStreamWidth( platform );
   _verletParticleStreamHeight   = getParticleStreamHeight( platform );

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

int BrookVerletDynamics::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVerletDynamics::_initializeStreams";

   BrookOpenMMFloat dangleValue         = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int sdParticleStreamSize             = getVerletDynamicsParticleStreamSize();
   int sdParticleStreamWidth            = getVerletDynamicsParticleStreamWidth();

    _verletStreams[VPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::VPrimeStream,
                                                                        sdParticleStreamSize, sdParticleStreamWidth,
                                                                        BrookStreamInternal::Float3, dangleValue );

    _verletStreams[XPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::XPrimeStream,
                                                                        sdParticleStreamSize, sdParticleStreamWidth,
                                                                        BrookStreamInternal::Float3, dangleValue );

    _verletStreams[InverseMassStream]   = new BrookFloatStreamInternal( BrookCommon::InverseMassStream,
                                                                        sdParticleStreamSize, sdParticleStreamWidth,
                                                                        BrookStreamInternal::Float, dangleValue );

   return DefaultReturnValue;
}

/** 
 * Update sd streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookVerletDynamics::_updateVerletStreams( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName       = "BrookVerletDynamics::_updateVerletStreams";

// ---------------------------------------------------------------------------------------

   int particleStreamSize                    = getVerletDynamicsParticleStreamSize();

   BrookOpenMMFloat* inverseMass = new BrookOpenMMFloat[particleStreamSize];
   memset( inverseMass, 0, particleStreamSize*sizeof( BrookOpenMMFloat ) ); 

   int numberOfParticles                         = getNumberOfParticles();
   for( int ii = 0; ii < numberOfParticles; ii++ ){
      inverseMass[ii]  = _inverseMasses[ii];
   }

   _verletStreams[InverseMassStream]->loadFromArray( inverseMass );

   delete[] inverseMass;

   return DefaultReturnValue;

}

/** 
 * Set masses 
 * 
 * @param masses             particle masses
 *
 */

int BrookVerletDynamics::_setInverseMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVerletDynamics::_setInverseSqrtMasses";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   // setup inverse sqrt masses

   _inverseMasses = new BrookOpenMMFloat[masses.size()];
   int index          = 0;
   for( std::vector<double>::const_iterator ii = masses.begin(); ii != masses.end(); ii++, index++ ){
      if( *ii != 0.0 ){
         BrookOpenMMFloat value    = static_cast<BrookOpenMMFloat>(*ii);
         _inverseMasses[index] = ( one/value );
      } else {
         _inverseMasses[index] = zero;
      }
   }

   return DefaultReturnValue;
}   
 
/*  
 * Setup of VerletDynamics parameters
 *
 * @param masses                masses
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookVerletDynamics::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVerletDynamics::setup";

// ---------------------------------------------------------------------------------------

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfParticles  = (int) masses.size();
   setNumberOfParticles( numberOfParticles );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, platform );
   _initializeStreams( platform );

   _setInverseMasses( masses );

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

std::string BrookVerletDynamics::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVerletDynamics::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%.5f", getStepSize() );
   message << _getLine( tab, "Step size:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _verletStreams[ii] ){
         message << _verletStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
/** 
 * Update
 * 
 * @param  positions                   particle positions
 * @param  velocities                  particle velocities
 * @param  forces                      particle forces
 * @param  brookShakeAlgorithm         BrookShakeAlgorithm reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookVerletDynamics::update( BrookStreamImpl& positionStream, BrookStreamImpl& velocityStream,
                                 const BrookStreamImpl& forceStreamC,
                                 BrookShakeAlgorithm& brookShakeAlgorithm ){

// ---------------------------------------------------------------------------------------

   static std::string methodName  = "\nBrookVerletDynamics::update";
   static int printOn             = 0;
   FILE* log;

// ---------------------------------------------------------------------------------------

   _internalStepCount++;

//setLog( stderr );
   printOn = (printOn && getLog()) ? printOn : 0; 

   BrookStreamImpl& forceStream = const_cast<BrookStreamImpl&> (forceStreamC);

   if( printOn ){

      static int showAux = 1;
      log                = getLog();
      if( showAux ){
         showAux = 0; 

/*
         std::string contents = _brookVelocityCenterOfMassRemoval->getContentsString( );
         (void) fprintf( log, "%s VelocityCenterOfMassRemoval contents\n%s", methodName, contents.c_str() );
*/
         (void) fprintf( log, "%s step=%d Shake contents\n%s", methodName.c_str(), _internalStepCount, brookShakeAlgorithm.getContentsString().c_str() );
         (void) fflush( log );
      }    

   }

   // To Shake or not to Shake

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){

      // integration step
   
      kupdate_md1( (float) getStepSize(),
                   positionStream.getBrookStream(),
                   velocityStream.getBrookStream(),
                   forceStream.getBrookStream(),
                   getInverseMassStream()->getBrookStream(),
                   getXPrimeStream()->getBrookStream() 
                 );
    
      // diagnostics
   
      if( printOn ){
         (void) fprintf( log, "\n%s step=%d Post kupdate_md1: particleStrW=%3d step=%.5f",
                                   methodName.c_str(), _internalStepCount, getVerletDynamicsParticleStreamWidth(), getStepSize() );
   
         (void) fprintf( log, "\nInverseMassStream %d\n", _internalStepCount );
         getInverseMassStream()->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nForceStream %d\n", _internalStepCount );
         BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
         std::vector<std::vector<double> > forceStatistics;
         brookStreamInternalF->getStatistics( forceStatistics, getNumberOfParticles() );

         std::stringstream tag;
         tag << _internalStepCount << " Fxx ";
         std::string stats = brookStreamInternalF->printStatistics( tag.str(), forceStatistics );
         (void) fprintf( log, "\nStep %d Force stats:\n%s", _internalStepCount, stats.c_str() );

         brookStreamInternalF->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalV   = velocityStream.getBrookStreamImpl();
         std::vector<std::vector<double> > velocityStatistics;
         brookStreamInternalV->getStatistics( velocityStatistics, getNumberOfParticles() );
         std::stringstream tagV;
         tagV << _internalStepCount << " Vxx ";
         stats = brookStreamInternalPos->printStatistics( tagV.str(), velocityStatistics );
         (void) fprintf( log, "\nStep %d Velocity stats:\n%s", _internalStepCount, stats.c_str() );
         (void) fprintf( log, "\nVelocityStream %d\n", _internalStepCount );
         brookStreamInternalV->printToFile( log );
   
         (void) fprintf( log, "\nXPrimeStream %d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 
      }   

      // Shake

      kshakeh_fix1( 
                    (float) brookShakeAlgorithm.getMaxIterations(),
                    (float) getVerletDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getShakeTolerance(),
                    brookShakeAlgorithm.getShakeParticleIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeParticleParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      if( printOn ){

         (void) fprintf( log, "\n%s Post kshakeh_fix1: particleStrW=%3d", methodName.c_str(), getVerletDynamicsParticleStreamWidth() );
   
         (void) fprintf( log, "\nShakeParticleIndicesStream %d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeParticleIndicesStream()->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nXPrimeStream %d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeParticleParameterStream %d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeParticleParameterStream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons0\n" );
         brookShakeAlgorithm.getShakeXCons0Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons1\n" );
         brookShakeAlgorithm.getShakeXCons1Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons2\n" );
         brookShakeAlgorithm.getShakeXCons2Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons3\n" );
         brookShakeAlgorithm.getShakeXCons3Stream()->printToFile( log ); 

      }   

      // Shake gather
   
      kshakeh_update1_fix1( 
                    (float) brookShakeAlgorithm.getShakeConstraintStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream() );
                    //positionStream.getBrookStream() );
   
      if( printOn ){

         (void) fprintf( log, "\n%s Post kshakeh_update1_fix1: step=%d ShakeConstraintStreamWidth=%3d",
                                   methodName.c_str(), _internalStepCount, brookShakeAlgorithm.getShakeConstraintStreamWidth() );
   
         (void) fprintf( log, "\nShakeInverseMapStream %d\n", _internalStepCount );
         brookShakeAlgorithm.getShakeInverseMapStream()->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nXPrimeStream %d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons0\n" );
         brookShakeAlgorithm.getShakeXCons0Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons1\n" );
         brookShakeAlgorithm.getShakeXCons1Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons2\n" );
         brookShakeAlgorithm.getShakeXCons2Stream()->printToFile( log ); 

         (void) fprintf( log, "\nShakeXCons3\n" );
         brookShakeAlgorithm.getShakeXCons3Stream()->printToFile( log ); 

      }   

      // second integration step

      float inverseStepSize = 1.0f/getStepSize();
      kupdate_md2( inverseStepSize,
                   getXPrimeStream()->getBrookStream(), 
                   positionStream.getBrookStream(),
                   velocityStream.getBrookStream(),
                   positionStream.getBrookStream() );

      if( printOn ){

         (void) fprintf( log, "\n%s step=%d Post kupdate_md2: inverseStepSize=%3e",
                                   methodName.c_str(), _internalStepCount, inverseStepSize );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         brookShakeAlgorithm.checkConstraints( brookStreamInternalPos, log, 0.0001f );
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );

         (void) fprintf( log, "\nXPrimeStream %d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 

         brookStreamInternalPos = velocityStream.getBrookStreamImpl();
         (void) fprintf( log, "\nVelocityStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log ); 
      }   

   } else {

      kupdateMdNoShake( getStepSize(),
                        positionStream.getBrookStream(),
                        velocityStream.getBrookStream(),
                        forceStream.getBrookStream(),
                        getInverseMassStream()->getBrookStream(),
                        velocityStream.getBrookStream(),
                        positionStream.getBrookStream() );
   }

//_brookVelocityCenterOfMassRemoval->removeVelocityCenterOfMass( velocities );

   return DefaultReturnValue;

}
