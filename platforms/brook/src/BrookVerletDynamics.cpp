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

// use random number generator

#include "../../reference/src/SimTKUtilities/SimTKOpenMMUtilities.h"

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

   // mark stream dimension variables as unset

   _verletParticleStreamWidth         = -1;
   _verletParticleStreamHeight        = -1;
   _verletParticleStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _verletStreams[ii]   = NULL;
   }

   _stepSize                          = oneMinus;

   // setup inverse sqrt masses

   _inverseMasses             = NULL;

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

   static int showUpdate         = 1;
   static int maxShowUpdate      = 3;
   static const char* methodName = "\nBrookVerletDynamics::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize( (BrookOpenMMFloat) stepSize );
   _updateVerletStreams( );

   // show update

   if( showUpdate && getLog() && (showUpdate++ < maxShowUpdate) ){
      std::string contents = getContentsString( );
      (void) fprintf( getLog(), "%s contents\n%s", methodName, contents.c_str() );
      (void) fflush( getLog() );

   }

   return DefaultReturnValue;

};

/** 
 * Update
 * 
 * @param  positions                   particle positions
 * @param  velocities                  particle velocities
 * @param  forces                      particle forces
 * @param  brookShakeAlgorithm         BrookShakeAlgorithm reference
 * @param  brookRandomNumberGenerator  BrookRandomNumberGenerator reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookVerletDynamics::update( Stream& positions, Stream& velocities,
                                 const Stream& forces,
                                 BrookShakeAlgorithm& brookShakeAlgorithm ){

// ---------------------------------------------------------------------------------------

   // unused Shake parameter

   float omega                   = 1.0f;

   static const char* methodName = "\nBrookVerletDynamics::update";

   static const int PrintOn      = 0;

// ---------------------------------------------------------------------------------------

   BrookStreamImpl& positionStream                     = dynamic_cast<BrookStreamImpl&>       (positions.getImpl());
   BrookStreamImpl& velocityStream                     = dynamic_cast<BrookStreamImpl&>       (velocities.getImpl());
   const BrookStreamImpl& forceStreamC                 = dynamic_cast<const BrookStreamImpl&> (forces.getImpl());
   BrookStreamImpl& forceStream                        = const_cast<BrookStreamImpl&> (forceStreamC);

   if( (1 || PrintOn) && getLog() ){

      static int showAux = 1;

      (void) fprintf( getLog(), "%s shake=%d\n", methodName, brookShakeAlgorithm.getNumberOfConstraints() );
      (void) fflush( getLog() );

      if( showAux ){
         showAux = 0; 

/*
         std::string contents = _brookVelocityCenterOfMassRemoval->getContentsString( );
         (void) fprintf( getLog(), "%s VelocityCenterOfMassRemoval contents\n%s", methodName, contents.c_str() );
*/
         std::string contents   = brookShakeAlgorithm.getContentsString( );
         (void) fprintf( getLog(), "%s Shake contents\n%s", methodName, contents.c_str() );

         (void) fflush( getLog() );
      }    

   }

   // integration step

   kupdate_md_verlet( (float) getStepSize(),
                       positionStream.getBrookStream(),
                       velocityStream.getBrookStream(),
                       forceStream.getBrookStream(),
                       getInverseMassStream()->getBrookStream(),
                       velocityStream.getBrookStream(),
                       getXPrimeStream()->getBrookStream() 
                    );
 
   // diagnostics

   if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "\nPost kupdate_md_verlet: particleStrW=%3d step=%.5f",
                                getVerletDynamicsParticleStreamWidth(), getStepSize() );

      (void) fprintf( getLog(), "\nInverseMassStream\n" );
      getInverseMassStream()->printToFile( getLog() );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStream\n" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

      BrookStreamInternal* brookStreamInternalV   = velocityStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nVelocityStream\n" );
      brookStreamInternalV->printToFile( getLog() );

      (void) fprintf( getLog(), "\nXPrimeStream\n" );
      getXPrimeStream()->printToFile( getLog() ); 
   }   

   // second Shake step

   if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){

      kshakeh_fix1( 
                    10.0f,
                    (float) getVerletDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getInverseHydrogenMass(),
                    omega, 
                    brookShakeAlgorithm.getShakeParticleIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeParticleParameterStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream() );
   
      if( (1|| PrintOn) && getLog() ){

         (void) fprintf( getLog(), "\nPost kshakeh_fix1: particleStrW=%3d invMass_H=%.5f",
                                   getVerletDynamicsParticleStreamWidth(), brookShakeAlgorithm.getInverseHydrogenMass() );
   
         (void) fprintf( getLog(), "\nShakeParticleIndicesStream\n" );
         brookShakeAlgorithm.getShakeParticleIndicesStream()->printToFile( getLog() );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( getLog() );
   
         (void) fprintf( getLog(), "\nXPrimeStream\n" );
         getXPrimeStream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeParticleParameterStream\n" );
         brookShakeAlgorithm.getShakeParticleParameterStream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons0\n" );
         brookShakeAlgorithm.getShakeXCons0Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons1\n" );
         brookShakeAlgorithm.getShakeXCons1Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons2\n" );
         brookShakeAlgorithm.getShakeXCons2Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons3\n" );
         brookShakeAlgorithm.getShakeXCons3Stream()->printToFile( getLog() ); 

      }   

      // second Shake gather
   
      kshakeh_update2_fix1( 
                    (float) getVerletDynamicsParticleStreamWidth(),
                    brookShakeAlgorithm.getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm.getShakeXCons3Stream()->getBrookStream(),
                    positionStream.getBrookStream() );
   
      if( ( 1 || PrintOn) && getLog() ){

         (void) fprintf( getLog(), "\nPost kshakeh_update2_fix1: particleStrW=%3d",
                                   getVerletDynamicsParticleStreamWidth() );
   
         (void) fprintf( getLog(), "\nShakeInverseMapStream\n" );
         brookShakeAlgorithm.getShakeInverseMapStream()->printToFile( getLog() );
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
         (void) fprintf( getLog(), "\nPositionStream\n" );
         brookStreamInternalPos->printToFile( getLog() );
   
         (void) fprintf( getLog(), "\nXPrimeStream\n" );
         getXPrimeStream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons0\n" );
         brookShakeAlgorithm.getShakeXCons0Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons1\n" );
         brookShakeAlgorithm.getShakeXCons1Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons2\n" );
         brookShakeAlgorithm.getShakeXCons2Stream()->printToFile( getLog() ); 

         (void) fprintf( getLog(), "\nShakeXCons3\n" );
         brookShakeAlgorithm.getShakeXCons3Stream()->printToFile( getLog() ); 

      }   

   } else {
      //kadd3( getXPrimeStream()->getBrookStream(), positionStream.getBrookStream() );
      ksetStr3( getXPrimeStream()->getBrookStream(), positionStream.getBrookStream() );
   }

//_brookVelocityCenterOfMassRemoval->removeVelocityCenterOfMass( velocities );

   return DefaultReturnValue;

};

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
