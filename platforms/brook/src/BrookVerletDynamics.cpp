/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include <sstream>
#include "BrookVerletDynamics.h"
#include "BrookPlatform.h"
#include "openmm/OpenMMException.h"
#include "BrookStreamImpl.h"
#include "kernels/kshakeh.h"
#include "kernels/kupdatemd.h"
#include "kernels/kcommon.h"

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

   BrookOpenMMFloat zero                    = static_cast<BrookOpenMMFloat>(  0.0 );
   BrookOpenMMFloat one                     = static_cast<BrookOpenMMFloat>(  1.0 );
   BrookOpenMMFloat oneMinus                = static_cast<BrookOpenMMFloat>( -1.0 );

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

   const BrookPlatform& brookPlatform       = dynamic_cast<const BrookPlatform&> (platform);
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
   if( printOn && getLog() ){
       log    = getLog();; 
   } else {
      printOn = 0;
   }

   BrookStreamImpl& forceStream = const_cast<BrookStreamImpl&> (forceStreamC);

   if( printOn ){

      static int showAux = 1;
      if( showAux ){
         showAux = 0; 

/*
         std::string contents = _brookVelocityCenterOfMassRemoval->getContentsString( );
         (void) fprintf( log, "%s VelocityCenterOfMassRemoval contents\n%s", methodName, contents.c_str() );
*/
         (void) fprintf( log, "%s step=%d \n%s\n\nShake contents\n%s", methodName.c_str(), _internalStepCount, getContentsString().c_str(),
                         brookShakeAlgorithm.getContentsString().c_str() );
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
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );
   
         (void) fprintf( log, "\nForceStream %d\n", _internalStepCount );
         BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamInternal();
         std::vector<std::vector<double> > forceStatistics;
         brookStreamInternalF->getStatistics( forceStatistics, getNumberOfParticles() );

         std::stringstream tag;
         tag << _internalStepCount << " Fxx ";
         std::string stats = brookStreamInternalF->printStatistics( tag.str(), forceStatistics );
         (void) fprintf( log, "\nStep %d Force stats:\n%s", _internalStepCount, stats.c_str() );

         brookStreamInternalF->printToFile( log );
   
         BrookStreamInternal* brookStreamInternalV   = velocityStream.getBrookStreamInternal();
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
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
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
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
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
   
         BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
         std::string violationString;
         brookShakeAlgorithm.checkConstraints( brookStreamInternalPos, violationString, 0.0001f );
         (void) fprintf( log, "Shake: %s\n", violationString.c_str() );
         (void) fprintf( log, "\nPositionStream %d\n", _internalStepCount );
         brookStreamInternalPos->printToFile( log );

         (void) fprintf( log, "\nXPrimeStream %d\n", _internalStepCount );
         getXPrimeStream()->printToFile( log ); 

         brookStreamInternalPos = velocityStream.getBrookStreamInternal();
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

   // diagnostics

   if( 0 && (_internalStepCount % 100) == 0 ){
      FILE*	log1     = stderr;
      float  epsilon = 1.0e-01f;

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamInternal();
      BrookStreamInternal* brookStreamInternalVel  = velocityStream.getBrookStreamInternal();
      BrookStreamInternal* brookStreamInternalFrc  = forceStream.getBrookStreamInternal();

      // check for nan and infinities

      int coordinateNans                           = brookStreamInternalPos->checkForNans( );
      int velocityNans                             = brookStreamInternalVel->checkForNans( );
      int forceNans                                = brookStreamInternalFrc->checkForNans( );
      int abort                                    = abs( coordinateNans ) + abs( velocityNans ) + abs( forceNans );

      // Shake violations

      std::string violationString;
      if( brookShakeAlgorithm.getNumberOfConstraints() > 0 ){
         int constraintViolations                     = brookShakeAlgorithm.checkConstraints( brookStreamInternalPos, violationString, 0.0001f );
         abort                                       += abs( constraintViolations );
      } else {
         violationString                              = "Shake not active";
      }

      // force sums ~ 0?
      std::vector<float> sums;
      brookStreamInternalFrc->sumColumns( sums );

      // check if should abort

      (void) fprintf( log1, "%d Nans: x=%d v=%d f=%d ", _internalStepCount, coordinateNans, velocityNans, forceNans );
      (void) fprintf( log1, " Fsum[" );
      for( int ii = 0; ii < 3; ii++ ){
         if( fabsf( sums[ii] ) > epsilon ){
            abort++;
         }
         (void) fprintf( log1, "%12.4e ", sums[ii] );
      }
      (void) fprintf( log1, "] %s abort=%d\n", violationString.c_str(), abort );
      (void) fflush( log1 );
      if( abort ){

         (void) fprintf( log1, "Aborting:\n" );

         brookStreamInternalPos->printToFile( log1 );
         brookStreamInternalVel->printToFile( log1 );
         brookStreamInternalFrc->printToFile( log1 );

         exit(1);
      }
  
/*
      std::vector<std::vector<double> > velocityStatistics;
      brookStreamInternalPos->getStatistics( velocityStatistics, getNumberOfParticles() );
      std::stringstream tagV;
      tagV << _internalStepCount << " Vxx "; 
      std::string stats = brookStreamInternalPos->printStatistics( tagV.str(), velocityStatistics );
      (void) fprintf( log1, "\nStep %d Velocity stats:\n%s", _internalStepCount, stats.c_str() );
*/
      //removeCom( brookStreamInternalPos, getInverseMassStream() );
   }
//_brookVelocityCenterOfMassRemoval->removeVelocityCenterOfMass( velocities );

   return DefaultReturnValue;

}
