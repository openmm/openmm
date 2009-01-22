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
#include "BrookVelocityCenterOfMassRemoval.h"
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"
#include "kernels/kcom.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookVelocityCenterOfMassRemoval::BrookVelocityCenterOfMassRemoval( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::BrookVelocityCenterOfMassRemoval";

   BrookOpenMMFloat zero                    = static_cast<BrookOpenMMFloat>( 0.0 );

// ---------------------------------------------------------------------------------------

   _numberOfParticles           = -1;

   // mark stream dimension variables as unset

   _particleStreamWidth         = -1;
   _particleStreamHeight        = -1;
   _particleStreamSize          = -1;

   _totalInverseMass            = zero;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _streams[ii]   = NULL;
   }

}   
 
/** 
 * Destructor
 * 
 */

BrookVelocityCenterOfMassRemoval::~BrookVelocityCenterOfMassRemoval( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::~BrookVelocityCenterOfMassRemoval";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _streams[ii];
   }

}

/** 
 * Get inverse of total mass
 * 
 * @return  inverse of total mass
 *
 */

BrookOpenMMFloat BrookVelocityCenterOfMassRemoval::getTotalInverseMass( void ) const {
   return _totalInverseMass;
}

/** 
 * Get inverse mass stream 
 *
 * @return  inverse mass stream
 *
 */

BrookFloatStreamInternal* BrookVelocityCenterOfMassRemoval::getMassStream( void ) const {
   return _streams[MassStream];
}

/** 
 * Get work stream 
 *
 * @return  work stream
 *
 */

BrookFloatStreamInternal* BrookVelocityCenterOfMassRemoval::getWorkStream( void ) const {
   return _streams[WorkStream];
}

/** 
 * Get linear momentum stream
 *
 * @return   linear momentum stream
 *
 */

BrookFloatStreamInternal* BrookVelocityCenterOfMassRemoval::getLinearMomentumStream( void ) const {
   return _streams[LinearMomentumStream];
}

/**
 * Get velocity-COM
 *
 * @param velocities velocities
 * @param  velocityCom                 output velocity com
 *
 * @return DefaultReturnValue
 *
 */
   
int BrookVelocityCenterOfMassRemoval::getVelocityCenterOfMass( BrookStreamImpl& vStream, BrookOpenMMFloat velocityCom[3], BrookOpenMMFloat* ke ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName  = "\nBrookVelocityCenterOfMassRemoval::getVelocityCenterOfMass";

   BrookOpenMMFloat zero    = static_cast<BrookOpenMMFloat>( 0.0 );

   // ---------------------------------------------------------------------------------------

   // calculate linear momentum via reduction
   // subtract it (/totalMass) from velocities

   BrookStreamInternal* velocityStream  = dynamic_cast<BrookStreamInternal*> (vStream.getBrookStreamInternal());
   *ke                                  = zero;

   void* velV                           = velocityStream->getData( 1 );
   const float* vArray                  = static_cast<float*>( velV );

   void* massV                          = getMassStream()->getData( 1 );
   const float* mArray                  = static_cast<float*>( massV );

   int numberOfParticles                = getNumberOfParticles();
   int index                            = 0;

   velocityCom[0] = velocityCom[1] = velocityCom[2] = zero;

   for( int ii = 0; ii < numberOfParticles; ii++, index += 3 ){
      velocityCom[0] += mArray[ii]*vArray[index];
      velocityCom[1] += mArray[ii]*vArray[index+1];
      velocityCom[2] += mArray[ii]*vArray[index+2];
                 *ke += mArray[ii]*( vArray[index]*vArray[index] + vArray[index+1]*vArray[index+1] + vArray[index+2]*vArray[index+2] );
if( 0 && ii < 3 ){
fprintf( stderr, "getVelocityCenterOfMass %5d %5d m=%15.6e v[%15.6e %15.6e %15.6e ] %d\n", ii, index, mArray[ii], vArray[index], vArray[index+1], vArray[index+2], numberOfParticles );
}
   }

   return DefaultReturnValue;
}

/** 
 * Get Particle stream size
 *
 * @return  Particle stream size
 *
 */

int BrookVelocityCenterOfMassRemoval::getComParticleStreamSize( void ) const {
   return _particleStreamSize;
}

/** 
 * Get particle stream width
 *
 * @return  particle stream width
 *
 */

int BrookVelocityCenterOfMassRemoval::getComParticleStreamWidth( void ) const {
   return _particleStreamWidth;
}

/** 
 * Get particle stream height
 *
 * @return particle stream height
 */

int BrookVelocityCenterOfMassRemoval::getComParticleStreamHeight( void ) const {
   return _particleStreamHeight;
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

int BrookVelocityCenterOfMassRemoval::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _particleStreamSize     = getParticleStreamSize( platform );
   _particleStreamWidth    = getParticleStreamWidth( platform );
   _particleStreamHeight   = getParticleStreamHeight( platform );

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

int BrookVelocityCenterOfMassRemoval::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::_initializeStreams";

   BrookOpenMMFloat dangleValue     = static_cast<BrookOpenMMFloat>( 0.0 );

// ---------------------------------------------------------------------------------------

   int particleStreamSize           = getComParticleStreamSize();
   int particleStreamWidth          = getComParticleStreamWidth();

    _streams[WorkStream]            = new BrookFloatStreamInternal( BrookCommon::BrookVelocityCenterOfMassRemovalWorkStream,
                                                                    particleStreamSize, particleStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );


    _streams[MassStream]            = new BrookFloatStreamInternal( BrookCommon::BrookVelocityCenterOfMassRemovalMassStream,
                                                                    particleStreamSize, particleStreamWidth,
                                                                    BrookStreamInternal::Float, dangleValue );


    _streams[LinearMomentumStream]  = new BrookFloatStreamInternal( "LinearMomentumStream",
                                                                    1, 1, BrookStreamInternal::Float3, dangleValue );


   return DefaultReturnValue;
}

/** 
 * Set masses & _totalInverseMass 
 * 
 * @param masses             particle masses
 *
 */

int BrookVelocityCenterOfMassRemoval::_setMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::_setMasses";

   BrookOpenMMFloat zero                    = static_cast<BrookOpenMMFloat>( 0.0 );
   BrookOpenMMFloat one                     = static_cast<BrookOpenMMFloat>( 1.0 );

// ---------------------------------------------------------------------------------------

   // check that masses vector is not larger than expected

   if( (int) masses.size() > getComParticleStreamSize() ){
      std::stringstream message;
      message << methodName << " mass array size=" << masses.size() << " is larger than stream size=" << getComParticleStreamSize();
      throw OpenMMException( message.str() );
   }

   // setup masses

   BrookOpenMMFloat* localMasses = new BrookOpenMMFloat[getComParticleStreamSize()];
   memset( localMasses, 0, sizeof( BrookOpenMMFloat )*getComParticleStreamSize() );

   int index          = 0;
   _totalInverseMass  = zero;
   for( std::vector<double>::const_iterator ii = masses.begin(); ii != masses.end(); ii++, index++ ){
      if( *ii != 0.0 ){
         BrookOpenMMFloat value    = static_cast<BrookOpenMMFloat>(*ii);
         localMasses[index]        = value;
         _totalInverseMass        += value;
      }
   }

   // 1/Sum[masses]

   if( _totalInverseMass > zero ){
      _totalInverseMass = one/_totalInverseMass;
   }

   // write masses to board

   _streams[MassStream]->loadFromArray( localMasses );

   delete[] localMasses;

   return DefaultReturnValue;
}   
 
/*  
 * Setup of StochasticDynamics parameters
 *
 * @param masses                masses
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookVelocityCenterOfMassRemoval::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::setup";

// ---------------------------------------------------------------------------------------

   const BrookPlatform& brookPlatform        = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfParticles                    = (int) masses.size();
   setNumberOfParticles( numberOfParticles );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, platform );
   _initializeStreams( platform );

   _setMasses( masses );

   if( 1 && getLog() ){
      FILE* log             = getLog();
      std::string contents  = getContentsString( 0 );
      (void) fprintf( log, "%s contents:\n%s\n", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }
   
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

std::string BrookVelocityCenterOfMassRemoval::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::getContentsString";

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

   (void) LOCAL_SPRINTF( value, "%d", getComParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getComParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getComParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", getTotalInverseMass() );
   message << _getLine( tab, "TotalInverseMass:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                  ? Set : NotSet) ); 

   message << _getLine( tab, "Mass:",                 (getMassStream()           ? Set : NotSet) ); 
   message << _getLine( tab, "Work:",                 (getWorkStream()           ? Set : NotSet) ); 
   message << _getLine( tab, "LinearMomentum:",       (getLinearMomentumStream() ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _streams[ii] ){
         message << _streams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
/**
 * Remove velocity-COM
 *
 * @param velocities velocities
 *
 * @return DefaultReturnValue
 *
 */
   
int BrookVelocityCenterOfMassRemoval::removeVelocityCenterOfMass( BrookStreamImpl& velocityStream ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName  = "BrookVelocityCenterOfMassRemoval::removeVelocityCenterOfMass";
   static       int printOn             = 0;
   static const double boltz            = 8.3145112119486e-03;
   FILE* log;
   BrookOpenMMFloat ke;

   // ---------------------------------------------------------------------------------------

//setLog( stderr );
   if( printOn && getLog() ){
      log = getLog();
   } else {
      printOn = 0;
   }

   // diagnostics

   if( printOn ){

      BrookOpenMMFloat com[3];
      getVelocityCenterOfMass( velocityStream, com, &ke );

      BrookOpenMMFloat denomiator = static_cast<float>( 3*getNumberOfParticles())*static_cast<float>( boltz );
      (void) fprintf( log, "%s Pre removal com: [%12.5e %12.5e %12.5e] ke=%.3f ~T=%.3f\n", 
                      methodName.c_str(), com[0], com[1], com[2], ke, ke/denomiator );
      BrookStreamInternal* outputVelocityStream = dynamic_cast<BrookStreamInternal*> (velocityStream.getBrookStreamInternal());

      void* velV                       = outputVelocityStream->getData( 1 );
      const float* vArray              = static_cast<float*>( velV );

      int index                        = 0;
      if( printOn > 1 ){
         for( int ii = 0; ii < getNumberOfParticles(); ii++, index += 3 ){
            (void) fprintf( log, "V %d [%12.5e %12.5e %12.5e]\n", ii, vArray[index], vArray[index+1], vArray[index+2] );
         }
      }
      (void) fflush( log );
   }

   // calculate linear momentum via reduction
   // subtract it (/totalMass) from velocities

   kCalculateLinearMomentum( getMassStream()->getBrookStream(), velocityStream.getBrookStream(), getWorkStream()->getBrookStream() );
   kSumLinearMomentum( (float) getComParticleStreamWidth(), (float) getNumberOfParticles(), (float) getTotalInverseMass(),
                        getWorkStream()->getBrookStream(), getLinearMomentumStream()->getBrookStream() );
   kRemoveLinearMomentum( getLinearMomentumStream()->getBrookStream(), velocityStream.getBrookStream(), velocityStream.getBrookStream() );

   // diagnostics

   if( printOn ){

      BrookOpenMMFloat com[3];

      getVelocityCenterOfMass( velocityStream, com, &ke );
      BrookOpenMMFloat denomiator = static_cast<float>( 3*getNumberOfParticles() )*static_cast<float>( boltz );

      (void) fprintf( log, "%s strW=%d iatm=%d invMass=%.4e\n   Post removal com: [%12.5e %12.5e %12.5e]  ke=%.3f ~T=%.3f\n",
                      methodName.c_str(),
                      getComParticleStreamWidth(), getNumberOfParticles(), getTotalInverseMass(),  com[0], com[1], com[2],
                      ke, ke/denomiator );

      void* linMoV = getLinearMomentumStream()->getData( 1 );
      float* linMo = static_cast<float*>( linMoV );
      (void) fprintf( log, "   LM [%12.5e %12.5e %12.5e]\n", linMo[0], linMo[1], linMo[2] );

      BrookStreamInternal* outputVelocityStream = dynamic_cast<BrookStreamInternal*> (velocityStream.getBrookStreamInternal());

      void* velV                      = outputVelocityStream->getData( 1 );
      const float* vArray             = static_cast<float*>( velV );

      void* w1                        = getWorkStream()->getData( 1 );
      const float* w2                 = static_cast<float*>( w1 );

      if( printOn > 1 ){
         int index                        = 0;
         for( int ii = 0; ii < getNumberOfParticles(); ii++, index += 3 ){
            (void) fprintf( log, "V %d [%12.5e %12.5e %12.5e] [%12.5e %12.5e %12.5e]\n", ii,
                            vArray[index], vArray[index+1], vArray[index+2], w2[index], w2[index+1], w2[index+2] );
   
         }
      }
      (void) fflush( log );

//exit(0);
   }

   return DefaultReturnValue;
}

