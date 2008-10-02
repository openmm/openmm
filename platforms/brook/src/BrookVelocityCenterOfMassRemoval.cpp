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
#include "gpu/kcom.h"

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

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms           = -1;

   // mark stream dimension variables as unset

   _atomStreamWidth         = -1;
   _atomStreamHeight        = -1;
   _atomStreamSize          = -1;

   _totalInverseMass        = zero;

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
 * Remove velocity-COM
 *
 * @param velocities velocities
 *
 * @return DefaultReturnValue
 *
 */
   
int BrookVelocityCenterOfMassRemoval::removeVelocityCenterOfMass( Stream& velocities ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "BrookVelocityCenterOfMassRemoval::removeVelocityCenterOfMass";

   static int debug               = 0;
   float zero                     = 0.0f;

   // ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      BrookOpenMMFloat com[3];
      getVelocityCenterOfMass( velocities, com );
      (void) fprintf( getLog(), "\n%s Pre removal com: [%12.5e %12.5e %12.5e]\n", methodName, com[0], com[1], com[2] );
   }

   // calculate linear momentum via reduction
   // subtract it (/totalMass) from velocities

   BrookStreamImpl& velocityStream  = dynamic_cast<BrookStreamImpl&> (velocities.getImpl());

   kCalculateLinearMomentum( getMassStream()->getBrookStream(), velocityStream.getBrookStream(), getWorkStream()->getBrookStream() );
   kSumLinearMomentum( (float) getComAtomStreamWidth(), (float) getNumberOfAtoms(), getWorkStream()->getBrookStream(), getLinearMomentumStream()->getBrookStream() );
   kScale( (float) getTotalInverseMass(), getLinearMomentumStream()->getBrookStream(), getLinearMomentumStream()->getBrookStream() );
   kRemoveLinearMomentum( getLinearMomentumStream()->getBrookStream(), velocityStream.getBrookStream(), velocityStream.getBrookStream() );

   if( (0 || debug) && getLog() ){
      BrookOpenMMFloat com[3];
      getVelocityCenterOfMass( velocities, com );
      (void) fprintf( getLog(), "\n%s strW=%d iatm=%d Post removal com: [%12.5e %12.5e %12.5e]\n", methodName,
                      getComAtomStreamWidth(), getNumberOfAtoms(),  com[0], com[1], com[2] );
/*
      (void) fprintf( getLog(), "LM [%12.5e %12.5e %12.5e]\n", com[0], com[1], com[2] );
      BrookStreamImpl& v1                  = dynamic_cast<BrookStreamImpl&> (velocities.getImpl());
      BrookStreamInternal* velocityStream  = dynamic_cast<BrookStreamInternal*> (v1.getBrookStreamImpl());

      void* velV                       = velocityStream->getData( 1 );
      const float* vArray              = (float*) velV;

      void* w1                        = getWorkStream()->getData( 1 );
      const float* w2                 = (float*) w1;

      int index                        = 0;
      for( int ii = 0; ii < getNumberOfAtoms(); ii++, index += 3 ){
         (void) fprintf( getLog(), "V %d [%12.5e %12.5e %12.5e] [%12.5e %12.5e %12.5e]\n", ii,
                         vArray[index], vArray[index+1], vArray[index+2], w2[index], w2[index+1], w2[index+2] );
      }
      (void) fflush( getLog() );
*/
//exit(0);
   }

   return DefaultReturnValue;
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
   
int BrookVelocityCenterOfMassRemoval::getVelocityCenterOfMass( Stream& velocities, BrookOpenMMFloat velocityCom[3] ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookVelocityCenterOfMassRemoval::getVelocityCenterOfMass";

   static int debug         = 0;
   BrookOpenMMFloat zero    = (BrookOpenMMFloat) 0.0;

   // ---------------------------------------------------------------------------------------

   // calculate linear momentum via reduction
   // subtract it (/totalMass) from velocities

   BrookStreamImpl& v1                  = dynamic_cast<BrookStreamImpl&> (velocities.getImpl());
   BrookStreamInternal* velocityStream  = dynamic_cast<BrookStreamInternal*> (v1.getBrookStreamImpl());

   void* velV                           = velocityStream->getData( 1 );
   const float* vArray                  = (float*) velV;

   void* massV                          = getMassStream()->getData( 1);
   const float* mArray                  = (float*) massV;

   int numberOfAtoms                    = getNumberOfAtoms();
   int index                            = 0;

   velocityCom[0] = velocityCom[1] = velocityCom[2] = zero;

   for( int ii = 0; ii < numberOfAtoms; ii++, index += 3 ){
      velocityCom[0] += mArray[ii]*vArray[index];
      velocityCom[1] += mArray[ii]*vArray[index+1];
      velocityCom[2] += mArray[ii]*vArray[index+2];
   }

   return DefaultReturnValue;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookVelocityCenterOfMassRemoval::getComAtomStreamSize( void ) const {
   return _atomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookVelocityCenterOfMassRemoval::getComAtomStreamWidth( void ) const {
   return _atomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookVelocityCenterOfMassRemoval::getComAtomStreamHeight( void ) const {
   return _atomStreamHeight;
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

int BrookVelocityCenterOfMassRemoval::_initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _atomStreamSize     = getAtomStreamSize( platform );
   _atomStreamWidth    = getAtomStreamWidth( platform );
   _atomStreamHeight   = getAtomStreamHeight( platform );

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

   BrookOpenMMFloat dangleValue            = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int atomStreamSize               = getComAtomStreamSize();
   int atomStreamWidth              = getComAtomStreamWidth();

    _streams[WorkStream]            = new BrookFloatStreamInternal( BrookCommon::BrookVelocityCenterOfMassRemovalWorkStream,
                                                                    atomStreamSize, atomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );


    _streams[MassStream]            = new BrookFloatStreamInternal( BrookCommon::BrookVelocityCenterOfMassRemovalMassStream,
                                                                    atomStreamSize, atomStreamWidth,
                                                                    BrookStreamInternal::Float, dangleValue );


    _streams[LinearMomentumStream]  = new BrookFloatStreamInternal( "LinearMomentumStream",
                                                                    1, 1, BrookStreamInternal::Float3, dangleValue );


   return DefaultReturnValue;
}

/** 
 * Set masses & _totalInverseMass 
 * 
 * @param masses             atomic masses
 *
 */

int BrookVelocityCenterOfMassRemoval::_setMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookVelocityCenterOfMassRemoval::_setMasses";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   // check that masses vector is not larger than expected

   if( (int) masses.size() > getComAtomStreamSize() ){
      std::stringstream message;
      message << methodName << " mass array size=" << masses.size() << " is larger than stream size=" << getComAtomStreamSize();
      throw OpenMMException( message.str() );
   }

   // setup masses

   BrookOpenMMFloat* localMasses = new BrookOpenMMFloat[getComAtomStreamSize()];
   memset( localMasses, 0, sizeof( BrookOpenMMFloat )*getComAtomStreamSize() );

   int index          = 0;
   _totalInverseMass  = (BrookOpenMMFloat) 0.0;
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

   const BrookPlatform brookPlatform        = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   int numberOfAtoms  = (int) masses.size();
   setNumberOfAtoms( numberOfAtoms );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfAtoms, platform );
   _initializeStreams( platform );

   _setMasses( masses );

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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getComAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getComAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getComAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTotalInverseMass() );
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
