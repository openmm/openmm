#ifndef OPENMM_BROOK_COMMON_H_
#define OPENMM_BROOK_COMMON_H_

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

#include <vector>
#include <set>

#include "BrookFloatStreamInternal.h"
#include "BrookIntStreamInternal.h"
#include "BrookPlatform.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class BrookCommon {

   public:
  
      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      // ---------------------------------------------------------------------------------------

      // Stream names
     
      // bonded stream names

      static const std::string BondedParticleIndicesStream;
      static const std::string BondedParametersStream;
      static const std::string UnrolledForceStream;
      static const std::string BondedChargeStream;
      static const std::string BondedInverseMapStreams;

      // nonbonded stream names

      static const std::string NonBondedExclusionStream;
      static const std::string OuterVdwStream;
      static const std::string InnerSigmaStream;
      static const std::string InnerEpsilonStream;
      static const std::string NonBondedChargeStream;
      static const std::string PartialForceStream;

      // OBC Gbsa streams

      static const std::string ObcParticleRadiiStream;
      static const std::string ObcScaledParticleRadiiStream;
      static const std::string ObcParticleRadiiWithDielectricOffsetStream;
      static const std::string ObcBornRadiiStream;
      static const std::string ObcBornRadii2Stream;
      static const std::string ObcIntermediateForceStream;
      static const std::string ObcChainStream;

      // Stochastic Dynamics streams

      static const std::string SDPC1Stream;
      static const std::string SDPC2Stream;
      static const std::string SD2XStream;
      static const std::string SD1VStream;
      static const std::string VPrimeStream;
      static const std::string XPrimeStream;
      static const std::string InverseMassStream;

      // Shake streams

      static const std::string ShakeParticleIndicesStream;
      static const std::string ShakeParticleParameterStream;
      static const std::string ShakeXCons0Stream;
      static const std::string ShakeXCons1Stream;
      static const std::string ShakeXCons2Stream;
      static const std::string ShakeXCons3Stream;
      static const std::string ShakeInverseMapStream;

      // Random number generator streams

      static const std::string ShuffleStream;
      static const std::string RandomValuesStream;

      // BrookVelocityCenterOfMassRemoval streams

      static const std::string BrookVelocityCenterOfMassRemovalWorkStream;
      static const std::string BrookVelocityCenterOfMassRemovalMassStream;

      // ---------------------------------------------------------------------------------------

     /** 
      * Constructor
      * 
      */
     
      BrookCommon( );
  
     /** 
      * Destructor
      * 
      */
     
      ~BrookCommon();
  
      /**
       * Return number of particles
       * 
       * @return number of particles
       *
       */

      int getNumberOfParticles( void ) const; 

      /** 
       * Get particle ceiling parameter
       * 
       * @return particle ceiling parameter
       *
       */
      
      int getParticleSizeCeiling( void ) const;
      
      /**
       * Get particle stream width
       *
       * @param platform platform reference
       *
       * @return particle stream width
       */

      int getParticleStreamWidth( const Platform& platform ); 

      /**
       * Get particle stream width
       *
       * @return particle stream width
       */

      int getParticleStreamWidth( void ) const; 

      /**
       * Get particle stream height
       *
       * @param platform platform reference
       *
       * @return particle stream height
       */

      int getParticleStreamHeight( const Platform& platform ); 

      /**
       * Get particle stream height
       *
       * @return particle stream height
       */

      int getParticleStreamHeight( void ) const;

      /**
       * Get particle stream size
       * 
       * @param platform platform reference
       *
       * @return particle stream size
       */

      int getParticleStreamSize( const Platform& platform ); 

      /**
       * Get particle stream size
       * 
       * @return particle stream size
       */

      int getParticleStreamSize( void ) const; 

      /**
       * Get flag signalling whether active
       * 
       * @return flag signalling whether active
       */

      int isActive( void ) const; 

      /**
       * Set flag signalling whether active
       * 
       * @param isActive  flag signalling whether active
       *
       * @return DefaultReturnValue
       */

      int setIsActive( int isActive ); 

      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;

      /* 
       * Given number of stream elements and width, returns the appropriate
       * height of the stream
       *
       * @param streamSize   stream size 
       * @param width        stream width
       *
       * @return stream height
       *
       */
      
      static int getStreamHeight( int streamSize, int streamWidth );
      
      /* 
       * Given number of stream elements, get stream width & height
       *
       * @param streamSize    stream size 
       * @param streamWidth   output stream width
       * @param streamHeight  output stream height
       *
       * @return stream height
       *
       */
      
      static void getStreamDimensions( int streamSize, int *streamWidth, int *streamHeight );

      /* 
       * Allocate array
       *
       * @param length        length of array
       * @param width         width  of array
       *
       * @return ptr to array
       *
       */
      
      RealOpenMM** allocateRealArray( int length, int width ) const;
      
      /* 
       * Free array
       *
       * @param array         array to be freed (assumed allocated using BrookCommon::allocateRealArray
       *
       * @return DefaultReturnValue
       *
       */
      
      int disposeRealArray( RealOpenMM** array ) const;

      /* 
       * Copy 1D BrookOpenMMFloat* array to 2D array of RealOpenMM
       *
       * @param length        length of array
       * @param width         width  of array
       * @param array1D       array to copy
       *
       * @return ptr to array
       *
       */
      
      RealOpenMM** copy1DArrayTo2DArray( int length, int width, BrookOpenMMFloat* array1D ) const;
      
   protected:
   
      // number of particles

      int _numberOfParticles;

      // particle stream dimensions

      int _particleStreamWidth;
      int _particleStreamHeight;
      int _particleStreamSize;

      // particle size mod

      int _particleSizeModified;

      // log file reference

      FILE* _log;

      // active flag 

      int _isActive;

      /**
       * Set number of particles
       * 
       * @param numberOfParticles number of particles
       *
       */

      int setNumberOfParticles( int numberOfParticles ); 

      /** 
       * Get particle stream dimensions
       * 
       * @param platform                  platform
       *
       */
      
      void _getParticleStreamDimensions( const Platform& platform );
      
      /* 
       * Get line
       *
       * @param tab         tab
       * @param description description
       * @param value       value
       *
       * @return string containing contents
       *
       * */
      
      std::string _getLine( const std::string& tab, const std::string& description, const std::string& value ) const;
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_COMMON_H_ */
