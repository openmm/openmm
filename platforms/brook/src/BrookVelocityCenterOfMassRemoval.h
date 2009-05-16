#ifndef OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_
#define OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_

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

#include "BrookStreamImpl.h"
//#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates removal of center of mass
 *
 */

class BrookVelocityCenterOfMassRemoval : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookVelocityCenterOfMassRemoval( );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookVelocityCenterOfMassRemoval( );
  
      /**
       * Get  particle stream width
       *
       * @return particle stream width
       */

      int getComParticleStreamWidth( void ) const; 

      /**
       * Get  particle stream height
       *
       * @return particle stream height
       */

      int getComParticleStreamHeight( void ) const;

      /**
       * Get  particle stream size
       * 
       * @return particle stream size
       */

      int getComParticleStreamSize( void ) const; 

      /** 
       * Remove velocity center-of-mass
       * 
       * @param  velocities                  particle velocities
       *
       * @return  DefaultReturnValue
       *
       */
      
      int removeVelocityCenterOfMass( BrookStreamImpl& velocityStream );

      /** 
       * Get velocity center-of-mass and kinetic energy (used for diagnostics)
       * 
       * @param  velocities                  particle velocities
       * @param  velocityCom                 output velocity com
       * @param  ke                          output kinetic energy
       *
       * @return  DefaultReturnValue
       *
       */
      
      int getVelocityCenterOfMass( BrookStreamImpl& vStream, BrookOpenMMFloat velocityCom[3], BrookOpenMMFloat* ke );

      /* 
       * Setup of  parameters
       *
       * @param masses                particle masses
       * @param platform              Brook platform
       *
       * @return ErrorReturnValue value if error, else DefaultReturnValue
       *
       * */
      
      int setup( const std::vector<double>& masses, const Platform& platform  );
      
      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContentsString( int level = 0 ) const;

      /**
       * Get inverse of total mass
       *
       * @return inverse of total mass
       *
       */
         
      BrookOpenMMFloat getTotalInverseMass( void ) const;
      
      /** 
       * Get inverse mass stream 
       *
       * @return  inverse mass stream
       *
       */
      
      BrookFloatStreamInternal* getMassStream( void ) const;
      
      /** 
       * Get work stream 
       *
       * @return  work stream
       *
       */
      
      BrookFloatStreamInternal* getWorkStream( void ) const;
      
      /** 
       * Get linear momentum stream 
       *
       * @return  linear momentum stream
       *
       */
      
      BrookFloatStreamInternal* getLinearMomentumStream( void ) const;
      
   private:
   
      // streams indices

      enum BrookVelocityCenterOfMassRemovalStreams { 
              MassStream,
              WorkStream,
              LinearMomentumStream,
              LastStreamIndex
           };

      BrookOpenMMFloat _totalInverseMass;

      // Particle stream dimensions

      int _particleStreamWidth;
      int _particleStreamHeight;
      int _particleStreamSize;

      // internal streams

      BrookFloatStreamInternal* _streams[LastStreamIndex];

      /* 
       * Setup of stream dimensions
       *
       * @param particleStreamSize        particle stream size
       * @param particleStreamWidth       particle stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int _initializeStreamSizes( int particleStreamSize, int particleStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfParticles         number of particles
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       */
      
      int _initializeStreamSizes( int numberOfParticles, const Platform& platform );
      
      /** 
       * Initialize stream dimensions and streams
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int _initializeStreams( const Platform& platform );

      /** 
       * Set inverse masses 
       * 
       * @param masses             particle masses
       *
       */
      
      int _setMasses( const std::vector<double>& masses );
      
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_ */
