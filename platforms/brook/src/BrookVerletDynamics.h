#ifndef OPENMM_BROOK_VERLET_DYNAMCIS_H_
#define OPENMM_BROOK_VERLET_DYNAMCIS_H_

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
#include "BrookShakeAlgorithm.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates stochastic dynamics algorithm 
 *
 */

class BrookVerletDynamics : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookVerletDynamics(  );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookVerletDynamics();
  
      /**
       * Get step size
       *
       * @return step size
       */

      BrookOpenMMFloat getStepSize( void ) const; 

      /**
       *
       * Get array of derived parameters indexed by 'DerivedParameters' enums
       *
       * @return array
       *
       */
      
      const BrookOpenMMFloat* getDerivedParameters( void ) const;
      
      /**
       * Get VerletDynamics particle stream width
       *
       * @return particle stream width
       */

      int getVerletDynamicsParticleStreamWidth( void ) const; 

      /**
       * Get VerletDynamics particle stream height
       *
       * @return particle stream height
       */

      int getVerletDynamicsParticleStreamHeight( void ) const;

      /**
       * Get VerletDynamics particle stream size
       * 
       * @return particle stream size
       */

      int getVerletDynamicsParticleStreamSize( void ) const; 

      /** 
       * Update parameters
       * 
       * @param  step size       step size
       *
       * @return   DefaultReturnValue
       *
       */
      
      int updateParameters( double stepSize );
      
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
      
      int update( BrookStreamImpl& positions, BrookStreamImpl& velocities,
                  const BrookStreamImpl& forces, BrookShakeAlgorithm& brookShakeAlgorithm );

      /** 
       * Get array of VerletDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of VerletDynamics parameters
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
       * Get V-prime stream 
       *
       * @return  V-prime stream
       *
       */
      
      BrookFloatStreamInternal* getVPrimeStream( void ) const;

      /** 
       * Get X-prime stream 
       *
       * @return  X-prime stream
       *
       */
      
      BrookFloatStreamInternal* getXPrimeStream( void ) const;

      /** 
       * Get inverse sqrt masses 
       * 
       * @return   inverse sqrt masses stream
       *
       */
      
      BrookFloatStreamInternal* getInverseMassStream( void ) const;

   private:
   
      // streams indices

      enum BrookVerletDynamicsStreams { 
              VPrimeStream,
              XPrimeStream,
              InverseMassStream,
              LastStreamIndex
           };

      // track step (diagnostics)

      int _internalStepCount;

      BrookOpenMMFloat _stepSize;

      // Particle stream dimensions

      int _verletParticleStreamWidth;
      int _verletParticleStreamHeight;
      int _verletParticleStreamSize;

      /*
       * Update streams
       * 
       * @return  DefaultReturn
       *
       */
      
      int _updateVerletStreams( void );
      
      // inverse masses

      BrookOpenMMFloat* _inverseMasses;

      // internal streams

      BrookFloatStreamInternal* _verletStreams[LastStreamIndex];

      /** 
       * Set stepSize
       * 
       * @param   stepSize
       *
       * @return      DefaultReturn
       *
       */
      
      int _setStepSize( BrookOpenMMFloat stepSize );
      
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
       * Set masses 
       * 
       * @param masses             particle masses
       *
       */
      
      int _setInverseMasses( const std::vector<double>& masses );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_VERLET_DYNAMCIS_H_ */
