#ifndef OPENMM_BROOK_LANGEVIN_DYNAMICS_H_
#define OPENMM_BROOK_LANGEVIN_DYNAMICS_H_

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
#include "BrookRandomNumberGenerator.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates stochastic dynamics algorithm 
 *
 */

class BrookLangevinDynamics : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookLangevinDynamics(  );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookLangevinDynamics();
  
      /**
       * Get tau
       *
       * @return tau
       */

      BrookOpenMMFloat getTau( void ) const; 

      /**
       * Get friction
       *
       * @return friction
       */

      BrookOpenMMFloat getFriction( void ) const; 

      /**
       * Get temperature
       *
       * @return temperature
       */

      BrookOpenMMFloat getTemperature( void ) const; 

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
       * Get LangevinDynamics particle stream width
       *
       * @return particle stream width
       */

      int getLangevinDynamicsParticleStreamWidth( void ) const; 

      /**
       * Get LangevinDynamics particle stream height
       *
       * @return particle stream height
       */

      int getLangevinDynamicsParticleStreamHeight( void ) const;

      /**
       * Get LangevinDynamics particle stream size
       * 
       * @return particle stream size
       */

      int getLangevinDynamicsParticleStreamSize( void ) const; 

      /** 
       * Update parameters
       * 
       * @param  temperature     temperature
       * @param  friction        friction
       * @param  step size       step size
       *
       * @return   DefaultReturnValue
       *
       */
      
      int updateParameters( double temperature, double friction, double stepSize );
      
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
      
      int update( BrookStreamImpl& positionStream, BrookStreamImpl& velocityStream,
                  BrookStreamImpl& forceStream,
                  BrookShakeAlgorithm& brookShakeAlgorithm,
                  BrookRandomNumberGenerator& brookRandomNumberGenerator );
      /** 
       * Get array of LangevinDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of LangevinDynamics parameters
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
       * Get SDPC1 stream 
       *
       * @return  SDPC1 stream
       *
       */
      
      BrookFloatStreamInternal* getSDPC1Stream( void ) const;
      
      /** 
       * Get SDPC2 stream 
       *
       * @return  SDPC2 stream
       *
       */
      
      BrookFloatStreamInternal* getSDPC2Stream( void ) const;
      
      /** 
       * Get SD2X stream 
       *
       * @return  SD2X stream
       *
       */
      
      BrookFloatStreamInternal* getSD2XStream( void ) const;
      
      /** 
       * Get SD1V stream 
       *
       * @return  SD1V stream
       *
       */
      
      BrookFloatStreamInternal* getSD1VStream( void ) const;

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

      /** 
       * Get Temperature
       * 
       * @param velocities             velocities
       * @param inverseMassStream      inverse masses
       * @param numberOfConstraints    number of constraints
       *
       * @return temperature
       *
       */
      
      float getTemperature( BrookStreamInternal* velocities, BrookFloatStreamInternal* inverseMassStream, 
                            int numberOfConstraints ) const;

      /** 
       * Remove velocity com
       * 
       * @param velocities             velocities
       * @param inverseMassStream      inverse masses
       *
       * @return DefaultReturnValue  
       */
      
      int removeCom( BrookStreamInternal* velocities, BrookFloatStreamInternal* inverseMassStream ) const;
      
      /** 
       * Reset velocities (diagnostics)
       * 
       * @param velocities             velocities
       *
       * @return DefaultReturnValue  
       */
      
      int resetVelocities( BrookStreamInternal* velocities ) const;
      
   private:
   
      enum DerivedParameters { GDT, EPH, EMH, EP, EM, B, C, D, V, X, Yv, Yx,
                               Sd1pc1, Sd1pc2, Sd1pc3, Sd2pc1, Sd2pc2, MaxDerivedParameters };

      BrookOpenMMFloat _derivedParameters[MaxDerivedParameters];

      // streams indices

      enum BrookLangevinDynamicsStreams { 
              SDPC1Stream,
              SDPC2Stream,
              SD2XStream,
              SD1VStream,
              VPrimeStream,
              XPrimeStream,
              InverseMassStream,
              LastStreamIndex
           };

      // randomNumberSeed

      unsigned int _randomNumberSeed;

      BrookOpenMMFloat _tau;
      BrookOpenMMFloat _temperature;
      BrookOpenMMFloat _stepSize;

      // internal step count

      int _internalStepCount;

      // Particle stream dimensions

      int _sdParticleStreamWidth;
      int _sdParticleStreamHeight;
      int _sdParticleStreamSize;

      /** 
       * Get derived parameter string
       * 
       * @return  string
       *
       */

      //std::string _getDerivedParametersString( BrookLangevinDynamics::DerivedParameters ) const;
      std::string _getDerivedParametersString( int id ) const;

      /** 
       * Update derived parameters
       * 
       * @return  DefaultReturn
       *
       */

      int _updateDerivedParameters( void );

      /*
       * Update streams
       * 
       * @return  DefaultReturn
       *
       */
      
      int _updateSdStreams( void );
      
      // inverse sqrt masses

      BrookOpenMMFloat* _inverseSqrtMasses;

      // internal streams

      BrookFloatStreamInternal* _sdStreams[LastStreamIndex];

      /** 
       * Set tau
       * 
       * @param tau   new tau value
       *
       * @return      DefaultReturn
       *
       */
      
      int _setTau( BrookOpenMMFloat tau );
      
      /** 
       * Set friction = 1/tau
       * 
       * @param friction   new friction value
       *
       * @return      DefaultReturn
       *
       */
      
      int _setFriction( BrookOpenMMFloat friction );
      
      /** 
       * Set temperature
       * 
       * @parameter   temperature
       *
       * @return      DefaultReturn
       *
       */
      
      int _setTemperature( BrookOpenMMFloat temperature );
      
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
       * @param numberOfParticles             number of particles
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
      
      int _setInverseSqrtMasses( const std::vector<double>& masses );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_LANGEVIN_DYNAMICS_H_ */
