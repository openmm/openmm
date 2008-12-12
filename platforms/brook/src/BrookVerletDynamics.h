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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include <vector>
#include <set>

#include "BrookFloatStreamInternal.h"
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
      
      int update( Stream& positions, Stream& velocities,
                  const Stream& forces, BrookShakeAlgorithm& brookShakeAlgorithm );

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
