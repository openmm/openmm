#ifndef OPENMM_BROOK_BROWNIAN_DYNAMCIS_H_
#define OPENMM_BROOK_BROWNIAN_DYNAMCIS_H_

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
#include "BrookRandomNumberGenerator.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates stochastic dynamics algorithm 
 *
 */

class BrookBrownianDynamics : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookBrownianDynamics(  );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookBrownianDynamics();
  
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
       * Get noise amplitude
       *
       * @return noise amplitude
       */

      BrookOpenMMFloat getNoiseAmplitude( void ) const; 

      /**
       * Get force scale
       *
       * @return force scale
       */

      BrookOpenMMFloat getForceScale( void ) const; 

      /**
       * Get BrownianDynamics atom stream width
       *
       * @return atom stream width
       */

      int getBrownianDynamicsAtomStreamWidth( void ) const; 

      /**
       * Get BrownianDynamics atom stream height
       *
       * @return atom stream height
       */

      int getBrownianDynamicsAtomStreamHeight( void ) const;

      /**
       * Get BrownianDynamics atom stream size
       * 
       * @return atom stream size
       */

      int getBrownianDynamicsAtomStreamSize( void ) const; 

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
       * @param  positions                   atom positions
       * @param  velocities                  atom velocities
       * @param  forces                      atom forces
       * @param  brookShakeAlgorithm         BrookShakeAlgorithm reference
       * @param  brookRandomNumberGenerator  BrookRandomNumberGenerator reference
       *
       * @return  DefaultReturnValue
       *
       */
      
      int update( Stream& positions, Stream& velocities,
                  const Stream& forces, BrookShakeAlgorithm& brookShakeAlgorithm,
                  BrookRandomNumberGenerator& brookRandomNumberGenerator );

      /** 
       * Get array of BrownianDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of BrownianDynamics parameters
       *
       * @param masses                atom masses
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

      enum BrookBrownianDynamicsStreams { 
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
      BrookOpenMMFloat _forceScale;
      BrookOpenMMFloat _noiseAmplitude;

      // Atom stream dimensions

      int _bdAtomStreamWidth;
      int _bdAtomStreamHeight;
      int _bdAtomStreamSize;

      /*
       * Update streams
       * 
       * @return  DefaultReturn
       *
       */
      
      int _updateStreams( void );
      
      // inverse sqrt masses

      BrookOpenMMFloat* _inverseSqrtMasses;

      // internal streams

      BrookFloatStreamInternal* _streams[LastStreamIndex];

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
      
      /** 
       * Set force scale
       * 
       * @return DefaultReturn
       *
       */
      
      int _setForceScale( void );
      
      /** 
       * Set noise amplitude
       * 
       * @return DefaultReturn
       *
       */
      
      int _setNoiseAmplitude( void );
      
      /* 
       * Setup of stream dimensions
       *
       * @param atomStreamSize        atom stream size
       * @param atomStreamWidth       atom stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int _initializeStreamSizes( int atomStreamSize, int atomStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfAtoms             number of atoms
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       */
      
      int _initializeStreamSizes( int numberOfAtoms, const Platform& platform );
      
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
       * @param masses             atomic masses
       *
       */
      
      int _setInverseSqrtMasses( const std::vector<double>& masses );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_BROWNIAN_DYNAMCIS_H_ */
