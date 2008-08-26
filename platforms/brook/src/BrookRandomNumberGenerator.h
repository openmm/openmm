#ifndef OPENMM_BROOK_RANDOM_NUMBER_GENERATOR_H_
#define OPENMM_BROOK_RANDOM_NUMBER_GENERATOR_H_

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

class BrookRandomNumberGenerator : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookRandomNumberGenerator(  );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookRandomNumberGenerator();
  
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
       * Get StochasticDynamics atom stream width
       *
       * @return atom stream width
       */

      int getStochasticDynamicsAtomStreamWidth( void ) const; 

      /**
       * Get StochasticDynamics atom stream height
       *
       * @return atom stream height
       */

      int getStochasticDynamicsAtomStreamHeight( void ) const;

      /**
       * Get StochasticDynamics atom stream size
       * 
       * @return atom stream size
       */

      int getStochasticDynamicsAtomStreamSize( void ) const; 

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
       * @param  positions           atom positions
       * @param  velocities          atom velocities
       * @param  forces              atom forces
       * @param  brookShakeAlgorithm BrookShakeAlgorithm reference
       *
       * @return  DefaultReturnValue
       *
       */
      
      int update( Stream& positions, Stream& velocities,
                  const Stream& forces, BrookShakeAlgorithm& brookShakeAlgorithm );
      /** 
       * Get array of StochasticDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of StochasticDynamics parameters
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
       * Get shuffle stream 
       *
       * @return  Shuffle stream
       *
       */
      
      BrookFloatStreamInternal* getShuffleStream( void ) const;
      
      /** 
       * Generate a random number using algorithm in Gromacs
       * 
       * @param ig seed
       *
       * @return  random number
       *
       */
      
      BrookOpenMMFloat generateGromacsRandomNumber( int* ig );
      
      /** 
       * Generate a random number using algorithm in Nvidia code
       * http://www.helsbreth.org/random/rng_kiss.html
       * 
       * @param randomV1   output random value
       * @param randomV2   output random value
       * @param randomV3   output random value
       * @param state      state
       *
       */
      
      void generateRandomsAlaNvidia( float* randomV1, float* randomV2, float* randomV3, 
                                     unsigned int state[4] );

      /** 
       * Load random number streams using Nvidia algorithm
       * 
       *
       * @return DefaultReturnValue;
       */
      
      int loadRandomNumberStreamsNvidia( void );
            
      /** 
       * Get random number seed
       *
       * @return random number seed
       */
      
      unsigned long int getRandomNumberSeed( void ) const;
            
      /** 
       * Increment random number seed
       *
       * @param increment    amount to increment random number seed; default = 1
       *
       * @return updated random number seed
       */
      
      unsigned long int incrementRandomNumberSeed( unsigned long int  increment = 1 );
            
      /** 
       * Set random number seed
       *
       * @param new random number seed; default = 1
       *
       * @return random number seed
       */
      
      unsigned long int setRandomNumberSeed( unsigned long int seed = 1 );
            
   private:
   
      // streams indices

      enum BrookRandomNumberGeneratorStreams { 
              RandomNumberStream,
              ShuffleStream,
              LastStreamIndex
           };

      // randomNumberSeed

      unsigned long int _randomNumberSeed;

      // number of random number streams

      int _numberOfRandomNumberStreams;

      // random number stream dimensions

      int _randomNumberStreamWidth;
      int _randomNumberStreamHeight;
      int _randomNumberStreamSize;

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

#endif /* OPENMM_BROOK_RANDOM_NUMBER_GENERATOR_H_ */
