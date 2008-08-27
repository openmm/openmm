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

#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates stochastic dynamics algorithm 
 *
 */

class BrookRandomNumberGenerator : public BrookCommon {

   public:
  
      // toggle between original rng & Kiss (Nvidia) code

      static const int UseOriginalRng = 1;

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
       * Get number of random number streams
       * 
       * @return     number of random number streams 
       *
       */
      
      int getNumberOfRandomNumberStreams( void ) const;
      
      /**
       * Get stream width
       *
       * @return stream width
       */

      int getRandomNumberStreamWidth( void ) const; 

      /**
       * Get stream height
       *
       * @return stream height
       */

      int getRandomNumberStreamHeight( void ) const;

      /**
       * Get stream size
       * 
       * @return stream size
       */

      int getRandomNumberStreamSize( void ) const; 

      /** 
       * Get array of StochasticDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of RNG parameters
       *
       * @param numberOfAtoms        number of atoms
       * @param platform             Brook platform
       *
       * @return ErrorReturnValue value if error, else DefaultReturnValue
       *
       * */
      
      int setup( int numberOfAtoms,  const Platform& platform  );
      
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
       * Get random number stream 
       *
       * @param index random number stream index     
       *
       * @return  random number stream
       *
       */
      
      BrookFloatStreamInternal* getRandomNumberStream( int index ) const;
      
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
            
      /** 
       * Get index of rv texture
       *
       * @return index of rv texture
       */
      
      int getRvStreamIndex( void ) const;
            
      /** 
       * Get max shuffles
       *
       * @return  max shuffles
       *
       */
      
      int getMaxShuffles( void ) const;
      
      /** 
       * Advance random values stream index
       *
       * @param numberOfEntriesToAdvance number of entries consumed in previous iteration
       *
       * @return  DefaultReturnValue
       *
       */
      
      int advanceGVCursor( int numberOfEntriesToAdvance );
      
      /** 
       * Get random value stream offset
       *
       * @return  random value stream offset
       *
       */
      
      int getRvStreamOffset( void ) const;
      
   private:
   
      // streams indices

      enum BrookRandomNumberGeneratorStreams { 
              ShuffleStream,
              LastStreamIndex
           };

      BrookFloatStreamInternal*  _auxiliaryStreams[LastStreamIndex];
      BrookFloatStreamInternal** _randomNumberGeneratorStreams;

      // randomNumberSeed

      unsigned long int _randomNumberSeed;

      // number of random number streams

      int _numberOfRandomNumberStreams;

      // random number stream dimensions

      int _randomNumberStreamWidth;
      int _randomNumberStreamHeight;
      int _randomNumberStreamSize;

      // control variables

      int _rvStreamIndex;
      int _rvStreamOffset;
      int _numberOfShuffles;
      int _maxShuffles;

      float* _loadBuffer;
      int*   _shuffleIndices;

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
      
      int _initializeStreamSizes(  int numberOfAtoms, const Platform& platform );
      
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
       * Increment random number offset
       *
       * @param increment increment for offset
       *
       * @return random number offset
       */
      
      int _incrementRvOffset( int increment );
            
      /** 
       * Get shuffle stream 
       *
       * @return  Shuffle stream
       *
       */
      
      BrookFloatStreamInternal* _getShuffleStream( void ) const;
      
      /** 
       * Generate a random number using algorithm in Gromacs
       * 
       * @param ig seed
       *
       * @return  random number
       *
       */
      
      BrookOpenMMFloat _generateGromacsRandomNumber( unsigned long int* ig );
      
      /** 
       * Generate a random number using Kiss (algorithm in Kiss code)
       * http://www.helsbreth.org/random/rng_kiss.html
       * 
       * @param randomV1   output random value
       * @param randomV2   output random value
       * @param randomV3   output random value
       * @param state      state
       *
       */
      
      void _generateRandomsKiss( float* randomV1, float* randomV2, float* randomV3, 
                                 unsigned int state[4] );

      /** 
       * Load random number streams using Kiss algorithm
       * 
       *
       * @return DefaultReturnValue;
       */
      
      int _loadRandomNumberStreamsKiss( void );

      /** 
       * Load random number streams using original gpu algorithm
       * 
       *
       * @return DefaultReturnValue;
       */
      
      int _loadGVStreamsOriginal( void );
      
      /** 
       * Loads a permutation of indices from 0 to gvSize-1 in
       * sdp->strShuffle. To make sure that the order of the
       * permutation is atleast NGVSHUFFLE, we create the
       * permutation by introducing a random number of p-cycles
       * where p is randomly determined from 2,3,5,7 and 11.
       * The LCM of these numbers is 2310. 
       * Ofcourse the p-cycles are not necessarily disjoint
       * the way it's done here, but there's a good chance 
       * there will enough disjoint cycles to make the 
       * order of the permutation larger than NGVSHUFFLE
       *
       *
       * This function is only called once at startup
       *
       * @return DefaultReturnValue;
       **/
      
      int _loadGVShuffle( void );

      /** 
       * Get number of shuffles
       *
       * @return  number of shuffles
       *
       */
      
      int _getNumberOfShuffles( void ) const;

      /** 
       * Load buffer
       *
       * @return ptr to load buffer
       *
       * @throw OpenMMException if rv stream size is < 1
       *
       **/
      
      float* _getLoadBuffer( void );
      
      /** 
       * Get ptr to shuffle indices
       *
       * @return ptr to shuffle indices
       *
       * @throw OpenMMException if size is < 1
       *
       **/
      
      int* _getShuffleIndices( int size );
      
      /** 
       * Shuffle streams
       *
       * @return DefaultReturnValue;
       */
      
      int _shuffleGVStreams( void );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_RANDOM_NUMBER_GENERATOR_H_ */
