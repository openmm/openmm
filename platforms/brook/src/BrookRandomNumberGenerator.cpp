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
#include "BrookRandomNumberGenerator.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMUtilities.h"
#include "OpenMMException.h"
#include "kernels/kupdatesd.h"

using namespace OpenMM;
using namespace std;

/** 
 *
 * Constructor
 * 
 */

BrookRandomNumberGenerator::BrookRandomNumberGenerator( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::BrookRandomNumberGenerator";

// ---------------------------------------------------------------------------------------

   // fixed for now

   _numberOfRandomNumberStreams     = 2;
   _randomNumberGeneratorStreams    = NULL;

   // mark stream dimension variables as unset

   _randomNumberStreamWidth         = -1;
   _randomNumberStreamHeight        = -1;
   _randomNumberStreamSize          = -1;

   _rvStreamIndex                   = 0;
   _rvStreamOffset                  = 0;
   _numberOfShuffles                = 0;
   //_maxShuffles                     = 0;
   _maxShuffles                     = 100;

   _loadBuffer                      = NULL;
   _shuffleIndices                  = NULL;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _auxiliaryStreams[ii]   = NULL;
   }

   // set randomNumber seed & generator

   _randomNumberSeed      = 1393;
   //_randomNumberGenerator = Mersenne;
   _randomNumberGenerator = Kiss;
   //_randomNumberGenerator = FixedValue;

   //SimTKOpenMMUtilities::setRandomNumberSeed( randomNumberSeed );
}   
 
/** 
 * Destructor
 * 
 */

BrookRandomNumberGenerator::~BrookRandomNumberGenerator( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::~BrookRandomNumberGenerator";

// ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      delete _auxiliaryStreams[ii];
   }

   delete[] _randomNumberGeneratorStreams;

   delete[] _loadBuffer;
   delete[] _shuffleIndices;
}

/** 
 * Get number of random number streams
 * 
 * @return     number of random number streams 
 *
 */

int BrookRandomNumberGenerator::getNumberOfRandomNumberStreams( void ) const {
   return _numberOfRandomNumberStreams;
}

/** 
 * Get random number stream width
 * 
 * @return     nndom number stream width
 *
 */

int BrookRandomNumberGenerator::getRandomNumberStreamWidth( void ) const {
   return _randomNumberStreamWidth;
}

/** 
 * Get random number stream height
 * 
 * @return     nndom number stream height
 *
 */

int BrookRandomNumberGenerator::getRandomNumberStreamHeight( void ) const {
   return _randomNumberStreamHeight;
}

/** 
 * Get random number seed
 *
 * @return random number seed
 */
    
unsigned long int BrookRandomNumberGenerator::getRandomNumberSeed( void ) const {
   return _randomNumberSeed;
}
          
/** 
 * Increment random number seed
 *
 * @param increment    amount to increment random number seed; default = 1
 *
 * @return updated random number seed
 */
     
unsigned long int BrookRandomNumberGenerator::incrementRandomNumberSeed( unsigned long int  increment ){
   _randomNumberSeed += increment;
   return _randomNumberSeed;
}

/** 
 * Set random number seed
 *
 * @param new random number seed; default = 1
 *
 * @return random number seed
 */
    
unsigned long int BrookRandomNumberGenerator::setRandomNumberSeed( unsigned long int seed ){
   _randomNumberSeed  = seed;
   return _randomNumberSeed;
}

/** 
 * Generate a random number using algorithm in Gromacs
 * 
 * @param ig seed
 *
 * @return  random number
 *
 */

BrookOpenMMFloat BrookRandomNumberGenerator::_generateGromacsRandomNumber( unsigned long int* ig ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nBrookRandomNumberGenerator::_generateGromacsRandomNumber";

   int  irand;
 
   int  m     = 100000000;
   float rm   = 100000000.0;  /* same number as m, but real format */
   int  m1    = 10000;
   int  mult  = 31415821;
   
   BrookOpenMMFloat r;
   int  irandh,irandl,multh,multl;
 
   // ---------------------------------------------------------------------------------------
 
   unsigned long int igg = (*ig > 0) ? *ig : -1*(*ig);
   irand = igg % m; 
   
   /* multiply irand by mult, but take into account that overflow
    * must be discarded, and do not generate an error.
    */

   irandh = irand / m1;
   irandl = irand % m1;
   multh  = mult / m1;
   multl  = mult % m1;
   irand  = ((irandh*multl+irandl*multh) % m1) * m1 + irandl*multl;
   irand  = (irand + 1) % m; 
 
   /* convert irand to a real random number between 0 and 1. */

   r = (BrookOpenMMFloat) (irand / 10); 
   r = r * 10 / rm;
   if ((r <= 0) || (r > 1))
     r = 0.0; 
   *ig = irand;
   
   return r;
}     

inline int MaxInt( unsigned int x, unsigned int y ){ return x > y ? x : y; }

/** 
 * Generate a random number using algorithm in Kiss code
 * http://www.helsbreth.org/random/rng_kiss.html
 * 
 * @param randomV1   output random value
 * @param randomV2   output random value
 * @param randomV3   output random value
 * @param state      state
 *
 */

void BrookRandomNumberGenerator::_generateRandomsKiss( float* randomV1, float* randomV2, float* randomV3, 
                                                       unsigned int state[4] ){
    
   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nBrookRandomNumberGenerator::_generateRandomsKiss";

   unsigned int carry          = 0;

   // ---------------------------------------------------------------------------------------
 
   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   unsigned int k       = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   unsigned int m       = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;
   unsigned int z1      = MaxInt(state[0] + state[1] + state[3], 0x00000001);
   float x1             = (float) ( (double) z1 / (double)UINT_MAX );

/*
if( x1 < 0.0f ){
   unsigned int z1 = MaxInt(state[0] + state[1] + state[3], 0x00000001);
   double       z2 = (double) z1/((double) UINT_MAX);
   (void) fprintf( logFile, "x1=%.6e state[%u %u %u] sum %u %.6e den=%u z2=%.6e\n",
                   x1, state[0], state[1], state[3], 
                   (state[0] + state[1] + state[3]), (float) (state[0] + state[1] + state[3]), z1, z2 );
}
*/

   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   x1                   = sqrt(-2.0f * log(x1));
   k                    = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   m                    = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;
   float x2             = (float)(state[0] + state[1] + state[3]) / (float)UINT_MAX;
   
   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   *randomV1            = x1 * cos(2.0f * 3.14159265f * x2);

   k                    = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   m                    = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;

   unsigned int z3      = MaxInt(state[0] + state[1] + state[3], 0x00000001);
   float x3             = (float) ( (double) z3 / (double)UINT_MAX );
   //float x3             = (float)MaxInt(state[0] + state[1] + state[3], 0x00000001) / (float)UINT_MAX;

   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   x3                   = sqrt(-2.0f * log(x3));
   k                    = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   m                    = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;
   float x4             = (float)(state[0] + state[1] + state[3]) / (float)UINT_MAX;
   
   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   *randomV2            = x3 * cos(2.0f * 3.14159265f * x4);
   k                    = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   m                    = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;

   //float x5             = (float)MaxInt(state[0] + state[1] + state[3], 0x00000001) / (float)UINT_MAX;
   unsigned int z5      = MaxInt(state[0] + state[1] + state[3], 0x00000001);
   float x5             = (float) ( (double) z5 / (double)UINT_MAX );

   state[0]             = state[0] * 69069 + 1;
   state[1]            ^= state[1] << 13;
   state[1]            ^= state[1] >> 17;
   state[1]            ^= state[1] << 5;
   x5                   = sqrt(-2.0f * log(x5));
   k                    = (state[2] >> 2) + (state[3] >> 3) + (carry >> 2);
   m                    = state[3] + state[3] + state[2] + carry;
   state[2]             = state[3];
   state[3]             = m;
   carry                = k >> 30;
   float x6             = (float)(state[0] + state[1] + state[3]) / (float)UINT_MAX;
   *randomV3            = x5 * cos(2.0f * 3.14159265f * x6); 
   
//(void) fprintf( logFile, "rv=%.6e %.6e %.6e\n", *randomV1, *randomV2, *randomV3 );
//exit(0);

   return;
}

/** 
 * Load random number streams using Kiss algorithm
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::_loadRandomNumberStreamsKiss( void ){

   // ---------------------------------------------------------------------------------------

   static unsigned int state[4];
   static int stateInitialized    = 0;
   int printOn                    = 0;
   FILE* log;
   static const int reseed        = 10000;
 
   static std::string methodName  = "\nBrookRandomNumberGenerator::_loadRandomNumberStreamsKiss";

   // ---------------------------------------------------------------------------------------
   
   if( printOn && getLog() ){
       log = getLog();
   } else {
      printOn = 0;
   }   

   // periodically reset seeds

   if( !stateInitialized || !(stateInitialized % reseed) ){

      state[0] = rand();
      state[1] = rand();
      state[2] = rand();
      state[3] = rand();

      if( printOn ){
         (void) fprintf( log, "%s reset state seeds stateInitialized=%d reseed=%d [%u %u %u %u]\n",
                         methodName.c_str(), stateInitialized, reseed,  state[0], state[1], state[2], state[3] );
         (void) fflush( log );
      }
   }
   stateInitialized++;

   // allocate memory once for download of random nos.

   float* loadBuffer = _getLoadBuffer();

   if( printOn ){
	   static float count   = 0.0f;
      float block          = (float) (3*getRandomNumberStreamSize() );
      count               += 1.0f;
      (void) fprintf( log, "%s: count=%.1f ttl=%.3e no./count=%.1f %d %d\n", methodName.c_str(),
                      count, block*count, block, getRandomNumberStreamSize(), getNumberOfRandomNumberStreams() );
      (void) fflush( log );
   }
	
   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
      for( int ii = 0; ii < 3*getRandomNumberStreamSize(); ii += 3 ){
         float v1,v2,v3;
         _generateRandomsKiss( &v1, &v2, &v3, state );
         loadBuffer[ii]   = v1;
         loadBuffer[ii+1] = v2;
         loadBuffer[ii+2] = v3;
   	}
	   getRandomNumberStream( jj )->loadFromArray( loadBuffer );
   }

   if( printOn ){
      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s: stats\n%s\n", methodName.c_str(), getStatisticsString().c_str() );
      (void) fflush( log );
   }

   return DefaultReturnValue;
}

/** 
 * Load random number streams using Mersenne algorithm
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::_loadRandomNumberStreamsMersenne( void ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName  = "\nBrookRandomNumberGenerator::_loadRandomNumberStreamsMersenne";
   int printOn                          = 0;

   // ---------------------------------------------------------------------------------------
   
   // allocate memory once for download of random nos.

   float* loadBuffer = _getLoadBuffer();

   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
      for( int ii = 0; ii < 3*getRandomNumberStreamSize(); ii += 3 ){
         loadBuffer[ii]   = (float) SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
         loadBuffer[ii+1] = (float) SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
         loadBuffer[ii+2] = (float) SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
   	}
	   getRandomNumberStream( jj )->loadFromArray( loadBuffer );
   }

   if( printOn ){
      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s: stats\n%s\n", methodName.c_str(), getStatisticsString().c_str() );
      (void) fflush( log );
   }

   return DefaultReturnValue;
}

/** 
 * Load random number streams w/ fixed value
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::_loadRandomNumberStreamsFixedValue( void ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName  = "\nBrookRandomNumberGenerator::_loadRandomNumberStreamsFixedValue";
   int printOn                          = 1;

   // ---------------------------------------------------------------------------------------
   
   // load fixed value

   float fixedValue  = 0.1f;
   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
	   getRandomNumberStream( jj )->fillWithValue( &fixedValue );
   }

   if( printOn ){
      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s: stats\n%s\n", methodName.c_str(), getStatisticsString().c_str() );
      (void) fflush( log );
   }

   return DefaultReturnValue;
}

/** 
 * Load random number streams using original gpu algorithm
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::_loadGVStreamsOriginal( void ){

   // ---------------------------------------------------------------------------------------

	unsigned long int jran;

   // static const std::string methodName = "\nBrookRandomNumberGenerator::_loadGVStreamsOriginal";

   // ---------------------------------------------------------------------------------------
	
   float* loadBuffer = _getLoadBuffer();
	jran              = getRandomNumberSeed();
	
   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
      for( int ii = 0; ii < 3*getRandomNumberStreamSize(); ii++ ){
			loadBuffer[ii] = _generateGromacsRandomNumber( &jran );
		}
	   getRandomNumberStream( jj )->loadFromArray( loadBuffer );
	}
	
   incrementRandomNumberSeed( 1 );

   return DefaultReturnValue;
}

/** 
 * Get load buffer
 *
 * @return ptr to load buffer
 *
 * @throw OpenMMException if rv stream size is < 1
 *
 **/

float* BrookRandomNumberGenerator::_getLoadBuffer( void ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nBrookRandomNumberGenerator::_getLoadBuffer";

   // ---------------------------------------------------------------------------------------
	
   if( getRandomNumberStreamSize() < 1 ){
      std::stringstream message;
      message << methodName << " rv stream size=" << getRandomNumberStreamSize() << " is less than 1.";
      throw OpenMMException( message.str() );
      return NULL;
   }

   if( _loadBuffer == NULL ){
   	_loadBuffer = (float*) malloc( sizeof(float)*3*getRandomNumberStreamSize() );
   }

   return _loadBuffer;

}	

/** 
 * Get ptr to shuffle indices
 *
 * @return ptr to shuffle indices
 *
 * @throw OpenMMException if size is < 1
 *
 **/

int* BrookRandomNumberGenerator::_getShuffleIndices( int size ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nBrookRandomNumberGenerator::_getShuffleIndices";

   // ---------------------------------------------------------------------------------------
	
   if( size < 1 ){
      std::stringstream message;
      message << methodName << " size=" << size << " is less than 1.";
      throw OpenMMException( message.str() );
      return NULL;
   }

   if( _shuffleIndices == NULL ){
   	_shuffleIndices = (int*) malloc( sizeof(int)*size );
   }

   return _shuffleIndices;

}	

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

int BrookRandomNumberGenerator::_loadGVShuffle( void ){

   // ---------------------------------------------------------------------------------------

	const int p[] = { 2, 3, 5, 7, 11 };
	const int np  = sizeof(p) / sizeof(p[0]);
   const int pmax = p[np-1];

   static const std::string methodName  = "\nBrookRandomNumberGenerator::loadGVShuffle";
   int printOn                          = 0;

   // ---------------------------------------------------------------------------------------
	
	float* loadBuffer  = _getLoadBuffer();
   int* indices       = _getShuffleIndices( pmax );
	
   int rvSize = getRandomNumberStreamSize();
	for ( int ii = 0; ii < rvSize; ii++ ) {
		loadBuffer[ii] = (float) ii;
	}

	//How to come up with this number here?

   unsigned long int seed = getRandomNumberSeed();
	for( int iter = 0; iter < 1000000; iter++ ) {

		//for each p

		for( int ii = 0; ii < np; ii++ ){

			//Come up with p random indices
			//Not checking that they are distinct
			//because that should be fairly rare

			for ( int jj = 0; jj < p[ii]; jj++ ) {
				//indices[j] = (int) ( gmx_rando( &sdp->seed ) * sdp->gvSize );
				indices[jj] = (int) ( _generateGromacsRandomNumber( &seed )*rvSize );
			}

			// do a p-cycle

			float tmp = loadBuffer[ indices[0] ];
			for ( int jj = 0; jj < p[ii]-1; jj++ ) {
				loadBuffer[ indices[jj] ] = loadBuffer[ indices[jj+1] ];
			}
			loadBuffer[ indices[ p[ii]-1 ] ] = tmp;
		}
	}
   _getShuffleStream()->loadFromArray( loadBuffer );
	
   if( printOn ){
      FILE* log = getLog() ? getLog() : stderr;
      (void) fprintf( log, "%s: sz=%d pmax=%d np=%d sample indices:\n", methodName.c_str(), rvSize, pmax, np );
      float maxIdx = -1.0f;
      float minIdx =  (float) (rvSize)*2.0f;
      for( int ii = 0; ii < rvSize; ii++ ){
         if( ii < 30 || ii > (rvSize-30) ){
            (void) fprintf( log, "  %d %.8f\n", ii, loadBuffer[ii] );
         }
         if( loadBuffer[ii] < minIdx ){
            minIdx = loadBuffer[ii];
         }
         if( loadBuffer[ii] > maxIdx ){
            maxIdx = loadBuffer[ii];
         }
      }
      (void) fprintf( log, "%s: min-max indices: %.8f %.8f\n", methodName.c_str(), minIdx, maxIdx );
   }
      
   return DefaultReturnValue;
}

/** 
 * Shuffle streams
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::_shuffleGVStreams( void ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nBrookRandomNumberGenerator::_shuffleGVStreams";

   // ---------------------------------------------------------------------------------------
	
   int numberOfRvStreams = getNumberOfRandomNumberStreams();

	for( int ii = 0; ii < (numberOfRvStreams - 1); ii++ ){
		kpermute_vectors( (float) getRandomNumberStreamWidth(),
                         _getShuffleStream()->getBrookStream(),
                         getRandomNumberStream( ii + 1 )->getBrookStream(), 
                         getRandomNumberStream( ii  )->getBrookStream() );
   }

	kpermute_vectors( (float) getRandomNumberStreamWidth(),
                      _getShuffleStream()->getBrookStream(),
                      getRandomNumberStream( 0 )->getBrookStream(), 
                      getRandomNumberStream( numberOfRvStreams - 1 )->getBrookStream() );

   return DefaultReturnValue;
}

/** 
 * Advances the current position by 2*gpu->particles
 * If we run out of rand numbers, we shuffle and 
 * reuse a few times before reloading from the cpu
 *
 * @param numberOfRandomValuesConsumedPerIteration    number of random values consumed/iteration
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::advanceGVCursor( int numberOfRandomValuesConsumedPerIteration ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookRandomNumberGenerator::advanceGVCursor";
   int printOn                         = 0;
   FILE* log;

   // ---------------------------------------------------------------------------------------
	
//setLog( stderr );
   if( printOn && getLog() ){
      log = getLog();
   } else {
      printOn = 0;
   }   

   int rvStreamSize  = getRandomNumberStreamSize();

   // use 2 random values per sd 

	_rvStreamOffset  += numberOfRandomValuesConsumedPerIteration;

	//Check if we've used up this texture

   char* action = "none";
	if ( _rvStreamOffset > rvStreamSize - numberOfRandomValuesConsumedPerIteration ){

		// next one if available

      _rvStreamOffset = 0;

		if ( _rvStreamIndex < _numberOfRandomNumberStreams - 1 ){
         _rvStreamIndex++;
         action = "incremented stream index";
		} else {

			//No more textures, need to shuffle

			if( _numberOfShuffles < _maxShuffles ){

				_shuffleGVStreams( );
            _numberOfShuffles++;
            action = "shuffled streams";

			} else { //Need to refresh random numbers from cpu

            if( _randomNumberGenerator == Mersenne ){
               action = "loaded new values to GPU using Mersenne rng";
               _loadRandomNumberStreamsMersenne( );
            } else if( _randomNumberGenerator == Kiss ){
               action = "loaded new values to GPU using KISS rng";
               _loadRandomNumberStreamsKiss( );
            } else if( _randomNumberGenerator == Original ){
               action = "loaded new values to GPU using original Gromac's rng";
               _loadGVStreamsOriginal( );
            } else if( _randomNumberGenerator == FixedValue ){
               action = "loaded new fixed values to GPU";
               _loadRandomNumberStreamsFixedValue( );
            }
            _numberOfShuffles  = 0;
			}
         _rvStreamIndex = 0;
		}

      if( printOn ){
         (void) fprintf( log, "%s offset=%d consume/itr=%d StrmSz=%d idx=%d shffle=%d action=%s\n",
                         methodName.c_str(), _rvStreamOffset, numberOfRandomValuesConsumedPerIteration,
                         rvStreamSize, _rvStreamIndex, _numberOfShuffles, action );
         (void) fflush( log );
      }

/*
      // check rng distribution

      static int count = 0;
      if( count++ < 2 ){

         // accumulate rng -- stats will be in cumulative fields in stat string

         int testIterations = 1000;
         for( int ii = 0; ii < testIterations; ii++ ){
            //_loadRandomNumberStreamsKiss( );
            _loadRandomNumberStreamsMersenne( );
            getStatisticsString();
         }
         (void) fprintf( log, "%s: stats\n%s\n", methodName.c_str(), getStatisticsString().c_str() );
      }
*/
	}

   return DefaultReturnValue;
}

/** 
 * Get random number stream size
 *
 * @return  random number stream size
 *
 */

int BrookRandomNumberGenerator::getRandomNumberStreamSize( void ) const {
   return _randomNumberStreamSize;
}

/** 
 * Get Shuffle stream 
 *
 * @return  shuffle stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::_getShuffleStream( void ) const {
   return _auxiliaryStreams[ShuffleStream];
}

/** 
 * Get random number stream 
 *
 * @param index     random number stream index     
 *
 * @return  random number stream
 *
 * @throw OpenMMException if index is invalid or _randomNumberGeneratorStreams not set
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getRandomNumberStream( int index ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nBrookRandomNumberGenerator::getRandomNumberStream";

   // ---------------------------------------------------------------------------------------
	
   if( index < 0 || index >= _numberOfRandomNumberStreams ){
      std::stringstream message;
      message << methodName << " index=" << index << " is out of range: [0," << _numberOfRandomNumberStreams;
      throw OpenMMException( message.str() );
      return NULL;
   }

   if( _randomNumberGeneratorStreams == NULL ){
      std::stringstream message;
      message << methodName << " rv streams not initialized; input index=" << index;
      throw OpenMMException( message.str() );
      return NULL;
   }

   return _randomNumberGeneratorStreams[index];
}

/** 
 * Get random value stream index
 *
 * @return  random value stream index
 *
 */

int BrookRandomNumberGenerator::getRvStreamIndex( void ) const {
   return _rvStreamIndex;
}

/** 
 * Get random value stream offset
 *
 * @return  random value stream offset
 *
 */

int BrookRandomNumberGenerator::getRvStreamOffset( void ) const {
   return _rvStreamOffset;
}

/** 
 * Get max shuffles
 *
 * @return  max shuffles
 *
 */

int BrookRandomNumberGenerator::getMaxShuffles( void ) const {
   return _maxShuffles;
}

/** 
 * Get number of shuffles
 *
 * @return  number of shuffles
 *
 */

int BrookRandomNumberGenerator::_getNumberOfShuffles( void ) const {
   return _numberOfShuffles;
}

/** 
 * Initialize stream dimensions
 * 
 * @param numberOfParticles         number of particles
 * @param platform                  platform
 *      
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookRandomNumberGenerator::_initializeStreamSizes( int numberOfParticles, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   setNumberOfParticles( numberOfParticles );

   // get randomNumber stream dimensions

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory() );
   _randomNumberStreamWidth                     = brookStreamFactory.getDefaultRandomNumberStreamWidth();
   _randomNumberStreamSize                      = brookStreamFactory.getDefaultRandomNumberStreamSize();
   _randomNumberStreamHeight                    = (int) ( ((float) _randomNumberStreamSize)/( (float) _randomNumberStreamWidth) + 0.001);

   _randomNumberStreamSize                      =  _randomNumberStreamWidth*_randomNumberStreamHeight;

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

int BrookRandomNumberGenerator::_initializeStreams( const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::_initializeStreams";

   float dangleValue                          = 0.0f;

// ---------------------------------------------------------------------------------------

   int randomNumberStreamSize                = getRandomNumberStreamSize();
   int randomNumberStreamWidth               = getRandomNumberStreamWidth();

   _auxiliaryStreams[ShuffleStream]         = new BrookFloatStreamInternal( BrookCommon::ShuffleStream,
                                                                            randomNumberStreamSize, randomNumberStreamWidth,
                                                                            BrookStreamInternal::Float, dangleValue );

   // insure number of random number streams is > 0
   // delete if already allocated and then initialize
 
   if( _numberOfRandomNumberStreams < 1 ){
      _numberOfRandomNumberStreams     = 1;
   }

   if( _randomNumberGeneratorStreams ){
      delete[] _randomNumberGeneratorStreams;
   }

   _randomNumberGeneratorStreams = new BrookFloatStreamInternal*[_numberOfRandomNumberStreams];

   for( int ii = 0; ii < _numberOfRandomNumberStreams; ii++ ){
      _randomNumberGeneratorStreams[ii]    = new BrookFloatStreamInternal( BrookCommon::RandomValuesStream,
                                                                           randomNumberStreamSize, randomNumberStreamWidth,
                                                                           BrookStreamInternal::Float3, dangleValue );
   }

   // create shuffle stream

   _loadGVShuffle();

   return DefaultReturnValue;
}

/*  
 * Setup of streams, ... associated w/ random number generator
 *
 * @param numberOfParticles     number of particles
 * @param platform              Brook platform
 *
 * @return nonzero value if error
 *
 * */
    
int BrookRandomNumberGenerator::setup( int numberOfParticles, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::setup";

// ---------------------------------------------------------------------------------------

   const BrookPlatform brookPlatform            = dynamic_cast<const BrookPlatform&> (platform);
   setLog( brookPlatform.getLog() );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfParticles, platform );
   _initializeStreams( platform );

   if( _randomNumberGenerator == Mersenne ){
      _loadRandomNumberStreamsMersenne( );
   } else if( _randomNumberGenerator == Kiss ){
      _loadRandomNumberStreamsKiss( );
   } else if( _randomNumberGenerator == Original ){
      _loadGVStreamsOriginal( );
   } else if( _randomNumberGenerator == FixedValue ){
      _loadRandomNumberStreamsFixedValue( );
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

std::string BrookRandomNumberGenerator::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRandomNumberGenerator::getContentsString";

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

   
   if( _randomNumberGenerator == Mersenne ){
      (void) LOCAL_SPRINTF( value, "%s", "Mersenne Rng");
   } else if( _randomNumberGenerator == Kiss ){
      (void) LOCAL_SPRINTF( value, "%s", "Kiss Rng");
   } else if( _randomNumberGenerator == Original ){
      (void) LOCAL_SPRINTF( value, "%s", "Original Gromacs Rng");
   } else if( _randomNumberGenerator == FixedValue ){
      (void) LOCAL_SPRINTF( value, "%s", "Fixed value Rng");
   }
   message << _getLine( tab, "Random number generator:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getRandomNumberStreamSize() );
   message << _getLine( tab, "RandomNumber stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getRandomNumberStreamWidth() );
   message << _getLine( tab, "RandomNumber stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getRandomNumberStreamHeight() );
   message << _getLine( tab, "RandomNumber stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfRandomNumberStreams() );
   message << _getLine( tab, "Number of RandomNumber streams:", value ); 

   (void) LOCAL_SPRINTF( value, "%ld", getRandomNumberSeed() );
   message << _getLine( tab, "Seed:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getRvStreamIndex() );
   message << _getLine( tab, "Rv stream index:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getRvStreamOffset() );
   message << _getLine( tab, "Rv stream offset:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", _getNumberOfShuffles() );
   message << _getLine( tab, "Number of rv stream shuffles:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getMaxShuffles() );
   message << _getLine( tab, "Max number of rv stream shuffles:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                 ? Set : NotSet) ); 

   message << _getLine( tab, "Shuffle:",              (_getShuffleStream()      ? Set : NotSet) ); 
 
   // show stats

   message << getStatisticsString( );

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _auxiliaryStreams[ii] ){
         message << _auxiliaryStreams[ii]->getContentsString( );
      }
   }

   for( int ii = 0; ii < _numberOfRandomNumberStreams; ii++ ){
      message << std::endl;
      if( _randomNumberGeneratorStreams[ii] ){
         message <<  _randomNumberGeneratorStreams[ii]->getContentsString();
      }
   }
#undef LOCAL_SPRINTF

   return message.str();
}

/* 
 * Get statistics
 *
 * @return string containing contents
 *
 * */

std::string BrookRandomNumberGenerator::getStatisticsString( void ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRandomNumberGenerator::getStatisticsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   static const char* Set                   = "Set";
   static const char* NotSet                = "Not set";
   static double cumulativeStatistics[7]    = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e+99, -1.0e+99 };

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   std::string tab   = "   ";

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#endif

   message << "\n   Statistics:\n";

   double statistics[7];
   getStatistics( statistics, -1, cumulativeStatistics );

   (void) LOCAL_SPRINTF( value, "%.5e", statistics[4] );
   message << _getLine( tab, "Count:",  value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", statistics[0] );
   message << _getLine( tab, "Average:",  value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", statistics[1] );
   message << _getLine( tab, "StdDev:",  value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", statistics[3] );
   message << _getLine( tab, "Kurtosis:",  value ); 

   (void) LOCAL_SPRINTF( value, "%.5e", statistics[5] );
   message << _getLine( tab, "Min:",  value ); 

   (void) LOCAL_SPRINTF( value, "%.6e", statistics[6] );
   message << _getLine( tab, "Max:",  value ); 

   if( cumulativeStatistics[4] > 1000.0 ){ 

      for( int ii = 0; ii < 7; ii++ ){
         statistics[ii] = cumulativeStatistics[ii];
      }

      statistics[0] /= statistics[4];
      statistics[1]  = statistics[1] - statistics[4]*statistics[0]*statistics[0];
      if( statistics[4] > 1.0 ){
         statistics[1] = sqrt( statistics[1]/( statistics[4] - 1.0 ) );
      }
      statistics[3]  = (statistics[3]/(statistics[4]*statistics[1]*statistics[1]) ) - 3.0;
      (void) LOCAL_SPRINTF( value, "%.5e", statistics[4] );
      message << _getLine( tab, "Cumulative Count:",  value ); 
   
      (void) LOCAL_SPRINTF( value, "%.5e", statistics[0] );
      message << _getLine( tab, "Cumulative Average:",  value ); 
   
      (void) LOCAL_SPRINTF( value, "%.5e", statistics[1] );
      message << _getLine( tab, "Cumulative StdDev:",  value ); 
   
      (void) LOCAL_SPRINTF( value, "%.5e", statistics[3] );
      message << _getLine( tab, "Cumulative Kurtosis:",  value ); 
   
      (void) LOCAL_SPRINTF( value, "%.5e", statistics[5] );
      message << _getLine( tab, "Cumulative Min:",  value ); 
   
      (void) LOCAL_SPRINTF( value, "%.6e", statistics[6] );
      message << _getLine( tab, "Cumulative Max:",  value ); 
   }
   
#undef LOCAL_SPRINTF

   return message.str();
}

/* 
 * Get statistics
 *
 * @param statistics   array of size 5:
 *                       0: mean
 *                       1: std dev
 *                       2: 3rd moment (not normalized)
 *                       3: kurtosis
 *                       4: count 
 *                       5: min
 *                       6: max
 *
 * @param streamIndex  stream index to analyze
 *
 * @return DefaultReturnValue
 *
 * */

int BrookRandomNumberGenerator::getStatistics( double statistics[7], int streamIndex, double cumulativeStatistics[7] ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRandomNumberGenerator::";

// ---------------------------------------------------------------------------------------

   statistics[0] = 0.0;
   statistics[1] = 0.0;
   statistics[2] = 0.0;
   statistics[3] = 0.0;
   statistics[4] = 0.0;
   statistics[5] =  1.0e+99;
   statistics[6] = -1.0e+99;

   for( int ii = 0; ii < _numberOfRandomNumberStreams; ii++ ){
      if( streamIndex < 0 || ii == streamIndex ){
         void* dataArrayV       = _randomNumberGeneratorStreams[ii]->getData( 1 );
         int numberOfValues     = _randomNumberGeneratorStreams[ii]->getStreamSize()*_randomNumberGeneratorStreams[ii]->getWidth();
         const float* dataArray = (float*) dataArrayV;
         for( int ii = 0; ii < numberOfValues; ii++ ){

             statistics[0] += dataArray[ii];

             double rv2     = dataArray[ii]*dataArray[ii];
             statistics[1] += rv2;
             statistics[2] += rv2*dataArray[ii];
             statistics[3] += rv2*rv2;

             if( statistics[5] > dataArray[ii] ){
                statistics[5] = dataArray[ii];
             }
             if( statistics[6] < dataArray[ii] ){
                statistics[6] = dataArray[ii];
             }
         } 
         statistics[4] += (double) numberOfValues;
      }   
   }   

   // accumulate moments, ... in cumulativeStatistics array

   for( int ii = 0; ii < 5; ii++ ){
      cumulativeStatistics[ii] += statistics[ii];
   }
   if( statistics[5] < cumulativeStatistics[5] ){
      cumulativeStatistics[5] = statistics[5];
   }
   if( statistics[6] > cumulativeStatistics[6] ){
      cumulativeStatistics[6] = statistics[6];
   }

   if( statistics[4] > 0.0 ){ 
      statistics[0] /= statistics[4];
      statistics[1]  = statistics[1] - statistics[4]*statistics[0]*statistics[0];
      if( statistics[4] > 1.0 ){
         statistics[1] = sqrt( statistics[1]/( statistics[4] - 1.0 ) );
      }
      statistics[3]  = (statistics[3]/(statistics[4]*statistics[1]*statistics[1]) ) - 3.0;
   }

   return DefaultReturnValue;
}
