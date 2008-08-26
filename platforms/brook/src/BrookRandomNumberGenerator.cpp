working on shuffleGVStreams
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
#include "BrookPlatform.h"
#include "OpenMMException.h"
#include "BrookStreamImpl.h"

// use random number generator

#include "SimTKOpenMMUtilities.h"

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

   // mark stream dimension variables as unset

   _randomNumberStreamWidth         = -1;
   _randomNumberStreamHeight        = -1;
   _randomNumberStreamSize          = -1;

   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      _randomNumberStreams[ii]   = NULL;
   }

   // set randomNumber seed 

   _randomNumberSeed  = 1393;

   //_randomNumberSeed = randomNumberSeed ? randomNumberSeed : 1393;
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
      delete _sdStreams[ii];
   }

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
 * Get number of random number streams
 * 
 * @return     number of random number streams 
 *
 */

int BrookRandomNumberGenerator::getNumberOfRandomNumberStreams( void ) const {
   return _numberOfRandomNumberStreams;
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
    
unsigned long int setRandomNumberSeed( unsigned long int seed = 1 );
   _randomNumberSeed = seed;
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

BrookOpenMMFloat BrookRandomNumberGenerator::generateGromacsRandomNumber( int* ig ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookRandomNumberGenerator::generateGromacsRandomNumber";

   int  irand;
 
   int  m     = 100000000;
   float rm   = 100000000.0;  /* same number as m, but real format */
   int  m1    = 10000;
   int  mult  = 31415821;
   
   BrookOpenMMFloat r;
   int  irandh,irandl,multh,multl;
 
   // ---------------------------------------------------------------------------------------
 
   irand = abs(*ig) % m; 
   
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
 * Generate a random number using algorithm in Nvidia code
 * http://www.helsbreth.org/random/rng_kiss.html
 * 
 * @param randomV1   output random value
 * @param randomV2   output random value
 * @param randomV3   output random value
 * @param state      state
 *
 */

void BrookRandomNumberGenerator::generateRandomsAlaNvidia( float* randomV1, float* randomV2, float* randomV3, 
                                                           unsigned int state[4] ){
    
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookRandomNumberGenerator::generateRandomsAlaNvidia";

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
 * Load random number streams using Nvidia algorithm
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::loadRandomNumberStreamsNvidia( void ){

   // ---------------------------------------------------------------------------------------

	static float *buf              = NULL;
   static unsigned int state[4];
   static int stateInitialized    = 0;
   static const int reseed        = 5;
 
   // static const char* methodName = "\nBrookRandomNumberGenerator::loadGVStreamsNvidia";

   // ---------------------------------------------------------------------------------------
   
   // periodically reset seeds

   if( !stateInitialized || !(stateInitialized % reseed) ){

      state[0] = rand();
      state[1] = rand();
      state[2] = rand();
      state[3] = rand();

      if( getVerbosity() && getLog() ){
         (void) fprintf( getLog(), "LoadGVStreamsNvidia: reset state seeds stateInitialized=%d reseed=%d\n",
                         stateInitialized, reseed );
         (void) fflush( getLog() );
      }

/*
state[0] = 9578;
state[1] = 29245;
state[2] = 16266;
state[3] = 27587;
*/

   }
   stateInitialized++;

   // allocate memory once for download of random nos.

   if( buf == NULL ){	
	   buf = (float*) malloc( sizeof(float) * 3 * getRandomNumberStreamSize() );
   }

   if( getVerbosity() && getLog() ){
	   static float count   = 0.0f;
      float block          = (float) (3*sdp->gvSize);
      count               += 1.0f;
      (void) fprintf( getLog(), "LoadGVStreamsNvidia: count=%.1f ttl=%.3e no./count=%.1f %d %d\n",
                      count, block*count, block, sdp->gvSize, NGVSTREAMS );
      (void) fflush( getLog() );
   }
	
   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
      for( int ii = 0; ii < 3*getRandomNumberStreamSize(); ii += 3 ){
         float v1,v2,v3;
         generateRandomsAlaNvidia( &v1, &v2, &v3, state, NULL );
         buf[ii]   = v1;
         buf[ii+1] = v2;
         buf[ii+2] = v3;
   	}
	   getRandomNumberStream( jj )->loadFromArray( buf );
   }

   return DefaultReturnValue;
}

/** 
 * Load random number streams using original gpu algorithm
 * 
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::loadGVStreamsOriginal( void ){

   // ---------------------------------------------------------------------------------------

	int i, j;
	static float *buf = NULL;
	unsigned long int jran;

   // static const char* methodName = "\nBrookRandomNumberGenerator::LoadGVStreamsOriginal";

   // ---------------------------------------------------------------------------------------
	
   if( buf == NULL ){
      buf = (float*) malloc( sizeof(float) * 3 * sdp->gvSize );
   }
	
	jran = getRandomNumberSeed();
	
   for( int jj = 0; jj < getNumberOfRandomNumberStreams(); jj++ ){
      for( int ii = 0; ii < 3*getRandomNumberStreamSize(); ii += 3 ){
			//buf[i] = sdp->fgauss( &jran );
			buf[i] = generateGromacsRandomNumber( &jran );
		}
	   getRandomNumberStream( jj )->loadFromArray( buf );
	}
	
   incrementRandomNumberSeed( 1 );

   return DefaultReturnValue;
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

int BrookRandomNumberGenerator::loadGVShuffle( void ){

   // ---------------------------------------------------------------------------------------

	const int p[] = { 2, 3, 5, 7, 11 };
	const int np  = sizeof(p) / sizeof(p[0]);
   const int pmax = p[np-1];

	static float *buf = 0;
	int iter, i, j;
	static int *indices; 
	float tmp;

   // static const char* methodName = "\nBrookRandomNumberGenerator::loadGVShuffle";

   // ---------------------------------------------------------------------------------------
	
   if( buf == NULL ){
	   buf = (float*) malloc( sizeof(float) * sdp->gvSize );
	   indices = (int*) malloc( sizeof(int) * pmax );
   }
	
   int rvSize = getRandomNumberStreamSize();
	for ( i = 0; i < rvSize; i++ ) {
		buf[i] = (float) i;
	}

	//How to come up with this number here?

   unsigned long int seed = getRandomNumberSeed();
	for ( iter = 0; iter < 1000000; iter++ ) {
		//for each p
		for ( i = 0; i < np; i++ ){
			//Come up with p random indices
			//Not checking that they are distinct
			//because that should be fairly rare
			for ( j = 0; j < p[i]; j++ ) {
				//indices[j] = (int) ( gmx_rando( &sdp->seed ) * sdp->gvSize );
				indices[j] = (int) ( generateGromacsRandomNumber( &seed )*rvSize );
			}
			//do a p-cycle
			tmp = buf[ indices[0] ];
			for ( j = 0; j < p[i]-1; j++ ) {
				buf[ indices[j] ] = buf[ indices[j+1] ];
			}
			buf[ indices[ p[i]-1 ] ] = tmp;
		}
	}
   getShuffleStream()->loadFromArray( buf );
	
   return DefaultReturnValue;
}

/** 
 * Shuffle streams
 *
 * @return DefaultReturnValue;
 */

int BrookRandomNumberGenerator::shuffleGVStreams( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookRandomNumberGenerator::shuffleGVStreams";

   // ---------------------------------------------------------------------------------------
	
   int numberOfRvStreams = getNumberOfRandomNumberStreams();
	for( int ii = 0; ii < numberOfRvStreams - 1; ii++ ){
		kpermute_vectors( (float) , sdp->strShuffle, sdp->strVGauss[i+1], 
		                  sdp->strVGauss[i] );
   }

   for( int jj = 0; jj < numberOfRvStreams - 1; jj++ ){
	   kpermute_vectors( (float) sdp->gvWidth, sdp->strShuffle, sdp->strVGauss[0], 
	                      sdp->strVGauss[ NGVSTREAMS-1 ] );
   }

   return DefaultReturnValue;
}

static int
PrintGVStreams( gpuContext gpu )
{
	int ii;
	
	gpuSDParams sdp = gpu->sdparams;
   char testName[128];
   char fileName[128];
      
   /* dump inputs */

   (void) sprintf( testName, "ShuffleGVStreams" );
   (void) sprintf( fileName, "%s.32bits.in", testName );
   Save32Bits( fileName, "iiiii", sdp->gvWidth, sdp->seed, sdp->gvOffset, sdp->gvCurStream, sdp->gvNumShuffles );  

   (void) sprintf( fileName, "%s.strShuffle.in", testName );
   SaveStream( fileName, sdp->strShuffle );

	for ( ii = 0; ii < NGVSTREAMS; ii++ ){
      (void) sprintf( testName, "ShuffleGVStreams%d", ii );
      (void) sprintf( fileName, "%s.strVGaussIn.in", testName );
      SaveStream( fileName, sdp->strVGauss[ii] );

   }

	return 1;
}

/* Advances the current position by 2*gpu->natoms
 * If we run out of rand numbers, we shuffle and 
 * reuse a few times before reloading from the cpu
 * */

static int
AdvanceGVCursor( gpuSDParams sdp, gpuContext gpu, int step )
{
   static clock_t rvTime      = 0;
   static clock_t shuffleTime = 0;
   static const int trackTime = 0;
   clock_t time;

   if( trackTime ){
      time               = clock();
   }

	sdp->gvOffset += 2*gpu->natoms;
	//Check if we've used up this texture
	if ( sdp->gvOffset > sdp->gvSize - 2*gpu->natoms ) {
		//Next one if available
		sdp->gvOffset = 0;
		if ( sdp->gvCurStream < NGVSTREAMS - 1 ) {
			sdp->gvCurStream++;
			//printf( "Using stream %d\n", sdp->gvCurStream ); 
			//fflush( stdout );
		} else {
			//No more textures, need to shuffle
			if ( sdp->gvNumShuffles < NGVSHUFFLE ) {
				//printf( "Shuffling streams\n" );
				//fflush( stdout );

            clock_t localTime;
            if( trackTime ){
               localTime = clock();
            }

				ShuffleGVStreams( sdp, gpu );

            if( trackTime ){
               shuffleTime += clock() - localTime;
            }
				sdp->gvNumShuffles++;
			} else { //Need to refresh random numbers from cpu
            static int reducePrint = 0;
				fflush( gpu->log );
            if( UseOriginalRng ){
               LoadGVStreamsOriginal( sdp );
            } else {
               LoadGVStreamsNvidia( sdp, gpu->log );
            }
				sdp->gvNumShuffles = 0;

   if( trackTime && !(reducePrint++ % 20) ){
      rvTime += clock() - time;
		(void) fprintf( gpu->log, "Reloading random number streams from cpu: seed=%d at step=%d", sdp->seed, step );
      (void) fprintf( gpu->log, " cmltv times: RValue=%.5f Shuffle=%.5f\n", ((double) rvTime/(double)CLOCKS_PER_SEC), ((double) shuffleTime/(double)CLOCKS_PER_SEC) );
      (void) fflush( gpu->log );
   }

			}
			sdp->gvCurStream = 0;
		}
	}

/*
   if( trackTime ){
      rvTime += clock() - time;
      if( (!sdp->gvNumShuffles && rvTime > 1) || !(sdp->gvNumShuffles %20) && (rvTime > 2*CLOCKS_PER_SEC) ){
         (void) fprintf( gpu->log, "\n   Total RValue=%.5f Shuffle=%.5f", ((double) rvTime/(double)CLOCKS_PER_SEC), ((double) shuffleTime/(double)CLOCKS_PER_SEC) );
      }
   }
*/

	return 1;
}


extern "C" int
gpuResetSdSeed( gpuContext gpu, int step, int seed )
{
	gpuSDParams sdp = gpu->sdparams;
   sdp->seed       = seed;
   (void) fprintf( gpu->log, "gpuResetSdSeed: reset SD seed to %d at step=%d.\n", seed, step );
   (void) fflush( gpu->log );
   return 0;
}
/** 
 * Update parameters -- only way parameters can be set
 * 
 * @param  temperature     temperature
 * @param  friction        friction
 * @param  step size       step size
 *
 * @return   solute dielectric
 *
 */

int BrookRandomNumberGenerator::updateParameters( double temperature, double friction, double stepSize ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nBrookRandomNumberGenerator::updateParameters";

   // ---------------------------------------------------------------------------------------

   _setStepSize(    (BrookOpenMMFloat)  stepSize );
   _setFriction(    (BrookOpenMMFloat)  friction );
   _setTemperature( (BrookOpenMMFloat)  temperature );

   _updateDerivedParameters( );
   _updateSdStreams( );

   return DefaultReturnValue;

};

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

int BrookRandomNumberGenerator::update( Stream& positions, Stream& velocities,
                                     const Stream& forces,
                                     BrookShakeAlgorithm& brookShakeAlgorithm ){

   // ---------------------------------------------------------------------------------------

      // unused Shake parameter

      float omega                   = 1.0f;

   // static const char* methodName = "\nBrookRandomNumberGenerator::update";

   // ---------------------------------------------------------------------------------------

      BrookOpenMMFloat* derivedParameters                 = getDerivedParameters();

      BrookStreamImpl& positionStream                     = dynamic_cast<BrookStreamImpl&>       (positions.getImpl());
      BrookStreamImpl& velocityStream                     = dynamic_cast<BrookStreamImpl&>       (velocities.getImpl());
      BrookStreamImpl& forceStreamC                       = dynamic_cast<BrookStreamImpl&>       (forces.getImpl());
      const BrookStreamImpl& forceStream                  = dynamic_cast<const BrookStreamImpl&> (forceStream.getImpl());


      // first integration step

      kupdate_sd1_fix1(
            (float) getStochasticDynamicsAtomStreamWidth(),
(float) sdp->gvWidth,
(float) sdp->gvOffset,
            derivedParameters[EM],
            derivedParameters[Sd1pc1],
            derivedParameters[Sd1pc2],
            derivedParameters[Sd1pc3],
            getSDPC1Stream()->getBrookStream(),
sdp->strVGauss[ sdp->gvCurStream ],
            getSD2XStream()->getBrookStream(),
            positionStream.getBrookStream(),
            forceStream.getBrookStream(),
            velocityStream.getBrookStream(),
            getInverseMassStream()->getBrookStream(),
            getSD1VStream()->getBrookStream(),
            getVPrimeStream()->getBrookStream(),
            getXPrimeStream()->getBrookStream() 
            );
 
AdvanceGVCursor( sdp, gpu, step );

      // first Shake step

      kshakeh_fix1( 
                    10.0f,
                    (float) getStochasticDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm->getInverseHydrogenMass(),
                    omega, 
                    brookShakeAlgorithm->getShakeAtomIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeAtomParameterStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons3Stream()->getBrookStream() );

      // first Shake gather

      kshakeh_update1_fix1(
                    (float) getStochasticDynamicsAtomStreamWidth(),
                    derivedParameters[Sd2pc1],
                    brookShakeAlgorithm->getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    getVPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons3Stream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream() );

      // second integration step

      kupdate_sd2_fix1(
            (float) getStochasticDynamicsAtomStreamWidth(),
(float) sdp->gvWidth,
(float) sdp->gvOffset,
            derivedParameters[Sd2pc1],
            derivedParameters[Sd2pc2],
            getSDPC2Stream()->getBrookStream(),
sdp->strVGauss[ sdp->gvCurStream ],
            getSD1VStream()->getBrookStream(),
            positionStream.getBrookStream(),
            getXPrimeStream()->getBrookStream(), 
            getVPrimeStream()->getBrookStream(),
            getS2XStream()->getBrookStream(),
            velocityStream.getBrookStream(),
            getXPrimeStream()->getBrookStream() 
            );

AdvanceGVCursor( sdp, gpu, step );

      // second Shake step

      kshakeh_fix1( 
                    10.0f,
                    (float) getStochasticDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm->getInverseHydrogenMass(),
                    omega, 
                    brookShakeAlgorithm->getShakeAtomIndicesStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeAtomParameterStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons3Stream()->getBrookStream() );

      // second Shake gather

      kshakeh_update1_fix1(
                    (float) getStochasticDynamicsAtomStreamWidth(),
                    derivedParameters[Sd2pc1],
                    brookShakeAlgorithm->getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    getVPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons3Stream()->getBrookStream(),
                    getXPrimeStream()->getBrookStream() );

      kshakeh_update2_fix1( 
                    (float) getStochasticDynamicsAtomStreamWidth(),
                    brookShakeAlgorithm->getShakeInverseMapStream()->getBrookStream(),
                    positionStream.getBrookStream(),
                    getXPrimeStream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons0Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons1Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons2Stream()->getBrookStream(),
                    brookShakeAlgorithm->getShakeXCons3Stream()->getBrookStream(),
                    positionStream.getBrookStream() );

   return DefaultReturnValue;

};

/**
 *
 * Get array of derived parameters indexed by 'DerivedParameters' enums
 *
 * @return array
 *
 */
   
const BrookOpenMMFloat* BrookRandomNumberGenerator::getDerivedParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName  = "\nBrookRandomNumberGenerator::getDerivedParameters";

   // ---------------------------------------------------------------------------------------

   return _derivedParameters;
}

/** 
 * Get Atom stream size
 *
 * @return  Atom stream size
 *
 */

int BrookRandomNumberGenerator::getStochasticDynamicsAtomStreamSize( void ) const {
   return _sdAtomStreamSize;
}

/** 
 * Get atom stream width
 *
 * @return  atom stream width
 *
 */

int BrookRandomNumberGenerator::getStochasticDynamicsAtomStreamWidth( void ) const {
   return _sdAtomStreamWidth;
}

/** 
 * Get atom stream height
 *
 * @return atom stream height
 */

int BrookRandomNumberGenerator::getStochasticDynamicsAtomStreamHeight( void ) const {
   return _sdAtomStreamHeight;
}

/** 
 * Get SDPC1 stream 
 *
 * @return  SDPC1 stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getSDPC1Stream( void ) const {
   return _sdStreams[SDPC1Stream];
}

/** 
 * Get SDPC2 stream 
 *
 * @return  SDPC2 stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getSDPC2Stream( void ) const {
   return _sdStreams[SDPC2Stream];
}

/** 
 * Get SD2X stream 
 *
 * @return  SD2X stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getSD2XStream( void ) const {
   return _sdStreams[SD2XStream];
}

/** 
 * Get SD1V stream 
 *
 * @return  SD1V stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getSD1VStream( void ) const {
   return _sdStreams[SD1VStream];
}

/** 
 * Get VPrime stream 
 *
 * @return  Vprime stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getVPrimeStream( void ) const {
   return _sdStreams[VPrimeStream];
}

/** 
 * Get XPrime stream 
 *
 * @return  Xprime stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getXPrimeStream( void ) const {
   return _sdStreams[XPrimeStream];
}

/** 
 * Get InverseMass stream 
 *
 * @return  inverse mass stream
 *
 */

BrookFloatStreamInternal* BrookRandomNumberGenerator::getInverseMassStream( void ) const {
   return _sdStreams[InverseMassStream];
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

int BrookRandomNumberGenerator::_initializeStreamSizes( int numberOfAtoms, const Platform& platform ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::_initializeStreamSizes";

// ---------------------------------------------------------------------------------------

   _sdAtomStreamSize     = getAtomStreamSize( platform );
   _sdAtomStreamWidth    = getAtomStreamWidth( platform );
   _sdAtomStreamHeight   = getAtomStreamHeight( platform );

   return DefaultReturnValue;
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

//std::string BrookRandomNumberGenerator::_getDerivedParametersString( BrookRandomNumberGenerator::DerivedParameters  derivedParametersIndex ) const {
std::string BrookRandomNumberGenerator::_getDerivedParametersString( int derivedParametersIndex ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::_getDerivedParametersString";

// ---------------------------------------------------------------------------------------

   std::string returnString;
 
   switch( derivedParametersIndex ){
 
      case GDT:
         returnString = "GDT";
         break;
 
      case EPH:
         returnString = "EPH";
         break;
 
      case EMH:
         returnString = "EMH";
         break;
 
      case EP:
         returnString = "EP";
         break;
 
      case EM:
         returnString = "EM";
         break;
 
      case B:
         returnString = "B";
         break;
 
      case C:
         returnString = "C";
         break;
 
      case D:
         returnString = "D";
         break;
 
      case V:
         returnString = "V";
         break;
 
      case X:
         returnString = "X";
         break;
 
      case Yv:
         returnString = "Yv";
         break;
 
      case Yx:
         returnString = "Yx";
         break;
 
      case Sd1pc1:
         returnString = "Sd1pc1";
         break;
 
      case Sd1pc2:
         returnString = "Sd1pc2";
         break;
 
      case Sd1pc3:
         returnString = "Sd1pc3";
         break;
 
      case Sd2pc1:
         returnString = "Sd2pc1";
         break;
 
      case Sd2pc2:
         returnString = "Sd2pc2";
         break;
 
      default:
         returnString = "Unknown";
         break;
    }
 
    return returnString;
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

   BrookOpenMMFloat dangleValue            = (BrookOpenMMFloat) 0.0;

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize             = getStochasticDynamicsAtomStreamSize();
   int sdAtomStreamWidth            = getStochasticDynamicsAtomStreamWidth();

    _sdStreams[SDPC1Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC1Stream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SDPC2Stream]         = new BrookFloatStreamInternal( BrookCommon::SDPC2Stream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float2, dangleValue );

    _sdStreams[SD2XStream]          = new BrookFloatStreamInternal( BrookCommon::SD2XStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[SD1VStream]          = new BrookFloatStreamInternal( BrookCommon::SD1VStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[VPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::VPrimeStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[XPrimeStream]        = new BrookFloatStreamInternal( BrookCommon::XPrimeStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float3, dangleValue );

    _sdStreams[InverseMassStream]   = new BrookFloatStreamInternal( BrookCommon::InverseMassStream,
                                                                    sdAtomStreamSize, sdAtomStreamWidth,
                                                                    BrookStreamInternal::Float, dangleValue );

   return DefaultReturnValue;
}

/** 
 * Update sd streams -- called after parameters change
 * 
 * @return ErrorReturnValue if error, else DefaultReturnValue
 *
 */

int BrookRandomNumberGenerator::_updateSdStreams( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName       = "BrookRandomNumberGenerator::_updateSdStreams";

// ---------------------------------------------------------------------------------------

   int sdAtomStreamSize                      = getStochasticDynamicsAtomStreamSize();

   BrookOpenMMFloat* sdpc[2];
   for( int ii = 0; ii < 2; ii++ ){
      sdpc[ii] = new BrookOpenMMFloat[2*sdAtomStreamSize];
      memset( sdpc[ii], 0, 2*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 
   }

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters( );
   int numberOfAtoms                         = getNumberOfAtoms();
   int index                                 = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++, index += 2 ){

      sdpc[0][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yv]) );
      sdpc[0][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[V])  );

      sdpc[1][index]      = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[Yx]) );
      sdpc[1][index+1]    = _inverseSqrtMasses[ii]*( static_cast<BrookOpenMMFloat> (derivedParameters[X])  );

   }

   _sdStreams[SDPC1Stream]->loadFromArray( sdpc[0] );
   _sdStreams[SDPC2Stream]->loadFromArray( sdpc[1] );

   for( int ii = 0; ii < 2; ii++ ){
      delete[] sdpc[ii];
   }

   // initialize SD2X

   BrookOpenMMFloat* sd2x = new BrookOpenMMFloat[3*sdAtomStreamSize];
   //SimTKOpenMMUtilities::setRandomNumberSeed( (uint32_t) getRandomNumberSeed() );

   memset( sd2x, 0, 3*sdAtomStreamSize*sizeof( BrookOpenMMFloat ) ); 

   index = 0;
   for( int ii = 0; ii < numberOfAtoms; ii++, index += 3 ){
      sd2x[index]        = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+1]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
      sd2x[index+2]      = _inverseSqrtMasses[ii]*derivedParameters[X]*( static_cast<BrookOpenMMFloat> (SimTKOpenMMUtilities::getNormallyDistributedRandomNumber()) );
   }
   
   _sdStreams[SD2XStream]->loadFromArray( sd2x );

   delete[] sd2x;

   return DefaultReturnValue;

}

/** 
 * Set masses 
 * 
 * @param masses             atomic masses
 *
 */

int BrookRandomNumberGenerator::_setInverseSqrtMasses( const std::vector<double>& masses ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookRandomNumberGenerator::_setInverseSqrtMasses";

   BrookOpenMMFloat zero                    = (BrookOpenMMFloat)  0.0;
   BrookOpenMMFloat one                     = (BrookOpenMMFloat)  1.0;

// ---------------------------------------------------------------------------------------

   // setup inverse sqrt masses

   _inverseSqrtMasses = new BrookOpenMMFloat[masses.size()];
   int index          = 0;
   for( std::vector<double>::const_iterator ii = masses.begin(); ii != masses.end(); ii++, index++ ){
      if( *ii != 0.0 ){
         BrookOpenMMFloat value    = static_cast<BrookOpenMMFloat>(*ii);
         _inverseSqrtMasses[index] = ( SQRT( one/value ) );
      } else {
         _inverseSqrtMasses[index] = zero;
      }
   }

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
    
int BrookRandomNumberGenerator::setup( const std::vector<double>& masses, const Platform& platform ){
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookRandomNumberGenerator::setup";

// ---------------------------------------------------------------------------------------

   int numberOfAtoms  = (int) masses.size();
   setNumberOfAtoms( numberOfAtoms );

   // set stream sizes and then create streams

   _initializeStreamSizes( numberOfAtoms, platform );
   _initializeStreams( platform );

   _setInverseSqrtMasses( masses );

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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTau() );
   message << _getLine( tab, "Tau:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getStepSize() );
   message << _getLine( tab, "Step size:", value ); 

   const BrookOpenMMFloat* derivedParameters = getDerivedParameters();
   for( int ii = 0; ii < MaxDerivedParameters; ii++ ){
      (void) LOCAL_SPRINTF( value, "%.5e", derivedParameters[ii] );
      message << _getLine( tab, _getDerivedParametersString( ii ), value ); 
   }

   (void) LOCAL_SPRINTF( value, "%.5f", getTemperature() );
   message << _getLine( tab, "Temperature:", value ); 

   message << _getLine( tab, "Log:",                  (getLog()                ? Set : NotSet) ); 

   message << _getLine( tab, "SDPC1:",                (getSDPC1Stream()        ? Set : NotSet) ); 
   message << _getLine( tab, "SDPC2:",                (getSDPC2Stream()        ? Set : NotSet) ); 
   message << _getLine( tab, "SD2X:",                 (getSD2XStream()         ? Set : NotSet) ); 
   message << _getLine( tab, "SD1V:",                 (getSD1VStream()         ? Set : NotSet) ); 
 
   for( int ii = 0; ii < LastStreamIndex; ii++ ){
      message << std::endl;
      if( _sdStreams[ii] ){
         message << _sdStreams[ii]->getContentsString( );
      }
   }

#undef LOCAL_SPRINTF

   return message.str();
}
