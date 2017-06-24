
/* Portions copyright (c) 2006-2009 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// class of shared, static utility methods

#include "openmm/internal/OSRngSeed.h"
#include "SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"

// fabs(), ...

#include <cmath>
#include <cstdio>
#include <string.h>
#include <iostream>

using namespace OpenMM;

uint32_t SimTKOpenMMUtilities::_randomNumberSeed = 0;
bool SimTKOpenMMUtilities::_randomInitialized = false;
bool SimTKOpenMMUtilities::nextGaussianIsValid = false;
double SimTKOpenMMUtilities::nextGaussian = 0;
OpenMM_SFMT::SFMT SimTKOpenMMUtilities::sfmt;

/* ---------------------------------------------------------------------------------------

   Allocate 1D double array (Simbios)

   array[i]

   @param iSize                i-dimension
   @param array1D              array (if null on entry allocated)
   @param initialize           if true, then initialize array
   @param initialValue         intitial value
   @param idString             id string

   @return allocated array

   --------------------------------------------------------------------------------------- */

double* SimTKOpenMMUtilities::allocateOneDRealOpenMMArray(int iSize, double* array1D, 
                                                              int initialize, double initialValue,
                                                              const std::string& idString) {
   if (array1D == NULL) {

      array1D = new double[iSize];

   }

   if (initialize) {
      if (initialValue == 0.0) {
         memset(array1D, 0, iSize*sizeof(double));
      } else {
         for (int ii = 0; ii < iSize; ii++) {
            array1D[ii] = initialValue;
         }
      }
   }

   return array1D;
}

/* ---------------------------------------------------------------------------------------

   Allocate 2D double array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialize           if true, then initialize array
   @param initialValue         intitial value
   @param idString             id string

   @return allocated array

   --------------------------------------------------------------------------------------- */

double** SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray(int iSize, int jSize, double** array2D, 
                                                               int initialize, double initialValue,
                                                               const std::string& idString) {
   if (array2D == NULL) {

      array2D = new double*[iSize];
      std::string blockString = idString;
      blockString.append("Block");

      double* block = new double[jSize*iSize];

      for (int ii = 0; ii < iSize; ii++) {
         array2D[ii]  = block;
         block       += jSize;
      }    
   }

   if (initialize) {
      initialize2DRealOpenMMArray(iSize, jSize, array2D, initialValue);
   }

   return array2D;
}

/* ---------------------------------------------------------------------------------------

   Free 2D double array (Simbios)

   array[i][j]

   @param array2D              array (if null on entry allocated)
   @param idString             id string

   --------------------------------------------------------------------------------------- */

void SimTKOpenMMUtilities::freeTwoDRealOpenMMArray(double** array2D, const std::string& idString) {
   if (array2D != NULL) {

      std::string blockString = idString;
      blockString.append("Block");

      delete[] array2D[0];
      delete[] array2D;
   }
}

/* ---------------------------------------------------------------------------------------

   Free 1D double array (Simbios)

   array[i]

   @param array1D              array (if null on entry allocated)
   @param idString             id string

   --------------------------------------------------------------------------------------- */

void SimTKOpenMMUtilities::freeOneDRealOpenMMArray(double* array1D, const std::string& idString) {
   if (array1D != NULL) {
      delete[] array1D;
   }
}

/* ---------------------------------------------------------------------------------------

   Initialize 2D double array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialValue         intitial value

   --------------------------------------------------------------------------------------- */

void SimTKOpenMMUtilities::initialize2DRealOpenMMArray(int iSize, int jSize,
                                                       double** array2D,
                                                       double initialValue) {
   bool useMemset;
   bool useMemsetSingleBlock;

   if (initialValue == 0.0f) {
      useMemset = true;
      if (jSize > 1 && (array2D[0] + jSize) == array2D[1]) {  
         useMemsetSingleBlock = true;
      } else {
         useMemsetSingleBlock = false;
      }

   } else {
      useMemset = false;
   }

   if (useMemset) {
      if (useMemsetSingleBlock) {
         memset(array2D[0], 0, iSize*jSize*sizeof(double));
      } else {
         for (int ii = 0; ii < iSize; ii++) {
            memset(array2D[ii], 0, jSize*sizeof(double));
         }
      }
   } else {
      for (int ii = 0; ii < iSize; ii++) {
         for (int jj = 0; jj < jSize; jj++) {
            array2D[ii][jj] = initialValue;
         }
      }
   }
}

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector  -- helper method

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
void SimTKOpenMMUtilities::crossProductVector3(double* vectorX,
                                               double* vectorY,
                                               double* vectorZ) {
   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

   return;
}

/**---------------------------------------------------------------------------------------

   Get normally distributed random number

   @return random value

   --------------------------------------------------------------------------------------- */

double SimTKOpenMMUtilities::getNormallyDistributedRandomNumber() {
    if (nextGaussianIsValid) {
        nextGaussianIsValid = false;
        return nextGaussian;
    }
    if (!_randomInitialized) {
        init_gen_rand(_randomNumberSeed, sfmt);
        _randomInitialized = true;
        nextGaussianIsValid = false;
    }
    
    // Use the polar form of the Box-Muller transformation to generate two Gaussian random numbers.
    
    double x, y, r2;
    do {
        x = 2.0 * genrand_real2(sfmt) - 1.0;
        y = 2.0 * genrand_real2(sfmt) - 1.0;
        r2 = x*x + y*y;
    } while (r2 >= 1.0 || r2 == 0.0);
    double multiplier = sqrt((-2.0*log(r2))/r2);
    nextGaussian = y*multiplier;
    nextGaussianIsValid = true;
    return x*multiplier;
}

/**---------------------------------------------------------------------------------------

   Get uniformly distributed random number in the range [0, 1)

   @return random value

   --------------------------------------------------------------------------------------- */

double SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() {
    if (!_randomInitialized) {
        init_gen_rand(_randomNumberSeed, sfmt);
        _randomInitialized = true;
        nextGaussianIsValid = false;
    }
    double value = genrand_real2(sfmt);
    return value;
}

/**---------------------------------------------------------------------------------------

   Get random number seed

   @return random number seed

   --------------------------------------------------------------------------------------- */

uint32_t SimTKOpenMMUtilities::getRandomNumberSeed() {
   return _randomNumberSeed;
}

/**---------------------------------------------------------------------------------------

   Set random number seed

   @param seed    new seed value

   --------------------------------------------------------------------------------------- */

void SimTKOpenMMUtilities::setRandomNumberSeed(uint32_t seed) {
    // If the seed is 0, use a unique seed
    if (seed == 0)
        _randomNumberSeed = (uint32_t) osrngseed();
    else
        _randomNumberSeed = seed;
   _randomInitialized = false;
}

void SimTKOpenMMUtilities::createCheckpoint(std::ostream& stream) {
    stream.write((char*) &_randomNumberSeed, sizeof(uint32_t));
    stream.write((char*) &_randomInitialized, sizeof(bool));
    if (_randomInitialized) {
        stream.write((char*) &nextGaussianIsValid, sizeof(bool));
        stream.write((char*) &nextGaussian, sizeof(double));
        sfmt.createCheckpoint(stream);
    }
}

void SimTKOpenMMUtilities::loadCheckpoint(std::istream& stream) {
    stream.read((char*) &_randomNumberSeed, sizeof(uint32_t));
    bool prevInitialized = _randomInitialized;
    stream.read((char*) &_randomInitialized, sizeof(bool));
    if (_randomInitialized) {
        if (!prevInitialized)
            init_gen_rand(0, sfmt);
        stream.read((char*) &nextGaussianIsValid, sizeof(bool));
        stream.read((char*) &nextGaussian, sizeof(double));
        sfmt.loadCheckpoint(stream);
    }
}
