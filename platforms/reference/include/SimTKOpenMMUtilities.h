
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

#ifndef __SimTKOpenMMUtilities_H_
#define __SimTKOpenMMUtilities_H_

// class of shared, static utility methods

#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include "openmm/internal/windowsExport.h"

#include <string>

namespace OpenMM {

/**---------------------------------------------------------------------------------------

   Class of static methods to be shared
   Most methods are standalone 'utility' methods

--------------------------------------------------------------------------------------- */

class OPENMM_EXPORT SimTKOpenMMUtilities {

   private:

       static uint32_t _randomNumberSeed;
       static bool _randomInitialized;
       static bool nextGaussianIsValid;
       static double nextGaussian;
       static OpenMM_SFMT::SFMT sfmt;

   public:

      // dummy constructor/destructor

       SimTKOpenMMUtilities() {};
      ~SimTKOpenMMUtilities() {};
      
      /**---------------------------------------------------------------------------------------
      
         Allocate 1D double array (Simbios)
      
         array[i]
      
         @param iSize                i-dimension
         @param array1D              array (if null on entry allocated)
         @param initialize           if true, then initialize array
         @param initialValue         intitial value
         @param idString             id string
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static double* allocateOneDRealOpenMMArray(int iSize, double* array1D, int initialize,
                                                      double initialValue,
                                                      const std::string& idString = std::string("1DArray"));
      
      /**---------------------------------------------------------------------------------------
      
         Allocate 2D double array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialize           if true, then initialize array
         @param initialValue         intitial value
         @param idString             id string
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static double** allocateTwoDRealOpenMMArray(int iSize, int jSize,
                                                       double** array2D, int initialize,
                                                       double initialValue,
                                                       const std::string& idString = std::string("2DArray"));
      
      /* ---------------------------------------------------------------------------------------
      
         Free 2D double array (Simbios)
      
         array[i][j]
      
         @param array2D              array (if null on entry allocated)
         @param idString             id string
      
         --------------------------------------------------------------------------------------- */
      
      static void freeOneDRealOpenMMArray(double* array1D,
                                          const std::string& idString = std::string("1DArray"));
      
      /* ---------------------------------------------------------------------------------------
      
         Free 2D double array (Simbios)
      
         array[i][j]
      
         @param array2D              array (if null on entry allocated)
         @param idString             id string
      
         --------------------------------------------------------------------------------------- */
      
      static void freeTwoDRealOpenMMArray(double** array2D,
                                          const std::string& idString = std::string("2DArray"));
      
      /**---------------------------------------------------------------------------------------
      
         Initialize 2D double array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialValue         intitial value
      
         --------------------------------------------------------------------------------------- */
      
      static void initialize2DRealOpenMMArray(int iSize, int jSize,
                                              double** array2D, double initialValue);

      /**---------------------------------------------------------------------------------------
      
         Compute cross product of two 3-vectors and place in 3rd vector  -- helper method
      
         @param vectorZ = vectorX x vectorY
      
         @param vectorX             x-vector
         @param vectorY             y-vector
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void crossProductVector3(double* vectorX, double* vectorY, double* vectorZ);
      
      /**---------------------------------------------------------------------------------------
      
         Get normally distributed random number
      
         @return random value
      
         --------------------------------------------------------------------------------------- */
      
      static double getNormallyDistributedRandomNumber();
      
      /**---------------------------------------------------------------------------------------
      
         Get uniformly distributed random number in the range [0, 1)
      
         @return random value
      
         --------------------------------------------------------------------------------------- */
      
      static double getUniformlyDistributedRandomNumber();

      /**---------------------------------------------------------------------------------------
      
         Get random number seed
      
         @return random number seed
      
         --------------------------------------------------------------------------------------- */
      
      static uint32_t getRandomNumberSeed();
      
      /**---------------------------------------------------------------------------------------
      
         Set random number seed
      
         @param seed    new seed value
      
         --------------------------------------------------------------------------------------- */
      
      static void setRandomNumberSeed(uint32_t seed);
      
      /**---------------------------------------------------------------------------------------
      
         Write out the internal state of the random number generator.
      
         @param stream  a stream to write the checkpoint to
      
         --------------------------------------------------------------------------------------- */

      static void createCheckpoint(std::ostream& stream);
      
      /**---------------------------------------------------------------------------------------
      
         Load a checkpoint created by createCheckpoint().
      
         @param stream  a stream to load the checkpoint from
      
         --------------------------------------------------------------------------------------- */

      static void loadCheckpoint(std::istream& stream);
};

} // namespace OpenMM

#endif // __SimTKOpenMMUtilities_H__
