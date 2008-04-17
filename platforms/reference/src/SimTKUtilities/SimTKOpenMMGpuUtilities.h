
/* Portions copyright (c) 2006 Stanford University and Simbios.
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

#ifndef __SimTKOpenMMGpuUtilities_H_
#define __SimTKOpenMMGpuUtilities_H_

// ---------------------------------------------------------------------------------------

// reserve the option to change between Real and double

#include "SimTKOpenMMRealType.h" 
#include "SimTKOpenMMCommon.h" 

// ---------------------------------------------------------------------------------------

// class of shared, static utility methods

#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>

// #include <string>

// ---------------------------------------------------------------------------------------

/**---------------------------------------------------------------------------------------

   Class of static methods to be shared
   Most methods are standalone 'utility' methods

--------------------------------------------------------------------------------------- */

class SimTKOpenMMGpuUtilities {

   public:

      // dummy constructor/destructor

       SimTKOpenMMGpuUtilities(){};
      ~SimTKOpenMMGpuUtilities(){};

      /**---------------------------------------------------------------------------------------
      
         Helper method to repack RealOpenMM arrays (Simbios)
      
         @param numberOfEntries     entries/sub-array
         @param subarraySize        number of subarrays
         @param array               array
      
         Input array = [subArray_1 subArray_2 subArray_3 ... subArray_Stacked] \n
                where each subArray_i is 1 x numberOfEntries
      
         Output array = [ subArray_1_1 subArray_2_1 subArray_3_1 ... subArray_Stacked_1 \n
                           subArray_1_2 subArray_2_2 subArray_3_2 ... subArray_Stacked_2 \n
                           subArray_1_3 subArray_2_3 subArray_3_3 ... subArray_Stacked_3 \n
                           ... \n
                           subArray_1_N subArray_2_N subArray_3_N ... subArray_Stacked_N ] \n

         where N = numberOfEntries

         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArray1( int numberOfEntries, int subarraySize, RealOpenMM* array );
      
      /**---------------------------------------------------------------------------------------
      
         Copy an array into a packed (e.g., RealOpenMM4 array) (Simbios)
      
         Example copy Born radii into last slot of force array
      
         @param numberOfEntries     entries/sub-array (no. atoms)
         @param subarraySize        number of subarrays (4 for RealOpenMM4)
         @param fullArray           full array (force array in example)
         @param copySlot            index of slot to copied into (3 in example, since want Born radius in .w slot)
         @param arrayToCopy         array to copy (Born array)
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int copySubArray( int numberOfEntries, int subarraySize, RealOpenMM* fullArray,
                               int copySlot, RealOpenMM* arrayToCopy );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to repack RealOpenMM arrays (Simbios)
      
         @param numberOfEntries     entries/sub-array
         @param subarraySize        number of subarrays
         @param array               array
      
         Input array = [subArray_1 subArray_2 subArray_3 ... subArray_N] \n
       
                where each subArray_i is vector of length M=numberOfEntries \n
       
         Output array = [ subArray_1_1 subArray_2_1 subArray_3_1 ... subArray_N_1 \n
                           subArray_1_2 subArray_2_2 subArray_3_2 ... subArray_N_2 \n
                           subArray_1_3 subArray_2_3 subArray_3_3 ... subArray_N_3 \n
                           ... \n
                           subArray_1_M subArray_2_M subArray_3_M ... subArray_M ] \n
         where N = numberOfEntries \n

         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArray( int numberOfEntries, int subarraySize, RealOpenMM* array );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to repack RealOpenMM arrays (Simbios)
      
         @param numberOfEntries     entries/sub-array
         @param subarraySize        number of subarrays
         @param array               repacked output array
         @param inputArrays         inputArrays[subarraySize][numberOfEntries]
      
         Output array = [ inputArrays[0][0]     inputArrays[1][0]    inputArrays[2][0] inputArrays[3][0] \n
                           inputArrays[0][1]     inputArrays[1][1]    inputArrays[2][1] inputArrays[3][1] \n
                           inputArrays[0][2]     inputArrays[1][2]    inputArrays[2][2] inputArrays[3][2] \n
       
                           ... \n
                           inputArrays[0][numberOfEntries] ... inputArrays[subarraySize-1][numberOfEntries] \n
                         ] \n
         where N = numberOfEntries \n
      
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArrayOfArrays( int numberOfEntries, int subarraySize,
                                      RealOpenMM* array, RealOpenMM** inputArrays );
      
      /**---------------------------------------------------------------------------------------
      
         Collapse 2D array into packed single array (Simbios)
         Example forces[3][N] -> array of size N containing RealOpenMM3 values 
      
         @param numberOfEntries     entries (no. atoms)
         @param iUnroll             iUnroll
         @param jUnroll             jUnroll
         @param arrays              arrays to be merged (dimension is [iUnroll][numberOfEntries/iUnroll]
         @param mergeArray          output array (if null, then allocated)
         @param log                 logging file descriptor
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM* collapseArrays( int numberOfEntries, int iUnroll, int jUnroll,
                                         RealOpenMM** arrays, RealOpenMM* mergeArray, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Merge 2 arrays based on sentinel value  (Simbios)
         Overflow array, if present, signals which entries had nonsentinel in both arrays
      
         Example   \n
                   array1 = { 1, 2, 3, s, s, s, 4, 5, s } \n
                   array2 = { s, s, 6, 8, 9, 3, s, 7, s } \n
             merge array  = { 1, 2, 3, 8, 9, 3, 4, 5, s } \n
          overflow array  = { s, s, 6, s, s, s, s, 7, s } \n
      
         @param numberOfEntries     entries (no. atoms)
         @param sentinelValue       sentinel value
         @param array1              first array
         @param array2              second array
         @param mergeArray          output merge array
         @param overflowArray       output overflow array  (if null, then ignored)
         @param log                 logging file descriptor
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int mergeArraysBasedSentinelValue( int numberOfEntries,
                                                RealOpenMM sentinelValue,
                                                RealOpenMM* array1, RealOpenMM* array2, 
                                                RealOpenMM* mergeArray, RealOpenMM* overflowArray,
                                                FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to store values from CPU loop in position seen in GPU output (Simbios)
      
         --------------------------------------------------------------------------------------- */
      
      static int storeInGpuFormat( int atomI, int atomJ, int debugAtomJ,
                                   int jUnroll, int numberOfStreams,
                                   RealOpenMM** debugStreams, RealOpenMM* RealOpenMMValues[2],
                                   RealOpenMM unsetDebugValue, bool useUpper, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Helper method to compare cpu and gpu computed arrays
      
         @param numberOfAtoms       entries (no. atoms)
         @param chunkSize           chunk size (usually 3 or 4)
         @param cpuArray            cpuArray[0-(chunkSize-1)][0,entries-1]
         @param gpuArray            gpuArray[chunkSize*entries]
         @param tolerance           check if relative difference is greater
                                    than this tolerance
         @param compareColumn       array of size [chunkSize] signalling whether
                                    column i is to be compared; it may be NULL
         @param absoluteMin         error if abs(cpu) + abs(gpu) > absoluteMin
                                    set negative to ignore this condition
         @param printOn             if not set, then no printing
         @param header              id header -- optional
         @param log                 logging file descriptor
      
         @return number of failed entries 
      
         --------------------------------------------------------------------------------------- */
      
      static int compareCpuGpuArrays( int numberOfAtoms, int chunkSize,
                                      RealOpenMM** cpuArray, RealOpenMM* gpuArray,
                                      RealOpenMM tolerance,
                                      int* compareColumn, RealOpenMM absoluteMin,
                                      int printOn, const char* header, FILE* log );

   // ---------------------------------------------------------------------------------------

};
   
// ---------------------------------------------------------------------------------------

#endif // __SimTKOpenMMGpuUtilities_H__
