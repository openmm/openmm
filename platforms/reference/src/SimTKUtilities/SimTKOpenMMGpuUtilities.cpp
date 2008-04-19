
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

// class of shared, static utility methods

#include "SimTKOpenMMGpuUtilities.h"
#include "SimTKOpenMMUtilities.h"

// fabs(), ...

#include <math.h>

//#define UseGromacsMalloc 1
//#ifdef UseGromacsMalloc
//extern "C" {
//#include "smalloc.h" 
//}
//#endif

/* ---------------------------------------------------------------------------------------

   Helper method to repack RealOpenMM arrays (Simbios)

   @param numberOfEntries     entries/sub-array
   @param subarraySize        number of subarrays
   @param array               array

   Input: array = [subArray_1 subArray_2 subArray_3 ... subArray_Stacked]
          where each subArray_i is 1 x numberOfEntries

   Output: array = [ subArray_1_1 subArray_2_1 subArray_3_1 ... subArray_Stacked_1
                     subArray_1_2 subArray_2_2 subArray_3_2 ... subArray_Stacked_2
                     subArray_1_3 subArray_2_3 subArray_3_3 ... subArray_Stacked_3
                     ...
                     subArray_1_N subArray_2_N subArray_3_N ... subArray_Stacked_N ]
   where N = numberOfEntries

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::repackArray1( int numberOfEntries, int subarraySize,
                                           RealOpenMM* array ){

   // ---------------------------------------------------------------------------------------

   RealOpenMM* tempArray;
   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::rePackArray";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfArray = sizeof( RealOpenMM )*numberOfEntries*subarraySize; 
   tempArray                = (RealOpenMM*) SimTKOpenMMUtilities::Xmalloc( "tempArray",
                                                                           __FILE__, __LINE__,
                                                                          sizeOfArray );

   memcpy( tempArray, array, sizeOfArray );

   int arrayIndex           = 0;
   for( int jj = 0; jj < subarraySize; jj++ ){
      for( int ii = 0; ii < numberOfEntries; ii++ ){
         array[arrayIndex++] = tempArray[jj*numberOfEntries+ii];
      }
   }

   SimTKOpenMMUtilities::Xfree( "tempArray", __FILE__, __LINE__, tempArray );

   return SimTKOpenMMCommon::DefaultReturn;
}      

/* ---------------------------------------------------------------------------------------

   Copy an array into a packed (e.g., RealOpenMM4 array) (Simbios)

   Example: copy Born radii into last slot of force array \n

   @param numberOfEntries     entries/sub-array (no. atoms)
   @param subarraySize        number of subarrays (4 for RealOpenMM4)
   @param fullArray           full array (force array in example)
   @param copySlot            index of slot to copied into (3 in example, since want Born radius in .w slot)
   @param arrayToCopy         array to copy (Born array)

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::copySubArray( int numberOfEntries, int subarraySize, RealOpenMM* fullArray,
                                           int copySlot, RealOpenMM* arrayToCopy ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::copySubArray";

   // ---------------------------------------------------------------------------------------

   int arrayIndex           = copySlot;
   for( int ii = 0; ii < numberOfEntries; ii++ ){
      fullArray[arrayIndex]  = arrayToCopy[ii];
      arrayIndex            += subarraySize;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to repack RealOpenMM arrays (Simbios)

   @param numberOfEntries     entries/sub-array
   @param subarraySize        number of subarrays
   @param array               array

   Input: array = [subArray_1 subArray_2 subArray_3 ... subArray_N] \n

          where each subArray_i is vector of length M=numberOfEntries \n
 
   Output: array = [ subArray_1_1 subArray_2_1 subArray_3_1 ... subArray_N_1 \n
                     subArray_1_2 subArray_2_2 subArray_3_2 ... subArray_N_2 \n
                     subArray_1_3 subArray_2_3 subArray_3_3 ... subArray_N_3 \n
                     ... \n
                     subArray_1_M subArray_2_M subArray_3_M ... subArray_M ] \n
   where N = numberOfEntries \n

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::repackArray( int numberOfEntries, int subarraySize, RealOpenMM* array ){

   // ---------------------------------------------------------------------------------------

   RealOpenMM* tempArray;
   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::rePackArray";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfArray = sizeof( RealOpenMM )*numberOfEntries*subarraySize; 
   tempArray                = (RealOpenMM*) SimTKOpenMMUtilities::Xmalloc( "tempArray", __FILE__, __LINE__,
                                                                           sizeOfArray );

   memcpy( tempArray, array, sizeOfArray );

   int arrayIndex           = 0;
   for( int ii = 0; ii < numberOfEntries; ii++ ){
      for( int jj = 0; jj < subarraySize; jj++ ){
         array[arrayIndex++] = tempArray[jj*numberOfEntries+ii];
      }
   }

   SimTKOpenMMUtilities::Xfree( "tempArray", __FILE__, __LINE__, tempArray );

   return SimTKOpenMMCommon::DefaultReturn;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to repack RealOpenMM arrays (Simbios)

   @param numberOfEntries     entries/sub-array
   @param subarraySize        number of subarrays
   @param array               repacked output array
   @param inputArrays         inputArrays[subarraySize][numberOfEntries]

   Output: array = [ inputArrays[0][0]     inputArrays[1][0]    inputArrays[2][0] inputArrays[3][0] \n
                     inputArrays[0][1]     inputArrays[1][1]    inputArrays[2][1] inputArrays[3][1] \n
                     inputArrays[0][2]     inputArrays[1][2]    inputArrays[2][2] inputArrays[3][2] \n
 
                     ... \n
                     inputArrays[0][numberOfEntries] ... inputArrays[subarraySize-1][numberOfEntries] \n
                   ] \n
   where N = numberOfEntries \n


   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::repackArrayOfArrays( int numberOfEntries, int subarraySize, 
                                                  RealOpenMM* array, RealOpenMM** inputArrays ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::repackArrayOfArrays";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfBlock = sizeof( RealOpenMM )*numberOfEntries;
   RealOpenMM* RealOpenMMPtr          = array;
   for( int ii = 0; ii < subarraySize; ii++ ){
      memcpy( RealOpenMMPtr, inputArrays[ii], sizeOfBlock );
      RealOpenMMPtr += numberOfEntries;
   }
   return repackArray( numberOfEntries, subarraySize, array );

}      

/* ---------------------------------------------------------------------------------------

   Collapse 2D array into packed single array (Simbios)
   Example: forces[3][N] -> array of size N containing RealOpenMM3 values 

   @param numberOfEntries     entries (no. atoms)
   @param iUnroll             iUnroll
   @param jUnroll             jUnroll
   @param arrays              arrays to be merged (dimension is [iUnroll][numberOfEntries/iUnroll]
   @param mergeArray          output array (if null, then allocated)
   @param log                 logging file descriptor

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

RealOpenMM* SimTKOpenMMGpuUtilities::collapseArrays( int numberOfEntries, int iUnroll, int jUnroll,
                                                     RealOpenMM** arrays, RealOpenMM* mergeArray,
                                                     FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName    = "\nSimTKOpenMMGpuUtilities::collapseArrays";
   static bool printOn              = false;

   // ---------------------------------------------------------------------------------------

   printOn = printOn && log != NULL;

   if( mergeArray == NULL ){
      unsigned int sizeOfArray = sizeof( RealOpenMM )*numberOfEntries*jUnroll; 
      mergeArray               = (RealOpenMM*) SimTKOpenMMUtilities::Xmalloc( "mergeArray", __FILE__, __LINE__,
                                                                              sizeOfArray );
   }

   int arrayIndex  = 0;
   int arrayOffset = 0;
   if( printOn ){
       (void) fprintf( log, "\n%s %d %d %d", methodName, numberOfEntries, iUnroll, jUnroll );
   }
   for( int ii = 0; ii < numberOfEntries/iUnroll; ii++ ){
      for( int kk = 0; kk < iUnroll; kk++ ){
         for( int jj = 0; jj < jUnroll; jj++ ){
            mergeArray[arrayIndex++] = arrays[kk][arrayOffset+jj];
         }

         if( printOn ){
            (void) fprintf( log, "\n%d %d [%.4f %.4f %.4f %.4f]", ii,kk,
                            arrays[kk][arrayOffset+0],
                            arrays[kk][arrayOffset+1],
                            arrays[kk][arrayOffset+2],
                            arrays[kk][arrayOffset+3] );
         }

      }
      arrayOffset += jUnroll;
   }

   return mergeArray;
}      

/* ---------------------------------------------------------------------------------------

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

int SimTKOpenMMGpuUtilities::mergeArraysBasedSentinelValue( int numberOfEntries,
                                                            RealOpenMM sentinelValue,
                                                            RealOpenMM* array1, RealOpenMM* array2, 
                                                            RealOpenMM* mergeArray, RealOpenMM* overflowArray,
                                                            FILE* log ){

   // ---------------------------------------------------------------------------------------

   int hits[3];
   RealOpenMM* arrays[4];

   const RealOpenMM tolerance             = (RealOpenMM) 0.00001;
   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::sentinelArraysBasedInitializedValue";

   // ---------------------------------------------------------------------------------------

   arrays[0] = array1;
   arrays[1] = array2;
   arrays[2] = mergeArray;
   arrays[3] = overflowArray;

   for( int ii = 0; ii < numberOfEntries; ii++ ){

      RealOpenMM value;
      RealOpenMM overflowValue = sentinelValue;

      for( int jj = 0; jj < 2; jj++ ){
         hits[jj] = fabs( arrays[jj][ii] - sentinelValue ) < tolerance ? 1 : 0;
      }

      // both missing
     
      if( hits[0] && hits[1] ){
         value = sentinelValue;   
      
      // single hit
 
      } else if( !hits[0] && hits[1] ){
         value = arrays[0][ii];   
      } else if( hits[0] && !hits[1] ){
         value = arrays[1][ii];   
      } else {

      // both present -- add to overflow array if available

         value         = arrays[0][ii];   
         overflowValue = arrays[1][ii];
      }
      arrays[2][ii] = value;
      if( arrays[3] != NULL ){
         arrays[3][ii] = overflowValue;
      }
   }

   return SimTKOpenMMCommon::DefaultReturn;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to store values from CPU loop in position seen in GPU output (Simbios)

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::storeInGpuFormat( int atomI, int atomJ, int debugAtomJ,
                                               int jUnroll, int numberOfStreams,
                                               RealOpenMM** debugStreams,
                                               RealOpenMM* RealOpenMMValues[2],
                                               RealOpenMM unsetDebugValue, bool useUpper,
                                               FILE* log ){

   // ---------------------------------------------------------------------------------------

   int diffIndices[2][2];
   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::storeInGpuFormat";

   // ---------------------------------------------------------------------------------------

   diffIndices[0][0]  = atomI - debugAtomJ;
   diffIndices[0][1]  = atomJ;
   diffIndices[1][0]  = atomJ - debugAtomJ;
   diffIndices[1][1]  = atomI;

   for( int ii = 0; ii < 2; ii++ ){
      if( diffIndices[ii][0] >= 0 && diffIndices[ii][0] < jUnroll ){
         for( int jj = 0; jj < numberOfStreams; jj++ ){
            RealOpenMM* RealOpenMMPtr  = &(debugStreams[jj][jUnroll*diffIndices[ii][1]]);
            RealOpenMMPtr             += diffIndices[ii][0];

/*
if( log && atomI < 10 && atomJ < 10 ){
(void) fprintf( log, "%s i=%d j=%d dbg=%d dff=%d crnt=%.4f new=%.4f ii=%d jj=%d", methodName,
                atomI, atomJ, debugAtomJ, diffIndices[ii][0], *RealOpenMMPtr,
                RealOpenMMValues[jj][diffIndices[ii][0]], ii, jj );
}
*/

/*
            if( fabs( *RealOpenMMPtr - unsetDebugValue ) < 1.0e-04 || 
                (atomI > atomJ && useUpper) ){
               *RealOpenMMPtr = RealOpenMMValues[jj][diffIndices[ii][0]];
            }
*/
            *RealOpenMMPtr = RealOpenMMValues[jj][diffIndices[ii][0]];
         }
      }
   }

   return SimTKOpenMMCommon::DefaultReturn;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to compare cpu and gpu computed arrays

   @param numberOfAtoms       entries (no. atoms)
   @param chunkSize           chunk size (usually 3 or 4)
   @param cpuArray            cpuArray[0-(chunkSize-1)][0,entries-1]
   @param gpuArray            gpuArray[chunkSize*entries]
   @param tolerance           check if relative difference is greater \n
                              than this tolerance
   @param compareColumn       array of size [chunkSize] signalling whether \n
                              column i is to be compared; it may be NULL
   @param absoluteMin         error if abs(cpu) + abs(gpu) > absoluteMin \n
                              set negative to ignore this condition
   @param printOn             if not set, then no printing
   @param header              id header -- optional
   @param log                 logging file descriptor

   @return number of failed entries 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGpuUtilities::compareCpuGpuArrays( int numberOfAtoms, int chunkSize,
                                                  RealOpenMM** cpuArray, RealOpenMM* gpuArray, RealOpenMM tolerance,
                                                  int* compareColumn, RealOpenMM absoluteMin,
                                                  int printOn, const char* header, FILE* log ){

   // ---------------------------------------------------------------------------------------

   char failedString[20];
   // static const char* methodName    = "\nSimTKOpenMMGpuUtilities::compareCpuGpuArrays";

   // ---------------------------------------------------------------------------------------

   printOn = printOn && log != NULL;

   // print header

   failedString[chunkSize] = '\0';

   if( printOn ){
      (void) fprintf( log, "\n" );
      if( header ){
         (void) fprintf( log, "%s", header );
      }
      (void) fprintf( log, " atoms=%d tol=%.3e absCutoff=%.3e printOn=%d",
                      numberOfAtoms, tolerance, absoluteMin, printOn );
   }

   // look for differences and print if flags set appropriately

   int returnFailed = 0;
   int offset       = 0;
   int printedOnce  = 0;

   for( int ii = 0; ii < numberOfAtoms; ii++ ){
         
      int localFailed = 0;
      for( int jj = 0; jj < chunkSize; jj++ ){
         if( compareColumn == NULL || compareColumn[jj] ){
            RealOpenMM   f1     = fabs( cpuArray[jj][ii] ) + fabs( gpuArray[offset+jj] );
            if( f1 > 0.0f ){
               RealOpenMM diff  = fabs( (cpuArray[jj][ii] - gpuArray[offset+jj]) )/f1;  
               if( diff > tolerance && f1 > absoluteMin ){ 
                  localFailed      = 1; 
                  returnFailed    += 1;
                  failedString[jj] = 'X'; 
               } else {
                  failedString[jj] = ' '; 
               }
            }
         } else {
            failedString[jj] = ' ';
         }
      }

      // print

      if( log && (printOn || localFailed) ){

         const char* numberFormat = "%9.4f ";
         if( !printedOnce ){
            (void) fprintf( log, "\n   CpuF                GpuF" );
            printedOnce = 1;
         }

         (void) fprintf( log, "\n%d [", ii );

         for( int jj = 0; jj < chunkSize; jj++ ){
            if( compareColumn == NULL || compareColumn[jj] ){
               (void) fprintf( log, numberFormat, cpuArray[jj][ii] );
            }
         }

         (void) fprintf( log, "] [" );

         for( int jj = 0; jj < chunkSize; jj++ ){
            if( compareColumn == NULL || compareColumn[jj] ){
               (void) fprintf( log, numberFormat, gpuArray[offset + jj] );
            }
         }

         (void) fprintf( log, "] [%s]", failedString ); 
         
      }
      offset += chunkSize;
   }

   // flush buffer

   if( log ){
      (void) fflush( log );
   }

   return returnFailed;
}      

