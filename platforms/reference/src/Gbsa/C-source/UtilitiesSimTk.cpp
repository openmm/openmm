
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

#include "UtilitiesSimTk.h"

// fabs(), ...

#include <math.h>

#define UseGromacsMalloc 1
#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given atom (Simbios)

   @param atomCoordinates     atom coordinates
   @param atomIndex           atom index to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least
                              numberOfAtoms
   @param log                 if set, then print error messages to log file

   @return distances

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::getDistanceSquaredFromSpecifiedAtom( const rvec* atomCoordinates, int atomIndex,
                                                         int numberOfAtoms, float* distances,
                                                         FILE* log ){

   // ---------------------------------------------------------------------------------------

   float atomXyz[3];
   // static const char* methodName    = "\nUtilitiesSimTk::getDistanceSquaredFromSpecifiedAtom";

   // ---------------------------------------------------------------------------------------

   for( int jj = 0; jj < 3; jj++ ){
      atomXyz[jj] = atomCoordinates[atomIndex][jj];
   }
      
   return getDistanceSquaredFromSpecifiedPoint( atomCoordinates, atomXyz,
                                                numberOfAtoms, distances, log );
}

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given point (Simbios)

   @param atomCoordinates     atom coordinates
   @param point               point to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least \n
                              numberOfAtoms
   @param log                 if set, then print error messages to log file

   @return distances

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::getDistanceSquaredFromSpecifiedPoint( const rvec* atomCoordinates, float* point,
                                                          int numberOfAtoms,
                                                          float* distances, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nUtilitiesSimTk::getDistanceSquaredFromSpecifiedPoint";

   // ---------------------------------------------------------------------------------------

   memset( distances, 0, sizeof( float )*numberOfAtoms );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      for( int jj = 0; jj < 3; jj++ ){
         float diff = (point[jj] - atomCoordinates[ii][jj]);
         distances[ii] += diff*diff;
      }
   }
      
   return 0;
}

/* ---------------------------------------------------------------------------------------

   Helper method to allocate float arrays (Simbios)

   @param bufferIndex         buffer index
   @param allocatedSz         array of allocated sizes
   @param bufferArray         array of allocated float arrays
   @param requestedSize       requested size
   @param dataAction          action flag: -1 = free memory \n
                                      1 = zero memory

   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::allocateFloatBufferArray( int bufferIndex, int* allocatedSz, float** bufferArray,
                                              int requestedSize, int dataAction ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nUtilitiesSimTk::allocateFloatBufferArray";

   // ---------------------------------------------------------------------------------------

   // clear data

   if( dataAction == -1 ){
      if( allocatedSz[bufferIndex] && bufferArray[bufferIndex] ){
         UtilitiesSimTk::Xfree( "bufferArray", __FILE__, __LINE__, bufferArray[bufferIndex] );
         allocatedSz[bufferIndex] = 0;
         bufferArray[bufferIndex] = NULL;
      }
   }
 
   // return if requested size is less than allcated

   if( allocatedSz[bufferIndex] > requestedSize ){
      return 0;
   }

   // free space if currently allocated

   if(  allocatedSz[bufferIndex] && bufferArray[bufferIndex] ){
      UtilitiesSimTk::Xfree( "bufferArray", __FILE__, __LINE__, bufferArray[bufferIndex] );
   }

   // allocate

   // bufferArray[bufferIndex] = (float*) malloc( requestedSize*sizeof( float ) );
   bufferArray[bufferIndex] = (float*) UtilitiesSimTk::Xmalloc( "bufferArray", __FILE__, __LINE__, requestedSize*sizeof( float ) );
   allocatedSz[bufferIndex] = requestedSize;

   // zero?

   if( dataAction == 1 ){
      memset( bufferArray[bufferIndex], 0, requestedSize*sizeof( float ) );
   }

   return 0;
}

/* ---------------------------------------------------------------------------------------

   Helper method to repack float arrays (Simbios)

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

   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::repackArray1( int numberOfEntries, int subarraySize, float* array ){

   // ---------------------------------------------------------------------------------------

   float* tempArray;
   // static const char* methodName    = "\nUtilitiesSimTk::rePackArray";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfArray = sizeof( float )*numberOfEntries*subarraySize; 
   tempArray = (float*) UtilitiesSimTk::Xmalloc( "tempArray", __FILE__, __LINE__, sizeOfArray );

   memcpy( tempArray, array, sizeOfArray );

   int arrayIndex           = 0;
   for( int jj = 0; jj < subarraySize; jj++ ){
      for( int ii = 0; ii < numberOfEntries; ii++ ){
         array[arrayIndex++] = tempArray[jj*numberOfEntries+ii];
      }
   }

   UtilitiesSimTk::Xfree( "tempArray", __FILE__, __LINE__, tempArray );

   return 0;
}      

/* ---------------------------------------------------------------------------------------

   Copy an array into a packed (e.g., float4 array) (Simbios)

   Example: copy Born radii into last slot of force array \n

   @param numberOfEntries     entries/sub-array (no. atoms)
   @param subarraySize        number of subarrays (4 for float4)
   @param fullArray           full array (force array in example)
   @param copySlot            index of slot to copied into (3 in example, since want Born radius in .w slot)
   @param arrayToCopy         array to copy (Born array)

   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::copySubArray( int numberOfEntries, int subarraySize, float* fullArray,
                                  int copySlot, float* arrayToCopy ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nUtilitiesSimTk::copySubArray";

   // ---------------------------------------------------------------------------------------

   int arrayIndex           = copySlot;
   for( int ii = 0; ii < numberOfEntries; ii++ ){
      fullArray[arrayIndex]  = arrayToCopy[ii];
      arrayIndex            += subarraySize;
   }

   return 0;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to repack float arrays (Simbios)

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

   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::repackArray( int numberOfEntries, int subarraySize, float* array ){

   // ---------------------------------------------------------------------------------------

   float* tempArray;
   // static const char* methodName    = "\nUtilitiesSimTk::rePackArray";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfArray = sizeof( float )*numberOfEntries*subarraySize; 
   tempArray                = (float*) UtilitiesSimTk::Xmalloc( "tempArray", __FILE__, __LINE__, sizeOfArray );

   memcpy( tempArray, array, sizeOfArray );

   int arrayIndex           = 0;
   for( int ii = 0; ii < numberOfEntries; ii++ ){
      for( int jj = 0; jj < subarraySize; jj++ ){
         array[arrayIndex++] = tempArray[jj*numberOfEntries+ii];
      }
   }

   UtilitiesSimTk::Xfree( "tempArray", __FILE__, __LINE__, tempArray );

   return 0;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to repack float arrays (Simbios)

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


   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::repackArrayOfArrays( int numberOfEntries, int subarraySize, 
                                         float* array, float** inputArrays ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nUtilitiesSimTk::repackArrayOfArrays";

   // ---------------------------------------------------------------------------------------

   unsigned int sizeOfBlock = sizeof( float )*numberOfEntries;
   float* floatPtr          = array;
   for( int ii = 0; ii < subarraySize; ii++ ){
      memcpy( floatPtr, inputArrays[ii], sizeOfBlock );
      floatPtr += numberOfEntries;
   }
   return repackArray( numberOfEntries, subarraySize, array );

}      

/* ---------------------------------------------------------------------------------------

   Collapse 2D array into packed single array (Simbios)
   Example: forces[3][N] -> array of size N containing float3 values 

   @param numberOfEntries     entries (no. atoms)
   @param iUnroll             iUnroll
   @param jUnroll             jUnroll
   @param arrays              arrays to be merged (dimension is [iUnroll][numberOfEntries/iUnroll]
   @param mergeArray          output array (if null, then allocated)
   @param log                 logging file descriptor

   @return 0 

   --------------------------------------------------------------------------------------- */

float* UtilitiesSimTk::collapseArrays( int numberOfEntries, int iUnroll, int jUnroll,
                                       float** arrays, float* mergeArray, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName    = "\nUtilitiesSimTk::collapseArrays";
   static bool printOn              = false;

   // ---------------------------------------------------------------------------------------

   printOn = printOn && log != NULL;

   if( mergeArray == NULL ){
      unsigned int sizeOfArray = sizeof( float )*numberOfEntries*jUnroll; 
      mergeArray                = (float*) UtilitiesSimTk::Xmalloc( "mergeArray", __FILE__, __LINE__, sizeOfArray );
   }

   int arrayIndex           = 0;
   int arrayOffset          = 0;
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

   @return 0 

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::mergeArraysBasedSentinelValue( int numberOfEntries,
                                                   float sentinelValue,
                                                   float* array1, float* array2, 
                                                   float* mergeArray, float* overflowArray,
                                                   FILE* log ){

   // ---------------------------------------------------------------------------------------

   int hits[3];
   float* arrays[4];

   const float tolerance             = 0.00001f;
   // static const char* methodName    = "\nUtilitiesSimTk::sentinelArraysBasedInitializedValue";

   // ---------------------------------------------------------------------------------------

   arrays[0] = array1;
   arrays[1] = array2;
   arrays[2] = mergeArray;
   arrays[3] = overflowArray;

   for( int ii = 0; ii < numberOfEntries; ii++ ){

      float value;
      float overflowValue = sentinelValue;

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

   return 0;
}      

/* ---------------------------------------------------------------------------------------

   Helper method to store values from CPU loop in position seen in GPU output (Simbios)

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::storeInGpuFormat( int atomI, int atomJ, int debugAtomJ,
                                      int jUnroll, int numberOfStreams,
                                      float** debugStreams, float* floatValues[2],
                                      float unsetDebugValue, bool useUpper, FILE* log ){

   // ---------------------------------------------------------------------------------------

   int diffIndices[2][2];
   // static const char* methodName    = "\nUtilitiesSimTk::storeInGpuFormat";

   // ---------------------------------------------------------------------------------------

   diffIndices[0][0]  = atomI - debugAtomJ;
   diffIndices[0][1]  = atomJ;
   diffIndices[1][0]  = atomJ - debugAtomJ;
   diffIndices[1][1]  = atomI;

   for( int ii = 0; ii < 2; ii++ ){
      if( diffIndices[ii][0] >= 0 && diffIndices[ii][0] < jUnroll ){
         for( int jj = 0; jj < numberOfStreams; jj++ ){
            float* floatPtr  = &(debugStreams[jj][jUnroll*diffIndices[ii][1]]);
            floatPtr        += diffIndices[ii][0];

/*
if( log && atomI < 10 && atomJ < 10 ){
(void) fprintf( log, "%s i=%d j=%d dbg=%d dff=%d crnt=%.4f new=%.4f ii=%d jj=%d", methodName,
                atomI, atomJ, debugAtomJ, diffIndices[ii][0], *floatPtr,
                floatValues[jj][diffIndices[ii][0]], ii, jj );
}
*/

/*
            if( fabs( *floatPtr - unsetDebugValue ) < 1.0e-04 || 
                (atomI > atomJ && useUpper) ){
               *floatPtr = floatValues[jj][diffIndices[ii][0]];
            }
*/
            *floatPtr = floatValues[jj][diffIndices[ii][0]];
         }
      }
   }

   return 0;
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

int UtilitiesSimTk::compareCpuGpuArrays( int numberOfAtoms, int chunkSize,
                                         float** cpuArray, float* gpuArray, float tolerance,
                                         bool* compareColumn, float absoluteMin,
                                         bool printOn, const char* header, FILE* log ){

   // ---------------------------------------------------------------------------------------

   char failedString[20];
   // static const char* methodName    = "\nUtilitiesSimTk::compareCpuGpuArrays";

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
            float   f1 = fabs( cpuArray[jj][ii] ) + fabs( gpuArray[offset+jj] );
            if( f1 > 0.0f ){
               float diff = fabs( (cpuArray[jj][ii] - gpuArray[offset+jj]) )/f1;  
               if( diff > tolerance && f1 > absoluteMin ){ 
                  localFailed   = 1; 
                  returnFailed += 1;
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

/* ---------------------------------------------------------------------------------------

   Print atom coordinates, ...

   @param numberAtoms         numberAtoms
   @param atomCoordinates     atomCoordinates (may be NULL)
   @param numberOf1Darrays    number of 1-d arrays (may be 0)
   @param oneDArrays          1-d arrays 
   @param idString            id string to be printed if set
   @param log                 print messages to log file

   @return 0

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::printCoordinateAnd1DArrays( int numberAtoms, rvec* atomCoordinates,
                                                int numberOf1Darrays,
                                                float** oneDArrays, const char* idString, FILE* log ){

   // ---------------------------------------------------------------------------------------

   if( log == NULL ){
      return 0;
   }

   if( idString ){
      (void) fprintf( log, "\n%s", idString );
   }
   for( int ii = 0; ii < numberAtoms; ii++ ){
      if( atomCoordinates != NULL ){
         (void) fprintf( log, "\n%d %12.4f %12.4f %12.4f",
                         ii + 1, atomCoordinates[ii][0], atomCoordinates[ii][1], atomCoordinates[ii][2] );
      } else {
         (void) fprintf( log, "\n" );
      }
      for( int jj = 0; jj < numberOf1Darrays; jj++ ){
         (void) fprintf( log, " %12.4f", oneDArrays[jj][ii] );
      }
   }

   (void) fflush(log);

   // ---------------------------------------------------------------------------------------
	
	return 0;	
}

/* ---------------------------------------------------------------------------------------

   Get atom name from top data struct

   @param atomIndex           atom index
   @param outputAtomName      output atom name
   @param top                 GMX t_topology struct

   @return 0

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::getAtomNameGivenAtomIndex( int atomIndex, char* outputAtomName,
                                               const t_topology* top ){

   // ---------------------------------------------------------------------------------------
   // ---------------------------------------------------------------------------------------

   char*** atomNames       = top->atoms.atomname;
   const char* atomName    =  (*(atomNames[atomIndex])) == NULL || strlen( (*(atomNames[atomIndex])) ) < 1 ||
                              strlen( (*(atomNames[atomIndex])) ) > 100 ? "NA" : (*(atomNames[atomIndex]));

   (void) strcpy( outputAtomName, atomName );

   return 0;

}

/* ---------------------------------------------------------------------------------------

   Get residue name from top data struct given atom index

   @param atomIndex           atom index
   @param top                 GMX t_topology struct
   @param outputResidueName   output residue name (assume enough memory has been allocated)
   @param outputResidueIndex  if not null, then *outputResidueIndex is residue index

   @return 0

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::getResidueNameGivenAtomIndex( int atomIndex, const t_topology* top,
                                                  char* outputResidueName, int* outputResidueIndex ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::getResidueNameGivenAtomIndex";

   // ---------------------------------------------------------------------------------------

   char*** residueNames = top->atoms.resname;
   int residueIndex     = top->atoms.atom[atomIndex].resnr;

   if( outputResidueIndex != NULL ){
      *outputResidueIndex = residueIndex;
   }

   const char* residueName =  (*(residueNames[residueIndex])) == NULL          ||
                              strlen( (*(residueNames[residueIndex])) ) < 1 ||
                              strlen( (*(residueNames[residueIndex])) ) > 100 ? "NA" : (*(residueNames[residueIndex]));

   (void) strcpy( outputResidueName, residueName );

   return 0;

   // ---------------------------------------------------------------------------------------
}

/* ---------------------------------------------------------------------------------------

   Get atom name from top data struct

   @param atomIndex           atom index
   @param top                 GMX t_topology struct
   @param buffer              output buffer (enough space should have been reserved) 
   @param maxAtoms            max number of atoms for this run (may change -- used mainly
                              to keep from reallocating cache array)
   @param tab                 tab spacing

   @return 0

   Atom info is cached in atomIdStrings to speed things up.
   The cache memory may be freed by calling the method w/
   atomIndex == -1

   @see freeArrayOfStrings()

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::getAtomIdStringGivenAtomIndex( int atomIndex, const t_topology* top,
                                                   char* buffer, 
                                                   int maxAtoms, unsigned int tab ){

   // ---------------------------------------------------------------------------------------

   static int maxAtomIndex         = -1;
   static char** atomIdStrings     = NULL;

   char atomName[32];
   char residueName[32];

   int residueIndex;

   // static const char* methodName = "\nUtilitiesSimTk::getAtomIdStringGivenAtomIndex";

   // ---------------------------------------------------------------------------------------

   // free cache memory if allocated

   if( atomIndex == -1 ){
      // (void) fprintf( stdout, "UtilitiesSimTk: called getAtomIdStringGivenAtomIndex to delete cached strings %d.", maxAtomIndex );
      if( maxAtomIndex > 0 && atomIdStrings != NULL ){
         freeArrayOfStrings( maxAtomIndex, atomIdStrings );
         atomIdStrings = NULL;
         maxAtomIndex  = -1;
      }         
      return 0;
   }         

   // allocate cache memory?

   if( maxAtoms > 0 && (maxAtomIndex == -1 || maxAtoms > maxAtomIndex) ){
      if(  maxAtoms > maxAtomIndex && maxAtomIndex > 0 ){
         freeArrayOfStrings( maxAtomIndex, atomIdStrings );
      }
      maxAtomIndex = maxAtoms + 1;

//    atomIdStrings = (char**) malloc( maxAtomIndex*sizeof( char* ) ); 
      atomIdStrings = (char**) UtilitiesSimTk::Xmalloc( "atomIdStrings", __FILE__, __LINE__, maxAtomIndex*sizeof( char* ) );
      memset( atomIdStrings, 0, maxAtomIndex*sizeof( char* ) );

   }

   // if id is cached, return it

   if( atomIndex < maxAtomIndex && atomIdStrings[atomIndex] != NULL ){ 
      (void) strcpy( buffer, atomIdStrings[atomIndex] );
      return 0;
   }

   // not cached -- assemble info

   getAtomNameGivenAtomIndex( atomIndex, atomName, top );
   getResidueNameGivenAtomIndex( atomIndex, top, residueName, &residueIndex );

   (void) sprintf( buffer, "%s_%d %s", residueName, residueIndex, atomName );

   // tab string

   if( tab > 0 && strlen( buffer ) < tab ){
      tabStringInPlace( buffer, tab );
   }

   // cache info if atomIdStrings array is allocated

   if( atomIndex < maxAtomIndex && atomIdStrings ){

      if( atomIdStrings[atomIndex] != NULL ){
         UtilitiesSimTk::Xfree( "atomIdStrings", __FILE__, __LINE__,  atomIdStrings[atomIndex] );
      }
      atomIdStrings[atomIndex] = (char*) UtilitiesSimTk::Xmalloc( "atomIdStrings[atomIndex]", __FILE__, __LINE__, (strlen( buffer ) + 1)*sizeof( char ) );
      (void) strcpy( atomIdStrings[atomIndex], buffer );
   }

   return 0;

}

/* ---------------------------------------------------------------------------------------

   Free array of strings

   @param arraySz             atom index
   @param arrayOfStrings      array of strings

   @return 0

   Used for freeing an array of strings

   @see getAtomIdStringGivenAtomIndex()

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::freeArrayOfStrings( int arraySz, char** arrayOfStrings ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::freeArrayOfStrings";

   // ---------------------------------------------------------------------------------------

   // free memory if allocated

   for( int ii = 0; ii < arraySz; ii++ ){
      if( arrayOfStrings[ii] ){
         UtilitiesSimTk::Xfree( "atomIdStrings", __FILE__, __LINE__,  arrayOfStrings[ii] );
      }
   }
   UtilitiesSimTk::Xfree( "arrayOfStrings", __FILE__, __LINE__,  arrayOfStrings );

   return 0;
}

/* ---------------------------------------------------------------------------------------

   Tab string in place

   @param string   string to tab; assume string is of at least length=tab + 1
   @param tab      tab spacing length

   @return 0

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::tabStringInPlace( char* string, int tab ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::tabStringInPlace";

   // ---------------------------------------------------------------------------------------

   for( int ii = strlen( string ); ii < tab; ii++ ){
      string[ii] = ' ';
   }
   string[tab] = '\0';

   return 0;
}

/* ---------------------------------------------------------------------------------------

   Write debug fields (Simbios)

   @param numberOfFields       number of fields to print
   @param fields               fields
   @param numberOfStringFields number of string fields to print
   @param stringFields         string fields
   @param comment              comment (optinal -- ignored if NULL)
   @param debugFileName        output debug file name
   @param action               0 open file and return w/o printing \n
                               1 open file and print \n
                               2 close file (no print) \n
   @param log                  if set, then print error messages to log file

   @return debugFile unless file coud not be opened (or other errors )
   or file is closed -- for these cases return NULL

   stringFields printed after float fields

   --------------------------------------------------------------------------------------- */

FILE* UtilitiesSimTk::writeDebugFile( int numberOfFields, float* fields,
                                      int numberOfStringFields,
                                      char* stringFields[MAX_DEBUG_FIELDS],
                                      char* comment, const char* debugFileName, int action,
                                      FILE* debugFile, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nUtilitiesSimTk::writeDebugFile";

   // ---------------------------------------------------------------------------------------

   if( debugFileName != NULL && (action <= WriteDebugFile || debugFile == NULL) ){

   // open file

      debugFile = fopen( debugFileName, "w" );
      if( debugFile != NULL ){
         if( log != NULL ){
            (void) fprintf( log, "%s opened file=<%s>.", methodName, debugFileName );
            (void) fflush( log );
         }
      } else {
         if( log != NULL ){
            (void) fprintf( log, "%s could not open file=<%s> -- abort output.",
                            methodName, debugFileName );
            (void) fflush( log );
         }
         return NULL;
      }

      if( action == OpenDebugFile ){
         return debugFile;
      }

   } else if( action == CloseDebugFile ){

   // close file

      if( debugFile ){
         if( log != NULL ){
            (void) fprintf( log, "%s closed debug file=<%s>.", 
                            methodName, debugFileName == NULL ? "??" : debugFileName );
            (void) fflush( log );
         }
         (void) fclose( debugFile );
      }
      return NULL;
   }   

   if( comment != NULL ){
      (void) fprintf( debugFile, "%s", comment );
   }

   if( numberOfFields > 0 || numberOfStringFields > 0 ){
      (void) fprintf( debugFile, "\n" );
      if( numberOfFields > 0 || fields != NULL ){
         for( int ii = 0; ii < numberOfFields; ii++ ){
            (void) fprintf( debugFile, "%.5e ", fields[ii] );
         }
      }
      if( numberOfStringFields > 0 && stringFields != NULL ){
         for( int ii = 0; ii < numberOfStringFields; ii++ ){
            (void) fprintf( debugFile, "%s ", stringFields[ii] );
         }
      }
   }

   (void) fflush( debugFile );

   return debugFile;
}

/* ---------------------------------------------------------------------------------------

   Allocate 2D float array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialize           if true, then initialize array
   @param initialValue         intitial value
   @param log                  if set, then print error messages to log file

   @return allocated array

   --------------------------------------------------------------------------------------- */

float** UtilitiesSimTk::allocate2DFloatArray( int iSize, int jSize, float** array2D, 
                                              bool initialize, float initialValue, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::allocate2DArray";

   // ---------------------------------------------------------------------------------------

   if( array2D == NULL ){

      array2D      = (float**) UtilitiesSimTk::Xmalloc( "array2D", __FILE__, __LINE__, iSize*sizeof( float* ) );
      float* block = (float*)  UtilitiesSimTk::Xmalloc( "block",   __FILE__, __LINE__, jSize*iSize*sizeof( float ) );

      for( int ii = 0; ii < iSize; ii++ ){
         array2D[ii]  = block;
         block       += jSize;
      }    
   }

   if( initialize ){
      initialize2DFloatArray( iSize, jSize, array2D, initialValue, log );
   }

   return array2D;
}

/* ---------------------------------------------------------------------------------------

   Initialize 2D float array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialValue         intitial value
   @param log                  if set, then print error messages to log file

   @return array

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::initialize2DFloatArray( int iSize, int jSize,
                                            float** array2D, float initialValue, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::initialize2DFloatArray";

   // ---------------------------------------------------------------------------------------

   bool useMemset;
   bool useMemsetSingleBlock;

   if( initialValue == 0.0f ){
      useMemset = true;
      if( jSize > 1 && (array2D[0] + jSize) == array2D[1] ){  
         useMemsetSingleBlock = true;
      } else {
         useMemsetSingleBlock = false;
      }

   } else {
      useMemset = false;
   }

   if( useMemset ){
      if( useMemsetSingleBlock ){
         memset( array2D[0], 0, iSize*jSize*sizeof( float ) );
      } else {
         for( int ii = 0; ii < iSize; ii++ ){
            memset( array2D[ii], 0, jSize*sizeof( float ) );
         }
      }
   } else {
      for( int ii = 0; ii < iSize; ii++ ){
         for( int jj = 0; jj < jSize; jj++ ){
            array2D[ii][jj] = initialValue;
         }
      }
   }

   return 0;
}

/* ---------------------------------------------------------------------------------------
      
   Malloc memory of size bytesToAllocate and zero
      
   @param bytesToAllocate      bytes to allocate
      
   @return ptr to allocated memory; NULL if bytesToAllocate <= 0
      
   --------------------------------------------------------------------------------------- */
          
char* UtilitiesSimTk::allocateAndZero( unsigned int bytesToAllocate ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::allocateAndZero";

   // ---------------------------------------------------------------------------------------

   if( bytesToAllocate <= 0 ){
      return NULL;
   }

   char* ptrToMemory = (char*) UtilitiesSimTk::Xmalloc( "ptrToMemory", __FILE__, __LINE__, bytesToAllocate*sizeof( char ) );
   memset( ptrToMemory, 0, bytesToAllocate );

   return ptrToMemory;
}

/* ---------------------------------------------------------------------------------------

   Normalize 3-vector -- helper method

   @param vector 3-vector to normalize

   --------------------------------------------------------------------------------------- */
     
void UtilitiesSimTk::normalizeVector3( Real* vector ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::normalizeVector3";

   // ---------------------------------------------------------------------------------------

   Real sum   = vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
   sum        = sum > 0.0 ? (Real) (1.0/sqrt( (double) sum )) : (Real) 0.0;

   vector[0] *= sum;
   vector[1] *= sum;
   vector[2] *= sum;

   return;
}

/* ---------------------------------------------------------------------------------------

   Remove 3-vector -- helper method

   @param vectorToRemove      vector to remove
   @param vector              vector to from which 'vectorToRemove' is to be removed
                              vector is normalized after the component is subtracted out

   --------------------------------------------------------------------------------------- */
     
void UtilitiesSimTk::removeVector3( Real* vectorToRemove, Real* vector ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::removeVector3";

   // ---------------------------------------------------------------------------------------

   Real dot   = vectorToRemove[0]*vector[0] + vectorToRemove[1]*vector[1] + vectorToRemove[2]*vector[2];

   vector[0] -= dot*vectorToRemove[0];
   vector[1] -= dot*vectorToRemove[1];
   vector[2] -= dot*vectorToRemove[2];

   normalizeVector3( vector );
}

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector  -- helper method

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
void UtilitiesSimTk::crossProductVector3( Real* vectorX, Real* vectorY, Real* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::crossProductVector3";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

   return;
}

/**---------------------------------------------------------------------------------------

   Compute matrix product of 3x3 matrix and 3-vector and place in 3rd vector  -- helper method

   vectorZ = matrixX . vectorY

   @param matrixX             matrixX
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
    
void UtilitiesSimTk::matrixProductVector3( Real* matrixX, Real* vectorY, Real* vectorZ ){
     
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::matrixProductVector3";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = matrixX[0]*vectorY[0] + matrixX[3]*vectorY[1] + matrixX[6]*vectorY[2];
   vectorZ[1]  = matrixX[1]*vectorY[0] + matrixX[4]*vectorY[1] + matrixX[7]*vectorY[2];
   vectorZ[2]  = matrixX[2]*vectorY[0] + matrixX[5]*vectorY[1] + matrixX[8]*vectorY[2];

   return;
}

/**---------------------------------------------------------------------------------------

   Compute cross product between two 3x3 matrices

   @param vectorZ = matrixX x matrixY

   @param matrixX             matrixX
   @param matrixY             matrixY
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
    
void UtilitiesSimTk::matrixCrossProductMatrix3( Real* matrixX, Real* matrixY, Real* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const int indices[3][2] = { { 3, 6 }, { 6, 0 }, { 0 , 3 } };
   // static const char* methodName = "\nUtilitiesSimTk::matrixCrossProductMatrix3";
   Real* xPtr[3];
   Real* yPtr[3];

   // ---------------------------------------------------------------------------------------

   xPtr[0] = matrixX;
   xPtr[1] = matrixX + 3;
   xPtr[2] = matrixX + 6;

   yPtr[0] = matrixY;
   yPtr[1] = matrixY + 3;
   yPtr[2] = matrixY + 6;

   vectorZ[0] = DOT3( xPtr[1], yPtr[2] ) - DOT3( xPtr[2], yPtr[1] );
   vectorZ[1] = DOT3( xPtr[2], yPtr[0] ) - DOT3( xPtr[0], yPtr[2] );
   vectorZ[2] = DOT3( xPtr[0], yPtr[1] ) - DOT3( xPtr[1], yPtr[0] );

   return;
}

/* ---------------------------------------------------------------------------------------

   Centralized malloc/new

   @param name                ptr name
   @param fileName            file name
   @param line                file line no.
   @param file line           size in bytes to be allocated

   @return ptr to allocated object

   --------------------------------------------------------------------------------------- */
     
void* UtilitiesSimTk::Xmalloc( char* name, char* fileName, int line, unsigned int size ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::Xmalloc";

   // ---------------------------------------------------------------------------------------

#ifdef UseGromacsMalloc
   return save_malloc( name, fileName, line, size );
#else
   return (void*) new char[size];
#endif

}

/* ---------------------------------------------------------------------------------------

   Centralized free/delete

   @param name                ptr name
   @param fileName            file name
   @param line                file line no.
   @param ptr                 ptr to be freed

   --------------------------------------------------------------------------------------- */
     
void UtilitiesSimTk::Xfree( char* name, char* fileName, int line, void* ptr ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::Xfree";

   // ---------------------------------------------------------------------------------------

#ifdef UseGromacsMalloc
   return save_free( name, fileName, line, ptr );
#else
   delete ptr;
   return;
#endif
}

/* ---------------------------------------------------------------------------------------
      
   Format array of reals
      
   @param message             input string stream
   @param realArray           array of Reals
   @param numberOfFields      number of fields (optional - defaults to 3)
   @param factor					scale factor
      
   @return 0

--------------------------------------------------------------------------------------- */
          
int UtilitiesSimTk::formatRealStringStream( std::stringstream& message, const Real* realArray,
		                                      int numberOfFields, Real factor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nUtilitiesSimTk::Xfree";

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < numberOfFields; ii++ ){ 
      message << factor*realArray[ii] << " ";
   }
   return 0;

}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::tokenizeString( char* lineBuffer, StringVector& tokenArray,
                                    const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::tokenizeString";

   // ---------------------------------------------------------------------------------------

// (void) fprintf( stdout, "\nIn UtilitiesSimTk::tokenizeString <%s>", lineBuffer );
// (void) fflush( stdout );

   char *ptr_c = NULL;

   for( ; (ptr_c = UtilitiesSimTk::strsep( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse lines

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return token

   --------------------------------------------------------------------------------------- */

char* UtilitiesSimTk::strsep( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::strsep"

	char *s;
	const char *spanp;
	int c, sc;
	char *token;

   // ---------------------------------------------------------------------------------------

	s = *lineBuffer;
	if( s == NULL ){
		return NULL;
   }

   for( token = s;; ){
      c = *s++;
		spanp = delimiter;
		do {
		   if( (sc = *spanp++) == c ){
			   if( c == 0 ){
				   s = NULL;
				} else {
					s[-1] = 0;
            }
				*lineBuffer = s;
				return token;
			}
		} while( sc != 0 );
	}

   return NULL;
}

/**---------------------------------------------------------------------------------------

   Local version of strncasecmp (missing in Windows) (static method) (Simbios)

   @param string1                 first string
   @param string2                 second string
   @param matchLength             match length

   @return 0

   --------------------------------------------------------------------------------------- */

int UtilitiesSimTk::localStrncasecmp( const char *string1, const char *string2, int matchLength ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::localStrncasecmp"

   char ch1,ch2;

   // ---------------------------------------------------------------------------------------

   if( matchLength == 0 ){ 
      return 0;
   }

   do {   
      ch1 = toupper(*(string1++));
      ch2 = toupper(*(string2++));
      if( ch1 != ch2 )return (ch1-ch2);
      matchLength--;
   } while( ch1 && matchLength ); 

   return 0;  

}

/**---------------------------------------------------------------------------------------

   Check that string is valid integer

   @param stringToCheck string to check

   @return true if string is a valid integer

   --------------------------------------------------------------------------------------- */

bool UtilitiesSimTk::isValidInteger( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::isValidInteger";

   // ---------------------------------------------------------------------------------------

   int ii;

   return checkString<int>(ii, stringToCheck, std::dec );
}

/**---------------------------------------------------------------------------------------

   Check that string is valid Real

   @param stringToCheck string to check

   @return true if string is a valid Real

   --------------------------------------------------------------------------------------- */

bool UtilitiesSimTk::isValidReal( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::isValidReal";

   // ---------------------------------------------------------------------------------------

   Real ii;

   return checkString<Real>(ii, stringToCheck, std::dec );
}

/**---------------------------------------------------------------------------------------

   Return lower case copy of string

   @param string                  string

   @return lower cased string

   --------------------------------------------------------------------------------------- */

//int UtilitiesSimTk::toLower( std::string& string ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nUtilitiesSimTk::toLower"

   // ---------------------------------------------------------------------------------------

	// transform string to lower case
	    
	// std::transform( string.begin(), string.end(), string.begin(), (int(*)(int)) std::tolower);
	//return 0;
		
//}
