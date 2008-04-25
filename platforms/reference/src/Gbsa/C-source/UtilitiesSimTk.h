
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

#ifndef __UtilitiesSimTk_H_
#define __UtilitiesSimTk_H_

// ---------------------------------------------------------------------------------------

// reserve the option to change between float and double

#include "RealTypeSimTk.h" 
#include "CommonSimTk.h" 

// ---------------------------------------------------------------------------------------

// class of shared, static utility methods

#include "CommonSimTk.h" 

#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>

// #include <string>

// ---------------------------------------------------------------------------------------

// Gromac's types

#include "typedefs.h" 

// ---------------------------------------------------------------------------------------

// reserve the option to change between float and double

#include "RealTypeSimTk.h" 

// degrees to radians conversion factor

#define DEG2RAD 0.017453292f

// ---------------------------------------------------------------------------------------

#define ATOM_ID_STRING_TAB 12
#define OpenDebugFile       0
#define WriteDebugFile      1
#define CloseDebugFile      2
#define MAX_DEBUG_FIELDS   20

// template is used to check if a string is integer, real, ...

template <class T>
bool checkString( T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&) ){
   std::istringstream iss(s); 
   return !(iss >> f >> t).fail();
}

/**---------------------------------------------------------------------------------------

   Class of static methods to be shared
   Most methods are standalone 'utility' methods

--------------------------------------------------------------------------------------- */

class UtilitiesSimTk {

   public:

      // dummy constructor/destructor

       UtilitiesSimTk(){};
      ~UtilitiesSimTk(){};

      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given atom (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param atomIndex           atom index to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
         @param log                 if set, then print error messages to log file
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedAtom( const rvec* atomCoordinates, int atomIndex,
                                                      int numberOfAtoms, float* distances, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given point (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param point               point to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
         @param log                 if set, then print error messages to log file
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedPoint( const rvec* atomCoordinates, float* point, 
                                                       int numberOfAtoms, float* distances, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Helper method to allocate float arrays (Simbios)
      
         @param bufferIndex         buffer index
         @param allocatedSz         array of allocated sizes
         @param bufferArray         array of allocated float arrays
         @param requestedSize       requested size
         @param dataAction          action flag: -1 = free memory \n
                                                  1 = zero memory
      
         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int allocateFloatBufferArray( int bufferIndex, int* allocatedSz, float** bufferArray,
                                           int requestedSize, int dataAction );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to repack float arrays (Simbios)
      
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

         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArray1( int numberOfEntries, int subarraySize, float* array );
      
      /**---------------------------------------------------------------------------------------
      
         Copy an array into a packed (e.g., float4 array) (Simbios)
      
         Example copy Born radii into last slot of force array
      
         @param numberOfEntries     entries/sub-array (no. atoms)
         @param subarraySize        number of subarrays (4 for float4)
         @param fullArray           full array (force array in example)
         @param copySlot            index of slot to copied into (3 in example, since want Born radius in .w slot)
         @param arrayToCopy         array to copy (Born array)
      
         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int copySubArray( int numberOfEntries, int subarraySize, float* fullArray,
                               int copySlot, float* arrayToCopy );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to repack float arrays (Simbios)
      
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

         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArray( int numberOfEntries, int subarraySize, float* array );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to repack float arrays (Simbios)
      
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
      
      
         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int repackArrayOfArrays( int numberOfEntries, int subarraySize,
                                      float* array, float** inputArrays );
      
      /**---------------------------------------------------------------------------------------
      
         Collapse 2D array into packed single array (Simbios)
         Example forces[3][N] -> array of size N containing float3 values 
      
         @param numberOfEntries     entries (no. atoms)
         @param iUnroll             iUnroll
         @param jUnroll             jUnroll
         @param arrays              arrays to be merged (dimension is [iUnroll][numberOfEntries/iUnroll]
         @param mergeArray          output array (if null, then allocated)
         @param log                 logging file descriptor
      
         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static float* collapseArrays( int numberOfEntries, int iUnroll, int jUnroll,
                                    float** arrays, float* mergeArray, FILE* log );
      
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
      
         @return 0 
      
         --------------------------------------------------------------------------------------- */
      
      static int mergeArraysBasedSentinelValue( int numberOfEntries,
                                                float sentinelValue,
                                                float* array1, float* array2, 
                                                float* mergeArray, float* overflowArray, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Helper method to store values from CPU loop in position seen in GPU output (Simbios)
      
         --------------------------------------------------------------------------------------- */
      
      static int storeInGpuFormat( int atomI, int atomJ, int debugAtomJ,
                                   int jUnroll, int numberOfStreams,
                                   float** debugStreams, float* floatValues[2],
                                   float unsetDebugValue, bool useUpper, FILE* log );

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
                                      float** cpuArray, float* gpuArray, float tolerance,
                                      bool* compareColumn, float absoluteMin,
                                      bool printOn, const char* header, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Print atom coordinates, ...
      
         @param numberAtoms         numberAtoms
         @param atomCoordinates     atomCoordinates (may be NULL)
         @param numberOf1Darrays     number of 1-d arrays (may be 0)
         @param oneDArrays          1-d arrays 
         @param idString            id string to be printed if set
         @param log                 print messages to log file
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      static int printCoordinateAnd1DArrays( int numberAtoms, rvec* atomCoordinates,
                                             int numberOf1Darrays, float** oneDArrays,
                                             const char* idString, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Get atom name from top data struct
      
         @param atomIndex           atom index
         @param outputAtomName      output atom name
         @param top                 GMX t_topology struct
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      static int getAtomNameGivenAtomIndex( int atomIndex, char* outputAtomName, const t_topology* top );

      /**---------------------------------------------------------------------------------------
      
         Get residue name from top data struct given atom index
      
         @param atomIndex           atom index
         @param top                 GMX t_topology struct
         @param outputResidueName   output residue name
         @param outputResidueIndex  if not null, then *outputResidueIndex is residue index
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      static int getResidueNameGivenAtomIndex( int atomIndex, const t_topology* top,
                                               char* outputResidueName, int* outputResidueIndex );

      /**---------------------------------------------------------------------------------------
      
         Get atom name from top data struct
      
         @param atomIndex           atom index
         @param top                 GMX t_topology struct
         @param buffer              output buffer (enough space should have been reserved) 
         @param maxAtoms            max number of atoms for this run (may change -- used mainly
                                    to keep from reallocating cache array)
         @param tab                 tab spacing

         @return 0

         --------------------------------------------------------------------------------------- */
      
      static int getAtomIdStringGivenAtomIndex( int atomIndex, const t_topology* top, char* buffer,
                                                int maxAtoms, unsigned int tab );

      /**---------------------------------------------------------------------------------------
      
         Free array of strings
      
         @param arraySz             atom index
         @param arrayOfStrings      array of strings
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      static int freeArrayOfStrings( int arraySz, char** arrayOfStrings );

      /**---------------------------------------------------------------------------------------
      
         Tab string in place
      
         @param string              string to tab; assume string is of at least length=tab + 1
         @param tab                 tab length
      
         --------------------------------------------------------------------------------------- */
      
      static int tabStringInPlace( char* string, int tab  );
      
      /**---------------------------------------------------------------------------------------
      
         Write debug fields (Simbios)
      
         @param numberOfFields       number of fields to print
         @param fields               fields
         @param numberOfStringFields number of string fields to print
         @param stringFields         string fields
         @param comment              comment (optinal -- ignored if NULL)
         @param debugFileName        output debug file name
         @param action               0 open file and @return w/o printing
                                     1 open file and print
                                     2 close file (no print)
         @param debugFile            debug file reference
         @param log                  if set, then print error messages to log file
      
         @return debugFile unless file is closed
      
         stringFields printed after float fields
      
         --------------------------------------------------------------------------------------- */
      
      static FILE* writeDebugFile( int numberOfFields, float* fields,
                                   int numberOfStringFields, char* stringFields[MAX_DEBUG_FIELDS],
                                   char* comment, const char* debugFileName, int action,
                                   FILE* debugFile, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Allocate 2D float array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialize           if true, then initialize array
         @param initialValue         intitial value
         @param log                  if set, then print error messages to log file
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static float** allocate2DFloatArray( int iSize, int jSize,
                                           float** array2D, bool initialize, float initialValue, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Initialize 2D float array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialValue         intitial value
         @param log                  if set, then print error messages to log file
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static int initialize2DFloatArray( int iSize, int jSize,
                                         float** array2D, float initialValue, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Malloc memory of size bytesToAllocate and zero
      
         @param bytesToAllocate      bytes to allocate
      
         @return ptr to allocated memory; NULL if bytesToAllocate <= 0
      
         --------------------------------------------------------------------------------------- */
      
      static char* allocateAndZero( unsigned int bytesToAllocate );

      /**---------------------------------------------------------------------------------------
      
         Normalize 3-vector -- helper method
      
         @param vector              vector to normalize
      
         --------------------------------------------------------------------------------------- */
            
      static void normalizeVector3( Real* vector );
      
      /**---------------------------------------------------------------------------------------
      
         Remove 3-vector -- helper method
      
         @param vectorToRemove      vector to remove
         @param vector              vector to from which 'vectorToRemove' is to be removed \n
                                    vector is normalized after the component is subtracted out
      
         --------------------------------------------------------------------------------------- */
      
      static void removeVector3( Real* vectorToRemove, Real* vector );
      
      /**---------------------------------------------------------------------------------------
      
         Compute cross product of two 3-vectors and place in 3rd vector  -- helper method
      
         @param vectorZ = vectorX x vectorY
      
         @param vectorX             x-vector
         @param vectorY             y-vector
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void crossProductVector3( Real* vectorX, Real* vectorY, Real* vectorZ );

      /**---------------------------------------------------------------------------------------
      
         Compute matrix product of 3x3 matrix and 3-vector and place in 3rd vector  -- helper method
      
         @param vectorZ = matrixX . vectorY
      
         @param matrixX             matrixX
         @param vectorY             y-vector
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void matrixProductVector3( Real* matrixX, Real* vectorY, Real* vectorZ );

      /**---------------------------------------------------------------------------------------
      
         Compute cross product between two 3x3 matrices
      
         @param vectorZ = matrixX . matrixY
      
         @param matrixX             matrixX
         @param matrixY             matrixY
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void matrixCrossProductMatrix3( Real* matrixX, Real* matrixY, Real* vectorZ );

      /* ---------------------------------------------------------------------------------------
      
         Centralized malloc/new
      
         @param name                ptr name
         @param fileName            file name
         @param line                file line no.
         @param file line           size in bytes to be allocated
      
         @return ptr to allocated object
      
         --------------------------------------------------------------------------------------- */
               
      static void* Xmalloc( char* name, char* fileName, int line, unsigned int size );
      
      /* ---------------------------------------------------------------------------------------
      
         Centralized free/delete
      
         @param name                ptr name
         @param fileName            file name
         @param line                file line no.
         @param ptr                 ptr to be freed
      
         --------------------------------------------------------------------------------------- */
      
      static void Xfree( char* name, char* fileName, int line, void* ptr );
      
      /* ---------------------------------------------------------------------------------------
      
         Format array of reals
      
         @param message             input string stream
         @param realArray				array of Reals
         @param numberOfFields      number of fields (optional - defaults to 3)
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      static int formatRealStringStream( std::stringstream& message, const Real* realArray, 
                                         int numberOfFields = 3, Real factor = (Real) 1.0f );
      
      /**---------------------------------------------------------------------------------------
      
         Tokenize a string (static method) (Simbios)
      
         @param lineBuffer           string to tokenize
         @param tokenArray           upon return vectory of tokens
         @param delimiter            token delimter
      
         @return number of args
      
         --------------------------------------------------------------------------------------- */
      
      static int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter = "\t\n " );
      
      /**---------------------------------------------------------------------------------------
      
         Replacement of sorts for strtok() (static method) (Simbios)
         Used to parse lines
      
         @param lineBuffer           string to tokenize
         @param delimiter            token delimter
      
         @return token
      
         --------------------------------------------------------------------------------------- */
      
      static char* strsep( char** lineBuffer, const char* delimiter );
      
      /**---------------------------------------------------------------------------------------
      
         Local version of strncasecmp (missing in Windows) (static method) (Simbios)
      
         @param string1                 first string
         @param string2                 second string
         @param matchLength             match length
         
         @return 0
         
         --------------------------------------------------------------------------------------- */
      
      static int localStrncasecmp( const char *string1, const char *string2, int matchLength );

      /**---------------------------------------------------------------------------------------
      
         Check that string is valid integer
      
         @param stringToCheck string to check
      
         @return true if string is a valid integer
      
         --------------------------------------------------------------------------------------- */
      
      static bool isValidInteger( std::string stringToCheck );

      /**---------------------------------------------------------------------------------------
      
         Check that string is valid Real
      
         @param stringToCheck string to check
      
         @return true if string is a valid Real
      
         --------------------------------------------------------------------------------------- */
      
      static bool isValidReal( std::string stringToCheck );

   // ---------------------------------------------------------------------------------------

};
   
// ---------------------------------------------------------------------------------------

#endif // __UtilitiesSimTk_H__
