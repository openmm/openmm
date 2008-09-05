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

/**
 * This tests the Brook stream.
 */

// #define ASSERT_EQUAL_RELTOL(expected, found, tol) {double _scale_ = std::fabs(expected) > 1.0 ? std::fabs(expected) : 1.0; if (std::fabs((expected)-(found))/_scale_ > (tol)) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

//#define ASSERT_EQUAL_RELVEC(expected, found, tol) {ASSERT_EQUAL_RELTOL((expected)[0], (found)[0], (tol)); ASSERT_EQUAL_RELTOL((expected)[1], (found)[1], (tol)); ASSERT_EQUAL_RELTOL((expected)[2], (found)[2], (tol));};

#include "../../../tests/AssertionUtilities.h"
//#include "OpenMMContext.h"
#include "BrookPlatform.h"
#include "ReferencePlatform.h"
#include "BrookStreamFactory.h"
#include "BrookBonded.h"
#include "BrookNonBonded.h"

#include "OpenMMContext.h"
#include "StandardMMForceField.h"
#include "System.h"
#include "LangevinIntegrator.h"
#include "BrookRandomNumberGenerator.h"
#include "BrookShakeAlgorithm.h"
#include "BrookStochasticDynamics.h"
#include "../src/sfmt/SFMT.h"
#include <iostream>
#include <vector>
#include <fstream>

typedef std::vector<std::string> StringVector;

using namespace OpenMM;
using namespace std;

const double TOL = 5e-4;

void testWriteRead( void ){

/*
   static const int ArraySz = 3000;
   
   BrookPlatform platform;

   // create and initialize arrays

   float* array      = new float[ArraySz];
   float* saveArray  = new float[ArraySz];
   for( int ii = 0; ii < ArraySz; ii++ ){
      array[ii] = (float) ii;
   }
   memset( saveArray, 0, sizeof( float )*ArraySz );

   // get factory & create stream

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());
   StreamImpl* testStream                       = brookStreamFactory.createStreamImpl( OpenMM::BrookStreamFactory::BondedAtomIndicesStream, ArraySz, Stream::Float, platform );

   // load & retreive data

   testStream->loadFromArray( array );
   testStream->saveToArray( saveArray );

   // test for equality

   for( int ii = 0; ii < ArraySz; ii++ ){
      ASSERT_EQUAL( array[ii], saveArray[ii] );
   }

   delete saveArray;
   delete array;
 */  
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse parameter file lines

   Should be moved to Utilities file

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args; if return value equals maxTokens, then more tokens than allocated

   --------------------------------------------------------------------------------------- */

char* strsepLocal( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nTinkerParameterSet::strsep"

   char *s;
   const char *spanp;
   int c, sc;
   char *tok;

   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }

   for( tok = s;; ){
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
            return( tok );
         }
      } while( sc != 0 );
   }
}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   // ---------------------------------------------------------------------------------------

// (void) fprintf( stdout, "\nIn SimTKOpenMMUtilities::tokenizeString <%s>", lineBuffer );
// (void) fflush( stdout );

   char *ptr_c = NULL;

   for( ; (ptr_c = strsepLocal( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

   Local version of strncasecmp (missing in Windows) (static method) (Simbios)

   @param string1                 first string
   @param string2                 second string
   @param matchLength             match length

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

/*
int SimTKOpenMMUtilities::localStrncasecmp( const char *string1, const char *string2, int matchLength ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::localStrncasecmp"

   char ch1,ch2;

   // ---------------------------------------------------------------------------------------

   if( matchLength == 0 ){ 
      return SimTKOpenMMCommon::DefaultReturn;
   }

   do {   
      ch1 = toupper(*(string1++));
      ch2 = toupper(*(string2++));
      if( ch1 != ch2 )return (ch1-ch2);
      matchLength--;
   } while( ch1 && matchLength ); 

   return SimTKOpenMMCommon::DefaultReturn;  

}
*/

/**---------------------------------------------------------------------------------------

   Check that string is valid integer

   @param stringToCheck string to check

   @return true if string is a valid integer

   --------------------------------------------------------------------------------------- */

/*
bool SimTKOpenMMUtilities::isValidInteger( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::isValidInteger";

   // ---------------------------------------------------------------------------------------

   int ii;

   return checkString<int>(ii, stringToCheck, std::dec );
}
*/

/**---------------------------------------------------------------------------------------

   Check that string is valid RealOpenMM

   @param stringToCheck string to check

   @return true if string is a valid RealOpenMM

   --------------------------------------------------------------------------------------- */

/*
bool SimTKOpenMMUtilities::isValidRealOpenMM( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::isValidRealOpenMM";

   // ---------------------------------------------------------------------------------------

   RealOpenMM ii;

   return checkString<RealOpenMM>(ii, stringToCheck, std::dec );
}
*/

/**---------------------------------------------------------------------------------------

   Read file into string vector (Simbios) 

   @param fileName       file name
   @param fileContents   string vector containing file contents upon return
                         one string per line

   @return SimTKOpenMMCommon::DefaultReturn unless file could not be opened  

   --------------------------------------------------------------------------------------- */

/*
int SimTKOpenMMUtilities::readFileIntoStringVector( const std::string& fileName,
                                                    StringVector& fileContents ){

   // ---------------------------------------------------------------------------------------

   static int bufferSize = 2048;
   static char buffer[2048];
   static const std::string methodName = "\nSimTKOpenMMUtilities::readFileIntoStringVector";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* file = NULL;
#ifdef WIN32
   fopen_s( &file, fileName.c_str(), "r" );
#else
   file = fopen( fileName.c_str(), "r" );
#endif

   if( file != NULL ){
      std::stringstream message;
      message << methodName.c_str() << " opened file=<" <<  fileName.c_str() << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName.c_str() << " could not open file=<" <<  fileName.c_str() << ">.";
//(void) fprintf( stderr, "\n%s\n", message.str().c_str() ); 
//(void) fflush( stderr );
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;  
   }

   // loop over all lines in file, checking for end-of-file

   int lineNumber = 0;

   while( !feof( file ) ){ 

      // read next line

      int bufferLen;
      lineNumber++;
      if( fgets( buffer, bufferSize, file ) != NULL ){
         bufferLen = (int) strlen( buffer );
         if( bufferLen > 0 ){
            buffer[bufferLen-1] = '\0';
         }
         // (void) fprintf( log, "%s", buffer );
         // (void) fflush( log );
         fileContents.push_back( buffer );
      }
   }

   // done

   (void) fclose( file );

   if( file ){
      //std::stringstream message;
      //message << methodName.c_str() << " read " << lineNumber << " lines from file=<" << fileName->c_str() << ">.";
      // AmoebaLog::printMessage( message );
   }

   return SimTKOpenMMCommon::DefaultReturn;  
}
*/

/**---------------------------------------------------------------------------------------

   Read file into string vector (Simbios) 

   @param charArray      character array 
   @param arrayLength    array length
   @param arrayContents  string vector containing array contents upon return
                         one string per line

   @return SimTKOpenMMCommon::DefaultReturn unless file could not be opened  

   --------------------------------------------------------------------------------------- */

/*
int SimTKOpenMMUtilities::readCharacterArrayIntoStringVector( const char* charArray, int arrayLength,
                                                              StringVector& fileContents ){

   // ---------------------------------------------------------------------------------------

   static int bufferSize = 2048;
   static char buffer[2048];
   static const std::string methodName = "\nSimTKOpenMMUtilities::readCharacterArrayIntoStringVector";

   // ---------------------------------------------------------------------------------------

   // loop over all lines in file, checking for end-of-file

   int byteIndex = 0;
   //int lineIndex = 0;
   while( byteIndex < arrayLength ){ 

      // get next line

      int lineLength = 0;
      int done       = 0;
      while( byteIndex < arrayLength && lineLength < bufferSize && !done ){
         buffer[lineLength++] = charArray[byteIndex++];
         if( (int) charArray[byteIndex] == 10 || (int) charArray[byteIndex] == 13 ){
            while( (int) charArray[byteIndex] < 32 && byteIndex < arrayLength )byteIndex++;
            done = 1;
         }
      }
      buffer[lineLength] = '\0';
      //lineIndex++;

      // (void) fprintf( stdout, "%s readCharacterArrayIntoStringVector: %d %d <%s>", methodName.c_str(), lineIndex, lineLength, buffer );
      // (void) fflush( stdout );
      fileContents.push_back( buffer );
   }

   return SimTKOpenMMCommon::DefaultReturn;  
}
*/

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param line                 string to tokenize
   @param tokenVector          upon return vector of tokens
   @param delimiter            token delimter
   @param clearTokenVector     if true, clear tokenVector

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

/*
int SimTKOpenMMUtilities::tokenizeString( const std::string& line, StringVector& tokenVector,
                                          const std::string& delimiter, int clearTokenVector ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   static int bufferSz       = 8192;
   static char* lineBuffer   = NULL;

   // ---------------------------------------------------------------------------------------

   char *ptr_c;

   // clear token vector

   if( clearTokenVector ){
      tokenVector.clear();
   }   

   // allocate space for line buffer and copy via sprintf()

   if( lineBuffer == NULL || bufferSz < (int) line.size() ){
      if( lineBuffer != NULL ){
         free( lineBuffer );
      }
      if( bufferSz < (int) line.size() ){
         bufferSz = (int) (2*line.size());
      } 
      lineBuffer = (char*) malloc( bufferSz*sizeof( char ) );
   }

#ifdef WIN32
   (void) sprintf_s( lineBuffer, bufferSz, "%s", line.c_str() );
#else
   (void) sprintf( lineBuffer, "%s", line.c_str() );
#endif

   // parse

   while( (ptr_c = SimTKOpenMMUtilities::strsep( &lineBuffer, delimiter.c_str() )) != NULL ){
      if( *ptr_c ){
         tokenVector.push_back( std::string( ptr_c ) );
      }   
   }

   return SimTKOpenMMCommon::DefaultReturn;  
}
*/

/**---------------------------------------------------------------------------------------

   Return lower case copy of string

   @param string                  string

   @return lower cased string

   --------------------------------------------------------------------------------------- */

//int SimTKOpenMMUtilities::toLower( std::string& string ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::toLower"

   // ---------------------------------------------------------------------------------------

	// transform string to lower case
	    
	// std::transform( string.begin(), string.end(), string.begin(), (int(*)(int)) std::tolower);
	//return SimTKOpenMMCommon::DefaultReturn;
		
//}

/**---------------------------------------------------------------------------------------

   Read parameter file

   @param	parameterFileName	parameter file name

   @return AmoebaCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int readParameterFile( int numberOfAtomIndices, int numberOfParameters,
                       vector<vector<int> >& atomIndices,
                       vector<vector<double> >& parameters,
                       const std::string& parameterFileName ){

   // ---------------------------------------------------------------------------------------

/*
   const int bufferSize = 1024;
   char buffer[bufferSize];
   static const std::string methodName              = "\nreadParameterFile";

   // ---------------------------------------------------------------------------------------

   std::ifstream parameterFileStream( parameterFileName.c_str() );

   if( parameterFileStream == NULL ){
      std::stringstream message;
      message << methodName.c_str() << " Could not open file=<" << parameterFileName.c_str() << ">.";
      cout << message.str().c_str();
      return -1;
   } else {
      std::stringstream message;
      message << methodName.c_str() << " Opened file=<" << parameterFileName.c_str() << ">.";
      //AmoebaLog::printMessage( message );
   }

   int lineCount     = 1;

   // skip first line

   parameterFileStream.getline( buffer, bufferSize );
   lineCount++;

   int startColumnIndex = 1;
   const std::string delimters = " ";
   while( parameterFileStream.getline( buffer, bufferSize ) ){

      // std::stringstream message;
      // message << "\n<" << buffer << ">";

      StringVector tokens;
      tokenizeString( buffer, tokens, delimters );
      if( tokens.size() < 4 ){
         std::stringstream message;
         message << "\nProblem w/ line=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      } else {
         std::stringstream message;
         message << "\nline=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      }

      vector<int> entry;
      atomIndices.push_back( entry );

      for( int ii = 0; ii < numberOfAtomIndices; ii++ ){ 
         int atomIndex = (int) atoi( tokens[startColumnIndex+ii].c_str() );
         entry.push_back( atomIndex );
      }

      vector<double> entryP;
      parameters.push_back( entryP );
      int startParameterIndex = startColumnIndex + numberOfAtomIndices;
      for( int ii = 0; ii < numberOfParameters; ii++ ){ 
         double parameter = (double) atof( tokens[startParameterIndex+ii].c_str() );
         entryP.push_back( parameter );
      }
      lineCount++;
   }

*/
   return 0;
}

/**---------------------------------------------------------------------------------------

   Read parameter file

   @param	parameterFileName	parameter file name

   @return AmoebaCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int readLJCoulombParameterFile( std::vector<std::vector<double> >& parameters,
                                std::vector<std::string>& atomType,
                                std::vector<std::vector<int> >& exclusions,
                                const std::string& parameterFileName ){

   // ---------------------------------------------------------------------------------------

/*
   const int bufferSize = 1024;
   char buffer[bufferSize];
   static const std::string methodName   = "\nreadLJCoulombParameterFile";

   // ---------------------------------------------------------------------------------------

   std::ifstream parameterFileStream( parameterFileName.c_str() );

   if( parameterFileStream == NULL ){
      std::stringstream message;
      message << methodName.c_str() << " Could not open file=<" << parameterFileName.c_str() << ">.";
      cout << message.str().c_str();
      return -1;
   } else {
      std::stringstream message;
      message << methodName.c_str() << " Opened file=<" << parameterFileName.c_str() << ">.";
      //AmoebaLog::printMessage( message );
   }

   int lineCount     = 1;

   // skip first line

   parameterFileStream.getline( buffer, bufferSize );
   lineCount++;

   int startColumnIndex           = 1;
   const std::string delimters    = " ";
   int startParameterIndex        = 1;
   int numberOfParameters         = 5;
   while( parameterFileStream.getline( buffer, bufferSize ) ){

      // std::stringstream message;
      // message << "\n<" << buffer << ">";

      StringVector tokens;
      tokenizeString( buffer, tokens, delimters );
      if( tokens.size() < 7 ){
         std::stringstream message;
         message << "\nProblem w/ line=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      } else {
         std::stringstream message;
         message << "\nline=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      }

      // float & atomtype parameters

      vector<double> entryP;
      parameters.push_back( entryP );
      for( int ii = 0; ii < numberOfParameters; ii++ ){ 
         if( ii == 3 ){
            atomType.push_back( tokens[startParameterIndex+ii] );
         } else {
            double parameter = (double) atof( tokens[startParameterIndex+ii].c_str() );
            entryP.push_back( parameter );
         }
      }

      // exclusions

      vector<int> entry;
      exclusions.push_back( entry );
      int numberOfExclusions = atoi( tokens[6].c_str() );
      for( int ii = 0; ii < numberOfExclusions; ii++ ){ 
         int atomIndex = (int) atoi( tokens[startColumnIndex+ii].c_str() );
         entry.push_back( atomIndex );
      }

      lineCount++;
   }
*/

   return 0;
}

/**---------------------------------------------------------------------------------------

   Read output file

   @param	coordinates output vector of coordinates
   @param	forces      output vector of forces
   @param	fileName    input file name

   @return 0

   --------------------------------------------------------------------------------------- */

int readCoordinateForceOutputFile( std::vector<std::vector<double> >& coordinates,
                                   std::vector<std::vector<double> >& forces,
                                   const std::string& fileName ){

   // ---------------------------------------------------------------------------------------

/*
   const int bufferSize = 1024;
   char buffer[bufferSize];
   static const std::string methodName   = "\nreadCoordinateForceOutputFile";

   // ---------------------------------------------------------------------------------------

   // open file

   std::ifstream parameterFileStream( fileName.c_str() );

   if( outputFileStream == NULL ){
      std::stringstream message;
      message << methodName.c_str() << " Could not open file=<" << fileName.c_str() << ">.";
      cout << message.str().c_str();
      return -1;
   } else {
      std::stringstream message;
      message << methodName.c_str() << " Opened file=<" << fileName.c_str() << ">.";
      //AmoebaLog::printMessage( message );
   }

   int lineCount     = 1;

   // skip first line

   outputFileStream.getline( buffer, bufferSize );
   lineCount++;

   int startCoordinateIndex       = 1;
   int startForceIndex            = 4;
   const std::string delimters    = " ";
   while( outputFileStream.getline( buffer, bufferSize ) ){

      // std::stringstream message;
      // message << "\n<" << buffer << ">";

      StringVector tokens;
      tokenizeString( buffer, tokens, delimters );
      if( tokens.size() < 7 ){
         std::stringstream message;
         message << "\nProblem w/ line=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      } else {
         std::stringstream message;
         message << "\nline=" << lineCount << " <" << buffer << "> size=" << tokens.size();
         cout << message.str().c_str();
*/
/*
      }

      // coordinates

      vector<double> entryC;
      coordinates.push_back( entryC );
      for( int ii = 0; ii < 3; ii++ ){ 
         double coordinate = (double) atof( tokens[startCoordinateIndex+ii].c_str() );
         entryC.push_back( parameter );
      }

      // forces

      vector<double> entryF;
      coordinates.push_back( entryF );
      for( int ii = 0; ii < 3; ii++ ){ 
         double parameter = (double) atof( tokens[startForceIndex+ii].c_str() );
         entryF.push_back( parameter );
      }

      lineCount++;
   }
*/
   return 0;
}
/*
class ParameterInfo  {

   private:

      int _numberOfAtomIndices;
      int _numberOfParameterIndices;
      std::string _fileName;
      std::string _directoryName;
      vector<vector<int> > _atomIndices;
      vector<vector<double> > _parameters;

   public:

      ParameterInfo( int numberOfAtomIndices, int numberOfParameterIndices, std::string fileName, std::string directoryName ){
         _numberOfAtomIndices        = numberOfAtomIndices; 
         _numberOfParameterIndices   = numberOfParameterIndices;
         _fileName                   = fileName;        
         _directoryName              = directoryName;
      }
      ~ParameterInfo( void ){};
      std::string getFullFileName( void ){ std::string name; name = _directoryName; name.append( _fileName ); return name; }
      int getNumberOfAtomIndices(  void ){ return _numberOfAtomIndices; };
      int getNumberOfParameters(   void ){ return _numberOfParameterIndices; };

      vector<vector<int> >& getAtomIndices( void ){ return _atomIndices; };
      vector<vector<double> >& getParameters( void ){ return _parameters; };
};
*/

void testBrookBonded( void ){

/*
   static const int debug = 1;
   FILE* log              = stdout;
   
   // directory & file names

   std::string parameterDirectoryName    = "C:\\cygwin\\home\\friedrim\\nvidia\\ParameterFiles\\villin\\"; 
   std::string resultsDirectoryName      = "C:\\cygwin\\home\\friedrim\\nvidia\\OutputFiles\\villin\\"; 

   std::string harmonicFileName          = "GromacsHarmonicBondParameter.txt"; 
   std::string angleFileName             = "GromacsAngleBondParameter.txt"; 
   std::string properDihedralFileName    = "GromacsProperDihedralParameter.txt"; 
   std::string rbDihedralFileName        = "GromacsRbDihedralParameter.txt"; 
   std::string lj14FileName              = "GromacsLJ14Parameter.txt"; 
   std::string ljCoulombFileName         = "GromacsLJCoulombParameter.txt"; 

   std::string bondedReferenceFileName   = "ReferenceBonded.txt"; 

   vector<ParameterInfo> parameterInfoVector;
   parameterInfoVector.push_back( ParameterInfo( 2, 2, harmonicFileName,       directoryName ) );
   parameterInfoVector.push_back( ParameterInfo( 3, 2, angleFileName,          directoryName ) );
   parameterInfoVector.push_back( ParameterInfo( 4, 3, properDihedralFileName, directoryName ) );
   parameterInfoVector.push_back( ParameterInfo( 4, 6, rbDihedralFileName,     directoryName ) );
   parameterInfoVector.push_back( ParameterInfo( 2, 4, lj14FileName,           directoryName ) );

   // read in parameter files

   for( unsigned int ii = 0; ii < parameterInfoVector.size(); ii++ ){
      ParameterInfo parameterInfo = parameterInfoVector[ii];
      readParameterFile( parameterInfo.getNumberOfAtomIndices(), parameterInfo.getNumberOfParameters(),
                         parameterInfo.getAtomIndices(), parameterInfo.getParameters(), parameterInfo.getFullFileName() );
      if( debug ){
         (void) fprintf( log, "%s %d\n", parameterInfo.getFullFileName().c_str(), parameterInfo.getParameters().size() ); 
      }
   }

   // LJ/Coulomb parameter file handled separately

   ParameterInfo ljCoulombParameterInfo( 0, 4, ljCoulombFileName, directoryName );
   vector<std::string> atomTypes;
   vector<vector<int>> exclusions;
   readLJCoulombParameterFile( ljCoulombParameterInfo.getParameters(), atomTypes, exclusions, ljCoulombParameterInfo.getFullFileName() );
   if( debug ){
      (void) fprintf( log, "%s param=%d atomTypes=%d exclusions=%d\n", ljCoulombParameterInfo.getFullFileName().c_str(), ljCoulombParameterInfo.getParameters().size(), atomTypes.size(), exclusions.size() ); 
   }

   BrookBonded brookBonded;
   BrookPlatform brookPlatform( 32 );
   double lj14Scale          = 8.33300e-001;
   double coulombScale       = 1.0;
   brookBonded.setup( atomTypes.size(),
                      parameterInfoVector[0].getAtomIndices(), parameterInfoVector[0].getParameters(),
                      parameterInfoVector[1].getAtomIndices(), parameterInfoVector[1].getParameters(),
                      parameterInfoVector[2].getAtomIndices(), parameterInfoVector[2].getParameters(),
                      parameterInfoVector[3].getAtomIndices(), parameterInfoVector[3].getParameters(),
                      parameterInfoVector[4].getAtomIndices(), ljCoulombParameterInfo.getParameters(),
                      lj14Scale, coulombScale,  brookPlatform, log );

   // read in coordinates and forces

   vector<std::double> coordinates;
   vector<std::double> forces;
   std::string referenceFileName      = resultsDirectoryName + bondedReferenceFileName; 
   readCoordinateForceOutputFile( coordinates, forces, referenceFileName );
exit(0);

   //brookBonded.calculateBondedForce( StreamImpl& positionsStreamImpl, StreamImpl& forcesImpl )
*/
/*
   // read bonded file
   BrookPlatform platform;

   // create and initialize arrays

   float* array      = new float[ArraySz];
   float* saveArray  = new float[ArraySz];
   for( int ii = 0; ii < ArraySz; ii++ ){
      array[ii] = (float) ii;
   }
   memset( saveArray, 0, sizeof( float )*ArraySz );

   // get factory & create stream

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());
   StreamImpl* testStream                       = brookStreamFactory.createStreamImpl( OpenMM::BrookStreamFactory::BondedAtomIndicesStream, ArraySz, Stream::Float, platform );

   // load & retreive data

   testStream->loadFromArray( array );
   testStream->saveToArray( saveArray );

   // test for equality

   for( int ii = 0; ii < ArraySz; ii++ ){
      ASSERT_EQUAL( array[ii], saveArray[ii] );
   }

   delete saveArray;
   delete array;
*/
   
}

void testBrookBondedHarmonicBond( void ){

   static const int debug = 1;
   FILE* log              = stdout;
   
   int numberOfAtoms      = 2;
   RealOpenMM mass        = 2.0;

   System system( numberOfAtoms, 0);
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      system.setAtomMass( ii, mass );
   }

   LangevinIntegrator integrator(0, 0.1, 0.01);

   BrookBonded brookBonded;
   brookBonded.setLog( log );
   BrookPlatform brookPlatform( 32, "cal", log );

   double lj14Scale          = 8.33300e-001;
   double coulombScale       = 1.0;

   std::vector<std::vector<int> >    harmonicBondsIndices;
   std::vector<int>  harmonicBonds;
   harmonicBonds.push_back( 0 );
   harmonicBonds.push_back( 1 );
   harmonicBondsIndices.push_back( harmonicBonds ); 

   std::vector<double>  parameters;
   parameters.push_back( 1.0 );
   parameters.push_back( 1.0 );
   std::vector<std::vector<double> > harmonicBondsParameters;
   harmonicBondsParameters.push_back( parameters );

   std::vector<std::vector<int> >    angleBondsIndices;
   std::vector<std::vector<double> > angleBondsParameters;

   std::vector<std::vector<int> >    rbTorsionBondsIndices;
   std::vector<std::vector<double> > rbTorsionBondsParameters;

   std::vector<std::vector<int> >    properTorsionBondsIndices;
   std::vector<std::vector<double> > properTorsionBondsParameters;

   std::vector<std::vector<int> >    lj14Indices;
   std::vector<std::vector<double> > lj14Parameters;

   (void) fprintf( log, "testBrookBondedHarmonicBond: Calling brookBonded.setup\n" );
   brookBonded.setup( numberOfAtoms,
                      harmonicBondsIndices,         harmonicBondsParameters,
                      angleBondsIndices,            angleBondsParameters,
                      rbTorsionBondsIndices,        rbTorsionBondsParameters,
                      properTorsionBondsIndices,    properTorsionBondsParameters,
                      lj14Indices,                  lj14Parameters,
                      lj14Scale, coulombScale,  brookPlatform );
  
   std::string contents = brookBonded.getContentsString( 0 );
   (void) fprintf( log, "testBrookBondedHarmonicBond: brookBonded::contents\n%s", contents.c_str() );
   (void) fflush( log );

   BrookStochasticDynamics    brookStochasticDynamics;
   BrookShakeAlgorithm        brookShakeAlgorithm;
   BrookRandomNumberGenerator brookRandomNumberGenerator;

   contents = brookStochasticDynamics.getContentsString( );
   (void) fprintf( log, "testBrookBondedHarmonicBond: brookStochasticDynamics::contents\n%s", contents.c_str() );

   contents = brookShakeAlgorithm.getContentsString( );
   (void) fprintf( log, "testBrookBondedHarmonicBond: brookShakeAlgorithm::contents\n%s", contents.c_str() );

   contents = brookRandomNumberGenerator.getContentsString( );
   (void) fprintf( log, "testBrookBondedHarmonicBond: brookRandomNumberGenerator::contents\n%s", contents.c_str() );

}

void testBrookNonBonded( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookNonBonded";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   System system( numberOfAtoms, 0);
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      system.setAtomMass( ii, mass );
   }

   LangevinIntegrator integrator(0, 0.1, 0.01);

   BrookNonBonded brookNonBonded;
   brookNonBonded.setLog( log );
   BrookPlatform brookPlatform( 32, "cal", log );

   // load nonbonded parameters

   std::vector<std::vector<double> > nonbondedParameters;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      std::vector<double>  parameters;
      parameters.push_back( 1.0 );
      parameters.push_back( 1.0 );
      parameters.push_back( 1.0 );
      nonbondedParameters.push_back( parameters );
   }

   // exclusions

   std::vector<std::set<int> > exclusions;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      std::set<int>  exclusion;
      exclusions.push_back( exclusion );
   }

   (void) fprintf( log, "testBrookBondedHarmonicBond: Calling brookNonBonded::setup\n" );
   (void) fflush( log );
   brookNonBonded.setup( numberOfAtoms, nonbondedParameters, exclusions, brookPlatform );
   std::string contents = brookNonBonded.getContentsString( );
   (void) fprintf( log, "testBrookBondedHarmonicBond: Called brookNonBonded.getContentsString \n" );
   (void) fflush( log );
   (void) fprintf( log, "testBrookBondedHarmonicBond: brookNonBonded::contents\n%s", contents.c_str() );
   (void) fprintf( log, "testBrookBondedHarmonicBond: exiting\n" );
   (void) fflush( log );
}

void testBrookBonds( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookBonds";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 3;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator(0, 0.1, 0.01);

   // int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   StandardMMForceField* forceField = new StandardMMForceField( 3, 2, 0, 0, 0 ); 

   // ( index, atom1, atom2, length, k )
   forceField->setBondParameters(0, 0, 1, 1.5, 0.8);
   forceField->setBondParameters(1, 1, 2, 1.2, 0.7);
   system.addForce(forceField);

   //(void) fprintf( log, "testBrookBonds:Calling context\n");
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(3);

   positions[0] = Vec3(0, 2, 0); 
   positions[1] = Vec3(0, 0, 0); 
   positions[2] = Vec3(1, 0, 0); 

   context.setPositions(positions);

   //(void) fprintf( log, "testBrookBonds:Calling getState\n");
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "Harmonic bond forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );

   ASSERT_EQUAL_VEC(Vec3(0, -0.8*0.5, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(0.7*0.2, 0, 0), forces[2], TOL);
   ASSERT_EQUAL_VEC(Vec3(-forces[0][0]-forces[2][0], -forces[0][1]-forces[2][1], -forces[0][2]-forces[2][2]), forces[1], TOL);

   (void) fprintf( log, "Harmonic bond forces ok -- skipping energy test\n");

   // ASSERT_EQUAL_TOL(0.5*0.8*0.5*0.5 + 0.5*0.7*0.2*0.2, state.getPotentialEnergy(), TOL);
}

void testBrookAngles( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookAngles";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 4;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, 0, 2, 0, 0 ); 

   // int index, int atom1, int atom2, int atom3, double angle, double k
   forceField->setAngleParameters(0, 0, 1, 2, PI_M/3, 1.1);
   forceField->setAngleParameters(1, 1, 2, 3, PI_M/2, 1.2);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(2, 1, 0);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "Angle bond forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );

   double torque1 = 1.1*PI_M/6;
   double torque2 = 1.2*PI_M/4;
   ASSERT_EQUAL_VEC(Vec3(torque1, 0, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(-0.5*torque2, 0.5*torque2, 0), forces[3], TOL); // reduced by sqrt(2) due to the bond length, another sqrt(2) due to the angle
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), TOL);

   (void) fprintf( log, "Angle bond forces ok -- skipping energy test\n");

   // ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), TOL);
}

void testBrookPeriodicTorsions( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "PeriodicTorsions";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 4;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, 0, 0, 1, 0 ); 

   // int index, int atom1, int atom2, int atom3, double angle, double k
   forceField->setPeriodicTorsionParameters(0, 0, 1, 2, 3, 2, PI_M/3, 1.1);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(1, 0, 2);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "Periodic torsion bond forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );


   double torque = -2*1.1*std::sin(2*PI_M/3);
   ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, 0), forces[3], TOL);
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), TOL);
   // ASSERT_EQUAL_TOL(1.1*(1+std::cos(2*PI_M/3)), state.getPotentialEnergy(), TOL);

   (void) fprintf( log, "Periodic torsion bond forces ok -- skipping energy test\n");

   // ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), TOL);
}

void testBrookRBTorsions( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "RBTorsions";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 4;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int numAtoms, int numBonds, int numAngles, int numPeriodicTorsions, int numRBTorsions

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, 0, 0, 0, 1 ); 

   forceField->setRBTorsionParameters(0, 0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   positions[0] = Vec3(0, 1, 0);
   positions[1] = Vec3(0, 0, 0);
   positions[2] = Vec3(1, 0, 0);
   positions[3] = Vec3(1, 1, 1);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "RB torsion bond forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );

   double psi    = 0.25*PI_M - PI_M;
   double torque = 0.0;
   for (int i = 1; i < 6; ++i) {
      double c  = 0.1*(i+1);
      torque   += -c*i*std::pow(std::cos(psi), i-1)*std::sin(psi);
   }
   ASSERT_EQUAL_VEC(Vec3(0, 0, torque), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(0, 0.5*torque, -0.5*torque), forces[3], TOL);
   ASSERT_EQUAL_VEC(Vec3(forces[0][0]+forces[1][0]+forces[2][0]+forces[3][0],
                         forces[0][1]+forces[1][1]+forces[2][1]+forces[3][1],
                         forces[0][2]+forces[1][2]+forces[2][2]+forces[3][2]),
                         Vec3(0, 0, 0), TOL);
/*
    double energy = 0.0;
    for (int i = 0; i < 6; ++i) {
        double c = 0.1*(i+1);
        energy += c*std::pow(std::cos(psi), i);
    }
    ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
*/


   (void) fprintf( log, "RB torsion bond forces ok -- skipping energy test\n");

   // ASSERT_EQUAL_TOL(0.5*1.1*(PI_M/6)*(PI_M/6) + 0.5*1.2*(PI_M/4)*(PI_M/4), state.getPotentialEnergy(), TOL);
}

void testBrookCoulomb( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "Coulomb";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, 0, 0, 0, 0 ); 
   forceField->setAtomParameters(0, 0.5, 1, 0);
   forceField->setAtomParameters(1, -1.5, 1, 0);
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   positions[0] = Vec3(0, 0, 0);
   positions[1] = Vec3(2, 0, 0);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "Coulomb forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );

   double force = 138.935485*(-0.75)/4.0;
   ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
   //ASSERT_EQUAL_TOL(138.935485*(-0.75)/2.0, state.getPotentialEnergy(), TOL);

   (void) fprintf( log, "Coulomb forces ok -- skipping energy test\n");

}

void testBrookLJ( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LJ";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   // ReferencePlatform platform;
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, 0, 0, 0, 0 ); 
   forceField->setAtomParameters(0, 0, 1.2, 1); 
   forceField->setAtomParameters(1, 0, 1.4, 2); 
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   positions[0] = Vec3(0, 0, 0);
   positions[1] = Vec3(2, 0, 0);

   context.setPositions(positions);

   //(void) fprintf( log, "%s :Calling getState\n", methodName.c_str() );
   //(void) fflush( log );

   State state = context.getState( State::Forces );

   const vector<Vec3>& forces = state.getForces();
   (void) fprintf( log, "LJ forces\n");
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
   }
   (void) fflush( log );

   double x = 1.3/2.0;
   double eps = SQRT_TWO;
   double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/2.0;
   ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
   ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[1], TOL);
   //ASSERT_EQUAL_TOL(4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0)), state.getPotentialEnergy(), TOL);

   (void) fprintf( log, "LJ forces ok -- skipping energy test\n");

}

void testBrookExclusionsAnd14( void ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "ExclusionsAnd14";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 5;
   RealOpenMM mass                          = 2.0;

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cpu", log );
   //BrookPlatform platform( 32, "cal", log );
   //ReferencePlatform platform;
   System system( numberOfAtoms, 0 ); 
   LangevinIntegrator integrator( 0, 0.1, 0.01 );

   // int index, double charge, double radius, double depth

   StandardMMForceField* forceField = new StandardMMForceField( numberOfAtoms, numberOfAtoms-1, 0, 0, 0 ); 
   for( int ii = 1; ii < numberOfAtoms; ii++ ){
      forceField->setBondParameters(ii-1, ii-1, ii, 1, 0);
   }
   system.addForce(forceField);

   //(void) fprintf( log, "%s: Calling context\n",  methodName.c_str() );
   //(void) fflush( log );

   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numberOfAtoms);

   const double r = 1.0;
   positions[0] = Vec3(0, 0, 0);
   for( int ii = 1; ii < numberOfAtoms; ii++ ){
      positions[ii] = Vec3(r, 0, 0);
   }
   for( int ii = 1; ii < numberOfAtoms; ii++ ){

      // Test LJ forces

     forceField->setAtomParameters(0, 0, 1.5, 1);
      for (int jj = 1; jj < numberOfAtoms; ++jj){
         forceField->setAtomParameters(jj, 0, 1.5, 0);
      }
      forceField->setAtomParameters(ii, 0, 1.5, 1);
      context.reinitialize();
      context.setPositions(positions);
      State state = context.getState( State::Forces );
      const vector<Vec3>& forces = state.getForces();
      double x = 1.5/r;
      double eps = 1.0;
      double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/r;
      double energy = 4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0));
      if( ii == 3 ){
         force *= 0.5;
         energy *= 0.5;
      }
      if( ii < 3 ){
         force = 0;
         energy = 0;
      }
      (void) fprintf( log, "14 LJ forces ii=%d F=%.6e\n", ii, force );
      for( int jj = 0; jj < numberOfAtoms; jj++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", jj, forces[jj][0], forces[jj][1], forces[jj][2] );
      }
      (void) fflush( log );
      ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
      ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[ii], TOL);
      (void) fprintf( log, "14 LJ forces ok for index=%d\n\n", ii );
      (void) fflush( log );
      //ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);

      // Test Coulomb forces

      forceField->setAtomParameters(0, 2, 1.5, 0);
      forceField->setAtomParameters(ii, 2, 1.5, 0);
      context.reinitialize();
      context.setPositions(positions);
      state = context.getState( State::Forces );
      const vector<Vec3>& forces2 = state.getForces();
      force = 138.935485*4/(r*r);
      energy = 138.935485*4/r;
      if( ii == 3 ){
         force /= 1.2;
         energy /= 1.2;
      }
      if( ii < 3 ){
          force = 0;
          energy = 0;
      }
      (void) fprintf( log, "14 Coulomb forces ii=%d F=%.6e\n", ii, force );
      for( int jj = 0; jj < numberOfAtoms; jj++ ){
         (void) fprintf( log, "%d [%.5e %.5e %.5e]\n", jj, forces[jj][0], forces[jj][1], forces[jj][2] );
      }
      (void) fflush( log );
      ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces2[0], TOL);
      ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces2[ii], TOL);
      (void) fprintf( log, "14 Coulomb forces ok for index=%d\n\n", ii );
      (void) fflush( log );
       // ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
   }
   (void) fprintf( log, "ExclusionsAnd14 ok -- skipping energy test\n");

}

int testBrookStreams( int streamSize ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "testBrookStreams";
   static const int debug                   = 1;
   static const int maxPrint                = 10;
   float epsilon                            = 1.0e-04f;
   FILE* log                                = stdout;
   int totalErrors                          = 0;
   
// ---------------------------------------------------------------------------------------

   float* testArray                          = new float[streamSize];
   float* saveArray                          = new float[streamSize];
   memset( testArray, 0 , streamSize*sizeof( float ) );
   memset( saveArray, 0 , streamSize*sizeof( float ) );

   float value                               = 11.0f;
   BrookFloatStreamInternal* brookStream     = new BrookFloatStreamInternal( "testFloatStream", streamSize, 32, BrookStreamInternal::Float, -1.0 );
   brookStream->fillWithValue( &value );
   brookStream->saveToArray( saveArray );

   (void) fprintf( log, "\n%s\n", methodName.c_str() );
   int errors = 0;
   (void) fprintf( log, "Fill test for float stream of size=%d\n", streamSize );
   for( int ii = 0; ii < streamSize; ii++ ){
      if( fabsf( value - saveArray[ii] ) > epsilon ){
         if( errors < maxPrint ){
            (void) fprintf( log, "Error entry %d [%.2f %.2f]\n", ii, testArray[ii], saveArray[ii] );
         }
         errors++;
      }
   }

   if( errors ){
      (void) fprintf( log, "%d errors detected for fill-value test.\n", errors );
      totalErrors += errors;
   } else {
      (void) fprintf( log, "No errors detected for fill-value test.\n" );
   }
   (void) fflush( log );

   delete brookStream;
   delete[] testArray;
   delete[] saveArray;

// ---------------------------------------------------------------------------------------

   // test 4 float types

   for( int jj = 1; jj <= 4 ; jj++ ){

      int totalSize                             = streamSize*jj;
      BrookStreamInternal::DataType dataType;
      switch( jj ){
         case 1:

            dataType    = BrookStreamInternal::Float;
            break;

         case 2:

            dataType    = BrookStreamInternal::Float2;
            break;

         case 3:

            dataType    = BrookStreamInternal::Float3;
            break;

         case 4:
            dataType    = BrookStreamInternal::Float4;
            break;

         default:
            dataType    = BrookStreamInternal::Float4;
            break;
      }


      char name[128];
      (void) sprintf( name, "testFloatStream%d", jj );
      BrookFloatStreamInternal* brookStream     = new BrookFloatStreamInternal( name, streamSize, 32, dataType, -1.0 );
      float* testArray                          = new float[totalSize];
      float* saveArray                          = new float[totalSize];

      for( int ii = 0; ii < totalSize; ii++ ){
         testArray[ii] = (float) ii;
      }
      brookStream->loadFromArray( testArray );
      brookStream->saveToArray( saveArray );
   
      int errors = 0;
      (void) fprintf( log, "Comparison test for float%d stream of size=%d\n", jj, totalSize );
      for( int ii = 0; ii < totalSize; ii++ ){
         if( fabsf( testArray[ii] - saveArray[ii] ) > epsilon ){
            if( errors < maxPrint ){
               (void) fprintf( log, "Error entry %d [%.2f %.2f]\n", ii, testArray[ii], saveArray[ii] );
            }
            errors++;
         }
      }
   
      if( errors ){
         (void) fprintf( log, "%d errors detected for comparison float%d test.\n", jj, errors );
         totalErrors += errors;
      } else {
         (void) fprintf( log, "No errors detected for comparison float%d test.\n", jj );
   
      }

//std::string message = brookStream->getContentsString();
//(void) fprintf( log, "Info %s.\n", message.c_str() );

      (void) fflush( log );
   
      delete[] testArray;
      delete[] saveArray;
      delete brookStream;
   }
      
   // ---------------------------------------------------------------------------------------
 
   return totalErrors;
}

void testLangevinSingleBond() {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinSingleBond";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0; 

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );

   System system(numberOfAtoms, 0);
   system.setAtomMass(0, 2.0);
   system.setAtomMass(1, 2.0);
   LangevinIntegrator integrator(0, 0.1, 0.01);
   StandardMMForceField* forceField = new StandardMMForceField(2, 1, 0, 0, 0);
   forceField->setBondParameters(0, 0, 1, 1.5, 1);
   system.addForce(forceField);
   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(2);
   positions[0] = Vec3(-1, 0, 0);
   positions[1] = Vec3(1, 0, 0);
   context.setPositions(positions);
   
   // This is simply a damped harmonic oscillator, so compare it to the analytical solution.
   
   double freq = std::sqrt(1-0.05*0.05);
   for (int i = 0; i < 1000; ++i) {
       State state = context.getState(State::Positions | State::Velocities);
       double time = state.getTime();
       double expectedDist = 1.5+0.5*std::exp(-0.05*time)*std::cos(freq*time);
       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedDist, 0, 0), state.getPositions()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedDist, 0, 0), state.getPositions()[1], 0.02);
       double expectedSpeed = -0.5*std::exp(-0.05*time)*(0.05*std::cos(freq*time)+freq*std::sin(freq*time));
       ASSERT_EQUAL_VEC(Vec3(-0.5*expectedSpeed, 0, 0), state.getVelocities()[0], 0.02);
       ASSERT_EQUAL_VEC(Vec3(0.5*expectedSpeed, 0, 0), state.getVelocities()[1], 0.02);
       integrator.step(1);
   }
   
   // Not set the friction to a tiny value and see if it conserves energy.
   
   integrator.setFriction(5e-5);
   context.setPositions(positions);
   State state = context.getState(State::Energy);
   double initialEnergy = state.getKineticEnergy()+state.getPotentialEnergy();
   for (int i = 0; i < 1000; ++i) {
       state = context.getState(State::Energy);
       double energy = state.getKineticEnergy()+state.getPotentialEnergy();
       ASSERT_EQUAL_TOL(initialEnergy, energy, 0.01);
       integrator.step(1);
   }
}

void testLangevinTemperature() {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinTemperature";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0; 

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );
   // ReferencePlatform platform;

    const int numAtoms = 8;
    const double temp = 100.0;
    //ReferencePlatform platform;
    System system(numAtoms, 0);
    LangevinIntegrator integrator(temp, 2.0, 0.01);
    StandardMMForceField* forceField = new StandardMMForceField(numAtoms, 0, 0, 0, 0);
    for (int i = 0; i < numAtoms; ++i) {
        system.setAtomMass(i, 2.0);
        forceField->setAtomParameters(i, (i%2 == 0 ? 1.0 : -1.0), 1.0, 5.0);
    }
    system.addForce(forceField);
    OpenMMContext context(system, integrator, platform);
    vector<Vec3> positions(numAtoms);
    for (int i = 0; i < numAtoms; ++i)
        positions[i] = Vec3((i%2 == 0 ? 2 : -2), (i%4 < 2 ? 2 : -2), (i < 4 ? 2 : -2));
    context.setPositions(positions);
    
    // Let it equilibrate.
    
    integrator.step(10000);
    
    // Now run it for a while and see if the temperature is correct.
    
    double ke = 0.0;
    for (int i = 0; i < 1000; ++i) {
        State state = context.getState(State::Energy);
        ke += state.getKineticEnergy();
        integrator.step(1);
    }
    ke /= 1000;
    double expected = 0.5*numAtoms*3*BOLTZ*temp;
    ASSERT_EQUAL_TOL(expected, ke, 3*expected/std::sqrt(1000.0));
}

void testLangevinConstraints() {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "LangevinConstraints";
   static const int debug                   = 1;
   FILE* log                                = stdout;
   
   int numberOfAtoms                        = 2;
   RealOpenMM mass                          = 2.0; 

// ---------------------------------------------------------------------------------------

   (void) fprintf( log, "%s\n", methodName.c_str() );
   (void) fflush( log );

   BrookPlatform platform( 32, "cal", log );

   const int numAtoms = 8;
   const double temp = 100.0;
   //ReferencePlatform platform;
   System system(numAtoms, numAtoms-1);
   LangevinIntegrator integrator(temp, 2.0, 0.01);
   StandardMMForceField* forceField = new StandardMMForceField(numAtoms, 0, 0, 0, 0);
   for (int i = 0; i < numAtoms; ++i) {
       system.setAtomMass(i, 10.0);
       forceField->setAtomParameters(i, (i%2 == 0 ? 0.2 : -0.2), 0.5, 5.0);
   }
   for (int i = 0; i < numAtoms-1; ++i)
       system.setConstraintParameters(i, i, i+1, 1.0);
   system.addForce(forceField);
   OpenMMContext context(system, integrator, platform);
   vector<Vec3> positions(numAtoms);
   vector<Vec3> velocities(numAtoms);
   init_gen_rand(0);
   for (int i = 0; i < numAtoms; ++i) {
       positions[i] = Vec3(i/2, (i+1)/2, 0);
       velocities[i] = Vec3(genrand_real2()-0.5, genrand_real2()-0.5, genrand_real2()-0.5);
   }
   context.setPositions(positions);
   context.setVelocities(velocities);
   
   // Simulate it and see whether the constraints remain satisfied.
   
   for (int i = 0; i < 1000; ++i) {
       State state = context.getState(State::Positions);
       for (int j = 0; j < numAtoms-1; ++j) {
           Vec3 p1 = state.getPositions()[j];
           Vec3 p2 = state.getPositions()[j+1];
           double dist = std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
           ASSERT_EQUAL_TOL(1.0, dist, 2e-4);
       }
       integrator.step(1);
   }
}

int main( ){

   (void) fflush( stdout );
   (void) fflush( stderr );
   try {
//      testBrookBondedHarmonicBond( );

//      testBrookNonBonded( );

/*
      testBrookStreams( 50 );
      testBrookStreams( 63 );
      testBrookStreams( 64 );
      testBrookStreams( 65 );

      testBrookBonds();
      testBrookAngles();
      testBrookPeriodicTorsions();
//      testBrookRBTorsions();
      testBrookCoulomb();
      testBrookLJ();

      testBrookExclusionsAnd14();
*/
      testLangevinSingleBond();
      //testForce(platform);
    }  catch( const exception& e ){
      cout << "exception: " << e.what() << endl;
      return 1;
   }   
   cout << "Done" << endl;

/*
   try {
      //testWriteRead();
      testBrookBonded();
   } catch( const exception& e ){
      cout << "exception: " << e.what() << endl;
      return 1;
   }
   cout << "Done" << endl;
*/
   return 0;
}
