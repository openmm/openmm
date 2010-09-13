/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "AmoebaTinkerParameterFile.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/internal/ContextImpl.h"
#include "kernels/amoebaGpuTypes.h"
#include "AmoebaCudaData.h"
#include "openmm/LocalEnergyMinimizer.h"
#include "../../../../../platforms//reference/src//SimTKUtilities/SimTKOpenMMUtilities.h"

#include <exception>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

extern int isNanOrInfinity( double number );

using namespace std;


extern "C" void* getAmoebaCudaData( ContextImpl& context );

/**---------------------------------------------------------------------------------------

    Replacement of sorts for strtok()
    Used to parse parameter file lines

    @param lineBuffer           string to tokenize
    @param delimiter            token delimter

    @return number of args

    --------------------------------------------------------------------------------------- */

static char* strsepLocal( char** lineBuffer, const char* delimiter ){

    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "strsepLocal";

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
       c     = *s++;
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

    Tokenize a string

    @param lineBuffer           string to tokenize
    @param tokenArray           upon return vector of tokens
    @param delimiter            token delimter

    @return number of tokens

    --------------------------------------------------------------------------------------- */

int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter ){

    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "tokenizeString";

    // ---------------------------------------------------------------------------------------

    char *ptr_c;
    while( (ptr_c = strsepLocal( &lineBuffer, delimiter.c_str() )) != NULL ){
        if( *ptr_c ){
            tokenArray.push_back( std::string( ptr_c ) );
        }
    }
    return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

    Open file

    @param fileName             file name
    @param mode                 file mode: "r", "w", "a"
    @param log                  optional logging file reference

    @return file pttr or NULL if file not opened

    --------------------------------------------------------------------------------------- */

static FILE* openFile( const std::string& fileName, const std::string& mode, FILE* log ){

    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "openFile";

    // ---------------------------------------------------------------------------------------

    FILE* filePtr;

#ifdef _MSC_VER
    fopen_s( &filePtr, fileName.c_str(), mode.c_str() );
#else
    filePtr = fopen( fileName.c_str(), mode.c_str() );
#endif

    if( log ){
        (void) fprintf( log, "openFile: file=<%s> %sopened w/ mode=%s.\n", fileName.c_str(), (filePtr == NULL ? "not " : ""), mode.c_str() );
        (void) fflush( log );
    }
    return filePtr;
}

/**---------------------------------------------------------------------------------------

    Tokenize a line into strings

    @param line                 line to tokenize
    @param tokenArray           upon return vector of tokens
    @param delimiter            token delimter

    @return number of tokens

    --------------------------------------------------------------------------------------- */

int tokenizeStringFromLineString( std::string& line, StringVector& tokenArray, const std::string delimiter ){

    // ---------------------------------------------------------------------------------------

    // static const std::string methodName = "tokenizeStringFromLineString";

    // ---------------------------------------------------------------------------------------

    char buffer[4096];
    (void) strcpy( buffer, line.c_str() );
    return tokenizeString( buffer, tokenArray, delimiter );
}

/**---------------------------------------------------------------------------------------
 *
 * Set string field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setStringFromMap( MapStringString& argumentMap, std::string fieldToCheck, std::string& fieldToSet ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName             = "setStringFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = (*check).second; 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set int field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setIntFromMap( MapStringString& argumentMap, std::string fieldToCheck, int& fieldToSet ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName             = "setIntFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atoi( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------

 * Set float field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setFloatFromMap( MapStringString& argumentMap, std::string fieldToCheck, float& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setFloatFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = static_cast<float>(atof( (*check).second.c_str() )); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set double field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setDoubleFromMap( MapStringString& argumentMap, std::string fieldToCheck, double& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setDoubleFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atof( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------

    Read a line from a file and tokenize into an array of strings

    @param filePtr              file to read from
    @param tokens               array of token strings
    @param lineCount            line count
    @param log                  optional file ptr for logging

    @return ptr to string containing line

    --------------------------------------------------------------------------------------- */

static int readLine( FILE* filePtr, StringVector& tokens, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "readLine";
    
    std::string delimiter                    = " \r\n";
    const int bufferSize                     = 4096;
    char buffer[bufferSize];

// ---------------------------------------------------------------------------------------

    char* isNotEof = fgets( buffer, bufferSize, filePtr );
    if( isNotEof ){
       (*lineCount)++;
       tokenizeString( buffer, tokens, delimiter );
       return 1;
    } else {
       return 0;
    }

}

/**---------------------------------------------------------------------------------------

    Read a file

    @param fileName             file name
    @param fileContents         output file contents
    @param log                  log

    @return 1 if file not opened; else return 0

    --------------------------------------------------------------------------------------- */

static int readFile( std::string fileName, StringVectorVector& fileContents, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readFile";

// ---------------------------------------------------------------------------------------

     fileContents.resize(0);

     // open file

     FILE* filePtr = openFile( fileName, "r", log );
     if( filePtr == NULL ){
        if( log ){
            (void) fprintf( log, "%s: file=<%s> not found.\n", methodName.c_str(), fileName.c_str() );
            (void) fflush( log );
        }
        return 1;
     } else if( log ){
         (void) fprintf( log, "%s: file=<%s> found.\n", methodName.c_str(), fileName.c_str() );
         (void) fflush( log );
     }

     // read contents

     StringVector firstLine;
     int lineCount  = 0;
     int isNotEof   = readLine( filePtr, firstLine, &lineCount, log );
     fileContents.push_back( firstLine );
     while( isNotEof ){
         StringVector lineTokens;
         isNotEof = readLine( filePtr, lineTokens, &lineCount, log );
         fileContents.push_back( lineTokens );
     }
     (void) fclose( filePtr );

     return 0;
}

/**---------------------------------------------------------------------------------------

    Read vector of double vectors

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param vectorOfVectors      output of vector of vectors
    @param lineCount            used to track line entries read from parameter file
    @param typeName             id of entries being read
    @param log                  log file pointer -- may be NULL

    @return number of entries read

    --------------------------------------------------------------------------------------- */

static int readVectorOfDoubleVectors( FILE* filePtr, const StringVector& tokens, std::vector< std::vector<double> >& vectorOfVectors, 
                                       int* lineCount, const std::string& typeName, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readVectorOfDoubleVectors";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no %s entries?\n", methodName.c_str(), typeName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    int numberToRead = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of %s to read: %d\n", methodName.c_str(), typeName.c_str(), numberToRead );
       (void) fflush( log );
    }

    for( int ii = 0; ii < numberToRead; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 1 ){
          int index            = atoi( lineTokens[0].c_str() );
          std::vector<double> nextEntry;
          for( unsigned int jj = 1; jj < lineTokens.size(); jj++ ){
              double value = atof( lineTokens[jj].c_str() );
              nextEntry.push_back( value );
          }
          vectorOfVectors.push_back( nextEntry );
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s %s tokens incomplete at line=%d\n", methodName.c_str(), typeName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint = MAX_PRINT;
       unsigned int   arraySize           = vectorOfVectors.size();
       (void) fprintf( log, "%s: sample of %s size=%u\n", methodName.c_str(), typeName.c_str(), arraySize );
       for( unsigned int ii = 0; ii < vectorOfVectors.size(); ii++ ){
          (void) fprintf( log, "%6u [", ii );
          for( unsigned int jj = 0; jj < vectorOfVectors[ii].size(); jj++ ){
             (void) fprintf( log, "%15.7e ", vectorOfVectors[ii][jj] );
          }
          (void) fprintf( log, "]\n" );

          // skip to end

          if( ii == maxPrint && (arraySize - maxPrint) > ii ){
             ii = arraySize - maxPrint - 1;
          } 
       }
    }

    return static_cast<int>(vectorOfVectors.size());
}

/**---------------------------------------------------------------------------------------

    Read vector of double vectors

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param vectorOfVectors      output of vector of vectors
    @param lineCount            used to track line entries read from parameter file
    @param typeName             id of entries being read
    @param log                  log file pointer -- may be NULL

    @return number of entries read

    --------------------------------------------------------------------------------------- */

static int readVectorOfIntVectors( FILE* filePtr, const StringVector& tokens, std::vector< std::vector<int> >& vectorOfVectors, 
                                   int* lineCount, std::string typeName, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readVectorOfIntVectors";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no %s entries?\n", methodName.c_str(), typeName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    int numberToRead = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of %s to read: %d\n", methodName.c_str(), typeName.c_str(), numberToRead );
       (void) fflush( log );
    }
    for( int ii = 0; ii < numberToRead; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 1 ){
          int index            = atoi( lineTokens[0].c_str() );
          std::vector<int> nextEntry;
          for( unsigned int jj = 1; jj < lineTokens.size(); jj++ ){
              int value = atoi( lineTokens[jj].c_str() );
              nextEntry.push_back( value );
          }
          vectorOfVectors.push_back( nextEntry );
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s %s tokens incomplete at line=%d\n", methodName.c_str(), typeName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint = MAX_PRINT;
       unsigned int   arraySize           = vectorOfVectors.size();
       (void) fprintf( log, "%s: sample of %s size=%u\n", methodName.c_str(), typeName.c_str(), arraySize );
       for( unsigned int ii = 0; ii < vectorOfVectors.size(); ii++ ){
          (void) fprintf( log, "%6u [", ii );
          for( unsigned int jj = 0; jj < vectorOfVectors[ii].size(); jj++ ){
             (void) fprintf( log, "%6d ", vectorOfVectors[ii][jj] );
          }
          (void) fprintf( log, "]\n" );

          // skip to end

          if( ii == maxPrint && (arraySize - maxPrint) > ii ){
             ii = arraySize - maxPrint - 1;
          } 
       }
    }

    return static_cast<int>(vectorOfVectors.size());
}

/**---------------------------------------------------------------------------------------

    Read vector of ints

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
                                line format: ... <numberOfTokens> <int value1> <int value2> ...<int valueN>,
                                where N=<numberOfTokens> 
    @param numberTokenIndex     index of entry <numberOfTokens> in token array
    @param intVector            output vector of ints
    @param lineCount            used to track line entries read from parameter file
    @param typeName             id of entries being read
    @param log                  log file pointer -- may be NULL

    @return number of entries read

    --------------------------------------------------------------------------------------- */

static int readIntVector( FILE* filePtr, const StringVector& tokens, int numberTokenIndex,
                           std::vector<int>& intVector, int* lineCount,
                           std::string typeName, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readIntVector";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no %s entries?\n", methodName.c_str(), typeName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }
 
    int numberToRead = atoi( tokens[numberTokenIndex].c_str() );
    intVector.resize( numberToRead );
    if( log ){
        (void) fprintf( log, "%s number of %s to read: %d\n", methodName.c_str(), typeName.c_str(), numberToRead );
        (void) fflush( log );
    }
    int startIndex = numberTokenIndex+1;
    for( int ii = startIndex; ii < numberToRead + startIndex; ii++ ){
        intVector[ii-3] = atoi( tokens[ii].c_str() );
    }
 
    return static_cast<int>(intVector.size());
}

/**---------------------------------------------------------------------------------------

    Read particles count

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of parameters read

    --------------------------------------------------------------------------------------- */

static int readParticles( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readParticles";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no particles number entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    int numberOfParticles = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s particles=%d\n", methodName.c_str(), numberOfParticles );
    }

    return numberOfParticles;
}

/**---------------------------------------------------------------------------------------

    Read particle masses

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of masses read

    --------------------------------------------------------------------------------------- */

static int readMasses( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readMasses";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no particle masses?\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    int numberOfParticles = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s particle masses=%d\n", methodName.c_str(), numberOfParticles );
    }
    for( int ii = 0; ii < numberOfParticles; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() >= 1 ){
          int tokenIndex       = 0;
          int index            = atoi( lineTokens[tokenIndex++].c_str() );
          double mass          = atof( lineTokens[tokenIndex++].c_str() );
          system.addParticle( mass );
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s particle tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint   = MAX_PRINT;
       unsigned int arraySize               = static_cast<unsigned int>(system.getNumParticles());
       (void) fprintf( log, "%s: sample of masses\n", methodName.c_str() );
       for( unsigned int ii = 0; ii < arraySize; ii++ ){
          (void) fprintf( log, "%6u %15.7e \n", ii, system.getParticleMass( ii ) );

          // skip to end

          if( ii == maxPrint && (arraySize - maxPrint) > ii ){
             ii = arraySize - maxPrint - 1;
          } 
       }
    }

    return system.getNumParticles();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba harmonic bond parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of bonds

    --------------------------------------------------------------------------------------- */

static int readAmoebaHarmonicBondParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                             System& system, int useOpenMMUnits,
                                             MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaHarmonicBondParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no bonds number entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    AmoebaHarmonicBondForce* bondForce = new AmoebaHarmonicBondForce();
    MapStringIntI forceActive          = forceMap.find( AMOEBA_HARMONIC_BOND_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( bondForce );
        if( log ){
            (void) fprintf( log, "Amoeba harmonic bond force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba harmonic bond force is not being included.\n" );
    }

    int numberOfBonds            = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of HarmonicBondForce terms=%d\n", methodName.c_str(), numberOfBonds );
    }
    for( int ii = 0; ii < numberOfBonds; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 4 ){
          int tokenIndex       = 0;
          int index            = atoi( lineTokens[tokenIndex++].c_str() );
          int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
          int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
          double length        = atof( lineTokens[tokenIndex++].c_str() );
          double k             = atof( lineTokens[tokenIndex++].c_str() );
          bondForce->addBond( particle1, particle2, length, k );
       } else {
          (void) fprintf( log, "%s AmoebaHarmonicBondForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          exit(-1);
       }
    }

    // get cubic and quartic factors

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 2 ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
 
            std::string field       = tokens[0];
            if( field == "#" ){
  
               // skip
               if( log ){
                     (void) fprintf( log, "skip <%s>\n", field.c_str());
               }
  
           } else if( field == "AmoebaHarmonicBondCubic" ){
               double cubicParameter = atof( tokens[1].c_str() );
               bondForce->setAmoebaGlobalHarmonicBondCubic( cubicParameter );
               hits++;
           } else if( field == "AmoebaHarmonicBondQuartic" ){
               double quarticParameter = atof( tokens[1].c_str() );
               bondForce->setAmoebaGlobalHarmonicBondQuartic( quarticParameter );
               hits++;
           }
        }
    }

    // convert units to kJ-nm from kCal-Angstrom?

    if( useOpenMMUnits ){

        double cubic         = bondForce->getAmoebaGlobalHarmonicBondCubic()/AngstromToNm;
        double quartic       = bondForce->getAmoebaGlobalHarmonicBondQuartic()/(AngstromToNm*AngstromToNm);

        bondForce->setAmoebaGlobalHarmonicBondCubic( cubic );
        bondForce->setAmoebaGlobalHarmonicBondQuartic( quartic );

        // scale equilibrium bond lengths/force prefactor k

        for( int ii = 0; ii < bondForce->getNumBonds(); ii++ ){
            int particle1, particle2;
            double length, k;
            bondForce->getBondParameters( ii, particle1, particle2, length, k );
            length           *= AngstromToNm;
            k                *= CalToJoule/(AngstromToNm*AngstromToNm);
            bondForce->setBondParameters( ii, particle1, particle2, length, k ); 
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(bondForce->getNumBonds());
        (void) fprintf( log, "%s: %u sample of AmoebaHarmonicBondForce parameters in %s units; cubic=%15.7e quartic=%15.7e\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba"),
                        bondForce->getAmoebaGlobalHarmonicBondCubic(), bondForce->getAmoebaGlobalHarmonicBondQuartic() );
        for( unsigned int ii = 0; ii < arraySize; ii++ ){
            int particle1, particle2;
            double length, k;
            bondForce->getBondParameters( ii, particle1, particle2, length, k );
            (void) fprintf( log, "%8d %8d %8d %15.7e %15.7e\n", ii, particle1, particle2, length, k );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
               ii = arraySize - maxPrint - 1;
            } 
        }
    }

    return bondForce->getNumBonds();
}

/**---------------------------------------------------------------------------------------
    Read Amoeba harmonic angle parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of angles

    --------------------------------------------------------------------------------------- */

static int readAmoebaHarmonicAngleParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                              System& system, int useOpenMMUnits,
                                              MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaHarmonicAngleParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no angles number entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    AmoebaHarmonicAngleForce* angleForce = new AmoebaHarmonicAngleForce();
    MapStringIntI forceActive            = forceMap.find( AMOEBA_HARMONIC_ANGLE_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( angleForce );
        if( log ){
            (void) fprintf( log, "Amoeba harmonic angle force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba harmonic angle force is not being included.\n" );
    }

    int numberOfAngles            = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of HarmonicAngleForce terms=%d\n", methodName.c_str(), numberOfAngles );
    }
    for( int ii = 0; ii < numberOfAngles; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 5 ){
          int tokenIndex       = 0;
          int index            = atoi( lineTokens[tokenIndex++].c_str() );
          int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
          int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
          int particle3        = atoi( lineTokens[tokenIndex++].c_str() );
          //double angle   = DegreesToRadians*atof( lineTokens[tokenIndex++].c_str() );
          double angle         = atof( lineTokens[tokenIndex++].c_str() );
          double k             = atof( lineTokens[tokenIndex++].c_str() );
          angleForce->addAngle( particle1, particle2, particle3, angle, k );
       } else {
          (void) fprintf( log, "%s AmoebaHarmonicAngleForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          exit(-1);
       }
    }

    // get cubic, quartic, pentic, sextic factors

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 4 ){
       StringVector tokens;
       isNotEof = readLine( filePtr, tokens, lineCount, log );
       if( isNotEof && tokens.size() > 0 ){

          std::string field       = tokens[0];
          if( field == "#" ){

             // skip
             if( log ){
                   (void) fprintf( log, "skip <%s>\n", field.c_str());
             }

          } else if( field == "AmoebaHarmonicAngleCubic" ){
             angleForce->setAmoebaGlobalHarmonicAngleCubic( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaHarmonicAngleQuartic" ){
             angleForce->setAmoebaGlobalHarmonicAngleQuartic( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaHarmonicAnglePentic" ){
             angleForce->setAmoebaGlobalHarmonicAnglePentic( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaHarmonicAngleSextic" ){
             angleForce->setAmoebaGlobalHarmonicAngleSextic( atof( tokens[1].c_str() ) );
             hits++;
          }
       }
    }

    // convert units to kJ-nm from kCal-Angstrom?

    if( useOpenMMUnits ){
        for( int ii = 0; ii < angleForce->getNumAngles(); ii++ ){
            int particle1, particle2, particle3;
            double length, k;
            angleForce->getAngleParameters( ii, particle1, particle2, particle3, length, k ); 
            k                *= CalToJoule;
            angleForce->setAngleParameters( ii, particle1, particle2, particle3, length, k ); 
        }
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint   = MAX_PRINT;
       unsigned int arraySize               = static_cast<unsigned int>(angleForce->getNumAngles());
       (void) fprintf( log, "%s: %u sample of AmoebaHarmonicAngleForce parameters in %s units; cubic=%15.7e quartic=%15.7e pentic=%15.7e sextic=%15.7e\n",
                       methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba"), 
                       angleForce->getAmoebaGlobalHarmonicAngleCubic(), 
                       angleForce->getAmoebaGlobalHarmonicAngleQuartic(), angleForce->getAmoebaGlobalHarmonicAnglePentic(),
                       angleForce->getAmoebaGlobalHarmonicAngleSextic() );

       for( unsigned int ii = 0; ii < arraySize; ii++ ){
          int particle1, particle2, particle3;
          double length, k;
          angleForce->getAngleParameters( ii, particle1, particle2, particle3, length, k ); 
          (void) fprintf( log, "%8d %8d %8d %8d %15.7e %15.7e\n",
                          ii, particle1, particle2, particle3, length, k );

          // skip to end

          if( ii == maxPrint && (arraySize - maxPrint) > ii ){
             ii = arraySize - maxPrint - 1;
          } 
       }
    }
    return angleForce->getNumAngles();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba harmonic angle parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of angles

    --------------------------------------------------------------------------------------- */

static int readAmoebaHarmonicInPlaneAngleParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                                      System& system, int useOpenMMUnits,
                                                      MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaHarmonicInPlaneAngleParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no angles number entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    AmoebaHarmonicInPlaneAngleForce* angleForce = new AmoebaHarmonicInPlaneAngleForce();
    MapStringIntI forceActive                   = forceMap.find( AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( angleForce );
        if( log ){
            (void) fprintf( log, "Amoeba harmonic in-plane angle force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba harmonic in-plane angle force is not being included.\n" );
    }

    int numberOfAngles            = atoi( tokens[1].c_str() );
    if( log ){
        (void) fprintf( log, "%s number of HarmonicInPlaneAngleForce terms=%d\n", methodName.c_str(), numberOfAngles );
    }
    for( int ii = 0; ii < numberOfAngles; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 6 ){
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle3        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle4        = atoi( lineTokens[tokenIndex++].c_str() );
            double angle         = atof( lineTokens[tokenIndex++].c_str() );
            double k             = atof( lineTokens[tokenIndex++].c_str() );
            angleForce->addAngle( particle1, particle2, particle3, particle4, angle, k );
        } else {
            (void) fprintf( log, "%s AmoebaHarmonicInPlaneAngleForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            exit(-1);
        }
    }

    // get cubic, quartic, pentic, sextic factors

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 4 ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
            std::string field       = tokens[0];
            if( field == "#" ){
  
               // skip
               if( log ){
                   (void) fprintf( log, "skip <%s>\n", field.c_str());
               }
  
            } else if( field == "AmoebaHarmonicInPlaneAngleCubic" ){
                angleForce->setAmoebaGlobalHarmonicInPlaneAngleCubic( atof( tokens[1].c_str() ) );
                hits++;
            } else if( field == "AmoebaHarmonicInPlaneAngleQuartic" ){
                angleForce->setAmoebaGlobalHarmonicInPlaneAngleQuartic( atof( tokens[1].c_str() ) );
                hits++;
            } else if( field == "AmoebaHarmonicInPlaneAnglePentic" ){
                angleForce->setAmoebaGlobalHarmonicInPlaneAnglePentic( atof( tokens[1].c_str() ) );
                hits++;
            } else if( field == "AmoebaHarmonicInPlaneAngleSextic" ){
                angleForce->setAmoebaGlobalHarmonicInPlaneAngleSextic( atof( tokens[1].c_str() ) );
                hits++;
            }
        }
    }

    // convert units to kJ-nm from kCal-Angstrom?

    if( useOpenMMUnits ){
        for( int ii = 0; ii < angleForce->getNumAngles(); ii++ ){
            int particle1, particle2, particle3, particle4;
            double length, k;
            angleForce->getAngleParameters( ii, particle1, particle2, particle3, particle4, length, k );
            k                *= CalToJoule;
            angleForce->setAngleParameters( ii, particle1, particle2, particle3, particle4, length, k );
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(angleForce->getNumAngles());
        (void) fprintf( log, "%s: %u sample of AmoebaHarmonicInPlaneAngleForce parameters in %s units; cubic=%15.7e quartic=%15.7e pentic=%15.7e sextic=%15.7e\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba"),
                        angleForce->getAmoebaGlobalHarmonicInPlaneAngleCubic(), 
                        angleForce->getAmoebaGlobalHarmonicInPlaneAngleQuartic(),
                        angleForce->getAmoebaGlobalHarmonicInPlaneAnglePentic(),
                        angleForce->getAmoebaGlobalHarmonicInPlaneAngleSextic() );
 
        for( unsigned int ii = 0; ii < arraySize; ii++ ){
            int particle1, particle2, particle3, particle4;
            double length, k;
            angleForce->getAngleParameters( ii, particle1, particle2, particle3, particle4, length, k );
            (void) fprintf( log, "%8d %8d %8d %8d %8d %15.7e %15.7e\n",
                            ii, particle1, particle2, particle3, particle4, length, k );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
    }

    return angleForce->getNumAngles();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba harmonic angle parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of angles

    --------------------------------------------------------------------------------------- */

static int readAmoebaTorsionParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                        System& system, int useOpenMMUnits,
                                        MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaTorsionParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no torsions number entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    AmoebaTorsionForce* torsionForce   = new AmoebaTorsionForce();
    MapStringIntI forceActive          = forceMap.find( AMOEBA_TORSION_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( torsionForce );
        if( log ){
            (void) fprintf( log, "Amoeba torsion force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba torsion force is not being included.\n" );
    }

    int numberOfTorsions            = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of TorsionForce terms=%d\n", methodName.c_str(), numberOfTorsions );
    }
    for( int ii = 0; ii < numberOfTorsions; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 10 ){
          std::vector<double> torsion1;
          std::vector<double> torsion2;
          std::vector<double> torsion3;
          int index            = 0;
          int torsionIndex     = atoi( lineTokens[index++].c_str() );
          int particle1        = atoi( lineTokens[index++].c_str() );
          int particle2        = atoi( lineTokens[index++].c_str() );
          int particle3        = atoi( lineTokens[index++].c_str() );
          int particle4        = atoi( lineTokens[index++].c_str() );

          double a1            = atof( lineTokens[index++].c_str() );
          double p1            = DegreesToRadians*atof( lineTokens[index++].c_str() );
          torsion1.push_back( a1 );
          torsion1.push_back( p1 );

          double a2            = atof( lineTokens[index++].c_str() );
          double p2            = DegreesToRadians*atof( lineTokens[index++].c_str() );
          torsion2.push_back( a2 );
          torsion2.push_back( p2 );


          double a3            = atof( lineTokens[index++].c_str() );
          double p3            = DegreesToRadians*atof( lineTokens[index++].c_str() );
          torsion3.push_back( a3 );
          torsion3.push_back( p3 );

          torsionForce->addTorsion( particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3 );
       } else {
          (void) fprintf( log, "%s AmoebaTorsionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          exit(-1);
       }
    }

    // convert units to kJ-nm from kCal-Angstrom?

    if( useOpenMMUnits ){
        for( int ii = 0; ii < torsionForce->getNumTorsions(); ii++ ){
            int particle1, particle2, particle3, particle4;
            std::vector<double> torsion1;
            std::vector<double> torsion2;
            std::vector<double> torsion3;
            torsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3 ); 
            torsion1[0] *= CalToJoule;
            torsion2[0] *= CalToJoule;
            torsion3[0] *= CalToJoule;
            torsionForce->setTorsionParameters( ii, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3 ); 
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(torsionForce->getNumTorsions());
        (void) fprintf( log, "%s: %u sample of AmoebaTorsionForce parameters in %s units.\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba") );
    
        for( unsigned int ii = 0; ii < arraySize; ii++ ){
            int particle1, particle2, particle3, particle4;
            std::vector<double> torsion1;
            std::vector<double> torsion2;
            std::vector<double> torsion3;
            torsionForce->getTorsionParameters( ii, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3 ); 
            (void) fprintf( log, "%8d %8d %8d %8d %8d [%15.7e %15.7e] [%15.7e %15.7e] [%15.7e %15.7e]\n",
                            ii, particle1, particle2, particle3, particle4, 
                            torsion1[0], torsion1[1]/DegreesToRadians, torsion2[0], torsion2[1]/DegreesToRadians, torsion3[0], torsion3[1]/DegreesToRadians );
    
            // skip to end
    
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
               ii = arraySize - maxPrint - 1;
            } 
        }
    } 

    return torsionForce->getNumTorsions();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba pi torsion parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of angles

    --------------------------------------------------------------------------------------- */

static int readAmoebaPiTorsionParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                          System& system, int useOpenMMUnits,
                                          MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaPiTorsionParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no pi torsions number entry???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    AmoebaPiTorsionForce* piTorsionForce = new AmoebaPiTorsionForce();
    MapStringIntI forceActive            = forceMap.find( AMOEBA_PI_TORSION_FORCE );

    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( piTorsionForce );
        if( log ){
            (void) fprintf( log, "Amoeba pi torsion force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba pi torsion force is not being included.\n" );
    }

    int numberOfPiTorsions            = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of PiTorsionForce terms=%d\n", methodName.c_str(), numberOfPiTorsions );
    }
    for( int ii = 0; ii < numberOfPiTorsions; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 7 ){
           std::vector<double> torsionK;
           int index            = 0;
           int torsionIndex     = atoi( lineTokens[index++].c_str() );
           int particle1        = atoi( lineTokens[index++].c_str() );
           int particle2        = atoi( lineTokens[index++].c_str() );
           int particle3        = atoi( lineTokens[index++].c_str() );
           int particle4        = atoi( lineTokens[index++].c_str() );
           int particle5        = atoi( lineTokens[index++].c_str() );
           int particle6        = atoi( lineTokens[index++].c_str() );
           double k             = atof( lineTokens[index++].c_str() );
 
           piTorsionForce->addPiTorsion( particle1, particle2, particle3, particle4, particle5, particle6, k );
        } else {
           (void) fprintf( log, "%s AmoebaPiTorsionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
           exit(-1);
        }
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){
        for( int ii = 0; ii < piTorsionForce->getNumPiTorsions(); ii++ ){
            int particle1, particle2, particle3, particle4, particle5, particle6;
            double torsionK;
            piTorsionForce->getPiTorsionParameters( ii, particle1, particle2, particle3, particle4, particle5, particle6, torsionK ); 
            torsionK *= CalToJoule;
            piTorsionForce->setPiTorsionParameters( ii, particle1, particle2, particle3, particle4, particle5, particle6, torsionK ); 
        }
    } 

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(piTorsionForce->getNumPiTorsions());
        (void) fprintf( log, "%s: %u sample of AmoebaPiTorsionForce parameters in %s units.\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba") );
 
        for( unsigned int ii = 0; ii < arraySize; ii++ ){
            int particle1, particle2, particle3, particle4, particle5, particle6;
            double torsionK;
            piTorsionForce->getPiTorsionParameters( ii, particle1, particle2, particle3, particle4, particle5, particle6, torsionK ); 
            (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %8d k=%15.7e\n",
                            ii, particle1, particle2, particle3, particle4,  particle5, particle6, torsionK );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
               ii = arraySize - maxPrint - 1;
            } 
        }
    } 

    return piTorsionForce->getNumPiTorsions();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba stretchBend parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of stretchBends

    --------------------------------------------------------------------------------------- */

static int readAmoebaStretchBendParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                            System& system, int useOpenMMUnits,
                                            MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaStretchBendParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no stretchBends number entry???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    // create force

    AmoebaStretchBendForce* stretchBendForce = new AmoebaStretchBendForce();
    MapStringIntI forceActive                = forceMap.find( AMOEBA_STRETCH_BEND_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( stretchBendForce );
        if( log ){
            (void) fprintf( log, "Amoeba stretchBend force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba stretchBend force is not being included.\n" );
    }

    // load in parameters

    int numberOfStretchBends            = atoi( tokens[1].c_str() );
    if( log ){
        (void) fprintf( log, "%s number of StretchBendForce terms=%d\n", methodName.c_str(), numberOfStretchBends );
    }
    for( int ii = 0; ii < numberOfStretchBends; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 7 ){
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle3        = atoi( lineTokens[tokenIndex++].c_str() );
            double lengthAB      = atof( lineTokens[tokenIndex++].c_str() );
            double lengthCB      = atof( lineTokens[tokenIndex++].c_str() );
            double angle         = atof( lineTokens[tokenIndex++].c_str() )*DegreesToRadians;
            double k             = atof( lineTokens[tokenIndex++].c_str() );
            stretchBendForce->addStretchBend( particle1, particle2, particle3, lengthAB, lengthCB, angle, k );
        } else {
            (void) fprintf( log, "%s AmoebaStretchBendForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            exit(-1);
        }
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){
        for( int ii = 0; ii < stretchBendForce->getNumStretchBends();  ii++ ){
            int particle1, particle2, particle3;
            double lengthAB, lengthCB, angle, k;
            stretchBendForce->getStretchBendParameters( ii, particle1, particle2, particle3, lengthAB, lengthCB, angle, k );
            lengthAB     *= AngstromToNm;
            lengthCB     *= AngstromToNm;
            k            *= CalToJoule/AngstromToNm;
            stretchBendForce->setStretchBendParameters( ii, particle1, particle2, particle3, lengthAB, lengthCB, angle, k );
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(stretchBendForce->getNumStretchBends());
        (void) fprintf( log, "%s: %u sample of AmoebaStretchBendForce parameters in %s units.\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba") );
        for( unsigned int ii = 0; ii < arraySize;  ii++ ){
            int particle1, particle2, particle3;
            double lengthAB, lengthCB, angle, k;
            stretchBendForce->getStretchBendParameters( ii, particle1, particle2, particle3, lengthAB, lengthCB, angle, k );
            (void) fprintf( log, "%8d %8d %8d %8d %15.7e %15.7e %15.7e %15.7e\n",
                            ii, particle1, particle2, particle3, lengthAB, lengthCB, angle/DegreesToRadians, k );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
               ii = arraySize - maxPrint - 1;
            } 
        }
    }

    return stretchBendForce->getNumStretchBends();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba stretchBend parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of stretchBends

    --------------------------------------------------------------------------------------- */

static int readAmoebaOutOfPlaneBendParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                               System& system, int useOpenMMUnits,
                                               MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaOutOfPlaneBendParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no outOfPlaneBends number entry???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    // create force

    AmoebaOutOfPlaneBendForce* outOfPlaneBendForce = new AmoebaOutOfPlaneBendForce();
    MapStringIntI forceActive                      = forceMap.find( AMOEBA_OUT_OF_PLANE_BEND_FORCE );

    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( outOfPlaneBendForce );
        if( log ){
            (void) fprintf( log, "Amoeba outOfPlaneBend force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba outOfPlaneBend force is not being included.\n" );
    }

    // load in parameters

    int numberOfOutOfPlaneBends            = atoi( tokens[1].c_str() );
    if( log ){
        (void) fprintf( log, "%s number of OutOfPlaneBendForce terms=%d\n", methodName.c_str(), numberOfOutOfPlaneBends );
    }
    for( int ii = 0; ii < numberOfOutOfPlaneBends; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 5 ){
           int tokenIndex       = 0;
           int index            = atoi( lineTokens[tokenIndex++].c_str() );
           int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
           int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
           int particle3        = atoi( lineTokens[tokenIndex++].c_str() );
           int particle4        = atoi( lineTokens[tokenIndex++].c_str() );
           double k             = atof( lineTokens[tokenIndex++].c_str() );
           outOfPlaneBendForce->addOutOfPlaneBend( particle1, particle2, particle3, particle4, k );
        } else {
           (void) fprintf( log, "%s AmoebaOutOfPlaneBendForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
           exit(-1);
        }
    }

    // get cubic, quartic, pentic, sextic factors

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 4 ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
           std::string field       = tokens[0];
           if( field == "#" ){
 
              // skip
              if( log ){
                  (void) fprintf( log, "skip <%s>\n", field.c_str());
              }
 
           } else if( field == "AmoebaOutOfPlaneBendCubic" ){
               outOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendCubic( atof( tokens[1].c_str() ) );
               hits++;
           } else if( field == "AmoebaOutOfPlaneBendQuartic" ){
               outOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendQuartic( atof( tokens[1].c_str() ) );
               hits++;
           } else if( field == "AmoebaOutOfPlaneBendPentic" ){
               outOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendPentic( atof( tokens[1].c_str() ) );
               hits++;
           } else if( field == "AmoebaOutOfPlaneBendSextic" ){
               outOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendSextic( atof( tokens[1].c_str() ) );
               hits++;
           }
        }
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){
        for( int ii = 0; ii < outOfPlaneBendForce->getNumOutOfPlaneBends();  ii++ ){
            int particle1, particle2, particle3, particle4;
            double k;
            outOfPlaneBendForce->getOutOfPlaneBendParameters( ii, particle1, particle2, particle3, particle4, k );
            //k *= CalToJoule/(AngstromToNm*AngstromToNm);
            k *= CalToJoule;
            outOfPlaneBendForce->setOutOfPlaneBendParameters( ii, particle1, particle2, particle3, particle4, k );
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(outOfPlaneBendForce->getNumOutOfPlaneBends());
        (void) fprintf( log, "%s: %u sample of AmoebaOutOfPlaneBendForce parameters in %s units.\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba") );
        (void) fprintf( log, "%s: %u sample of AmoebaOutOfPlaneBendForce parameters; cubic=%15.7e quartic=%15.7e pentic=%15.7e sextic=%15.7e\n",
                        methodName.c_str(), arraySize,
                        outOfPlaneBendForce->getAmoebaGlobalOutOfPlaneBendCubic(), 
                        outOfPlaneBendForce->getAmoebaGlobalOutOfPlaneBendQuartic(),
                        outOfPlaneBendForce->getAmoebaGlobalOutOfPlaneBendPentic(),
                        outOfPlaneBendForce->getAmoebaGlobalOutOfPlaneBendSextic() );
 
        for( unsigned int ii = 0; ii < arraySize;  ii++ ){
            int particle1, particle2, particle3, particle4;
            double k;
            outOfPlaneBendForce->getOutOfPlaneBendParameters( ii, particle1, particle2, particle3, particle4, k );
            (void) fprintf( log, "%8d %8d %8d %8d %8d %15.7e\n",
                            ii, particle1, particle2, particle3, particle4, k );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
    }

    return outOfPlaneBendForce->getNumOutOfPlaneBends();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba torsion-torsion parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of torsion-torsion parameters

    --------------------------------------------------------------------------------------- */

static int readAmoebaTorsionTorsionGrid( FILE* filePtr, int numX, int numY, TorsionTorsionGrid& grid, 
                                         int useOpenMMUnits,
                                         MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaTorsionTorsionForce";
    
// ---------------------------------------------------------------------------------------

    int gridCount = numX*numY;
    if( log ){
        (void) fprintf( log, "%s number of TorsionTorsion grid entries: %d %d %d\n", methodName.c_str(), numX, numY, gridCount);
    }

    grid.resize( numX );
    for( int ii = 0; ii < numX; ii++ ){
        grid[ii].resize( numY );
    }
    for( int ii = 0; ii < gridCount; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 8 ){
           int tokenIndex             = 0;
           int index                  = atoi( lineTokens[tokenIndex++].c_str() );
           int xIndex                 = atoi( lineTokens[tokenIndex++].c_str() );
           int yIndex                 = atoi( lineTokens[tokenIndex++].c_str() );
 
           std::vector<double> values;
           values.resize( 6 ); 
           int vIndex                 = 0;
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // xValue
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // yValue
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // fValue
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // dfdxValue
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // dfdyValue
           values[vIndex++]           = atof( lineTokens[tokenIndex++].c_str() ); // dfdxyValue
           if( useOpenMMUnits ){
               values[2] *= CalToJoule;
               values[3] *= CalToJoule;
               values[4] *= CalToJoule;
               values[5] *= CalToJoule;
           }
           grid[xIndex][yIndex]       = values;
       } else {
           (void) fprintf( log, "%s grid tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
           exit(-1);
       }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = 20;
        (void) fprintf( log, "%s: %dx%d sample of grid values; using %s units\n", methodName.c_str(), numX, numY,
                        (useOpenMMUnits ? "OpenMM" : "Amoeba") );
 
        int xI = 0;
        int yI = 0;
 
        // 'top' of grid
 
        for( int ii = 0; ii < gridCount;  ii++ ){
            std::vector<double> values = grid[xI][yI];
            (void) fprintf( log, "%4d %4d %4d ", ii, xI, yI );
            for( unsigned int jj = 0; jj < values.size(); jj++ ){
                (void) fprintf( log, "%15.7e ", values[jj] );
            }
            fprintf( log, "\n" );
  
            // increment or quit
  
            xI++; 
            if( xI == numX ){
               xI = 0;
               yI++;
               if( yI == 3 )ii = gridCount;
            }
        }
 
        // 'bottom' of grid
 
        xI = 0;
        yI = (gridCount/numX) - 3;
        for( int ii = 0; ii < gridCount;  ii++ ){
            std::vector<double> values = grid[xI][yI];
            (void) fprintf( log, "%4d %4d %4d ",
                            ii, xI, yI );
            for( unsigned int jj = 0; jj < values.size(); jj++ ){
                (void) fprintf( log, "%15.7e ", values[jj] );
            }
            fprintf( log, "\n" );
            xI++; 
            if( xI == numX ){
               xI = 0;
               yI++;
               if( yI == numY )ii = gridCount;
            }
        }
    }

    return 0;
}

/**---------------------------------------------------------------------------------------

    Read Amoeba torsion-torsion parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of torsion-torsion parameters

    --------------------------------------------------------------------------------------- */

static int readAmoebaTorsionTorsionParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                               System& system, int useOpenMMUnits,
                                               MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaTorsionTorsionParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no torsionTorsions number entry???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    // create force

    AmoebaTorsionTorsionForce* torsionTorsionForce = new AmoebaTorsionTorsionForce();
    MapStringIntI forceActive                      = forceMap.find( AMOEBA_TORSION_TORSION_FORCE );

    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( torsionTorsionForce );
        if( log ){
            (void) fprintf( log, "Amoeba torsionTorsion force is being included.\n" );
        }
    } else if( log ){
        (void) fprintf( log, "Amoeba torsionTorsion force is not being included.\n" );
    }

    // load in parameters

    int numberOfTorsionTorsions            = atoi( tokens[1].c_str() );
    if( log ){
        (void) fprintf( log, "%s number of TorsionTorsionForce terms=%d\n", methodName.c_str(), numberOfTorsionTorsions );
    }
    for( int ii = 0; ii < numberOfTorsionTorsions; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 7 ){
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            int particle1        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle2        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle3        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle4        = atoi( lineTokens[tokenIndex++].c_str() );
            int particle5        = atoi( lineTokens[tokenIndex++].c_str() );
            int chiralAtomIndex  = atoi( lineTokens[tokenIndex++].c_str() );
            int gridIndex        = atoi( lineTokens[tokenIndex++].c_str() );
            torsionTorsionForce->addTorsionTorsion( particle1, particle2, particle3, particle4, particle5, chiralAtomIndex, gridIndex );
        } else {
            (void) fprintf( log, "%s AmoebaTorsionTorsionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            exit(-1);
        }
    }

    // get grid

    int isNotEof                   = 1;
    int totalNumberOfGrids         = 1;
    int gridCount                  = 0;
    while( gridCount < totalNumberOfGrids ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
            std::string field       = tokens[0];
            if( field == "#" ){
  
               // skip
               if( log ){
                   (void) fprintf( log, "skip <%s>\n", field.c_str());
               }
  
            } else if( field == "AmoebaTorsionTorsionGrids" ){
                totalNumberOfGrids = atoi( tokens[1].c_str() );
            } else if( field == "AmoebaTorsionTorsionGridPoints" ){
                int gridIndex = atoi( tokens[1].c_str() );
                int numX      = atoi( tokens[2].c_str() );
                int numY      = atoi( tokens[3].c_str() );
                TorsionTorsionGrid torsionTorsionGrid;
                readAmoebaTorsionTorsionGrid( filePtr, numX, numY, torsionTorsionGrid, useOpenMMUnits, inputArgumentMap, lineCount, log );
                torsionTorsionForce->setTorsionTorsionGrid( gridIndex, torsionTorsionGrid );
                gridCount++;
            }
        }
    }

    // diagnostics

    if( log ){
        static const unsigned int maxPrint   = MAX_PRINT;
        unsigned int arraySize               = static_cast<unsigned int>(torsionTorsionForce->getNumTorsionTorsions());
        (void) fprintf( log, "%s: %u sample of AmoebaTorsionTorsionForce parameters\n",
                        methodName.c_str(), arraySize );
        (void) fprintf( log, "%s: %u sample of AmoebaTorsionTorsionForce parameters\n",
                        methodName.c_str(), arraySize );
 
        for( unsigned int ii = 0; ii < arraySize;  ii++ ){
            int particle1, particle2, particle3, particle4, particle5, chiralAtomIndex, gridIndex;
            torsionTorsionForce->getTorsionTorsionParameters( ii, particle1, particle2, particle3, particle4, particle5, chiralAtomIndex, gridIndex );
            (void) fprintf( log, "%8d %8d %8d %8d %8d %8d %8d %8d\n",
                            ii, particle1, particle2, particle3, particle4, particle5, chiralAtomIndex, gridIndex );
  
            // skip to end
  
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
    }

    return torsionTorsionForce->getNumTorsionTorsions();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba multipole parameters

    @param filePtr              file pointer to parameter file
    @param multipoleForce       AmoebaMultipoleForce reference
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of multipole parameters

    --------------------------------------------------------------------------------------- */

static void readAmoebaMultipoleCovalent( FILE* filePtr, AmoebaMultipoleForce* multipoleForce,
                                          int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

     static const std::string methodName      = "readAmoebaMultipoleCovalent";
    
// ---------------------------------------------------------------------------------------

     // load in parameters
 
     int numberOfMultipoles            = multipoleForce->getNumMultipoles();
     if( log ){
         (void) fprintf( log, "%s number of multipoles=%d\n", methodName.c_str(), numberOfMultipoles );
     }
     unsigned int numberOfCovalentTypes = AmoebaMultipoleForce::CovalentEnd;
     AmoebaMultipoleForce::CovalentType covalentTypes[AmoebaMultipoleForce::CovalentEnd] = {
               AmoebaMultipoleForce::Covalent12,
               AmoebaMultipoleForce::Covalent13,
               AmoebaMultipoleForce::Covalent14,
               AmoebaMultipoleForce::Covalent15,
               AmoebaMultipoleForce::PolarizationCovalent11,
               AmoebaMultipoleForce::PolarizationCovalent12,
               AmoebaMultipoleForce::PolarizationCovalent13,
               AmoebaMultipoleForce::PolarizationCovalent14 };
 
     int isNotEof = 1;
     for( int ii = 0; ii < numberOfMultipoles && isNotEof; ii++ ){
         for( unsigned int jj = 0; jj < numberOfCovalentTypes && isNotEof; jj++ ){
             StringVector lineTokens;
             isNotEof = readLine( filePtr, lineTokens, lineCount, log );
             std::vector<int> intVector;
             readIntVector( filePtr, lineTokens, 2, intVector, lineCount, lineTokens[0], (ii < 10 ? log: NULL) );
             multipoleForce->setCovalentMap( ii, covalentTypes[jj], intVector ); 
         }
     }
}

/**---------------------------------------------------------------------------------------

    Read Amoeba multipole parameters

    @param filePtr              file pointer to parameter file
    @param version              version id for file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of multipole parameters

    --------------------------------------------------------------------------------------- */

static int readAmoebaMultipoleParameters( FILE* filePtr, int version, MapStringInt& forceMap, const StringVector& tokens,
                                           System& system, int useOpenMMUnits, MapStringVectorOfVectors& supplementary,
                                           MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaMultipoleParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no multipoles number entry???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    // create force

    AmoebaMultipoleForce* multipoleForce = new AmoebaMultipoleForce();
    MapStringIntI forceActive            = forceMap.find( AMOEBA_MULTIPOLE_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( multipoleForce );
        if( log ){
            (void) fprintf( log, "Amoeba multipole force is being included.\n" );  (void) fflush( log );
        }    
    } else if( log ){
        (void) fprintf( log, "Amoeba multipole force is not being included.\n" ); (void) fflush( log );
    }

    // load in parameters
 
    int numberOfMultipoles            = atoi( tokens[1].c_str() );
    int usePme                        = 0;
    double aewald                     = 0.0;
    double cutoffDistance             = 0.0;
    double box[3]                     = { 10.0, 10.0, 10.0 };

    // usePme, aewald, cutoffDistance added w/ Version 1

    if( version > 0 ){
        usePme                        = atoi( tokens[2].c_str() );
        aewald                        = atof( tokens[3].c_str() );
        cutoffDistance                = atof( tokens[4].c_str() );
        box[0]                        = atof( tokens[5].c_str() );
        box[1]                        = atof( tokens[6].c_str() );
        box[2]                        = atof( tokens[7].c_str() );
    }
  
    if( usePme ){
        multipoleForce->setNonbondedMethod( AmoebaMultipoleForce::PME );
    } else {
        multipoleForce->setNonbondedMethod( AmoebaMultipoleForce::NoCutoff );
    }
    multipoleForce->setCutoffDistance( cutoffDistance );
    multipoleForce->setAEwald( aewald );
    system.setDefaultPeriodicBoxVectors( Vec3(box[0], 0.0, 0.0), Vec3(0.0, box[1], 0.0), Vec3(0.0, 0.0, box[2]) );
    if( log ){
        (void) fprintf( log, "%s number of MultipoleParameter terms=%d usePme=%d aewald=%15.7e cutoffDistance=%12.4f\n",
                        methodName.c_str(), numberOfMultipoles, usePme, aewald, cutoffDistance );
        (void) fflush( log );
    }
    for( int ii = 0; ii < numberOfMultipoles; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 7 ){
            std::vector<double> dipole;
            std::vector<double> quadrupole;
            dipole.resize( 3 );
            quadrupole.resize( 9 );
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            int axisType         = atoi( lineTokens[tokenIndex++].c_str() );
            int zAxis            = atoi( lineTokens[tokenIndex++].c_str() );
            int xAxis            = atoi( lineTokens[tokenIndex++].c_str() );
            double pdamp         = atof( lineTokens[tokenIndex++].c_str() );
            double tholeDamp     = atof( lineTokens[tokenIndex++].c_str() );
            double polarity      = atof( lineTokens[tokenIndex++].c_str() );
            double charge        = atof( lineTokens[tokenIndex++].c_str() );
            dipole[0]            = atof( lineTokens[tokenIndex++].c_str() );
            dipole[1]            = atof( lineTokens[tokenIndex++].c_str() );
            dipole[2]            = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[0]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[1]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[2]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[3]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[4]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[5]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[6]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[7]        = atof( lineTokens[tokenIndex++].c_str() );
            quadrupole[8]        = atof( lineTokens[tokenIndex++].c_str() );
            multipoleForce->addParticle( charge, dipole, quadrupole, axisType, zAxis, xAxis, tholeDamp, pdamp, polarity );
        } else {
            (void) fprintf( log, "%s AmoebaMultipoleForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            (void) fflush( log );
            exit(-1);
        }
    }

    // get supplementary fields

    int isNotEof                = 1;
    int totalFields                = 2;
    int fieldCount                 = 0;
    int done                       = 0;
    while( done == 0 && isNotEof ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
 
            std::string field        = tokens[0];
            if( field == "#" ){
   
               // skip
               if( log ){
                   (void) fprintf( log, "skip <%s>\n", field.c_str());
                   (void) fflush( log );
               }
   
            } else if( field == "AmoebaMultipoleEnd" ){
                done++;
            } else if( field == AMOEBA_MULTIPOLE_ROTATION_MATRICES  || 
                       field == AMOEBA_MULTIPOLE_ROTATED            ||
                       field == AMOEBA_FIXED_E                      ||
                       field == AMOEBA_FIXED_E_GK                   ||
                       field == AMOEBA_INDUCDED_DIPOLES             ||
                       field == AMOEBA_INDUCDED_DIPOLES_GK          ){
                fieldCount++;
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                readVectorOfDoubleVectors( filePtr, tokens, vectorOfDoubleVectors, lineCount, field, log );
                supplementary[field] = vectorOfDoubleVectors;
            } else if( field == "AmoebaMultipoleCovalent" ){
                fieldCount++;
                readAmoebaMultipoleCovalent( filePtr, multipoleForce, lineCount, log );
            }
        }
    }

    // set parameters if available

    MapStringStringI isPresent = inputArgumentMap.find( MUTUAL_INDUCED_MAX_ITERATIONS );
    if( isPresent != inputArgumentMap.end() ){
         multipoleForce->setMutualInducedMaxIterations( atoi( isPresent->second.c_str() ) );
    }

    isPresent = inputArgumentMap.find( MUTUAL_INDUCED_TARGET_EPSILON );
    if( isPresent != inputArgumentMap.end() ){
         multipoleForce->setMutualInducedTargetEpsilon( atof( isPresent->second.c_str() ) );
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){

        double dipoleConversion        = AngstromToNm;
        double quadrupoleConversion    = AngstromToNm*AngstromToNm;
        double polarityConversion      = AngstromToNm*AngstromToNm*AngstromToNm;
        double dampingFactorConversion = sqrt( AngstromToNm );
      
        multipoleForce->setAEwald(                multipoleForce->getAEwald()/AngstromToNm );
        multipoleForce->setCutoffDistance(        multipoleForce->getCutoffDistance()*AngstromToNm );
        multipoleForce->setScalingDistanceCutoff( multipoleForce->getScalingDistanceCutoff()*AngstromToNm );

        Vec3 a,b,c;
        system.getDefaultPeriodicBoxVectors( a, b, c);

        a[0] *= AngstromToNm;
        a[1] *= AngstromToNm;
        a[2] *= AngstromToNm;

        b[0] *= AngstromToNm;
        b[1] *= AngstromToNm;
        b[2] *= AngstromToNm;

        c[0] *= AngstromToNm;
        c[1] *= AngstromToNm;
        c[2] *= AngstromToNm;

        system.setDefaultPeriodicBoxVectors( a, b, c);

        for( int ii = 0; ii < multipoleForce->getNumMultipoles();  ii++ ){

            int axisType, zAxis, xAxis;
            std::vector<double> dipole;
            std::vector<double> quadrupole;
            double charge, thole, dampingFactor, polarity;
            multipoleForce->getMultipoleParameters( ii, charge, dipole, quadrupole, axisType, zAxis, xAxis, thole, dampingFactor, polarity );

            for( unsigned int jj = 0; jj < dipole.size(); jj++ ){
                dipole[jj] *= dipoleConversion;
            }
            for( unsigned int jj = 0; jj < quadrupole.size(); jj++ ){
                quadrupole[jj] *= quadrupoleConversion;
            }
            polarity          *= polarityConversion;
            dampingFactor     *= dampingFactorConversion;

            multipoleForce->setMultipoleParameters( ii, charge, dipole, quadrupole, axisType, zAxis, xAxis, thole, dampingFactor, polarity );

        }
    } else {
        float electricConstant         = static_cast<float>(multipoleForce->getElectricConstant());
              electricConstant        /= static_cast<float>(AngstromToNm*CalToJoule);
        multipoleForce->setElectricConstant( electricConstant );
    }

    // diagnostics

    if( log ){

        (void) fprintf( log, "%s Sample of parameters using %s units.\n", methodName.c_str(),
                        (useOpenMMUnits ? "OpenMM" : "Amoeba") );

        std::string nonbondedMethod = multipoleForce->getNonbondedMethod( ) == AmoebaMultipoleForce::PME ? "PME" : "NoCutoff";
        (void) fprintf( log, "NonbondedMethod=%s aEwald=%15.7e cutoff=%15.7e.\n", nonbondedMethod.c_str(),
                        multipoleForce->getAEwald(), multipoleForce->getCutoffDistance() );

        Vec3 a,b,c;
        system.getDefaultPeriodicBoxVectors( a, b, c );
        (void) fprintf( log, "Box=[%12.3f %12.3f %12.3f] [%12.3f %12.3f %12.3f] [%12.3f %12.3f %12.3f]\n",
                        a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );

        (void) fprintf( log, "Supplementary fields %u: ",  static_cast<unsigned int>(supplementary.size())  );
        for( MapStringVectorOfVectorsCI ii = supplementary.begin(); ii != supplementary.end();  ii++ ){
            (void) fprintf( log, "%s ", (*ii).first.c_str() );
        }
        (void) fprintf( log, "\n" );
        (void) fflush( log );
        AmoebaMultipoleForce::CovalentType covalentTypes[AmoebaMultipoleForce::CovalentEnd] = {
                  AmoebaMultipoleForce::Covalent12,
                  AmoebaMultipoleForce::Covalent13,
                  AmoebaMultipoleForce::Covalent14,
                  AmoebaMultipoleForce::Covalent15,
                  AmoebaMultipoleForce::PolarizationCovalent11,
                  AmoebaMultipoleForce::PolarizationCovalent12,
                  AmoebaMultipoleForce::PolarizationCovalent13,
                  AmoebaMultipoleForce::PolarizationCovalent14 };
     
        //static const unsigned int maxPrint   = MAX_PRINT;
        static const unsigned int maxPrint   = 15;
        unsigned int arraySize               = static_cast<unsigned int>(multipoleForce->getNumMultipoles());
        (void) fprintf( log, "%u maxIter=%d targetEps=%15.7e\n",
                        arraySize,
                        multipoleForce->getMutualInducedMaxIterations(),
                        multipoleForce->getMutualInducedTargetEpsilon() );
        (void) fprintf( log, "Sample of AmoebaMultipoleForce parameters\n" );
        (void) fflush( log );
 
        for( unsigned int ii = 0; ii < arraySize;  ii++ ){
            int axisType, zAxis, xAxis;
            std::vector<double> dipole;
            std::vector<double> quadrupole;
            double charge, thole, dampingFactor, polarity;
            multipoleForce->getMultipoleParameters( ii, charge, dipole, quadrupole, axisType, zAxis, xAxis, thole, dampingFactor, polarity );
            (void) fprintf( log, "%8d %8d %8d %8d q %10.4f thl %10.4f pgm %10.4f pol %10.4f d[%10.4f %10.4f %10.4f]\n",
                            ii, axisType, zAxis, xAxis, charge, thole, dampingFactor, polarity, dipole[0], dipole[1], dipole[2] );
            (void) fprintf( log, "   q[%10.4f %10.4f %10.4f] [%10.4f %10.4f %10.4f] [%10.4f %10.4f %10.4f]\n",
                            quadrupole[0], quadrupole[1], quadrupole[2],
                            quadrupole[3], quadrupole[4], quadrupole[5],
                            quadrupole[6], quadrupole[7], quadrupole[8] );
   
            for( int jj = 0; jj < AmoebaMultipoleForce::CovalentEnd; jj++ ){
                std::vector<int> covalentAtoms;
                multipoleForce->getCovalentMap( ii, covalentTypes[jj], covalentAtoms );
                (void) fprintf( log, "   CovTypeId=%d %u [", jj, static_cast<unsigned int>(covalentAtoms.size()) );
                for( unsigned int kk = 0; kk < covalentAtoms.size(); kk++ ){
                    (void) fprintf( log, "%5d ", covalentAtoms[kk] );
                }
                (void) fprintf( log, "]\n" );
                (void) fflush( log );
            }

            // skip to end
   
            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
        (void) fflush( log );
    }

    return multipoleForce->getNumMultipoles();
}

/**---------------------------------------------------------------------------------------

    Read GK Force parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of parameters read

    --------------------------------------------------------------------------------------- */

static int readAmoebaGeneralizedKirkwoodParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                                    System& system, int useOpenMMUnits,
                                                    MapStringVectorOfVectors& supplementary,
                                                    MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaGeneralizedKirkwoodParameters";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no GK terms entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    AmoebaGeneralizedKirkwoodForce* gbsaObcForce   = new AmoebaGeneralizedKirkwoodForce();
    MapStringIntI forceActive                      = forceMap.find( AMOEBA_GK_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
       system.addForce( gbsaObcForce );
       if( log ){
          (void) fprintf( log, "GK force is being included.\n" );
       }
    } else if( log ){
       (void) fprintf( log, "GK force is not being included.\n" );
    }

    int numberOfParticles           = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of GK force terms=%d\n", methodName.c_str(), numberOfParticles );
       (void) fflush( log );
    }
    for( int ii = 0; ii < numberOfParticles; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       int tokenIndex = 0;
       if( lineTokens.size() > 3 ){
          int index            = atoi( lineTokens[tokenIndex++].c_str() );
          double charge        = atof( lineTokens[tokenIndex++].c_str() );
          double radius        = atof( lineTokens[tokenIndex++].c_str() );
          double scalingFactor = atof( lineTokens[tokenIndex++].c_str() );
          gbsaObcForce->addParticle( charge, radius, scalingFactor );
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s GK force tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 2 ){
       StringVector tokens;
       isNotEof = readLine( filePtr, tokens, lineCount, log );
       if( isNotEof && tokens.size() > 0 ){

          std::string field       = tokens[0];
          if( field == "SoluteDielectric" ){
             gbsaObcForce->setSoluteDielectric( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "SolventDielectric" ){
             gbsaObcForce->setSolventDielectric( atof( tokens[1].c_str() ) );
             hits++;
          } else {
                char buffer[1024];
                (void) sprintf( buffer, "%s read past GK block at line=%d\n", methodName.c_str(), *lineCount );
                throwException(__FILE__, __LINE__, buffer );
                exit(-1);
          }
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s invalid token count at line=%d?\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // get supplementart fields

    isNotEof                       = 1;
    int totalFields                = 2;
    int fieldCount                 = 0;
    int done                       = 0;
    while( done == 0 && isNotEof ){
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, lineCount, log );
        if( isNotEof && tokens.size() > 0 ){
 
            std::string field        = tokens[0];
            if( field == "#" ){
   
               // skip
               if( log ){
                   (void) fprintf( log, "skip <%s>\n", field.c_str());
               }
   
            } else if( field == "AmoebaGeneralizedKirkwoodEnd" ){
                done++;
            } else if( field == "AmoebaGeneralizedKirkwoodBornRadii" ){
                fieldCount++;
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                readVectorOfDoubleVectors( filePtr, tokens, vectorOfDoubleVectors, lineCount, field, log );
                supplementary[field] = vectorOfDoubleVectors;
                done++;
            }
        }
    }

    // check if cavity term is to be included

    MapStringStringI isPresent = inputArgumentMap.find( INCLUDE_OBC_CAVITY_TERM );
    if( isPresent != inputArgumentMap.end() ){
          gbsaObcForce->setIncludeCavityTerm( atoi( isPresent->second.c_str() ) );
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){
        gbsaObcForce->setDielectricOffset( 0.009 );
        for( int ii = 0; ii < gbsaObcForce->getNumParticles(); ii++ ){
            double charge, radius, scalingFactor;
            gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor );
            radius      *= AngstromToNm;
            gbsaObcForce->setParticleParameters( ii, charge, radius, scalingFactor );
        }
    } else {
        gbsaObcForce->setDielectricOffset( 0.09 );
        gbsaObcForce->setProbeRadius( 1.4 );
        double surfaceAreaFactor = gbsaObcForce->getSurfaceAreaFactor( );
        surfaceAreaFactor       *= (AngstromToNm*AngstromToNm)/CalToJoule;
        gbsaObcForce->setSurfaceAreaFactor( surfaceAreaFactor );
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint   = MAX_PRINT;
       unsigned int arraySize               = static_cast<unsigned int>(gbsaObcForce->getNumParticles());
       (void) fprintf( log, "%s: sample of GK force parameters; no. of particles=%d using %s units.\n",
                       methodName.c_str(), gbsaObcForce->getNumParticles(), 
                       (useOpenMMUnits ? "OpenMM" : "Amoeba") );
       (void) fprintf( log, "solute/solvent dielectrics: [%10.4f %10.4f] includeCavityTerm=%1d probeRadius=%15.7e SA prefactor=%15.7e\n",
                       gbsaObcForce->getSoluteDielectric(),  gbsaObcForce->getSolventDielectric(),
                       gbsaObcForce->getIncludeCavityTerm(), gbsaObcForce->getProbeRadius( ), gbsaObcForce->getSurfaceAreaFactor( ) );

       for( unsigned int ii = 0; ii < arraySize; ii++ ){
          double charge, radius, scalingFactor;
          gbsaObcForce->getParticleParameters( ii, charge, radius, scalingFactor );
          (void) fprintf( log, "%8d  %15.7e %15.7e %15.7e\n", ii, charge, radius, scalingFactor );
          if( ii == maxPrint ){
             ii = arraySize - maxPrint;
             if( ii < maxPrint )ii = maxPrint;
          }
       }
    }

    return gbsaObcForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba vdw parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of multipole parameters

    --------------------------------------------------------------------------------------- */

static int readAmoebaVdwParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                    System& system,  int useOpenMMUnits,
                                    MapStringVectorOfVectors& supplementary,
                                    MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

     static const std::string methodName      = "readAmoebaVdwParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens

    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no vdw entries???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }

    // create force

    AmoebaVdwForce* vdwForce   = new AmoebaVdwForce();
    MapStringIntI forceActive  = forceMap.find( AMOEBA_VDW_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( vdwForce );
        if( log ){
            (void) fprintf( log, "Amoeba Vdw force is being included.\n" );
        }    
    } else if( log ){
        (void) fprintf( log, "Amoeba Vdw force is not being included.\n" );
    }

    int numberOfParticles = atoi( tokens[1].c_str() );

    // read in parameters
 
    if( log ){
        (void) fprintf( log, "%s number of vdwForce terms=%d\n", methodName.c_str(), numberOfParticles );
    }
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 2 ){
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            //int indexI           = atoi( lineTokens[tokenIndex++].c_str() );
            int indexIV          = atoi( lineTokens[tokenIndex++].c_str() );
            int indexClass       = atoi( lineTokens[tokenIndex++].c_str() );
            double sigma         = atof( lineTokens[tokenIndex++].c_str() );
            double epsilon       = atof( lineTokens[tokenIndex++].c_str() );
            double reduction     = atof( lineTokens[tokenIndex++].c_str() );
            vdwForce->addParticle( indexIV, indexClass, sigma, epsilon, reduction );
        } else {
            (void) fprintf( log, "%s AmoebaVdwForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            exit(-1);
        }
    }

    // scale factors -- just check that have not changed from assummed values
    // abort if have changed

    StringVector lineTokensT;
    int isNotEof = readLine( filePtr, lineTokensT, lineCount, log );
    if( lineTokensT[0] == "AmoebaVdw14_7Scales" ){
        int tokenIndex    = 1;
        double scale2     = atof( lineTokensT[tokenIndex++].c_str() );
        double scale3     = atof( lineTokensT[tokenIndex++].c_str() );
        double scale4     = atof( lineTokensT[tokenIndex++].c_str() );
        double scale5     = atof( lineTokensT[tokenIndex++].c_str() );
        if( fabs( scale2 )       > 0.0 || 
            fabs( scale3 )       > 0.0 || 
            fabs( 1.0 - scale4 ) > 0.0 || 
            fabs( 1.0 - scale5 ) > 0.0 ){
            char buffer[1024];
            (void) sprintf( buffer, "Vdw scaling factors different from assummed values [0.0 0.0 1.0 1.0] [%12.5e %12.5e %12.5e %12.5e]\n",
                            scale2, scale3, scale4, scale5 );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
        }
    }

    lineTokensT.resize(0);
    isNotEof = readLine( filePtr, lineTokensT, lineCount, log );
    if( lineTokensT[0] == "AmoebaVdw14_7Exclusion" ){
        int numberOfParticles =  atoi( lineTokensT[1].c_str() );
        for( int ii = 0; ii < numberOfParticles; ii++ ){
            StringVector lineTokens;
            int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
            std::vector< int > exclusions;
            if( lineTokens.size() > 1 ){
                int tokenIndex       = 0;
                int index            = atoi( lineTokens[tokenIndex++].c_str() );
                int exclusionCount   = atoi( lineTokens[tokenIndex++].c_str() );
                for( int jj = 0; jj < exclusionCount; jj++ ){
                    int atomIndex = atoi( lineTokens[tokenIndex++].c_str() );
                    exclusions.push_back( atomIndex );
                }
                vdwForce->setParticleExclusions( ii, exclusions );
            } else {
                (void) fprintf( log, "%s AmoebaVdwForce exclusion tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
                exit(-1);
            }
        }
    }

    // 14-7 factors -- just check that have not changed from assummed values
    // abort if have changed

    lineTokensT.resize(0);
    isNotEof = readLine( filePtr, lineTokensT, lineCount, log );
    if( lineTokensT[0] == "AmoebaVdw14_7Hal" ){
        int tokenIndex    = 1;
        double hal1       = atof( lineTokensT[tokenIndex++].c_str() );
        double hal2       = atof( lineTokensT[tokenIndex++].c_str() );
        if( fabs( hal1 - 0.07 )  > 0.0 || 
            fabs( hal2 - 0.12 )  > 0.0 ){ 
            char buffer[1024];
            (void) sprintf( buffer, "Vdw hal values different from assummed values [0.07 0.12 ] [%12.5e %12.5e]\n",
                            hal1, hal2 );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
        }
    }

    // get combining rule

    lineTokensT.resize(0);
    isNotEof = readLine( filePtr, lineTokensT, lineCount, log );
    if( lineTokensT[0] == "AmoebaVdw14_CombiningRule" ){
        int tokenIndex    = 1;
        std::string sigmaCombiningRule   = lineTokensT[tokenIndex++].c_str();
        std::string epsilonCombiningRule = lineTokensT[tokenIndex++].c_str();
        vdwForce->setSigmaCombiningRule( sigmaCombiningRule );
        vdwForce->setEpsilonCombiningRule( epsilonCombiningRule );
    }

    // convert units to kJ-nm from kCal-Angstrom?

    if( useOpenMMUnits ){
        for( int ii = 0; ii < vdwForce->getNumParticles();  ii++ ){
            int indexIV, indexClass;
            double sigma, epsilon, reduction;
            vdwForce->getParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
            sigma        *= AngstromToNm;
            epsilon      *= CalToJoule;
            vdwForce->setParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
        }
    }

    // diagnostics
 
    if( log ){

        //static const unsigned int maxPrint   = MAX_PRINT;
        static const int maxPrint            = 15;
        unsigned int arraySize               = static_cast<unsigned int>(vdwForce->getNumParticles());
        (void) fprintf( log, "%s: %u sample of AmoebaVdwForce parameters using %s units; combining rules=[sig=%s eps=%s]\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba"),
                        vdwForce->getSigmaCombiningRule().c_str(), vdwForce->getEpsilonCombiningRule().c_str() );
  
        for( int ii = 0; ii < vdwForce->getNumParticles();  ii++ ){
            int indexIV, indexClass;
            double sigma, epsilon, reduction;
            std::vector< int > exclusions;
            vdwForce->getParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
            vdwForce->getParticleExclusions( ii, exclusions );
            (void) fprintf( log, "%8d %8d %8d sig=%10.4f eps=%10.4f redct=%10.4f ",
                            ii, indexIV, indexClass, sigma, epsilon, reduction );

            (void) fprintf( log, "Excl=%3u [", static_cast<unsigned int>(exclusions.size()) );
            for( unsigned int jj = 0; jj < exclusions.size();  jj++ ){
                (void) fprintf( log, "%5d ", exclusions[jj] );
            }
            (void) fprintf( log, "]\n", exclusions.size() );

            // skip to end
   
            if( ii == maxPrint && static_cast<int>(arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
        (void) fflush( log );
    }

    return vdwForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

    Read Amoeba WCA dispersion parameters

    @param filePtr              file pointer to parameter file
    @param forceMap             map of forces to be included
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param useOpenMMUnits       if set, use OpenMM units (override input (kcal/A) units)
    @param inputArgumentMap     supplementary arguments
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of particles

    --------------------------------------------------------------------------------------- */

static int readAmoebaWcaDispersionParameters( FILE* filePtr, MapStringInt& forceMap, const StringVector& tokens,
                                              System& system, int useOpenMMUnits, MapStringVectorOfVectors& supplementary,
                                              MapStringString& inputArgumentMap, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

     static const std::string methodName      = "readAmoebaWcaDispersionParameters";
    
// ---------------------------------------------------------------------------------------

    // validate number of tokens
 
    if( tokens.size() < 1 ){
        char buffer[1024];
        (void) sprintf( buffer, "%s no wca entries???\n", methodName.c_str() );
        throwException(__FILE__, __LINE__, buffer );
        exit(-1);
    }
 
   // create force

    AmoebaWcaDispersionForce* wcaDispersionForce   = new AmoebaWcaDispersionForce();
    MapStringIntI forceActive  = forceMap.find( AMOEBA_WCA_DISPERSION_FORCE );
    if( forceActive != forceMap.end() && (*forceActive).second ){
        system.addForce( wcaDispersionForce );
        if( log ){
            (void) fprintf( log, "Amoeba WcaDispersion force is being included.\n" );
        }    
    } else if( log ){
        (void) fprintf( log, "Amoeba WcaDispersion force is not being included.\n" );
    }

    int numberOfParticles = atoi( tokens[1].c_str() );

    // read in parameters
 
    if( log ){
        (void) fprintf( log, "%s number of wcaDispersionForce terms=%d\n", methodName.c_str(), numberOfParticles );
    }

    std::vector<double> maxDispersionEnergyVector;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        StringVector lineTokens;
        int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
        if( lineTokens.size() > 2 ){
            int tokenIndex       = 0;
            int index            = atoi( lineTokens[tokenIndex++].c_str() );
            double radius        = atof( lineTokens[tokenIndex++].c_str() );
            double epsilon       = atof( lineTokens[tokenIndex++].c_str() );
            wcaDispersionForce->addParticle( radius, epsilon );

            if( tokenIndex < static_cast<int>(lineTokens.size()) ){
                double cdisp         = atof( lineTokens[tokenIndex++].c_str() );
                maxDispersionEnergyVector.push_back( cdisp );
            }

        } else {
            (void) fprintf( log, "%s AmoebaWcaDispersionForce tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
            exit(-1);
        }
    }

    int isNotEof                   = 1;
    int hits                       = 0;
    while( hits < 6 ){
       StringVector tokens;
       isNotEof = readLine( filePtr, tokens, lineCount, log );
       if( isNotEof && tokens.size() > 0 ){
          std::string field       = tokens[0];
          if( field == "AmoebaWcaDispersionAwater" ){
             wcaDispersionForce->setAwater( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaWcaDispersionSlevy" ){
             wcaDispersionForce->setSlevy( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaWcaDispersionShctd" ){
             wcaDispersionForce->setShctd( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaWcaDispersionDispoff" ){
             wcaDispersionForce->setDispoff( atof( tokens[1].c_str() ) );
             hits++;
          } else if( field == "AmoebaWcaDispersionEps" ){
             wcaDispersionForce->setEpso( atof( tokens[1].c_str() ) );
             wcaDispersionForce->setEpsh( atof( tokens[2].c_str() ) );
             hits++;
          } else if( field == "AmoebaWcaDispersionRmin" ){
             wcaDispersionForce->setRmino( atof( tokens[1].c_str() ) );
             wcaDispersionForce->setRminh( atof( tokens[2].c_str() ) );
             hits++;
          } else {
                char buffer[1024];
                (void) sprintf( buffer, "%s read past WcaDispersion block at line=%d\n", methodName.c_str(), *lineCount );
                throwException(__FILE__, __LINE__, buffer );
                exit(-1);
          }
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s invalid token count at line=%d?\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // convert to OpenMM units

    if( useOpenMMUnits ){

        // slevy enthalpy-to-free energy scale factor for dispersion

        // awater is number density at STP

        double aWater       = wcaDispersionForce->getAwater( );
        aWater             /= (AngstromToNm*AngstromToNm*AngstromToNm);
        wcaDispersionForce->setAwater( aWater );

        double dispoff      = wcaDispersionForce->getDispoff( );
        dispoff            *= AngstromToNm;
        wcaDispersionForce->setDispoff( dispoff );

        // rmino water-oxygen Rmin for implicit dispersion term
        // rminh water-hydrogen Rmin for implicit dispersion term

        double rmino        = wcaDispersionForce->getRmino( );
        double rminh        = wcaDispersionForce->getRminh( );
        rmino              *= AngstromToNm;
        rminh              *= AngstromToNm;
        wcaDispersionForce->setRmino( rmino );
        wcaDispersionForce->setRminh( rminh );

        // epso water-oxygen epsilon for implicit dispersion term
        // epsh water-hydrogen epsilon for implicit dispersion term

        double epso         = wcaDispersionForce->getEpso( );
        double epsh         = wcaDispersionForce->getEpsh( );
        epso               *= CalToJoule;
        epsh               *= CalToJoule;
        wcaDispersionForce->setEpso( epso );
        wcaDispersionForce->setEpsh( epsh );

        for( int ii = 0; ii < wcaDispersionForce->getNumParticles(); ii++ ){

            double radius, epsilon, maxDispersionEnergy;
            wcaDispersionForce->getParticleParameters( ii, radius, epsilon );
            radius           *= AngstromToNm;
            epsilon          *= CalToJoule;
            wcaDispersionForce->setParticleParameters( ii, radius, epsilon );

            if( ii < static_cast<int>(maxDispersionEnergyVector.size()) ){
                AmoebaWcaDispersionForceImpl::getMaximumDispersionEnergy( *wcaDispersionForce, ii, maxDispersionEnergy );
                double tinkerValue  = maxDispersionEnergyVector[ii];
                       tinkerValue *= CalToJoule;
                double delta        = fabs( maxDispersionEnergy - tinkerValue );
                const char* error   = (delta > 1.0e-05) ? "XXX" : "";
                if( delta > 1.0e-05 && log ){
                    (void) fprintf( log, "useOpenMMUnits: maxDispEDiff=%12.5e %14.7f %14.7f  %s\n",
                                    delta,  maxDispersionEnergy, tinkerValue, error );
                }
            }
        }
    }
 
    // diagnostics

    if( log ){

        //static const unsigned int maxPrint   = MAX_PRINT;
        static const unsigned int maxPrint   = 15;
        unsigned int arraySize               = static_cast<unsigned int>(wcaDispersionForce->getNumParticles());
        (void) fprintf( log, "%s: %u sample of AmoebaVdwForce parameters in %s units.\n",
                        methodName.c_str(), arraySize, (useOpenMMUnits ? "OpenMM" : "Amoeba") );

        (void) fprintf( log, "Eps[%14.7f %14.7f] Rmin[%14.7f %14.7f]\nAwater %14.7f Shctd %14.7f Dispoff %14.7f Slevy %14.7f\n",
                        wcaDispersionForce->getEpso( ), wcaDispersionForce->getEpsh( ),
                        wcaDispersionForce->getRmino( ), wcaDispersionForce->getRminh( ),
                        wcaDispersionForce->getAwater( ), wcaDispersionForce->getShctd( ), wcaDispersionForce->getDispoff( ), wcaDispersionForce->getSlevy( ) );

        for( unsigned int ii = 0; ii < arraySize;  ii++ ){
            double radius, epsilon, maxDispersionEnergy;
            wcaDispersionForce->getParticleParameters( ii, radius, epsilon );
            (void) fprintf( log, "%8d %10.4f %10.4f", ii, radius, epsilon );
            if( ii < maxDispersionEnergyVector.size() ){
                AmoebaWcaDispersionForceImpl::getMaximumDispersionEnergy( *wcaDispersionForce, ii, maxDispersionEnergy );
                if( useOpenMMUnits )maxDispersionEnergy /= CalToJoule;
                double delta = fabs( maxDispersionEnergy - maxDispersionEnergyVector[ii] );
                const char* error  = (delta > 1.0e-05) ? "XXX" : "";
                (void) fprintf( log, " maxDispEDiff=%12.5e %14.7f %14.7f  %s",
                                delta,  maxDispersionEnergy, maxDispersionEnergyVector[ii], error );
            }
            (void) fprintf( log, "\n" );

            // skip to end

            if( ii == maxPrint && (arraySize - maxPrint) > ii ){
                ii = arraySize - maxPrint - 1;
            } 
        }
        (void) fflush( log );

        // check max dispersion energy for all particles

        int errors = 0;
        for( unsigned int ii = 0; ii < arraySize && ii < maxDispersionEnergyVector.size(); ii++ ){
            double maxDispersionEnergy;
            AmoebaWcaDispersionForceImpl::getMaximumDispersionEnergy( *wcaDispersionForce, ii, maxDispersionEnergy );
            if( useOpenMMUnits )maxDispersionEnergy /= CalToJoule;
            double delta = fabs( maxDispersionEnergy - maxDispersionEnergyVector[ii] );
            if( delta > 1.0e-05 ){
                (void) fprintf( log, " maxDispEDiff=%12.5e %14.7f %14.7f  XXX\n", delta, maxDispersionEnergy, maxDispersionEnergyVector[ii] );
                errors++;
            } 
        }
        if( errors ){
            exit(-1);
        } else {
            (void) fprintf( log, "No errors detected in maxDispEnergy!\n" );
        }
        (void) fflush( log );
    }

    return wcaDispersionForce->getNumParticles();
}

/**---------------------------------------------------------------------------------------

    Read Constraints

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of parameters read

    --------------------------------------------------------------------------------------- */

static int readConstraints( FILE* filePtr, const StringVector& tokens, System& system, int useOpenMMUnits,
                            MapStringVectorOfVectors& supplementary, MapStringString& inputArgumentMap,
                            int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readConstraints";
    int applyConstraints                     = 1;
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no constraints terms entry???\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    setIntFromMap( inputArgumentMap, "applyConstraints",  applyConstraints);
    if( log ){
       (void) fprintf( log, "%s: constraints are %sbeing applied.\n", methodName.c_str(), (applyConstraints ? "" : "not ") );
    }

    int numberOfConstraints = atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of constraints=%d\n", methodName.c_str(), numberOfConstraints );
       (void) fflush( log );
    }
    for( int ii = 0; ii < numberOfConstraints; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 3 ){
          int index            = atoi( lineTokens[0].c_str() );
          int particle1        = atoi( lineTokens[1].c_str() );
          int particle2        = atoi( lineTokens[2].c_str() );
          double distance      = atof( lineTokens[3].c_str() );
          if( applyConstraints ){
              system.addConstraint( particle1, particle2, distance );
          }
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s constraint tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // convert to OpenMM units
    // scale constraint distances

    if( useOpenMMUnits ){
       for( int ii = 0; ii < system.getNumConstraints(); ii++ ){
          int particle1, particle2;
          double distance;
          system.getConstraintParameters( ii, particle1, particle2, distance ); 
          distance *= AngstromToNm;
          system.setConstraintParameters( ii, particle1, particle2, distance ); 
       }
    }

    // diagnostics

    if( log && system.getNumConstraints() ){
       int maxPrint = 10;
       (void) fprintf( log, "%s: sample of %d constraints using %s units.\n", methodName.c_str(),
                       system.getNumConstraints(), (useOpenMMUnits ? "OpenMM" : "Amoeba") );
       for( int ii = 0; ii < system.getNumConstraints(); ii++ ){
          int particle1, particle2;
          double distance;
          system.getConstraintParameters( ii, particle1, particle2, distance ); 
          (void) fprintf( log, "%8u %8d %8d %15.7e\n", ii, particle1, particle2, distance );
          if( ii == maxPrint && (system.getNumConstraints() - maxPrint) > ii ){
             ii = system.getNumConstraints() - maxPrint - 1;
          } 
       }
    }

    return system.getNumConstraints();
}

/**---------------------------------------------------------------------------------------

    Read integrator

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param system               System reference
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return integrator

    --------------------------------------------------------------------------------------- */

static Integrator* readIntegrator( FILE* filePtr, const StringVector& tokens, System& system, int* lineCount, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readIntegrator";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s integrator name missing?\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    std::string integratorName = tokens[1];
    if( log ){
       (void) fprintf( log, "%s integrator=%s\n", methodName.c_str(), integratorName.c_str() );
       (void) fflush( log );
    }

    // set number of parameters (lines to read)

    int readLines;
    if( integratorName == "LangevinIntegrator" ){
       readLines = 5;
    } else if( integratorName == "VariableLangevinIntegrator" ){
       readLines = 6;
    } else if( integratorName == "VerletIntegrator" ){
       readLines = 2;
    } else if( integratorName == "VariableVerletIntegrator" ){
       readLines = 3;
    } else if( integratorName == "BrownianIntegrator" ){
       readLines = 5;
    } else {
       (void) fprintf( log, "%s integrator=%s not recognized.\n", methodName.c_str(), integratorName.c_str() );
       (void) fflush( log );
       exit(-1);
    }
    
    // read in parameters

    double stepSize               = 0.001;
    double constraintTolerance    = 1.0e-05;
    double temperature            = 300.0;
    double friction               = 0.01099;
    double errorTolerance         = 1.0e-05;
    int randomNumberSeed          = 1993;

    for( int ii = 0; ii < readLines; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 1 ){
          if( lineTokens[0] == "StepSize" ){
             stepSize            =  atof( lineTokens[1].c_str() );
          } else if( lineTokens[0] == "ConstraintTolerance" ){
             constraintTolerance =  atof( lineTokens[1].c_str() );
          } else if( lineTokens[0] == "Temperature" ){
             temperature         =  atof( lineTokens[1].c_str() );
          } else if( lineTokens[0] == "Friction" ){
             friction            =  atof( lineTokens[1].c_str() );
          } else if( lineTokens[0] == "ErrorTolerance" ){
             errorTolerance      =  atof( lineTokens[1].c_str() );
          } else if( lineTokens[0] == "RandomNumberSeed" ){
             randomNumberSeed    =  atoi( lineTokens[1].c_str() );
          } else {
             (void) fprintf( log, "%s integrator field=%s not recognized.\n", methodName.c_str(), lineTokens[0].c_str() );
             (void) fflush( log );
             exit(-1);
          }
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s integrator parameters incomplete at line=%d\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // build integrator

    Integrator* returnIntegrator = NULL;

    if( integratorName == "LangevinIntegrator" ){

        LangevinIntegrator* langevinIntegrator = new LangevinIntegrator( temperature, friction, stepSize );
        langevinIntegrator->setRandomNumberSeed( randomNumberSeed );
        returnIntegrator = langevinIntegrator;

    } else if( integratorName == "VariableLangevinIntegrator" ){

        VariableLangevinIntegrator* variableLangevinIntegrator = new VariableLangevinIntegrator( temperature, friction, errorTolerance );
        variableLangevinIntegrator->setStepSize( stepSize );
        variableLangevinIntegrator->setRandomNumberSeed( randomNumberSeed );
        returnIntegrator = variableLangevinIntegrator;

    } else if( integratorName == "VerletIntegrator" ){
       returnIntegrator = new VerletIntegrator( stepSize );
    } else if( integratorName == "VariableVerletIntegrator" ){
       returnIntegrator = new VariableVerletIntegrator( errorTolerance );
       returnIntegrator->setStepSize( stepSize );
    } else if( integratorName == "BrownianIntegrator" ){
        BrownianIntegrator* brownianIntegrator = new BrownianIntegrator( temperature, friction, stepSize );
        brownianIntegrator->setRandomNumberSeed( randomNumberSeed );
        returnIntegrator = brownianIntegrator;
    }
    returnIntegrator->setConstraintTolerance( constraintTolerance );
    
    if( log ){
       static const unsigned int maxPrint   = MAX_PRINT;
       (void) fprintf( log, "%s: ", methodName.c_str() );
       (void) fprintf( log, "stepSize=%12.3f constraint tolerance=%12.3e ", stepSize, constraintTolerance );
       if( integratorName == "LangevinIntegrator"  || integratorName == "BrownianIntegrator"  ||  integratorName == "VariableLangevinIntegrator" ){
           (void) fprintf( log, "temperature=%12.3f friction=%12.3f seed=%d ", temperature, friction, randomNumberSeed );
       }
       if( integratorName == "VariableLangevinIntegrator"  || integratorName == "VariableVerletIntegrator" ){
           (void) fprintf( log, "error tolerance=%12.3e", errorTolerance);
       }
       (void) fprintf( log, "\n" );
    }

    return returnIntegrator;
}

/**---------------------------------------------------------------------------------------

    Read arrays of Vec3s (coordinates/velocities/forces/...)

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param coordinates          Vec3 array
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of entries read

    --------------------------------------------------------------------------------------- */

static int readVec3( FILE* filePtr, const StringVector& tokens, std::vector<Vec3>& coordinates, int* lineCount,
                      std::string typeName, FILE* log ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readVec3";
    
// ---------------------------------------------------------------------------------------

    if( tokens.size() < 1 ){
       char buffer[1024];
       (void) sprintf( buffer, "%s no entries?\n", methodName.c_str() );
       throwException(__FILE__, __LINE__, buffer );
       exit(-1);
    }

    int numberOfCoordinates= atoi( tokens[1].c_str() );
    if( log ){
       (void) fprintf( log, "%s number of %s=%d\n", methodName.c_str(), typeName.c_str(), numberOfCoordinates );
       (void) fflush( log );
    }
    for( int ii = 0; ii < numberOfCoordinates; ii++ ){
       StringVector lineTokens;
       int isNotEof = readLine( filePtr, lineTokens, lineCount, log );
       if( lineTokens.size() > 3 ){
          int index            = atoi( lineTokens[0].c_str() );
          double xCoord        = atof( lineTokens[1].c_str() );
          double yCoord        = atof( lineTokens[2].c_str() );
          double zCoord        = atof( lineTokens[3].c_str() );
          coordinates.push_back( Vec3( xCoord, yCoord, zCoord ) );
       } else {
          char buffer[1024];
          (void) sprintf( buffer, "%s coordinates tokens incomplete at line=%d\n", methodName.c_str(), *lineCount );
          throwException(__FILE__, __LINE__, buffer );
          exit(-1);
       }
    }

    // diagnostics

    if( log ){
       static const unsigned int maxPrint = MAX_PRINT;
       unsigned int arraySize             = coordinates.size();
       (void) fprintf( log, "%s: sample of vec3 (raw values): %u\n", methodName.c_str(), arraySize );
       for( unsigned int ii = 0; ii < arraySize; ii++ ){
          (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]\n", ii,
                          coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] );
          // skip to end

          if( ii == maxPrint && (arraySize - maxPrint) > ii ){
             ii = arraySize - maxPrint - 1;
          } 
       }
    }

    return static_cast<int>(coordinates.size());
}

/**---------------------------------------------------------------------------------------

    Read arrays of Vec3s (coordinates/velocities/forces/...)

    @param filePtr              file pointer to parameter file
    @param tokens               array of strings from first line of parameter file for this block of parameters
    @param coordinates          Vec3 array
    @param lineCount            used to track line entries read from parameter file
    @param log                  log file pointer -- may be NULL

    @return number of entries read

    --------------------------------------------------------------------------------------- */

static int addForces( std::vector<Vec3>& forceToAdd, std::vector<Vec3>& forceAccumulator ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "addForces";
    
// ---------------------------------------------------------------------------------------

    if( forceAccumulator.size() < forceToAdd.size() ){
       forceAccumulator.resize( forceToAdd.size() );
    }
    for( unsigned int ii = 0; ii < forceToAdd.size(); ii++ ){
       forceAccumulator[ii][0] += forceToAdd[ii][0];
       forceAccumulator[ii][1] += forceToAdd[ii][1];
       forceAccumulator[ii][2] += forceToAdd[ii][2];
    }
    return static_cast<int>(forceAccumulator.size());
}

static void getStringForceMap( System& system, MapStringForce& forceMap, FILE* log ){

    // print active forces and relevant parameters

    for( int ii = 0; ii < system.getNumForces(); ii++ ) {

        int hit                 = 0;
        Force& force            = system.getForce(ii);

        // bond

        if( !hit ){

            try {
               AmoebaHarmonicBondForce& harmonicBondForce = dynamic_cast<AmoebaHarmonicBondForce&>(force);
               forceMap[AMOEBA_HARMONIC_BOND_FORCE]       = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // multipole

        if( !hit ){

            try {
               AmoebaMultipoleForce& multipoleForce = dynamic_cast<AmoebaMultipoleForce&>(force);
               forceMap[AMOEBA_MULTIPOLE_FORCE]     = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // out-of-plane-bend Force

        if( !hit ){

            try {
               AmoebaOutOfPlaneBendForce& outOfPlaneBend = dynamic_cast<AmoebaOutOfPlaneBendForce&>(force);
               forceMap[AMOEBA_OUT_OF_PLANE_BEND_FORCE]  = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // Pi-torsion Force

        if( !hit ){

            try {
               AmoebaPiTorsionForce& piTorsion         = dynamic_cast<AmoebaPiTorsionForce&>(force);
               forceMap[AMOEBA_PI_TORSION_FORCE]       = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // Torsion-torsion Force

        if( !hit ){

            try {
               AmoebaTorsionTorsionForce& torsionTorsion = dynamic_cast<AmoebaTorsionTorsionForce&>(force);
               forceMap[AMOEBA_TORSION_TORSION_FORCE]    = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // Torsion Force

        if( !hit ){

            try {
               AmoebaTorsionForce& torsion          = dynamic_cast<AmoebaTorsionForce&>(force);
               forceMap[AMOEBA_TORSION_FORCE]       = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // stretch bend force

        if( !hit ){

            try {
               AmoebaStretchBendForce& stretchBend = dynamic_cast<AmoebaStretchBendForce&>(force);
               forceMap[AMOEBA_STRETCH_BEND_FORCE] = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // vdw force

        if( !hit ){

            try {
               AmoebaVdwForce& vdw              = dynamic_cast<AmoebaVdwForce&>(force);
               forceMap[AMOEBA_VDW_FORCE]       = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // WCA dspersion force

        if( !hit ){

            try {
               AmoebaWcaDispersionForce& wcaDispersionForce = dynamic_cast<AmoebaWcaDispersionForce&>(force);
               forceMap[AMOEBA_WCA_DISPERSION_FORCE]        = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // angle

        if( !hit ){
    
            try {
               AmoebaHarmonicAngleForce & harmonicAngleForce = dynamic_cast<AmoebaHarmonicAngleForce&>(force);
               forceMap[AMOEBA_HARMONIC_ANGLE_FORCE]         = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // in-plane angle

        if( !hit ){
    
            try {
               AmoebaHarmonicInPlaneAngleForce & harmonicAngleForce = dynamic_cast<AmoebaHarmonicInPlaneAngleForce&>(force);
               forceMap[AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE]       = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // Kirkwood
    
        if( !hit ){
            try {
               AmoebaGeneralizedKirkwoodForce& kirkwoodForce = dynamic_cast<AmoebaGeneralizedKirkwoodForce&>(force);
               forceMap[AMOEBA_GK_FORCE]                     = &force;
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // COM

        if( !hit ){
    
            try {
               CMMotionRemover& cMMotionRemover = dynamic_cast<CMMotionRemover&>(force);
               hit++;
            } catch( std::bad_cast ){
            }
        }

        if( !hit && log ){
           (void) fprintf( log, "   entry=%2d force not recognized.\n", ii );
        }

    }

    return;
}

static void getForceStrings( System& system, StringVector& forceStringArray, FILE* log ){

    MapStringForce forceMap;
    getStringForceMap( system, forceMap, log );
    for( MapStringForceI ii = forceMap.begin(); ii != forceMap.end(); ii++ ) {
        forceStringArray.push_back( ii->first );
    }
}

/**---------------------------------------------------------------------------------------

    Read parameter file

    @param inputParameterFile   input parameter file name
    @param system               system to which forces based on parameters are to be added
    @param coordinates          Vec3 array containing coordinates on output
    @param velocities           Vec3 array containing velocities on output
    @param inputLog             log file pointer -- may be NULL

    @return number of lines read

    --------------------------------------------------------------------------------------- */

Integrator* readAmoebaParameterFile( const std::string& inputParameterFile, MapStringInt& forceMap, System& system,
                                      std::vector<Vec3>& coordinates, 
                                      std::vector<Vec3>& velocities,
                                      MapStringVec3& forces, MapStringDouble& potentialEnergy,
                                      MapStringVectorOfVectors& supplementary, int useOpenMMUnits,
                                      MapStringString& inputArgumentMap, FILE* inputLog ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "readAmoebaParameterFile";
    int PrintOn                              = 1; 
     
// ---------------------------------------------------------------------------------------
 
    FILE* log;
    if( PrintOn == 0 && inputLog ){
        log = NULL;
    } else {
        log = inputLog;
    }
 
    if( log ){
        (void) fprintf( log, "%s\n", methodName.c_str() );
        (void) fflush( log );
    }   
 
    // open parameter file
 
    FILE* filePtr = openFile( inputParameterFile, "r", log );
    if( filePtr == NULL ){
        char buffer[1024];
        (void) sprintf( buffer, "Input parameter file=<%s> could not be opened -- aborting.\n", inputParameterFile.c_str() );
        throwException(__FILE__, __LINE__, buffer );
    } else if( log ){
        (void) fprintf( log, "Input parameter file=<%s> opened.\n", inputParameterFile.c_str() );
        (void) fflush( log );
    }

    int lineCount                  = 0;
    int version                    = 0; 
    int isNotEof                   = 1;
    Integrator* returnIntegrator   = NULL;

    // loop over lines in file

    while( isNotEof ){
 
        // read line and continue if not EOF and tokens found on line
 
        StringVector tokens;
        isNotEof = readLine( filePtr, tokens, &lineCount, log );
 
        if( isNotEof && tokens.size() > 0 ){
 
            std::string field       = tokens[0];
  
            if( log ){
                (void) fprintf( log, "Field=<%s> at line=%d\n", field.c_str(), lineCount );
                (void) fflush( log );
            }
   
            if( field == "Version" ){
                if( tokens.size() > 1 ){
                   version = atoi( tokens[1].c_str() );
                   if( log ){
                       (void) fprintf( log, "Version=%d at line=%d\n", version, lineCount );
                   }
                }
            } else if( field == "Masses" ){
                readMasses( filePtr, tokens, system, &lineCount, log );
            } else if( field == "CMMotionRemover" ){
                int frequency = atoi( tokens[1].c_str() );
                system.addForce( new CMMotionRemover( frequency ) );
                if( log ){
                    (void) fprintf( log, "CMMotionRemover added w/ frequency=%d at line=%d\n", frequency, lineCount );
                }
   
            // not used any longer -- was used in SASA force

            } else if( field == "AmoebaSurfaceProbe" ){
           
            // All forces/energy
  
            } else if( field == ALL_FORCES ){
                readVec3( filePtr, tokens, forces[ALL_FORCES], &lineCount, field, log );
            } else if( field == "AllEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[ALL_FORCES] = atof( tokens[1].c_str() ); 
                }
   
            // AmoebaHarmonicBond
  
            } else if( field == "AmoebaHarmonicBondParameters" ){
                readAmoebaHarmonicBondParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap,&lineCount, log );
            } else if( field == "AmoebaHarmonicBondForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_HARMONIC_BOND_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaHarmonicBondEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_HARMONIC_BOND_FORCE] = atof( tokens[1].c_str() ); 
                }
   
            // AmoebaHarmonicAngle
  
            } else if( field == "AmoebaHarmonicAngleParameters" ){
                readAmoebaHarmonicAngleParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaHarmonicAngleForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_HARMONIC_ANGLE_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaHarmonicAngleEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_HARMONIC_ANGLE_FORCE] = atof( tokens[1].c_str() );
                }
   
            // AmoebaHarmonicInPlaneAngle
  
            } else if( field == "AmoebaHarmonicInPlaneAngleParameters" ){
                readAmoebaHarmonicInPlaneAngleParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaHarmonicInPlaneAngleForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaHarmonicInPlaneAngleEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE] = atof( tokens[1].c_str() );
                }
   
            // AmoebaTorsion
  
            } else if( field == "AmoebaTorsionParameters" ){
                readAmoebaTorsionParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaTorsionForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_TORSION_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaTorsionEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_TORSION_FORCE] = atof( tokens[1].c_str() ); 
                }
   
            // AmoebaPiTorsion
  
            } else if( field == "AmoebaPiTorsionParameters" ){
               readAmoebaPiTorsionParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaPiTorsionForce" ){
               readVec3( filePtr, tokens, forces[AMOEBA_PI_TORSION_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaPiTorsionEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_PI_TORSION_FORCE] = atof( tokens[1].c_str() );
                }
  
            // AmoebaStretchBend
  
            } else if( field == "AmoebaStretchBendParameters" ){
                readAmoebaStretchBendParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaStretchBendForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_STRETCH_BEND_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaStretchBendEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_STRETCH_BEND_FORCE] = atof( tokens[1].c_str() );
                }
  
            // AmoebaOutOfPlaneBend
  
            } else if( field == "AmoebaOutOfPlaneBendParameters" ){
                readAmoebaOutOfPlaneBendParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaOutOfPlaneBendForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_OUT_OF_PLANE_BEND_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaOutOfPlaneBendEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_OUT_OF_PLANE_BEND_FORCE] = atof( tokens[1].c_str() ); 
                }
  
            // AmoebaTorsionTorsion
  
            } else if( field == "AmoebaTorsionTorsionParameters" ){
                readAmoebaTorsionTorsionParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaTorsionTorsionForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_TORSION_TORSION_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaTorsionTorsionEnergy" ){
                if( tokens.size() > 1 ){
                   potentialEnergy[AMOEBA_TORSION_TORSION_FORCE] = atof( tokens[1].c_str() ); 
                }
  
            // AmoebaMultipole
  
            } else if( field == "AmoebaMultipoleParameters" ){
                readAmoebaMultipoleParameters( filePtr, version, forceMap, tokens, system, useOpenMMUnits, supplementary, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaMultipoleForce" || field == "AmoebaPmeForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_MULTIPOLE_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaMultipoleEnergy" || field == "AmoebaPmeEnergy" ){
               if( tokens.size() > 1 ){
                 potentialEnergy[AMOEBA_MULTIPOLE_FORCE] = atof( tokens[1].c_str() ); 
               }
            } else if( field == "AmoebaRealPmeForce"                   || 
                       field == "AmoebaKSpacePmeForce"                 ||
                       field == "AmoebaSelfPmeForce"                   ){
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                readVectorOfDoubleVectors( filePtr, tokens, vectorOfDoubleVectors, &lineCount, field, log );
                supplementary[field] = vectorOfDoubleVectors;
            } else if( field == "AmoebaRealPmeEnergy"                  || 
                       field == "AmoebaKSpacePmeEnergy"                ||
                       field == "AmoebaSelfPmeEnergy"                  ){
                double value = atof( tokens[1].c_str() );
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                std::vector<double> doubleVectors;
                doubleVectors.push_back( value );
                vectorOfDoubleVectors.push_back( doubleVectors ); 
                supplementary[field] =  vectorOfDoubleVectors;

            // Amoeba GK

            } else if( field == "AmoebaGeneralizedKirkwoodParameters" ){
                readAmoebaGeneralizedKirkwoodParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, supplementary, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaGkForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_GK_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaGkAndCavityForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_GK_CAVITY_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaGk_A_ForceAndTorque"            || 
                       field == "AmoebaGk_A_Force"                     ||
                       field == "AmoebaSurfaceParameters"              ||
                       field == "AmoebaGk_A_DrB"                       ||
                       field == "AmoebaDBorn"                          ||
                       field == "AmoebaBorn1Force"                     ||
                       field == "AmoebaBornForce"                      ||
                       field == "AmoebaGkEdiffForceAndTorque"          ||
                       field == "AmoebaGkEdiffForce"                   ){
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                readVectorOfDoubleVectors( filePtr, tokens, vectorOfDoubleVectors, &lineCount, field, log );
                supplementary[field] = vectorOfDoubleVectors;
            } else if( field == "AmoebaGkEnergy"            || 
                       field == "AmoebaGkEdiffEnergy"       ||
                       field == "AmoebaGk_A_Energy"         ||
                       field == "AmoebaBorn1Energy"         ||
                       field == "AmoebaBornEnergy"          ||
                       field == "AmoebaGkAndCavityEnergy"      ){
                double value = atof( tokens[1].c_str() );
                std::vector< std::vector<double> > vectorOfDoubleVectors;
                std::vector<double> doubleVectors;
                doubleVectors.push_back( value );
                vectorOfDoubleVectors.push_back( doubleVectors ); 
                supplementary[field] =  vectorOfDoubleVectors;
                if( field == "AmoebaGkEnergy" ){
                    potentialEnergy[AMOEBA_GK_FORCE] = value; 
                } else if( field == "AmoebaGkAndCavityEnergy" ){
                    potentialEnergy[AMOEBA_GK_CAVITY_FORCE] = value; 
                }
 
            // Amoeba Vdw
 
            } else if( field == "AmoebaVdw14_7SigEpsTable"  || field == "AmoebaVdw14_7Reduction" ){
                readAmoebaVdwParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, supplementary, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaVdwForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_VDW_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaVdwEnergy" ){
                if( tokens.size() > 1 ){
                    potentialEnergy[AMOEBA_VDW_FORCE] = atof( tokens[1].c_str() ); 
                }
            } else if( field == "AmoebaWcaDispersionParameters" ){
                readAmoebaWcaDispersionParameters( filePtr, forceMap, tokens, system, useOpenMMUnits, supplementary, inputArgumentMap, &lineCount, log );
            } else if( field == "AmoebaWcaDispersionForce" ){
                readVec3( filePtr, tokens, forces[AMOEBA_WCA_DISPERSION_FORCE], &lineCount, field, log );
            } else if( field == "AmoebaWcaDispersionEnergy" ){
                if( tokens.size() > 1 ){
                    potentialEnergy[AMOEBA_WCA_DISPERSION_FORCE] = atof( tokens[1].c_str() ); 
                }
            } else if( field == "Constraints" ){
                readConstraints( filePtr, tokens, system, useOpenMMUnits, supplementary, inputArgumentMap, &lineCount, log );
            } else if( field == "Integrator" ){
                returnIntegrator = readIntegrator( filePtr, tokens, system, &lineCount, log );
            } else if( field == "Positions"  || field == "Coordinates" ){
                readVec3( filePtr, tokens, coordinates, &lineCount, field, log );
                if( useOpenMMUnits ){
                    for( unsigned int ii = 0; ii < coordinates.size(); ii++ ){
                        coordinates[ii][0] *= AngstromToNm;
                        coordinates[ii][1] *= AngstromToNm;
                        coordinates[ii][2] *= AngstromToNm;
                    }
                }
            } else if( field == "Velocities" ){
                readVec3( filePtr, tokens, velocities, &lineCount, field, log );
                if( useOpenMMUnits ){
                    for( unsigned int ii = 0; ii < velocities.size(); ii++ ){
                        velocities[ii][0] *= AngstromToNm;
                        velocities[ii][1] *= AngstromToNm;
                        velocities[ii][2] *= AngstromToNm;
                    }
                }
            } else {
                char buffer[1024];
                (void) sprintf( buffer, "Field=<%s> not recognized at line=%d.\n", field.c_str(), lineCount );
                throwException(__FILE__, __LINE__, buffer );
                exit(-1);
            }
       }
    }

    // if integrator not set, default to Verlet integrator

    if( returnIntegrator == NULL ){
       returnIntegrator = new VerletIntegrator(0.001);
    }
 
    // sum energies

    double totalPotentialEnergy = 0.0;
    if( log )(void) fprintf( log, "Potential energies\n" );
   
    double allEnergy = 0.0;
    for( MapStringDoubleI ii = potentialEnergy.begin(); ii != potentialEnergy.end(); ii++ ){
        if( ii->first == ALL_FORCES ){
            allEnergy = ii->second;
        } else if( ii->first != AMOEBA_GK_CAVITY_FORCE ){
            totalPotentialEnergy += ii->second;
        }
        if( log )(void) fprintf( log, "%30s %15.7e\n", ii->first.c_str(), ii->second );
    }
    potentialEnergy["SumOfInputEnergies"] = totalPotentialEnergy;
       
    if( log ){
       MapStringDoubleI isPresent = potentialEnergy.find( AMOEBA_GK_CAVITY_FORCE );
       if( isPresent != potentialEnergy.end() ){
           double cavityEnergy = potentialEnergy[AMOEBA_GK_CAVITY_FORCE] - potentialEnergy[AMOEBA_GK_FORCE];
           (void) fprintf( log, "Cavity energy %15.7e\n", cavityEnergy );
       }
       (void) fprintf( log, "Total PE %15.7e  %15.7e\n", totalPotentialEnergy, allEnergy );
       (void) fprintf( log, "Read %d lines from file=<%s>\n", lineCount, inputParameterFile.c_str() );
       (void) fflush( log );
    }

   return returnIntegrator;
}

/**---------------------------------------------------------------------------------------
 * Initialize forceMap
 *
 * @param forceMap        has w/ force name as key and int as value
 * @param initialValue    initial value
 *
 *
    --------------------------------------------------------------------------------------- */

void initializeForceMap( MapStringInt& forceMap, int initialValue ){

     forceMap[AMOEBA_HARMONIC_BOND_FORCE]            = initialValue;
     forceMap[AMOEBA_HARMONIC_ANGLE_FORCE]           = initialValue;
     forceMap[AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE]  = initialValue;
     forceMap[AMOEBA_TORSION_FORCE]                  = initialValue;
     forceMap[AMOEBA_PI_TORSION_FORCE]               = initialValue;
     forceMap[AMOEBA_STRETCH_BEND_FORCE]             = initialValue;
     forceMap[AMOEBA_OUT_OF_PLANE_BEND_FORCE]        = initialValue;
     forceMap[AMOEBA_TORSION_TORSION_FORCE]          = initialValue;

     forceMap[AMOEBA_MULTIPOLE_FORCE]                = initialValue;
     forceMap[AMOEBA_GK_FORCE]                       = initialValue;
     forceMap[AMOEBA_VDW_FORCE]                      = initialValue;
     forceMap[AMOEBA_WCA_DISPERSION_FORCE]           = initialValue;
 
     return;

}

void checkIntermediateMultipoleQuantities( Context* context, MapStringVectorOfVectors& supplementary,
                                           int useOpenMMUnits, FILE* log ) {

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "checkIntermediateMultipoleQuantities";
 
// ---------------------------------------------------------------------------------------

    // get pointer to AmoebaCudaData for this context

    ContextImpl* contextImpl          = *reinterpret_cast<ContextImpl**>(context);
    void* amoebaCudaDataV             = getAmoebaCudaData( *contextImpl );
    if( log == NULL ){
        return;
    }

    (void) fprintf( log, "%s amoebaCudaDataV=%p\n", methodName.c_str(), amoebaCudaDataV );
    AmoebaCudaData* amoebaCudaData = static_cast<AmoebaCudaData*>(amoebaCudaDataV);
    amoebaGpuContext amoebaGpu     = amoebaCudaData->getAmoebaGpu();

    unsigned int totalMisses       = 0;

    // compare rotation matrix 

    try {
        amoebaGpu->psRotationMatrix->Download();
        unsigned int numberOfEntries                                  = amoebaGpu->psRotationMatrix->_length/9;
        float* rotationMatrices                                       = reinterpret_cast<float*>(amoebaGpu->psRotationMatrix->_pSysData);
        std::vector< std::vector<double> > expectedRotationMatrices   =  supplementary[AMOEBA_MULTIPOLE_ROTATION_MATRICES];
    
        unsigned int index  = 0;
        unsigned int misses = 0;
        double tolerance    = 5.0e-03;
        numberOfEntries     = expectedRotationMatrices.size() < numberOfEntries ? expectedRotationMatrices.size() : numberOfEntries;
        for( unsigned int ii = 0; ii < numberOfEntries; ii++ ){
            std::vector<double> expectedRotationMatrix = expectedRotationMatrices[ii];
            int rowHit = 0;
            for( unsigned int jj = 0; jj < expectedRotationMatrix.size(); jj++ ){
                double diff = fabs( rotationMatrices[index] - expectedRotationMatrix[jj] );
                if( diff > 0.0 ){
                    diff = 2.0*diff/(fabs( rotationMatrices[index] ) + fabs( expectedRotationMatrix[jj] ) );
                }
                if( diff > tolerance ){
                    misses++;
                    if( misses == 1 ){
                        (void) fprintf( log, "%s: RotationMatrix tolerance=%10.3e\n", methodName.c_str(), tolerance );
                    }
                    if( !rowHit ){
                        (void) fprintf( log, "%5u ", ii );
                        rowHit = 1;
                    }
                    (void) fprintf( log, "%5u [%15.7e %15.7e %15.7e]\n", jj,
                                    diff, rotationMatrices[index], expectedRotationMatrix[jj] );
                }
                index++;
            }
        }
        if( misses == 0 ){
            (void) fprintf( log, "%u rotation matricies agree to relative tolerance of %10.3e\n", numberOfEntries, tolerance );
            (void) fflush( log );
        } else {
            totalMisses += misses;
        }
    } catch( exception& e ){
        (void) fprintf( log, "Rotation matricies not available %s\n", e.what() );
        (void) fflush( log );
    }
    
    // compare fixed E field
    
    try {
        amoebaGpu->psE_Field->Download();
        amoebaGpu->psE_FieldPolar->Download();
        unsigned int numberOfEntries                                  = amoebaGpu->psE_Field->_length/3;
        float* E_Field                                                = reinterpret_cast<float*>(amoebaGpu->psE_Field->_pSysData);
        float* E_FieldPolar                                           = reinterpret_cast<float*>(amoebaGpu->psE_FieldPolar->_pSysData);
        std::vector< std::vector<double> > expectedEFields            = supplementary[AMOEBA_FIXED_E];
    
        double dipoleConversion       = useOpenMMUnits ? 1.0/AngstromToNm : 1.0;
        unsigned int misses           = 0;
        double tolerance              = 1.0e-03;
        numberOfEntries               = expectedEFields.size() < numberOfEntries ? expectedEFields.size() : numberOfEntries;
        for( unsigned int ii = 0; ii < numberOfEntries; ii++ ){
            std::vector<double> expectedEField = expectedEFields[ii];
            int rowHit = 0;
            for( unsigned int jj = 0; jj < expectedEField.size(); jj++ ){
                double eFieldValue  = (jj < 3) ? E_Field[ii*3+jj] : E_FieldPolar[ii*3+jj-3];
                       eFieldValue *= dipoleConversion;
                double diff         = fabs( eFieldValue - expectedEField[jj] );
                if( diff > 1.0e-04 ){
                    diff = 2.0*diff/(fabs( eFieldValue ) + fabs( expectedEField[jj] ) );
                }
                if( diff > tolerance ){
                    misses++;
                    if( misses == 1 ){
                        (void) fprintf( log, "%s: EField\n", methodName.c_str() );
                    }
                    if( !rowHit ){
                        (void) fprintf( log, "     Row %5u\n", ii );
                        rowHit = 1;
                    }
                    (void) fprintf( log, "         %5u [%15.7e %15.7e %15.7e]\n", jj, diff, eFieldValue, expectedEField[jj] );
                }
            }
        }
        if( misses == 0 ){
            (void) fprintf( log, "%u fixed-E fields agree to relative tolerance of %10.3e\n", numberOfEntries, tolerance); 
            (void) fflush( log );
        } else {
            totalMisses += misses;
        }
    } catch( exception& e ){
        (void) fprintf( log, "Fixed-E fields not available %s\n", e.what() );
        (void) fflush( log );
    }
    
    try {
        // compare induced dipoles
    
        amoebaGpu->psInducedDipole->Download();
        amoebaGpu->psInducedDipolePolar->Download();
        unsigned int numberOfEntries                                  = amoebaGpu->psInducedDipole->_length/3;
        float* inducedDipole                                          = reinterpret_cast<float*>(amoebaGpu->psInducedDipole->_pSysData);
        float* inducedDipolePolar                                     = reinterpret_cast<float*>(amoebaGpu->psInducedDipolePolar->_pSysData);
        std::vector< std::vector<double> > expectedInducedDipoles     = supplementary[AMOEBA_INDUCDED_DIPOLES];
    
        unsigned int misses           = 0;
        double tolerance              = 1.0e-03;
        double dipoleConversion       = useOpenMMUnits ? 1.0/AngstromToNm : 1.0;
        numberOfEntries               = expectedInducedDipoles.size() < numberOfEntries ? expectedInducedDipoles.size() : numberOfEntries;
        for( unsigned int ii = 0; ii < numberOfEntries; ii++ ){

            std::vector<double> expectedInducedDipole = expectedInducedDipoles[ii];
            int rowHit                                = 0;

            for( unsigned int jj = 0; jj < expectedInducedDipole.size(); jj++ ){
                double inducedDipoleValue  = (jj < 3) ? inducedDipole[ii*3+jj] : inducedDipolePolar[ii*3+jj-3];
                       inducedDipoleValue *= dipoleConversion;
                double diff         = fabs( inducedDipoleValue - expectedInducedDipole[jj] );
                if( diff > 1.0e-04 ){
                    diff = 2.0*diff/(fabs( inducedDipoleValue ) + fabs( expectedInducedDipole[jj] ) );
                }

                int printDipole = 0;
                if( diff > tolerance ){
                    misses++;
                    printDipole = 1;
                }
                if( misses == 1 && printDipole ){
                    (void) fprintf( log, "%s: induced dipoles\n", methodName.c_str() );
                }
                if( printDipole ){
                    if( !rowHit ){
                        (void) fprintf( log, "     Row %5u\n", ii );
                        rowHit = 1;
                    }
                    (void) fprintf( log, "         %5u [%15.7e %15.7e %15.7e]\n", jj, diff, inducedDipoleValue, expectedInducedDipole[jj] );
                }
            }
        }
        if( misses == 0 ){
            (void) fprintf( log, "%u induced dipoles agree to relative tolerance of %10.3e\n", numberOfEntries, tolerance);
            (void) fflush( log );
        } else {
            totalMisses += misses;
        }
    } catch( exception& e ){
        (void) fprintf( log, "Induced dipoles not available %s\n", e.what() ); 
        (void) fflush( log );
    }

}

void calculateBorn1( System& amoebaSystem, std::vector<Vec3>& tinkerCoordinates, FILE* log ) {
/*
// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "calculateBorn1";
 
// ---------------------------------------------------------------------------------------
 
    System* system  = new System();
    int hit         = 0;
    for( int ii = 0; ii < amoebaSystem.getNumForces() && hit == 0; ii++ ){
        const Force& force = amoebaSystem.getForce(ii);
        try {
            const AmoebaGeneralizedKirkwoodForce& amoebaGeneralizedKirkwoodForce = dynamic_cast<const AmoebaGeneralizedKirkwoodForce&>(force);
            hit = 1;
            GBSAOBCForce* gbsa = new GBSAOBCForce();
            system->addForce( gbsa );
            for( int jj = 0; jj < amoebaGeneralizedKirkwoodForce.getNumParticles(); jj++ ){
                double charge, radius, scalingFactor;
                amoebaGeneralizedKirkwoodForce.getParticleParameters(jj, charge, radius, scalingFactor);
                radius *= 0.1;
                gbsa->addParticle( charge, radius, scalingFactor);
                system->addParticle( 1.0 );
            }
            gbsa->setSoluteDielectric( amoebaGeneralizedKirkwoodForce.getSoluteDielectric() );
            gbsa->setSolventDielectric( amoebaGeneralizedKirkwoodForce.getSolventDielectric() );
        } catch( std::bad_cast ){
        }    
    }
    if( hit == 0 ){
        (void) fprintf( log, "%s: AmoebaGeneralizedKirkwoodForce not found\n", methodName.c_str() );
    } else {
        (void) fprintf( log, "%s: AmoebaGeneralizedKirkwoodForce found\n", methodName.c_str() );
    }
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(*system, integrator, Platform::getPlatformByName( "Reference"));
    std::vector<Vec3> coordinates;
    coordinates.resize( tinkerCoordinates.size() );
    for( unsigned int ii = 0; ii < tinkerCoordinates.size(); ii++ ){
        Vec3 coordinate = tinkerCoordinates[ii];
        coordinates[ii] = Vec3( coordinate[0]*0.1, coordinate[1]*0.1, coordinate[2]*0.1 );
    }
    context.setPositions(coordinates);

    State state                            = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3> forces         = state.getForces();

    if( log ){
        (void) fprintf( log, "%s: energy=%15.7e\n", methodName.c_str(), state.getPotentialEnergy() );
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e]\n", ii, forces[ii][0], forces[ii][1], forces[ii][2] );
        }
        (void) fflush( log );
    }
*/
}

/**---------------------------------------------------------------------------------------
 * Get integrator
 * 
 * @param  inputArgumentMap     StringString Map w/ argumement values
 * @param  log                  optional logging reference
 *
 * @return OpenMM integrator
 *
   --------------------------------------------------------------------------------------- */

Integrator* getIntegrator( MapStringString& inputArgumentMap, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "getIntegrator";

   std::string integratorName               = "LangevinIntegrator";
   double timeStep                          = 0.001;
   double friction                          = 91.0;
   double temperature                       = 300.0;
   double shakeTolerance                    = 1.0e-05;
   double errorTolerance                    = 1.0e-05;
   int randomNumberSeed                     = 1993;

// ---------------------------------------------------------------------------------------

    setStringFromMap( inputArgumentMap, "integrator",       integratorName );
    setDoubleFromMap( inputArgumentMap, "timeStep",         timeStep );
    setDoubleFromMap( inputArgumentMap, "friction",         friction );
    setDoubleFromMap( inputArgumentMap, "temperature",      temperature );
    setDoubleFromMap( inputArgumentMap, "shakeTolerance",   shakeTolerance );
    setDoubleFromMap( inputArgumentMap, "errorTolerance",   errorTolerance );
    setIntFromMap(    inputArgumentMap, "randomNumberSeed", randomNumberSeed );

    // create integrator 
    
    Integrator* integrator;

    if( integratorName == "VerletIntegrator" ){
        integrator = new VerletIntegrator( timeStep );
    } else if( integratorName == "VariableVerletIntegrator" ){
        integrator = new VariableVerletIntegrator( errorTolerance );
    } else if( integratorName == "BrownianIntegrator" ){
        integrator = new BrownianIntegrator( temperature, friction, timeStep );
    } else if( integratorName == "LangevinIntegrator" ){
        integrator                                = new LangevinIntegrator( temperature, friction, timeStep );
        LangevinIntegrator* langevinIntegrator    = dynamic_cast<LangevinIntegrator*>(integrator);
        langevinIntegrator->setRandomNumberSeed( randomNumberSeed );
    } else if( integratorName == "VariableLangevinIntegrator" ){
        integrator                                        = new VariableLangevinIntegrator( temperature, friction, errorTolerance );
        VariableLangevinIntegrator* langevinIntegrator    = dynamic_cast<VariableLangevinIntegrator*>(integrator);
        langevinIntegrator->setRandomNumberSeed( randomNumberSeed );
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "%s  integrator=<%s> not recognized.\n", methodName.c_str(), integratorName.c_str() );
       if( log ){
          (void) fprintf( log , "%s", buffer );
          (void) fflush( log );
       }
       throwException(__FILE__, __LINE__, buffer );
       return NULL;
    }    

    integrator->setConstraintTolerance( shakeTolerance );
    
   return integrator;
}

/**---------------------------------------------------------------------------------------
 * Get integrator type
 * 
 * @param  integrator 
 *
 * @return name or "NotFound"
 *
   --------------------------------------------------------------------------------------- */

static std::string getIntegratorName( Integrator* integrator ){

// ---------------------------------------------------------------------------------------

//   static const std::string methodName      = "getIntegratorName";

// ---------------------------------------------------------------------------------------

   // LangevinIntegrator

   try {
      LangevinIntegrator& langevinIntegrator = dynamic_cast<LangevinIntegrator&>(*integrator);
      return  "LangevinIntegrator";
   } catch( std::bad_cast ){
   }

   // VariableLangevinIntegrator

   try {
      VariableLangevinIntegrator& langevinIntegrator = dynamic_cast<VariableLangevinIntegrator&>(*integrator);
      return "VariableLangevinIntegrator";
   } catch( std::bad_cast ){
   }

   // VerletIntegrator

   try {
      VerletIntegrator& verletIntegrator = dynamic_cast<VerletIntegrator&>(*integrator);
      return "VerletIntegrator";
   } catch( std::bad_cast ){
   }
    
   // VariableVerletIntegrator

   try {
      VariableVerletIntegrator & variableVerletIntegrator = dynamic_cast<VariableVerletIntegrator&>(*integrator);
      return "VariableVerletIntegrator";
   } catch( std::bad_cast ){
   }
    
   // BrownianIntegrator

   try {
      BrownianIntegrator& brownianIntegrator = dynamic_cast<BrownianIntegrator&>(*integrator);
      return "BrownianIntegrator";
   } catch( std::bad_cast ){
   }
    
   return "NotFound";
}

/**---------------------------------------------------------------------------------------
 * Print Integrator info to log
 * 
 * @param integrator   integrator
 * @param log          optional log reference
 *
 * @return DefaultReturnValue
 *
   --------------------------------------------------------------------------------------- */

static void printIntegratorInfo( Integrator& integrator, FILE* log ){
    
// ---------------------------------------------------------------------------------------

   //static const std::string methodName    = "printIntegratorInfo";

// ---------------------------------------------------------------------------------------

   std::string integratorName           = getIntegratorName( &integrator );
   (void) fprintf( log, "Integrator=%s ShakeTol=%.3e ", 
                   integratorName.c_str(), integrator.getConstraintTolerance() );

   // stochastic integrators (seed, friction, temperature)

   if( integratorName == "LangevinIntegrator" || integratorName == "VariableLangevinIntegrator"  ||
       integratorName ==  "BrownianIntegrator" ){

      double temperature = 300.0;
      double friction    = 100.0;
      int seed           = 0;

      if( integratorName == "LangevinIntegrator" ){
         LangevinIntegrator&  langevinIntegrator         = dynamic_cast<LangevinIntegrator&>(integrator);
         temperature                                     = langevinIntegrator.getTemperature();
         friction                                        = langevinIntegrator.getFriction();
         seed                                            = langevinIntegrator.getRandomNumberSeed();
      } else if( integratorName == "VariableLangevinIntegrator" ){
         VariableLangevinIntegrator&  langevinIntegrator = dynamic_cast<VariableLangevinIntegrator&>(integrator);
         temperature                                     = langevinIntegrator.getTemperature();
         friction                                        = langevinIntegrator.getFriction();
         seed                                            = langevinIntegrator.getRandomNumberSeed();
      } else if( integratorName == "BrownianIntegrator" ){
         BrownianIntegrator& brownianIntegrator          = dynamic_cast<BrownianIntegrator&>(integrator);
         temperature                                     = brownianIntegrator.getTemperature();
         friction                                        = brownianIntegrator.getFriction();
         seed                                            = brownianIntegrator.getRandomNumberSeed();
      }
   
      (void) fprintf( log, "T=%.3f friction=%.3f seed=%d ", temperature, friction, seed );
   }

   // variable integrators -- error tolerance

   if( integratorName == "VariableLangevinIntegrator" || integratorName== "VariableVerletIntegrator" ){
      double errorTolerance = 0.0;
      if( integratorName == "VariableLangevinIntegrator" ){
         VariableLangevinIntegrator& langevinIntegrator          = dynamic_cast<VariableLangevinIntegrator&>(integrator);
         errorTolerance                                          = langevinIntegrator.getErrorTolerance();
      } else {
         VariableVerletIntegrator& verletIntegrator              = dynamic_cast<VariableVerletIntegrator&>(integrator);
         errorTolerance                                          = verletIntegrator.getErrorTolerance();
      }
      (void) fprintf( log, "Error tolerance=%.3e\n", errorTolerance );
   } else {
      (void) fprintf( log, "Step size=%12.3e\n", integrator.getStepSize() );
   }

   (void) fflush( log );

   return;
}

/**---------------------------------------------------------------------------------------

   Set the velocities/positions of context2 to those of context1

   @param context1                 context1 
   @param context2                 context2 

   @return 0

   --------------------------------------------------------------------------------------- */

static int synchContexts( const Context& context1, Context& context2 ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName  = "\n_synchContexts: ";

   // ---------------------------------------------------------------------------------------

   const State state                       = context1.getState(State::Positions | State::Velocities);
   const std::vector<Vec3>& positions      = state.getPositions();
   const std::vector<Vec3>& velocities     = state.getVelocities();

   context2.setPositions( positions );
   context2.setVelocities( velocities );

   return DefaultReturnValue;
}

/**---------------------------------------------------------------------------------------
      
   Get statistics of elements in array
   
   @param array               array to collect stats
   @param statistics          statistics of array
      index = 0   mean
      index = 1   stddev
      index = 2   min 
      index = 3   index of min value 
      index = 4   max
      index = 5   index of max value 
      index = 6   size of array

   @return DefaultReturnValue
      
   --------------------------------------------------------------------------------------- */

static void getStatistics( const std::vector<double> & array,  std::vector<double> & statistics ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "getStatistics";

   static const int mean         = 0;
   static const int stddev       = 1;
   static const int min          = 2;
   static const int minIndex     = 3;
   static const int max          = 4;
   static const int maxIndex     = 5;
   static const int size         = 6;

   // ---------------------------------------------------------------------------------------

   // initialize stat array

   statistics.resize( 10 );
   for( unsigned int jj = 0; jj < statistics.size(); jj++ ){
      statistics[jj] = 0.0;
   }
   statistics[min] =  1.0e+30;
   statistics[max] = -1.0e+30;

   // collect stats

   int index       = 0;
   for( std::vector<double>::const_iterator ii = array.begin(); ii != array.end(); ii++ ){

      // first/second moments

      statistics[mean]     += *ii;
      statistics[stddev]   += (*ii)*(*ii);

      // min/max

      if( *ii < statistics[min] ){
         statistics[min]      = *ii;
         statistics[minIndex] = index;
      }
      if( *ii > statistics[max] ){
         statistics[max]      = *ii;
         statistics[maxIndex] = index;
      }
      index++;
   }

   // compute mean & std dev

   double arraySz      = (double) index;
   statistics[size]    = arraySz;
   if( index ){
      statistics[mean]   /= arraySz;
      statistics[stddev]  = statistics[stddev] - arraySz*statistics[mean]*statistics[mean];
      if( index > 1 ){
         statistics[stddev]  = std::sqrt( statistics[stddev] / ( arraySz - 1.0 ) );
      }
   }

   return;
}

/**---------------------------------------------------------------------------------------
      
   Cret a OpenMM context
   
   @param amoebaTinkerParameterFileName   parameter file name
   @param forceMap                        StringInt map[Force] = 1 or 0 (include/not include)
   @param useOpenMMUnits                  if set, convert to OpenMM units (kJ/nm)
   @param inputArgumentMap                StringString map w/ command-line arguments/values
   @param supplementary                   output of supplementary info (rotation matrices, ...)
   @param tinkerForces                    Tinker calculated forces
   @param tinkerEnergies                  Tinker calculated energies
   @param log                             optional file logging reference

   @return OpenMM context 
      
   --------------------------------------------------------------------------------------- */

Context* createContext( const std::string& amoebaTinkerParameterFileName, MapStringInt& forceMap,
                        int useOpenMMUnits, MapStringString& inputArgumentMap, MapStringVectorOfVectors& supplementary,
                        MapStringVec3& tinkerForces, MapStringDouble& tinkerEnergies, FILE* log ) {

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "createContext";
    int  cudaDevice                          = -1;
 
// ---------------------------------------------------------------------------------------
 
    setIntFromMap(    inputArgumentMap, "cudaDevice", cudaDevice );

    System* system  = new System();
 
    std::vector<Vec3> coordinates; 
    std::vector<Vec3> velocities; 
    MapStringIntI isPresent = forceMap.find( AMOEBA_GK_FORCE );
    bool gkIsActive;
    if( isPresent != forceMap.end() && isPresent->second != 0 ){
        forceMap[AMOEBA_MULTIPOLE_FORCE] = 1;
        gkIsActive                       = true;
    } else {
        gkIsActive                       = false;
    }

    // read parameters into system and coord/velocities into appropriate arrays
 
    readAmoebaParameterFile( amoebaTinkerParameterFileName, forceMap, *system, coordinates, velocities,
                             tinkerForces, tinkerEnergies, supplementary, useOpenMMUnits, inputArgumentMap, log );
 
    Integrator* integrator = getIntegrator( inputArgumentMap, log );

    Platform& platform = Platform::getPlatformByName("Cuda");
    map<string, string> properties;
    if( getenv("CudaDevice") || cudaDevice > -1 ){
        std::string cudaDeviceStr;
        if( getenv("CudaDevice") ){
            cudaDeviceStr = getenv("CudaDevice");
        } else {
            std::stringstream cudaDeviceStrStr;
            cudaDeviceStrStr << cudaDevice;
            cudaDeviceStr = cudaDeviceStrStr.str();
        }
        properties["CudaDevice"] = cudaDeviceStr;
        if( log ){
            (void) fprintf( log, "Setting Cuda device to %s.\n", cudaDeviceStr.c_str() );
        }
    }
    Context* context = new Context(*system, *integrator, platform, properties);
    context->setPositions(coordinates);

    return context;

}

void checkIntermediateStatesUsingAmoebaTinkerParameterFile( const std::string& amoebaTinkerParameterFileName, MapStringInt& forceMap,
                                                            int useOpenMMUnits, MapStringString& inputArgumentMap,
                                                            FILE* summaryFile,  FILE* log ) {

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "checkIntermediateStatesUsingAmoebaTinkerParameterFile";
    std::string statesFileName               = "states.txt";
 
// ---------------------------------------------------------------------------------------
 
    setStringFromMap( inputArgumentMap, "states",  statesFileName);

    StringVector forceList;
    std::string activeForceNames;
    for( MapStringInt::const_iterator ii = forceMap.begin(); ii != forceMap.end(); ii++ ){
        if( ii->second ){
            forceList.push_back( ii->first );
            activeForceNames += ii->first + ":";
        }
    }
    if( forceList.size() >= 11 ){
        activeForceNames =ALL_FORCES;
    }
  
    MapStringVec3 tinkerForces;
    MapStringDouble tinkerEnergies;
    MapStringVectorOfVectors supplementary;
 
    MapStringIntI isPresent = forceMap.find( AMOEBA_GK_FORCE );
    bool gkIsActive;
    if( isPresent != forceMap.end() && isPresent->second != 0 ){
        forceMap[AMOEBA_MULTIPOLE_FORCE] = 1;
        gkIsActive                       = true;
    } else {
        gkIsActive                       = false;
    }

    // read parameters into system and coord/velocities into appropriate arrays
    // and create context
 
    Context* context = createContext( amoebaTinkerParameterFileName, forceMap,
                                      useOpenMMUnits, inputArgumentMap, supplementary, tinkerForces, tinkerEnergies, log );

    StringVectorVector fileContents;
    if( readFile( statesFileName, fileContents, log ) ){
        (void) fprintf( stderr, "%s: File %s not read.\n", methodName.c_str(), statesFileName.c_str() );
        (void) fflush( stderr );
        exit(-1);
    }
    unsigned int lineIndex  = 0;
    unsigned int stateIndex = 0;

 
    while( lineIndex < (fileContents.size()-1) ){

        int numberOfAtoms = atoi( fileContents[lineIndex++][0].c_str() );

        stateIndex++;
        std::vector<Vec3> coordinates;
        coordinates.resize( numberOfAtoms );
        int skip = 0;
        for( int ii = 0; ii < numberOfAtoms; ii++ ){
            StringVector& stateTokenArray = fileContents[lineIndex++];
            if( stateTokenArray[1] == "nan" || stateTokenArray[2] == "nan" || stateTokenArray[3] == "nan" ){
                skip = 1;
            } else {
                coordinates[ii] = Vec3( atof( stateTokenArray[1].c_str() ), 
                                        atof( stateTokenArray[2].c_str() ), 
                                        atof( stateTokenArray[3].c_str() ) ); 
            }
        }
        if( skip && log ){
            (void) fprintf( log, "Skipping state=%u line=%u\n", stateIndex, lineIndex );
        } else if( !skip ){
            if( log ){
                (void) fprintf( log, "State=%u coordinates=%u\n", stateIndex, static_cast<unsigned int>(coordinates.size()) );
            }
            context->setPositions( coordinates );
            State state                            = context->getState(State::Forces | State::Energy);
            System& system                         = context->getSystem();
            double potentialEnergy                 = state.getPotentialEnergy();
           
            if( summaryFile ){
                int lastIndex = coordinates.size() - 1;
                FILE* filePtr = summaryFile;
                (void) fprintf( filePtr, "%8u %15.7e   %30s  [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e]\n",
                                stateIndex, potentialEnergy, activeForceNames.c_str(),
                                coordinates[0][0],         coordinates[0][1],         coordinates[0][2],
                                coordinates[lastIndex][0], coordinates[lastIndex][1], coordinates[lastIndex][2] );
                (void) fflush( filePtr );
            }
            if( log ){
                std::vector<Vec3> forces     = state.getForces();
                int lastIndex                = forces.size() - 1;
                FILE* filePtr                = log;
                (void) fprintf( filePtr, "%8u %15.7e   %30s  [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e]\n\nForces\n",
                                stateIndex, potentialEnergy, activeForceNames.c_str(),
                                coordinates[0][0],         coordinates[0][1],         coordinates[0][2],
                                coordinates[lastIndex][0], coordinates[lastIndex][1], coordinates[lastIndex][2] );

                for( int ii = 0; ii < numberOfAtoms; ii++ ){
                    double forceNorm = sqrt( forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2] );
                    (void) fprintf( filePtr, "%8d %15.7e [%15.7e %15.7e %15.7e] %s\n",
                                    ii, forceNorm, forces[ii][0], forces[ii][1], forces[ii][2],
                                    (forceNorm > 1.0e+06 ? "YYY" : "" ) );
                }
                (void) fflush( filePtr );
            }
        }
    }
}

void testUsingAmoebaTinkerParameterFile( const std::string& amoebaTinkerParameterFileName, MapStringInt& forceMap,
                                         int useOpenMMUnits, MapStringString& inputArgumentMap,
                                         FILE* summaryFile,  FILE* log ) {

// ---------------------------------------------------------------------------------------

    int applyAssert                          = 0;
    int includeCavityTerm                    = 0;
    double tolerance                         = 1.0e-02;
    static const std::string methodName      = "testUsingAmoebaTinkerParameterFile";
 
// ---------------------------------------------------------------------------------------
 
    setIntFromMap(    inputArgumentMap, "applyAssert",  applyAssert);
    setDoubleFromMap( inputArgumentMap, "tolerance",    tolerance );
    setIntFromMap(    inputArgumentMap, INCLUDE_OBC_CAVITY_TERM,  includeCavityTerm );

    MapStringVec3 tinkerForces;
    MapStringDouble tinkerEnergies;
    MapStringVectorOfVectors supplementary;
 

    MapStringIntI isPresent = forceMap.find( AMOEBA_GK_FORCE );
    bool gkIsActive;
    if( isPresent != forceMap.end() && isPresent->second != 0 ){
        forceMap[AMOEBA_MULTIPOLE_FORCE] = 1;
        gkIsActive                       = true;
    } else {
        gkIsActive                       = false;
    }

    // read parameters into system and coord/velocities into appropriate arrays
    // and create context
 
    Context* context = createContext( amoebaTinkerParameterFileName, forceMap,
                                      useOpenMMUnits, inputArgumentMap, supplementary, tinkerForces, tinkerEnergies, log );

    State state                            = context->getState( State::Positions | State::Forces | State::Energy);
    System& system                         = context->getSystem();
    std::vector<Vec3> coordinates          = state.getPositions();
    std::vector<Vec3> forces               = state.getForces();

    // get list of forces and then accumulate expected energies/forces

    StringVector forceList;
    std::string activeForceNames;
    for( MapStringInt::const_iterator ii = forceMap.begin(); ii != forceMap.end(); ii++ ){
        if( ii->second ){
            if( includeCavityTerm && ii->first == AMOEBA_GK_FORCE ){
                forceList.push_back( AMOEBA_GK_CAVITY_FORCE );
                activeForceNames += AMOEBA_GK_CAVITY_FORCE + ":";
            } else {
                forceList.push_back( ii->first );
                activeForceNames += ii->first + ":";
            }
        }
    }
    if( forceList.size() >= 11 ){
        activeForceNames = ALL_FORCES;
    }
  
    std::vector<Vec3> expectedForces;
    expectedForces.resize( system.getNumParticles() );
    for( int ii = 0; ii < system.getNumParticles(); ii++ ){
        expectedForces[ii][0] = 0.0;
        expectedForces[ii][1] = 0.0;
        expectedForces[ii][2] = 0.0;
    }
    double expectedEnergy = 0.0;

    for( unsigned int ii = 0; ii < forceList.size(); ii++ ){
        expectedEnergy                += tinkerEnergies[forceList[ii]];
        std::vector<Vec3> forces       = tinkerForces[forceList[ii]];
        for( int jj = 0; jj < system.getNumParticles(); jj++ ){
            expectedForces[jj][0] += forces[jj][0];
            expectedForces[jj][1] += forces[jj][1];
            expectedForces[jj][2] += forces[jj][2];
        }
    }

    int showAll                            = 1;
    double energyConversion;
    double forceConversion;
    double coordinateConversion;
    if( useOpenMMUnits ){
        energyConversion     = 1.0/CalToJoule;
        forceConversion      = -energyConversion*AngstromToNm;
        coordinateConversion = 1.0/AngstromToNm;
    } else {
        energyConversion     = 1.0;
        forceConversion      = -energyConversion;
        coordinateConversion = 1.0;
    }

    // output to log and/or summary file

    if( log ){
        std::vector<FILE*> fileList;
        if( log )fileList.push_back( log );
        double cutoffDelta = 0.02;
        for( unsigned int ii = 0; ii < fileList.size(); ii++ ){
            FILE* filePtr = fileList[ii];
            (void) fprintf( filePtr, "\n" );
            (void) fprintf( filePtr, "%s: conversion factors %15.7e %15.7e %12.3f tolerance=%15.7e %s\n",
                            methodName.c_str(), energyConversion, forceConversion, coordinateConversion, tolerance, amoebaTinkerParameterFileName.c_str() );
            double deltaE = fabs( expectedEnergy   -       energyConversion*state.getPotentialEnergy());
            double denom  = fabs( expectedEnergy ) + fabs( energyConversion*state.getPotentialEnergy());
            if( denom > 0.0 )deltaE *= 2.0/denom;
            (void) fprintf( filePtr, "expectedE        %10.3e %15.7e %15.7e %20s %30s\n",
                            deltaE, expectedEnergy, energyConversion*state.getPotentialEnergy(),
                            amoebaTinkerParameterFileName.c_str(), activeForceNames.c_str() );
            (void) fprintf( filePtr, "%s: %u %u Active forces: %s\n",
                            methodName.c_str(), static_cast<unsigned int>(expectedForces.size()), static_cast<unsigned int>(forces.size()), activeForceNames.c_str() );
            double maxRelativeDelta            = -1.0e+30;
            unsigned int maxRelativeDeltaIndex = -1;
            for( unsigned int ii = 0; ii < forces.size(); ii++ ){
                double normF1          = std::sqrt( (expectedForces[ii][0]*expectedForces[ii][0]) + (expectedForces[ii][1]*expectedForces[ii][1]) + (expectedForces[ii][2]*expectedForces[ii][2]) );
                double normF2          = std::sqrt( (forces[ii][0]*forces[ii][0]) + (forces[ii][1]*forces[ii][1]) + (forces[ii][2]*forces[ii][2]) );
                       normF2         *= fabs( forceConversion );
                double delta           = fabs( normF1 - normF2 );
                double sumNorms        = 0.5*(normF1 + normF2);
                double relativeDelta   = sumNorms > 0.0 ? fabs( normF1 - normF2 )/sumNorms : 0.0;
                bool badMatch          = (cutoffDelta < relativeDelta) && (sumNorms > 0.1) ? true : false;
                     badMatch          = badMatch || (normF1 == 0.0 && normF2 > 0.0) || (normF2 == 0.0 && normF1 > 0.0);
                if( badMatch || showAll ){
                    (void) fprintf( filePtr, "%6u %10.3e %10.3e [%15.7e %15.7e %15.7e]   [%15.7e %15.7e %15.7e]  %s\n", ii, relativeDelta, delta,
                                    expectedForces[ii][0], expectedForces[ii][1], expectedForces[ii][2],
                                    forceConversion*forces[ii][0], forceConversion*forces[ii][1], forceConversion*forces[ii][2],
                                    ( (showAll && badMatch) ? " XXX" : "") );
                    if( ( (maxRelativeDelta < relativeDelta) && (sumNorms > 0.1)) ){
                        maxRelativeDelta      =  relativeDelta;
                        maxRelativeDeltaIndex = ii;
                    }
                }
            }
            (void) fprintf( filePtr, "maxRelativeDelta %10.3e at %6u %20s %30s\n", maxRelativeDelta, maxRelativeDeltaIndex, amoebaTinkerParameterFileName.c_str(), activeForceNames.c_str() );

            // get box dimensions and bond distance for atom 0 

            double box[2][3];
            for( unsigned int jj = 0; jj < 3; jj++ ){
                box[0][jj] = coordinates[0][jj];
            }
            
            double minDistToAtom0     = 1.0e+30;
            double nextMinDistToAtom0 = 1.0e+30;
            for( unsigned int ii = 1; ii < coordinates.size(); ii++ ){
                double dist = 0.0;
                for( unsigned int jj = 0; jj < 3; jj++ ){
                    if( box[0][jj] > coordinates[ii][jj] ){
                        box[0][jj] = coordinates[ii][jj];
                    }
                    if( box[1][jj] < coordinates[ii][jj] ){
                        box[1][jj] = coordinates[ii][jj];
                    }
                    dist += (coordinates[ii][jj] - coordinates[0][jj])*(coordinates[ii][jj] - coordinates[0][jj]);
                }
                if( dist < minDistToAtom0 ){
                    nextMinDistToAtom0 = minDistToAtom0;
                    minDistToAtom0     = dist;
                }
            }
                
            (void) fprintf( filePtr, "Mindist atom 0 (in A) %10.3e %10.3e Box [%15.7e %15.7e] [%15.7e %15.7e] [%15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n",
                            sqrt( minDistToAtom0 )*coordinateConversion,
                            sqrt( nextMinDistToAtom0 )*coordinateConversion,
                            coordinateConversion*box[0][0], coordinateConversion*box[1][0],
                            coordinateConversion*box[0][1], coordinateConversion*box[1][1],
                            coordinateConversion*box[0][2], coordinateConversion*box[1][2],
                            coordinateConversion*(box[1][0] - box[0][0]),
                            coordinateConversion*(box[1][1] - box[0][1]),
                            coordinateConversion*(box[1][2] - box[0][2]) );

            unsigned int maxPrint = 5;
/*
            (void) fprintf( filePtr, "Sample raw coordinates (w/o conversion) %8u\n", static_cast<unsigned int>(coordinates.size()) );
            for( unsigned int ii = 0; ii < coordinates.size(); ii++ ){
                (void) fprintf( filePtr, "%8u [%16.7f %16.7f %16.7f]\n", ii,
                            coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] );
                if( ii == maxPrint && (coordinates.size()- maxPrint) > ii ){
                    ii = coordinates.size() - maxPrint - 1;
                }
            } 
*/
            (void) fflush( filePtr );
        }
    }

    if( summaryFile ){
        std::vector<FILE*> fileList;
        if( summaryFile )fileList.push_back( summaryFile );
        for( unsigned int ii = 0; ii < fileList.size(); ii++ ){

            FILE* filePtr = fileList[ii];

            double deltaE = fabs( expectedEnergy   -       energyConversion*state.getPotentialEnergy());
            double denom  = fabs( expectedEnergy ) + fabs( energyConversion*state.getPotentialEnergy());
            if( denom > 0.0 )deltaE *= 2.0/denom;

            double maxRelativeDelta            = -1.0e+30;
            unsigned int maxRelativeDeltaIndex = -1;
            for( unsigned int ii = 0; ii < forces.size(); ii++ ){
                double normF1          = std::sqrt( (expectedForces[ii][0]*expectedForces[ii][0]) + (expectedForces[ii][1]*expectedForces[ii][1]) + (expectedForces[ii][2]*expectedForces[ii][2]) );
                double normF2          = std::sqrt( (forces[ii][0]*forces[ii][0]) + (forces[ii][1]*forces[ii][1]) + (forces[ii][2]*forces[ii][2]) );
                       normF2         *= fabs( forceConversion );
                double delta           = fabs( normF1 - normF2 );
                double sumNorms        = 0.5*(normF1 + normF2);
                double relativeDelta   = sumNorms > 0.0 ? fabs( normF1 - normF2 )/sumNorms : 0.0;
                if( ( (maxRelativeDelta < relativeDelta) && (sumNorms > 0.1)) || showAll ){
                    if( ( (maxRelativeDelta < relativeDelta) && (sumNorms > 0.1)) ){
                        maxRelativeDelta      =  relativeDelta;
                        maxRelativeDeltaIndex = ii;
                    }
                }
            }
            (void) fprintf( filePtr, "%40s maxRelF/E %10.3e %10.3e E[%15.7e %15.7e] %20s %d\n", activeForceNames.c_str(), maxRelativeDelta,
                            deltaE, expectedEnergy, energyConversion*state.getPotentialEnergy(), amoebaTinkerParameterFileName.c_str(), useOpenMMUnits );
            (void) fflush( filePtr );
        }
    }

    if( gkIsActive == false ){
        isPresent = forceMap.find( AMOEBA_MULTIPOLE_FORCE );
        if( isPresent != forceMap.end() && isPresent->second != 0 ){
             checkIntermediateMultipoleQuantities( context, supplementary, useOpenMMUnits, log );
        }
    }

    if( applyAssert ){
        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            forces[ii][0] *= forceConversion;
            forces[ii][1] *= forceConversion;
            forces[ii][2] *= forceConversion;
            ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
        }
        ASSERT_EQUAL_TOL( expectedEnergy, energyConversion*state.getPotentialEnergy(), tolerance );
    }

    if( log ){
        (void) fprintf( log, "No issues w/ tolerance=%10.3e\n", tolerance );
        (void) fflush( log );
    }
}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

void testEnergyForcesConsistent( std::string parameterFileName, MapStringInt& forceMap, int useOpenMMUnits, 
                                 MapStringString& inputArgumentMap,
                                 FILE* log, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

   int applyAssertion                     = 1;
   double delta                           = 1.0e-04;
   double tolerance                       = 0.01;
  
   static const std::string methodName    = "checkEnergyForceConsistent";

// ---------------------------------------------------------------------------------------

    MapStringVectorOfVectors supplementary;
    MapStringVec3 tinkerForces;
    MapStringDouble tinkerEnergies;

    Context* context = createContext( parameterFileName, forceMap, useOpenMMUnits, inputArgumentMap, supplementary,
                                      tinkerForces, tinkerEnergies, log );

    setIntFromMap(    inputArgumentMap, "applyAssert",             applyAssertion  );
    setDoubleFromMap( inputArgumentMap, "energyForceDelta",        delta           );
    setDoubleFromMap( inputArgumentMap, "energyForceTolerance",    tolerance       );
 
    StringVector forceStringArray;
    System& system = context->getSystem();
    getForceStrings( system, forceStringArray, log );
 
    if( log ){
        (void) fprintf( log, "%s delta=%.3e tolerance=%.3e applyAssertion=%d\n", methodName.c_str(), delta, tolerance, applyAssertion );
        (void) fprintf( log, "\nForces:\n" );
        for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
           (void) fprintf( log, "   %s\n", (*ii).c_str() );
        }
        (void) fflush( log );
    }
 
    int returnStatus                       = 0;
 
    // get positions, forces and potential energy
 
    int types                              = State::Positions | State::Velocities | State::Forces | State::Energy;
 
    State state                            = context->getState( types );
 
    std::vector<Vec3> coordinates          = state.getPositions();
    std::vector<Vec3> velocities           = state.getVelocities();
    std::vector<Vec3> forces               = state.getForces();
    double kineticEnergy                   = state.getKineticEnergy();
    double potentialEnergy                 = state.getPotentialEnergy();
 
    // compute norm of force
 
    double forceNorm         = 0.0;
    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forceNorm += forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2];
    }
 
    // check norm is not nan
 
    if( isNanOrInfinity( forceNorm ) ){ 
        if( log ){
            (void) fprintf( log, "%s norm of force is nan -- aborting.\n", methodName.c_str() );
            unsigned int hitNan = 0;
            for( unsigned int ii = 0; (ii < forces.size()) && (hitNan < 10); ii++ ){
   
               if( isNanOrInfinity( forces[ii][0] ) ||
                   isNanOrInfinity( forces[ii][1] ) ||
                   isNanOrInfinity( forces[ii][2] ) )hitNan++;
   
                (void) fprintf( log, "%6u x[%15.7e %15.7e %15.7e] f[%15.7e %15.7e %15.7e]\n", ii,
                                coordinates[ii][0], coordinates[ii][1], coordinates[ii][2],
                                forces[ii][0], forces[ii][1], forces[ii][2] );
            }
            char buffer[1024];
            (void) sprintf( buffer, "%s : nans detected -- aborting.\n", methodName.c_str() );
            throwException(__FILE__, __LINE__, buffer );
        }    
    }
 
    forceNorm = std::sqrt( forceNorm );
 
     if( forceNorm <= 0.0 ){
        if( log ){
           (void) fprintf( log, "%s norm of force is <= 0 norm=%.3e\n", methodName.c_str(), forceNorm );
           (void) fflush( log );
        }
        return;
    }
  
    // take step in direction of energy gradient
 
    double step = delta/forceNorm;
    std::vector<Vec3> perturbedPositions; 
    perturbedPositions.resize( forces.size() );
    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        perturbedPositions[ii] = Vec3( coordinates[ii][0] - step*forces[ii][0], coordinates[ii][1] - step*forces[ii][1], coordinates[ii][2] - step*forces[ii][2] ); 
    }
 
    context->setPositions( perturbedPositions );
 
    // get new potential energy
 
    state    = context->getState( types );
 
    // report energies
 
    double perturbedPotentialEnergy = state.getPotentialEnergy();
    double deltaEnergy              = ( perturbedPotentialEnergy - potentialEnergy )/delta;
    double difference               = fabs( deltaEnergy - forceNorm );
    double denominator              = forceNorm; 
    if( denominator > 0.0 ){
        difference /= denominator;
    }
 
    if( log ){
        (void) fprintf( log, "%s  difference=%14.8e dE=%14.8e Pe2/1 [%16.10e %16.10e] delta=%10.4e nrm=%16.10e\n",
                        methodName.c_str(), difference, deltaEnergy, perturbedPotentialEnergy,
                        potentialEnergy, delta, forceNorm );
        (void) fflush( log );
    }
    if( summaryFile ){
        std::string forceString;
        if( forceStringArray.size() > 11 ){
           forceString = "All ";
        } else {
           for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
              forceString += (*ii) + " ";
           }
        }
        if( forceString.size() < 1 ){
           forceString = "NA";
        }
        (void) fprintf( summaryFile, "EFCnstnt %30s %15.7e dE[%14.6e %15.7e] E[%15.7e %15.7e FNorm %15.7e Delta %15.7e %20s %s\n",
                        forceString.c_str(), difference, deltaEnergy, forceNorm, 
                        potentialEnergy, perturbedPotentialEnergy, forceNorm, delta, parameterFileName.c_str(), context->getPlatform().getName().c_str() );
        (void) fflush( summaryFile );
    }
 
    if( applyAssertion ){
        ASSERT( difference < tolerance );
        if( log ){
           (void) fprintf( log, "\n%s passed\n", methodName.c_str() );
           (void) fflush( log );
        }
    }

    delete context;

    return;

}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

void testEnergyForceByFiniteDifference( std::string parameterFileName, MapStringInt& forceMap, int useOpenMMUnits, 
                                        MapStringString& inputArgumentMap,
                                        FILE* log, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

   int applyAssertion                     = 1;
   double energyForceDelta                = 1.0e-04;
   double tolerance                       = 0.01;
  
   static const std::string methodName    = "testEnergyForceByFiniteDifference";

// ---------------------------------------------------------------------------------------

    MapStringVectorOfVectors supplementary;
    MapStringVec3 tinkerForces;
    MapStringDouble tinkerEnergies;

    Context* context = createContext( parameterFileName, forceMap, useOpenMMUnits, inputArgumentMap, supplementary,
                                      tinkerForces, tinkerEnergies, log );

    setIntFromMap(    inputArgumentMap, "applyAssert",             applyAssertion   );
    setDoubleFromMap( inputArgumentMap, "energyForceDelta",        energyForceDelta );
    setDoubleFromMap( inputArgumentMap, "energyForceTolerance",    tolerance        );
 
    StringVector forceStringArray;
    System& system = context->getSystem();
    getForceStrings( system, forceStringArray, log );
 
    if( log ){
        (void) fprintf( log, "%s energyForceDelta=%.3e tolerance=%.3e applyAssertion=%d\n", methodName.c_str(), energyForceDelta, tolerance, applyAssertion );
        (void) fprintf( log, "\nForces:\n" );
        for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
           (void) fprintf( log, "   %s\n", (*ii).c_str() );
        }
        (void) fflush( log );
    }
 
    int returnStatus                       = 0;
 
    // get positions, forces and potential energy
 
    int types                              = State::Positions | State::Velocities | State::Forces | State::Energy;
 
    State state                            = context->getState( types );
 
    std::vector<Vec3> coordinates          = state.getPositions();
    std::vector<Vec3> velocities           = state.getVelocities();
    std::vector<Vec3> forces               = state.getForces();
    double kineticEnergy                   = state.getKineticEnergy();
    double potentialEnergy                 = state.getPotentialEnergy();
 
    // compute norm of force
 
    double forceNorm         = 0.0;
    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forceNorm += forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2];
    }
 
    // check norm is not nan
 
    if( isNanOrInfinity( forceNorm ) ){ 
        if( log ){
            (void) fprintf( log, "%s norm of force is nan -- aborting.\n", methodName.c_str() );
            unsigned int hitNan = 0;
            for( unsigned int ii = 0; (ii < forces.size()) && (hitNan < 10); ii++ ){
   
               if( isNanOrInfinity( forces[ii][0] ) ||
                   isNanOrInfinity( forces[ii][1] ) ||
                   isNanOrInfinity( forces[ii][2] ) )hitNan++;
   
                (void) fprintf( log, "%6u x[%15.7e %15.7e %15.7e] f[%15.7e %15.7e %15.7e]\n", ii,
                                coordinates[ii][0], coordinates[ii][1], coordinates[ii][2],
                                forces[ii][0], forces[ii][1], forces[ii][2] );
            }
            char buffer[1024];
            (void) sprintf( buffer, "%s : nans detected -- aborting.\n", methodName.c_str() );
            throwException(__FILE__, __LINE__, buffer );
        }    
    }
 
    std::vector<Vec3> perturbedPositions; 
    perturbedPositions.resize( forces.size() );
    for( unsigned int ii = 0; ii < coordinates.size(); ii++ ){
        perturbedPositions[ii] = Vec3( coordinates[ii][0], coordinates[ii][1], coordinates[ii][2] ); 
    }
    
    std::vector<double> energyForceDeltas;
    int scanEnergyForceDeltas = 0;
    if( scanEnergyForceDeltas ){
        energyForceDeltas.push_back( 1.0e-02 );
        energyForceDeltas.push_back( 5.0e-03 );
        energyForceDeltas.push_back( 1.0e-03 );
        energyForceDeltas.push_back( 5.0e-04 );
        energyForceDeltas.push_back( 1.0e-04 );
        energyForceDeltas.push_back( 5.0e-05 );
        energyForceDeltas.push_back( 1.0e-05 );
        energyForceDeltas.push_back( 5.0e-06 );
    } else {
        energyForceDeltas.push_back( energyForceDelta );
    }
    for(  unsigned int kk = 0; kk < energyForceDeltas.size(); kk++ ){
        energyForceDelta = energyForceDeltas[kk];
        std::vector<double> relativeDifferenceStatistics;
        for( unsigned int jj = 0; jj < coordinates.size(); jj++ ){
            perturbedPositions[jj][0] += energyForceDelta;
            context->setPositions( perturbedPositions );
     
            // get new potential energy
     
            state    = context->getState( types );
     
            // report energies
     
            double perturbedPotentialEnergy       = state.getPotentialEnergy();
            std::vector<Vec3> perturbedForces     = state.getForces();
            double deltaEnergy                    = ( potentialEnergy - perturbedPotentialEnergy )/energyForceDelta;
            double difference                     = fabs( deltaEnergy - perturbedForces[jj][0]);
            double denominator                    = 0.5*fabs( deltaEnergy ) + fabs( perturbedForces[jj][0] );
            double relativeDifference             = denominator > 0.0 ? difference/denominator : 0.0;
            if( log ){
                (void) fprintf( log, "   %5u fDiff=%14.8e %14.8e dE=[%16.9e %16.9e] delta=%12.1e\n",
                                jj, relativeDifference, difference, deltaEnergy, perturbedForces[jj][0], energyForceDelta);
                (void) fflush( log );
            }
            if( denominator > 1.0e-02 ){
                relativeDifferenceStatistics.push_back( relativeDifference );
            }
            perturbedPositions[jj][0] -= energyForceDelta;
        }
    
        std::vector<double> statistics;
        getStatistics( relativeDifferenceStatistics, statistics );
        if( log ){
            (void) fprintf( log, "Stats on relative diff average=%14.8e stddev=%14.8e max=%16.9e %8.1f %8.1f %12.3e\n",
                            statistics[0], statistics[1], statistics[4], statistics[5], statistics[6], energyForceDelta );
            (void) fflush( log );
        }
    
        if( summaryFile ){
            std::string forceString;
            if( forceStringArray.size() > 11 ){
               forceString = "All ";
            } else {
               for( StringVectorCI ii = forceStringArray.begin(); ii != forceStringArray.end(); ii++ ){
                  forceString += (*ii) + " ";
               }
            }
            if( forceString.size() < 1 ){
               forceString = "NA";
            }
            (void) fprintf( summaryFile, "FD %30s %15.7e %14.6e %15.7e at %18.1f %8.1f delta %15.7e %20s %s\n",
                            forceString.c_str(), statistics[0], statistics[1], statistics[4], statistics[5], statistics[6],
                            energyForceDelta, parameterFileName.c_str(), context->getPlatform().getName().c_str() );
            (void) fflush( summaryFile );
        }
    }
     
/*
    if( applyAssertion ){
        ASSERT( difference < tolerance );
        if( log ){
           (void) fprintf( log, "\n%s passed\n", methodName.c_str() );
           (void) fflush( log );
        }
    }
*/
    delete context;
    return;

}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

System* getCopyOfSystem( System& system, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "getCopyOfSystem";

// ---------------------------------------------------------------------------------------

    System* newSystem    = new System();
    for( int ii = 0; ii < system.getNumParticles(); ii++ ){
        newSystem->addParticle( system.getParticleMass( ii ) );
    }

    for( int ii = 0; ii < system.getNumConstraints(); ii++ ){
        int particle1, particle2;
        double distance;
        system.getConstraintParameters( ii, particle1, particle2, distance );
        newSystem->addConstraint( particle1, particle2, distance );
    }
    return newSystem;
}

/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

double getEnergyForceBreakdown( Context& context, MapStringDouble& mapEnergies, FILE* log ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "getEnergyForceBreakdown";

// ---------------------------------------------------------------------------------------

    int allTypes                       = State::Positions | State::Velocities | State::Forces | State::Energy;

    State state                        = context.getState( allTypes );
    std::vector<Vec3> coordinates      = state. getPositions();
    System& system                     = context.getSystem();

    MapStringForce forceMap;
    getStringForceMap( system, forceMap, log );

    MapStringForceI gkIsPresent        = forceMap.find( AMOEBA_GK_FORCE );
    bool gkIsActive                    = gkIsPresent == forceMap.end() ? false : true;

    double totalEnergy                 = 0.0;
    for( MapStringForceI ii = forceMap.begin(); ii != forceMap.end(); ii++ ){
        Force* force               = ii->second;
        int addForce               = 1;
        if( gkIsActive ){
            if( ii->first == AMOEBA_MULTIPOLE_FORCE ){
                addForce = 0; 
            } else if( ii->first == AMOEBA_GK_FORCE ){
                addForce = 2;
            }
        }
        if( addForce ){
            System* newSystem          = getCopyOfSystem( system, log );
            newSystem->addForce( force );
            if( addForce == 2 ){
                newSystem->addForce( forceMap[AMOEBA_MULTIPOLE_FORCE] );
            }
             
            Platform& platform          = Platform::getPlatformByName( "Cuda");
            platform.setPropertyDefaultValue( "CudaDevice",  "3");
            //Context newContext         = Context( *newSystem, context.getIntegrator(), Platform::getPlatformByName( "Cuda"));
            Context newContext         = Context( *newSystem, context.getIntegrator(), platform );
            newContext.setPositions(coordinates);
            State newState             = newContext.getState( allTypes );
            mapEnergies[ii->first]     = newState.getPotentialEnergy();
            totalEnergy               += newState.getPotentialEnergy();
        }
    }
    return totalEnergy;
}

/**---------------------------------------------------------------------------------------
 * Set velocities based on temperature
 * 
 * @param system       System reference -- retrieve particle masses
 * @param velocities   array of Vec3 for velocities (size must be set)
 * @param temperature  temperature
 * @param log          optional log reference
 *
 * @return DefaultReturnValue
 *
   --------------------------------------------------------------------------------------- */

static void setVelocitiesBasedOnTemperature( const System& system, int seed, std::vector<Vec3>& velocities, double temperature, FILE* log ) {
    
// ---------------------------------------------------------------------------------------

   static const std::string methodName    = "setVelocitiesBasedOnTemperature";
   double randomValues[3];

// ---------------------------------------------------------------------------------------

    if( seed ){
         SimTKOpenMMUtilities::setRandomNumberSeed( static_cast<uint32_t>(seed) );
         if( log ){
             (void) fprintf( log, "%s set random number seed to %d\n", methodName.c_str(), seed );
         }
    }

    // set velocities based on temperature

    double scaledTemperature     = temperature*2.0*BOLTZ;
    double kineticEnergy         = 0.0;
    for( unsigned int ii = 0; ii < velocities.size(); ii++ ){
        double mass               = system.getParticleMass(ii);
        double velocityScale      = std::sqrt( scaledTemperature/mass );
        randomValues[0]           = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        randomValues[1]           = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        randomValues[2]           = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
        velocities[ii]            = Vec3( randomValues[0]*velocityScale, randomValues[1]*velocityScale, randomValues[2]*velocityScale );
        kineticEnergy            += mass*(velocities[ii][0]*velocities[ii][0] + velocities[ii][1]*velocities[ii][1] + velocities[ii][2]*velocities[ii][2]);
    }
    kineticEnergy *= 0.5;
 
    //double degreesOfFreedom  = static_cast<double>(3*velocities.size() - system.getNumConstraints() - 3 );
    double degreesOfFreedom  = static_cast<double>(3*velocities.size() );
    double approximateT      = (kineticEnergy)/(degreesOfFreedom*BOLTZ);
     if( approximateT > 0.0 ){
        double scale             = approximateT > 0.0 ? std::sqrt(temperature/approximateT) : 1.0;
        for( unsigned int ii = 0; ii < velocities.size(); ii++ ){
            velocities[ii][0] *= scale;
            velocities[ii][1] *= scale;
            velocities[ii][2] *= scale;
        }
     }
 
    if( log ){
        double finalKineticEnergy  = 0.0;
        for( unsigned int ii = 0; ii < velocities.size(); ii++ ){
            double mass             = system.getParticleMass(ii);
            finalKineticEnergy     += mass*(velocities[ii][0]*velocities[ii][0] + velocities[ii][1]*velocities[ii][1] + velocities[ii][2]*velocities[ii][2]);
        }
        finalKineticEnergy *= 0.5;
        double finalT       = (finalKineticEnergy)/(degreesOfFreedom*BOLTZ);
    
        (void) fprintf( log, "%s KE=%15.7e ~T=%15.7e desiredT=%15.7e dof=%12.3f final KE=%12.3e T=%12.3e\n",
                        methodName.c_str(), kineticEnergy, approximateT, temperature,
                        degreesOfFreedom, finalKineticEnergy, finalT );
    }
 
 
    return;
}

/** 
 * Write intermediate state to file
 * 
 * @param context                OpenMM context
 * @param intermediateStateFile  file to write to
 * @param log                    optional logging reference
 *
 */

void writeIntermediateStateFile( Context& context, FILE* intermediateStateFile, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName             = "writeIntermediateStateFile";

// ---------------------------------------------------------------------------------------

    if( intermediateStateFile == NULL )return;

    int allTypes                                   = State::Positions | State::Velocities | State::Forces | State::Energy;
    State state                                    = context.getState( allTypes );

    const std::vector<Vec3> positions              = state.getPositions();
    const std::vector<Vec3> velocities             = state.getVelocities();
    const std::vector<Vec3> forces                 = state.getForces();

    (void) fprintf( intermediateStateFile, "%7u %12.3f %15.7e %15.7e %15.7e State (x,v,f)\n",
                    static_cast<unsigned int>(positions.size()), state.getTime(), state.getKineticEnergy(), state.getPotentialEnergy(),
                    state.getKineticEnergy() + state.getPotentialEnergy() );

    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        (void) fprintf( intermediateStateFile, "%7u %15.7e %15.7e %15.7e  %15.7e %15.7e %15.7e  %15.7e %15.7e %15.7e\n", ii,
                         positions[ii][0],  positions[ii][1],  positions[ii][2], 
                        velocities[ii][0], velocities[ii][1], velocities[ii][2], 
                            forces[ii][0],     forces[ii][1],     forces[ii][2] );
    }
    (void) fflush( intermediateStateFile );
    return;
}

/** 
 * Write intermediate state to file
 * 
 * @param context                OpenMM context
 * @param intermediateStateFile  file to write to
 * @param log                    optional logging reference
 *
 */

static void getVerletKineticEnergy( Context& context, double& currentTime, double& potentialEnergy, double& kineticEnergy, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName             = "getVerletKineticEnergy";

// ---------------------------------------------------------------------------------------

    int stateFieldsToRetreive                 = State::Energy | State::Velocities;
    State state                               = context.getState( stateFieldsToRetreive );

    const std::vector<Vec3>& velocitiesI      = state.getVelocities();

    context.getIntegrator().step( 1 );

    State statePlus1                          = context.getState( stateFieldsToRetreive );
    currentTime                               = statePlus1.getTime();
    potentialEnergy                           = statePlus1.getPotentialEnergy();
    const std::vector<Vec3>& velocitiesIPlus1 = statePlus1.getVelocities();

    System& system                            = context.getSystem();
    kineticEnergy                             = 0.0;
    for( unsigned int ii = 0; ii < velocitiesI.size(); ii++ ){
        double velocity = (velocitiesIPlus1[ii][0] + velocitiesI[ii][0])*(velocitiesIPlus1[ii][0] + velocitiesI[ii][0]) +
                          (velocitiesIPlus1[ii][1] + velocitiesI[ii][1])*(velocitiesIPlus1[ii][1] + velocitiesI[ii][1]) +
                          (velocitiesIPlus1[ii][2] + velocitiesI[ii][2])*(velocitiesIPlus1[ii][2] + velocitiesI[ii][2]);
        kineticEnergy  += velocity*system.getParticleMass(ii);
    }
    kineticEnergy *= 0.125;
    //kineticEnergy = statePlus1.getKineticEnergy();

    return;
}

/** 
 * Check for constraint violations 
 * 
 * @param context                OpenMM context
 * @param log                    optional logging reference
 *
 * @return number of violations
 *
 */

static int checkConstraints( Context& context, double shakeTolerance, double& maxViolation, int& maxViolationIndex, FILE* log ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName             = "getVerletKineticEnergy";

// ---------------------------------------------------------------------------------------

    int stateFieldsToRetreive                 = State::Positions;
    State state                               = context.getState( stateFieldsToRetreive );

    const std::vector<Vec3>& positions        = state.getPositions();

    System& system                            = context.getSystem();
    int violationCount                        = 0;
    maxViolation                              = 0.0;
    maxViolationIndex                         = 0;
    for( int ii = 0; ii < system.getNumConstraints(); ii++ ){
        int particle1, particle2;
        double constrainedDistance;
        system.getConstraintParameters( ii, particle1, particle2, constrainedDistance );

        double distance = (positions[particle2][0] - positions[particle1][0])*(positions[particle2][0] - positions[particle1][0]) +
                          (positions[particle2][1] - positions[particle1][1])*(positions[particle2][1] - positions[particle1][1]) +
                          (positions[particle2][2] - positions[particle1][2])*(positions[particle2][2] - positions[particle1][2]);

        double delta    = fabs( sqrt( distance ) - constrainedDistance );
        if( delta > shakeTolerance ){
            violationCount++;
            if( delta > maxViolation ){
                maxViolation      = delta;
                maxViolationIndex = ii;
            }
        }
    }
    return violationCount;
}

/** 
 * Get time of day (implementation different for Linux/Windows
 * 
 * @return time
 *
 */

double getTimeOfDay( void ){

#ifdef WIN32
    static double cycles_per_usec = 0;
    LARGE_INTEGER counter;
 
    if (cycles_per_usec == 0) {
       static LARGE_INTEGER lFreq;
       if (!QueryPerformanceFrequency(&lFreq)) {
          fprintf(stderr, "Unable to read the performance counter frquency!\n");
          return 0;
       }   
 
       cycles_per_usec = 1000000 / ((double) lFreq.QuadPart);
    }   
 
    if (!QueryPerformanceCounter(&counter)) {
       fprintf(stderr,"Unable to read the performance counter!\n");
       return 0;
    }   
 
    double time = ((((double) counter.QuadPart) * cycles_per_usec));
    return time*1.0e-06;
#else
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return static_cast<double>(tv.tv_sec) + 1.0e-06*static_cast<double>(tv.tv_usec);
#endif
}

double getEnergyDrift( std::vector<double>& totalEnergyArray, std::vector<double>& kineticEnergyArray, double degreesOfFreedom, double deltaTime, FILE* log ){

       // total energy constant
 
    std::vector<double> statistics;
    getStatistics( totalEnergyArray, statistics );
 
    std::vector<double> kineticEnergyStatistics;
    getStatistics( kineticEnergyArray, kineticEnergyStatistics );
    double temperature  = 2.0*kineticEnergyStatistics[0]/(degreesOfFreedom*BOLTZ);
    double kT           = temperature*BOLTZ;
 
    // compute stddev in units of kT/dof/ns
 
    double stddevE      = statistics[1]/kT;
           stddevE     /= degreesOfFreedom;
           stddevE     /= deltaTime*0.001;
 
    if( log ){
        (void) fprintf( log, "Simulation results: mean=%15.7e stddev=%15.7e  kT/dof/ns=%15.7e kT=%15.7e T=%12.3f  min=%15.7e  %d max=%15.7e %d\n",
                        statistics[0], statistics[1], stddevE, kT, temperature, statistics[2], (int) (statistics[3] + 0.001), statistics[4], (int) (statistics[5] + 0.001) );
    }

    return stddevE;
}
 
/** 
 * Check that energy and force are consistent
 * 
 * @return DefaultReturnValue or ErrorReturnValue
 *
 */

void testEnergyConservation( std::string parameterFileName, MapStringInt& forceMap, int useOpenMMUnits, 
                             MapStringString& inputArgumentMap,
                             FILE* log, FILE* summaryFile ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "testEnergyConservation";

   // tolerance for thermostat

   double temperatureTolerance                     = 3.0;
   double initialTemperature                       = 300.0;

   int energyMinimize                              = 1;
   int applyAssertion                              = 0;
   int randomNumberSeed                            = 0;

   // tolerance for energy conservation test

   double energyTolerance                          = 0.05;

   double equilibrationTime                        = 1000.0;
   double equilibrationTimeBetweenReportsRatio     = 0.1;

   double simulationTime                           = 1000.0;
   double simulationTimeBetweenReportsRatio        = 0.01;

   int allTypes                                    = State::Positions | State::Velocities | State::Forces | State::Energy;

// ---------------------------------------------------------------------------------------

    MapStringVectorOfVectors supplementary;
    MapStringVec3 tinkerForces;
    MapStringDouble tinkerEnergies;

    Context* context = createContext( parameterFileName, forceMap, useOpenMMUnits, inputArgumentMap, supplementary,
                                      tinkerForces, tinkerEnergies, log );

    setIntFromMap( inputArgumentMap, "applyAsser",          applyAssertion );


    setDoubleFromMap(     inputArgumentMap, "equilibrationTime",          equilibrationTime  );
    setDoubleFromMap(     inputArgumentMap, "simulationTime",             simulationTime     );
    const double totalTime  = equilibrationTime + simulationTime;

    std::string intermediateStateFileName = "NA";
    setStringFromMap( inputArgumentMap, "intermediateStateFileName",  intermediateStateFileName );

    FILE* intermediateStateFile = NULL;
    if( intermediateStateFileName != "NA" ){
        intermediateStateFile = openFile( intermediateStateFileName, "w", log );
        writeIntermediateStateFile( *context, intermediateStateFile, log );
    }
 
    System& system                     = context->getSystem();
    int numberOfAtoms                  = system.getNumParticles();
 
    std::vector<Vec3> velocities; 
    velocities.resize( numberOfAtoms );
    setIntFromMap(     inputArgumentMap, "randomNumberSeed",     randomNumberSeed );
    setDoubleFromMap(  inputArgumentMap, "temperature",          initialTemperature );
    setVelocitiesBasedOnTemperature( system, randomNumberSeed, velocities, initialTemperature, log );
    context->setVelocities(velocities);

    // energy minimize

    setIntFromMap( inputArgumentMap, "energyMinimize",          energyMinimize );
    if( log ){ 
        if( energyMinimize ){
            (void) fprintf( log, "Applying energy minimization before equilibration.\n" ); 
        } else {
            (void) fprintf( log, "Not applying energy minimization before equilibration.\n" ); 
        }
        (void) fflush( log );
    }
    if( energyMinimize ){
        State preState                            = context->getState( State::Energy );
        LocalEnergyMinimizer::minimize(*context);
        State postState                           = context->getState( State::Energy );
        if( log ){ 
            (void) fprintf( log, "Energy pre/post energies [%15.7e %15.7e] [%15.7e %15.7e].\n",
                            preState.getKineticEnergy(), preState.getPotentialEnergy(),
                            postState.getKineticEnergy(), postState.getPotentialEnergy() );
            (void) fflush( log );
        }
        if( intermediateStateFile ){
            writeIntermediateStateFile( *context, intermediateStateFile, log );
        }
    }

// ---------------------------------------------------------------------------------------

    int returnStatus                    = 0;
    double currentTime                  = 0.0;
 
    // set velocities based on temperature
 
    // get integrator 
 
    Integrator& integrator                                  = context->getIntegrator();
    std::string integratorName                              = getIntegratorName( &integrator );
    int  isVariableIntegrator                               = 0;
    int  isVerletIntegrator                                 = 0;
    VariableLangevinIntegrator*  variableLangevinIntegrator = NULL;
    VariableVerletIntegrator*  variableVerletIntegrator     = NULL;
    int stateFieldsToRetreive                               = State::Energy;

    if( integratorName == "VariableLangevinIntegrator" ){
        variableLangevinIntegrator = dynamic_cast<VariableLangevinIntegrator*>(&integrator);
        isVariableIntegrator       = 1;
    } else if( integratorName == "VariableVerletIntegrator" ){
        variableVerletIntegrator   = dynamic_cast<VariableVerletIntegrator*>(&integrator);
        isVariableIntegrator       = 2;
    } else if( integratorName == "VerletIntegrator" ){
        isVerletIntegrator         = 1;
        //stateFieldsToRetreive     |= State::Velocities;
    }

    if( log ){
        printIntegratorInfo( integrator, log );
    }
 
    // create/initialize arrays used to track energies
 
    std::vector<double> timeArray;
    std::vector<double> kineticEnergyArray;
    std::vector<double> potentialEnergyArray;
    std::vector<double> totalEnergyArray;
 
    State state                            = context->getState( State::Energy );
    double kineticEnergy                   = state.getKineticEnergy();
    double potentialEnergy                 = state.getPotentialEnergy();
    double totalEnergy                     = kineticEnergy + potentialEnergy;
 
    // log
 
    if( log ){
        (void) fprintf( log, "Initial energies: E=%15.7e [%15.7e %15.7e]\n",
                        (kineticEnergy + potentialEnergy), kineticEnergy, potentialEnergy );
        (void) fflush( log );
    }
 
    /* -------------------------------------------------------------------------------------------------------------- */
 
    setDoubleFromMap(    inputArgumentMap, "equilibrationTimeBetweenReportsRatio",      equilibrationTimeBetweenReportsRatio );
    double equilibrationTimeBetweenReports   = equilibrationTime*equilibrationTimeBetweenReportsRatio;
 
    setDoubleFromMap(    inputArgumentMap, "simulationTimeBetweenReportsRatio",         simulationTimeBetweenReportsRatio );
    double simulationTimeBetweenReports   = simulationTime*simulationTimeBetweenReportsRatio;

    // if equilibrationTimeBetweenReports || simulationTimeBetweenReports <= 0, take one step at a time
    if( equilibrationTimeBetweenReports <= 0.0 || simulationTimeBetweenReports <= 0.0 ){
         isVariableIntegrator  = -1;
    }
 
    if( log ){
        (void) fprintf( log, "Equilibration/simulation times [%12.3f %12.3f] timeBetweenReports [ %12.3e %12.3e] ratios [%12.4f %12.4f] variableIntegrator=%d VerletIntegrator=%d\n", 
                        equilibrationTime, simulationTime,
                        equilibrationTimeBetweenReports, simulationTimeBetweenReports,
                        equilibrationTimeBetweenReportsRatio, simulationTimeBetweenReportsRatio, isVariableIntegrator, isVerletIntegrator );
        (void) fflush( log );
    }
 
    // set dof
 
    double degreesOfFreedom  = static_cast<double>(3*numberOfAtoms - system.getNumConstraints() - 3 );

    // main simulation loop
 
    double timeBetweenReports;
    int    stepsBetweenReports;
    bool   equilibrating;
    if( equilibrationTime > 0.0 ){
        timeBetweenReports  = equilibrationTimeBetweenReports;
        stepsBetweenReports = isVariableIntegrator > 0 ? 1 : static_cast<int>(equilibrationTimeBetweenReports/integrator.getStepSize() + 1.0e-04);
        equilibrating       = true;
    } else {
        timeBetweenReports  = simulationTimeBetweenReports;
        stepsBetweenReports = isVariableIntegrator > 0 ? 1 : static_cast<int>(simulationTimeBetweenReports/integrator.getStepSize() + 1.0e-05);
        equilibrating       = false;
        if( isVerletIntegrator && stepsBetweenReports > 1 )stepsBetweenReports -= 1;
    }
    if( stepsBetweenReports < 1 )stepsBetweenReports = 1;

    double simulationStartTime = 0.0;
    double totalWallClockTime  = 0.0;
    double energyDrift         = 0.0;
    int totalShakeViolations   = 0;
    while( currentTime < totalTime ){
 
        double startTime        = getTimeOfDay();

        if( isVariableIntegrator <= 0 ){ 
            integrator.step( stepsBetweenReports );
        } else if( isVariableIntegrator == 1 ){ 
            variableLangevinIntegrator->stepTo( currentTime + timeBetweenReports);
        } else if( isVariableIntegrator == 2 ){ 
            variableVerletIntegrator->stepTo( currentTime + timeBetweenReports);
        }

        double elapsedTime                     = getTimeOfDay() - startTime;
        totalWallClockTime                    += elapsedTime;
  
        State state                            = context->getState( stateFieldsToRetreive );
        currentTime                            = state.getTime();
        double kineticEnergy                   = state.getKineticEnergy();
        double potentialEnergy                 = state.getPotentialEnergy();
        double totalEnergy                     = kineticEnergy + potentialEnergy;
  
        if( intermediateStateFile ){
            writeIntermediateStateFile( *context, intermediateStateFile, log );
        }
  
        if( equilibrating && currentTime >= equilibrationTime ){ 

            equilibrating       = false;
            simulationStartTime = state.getTime();
            timeBetweenReports  = simulationTimeBetweenReports;
            stepsBetweenReports = isVariableIntegrator != 0 ? 1 : static_cast<int>(simulationTimeBetweenReports/integrator.getStepSize() + 1.0e-04);
            if( stepsBetweenReports < 1 )stepsBetweenReports  = 1;
            if( isVerletIntegrator && stepsBetweenReports > 1 )stepsBetweenReports -= 1;

        } else if( !equilibrating ){ 
   
            if( isVerletIntegrator ){
                getVerletKineticEnergy( *context, currentTime, potentialEnergy, kineticEnergy, log );
            }

            // record energies
      
            timeArray.push_back( currentTime - simulationStartTime );
            kineticEnergyArray.push_back( kineticEnergy );
            potentialEnergyArray.push_back( potentialEnergy );
            totalEnergyArray.push_back( totalEnergy );
            energyDrift = getEnergyDrift( totalEnergyArray, kineticEnergyArray, degreesOfFreedom, (currentTime-simulationStartTime), NULL );
        } 

        // diagnostics & check for nans
  
        if( log ){
            double nsPerDay = 86.4*currentTime/totalWallClockTime;
            (void) fprintf( log, "%12.3f KE=%15.7e PE=%15.7e E=%15.7e wallClock=%12.3e %12.3e %12.3f ns/day", currentTime, kineticEnergy, potentialEnergy, totalEnergy,
                            elapsedTime, totalWallClockTime, nsPerDay );
            if( equilibrating ){
                (void) fprintf( log, " equilibrating" );
            } else if( isVerletIntegrator ){
                (void) fprintf( log, " drift=%12.3e", energyDrift );
            }
        }
  
        if( isNanOrInfinity( totalEnergy ) ){
            char buffer[1024];
            (void) sprintf( buffer, "%s nans detected at time %12.3f -- aborting.\n", methodName.c_str(), currentTime );
            throwException(__FILE__, __LINE__, buffer );
            exit(-1);
        }
 
        // check constraints

        if( system.getNumConstraints() > 0 ){
            double maxViolation;
            int maxViolationIndex;
            int violations        = checkConstraints( *context, integrator.getConstraintTolerance(), maxViolation, maxViolationIndex, log );
            totalShakeViolations += violations;
            if( violations && log ){
                (void) fprintf( log, " Shake violations %d max=%12.3f at index=%d", violations, maxViolation, maxViolationIndex );
            } 
        }
        if( log ){
            (void) fprintf( log, "\n" );
            (void) fflush( log );
        }
    }
 
    state                                  = context->getState( State::Energy );
    double simulationEndTime               = state.getTime();
    if( isVerletIntegrator ){
        getVerletKineticEnergy( *context, simulationEndTime, potentialEnergy, kineticEnergy, log );
    } else {
        kineticEnergy                      = state.getKineticEnergy();
        potentialEnergy                    = state.getPotentialEnergy();
    }
    totalEnergy                            = kineticEnergy + potentialEnergy;
 
    // log times and energies
 
    if( log ){
       double nsPerDay = 86.4*totalTime/totalWallClockTime;
       (void) fprintf( log, "Final Simulation: %12.3f  E=%15.7e [%15.7e %15.7e]  total wall time=%12.3e ns/day=%.3e Shake violations=%d\n",
                       currentTime, (kineticEnergy + potentialEnergy), kineticEnergy, potentialEnergy,
                       totalWallClockTime, nsPerDay, totalShakeViolations );
       (void) fprintf( log, "\n%8u Energies\n", static_cast<unsigned int>(kineticEnergyArray.size()) );
       for( unsigned int ii = 0; ii < kineticEnergyArray.size(); ii++ ){
           (void) fprintf( log, "%15.7e   %15.7e %15.7e  %15.7e    Energies\n",
            timeArray[ii], kineticEnergyArray[ii], potentialEnergyArray[ii], totalEnergyArray[ii] );
       }
       (void) fflush( log );
    }
 
    double conversionFactor  = degreesOfFreedom*0.5*BOLTZ;
           conversionFactor  = 1.0/conversionFactor;
 
    // if Langevin or Brownian integrator, then check that temperature constant
    // else (Verlet integrator) check that energy drift is acceptable
 
    if( (integratorName == "LangevinIntegrator"           ||
         integratorName == "VariableLangevinIntegrator"   ||
         integratorName == "BrownianIntegrator" ) && numberOfAtoms > 0 ){
 
       // check that temperature constant
       // convert KE to temperature
 
       std::vector<double> temperature;
       for( std::vector<double>::const_iterator ii = kineticEnergyArray.begin(); ii != kineticEnergyArray.end(); ii++ ){
          temperature.push_back( (*ii)*conversionFactor );
       }
 
       // get temperature stats
 
       std::vector<double> temperatureStatistics;
       getStatistics( temperature, temperatureStatistics );
       double initialTemperature = 300.0;
 
       if( integratorName == "LangevinIntegrator" ){
           LangevinIntegrator* langevinIntegrator          = dynamic_cast<LangevinIntegrator*>(&integrator);
           initialTemperature                              = langevinIntegrator->getTemperature();
       } else if( integratorName == "VariableLangevinIntegrator" ){
           VariableLangevinIntegrator* langevinIntegrator  = dynamic_cast<VariableLangevinIntegrator*>(&integrator);
           initialTemperature                              = langevinIntegrator->getTemperature();
       }
 
       if( log ){
          (void) fprintf( log, "Simulation temperature results: mean=%15.7e stddev=%15.7e   min=%15.7e   %d max=%15.7e %d\n",
                          temperatureStatistics[0], temperatureStatistics[1], temperatureStatistics[2],
                          (int) (temperatureStatistics[3] + 0.001), temperatureStatistics[4],
                          (int) (temperatureStatistics[5] + 0.001) );
 
       }
 
       // check that <temperature> is within tolerance
 
       if( applyAssertion ){
           ASSERT_EQUAL_TOL( temperatureStatistics[0], initialTemperature, temperatureTolerance );
       }
 
    } else {
 
        double stddevE = getEnergyDrift( totalEnergyArray, kineticEnergyArray, degreesOfFreedom, (simulationEndTime-simulationStartTime), log );
 
        // check that energy fluctuation is within tolerance
 
        if( applyAssertion ){
            ASSERT_EQUAL_TOL( stddevE, 0.0, energyTolerance );
        }
 
    }
 
    return;

}

// ---------------------------------------------------------------------------------------

int runTestsUsingAmoebaTinkerParameterFile( MapStringString& argumentMap ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName       = "runTestsUsingAmoebaTinkerParameterFile";
   MapStringString inputArgumentMap;
   std::string openmmPluginDirectory         = ".";

   FILE* log                                 = NULL;
   FILE* summaryFile                         = NULL;

// ---------------------------------------------------------------------------------------

/* command-line args

Int    "cudaDevice"                              cudaDevice
Int    "applyAssert"                             apply assertions 

                                                 // integrator args

String "integrator"                              integratorName
Double "timeStep"                                timeStep
Double "friction"                                friction
Double "temperature"                             temperature
Double "shakeTolerance"                          shakeTolerance
Double "errorTolerance"                          errorTolerance
Int    "randomNumberSeed"                        randomNumberSeed

                                                 // save states or read
String "states"                                  statesFileName

Double "tolerance"                               general tolerance
Double "energyForceDelta"                        delta energy/force consistency check          
Double "energyForceTolerance"                    tolerance nergy/force consistency check

Int    "energyMinimize"                          energyMinimize -- minimize structures before runs
Double "equilibrationTime"                       equilibrationTime 
Double "simulationTime"                          simulationTime    
String "intermediateStateFileName"               intermediateStateFileName -- name of file to write intermediate structures
Double "equilibrationTimeBetweenReportsRatio"    equilibrationTimeBetweenReportsRatio -- if for example set to 0.1, then output at t=0.1*equilibrationTime 
                                                                                                                                     0.2*equilibrationTime
                                                                                                                                     ...
                                                                                                                                     0.9*equilibrationTime
                                                                                                                                     1.0*equilibrationTime
Double "simulationTimeBetweenReportsRatio"       simulationTimeBetweenReportsRatio (see equilibrationTimeBetweenReportsRatio)

*/

// ---------------------------------------------------------------------------------------


    std::string parameterFileName            = "1UBQ.prm";
    MapStringInt forceMap;
    initializeForceMap( forceMap, 0 );

    int logFileNameIndex                     = 0;
    std::string logFileName;

    int summaryFileNameIndex                 = 0;
    std::string summaryFileName;

    int specifiedOpenmmPluginDirectory       = 0;
    int useOpenMMUnits                       = 1;
    int logControl                           = 0;

    int checkForces                          = 1;
    int checkEnergyForceConsistency          = 0;
    int checkEnergyForceByFiniteDifference   = 0;
    int checkEnergyConservation              = 0;
    int checkIntermediateStates              = 0;

    // parse arguments

    for( MapStringStringCI ii = argumentMap.begin(); ii != argumentMap.end(); ii++ ){
        std::string key                                = ii->first;
        std::string value                              = ii->second;
        if( key == "parameterFileName" ){
            parameterFileName                          = value;
        } else if( key == "logFileName" ){
            logFileNameIndex                           = 1;
            logFileName                                = value;
        } else if( key == "summaryFileName" ){
            summaryFileNameIndex                       = 1;
            summaryFileName                            = value;
        } else if( key == "openmmPluginDirectory" ){
            specifiedOpenmmPluginDirectory             = 1;
            openmmPluginDirectory                      = value;
        } else if( key == "useOpenMMUnits" ){
            useOpenMMUnits                             = atoi( value.c_str() );
        } else if( key == "checkEnergyForceConsistency" ){
            checkEnergyForceConsistency                = atoi( value.c_str() );
        } else if( key == "checkEnergyForceByFiniteDifference" ){
            checkEnergyForceByFiniteDifference         = atoi( value.c_str() );
        } else if( key == "checkEnergyConservation" ){
            checkEnergyConservation                    = atoi( value.c_str() );
        } else if( key == "checkIntermediateStates" ){
            checkIntermediateStates                    = atoi( value.c_str() );
        } else if( key == "log" ){
            logControl                                 = atoi( value.c_str() );
        } else if( key == ALL_FORCES ){
            initializeForceMap( forceMap, atoi( value.c_str() ) );
        } else if( key == AMOEBA_HARMONIC_BOND_FORCE              ||
                   key == AMOEBA_HARMONIC_ANGLE_FORCE             ||
                   key == AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE    ||
                   key == AMOEBA_TORSION_FORCE                    ||
                   key == AMOEBA_PI_TORSION_FORCE                 ||
                   key == AMOEBA_STRETCH_BEND_FORCE               ||
                   key == AMOEBA_OUT_OF_PLANE_BEND_FORCE          ||
                   key == AMOEBA_TORSION_TORSION_FORCE            ||
                   key == AMOEBA_MULTIPOLE_FORCE                  ||
                   key == AMOEBA_GK_FORCE                         ||
                   key == AMOEBA_VDW_FORCE                        ||
                   key == AMOEBA_WCA_DISPERSION_FORCE             ){
            forceMap[key]                              = atoi( value.c_str() );
        } else {
            inputArgumentMap[key]                      = value;
        }
    }

    // open log file

    if( logControl ){
        std::string mode = logControl == 1 ? "w" : "a";
        if( logFileNameIndex > -1 ){
            log = openFile( logFileName, mode, NULL );
        } else {
            log = stderr;
        }
    }
   
    // summary file

    if( summaryFileNameIndex > 0 ){
        std::string mode = "a";
        summaryFile = openFile( summaryFileName, mode, log );
    }

    // log info

    if( log ){
        (void) fprintf( log, "Input arguments:\n" );
        for( MapStringStringCI ii = argumentMap.begin(); ii != argumentMap.end(); ii++ ){
            std::string key   = ii->first;
            std::string value = ii->second;
            (void) fprintf( log, "      %30s %40s\n", key.c_str(), value.c_str() );
        }
        (void) fprintf( log, "\nParameter file=<%s>\n", parameterFileName.c_str() );

        (void) fprintf( log, "\nArgument map: %u\n", static_cast<unsigned int>(inputArgumentMap.size()) );
        for( MapStringStringCI ii = inputArgumentMap.begin(); ii != inputArgumentMap.end(); ii++ ){
            (void) fprintf( log, "   %s=%s\n", (*ii).first.c_str(), (*ii).second.c_str() );
        }
        (void) fprintf( log, "\nForce map: %u\n", static_cast<unsigned int>(forceMap.size()) );
        for( MapStringIntCI ii = forceMap.begin(); ii != forceMap.end(); ii++ ){
            (void) fprintf( log, "   %s=%d\n", (*ii).first.c_str(), (*ii).second );
        }
        (void) fflush( log );
    }

    // load plugins

    if( specifiedOpenmmPluginDirectory ){
        Platform::loadPluginsFromDirectory( openmmPluginDirectory );
    }

    if( checkEnergyForceConsistency ){
        // args:

        //     applyAssertion
        //     energyForceDelta 
        //     energyForceTolerance
        testEnergyForcesConsistent( parameterFileName, forceMap, useOpenMMUnits, 
                                    inputArgumentMap, log, summaryFile );

    } else if( checkEnergyForceByFiniteDifference ){
        // args:

        //     applyAssertion
        //     energyForceDelta 
        //     energyForceTolerance
        testEnergyForceByFiniteDifference( parameterFileName, forceMap, useOpenMMUnits, 
                                           inputArgumentMap, log, summaryFile );

    } else if( checkEnergyConservation ){
        // args:

        testEnergyConservation( parameterFileName, forceMap, useOpenMMUnits, 
                                inputArgumentMap, log, summaryFile );

    } else if( checkIntermediateStates ){
        // args:

        checkIntermediateStatesUsingAmoebaTinkerParameterFile( parameterFileName, forceMap, useOpenMMUnits, 
                                                               inputArgumentMap, summaryFile, log );

    } else {
        // args:
        //     tolerance
        testUsingAmoebaTinkerParameterFile( parameterFileName, forceMap,
                                            useOpenMMUnits, inputArgumentMap, summaryFile, log );
    }
    if( log ){
        (void) fprintf( log, "\n%s done\n", methodName.c_str() ); (void) fflush( log );
    }
    if( summaryFile ){
        (void) fclose( summaryFile );
    }

    return 0;
}

// ---------------------------------------------------------------------------------------

void appendInputArgumentsToArgumentMap( int numberOfArguments, char* argv[], MapStringString& argumentMap ){

// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------

    for( int ii = 1; ii < numberOfArguments; ii += 2 ){
        char* key        = argv[ii];
        if( *key == '-' )key++;
        argumentMap[key] = (ii+1) < numberOfArguments ? argv[ii+1] : "NA";
    }

    return;
}
