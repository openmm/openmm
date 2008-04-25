
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

#ifndef __GbsaParameters_H__
#define __GbsaParameters_H__

// Gromacs data structs

#ifndef M_PI
#define M_PI 3.1415926535
#endif

#include "typedefs.h"

// ---------------------------------------------------------------------------------------

// buffer size of atom id strings
#define MAX_BUFFER_SZ 128

#include "GbsaAtomParameter.h"

struct _gbsaStretchBondStruct {

   int atomI;
   int atomJ;
   Real bondLength;
   struct _gbsaStretchBondStruct* nextBond;
 
};
typedef struct _gbsaStretchBondStruct gbsaStretchBond;

struct _gbsaAngleBondStruct {

   gbsaStretchBond *stretchBondI;
   gbsaStretchBond *stretchBondJ;
   int pivotAtomIndex;
   int volumeIndex;
   int orderedIndices[3];
   Real harmonicAngleWidth;
   struct _gbsaAngleBondStruct* nextBond;
 
};
typedef struct _gbsaAngleBondStruct gbsaAngleBond;

#define MAX_BOND_EXCLUSIONS 20
struct _gbsaBondStruct {

   int stretchBondCount;
   gbsaStretchBond* stretchBonds;

   int angleBondCount;
   gbsaAngleBond*   angleBonds;

   int numberOfExclusions;
   int exclusions[MAX_BOND_EXCLUSIONS];

   int minExclusionIndex;
   int maxExclusionIndex;

   int startExclusionIndex;
   int stopExclusionIndex;

};
typedef struct _gbsaBondStruct gbsaBonds;

/*
class GbsaAtomInfo {

   private:

      static const int MaxBondExclusions = 20;

      int _stretchBondCount;
      gbsaStretchBond* _stretchBonds;
   
      int _angleBondCount;
      gbsaAngleBond* _angleBonds;
   
      int _numberOfExclusions;
      int _exclusions[MaxBondExclusions];
   
      int _minExclusionIndex;
      int _maxExclusionIndex;
   
      int _startExclusionIndex;
      int _stopExclusionIndex;

   public:

      // constructor/destructor

      GbsaAtomInfo( );
      ~GbsaAtomInfo( );
   
      // accessors

      int getStretchBondCount( void               ){ return _stretchBondCount;   };
      gbsaStretchBond* getGbsaStretchBond( void   ){ return _stretchBonds;       };

      int getAngleBondCount( void                 ){ return _angleBondCount;     };
      gbsaAngleBond* getGbsaAngleBond( void       ){ return _angleBonds;         };

      int  getNumberOfExclusions( void            ){ return _numberOfExclusions;  };
      int* get_exclusions( void                   ){ return _exclusions;          };

      int  getMinExclusionIndex( void             ){ return _minExclusionIndex;   };
      int  getMaxExclusionIndex( void             ){ return _maxExclusionIndex;   };

      int  getStartExclusionIndex( void           ){ return _startExclusionIndex; };
      int  getStopExclusionIndex( void            ){ return _stopExclusionIndex;  };
};
*/

class GbsaParameters {

   private:

      // GBSA constants & parameters; parameters from '97 paper
   
      Real _solventDielectric;
      Real _innerDielectric;
      Real _kcalA_To_kJNm;
      Real _phi;
      Real _P1;
      Real _P2;
      Real _P3;
      Real _P4;
      Real _P4_2;
      Real _P5;
      Real _electricConstant;
      Real _probeRadius;
   
      Real _preFactor;
      Real _P5Inverse;
      Real _piP5;
      Real _pi4Asolv;

      // ---------------------------------------------------------------------------------------

      int _numberOfAtoms;
      FILE* _log;

      bool _freeArrays;

      Real* _vdwRadii;
      Real* _volume;
      Real* _gPolFixed;
      Real* _bornRadii;

      // bond and exclusion info

      gbsaBonds** _gbsaBondsArray;

      // work array

      int* _exclusionWorkArray;

      // currently not exposed 

      int setGbsaBondsArray( gbsaBonds** gbsaBondsArray ){ _gbsaBondsArray = gbsaBondsArray; return 0; };

      /**---------------------------------------------------------------------------------------
      
         Reset prefactor (Simbios) 
      
         called when _electricConstant, _innerDielectric, or _solventDielectric are modified
      
         --------------------------------------------------------------------------------------- */
      
      void _resetPreFactor( void );

   public:

       // constructor/destructor

       GbsaParameters( int numberOfAtoms, FILE* log );
       ~GbsaParameters( );

      // override of new/delete

      static void* operator new( size_t size );
      static void  operator delete( void *p );

      static void* operator new[]( size_t size );
      static void  operator delete[]( void *p );

      // accessors

      int getNumberOfAtoms( void ){ return _numberOfAtoms; };

      FILE* getLog( void ){ return _log; };
      int setLog( FILE* log ){ _log = log; return 0; };

      // array accessors

      Real* getVdwRadii(  void ){ return _vdwRadii;  };
      Real* getVolume(    void ){ return _volume;    };
      Real* getGPolFixed( void ){ return _gPolFixed; };
      Real* getBornRadii( void ){ return _bornRadii; };

      int setVdwRadii(  Real* vdwRadii   ){ _vdwRadii   = vdwRadii;  return 0; };
      int setVolume(    Real* volume     ){ _volume     = volume;    return 0; };
      int setGPolFixed( Real* gPolFixed  ){ _gPolFixed  = gPolFixed; return 0; };
      int setBornRadii( Real* bornRadii  ){ _bornRadii  = bornRadii; return 0; };

      // flag signalling whether arrays should be deleted
      // if destructor is called

      int setFreeArrays( bool freeArrays  ){ _freeArrays = freeArrays; return 0; };

      // work array
      int* getExclusionWorkArray( void );

      // bond info

      gbsaBonds** getGbsaBondsArray( void ){ return _gbsaBondsArray; };

      // get accessor for constants

      Real getSolventDielectric(  void   ){ return _solventDielectric; };
      Real getInnerDielectric(    void   ){ return _innerDielectric;   };
      Real getKcalA_To_kJNm(      void   ){ return _kcalA_To_kJNm;     };
      Real getPhi(                void   ){ return _phi;               };
      Real getP1(                 void   ){ return _P1;                };
      Real getP2(                 void   ){ return _P2;                };
      Real getP3(                 void   ){ return _P3;                };
      Real getP4(                 void   ){ return _P4;                };
      Real getP4_2(               void   ){ return _P4_2;              };
      Real getP5(                 void   ){ return _P5;                };
      Real getElectricConstant(   void   ){ return _electricConstant;  };
      Real getProbeRadius(        void   ){ return _probeRadius;  };
   
      Real getPreFactor( void ){ return _preFactor; };
      Real getP5Inverse( void ){ return _P5Inverse; };
      Real getPiP5(      void ){ return _piP5;      };
      Real getPi4Asolv(  void ){ return _pi4Asolv;  };

      /**---------------------------------------------------------------------------------------
      
         Set solvent dielectric (Simbios) 
      
         @param solventDielectric         solvent dielectric
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      int setSolventDielectric( Real solventDielectric );
      
      /**---------------------------------------------------------------------------------------
      
         Set inner dielectric (Simbios) 
      
         @param innerDielectric    inner dielectric

         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int setInnerDielectric( Real innerDielectric );

      /**---------------------------------------------------------------------------------------
      
         Set electric constant (Simbios) 
      
         @param electricConstant   electric constant
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      int setElectricConstant( Real electricConstant );

      /**---------------------------------------------------------------------------------------
      
         Initialize Gbsa Constants (Simbios) 
      
         --------------------------------------------------------------------------------------- */

      void initializeConstants( void );

      /**---------------------------------------------------------------------------------------
      
         Initialize Gbsa Parameters (Simbios) 
      
         @param top                 GMX topology data struct 
      
         @return 0

         --------------------------------------------------------------------------------------- */

      int initializeParameters( const t_topology* top );

      /**---------------------------------------------------------------------------------------
      
         Set Gbsa Exclusions (Simbios) 
      
         @param top                 GMX topology data struct -- used to get harmonic angle
         @param printOn             print diagnostics 
         @param log                 file descriptor 
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int setGbsaExclusions( const t_topology* top, bool printOn, FILE* log );
      
      /**---------------------------------------------------------------------------------------

         Allocate memory for bond array (Simbios) 

         @param maxAtoms            max number of atoms

         array entries are intialized to zero

         free array[0] and then array when done

         @return ptr to allocated array or NULL if out of memory

         --------------------------------------------------------------------------------------- */

      gbsaBonds** allocateBondsArray( int maxAtoms );
      
      /**---------------------------------------------------------------------------------------
      
         Deallocate memory for bond array (Simbios) 
      
         @param maxAtoms            max number of atoms
         @param gbsaBondsArray      array to be freed
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      int freeBondArray( int maxAtoms, gbsaBonds** gbsaBondsArray );
      
      /**---------------------------------------------------------------------------------------
      
         Print bond array (Simbios) 
      
         @param numberOfAtoms       number of atoms
         @param top                 GMX topology data struct -- used to get harmonic angle
         @param gbsaBondsArray      bond array
         @param title               title string (optional)
         @param log                 file descriptor 
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int printBondsArray( int numberOfAtoms, const t_topology* top,
                           gbsaBonds** gbsaBondsArray,
                           const char* title, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Print bond (Simbios) 
      
         @param numberOfAtoms         number of atoms
         @param top                 GMX topology data struct -- used to get harmonic angle
         @param atomIndex           atom index
         @param gbsaBond             bond data struct
         @param log                 file descriptor 
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int printBond( int numberOfAtoms, const t_topology* top, int atomIndex, 
                     const gbsaBonds* gbsaBond, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Print parameter info (Simbios) 
      
         @param numberOfAtoms         number of atoms
         @param top                   GMX topology data struct -- used to get harmonic angle
         @param parameterInfoFileName parameterInfo FileName
         @param log                   file descriptor 
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int printParameterInfo( int numberOfAtoms, const t_topology* top,
                              const char* parameterInfoFileName, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Add angle bonds using stretch bond lists (Simbios) 
         Also load in angle width
      
         @param maxAtoms            max number of atoms
         @param gbsaBondsArray      array of gbsaBonds
         @param top                 GMX topology data struct -- used to get harmonic angle
         @param idefArrayIndex      parameter index for GMX iparams data struct
         @param log                 if set, then print error messages to log file
      
         @return 0 if no errors or
      
         --------------------------------------------------------------------------------------- */
      
      int addAngleBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                         const t_topology* top, int idefArrayIndex, FILE* log  );
                  
      /**---------------------------------------------------------------------------------------
      
         Add stretch bonds to angle data struct (Simbios) 
      
         @param stretchBondI        stretch bond I
         @param stretchBondJ        stretch bond J
         @param pivotAtomIndex      index of shared atom
         @param gbsaBonds           gbsaBond
         @param saveIndex           index to save????
         @param log                 if set, then print error messages to log file
      
         @return allocated gbsaAngleBond 
      
         --------------------------------------------------------------------------------------- */
      
      gbsaAngleBond* addStretchBondsToAngleBondList( gbsaStretchBond* stretchBondI, 
                                                     gbsaStretchBond* stretchBondJ,
                                                     int pivotAtomIndex,
                                                     gbsaBonds** gbsaBonds, int saveIndex, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Find angle bond with atom indices that match input angle indices (Simbios) 
      
         @param gbsaBond             gbsaBond
         @param atomI               index of atomI
         @param atomJ               index of atomJ
         @param atomK               index of atomK
      
         @return allocated gbsaAngleBond 
      
         --------------------------------------------------------------------------------------- */
      
      gbsaAngleBond* findAngleBond( const gbsaBonds* gbsaBond, int atomI, int atomJ, int atomK ) const;

      /**---------------------------------------------------------------------------------------
      
         Find angle bond with atom indices that match input ordered angle indices (Simbios) 
      
         gbsaBond             gbsaBond
         orderedAtomIndices  array of ordered indices
      
         return gbsaAngleBond if found
      
         --------------------------------------------------------------------------------------- */
      
      gbsaAngleBond* findAngleBondGivenOrderedArray( const gbsaBonds* gbsaBond, int* orderedAtomIndices ) const;

      /**---------------------------------------------------------------------------------------
      
         Sort atom indices
      
         @param atomI               index of atomI
         @param atomJ               index of atomJ
         @param atomK               index of atomK
         @param atomL               index of atomL
         @param orderedIndices      output array of ordered indices assumed to be of size 3
      
         @return count (should be 3, if 4, then all the atom indices were distinct)
      
         --------------------------------------------------------------------------------------- */
      
      int sortAtomAngleIndices( int atomI, int atomJ, int atomK, int atomL,
                                int* orderedIndices ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get bond counts (Simbios) 
      
         @param gbsaBond             gbsaBonds to check
         @param stretchCount        stretch count on return
         @param angleCount          angle count on return
      
         @return 0

         --------------------------------------------------------------------------------------- */
      
      int getBondCounts( const gbsaBonds* gbsaBond, int* stretchCount, int* angleCount ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get stretch bond count (Simbios) 
      
         @param gbsaBond              gbsaBonds to check
      
         @return stretch bond count
      
         --------------------------------------------------------------------------------------- */
      
      int getStretchBondCount( const gbsaBonds* gbsaBond ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get angle bond count (Simbios) 
      
         @param gbsaBond              gbsaBonds to check
      
         @return angle bond count
      
         --------------------------------------------------------------------------------------- */
      
      int getAngleBondCount( const gbsaBonds* gbsaBond ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Add stretch bonds (Simbios) 
      
         @param maxAtoms            max number of atoms
         @param gbsaBondsArray       array of gbsaBonds data structs
         @param top                 Gromacs t_topolgy struct
         @param idefArrayIndex      index to bond parameters (F_BONDS, F_SHAKE, ... )
         @param offset              number of entries for each bond block
         @param atomIndexOffset     offset into block for atom indices
         @param log                 if set, then print error messages to log file
      
         @return 0 if no errors or
         return x, where x is the number of errors encountered
      
         --------------------------------------------------------------------------------------- */
      
      int addStretchBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                           const t_topology* top, int idefArrayIndex, int offset,
                           int atomIndexOffset, FILE* log  );
      
      /**---------------------------------------------------------------------------------------
      
         Add SETTLE stretch bonds (Simbios) 
      
         @param maxAtoms            max number of atoms
         @param gbsaBondsArray       array of gbsaBonds data structs
         @param top                 Gromacs t_topolgy struct
         @param idefArrayIndex      index to bond parameters (F_BONDS, F_SHAKE, ... )
         @param offset              number of entries for each bond block
         @param atomIndexOffset     offset into block for atom indices
         @param log                 if set, then print error messages to log file
      
         @return 0 if no errors or
         return x, where x is the number of errors encountered
      
         --------------------------------------------------------------------------------------- */
      
      int addSettleStretchBonds( int maxAtoms, gbsaBonds** gbsaBondsArray,
                                 const t_topology* top, int idefArrayIndex, int offset,
                                 int atomIndexOffset, FILE* log  );
      
      /**---------------------------------------------------------------------------------------
      
         Check if atom is in bond list (as atom j) (Simbios)
      
         @param atomIndex           index of atom to be searched for
         @param bond                gbsaBond data struct
         @param log                 if set, then print error messages to log file
      
         @return 0 if atom is in StretchBondList; 1 otherwise
      
         --------------------------------------------------------------------------------------- */
      
      int isAtomInStretchBondList( int atomIndex, const gbsaBonds* bond, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Add a stretch bonds (Simbios) 
      
         @param atomIndexI          index of atom I
         @param atomIndexJ          index of atom J
         @param bonds                gpu bond data struct
         @param bondLength          bond length
         @param log                 if set, then print error messages to log file
      
         @return gbsaStretchBond
      
         --------------------------------------------------------------------------------------- */
            
      gbsaStretchBond* addStretchBond( int atomIndexJ, int atomIndexI,
                                       gbsaBonds* bonds, Real bondLength, FILE* log );
      
      /**--------------------------------------------------------------------------------------- 
      
         Assign standard radii for GB/SA methods other than ACE;
         taken from Macromodel and OPLS-AA, except for hydrogens (Simbios)
      
         Logic based on logic in Tinker's ksolv.f
      
         Currently only works for standard amino acid atoms
         If invalid atom name is encountered, a message is printed to log file and the
         radius for that atom is set to 1.0f
      
         @param numberOfAtoms       number of atoms
         @param gbsaBondsArray       array of gbsaBonds
         @param atomNames           array of atom names from GMX top data struct
         @param radii               array to store Macromodel radii for each atom
         @param log                 if set, then print error messages to log file
      
         @return 0 always
      
         --------------------------------------------------------------------------------------- */
      
      int getMacroModelAtomicRadii( int numberOfAtoms, gbsaBonds** gbsaBondsArray,
                                    char*** atomNames, Real* radii, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Map Gmx atom name to Tinker atom number (Simbios)
      
         @param atomName            atom name (CA, HA, ...); upper and lower case should both work
         @param log                 if set, then print error messages to log file
      
         return Tinker atom number if atom name is valid; else return -1
      
         --------------------------------------------------------------------------------------- */
            
      int mapGmxAtomNameToTinkerAtomNumber( const char* atomName, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Get atomic volumes minus subvolumes that lie inside directly bonded atoms
         Eqs 6 & 7 in J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         Use Macromodel radii
      
         @param numberOfAtoms       number of atoms
         @param top                 GMX t_topology struct
         @param gbsaBondsArray      array of gbsaBonds
         @param vdwRadius           array of Macromodel radii for each atom
         @param volume              output array of volumes
         @param log                 if set, then print error messages to log file
      
         @return 0 always; atomic volumes in array 'volume'
         
         --------------------------------------------------------------------------------------- */
      
      int getMacroModelAtomicVolumes( int numberOfAtoms, const t_topology* top,
                                      gbsaBonds** gbsaBondsArray, Real* vdwRadius,
                                      Real* volume, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Calculate first two terms in Gpol equation (nearest and next-nearest neighbors) 
         Eqs 6 & 7 in J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         Use Macromodel radii
      
         @param numberOfAtoms       number of atoms
         @param top                 GMX t_topology struct
         @param gbsaBondsArray      array of gpu bonds
         @param vdwRadius           array of Macromodel radii for each atom
         @param volume              array of atomic volumes
         @param gPolFixed           output array of fixed GBSA values for GPol
         @param log                 if set, then print error messages to log file
      
         return 0 always; fixed value for G_pol in array 'gPolFixed'
      
         --------------------------------------------------------------------------------------- */
            
      int getFixedGBSA_GPol( int numberOfAtoms, const t_topology* top,
                             gbsaBonds** gbsaBondsArray,
                             Real* vdwRadius, Real* volume, Real* gPolFixed, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Get exclusions for specific atom (Simbios)
      
         @param numberOfAtoms       number of atoms
         @param gbsaBondsArray      array of gpu bonds
         @param atomI               atom index of atom for which exclusions are to be set
         @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
                                    value may be null on input in which space is allocated
         @param previousIndex       previousIndex -- if < 0, then iniitialize all entries to
                                    0
         @param log                 if set, then print error messages to log file
      
         @return 0 if no problems
      
         --------------------------------------------------------------------------------------- */
      
      int getExclusionsForAtom( int numberOfAtoms, const gbsaBonds** gbsaBondsArray,
                                 int atomI, int* exclusionWorkArray, int previousIndex, FILE* log ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Set exclusions for specific atom (Simbios)
      
         @param gbsaBonds           gpu bond
         @param setValue            set value
         @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
                                    value may be null on input in which space is allocated
         @param log                 if set, then print error messages to log file
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      int setExclusionValue( const gbsaBonds* gbsaBonds, int setValue, int* exclusionWorkArray, FILE* log ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Collect Gbsa Bond Indices for a particular atom (Simbios) 
      
         @param gbsaBond            gbsaBond 
         @param collectionBinMax    size of collectionBin array
         @param collectionBin       array of indices to be excluded
         @param excludeIndex        exclude index
         @param log                 if set, then print error messages to log file
      
         @return number of collected indices (includes excludeIndex at end of array)
      
         --------------------------------------------------------------------------------------- */
      
      int collectGbsaBondIndices( gbsaBonds* gbsaBond, int collectionBinMax,
                                  int* collectionBin, int excludeIndex, FILE* log );

      /**---------------------------------------------------------------------------------------
            
         Print state to log file (Simbios)
         
         @param title               title (optional)
         @param log                 print state to log file
            
         @return 0 always;
            
         --------------------------------------------------------------------------------------- */
      
      int logState( const char* title, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Get Vdw radii from parameter file (Simbios) 
      
         @param numberOfAtoms       number of atoms
         @param parameterFileName   parameterFileName
         @param top                 Gromacs topology data struct
         @param atomNames           array of atom names from GMX top data struct
         @param radii               array to store Macromodel radii for each atom
         @param log                 if set, then print error messages to log file
      
         @return 0 always
      
         --------------------------------------------------------------------------------------- */
      
      int getMacroModelAtomicRadii( int numberOfAtoms, const std::string parameterFileName,
                                    const t_topology* top, char*** atomNames,
                                    Real* radii, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Read file (Simbios) 
      
         @param fileName       file name
      
         @return vector of strings, one string per line
      
         --------------------------------------------------------------------------------------- */
      
      static StringVector* readFile( const std::string& fileName );
      
      /**---------------------------------------------------------------------------------------
      
         Parse parameter file (Simbios) 
      
         @param fileContents  file contents
      
         @return vector of strings, one string per line
      
         --------------------------------------------------------------------------------------- */
      
      static GbsaAtomParameterMap* parseParameterFile( const StringVector& fileContents );
      
      /**---------------------------------------------------------------------------------------
      
         Tokenize a string (static method) (Simbios)
      
         @param line                 string to tokenize
         @param tokenVector          upon return vector of tokens
         @param delimiter            token delimter
         @param clearTokenVector     if true, clear tokenVector
      
         @return 1
      
         --------------------------------------------------------------------------------------- */
      
      static int tokenizeString( const std::string& line, StringVector& tokenVector,
                                 const std::string& delimiter, int clearTokenVector );
            
      /**---------------------------------------------------------------------------------------
      
         Replacement of sorts for strtok() (static method) (Simbios)
         Used to parse parameter file lines
      
         Should be moved to Utilities file
      
         @param lineBuffer           string to tokenize
         @param delimiter            token delimter
      
         @return number of args; if return value equals maxTokens, then more tokens than allocated
      
         --------------------------------------------------------------------------------------- */
      
      static char* strsep( char** lineBuffer, const char* delimiter );
      
};
   
/**---------------------------------------------------------------------------------------
      
   Qsort/heapsort integer comparison (Simbios) 
      
   @parma a first value to compare
   @param b second value to compare

   @return -1, 0, 1
      
--------------------------------------------------------------------------------------- */

int integerComparison( const void *a, const void *b);

// ---------------------------------------------------------------------------------------

#endif // __GbsaParameters_H__
