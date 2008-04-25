
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

#ifndef __CpuGbsa_H__
#define __CpuGbsa_H__

#include "GbsaParameters.h"

// ---------------------------------------------------------------------------------------

class CpuGbsa {

   private:

      // used for direct calls from Gromacs

      static CpuGbsa* _cpuGbsa;

      // GBSA parameters

      GbsaParameters* _gbsaParameters;

      // flag to signal whether ACE approximation
      // is to be included

      bool _includeAceApproximation;

      // work arrays and their accessors

      float* _gbsaBornForce;
      float* _gbsaBornRadiiTemp;

      float* getGbsaBornForce( void );
      float* getGbsaBornRadiiTemp( void );

      // initialize data members (more than
      // one constructor, so centralize intialization here)

      void initializeDataMembers( void );

      // convert units for energy/force

      float _forceConversionFactor;
      float _energyConversionFactor;

      // Ed, 2007-04-27: Store the energy internally
      float _gbsaEnergy; 

   public:

       // constructors/destructor

       CpuGbsa( int numberOfAtoms, const t_topology* top, FILE* log );
       CpuGbsa( GbsaParameters* gbsaParameters );

       ~CpuGbsa( );

      // override of new/delete

      static void* operator new( size_t size ); 
      static void  operator delete( void *p );

      static void* operator new[]( size_t size ); 
      static void  operator delete[]( void *p );

      // accessor for the energy

		float getEnergy( void ) const { return _gbsaEnergy; }

      // accessor for static member _cpuGbsa

      static int setCpuGbsa( CpuGbsa* cpuGbsa ){ _cpuGbsa = cpuGbsa; return 0; };
      static CpuGbsa* getCpuGbsa( void ){ return _cpuGbsa; };
      static int deleteCpuGbsa( void );

      // accessors for Gbsa parameters

      GbsaParameters* getGbsaParameters( void ) const { return _gbsaParameters; };
      int setGbsaParameters( GbsaParameters* gbsaParameters ){ _gbsaParameters = gbsaParameters; return 0; };
 
      // accessors ACE approximation

      bool includeAceApproximation( void ) const { return _includeAceApproximation; };
      int setIncludeAceApproximation( bool includeAceApproximation ){ _includeAceApproximation = includeAceApproximation; return 0; };

      // energy/force conversion factors

      float getForceConversionFactor(  void  ) const { return _forceConversionFactor;  };  
      float getEnergyConversionFactor( void  ) const { return _energyConversionFactor; };  

      int setForceConversionFactor(  float forceConversionFactor  ){ _forceConversionFactor  = forceConversionFactor;  return 0; };  
      int setEnergyConversionFactor( float energyConversionFactor ){ _energyConversionFactor = energyConversionFactor; return 0; };  

      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces (output); if not set on input, then memory is allocated
         @param log               if set, then print error messages to log file
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      static float** computeGbsaForces( const rvec *atomCoordinates, const float* partialCharges,
                                        float** forces, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadii( const rvec *atomCoordinates, float* bornRadii );
      
      /**---------------------------------------------------------------------------------------
      
         Get exclusions for specific atom (Simbios)
      
         @param numberOfAtoms         number of atoms
         @param gbsaBondsArray        array of bonds
         @param atomI                 atom index of atom for which exclusions are to be set
         @param exclusionWorkArray    exclusionWorkArray[j] = 1, if atom j is to be excluded
                                      value may be null on input in which space is allocated
         @param previousIndex         previousIndex -- if < 0, then iniitialize all entries to 0
         @param log                   if set, then print error messages to log file
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      int getExclusionsForAtom( int numberOfAtoms, gbsaBonds** gbsaBondsArray,
                                int atomI, int* exclusionWorkArray, int previousIndex,
                                FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Set exclusions for specific atom (Simbios)
      
         @param gbsaBonds           gbsa bond
         @param setValue            set value
         @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
                                    value may be null on input in which space is allocated
         @param log                 if set, then print error messages to log file
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      int setExclusionValue( gbsaBonds* gbsaBonds, int setValue, int* exclusionWorkArray, FILE* log ) const;
         
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      float** computeBornEnergyForces( float* bornRadii, const rvec *atomCoordinates,
                                       const float* partialCharges, float** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param top               GMX t_topology struct
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    atom partial charges
         @param energy            energy
         @param debugOn           enable debug
         @param debugAtomI        debug flag for atomI (outer loop atom)
         @param debugAtomJ        debug flag for atomJ (inner loop atom)
         @param debugReport       flag signalling what quantities are to be saved for comparisons
         @param saveLoopForces    flag signalling whether intermediate results for each loop
                                  are to be retained
         @param printOn           print flag 
         @param unsetDebugValue   ?
         @param numberOfDebugStreams 
                                  number of debug streams
         @param debugStreams      array of debug streams
         @param forces            force array (may be null)
         @param log               if set, then print error messages to log file
      
         @return 0 always; fixed value for G_pol in array 'gPolFixed'
      
         --------------------------------------------------------------------------------------- */
      
      float** computeBornEnergyForcesDebug( const t_topology* top, 
                                            float* bornRadii, const rvec *atomCoordinates,
                                            const float* partialCharges, float* energy,
                                            bool debugOn, int debugAtomI, int debugAtomJ,
                                            int debugReport, bool* saveLoopForces,
                                            bool printOn, float unsetDebugValue, int numberOfDebugStreams,
                                            float** debugStreams, float** forces, FILE* log );
 
      /**---------------------------------------------------------------------------------------
      
         Write Born energy and forces (Simbios)
      
         @param top               GMX t_topology struct
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial atom charges
         @param forces            force array
         @param gbsaResultsFileName  output file name
         @param log               if set, then print error messages to log file
      
         @return 0 always
      
         --------------------------------------------------------------------------------------- */
      
      int writeBornEnergyForces( const t_topology* top, const rvec *atomCoordinates,
                                 const float* partialCharges, float** forces,
                                 const char* gbsaResultsFileName, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Print state to log file (Simbios)
      
         title             title (optional)
         log               print state to log file
      
         @return 0 always;
      
         --------------------------------------------------------------------------------------- */
      
      int logState( const char* title, FILE* log ) const;

      /**---------------------------------------------------------------------------------------
      
         Write Tinker xyz file (Simbios)
      
         @param numberOfAtoms      number of atoms
         @param atomCoordinates    atom coordinates
         @param header             header
         @param xyzFileName        output file name
         @param gbsaBondsArray     bond array -- used to print 1-2 bonds
         @param top                GMX t_topology struct
         @param log                if set, then print error messages to log file
      
         @return 0 unless error detected
      
         Currently no attempt is made to get the atom name/type to accurately 
         reflect the Tinker names/types. Rather method is used to output atoms
         in Gromacs order and then reorder those in a corresponding xyz file
         w/ the correct atom names/types so that they match the Gromacs order
         This makes it easier to compare results between Gromacs and Tinker
      
         --------------------------------------------------------------------------------------- */
      
      static int writeXyzFile( int numberOfAtoms, const rvec* atomCoordinates, 
                               const char* header, const char* xyzFileName,
                               gbsaBonds** gbsaBondsArray, const t_topology* top, FILE* log );
      
};

// ---------------------------------------------------------------------------------------

#endif // __CpuGbsa_H__
