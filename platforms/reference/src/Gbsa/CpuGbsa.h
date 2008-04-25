
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

      // used for direct calls 

      static CpuGbsa* _cpuGbsa;

      // GBSA parameters

      GbsaParameters* _gbsaParameters;

      // flag to signal whether ACE approximation
      // is to be included

      bool _includeAceApproximation;

      // force index call 

      int _forceCallIndex;

      // work arrays

      RealOpenMM* _gbsaBornForce;
      RealOpenMM* _gbsaBornRadiiTemp;

      // convert units for energy/force

      RealOpenMM _forceConversionFactor;
      RealOpenMM _energyConversionFactor;

      // Ed, 2007-04-27: Store the energy internally

      RealOpenMM _gbsaEnergy; 

      /**---------------------------------------------------------------------------------------
      
         Initialize data members -- potentially more than
         one constructor, so centralize intialization here
      
         --------------------------------------------------------------------------------------- */

      void _initializeDataMembers( void );

      /**---------------------------------------------------------------------------------------
      
         Return gbsaBornForce, a work array of size _gbsaParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getGbsaBornForce( void );

      /**---------------------------------------------------------------------------------------
      
         Return gbsaBornRadiiTemp, a work array of size _gbsaParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getGbsaBornRadiiTemp( void );

      /**---------------------------------------------------------------------------------------
      
         Set energy 

         @param energy new energy
      
         ireturn SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

		int setEnergy( RealOpenMM energy );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param gbsaParameters    GbsaParameters reference
      
         @return CpuGbsa object
      
         --------------------------------------------------------------------------------------- */

       CpuGbsa( GbsaParameters* gbsaParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuGbsa( );

      // override of new/delete -- used when run in PS3 framework(?)

      // static void* operator new( size_t size ); 
      // static void  operator delete( void *p );

      // static void* operator new[]( size_t size ); 
      // static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Delete static _cpuGbsa object if set
      
         @return SimTKOpenMMCommon::DefaultReturn if _cpuGbsa was set; 
                 otherwise return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int deleteCpuGbsa( void );

      /**---------------------------------------------------------------------------------------
      
         Set static member _cpuGbsa
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      static int setCpuGbsa( CpuGbsa* cpuGbsa );

      /**---------------------------------------------------------------------------------------
      
         Get static member cpuGbsa
      
         @return static member cpuGbsa
      
         --------------------------------------------------------------------------------------- */
      
      static CpuGbsa* getCpuGbsa( void );

      /**---------------------------------------------------------------------------------------
      
         Get energy 
      
         @return energy
      
         --------------------------------------------------------------------------------------- */

		RealOpenMM getEnergy( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Return GbsaParameters
      
         @return GbsaParameters
      
         --------------------------------------------------------------------------------------- */
      
      GbsaParameters* getGbsaParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set GbsaParameters
      
         @param GbsaParameters
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setGbsaParameters( GbsaParameters* gbsaParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return flag signalling whether AceApproximation for nonpolar term is to be included
      
         @return flag
      
         --------------------------------------------------------------------------------------- */

      int includeAceApproximation( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether AceApproximation is to be included
      
         @param includeAceApproximation new includeAceApproximation value
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setIncludeAceApproximation( bool includeAceApproximation );

      /**---------------------------------------------------------------------------------------
      
         Return ForceConversionFactor for units
      
         @return ForceConversionFactor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getForceConversionFactor(  void  ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ForceConversionFactor
      
         @param ForceConversionFactor (units conversion)
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setForceConversionFactor(  RealOpenMM forceConversionFactor  );

      /**---------------------------------------------------------------------------------------
      
         Return EnergyConversionFactor for units
      
         @return EnergyConversionFactor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getEnergyConversionFactor( void  ) const;

      /**---------------------------------------------------------------------------------------
      
         Set EnergyConversionFactor
      
         @param EnergyConversionFactor (units conversion)
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setEnergyConversionFactor( RealOpenMM energyConversionFactor );

      /**---------------------------------------------------------------------------------------
      
         Return ForceCallIndex -- number of times forces have been calculated
      
         @return ForceCallIndex
      
         --------------------------------------------------------------------------------------- */
      
      int getForceCallIndex( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Increment ForceCallIndex
      
         @return incremented forceCallIndex
      
         --------------------------------------------------------------------------------------- */
      
      int incrementForceCallIndex( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param born radii        Born radii (may be NULL)
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces (output); if not set on input, then memory is allocated
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM** computeGbsaForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                             const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param atomCoordinates   atomic coordinates dimension: [0-numberAtoms-1][0-2]
         @param bornRadii         output array of Born radii
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii );
      
      /**---------------------------------------------------------------------------------------
      
         Get exclusions for specific atom (Simbios)
      
         @param numberOfAtoms         number of atoms
         @param gbsaBondsArray        array of bonds
         @param atomI                 atom index of atom for which exclusions are to be set
         @param exclusionWorkArray    exclusionWorkArray[j] = 1, if atom j is to be excluded
                                      value may be null on input in which space is allocated
         @param previousIndex         previousIndex -- if < 0, then iniitialize all entries to 0
      
         @return SimTKOpenMMCommon::DefaultReturn; abort if exclusionWorkArray is not set
      
         --------------------------------------------------------------------------------------- */
      
      int getExclusionsForAtom( int numberOfAtoms, gbsaBonds** gbsaBondsArray,
                                int atomI, int* exclusionWorkArray, int previousIndex ) const;

      /**---------------------------------------------------------------------------------------
      
         Set exclusions for specific atom (Simbios)
      
         @param gbsaBonds           gbsa bond
         @param setValue            set value
         @param exclusionWorkArray  exclusionWorkArray[j] = 1, if atom j is to be excluded
                                    value may be null on input in which space is allocated
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setExclusionValue( gbsaBonds* gbsaBonds, int setValue, int* exclusionWorkArray ) const;
         
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM** computeBornEnergyForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                            const RealOpenMM* partialCharges, RealOpenMM** forces );
      
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
      
/*
      RealOpenMM** computeBornEnergyForcesDebug( const t_topology* top, 
                                                 RealOpenMM* bornRadii, const rvec *atomCoordinates,
                                                 const RealOpenMM* partialCharges, RealOpenMM* energy,
                                                 bool debugOn, int debugAtomI, int debugAtomJ,
                                                 int debugReport, bool* saveLoopForces,
                                                 bool printOn, RealOpenMM unsetDebugValue, int numberOfDebugStreams,
                                                 RealOpenMM** debugStreams, RealOpenMM** forces ); */
 
      /**---------------------------------------------------------------------------------------
      
         Write Born energy and forces (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial atom charges
         @param forces            force array
         @param gbsaResultsFileName  output file name
         @param log               if set, then print error messages to log file
      
         @return 0 always
      
         --------------------------------------------------------------------------------------- */
      
      int writeBornEnergyForces( const RealOpenMM** atomCoordinates,
                                 const RealOpenMM* partialCharges, const RealOpenMM** forces,
                                 const std::string& gbsaResultsFileName ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state 
         
         @param title               title (optional)
            
         @return string containing state
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Write Tinker xyz file (Simbios)
      
         @param numberOfAtoms      number of atoms
         @param atomCoordinates    atom coordinates
         @param atomNames          atom names
         @param header             header
         @param xyzFileName        output file name
         @param bondsArray         bond array -- used to print 1-2 bonds
      
         @return 0 unless error detected
      
         Currently no attempt is made to get the atom name/type to accurately 
         reflect the Tinker names/types. Rather method is used to output atoms
         in Gromacs order and then reorder those in a corresponding xyz file
         w/ the correct atom names/types so that they match the Gromacs order
         This makes it easier to compare results between Gromacs and Tinker
      
         --------------------------------------------------------------------------------------- */
      
      static int writeXyzFile( int numberOfAtoms, RealOpenMM** atomCoordinates, 
                               char** atomNames,
                               const std::string& header, const std::string& xyzFileName,
                               gbsaBonds** bondsArray );
      
};

// ---------------------------------------------------------------------------------------

#endif // __CpuGbsa_H__
