
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

#ifndef __CpuImplicitSolvent_H__
#define __CpuImplicitSolvent_H__

#include "ImplicitSolventParameters.h"
#include <vector>

// ---------------------------------------------------------------------------------------

class OPENMM_EXPORT CpuImplicitSolvent {

   private:

      // used for direct calls 

      static CpuImplicitSolvent* _cpuImplicitSolvent;

      // parameters

      ImplicitSolventParameters* _implicitSolventParameters;

      // flag to signal whether ACE approximation
      // is to be included

      int _includeAceApproximation;

      // force index call 

      int _forceCallIndex;

      // work arrays

      std::vector<RealOpenMM> _bornForce;

      // Born radii and force

      std::vector<RealOpenMM> _bornRadii;
      std::vector<RealOpenMM> _tempBornRadii;

      // convert units for energy/force

      RealOpenMM _forceConversionFactor;
      RealOpenMM _energyConversionFactor;

      // Ed, 2007-04-27: Store the energy internally

      RealOpenMM _implicitSolventEnergy; 

      /**---------------------------------------------------------------------------------------
      
         Initialize data members -- potentially more than
         one constructor, so centralize intialization here
      
         --------------------------------------------------------------------------------------- */

      void _initializeDataMembers( void );

   protected:

      /**---------------------------------------------------------------------------------------
      
         Return implicitSolventBornForce, a work array of size _implicitSolventParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      std::vector<RealOpenMM>& getBornForce( void );

      /**---------------------------------------------------------------------------------------
      
         Return Born radii temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      std::vector<RealOpenMM>& getBornRadiiTemp( void );

      /**---------------------------------------------------------------------------------------
      
         Set energy 

         @param energy new energy
      
         --------------------------------------------------------------------------------------- */

      void setEnergy( RealOpenMM energy );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuImplicitSolvent( ImplicitSolventParameters* implicitSolventParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

      virtual ~CpuImplicitSolvent( );

      /**---------------------------------------------------------------------------------------
      
         Delete static _cpuImplicitSolvent object if set
      
         @return SimTKOpenMMCommon::DefaultReturn if _cpuImplicitSolvent was set; 
                 otherwise return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */

      static int deleteCpuImplicitSolvent( void );

      /**---------------------------------------------------------------------------------------
      
         Set static member _cpuImplicitSolvent
      
         --------------------------------------------------------------------------------------- */

      static void setCpuImplicitSolvent( CpuImplicitSolvent* cpuImplicitSolvent );

      /**---------------------------------------------------------------------------------------
      
         Get static member cpuImplicitSolvent
      
         @return static member cpuImplicitSolvent
      
         --------------------------------------------------------------------------------------- */
      
      static CpuImplicitSolvent* getCpuImplicitSolvent( void );

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */

		int getNumberOfAtoms( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get energy 
      
         @return energy
      
         --------------------------------------------------------------------------------------- */

		RealOpenMM getEnergy( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Return ImplicitSolventParameters
      
         @return ImplicitSolventParameters
      
         --------------------------------------------------------------------------------------- */
      
      ImplicitSolventParameters* getImplicitSolventParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         --------------------------------------------------------------------------------------- */
      
      void setImplicitSolventParameters( ImplicitSolventParameters* implicitSolventParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return flag signalling whether AceApproximation for nonpolar term is to be included
      
         @return flag
      
         --------------------------------------------------------------------------------------- */

      int includeAceApproximation( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether AceApproximation is to be included
      
         @param includeAceApproximation new includeAceApproximation value
      
         --------------------------------------------------------------------------------------- */
      
      void setIncludeAceApproximation( int includeAceApproximation );

      /**---------------------------------------------------------------------------------------
      
         Return ForceConversionFactor for units
      
         @return ForceConversionFactor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getForceConversionFactor(  void  ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ForceConversionFactor
      
         @param ForceConversionFactor (units conversion)
      
         --------------------------------------------------------------------------------------- */
      
      void setForceConversionFactor(  RealOpenMM forceConversionFactor  );

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
      
      void setEnergyConversionFactor( RealOpenMM energyConversionFactor );

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
      
         Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      std::vector<RealOpenMM>& getBornRadii( void );
      
      /**---------------------------------------------------------------------------------------
      
         Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const std::vector<RealOpenMM>& getBornRadiiConst( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces (output); if not set on input, then memory is allocated
         @param updateBornRadii   if set, then Born radii are updated for current configuration; 
                                  otherwise radii correspond to configuration from previous iteration
      
         --------------------------------------------------------------------------------------- */
      
      void computeImplicitSolventForces( std::vector<OpenMM::RealVec>& atomCoordinates,
                                               const RealOpenMM* partialCharges,
                                               std::vector<OpenMM::RealVec>& forces, int updateBornRadii = 0 );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param atomCoordinates   atomic coordinates dimension: [0-numberAtoms-1][0-2]
         @param bornRadii         output array of Born radii
         @param obcChain          output array of OBC chain derivative
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      virtual void computeBornRadii( std::vector<OpenMM::RealVec>& atomCoordinates, std::vector<RealOpenMM>& bornRadii ) = 0;
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      virtual void computeBornEnergyForces( std::vector<RealOpenMM>& bornRadii, std::vector<OpenMM::RealVec>& atomCoordinates,
                                           const RealOpenMM* partialCharges, std::vector<OpenMM::RealVec>& forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get nonpolar solvation force constribution via ACE approximation
      
         @param implicitSolventParameters parameters
         @param bornRadii                 Born radii
         @param energy                    energy (output): value is incremented from input value 
         @param forces                    forces: values are incremented from input values
      
         --------------------------------------------------------------------------------------- */
      
      void computeAceNonPolarForce( const ImplicitSolventParameters* implicitSolventParameters,
                                   const std::vector<RealOpenMM>& bornRadii, RealOpenMM* energy,
                                   std::vector<RealOpenMM>& forces ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state 
         
         @param title               title (optional)
            
         @return string containing state
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;
};

// ---------------------------------------------------------------------------------------

#endif // __CpuImplicitSolvent_H__
