
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

#ifndef __CpuObc_H__
#define __CpuObc_H__

#include "ObcParameters.h"

// ---------------------------------------------------------------------------------------

class CpuObc {

   private:

      // GBSA/OBC parameters

      ObcParameters* _obcParameters;

      // arrays containing OBC chain derivative 

      RealOpenMMVector _obcChain;

      // flag to signal whether ACE approximation
      // is to be included

      int _includeAceApproximation;


   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuObc( ObcParameters* obcParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuObc( );

      /**---------------------------------------------------------------------------------------
      
         Return ObcParameters
      
         @return ObcParameters
      
         --------------------------------------------------------------------------------------- */

      ObcParameters* getObcParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         --------------------------------------------------------------------------------------- */

      void setObcParameters( ObcParameters* obcParameters );
 
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
      
         Return OBC chain derivative: size = _implicitSolventParameters->getNumberOfAtoms()
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMMVector& getObcChain( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on OBC 

         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      void computeBornRadii( const std::vector<OpenMM::RealVec>& atomCoordinates, RealOpenMMVector& bornRadii );
      
      /**---------------------------------------------------------------------------------------
        
         Get nonpolar solvation force constribution via ACE approximation
        
         @param obcParameters parameters
         @param vdwRadii                  Vdw radii
         @param bornRadii                 Born radii
         @param energy                    energy (output): value is incremented from input value 
         @param forces                    forces: values are incremented from input values
        
            --------------------------------------------------------------------------------------- */
        
      void computeAceNonPolarForce( const ObcParameters* obcParameters, const RealOpenMMVector& bornRadii, 
                                    RealOpenMM* energy, RealOpenMMVector& forces ) const;
        
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on OBC 
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM computeBornEnergyForces( const std::vector<OpenMM::RealVec>& atomCoordinates,
                                          const RealOpenMMVector& partialCharges, std::vector<OpenMM::RealVec>& forces );
      
    /**---------------------------------------------------------------------------------------
    
        Print Obc parameters, radii, forces, ...
    
        @param atomCoordinates     atomic coordinates
        @param partialCharges      partial charges
        @param bornRadii           Born radii (may be empty)
        @param bornForces          Born forces (may be empty)
        @param forces              forces (may be empty)
        @param idString            id string (who is calling)
        @param log                 log file
    
        --------------------------------------------------------------------------------------- */
    
    void printObc( const std::vector<OpenMM::RealVec>& atomCoordinates,
                   const RealOpenMMVector& partialCharges,
                   const RealOpenMMVector& bornRadii,
                   const RealOpenMMVector& bornForces,
                   const std::vector<OpenMM::RealVec>& forces,
                   const std::string& idString, FILE* log );
    
};

// ---------------------------------------------------------------------------------------

#endif // __CpuObc_H__
