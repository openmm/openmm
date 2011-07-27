
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

#ifndef __CpuObcSoftcore_H__
#define __CpuObcSoftcore_H__

#include "ObcSoftcoreParameters.h"
#include "gbsa/CpuImplicitSolvent.h"

// ---------------------------------------------------------------------------------------

class CpuObcSoftcore : public CpuImplicitSolvent {

   private:

      // GBSA/OBC parameters

      ObcSoftcoreParameters* _obcSoftcoreParameters;

      // arrays containing OBC chain derivative 

      std::vector<RealOpenMM> _obcChain;
      std::vector<RealOpenMM> _obcChainTemp;

      // initialize data members (more than
      // one constructor, so centralize intialization here)

      void _initializeObcDataMembers( void );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuObcSoftcore( ImplicitSolventParameters* obcSoftcoreParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuObcSoftcore( );

      /**---------------------------------------------------------------------------------------
      
         Return ObcSoftcoreParameters
      
         @return ObcSoftcoreParameters
      
         --------------------------------------------------------------------------------------- */

      ObcSoftcoreParameters* getObcSoftcoreParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         --------------------------------------------------------------------------------------- */

      void setObcSoftcoreParameters( ObcSoftcoreParameters* obcSoftcoreParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return OBC chain derivative: size = _implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      std::vector<RealOpenMM>& getObcChain( void );
      const std::vector<RealOpenMM>& getObcChainConst( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Return OBC chain temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      std::vector<RealOpenMM>& getObcChainTemp( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on OBC 

         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
         @param obcChain          output array of OBC chain derivative factors; if NULL,
                                  then ignored
      
         --------------------------------------------------------------------------------------- */
      
      void computeBornRadii( std::vector<OpenMM::RealVec>& atomCoordinates,  std::vector<RealOpenMM>& bornRadii );
      
      /**---------------------------------------------------------------------------------------
      
         Get nonpolar solvation force constribution via ACE approximation
      
         @param obcSoftcoreParameters     parameters
         @param vdwRadii                  Vdw radii
         @param bornRadii                 Born radii
         @param energy                    energy (output): value is incremented from input value 
         @param forces                    forces: values are incremented from input values
      
         --------------------------------------------------------------------------------------- */
      
      void computeAceNonPolarForce( const ObcSoftcoreParameters* obcSoftcoreParameters,
                                   const std::vector<RealOpenMM>& bornRadii, RealOpenMM* energy,
                                   std::vector<RealOpenMM>& forces ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on OBC 
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         --------------------------------------------------------------------------------------- */
      
      void computeBornEnergyForces( std::vector<RealOpenMM>& bornRadii, std::vector<OpenMM::RealVec>& atomCoordinates,
                                   const RealOpenMM* partialCharges, std::vector<OpenMM::RealVec>& forces );
};

// ---------------------------------------------------------------------------------------

#endif // __CpuObcSoftcore_H__
