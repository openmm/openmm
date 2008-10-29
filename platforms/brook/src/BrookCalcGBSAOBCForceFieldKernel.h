#ifndef OPENMM_BROOK__CALCL_GBSAOBC_FORCEFIELD_KERNEL_H_
#define OPENMM_BROOK__CALCL_GBSAOBC_FORCEFIELD_KERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "kernels.h"
#include "../../reference/src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "BrookGbsa.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
 */
class BrookCalcGBSAOBCForceKernel : public CalcGBSAOBCForceKernel {

   public:
  
      /**
       * BrookCalcGBSAOBCForceKernel constructor
       */

      BrookCalcGBSAOBCForceKernel( std::string name, const Platform& platform );
  
      /**
       * BrookCalcGBSAOBCForceKernel destructor
       */

      ~BrookCalcGBSAOBCForceKernel();
  
      /**
       * Initialize the kernel, setting up the values of all the force field parameters.
       * 
       * @param atomParameters            vector containing atom index, charge, radius, scalingFactor
       * @param solventDielectric         solvent dielectric
       * @param soluteDielectric          solute dielectric
       *
       */

      void initialize( const std::vector<std::vector<double> >& atomParameters, double solventDielectric, double soluteDielectric );
  
      /**
       * Execute the kernel to calculate the forces.
       * 
       * @param positions   atom coordiantes
       * @param forces      output forces
       *
       */

      void executeForces( const Stream& positions, Stream& forces );
  
      /**
       * Execute the kernel to calculate the energy.
       * 
       * @param positions   atom positions
       *
       * @return  potential energy due to the NonbondedForce
       * Currently always return 0.0 since energies not calculated on gpu
       *
       */

      double executeEnergy( const Stream& positions );

      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
   private:
   
      // log file reference

      FILE* _log;

      // number of atoms

      int _numberOfAtoms;
   
      // Brook Gbsa

      BrookGbsa* _brookGbsa;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK__CALCL_GBSAOBC_FORCEFIELD_KERNEL_H_ */
