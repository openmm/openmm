#ifndef OPENMM_BROOK_CALC_STANDARD_MM_FORCEFIELD_KERNEL_H_
#define OPENMM_BROOK_CALC_STANDARD_MM_FORCEFIELD_KERNEL_H_

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
#include "BrookBonded.h"
#include "BrookNonBonded.h"
#include "StandardMMForceField.h"

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class BrookCalcStandardMMForceFieldKernel : public CalcStandardMMForceFieldKernel {

   public:
  
      BrookCalcStandardMMForceFieldKernel( std::string name, const Platform& platform );
  
      ~BrookCalcStandardMMForceFieldKernel();
  
      /**
       * Initialize the kernel, setting up the values of all the force field parameters.
       * 
       * @param bondIndices               the two atoms connected by each bond term
       * @param bondParameters            the force parameters (length, k) for each bond term
       * @param angleIndices              the three atoms connected by each angle term
       * @param angleParameters           the force parameters (angle, k) for each angle term
       * @param periodicTorsionIndices    the four atoms connected by each periodic torsion term
       * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
       * @param rbTorsionIndices          the four atoms connected by each Ryckaert-Bellemans torsion term
       * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
       * @param bonded14Indices           each element contains the indices of two atoms whose nonbonded interactions should be reduced since
       *                                  they form a bonded 1-4 pair
       * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
       * @param coulomb14Scale            the factor by which Coulomb interactions should be reduced for bonded 1-4 pairs
       * @param exclusions                the i'th element lists the indices of all atoms with which the i'th atom should not interact through
       *                                  nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
       *                                  the standard nonbonded calculation.
       * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
       * @param nonbondedMethod           the method to use for handling long range nonbonded interactions
       * @param nonbondedCutoff           the cutoff distance for nonbonded interactions (if nonbondedMethod involves a cutoff)
       * @param periodicBoxSize           the size of the periodic box (if nonbondedMethod involves a periodic boundary conditions)
       *
       */

      void initialize( const std::vector<std::vector<int> >& bondIndices,      const std::vector<std::vector<double> >& bondParameters,
                       const std::vector<std::vector<int> >& angleIndices,     const std::vector<std::vector<double> >& angleParameters,
                       const std::vector<std::vector<int> >& periodicTorsionIndices, const std::vector<std::vector<double> >& periodicTorsionParameters,
                       const std::vector<std::vector<int> >& rbTorsionIndices, const std::vector<std::vector<double> >& rbTorsionParameters,
                       const std::vector<std::vector<int> >& bonded14Indices,  double lj14Scale, double coulomb14Scale,
                       const std::vector<std::set   <int> >& exclusions,       const std::vector<std::vector<double> >& nonbondedParameters,
                       NonbondedMethod nonbondedMethod, double nonbondedCutoff, double periodicBoxSize[3] );
  
      /**
       * Execute the kernel to calculate the forces.
       * 
       * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
       * @param forces      a Stream of type Double3 containing the force (x, y, z) on each atom.  On entry, this contains the forces that
       *                    have been calculated so far.  The kernel should add its own forces to the values already in the stream.
       */

      void executeForces( const Stream& positions, Stream& forces );
  
      /**
       * Execute the kernel to calculate the energy.
       * 
       * @param positions   a Stream of type Double3 containing the position (x, y, z) of each atom
       *
       * @return the potential energy due to the StandardMMForceField
       *
       * Currently always return 0.0 since energies not calculated on gpu
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
   
       // Brook bonded & nonbonded

       BrookBonded* _brookBonded;
       BrookNonBonded* _brookNonBonded;

       // used to calculate energy

       StandardMMForceField* _refForceField;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_CALC_STANDARD_MM_FORCEFIELD_KERNEL_H_ */
