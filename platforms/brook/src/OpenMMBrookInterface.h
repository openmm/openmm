#ifndef OPENMM_BROOK_INTERFACE_H_
#define OPENMM_BROOK_INTERFACE_H_

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
#include "NonbondedForce.h"
#include "OpenMMContext.h"
#include "System.h"
#include "ReferencePlatform.h"
#include "VerletIntegrator.h"

class BrookBondParameters;

namespace OpenMM {

/**
 * OpenMM-Brook interface methods
 */
class OpenMMBrookInterface {

   public:
  
      OpenMMBrookInterface( void );
  
      ~OpenMMBrookInterface();
  
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
                       //NonbondedMethod nonbondedMethod,
                       double nonbondedCutoff, double periodicBoxSize[3] );

  
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
       * @return the potential energy due to the NonbondedForce
       *
       * Currently always return 0.0 since energies not calculated on gpu
       */

      double executeEnergy( const Stream& positions );
      double executeEnergyOld( const Stream& positions );

      /**
       * Get reference Context
       * 
       * @param numberOfAtoms  number of atoms
       *
       * @return  OpenMMContext
       *
       */
      
      OpenMMContext* getReferenceOpenMMContext( int numberOfAtoms );
      
      
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
      
      /** 
       * Compute forces
       *
       * @param context OpenMMContextImpl context 
       *
       */
      
      void computeForces( OpenMMContextImpl& context );
      
      /** 
       * Compute energy
       *
       * @param context OpenMMContextImpl context 
       *
       */
      
      float computeEnergy( OpenMMContextImpl& context );
      
      /** 
       * Set trigger Force Kernel
       *
       * @param triggerForceKernel kernel to calculate force
       *
       */
      
      void setTriggerForceKernel( KernelImpl* triggerForceKernel );
      
      /** 
       * Set trigger Energy Kernel
       *
       * @param triggerEnergyKernel kernel to calculate energy
       *
       */
      
      void setTriggerEnergyKernel( KernelImpl* triggerForceKernel );
      
      /** 
       * Get trigger Force Kernel
       *
       * @return triggerForceKernel kernel to calculate force
       *
       */
      
      KernelImpl* getTriggerForceKernel( void ) const;
      
      /** 
       * Get trigger Energy Kernel
       *
       * @return triggerEnergyKernel kernel to calculate energy
       *
       */
      
      KernelImpl* getTriggerEnergyKernel( void ) const;
      
      /** 
       * Set BrookBondParameters for harmonic angle force
       * 
       * @param brookBondParameters brookBondParameters for BrookBondParameters for harmonic angle force
       *
       * @return DefaultReturnValue
       *
       */
      
      int setHarmonicBondForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Set BrookBondParameters for harmonic angle force
       * 
       * @param brookBondParameters brookBondParameters for BrookBondParameters for harmonic angle force
       *
       * @return DefaultReturnValue
       *
       */
      
      int setHarmonicAngleForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Set BrookBondParameters for proper dihedral force
       * 
       * @param brookBondParameters brookBondParameters for proper dihedral force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setPeriodicTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Set BrookBondParameters for RB dihedral force
       * 
       * @param brookBondParameters brookBondParameters for RB force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setRBTorsionForceParameters( BrookBondParameters* brookBondParameters );
      
      /** 
       * Set BrookBondParameters for LJ 14 force
       * 
       * @param brookBondParameters brookBondParameters for LJ 14 force
       *
       * @return  DefaultReturnValue
       *
       */
      
      int setLJ14( BrookBondParameters* brookBondParameters );
      
   private:
   
      enum BondParameterIndices { HarmonicBondIndex, HarmonicAngleIndex, PeriodicTorsionForceIndex, RbTorsionForceIndex, LJ14Index, LastBondForce };

      // log file reference

      FILE* _log;

       // number of atoms

       int _numberOfAtoms;
   
       // Brook bonded & nonbonded

       BrookBonded* _brookBonded;
       BrookNonBonded* _brookNonBonded;

       BrookBondParameters* _bondParameters[LastBondForce];

       // used to calculate energy

       NonbondedForce*       _refForceField;
       System*               _refSystem;
       OpenMMContext*        _refOpenMMContext;
       ReferencePlatform*    _referencePlatform;
       VerletIntegrator*     _refVerletIntegrator;

      /** 
       * Set BrookBondParameters 
       * 
       * @param index
       * @param brookBondParameters brookBondParameters for BrookBondParameters
       *
       * @return DefaultReturnValue
       *
       */
      
      int _setBondParameters( BondParameterIndices index, BrookBondParameters* brookBondParameters );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INTERFACE_H_ */
