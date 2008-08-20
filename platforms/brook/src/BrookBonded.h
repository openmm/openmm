#ifndef OPENMM_BROOK_BONDED_H_
#define OPENMM_BROOK_BONDED_H_

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

#include <vector>

#include "BrookFloatStreamInternal.h"
#include "BrookIntStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"
#include "OpenMMContext.h"

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class BrookBonded : public BrookCommon {

   public:
  
      /** 
       *
       * BrookBonded constructor
       * 
       */
      
      BrookBonded( );
  
      /** 
       *
       * BrookBonded desstructor
       * 
       */
      
      ~BrookBonded( );
  
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
       * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
       * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
       * @param coulomb14Scale            the factor by which Coulomb interactions should be reduced for bonded 1-4 pairs
       * @param log                       log reference
       *
       * @return nonzero value if error
       *
       */

      int setup( int numberOfAtoms, 
                 const std::vector<std::vector<int> >& bondIndices,            const std::vector<std::vector<double> >& bondParameters,
                 const std::vector<std::vector<int> >& angleIndices,           const std::vector<std::vector<double> >& angleParameters,
                 const std::vector<std::vector<int> >& periodicTorsionIndices, const std::vector<std::vector<double> >& periodicTorsionParameters,
                 const std::vector<std::vector<int> >& rbTorsionIndices,       const std::vector<std::vector<double> >& rbTorsionParameters,
                 const std::vector<std::vector<int> >& bonded14Indices,        const std::vector<std::vector<double> >& nonbondedParameters,
                 double lj14Scale, double coulomb14Scale, const Platform& platform );

      /**
       * Get inverse map stream width
       * 
       * @return stream width
       *
       */

      int getInverseMapStreamWidth( void ) const;
      
      /**
       * Return number of parameter streams
       * 
       * @return number of parameter streams
       *
       */

      int getNumberOfParameterStreams( void ) const; 

      /**
       * Return number of force streams
       * 
       * @return number of force streams
       *
       */

      int getNumberOfForceStreams( void ) const; 

      /**
       * Return stream count for specifed index (i, j, k, l)
       * 
       * @return stream count for specifed index
       *
       * @throw throws OpenMMException if index out of range
       */

      int getInverseMapStreamCount( int index ) const; 

      /**
       * Return max stream count
       * 
       * @return max stream count
       */

      int getMaxInverseMapStreamCount( void ) const; 

      /**
       * Return max stream count for specified index
       * 
       * @param index index of force stream
       *
       * @return max stream count
       *
       */

      int getMaxInverseMapStreamCount( int index ) const; 

      /**
       * Return Brook stream handle
       * 
       * @return 
       */

      BrookFloatStreamInternal* getBrookAtomIndices( void ) const; 

      /**
       * Get LJ 14 scale factor
       *
       * @return  LJ 14 scaling (fudge) factor
       * 
       */

       BrookOpenMMFloat getLJ_14Scale( void ) const;

      /** 
       * Get Coulomb factor
       * 
       * @return Coulomb factor
       *
       */
      
      BrookOpenMMFloat getCoulombFactor( void ) const;

      /** 
       * Get bonded atom indices stream
       * 
       * @return  atom indices stream
       *
       */
      
      BrookFloatStreamInternal* getAtomIndicesStream( void ) const;
      
      /** 
       * Get array of bonded parameter streams
       * 
       * @return  array of bonded parameter streams
       *
       */
      
      BrookFloatStreamInternal** getBondedParameterStreams( void );
      
      /** 
       * Get array of force streams
       * 
       * @return  array
       *
       */
      
      BrookFloatStreamInternal** getBondedForceStreams( void );
      
      /** 
       * Get array of inverse map streams
       * 
       * @param index  array index 
       *
       * @return  array inverse map streams
       *
       */
      
      BrookFloatStreamInternal** getInverseStreamMapsStreams( int index );
      
      /** 
       * Return true if force[index] stream is set
       *
       * @param    index into force stream
       * @return   true if index is valid && force[index] stream is set; else false
       *
       */
      
      int isForceStreamSet( int index ) const;
       
      /** 
       * Return true if paramsterSet[index] stream is set
       *
       * @param    index into parameter stream
       *
       * @return   true if index is valid && paramsterSet[index] stream is set; else false
       *
       */
      
      int isParameterStreamSet( int index ) const;
          
      /** 
       * Return string showing if all inverse map streams are set
       *
       * @param    index into inverse map stream array
       *
       * @return   informative string
       *
       */
      
      std::string checkInverseMapStream( int index ) const;

      /* 
       * Get contents of object
       *
       *
       * @param level   level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContentsString( int level = 0 ) const;

   private:
   
      static const int NumberOfParameterStreams = 5;
      static const int NumberOfForceStreams     = 4;
      static const int MaxNumberOfInverseMaps   = 9;

      enum { StreamI, StreamJ, StreamK, StreamL, StreamMax };

      // inverse map stream width

      int _inverseMapStreamWidth;

      // actual max number of inverse maps 

      int _maxNumberOfInverseMaps;

      // scale factors for 1-4 ixn

      BrookOpenMMFloat _ljScale;
      BrookOpenMMFloat _coulombFactor;

      // streams

      BrookFloatStreamInternal* _atomIndicesStream;
      BrookFloatStreamInternal* _bondedParameters[NumberOfParameterStreams];
      BrookFloatStreamInternal* _bondedForceStreams[NumberOfForceStreams];
      BrookFloatStreamInternal* _chargeStream;
      BrookFloatStreamInternal* _inverseStreamMaps[NumberOfForceStreams][MaxNumberOfInverseMaps];

      int _maxInverseMapStreamCount[NumberOfForceStreams];
      int _inverseMapStreamCount[NumberOfForceStreams];

      // helper methods in setup of parameters

      void flipQuartet( int ibonded, int *atoms );
      int matchTorsion( int i, int j, int k, int l, int nbondeds, int *atoms );
      int matchAngle( int i, int j, int k, int nbondeds, int *atoms, int *flag );
      int matchBond( int i, int j, int nbondeds, int *atoms, int *flag );
      int matchPair( int i, int j, int nbondeds, int *atoms );

      /**
       * Setup Ryckaert-Bellemans parameters/atom indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param atoms                     array of atom indices
       * @param params                    arrays of bond parameters
       * @param rbTorsionIndices          the four atoms connected by each Ryckaert-Bellemans torsion term
       * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
       *
       * @return nonzero value if error
       *
       */

      int addRBDihedrals( int *nbondeds, int *atoms, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& rbTorsionIndices,  
                          const std::vector<std::vector<double> >& rbTorsionParameters );

      /**
       * Setup periodic torsion parameters/atom indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param atoms                     array of atom indices
       * @param params                    arrays of bond parameters
       * @param periodicTorsionIndices    the four atoms connected by each periodic torsion term
       * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
       *
       * @return nonzero value if error
       *
       */

      int addPDihedrals( int *nbondeds, int *atoms, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& periodicTorsionIndices,  
                          const std::vector<std::vector<double> >& periodicTorsionParameters );

      /**
       * Setup angle bond parameters/atom indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param atoms                     array of atom indices
       * @param params                    arrays of bond parameters
       * @param angleIndices              the angle bond atom indices
       * @param angleParameters           the angle parameters (angle in radians, force constant)
       *
       * @return nonzero value if error
       *
       */

      int addAngles( int *nbondeds, int *atoms, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& angleIndices,
                     const std::vector<std::vector<double> >& angleParameters );

      /**
       * Setup harmonic bond parameters/atom indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param atoms                     array of atom indices
       * @param params                    arrays of bond parameters
       * @param bondIndices               two harmonic bond atom indices
       * @param bondParameters            the force parameters (distance, k)
       *
       * @return nonzero value if error
       *
       */

      int addBonds( int *nbondeds, int *atoms, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& bondIndices,
                    const std::vector<std::vector<double> >& bondParameters );

      /**
       * Setup LJ/Coulomb 1-4 parameters/atom indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param atoms                     array of atom indices
       * @param params                    arrays of bond parameters
       * @param charges                   array of charges
       * @param bonded14Indices           each element contains the indices of two atoms whose nonbonded interactions should be reduced since
       *                                  they form a bonded 1-4 pair
       * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
       * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
       *
       * @return nonzero value if error
       *
       */

      int addPairs( int *nbondeds, int *atoms, BrookOpenMMFloat* params[], BrookOpenMMFloat* charges,
                    const std::vector<std::vector<int> >& bonded14Indices, const std::vector<std::vector<double> >& nonbondedParameters,
                    double lj14Scale, double coulombScale );
      
      /**
       * Create and load inverse maps for bonded ixns
       * 
       * @param nbondeds                  number of bonded entries
       * @param natoms                    number of atoms
       * @param atoms                     arrays of atom indices (atoms[numberOfBonds][4])
       * @param platform                  BrookPlatform reference
       * @param log                       log file reference (optional)
       *
       * @return nonzero value if error
       *
       */

      int loadInvMaps( int nbondeds, int natoms, int *atoms, const BrookPlatform& platform );
      
      /**
       * Validate inverse map count
       * 
       * @param index index to check (1-4) 
       * @param count count for index
       *
       * @return -1 if count invalid
       *
       * @tthrow OpenMMException exeception if count invalid
       *
       */

      int validateInverseMapStreamCount( int index, int count ) const;
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_BONDED_H_ */
