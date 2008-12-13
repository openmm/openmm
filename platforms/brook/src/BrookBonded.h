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

#include "BrookStreamImpl.h"
//#include "BrookPlatform.h"
#include "BrookCommon.h"
#include "OpenMMContext.h"
#include "BrookBondParameters.h"

namespace OpenMM {

/**
 * This kernel is invoked by NonbondedForce to calculate the forces acting on the system.
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
       * @param bondIndices               the two particles connected by each bond term
       * @param bondParameters            the force parameters (length, k) for each bond term
       * @param angleIndices              the three particles connected by each angle term
       * @param angleParameters           the force parameters (angle, k) for each angle term
       * @param periodicTorsionIndices    the four particles connected by each periodic torsion term
       * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
       * @param rbTorsionIndices          the four particles connected by each Ryckaert-Bellemans torsion term
       * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
       * @param bonded14Indices           each element contains the indices of two particles whose nonbonded interactions should be reduced since
       *                                  they form a bonded 1-4 pair
       * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each particle
       * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
       * @param coulomb14Scale            the factor by which Coulomb interactions should be reduced for bonded 1-4 pairs
       * @param log                       log reference
       *
       * @return nonzero value if error
       *
       */

/*
      int setup( int numberOfParticles, 
                 const std::vector<std::vector<int> >& bondIndices,            const std::vector<std::vector<double> >& bondParameters,
                 const std::vector<std::vector<int> >& angleIndices,           const std::vector<std::vector<double> >& angleParameters,
                 const std::vector<std::vector<int> >& periodicTorsionIndices, const std::vector<std::vector<double> >& periodicTorsionParameters,
                 const std::vector<std::vector<int> >& rbTorsionIndices,       const std::vector<std::vector<double> >& rbTorsionParameters,
                 const std::vector<std::vector<int> >& bonded14Indices,        const std::vector<std::vector<double> >& nonbondedParameters,
                 double lj14Scale, double coulomb14Scale, const Platform& platform );
*/

      int setup( int numberOfParticles,
                 BrookBondParameters* harmonicBondBrookBondParameters,
                 BrookBondParameters* harmonicAngleBrookBondParameters,
                 BrookBondParameters* periodicTorsionBrookBondParameters,
                 BrookBondParameters* rbTorsionBrookBondParameters,
                 BrookBondParameters* nonBonded14ForceParameters,  double lj14Scale, double coulombScale, int particleStreamWidth, int particleStreamSize );

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

      BrookFloatStreamInternal* getBrookParticleIndices( void ) const; 

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
       * Get bonded particle indices stream
       * 
       * @return  particle indices stream
       *
       */
      
      BrookFloatStreamInternal* getParticleIndicesStream( void ) const;
      
      /** 
       * Get bonded charge stream
       * 
       * @return  charge stream
       *
       */
      
      BrookFloatStreamInternal* getChargeStream( void ) const;
      
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

      /** 
       * Compute forces
       * 
       */
      
      void computeForces( BrookStreamImpl& positionStream, BrookStreamImpl& forceStream );
   
      /** 
       * Return SetupCompleted flag
       *
       * @return SetupCompleted flag
       */
      
      int isSetupCompleted( void ) const;
   
   private:
   
      static const int NumberOfParameterStreams = 5;
      static const int NumberOfForceStreams     = 4;
      static const int MaxNumberOfInverseMaps   = 9;

      enum { StreamI, StreamJ, StreamK, StreamL, StreamMax };

      int _setupCompleted;

      // inverse map stream width

      int _inverseMapStreamWidth;

      // actual max number of inverse maps 

      int _maxNumberOfInverseMaps;

      // scale factors for 1-4 ixn

      BrookOpenMMFloat _ljScale;
      BrookOpenMMFloat _coulombFactor;

      // streams

      BrookFloatStreamInternal* _particleIndicesStream;
      BrookFloatStreamInternal* _bondedParameters[NumberOfParameterStreams];
      BrookFloatStreamInternal* _bondedForceStreams[NumberOfForceStreams];
      BrookFloatStreamInternal* _chargeStream;
      BrookFloatStreamInternal* _inverseStreamMaps[NumberOfForceStreams][MaxNumberOfInverseMaps];

      int _maxInverseMapStreamCount[NumberOfForceStreams];
      int _inverseMapStreamCount[NumberOfForceStreams];

      // helper methods in setup of parameters

      void flipQuartet( int ibonded, int *particles );
      int matchTorsion( int i, int j, int k, int l, int nbondeds, int *particles );
      int matchAngle( int i, int j, int k, int nbondeds, int *particles, int *flag );
      int matchBond( int i, int j, int nbondeds, int *particles, int *flag );
      int matchPair( int i, int j, int nbondeds, int *particles );

      /**
       * Setup Ryckaert-Bellemans parameters/particle indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param particles                     array of particle indices
       * @param params                    arrays of bond parameters
       * @param rbTorsionIndices          the four particles connected by each Ryckaert-Bellemans torsion term
       * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
       *
       * @return nonzero value if error
       *
       */

      int addRBDihedrals( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& rbTorsionIndices,  
                          const std::vector<std::vector<double> >& rbTorsionParameters );

      /**
       * Setup periodic torsion parameters/particle indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param particles                     array of particle indices
       * @param params                    arrays of bond parameters
       * @param periodicTorsionIndices    the four particles connected by each periodic torsion term
       * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
       *
       * @return nonzero value if error
       *
       */

      int addPDihedrals( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& periodicTorsionIndices,  
                          const std::vector<std::vector<double> >& periodicTorsionParameters );

      /**
       * Setup angle bond parameters/particle indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param particles                     array of particle indices
       * @param params                    arrays of bond parameters
       * @param angleIndices              the angle bond particle indices
       * @param angleParameters           the angle parameters (angle in radians, force constant)
       *
       * @return nonzero value if error
       *
       */

      int addAngles( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& angleIndices,
                     const std::vector<std::vector<double> >& angleParameters );

      /**
       * Setup harmonic bond parameters/particle indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param particles                     array of particle indices
       * @param params                    arrays of bond parameters
       * @param bondIndices               two harmonic bond particle indices
       * @param bondParameters            the force parameters (distance, k)
       *
       * @return nonzero value if error
       *
       */

      int addBonds( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& bondIndices,
                    const std::vector<std::vector<double> >& bondParameters );

      /**
       * Setup LJ/Coulomb 1-4 parameters/particle indices
       * 
       * @param nbondeds                  number of bonded entries
       * @param particles                     array of particle indices
       * @param params                    arrays of bond parameters
       * @param charges                   array of charges
       * @param bonded14Indices           each element contains the indices of two particles whose nonbonded interactions should be reduced since
       *                                  they form a bonded 1-4 pair
       * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each particle
       * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
       *
       * @return nonzero value if error
       *
       */

      int addPairs( int *nbondeds, int *particles, BrookOpenMMFloat* params[], BrookOpenMMFloat* charges,
                    const std::vector<std::vector<int> >& bonded14Indices, const std::vector<std::vector<double> >& nonbondedParameters,
                    double lj14Scale, double coulombScale );
      
      /**
       * Create and load inverse maps for bonded ixns
       * 
       * @param nbondeds                  number of bonded entries
       * @param nparticles                    number of particles
       * @param particles                     arrays of particle indices (particles[numberOfBonds][4])
       * @param platform                  BrookPlatform reference
       * @param log                       log file reference (optional)
       *
       * @return nonzero value if error
       *
       */

      //int loadInvMaps( int nbondeds, int nparticles, int *particles, const BrookPlatform& platform );
      int loadInvMaps( int nbondeds, int nparticles, int *particles, int particleStreamWidth, int particleStreamSize );
      
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
      
      /*
       * Helper functions for building inverse maps for 
       * torsions, impropers and angles.
       * 
       * For each particle, calculates the positions at which it's
       * forces are to be picked up from and stores the position
       * in the appropriate index.
       *
       * Input: number of dihedrals, the particle indices, and a flag indicating
       *        whether we're doing i(0), j(1), k(2) or l(3)
       * Output: an array of counts per particle
       *         arrays of inversemaps
       *         nimaps - the number of invmaps actually used.
       *
       * @param posflag       0-niparticles-1
       * @param niparticles       3 for angles, 4 for torsions, impropers
       * @param nints         number of interactions
       * @param nparticles        number of particles
       * @param *particles        gromacs interaction list
       * @param nmaps         maximum number of inverse maps
       * @param   counts[]    output counts of how many places each particle occurs
       * @param *invmaps[]    output array of nmaps inverse maps
       * @param *nimaps,      output max number of inverse maps actually used
       *
       * @return DefaultReturnValue, unless error in which case exits w/ OpenMM exception
       *
       **/
      
      int gpuCalcInvMap( int posflag, int niparticles, int nints, int nparticles,
                          int *particles, int nmaps, int counts[], float4 *invmaps[],
                          int *nimaps );
      
      void gpuPrintInvMaps( int nmaps, int nparticles, int counts[], float4 *invmap[], FILE* logFile );
      
      /* We are still plagued by kernel call overheads. This is for a big fat
       * merged inverse gather kernel:
       * Since we have 32 bit floats, we have 23 bits of mantissa or the largest
       * integer we can represent is 2^23. So it should be quite safe to add 
       * 100000 * n to the index where n is the stream in which we should do the
       * lookup. This assumes that nints < 100000, preferably nints << 100000
       * which should always be true
       * */
      int gpuCalcInvMap_merged( int nints, int nparticles, int *particles, int nmaps, int counts[], float4 *invmaps[], int *nimaps );
      
      /* Repacks the invmap streams for more efficient access in the
       * merged inverse gather kernel
       *
       * buf should be nimaps * nparticles large.
       * */
      int gpuRepackInvMap_merged( int nparticles, int nmaps, int *counts, float4 *invmaps[], float4 *buf );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_BONDED_H_ */
