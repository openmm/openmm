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
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include <vector>

#include "BrookStreamImpl.h"
#include "BrookCommon.h"
#include "openmm/Context.h"
#include "BrookBondParameters.h"

namespace OpenMM {

/**
 * Calculate Brook bonded forces
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
       * @param numberOfParticles                   numberOfParticles
       * @param harmonicBondBrookBondParameters     force parameters (length, k) for each bond term
       * @param harmonicAngleBrookBondParameters    force parameters (angle, k) for each angle term
       * @param periodicTorsionBrookBondParameters  force parameters (k, phase, periodicity) for each periodic torsion term
       * @param rbTorsionBrookBondParameters        coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
       * @param nonBonded14ForceParameters          nonbonded 14 force parameters (charge, sigma, epsilon) for each particle
       * @param particleStreamWidth                 stream width
       * @param particleStreamSize                  stream size
       *
       * @return nonzero value if error
       *
       */

      int setup( int numberOfParticles,
                 BrookBondParameters* harmonicBondBrookBondParameters,
                 BrookBondParameters* harmonicAngleBrookBondParameters,
                 BrookBondParameters* periodicTorsionBrookBondParameters,
                 BrookBondParameters* rbTorsionBrookBondParameters,
                 BrookBondParameters* nonBonded14ForceParameters,  int particleStreamWidth, int particleStreamSize );

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
   
      /** 
       * Set SetupCompleted flag
       *
       * @param SetupCompleted flag
       *
       * @return SetupCompleted flag
       */
      
      int setupCompleted( int isSetupCompleted );
   
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

      void _flipQuartet( int ibonded, int *particles );

      int _matchTorsion( int i, int j, int k, int l, int nbondeds, int *particles );
      int _matchAngle(   int i, int j, int k, int nbondeds, int *particles, int *flag );
      int _matchBond(    int i, int j, int nbondeds, int *particles, int *flag );
      int _matchPair(    int i, int j, int nbondeds, int *particles );

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

      int _addRBTorsions( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& rbTorsionIndices,  
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

      int _addPTorsions( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& periodicTorsionIndices,  
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

      int _addAngles( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& angleIndices,
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

      int  _addBonds( int *nbondeds, int *particles, BrookOpenMMFloat* params[], const std::vector<std::vector<int> >& bondIndices,
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
       *
       * @return nonzero value if error
       *
       */

      int _addPairs( int *nbondeds, int *particles, BrookOpenMMFloat* params[], BrookOpenMMFloat* charges,
                     const std::vector<std::vector<int> >& bonded14Indices, const std::vector<std::vector<double> >& nonbondedParameters );
      
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
      int _loadInvMaps( int nbondeds, int nparticles, int *particles, int particleStreamWidth, int particleStreamSize );
      
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

      int _validateInverseMapStreamCount( int index, int count ) const;
      
      /*
       * Helper functions for building inverse maps for 
       * torsions, impropers and angles.
       * 
       * For each particle, calculates the positions at which it's
       * forces are to be picked up from and stores the position
       * in the appropriate index.
       *
       * Input: number of torsions, the particle indices, and a flag indicating
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
      
      int _gpuCalcInvMap( int posflag, int niparticles, int nints, int nparticles,
                          int *particles, int nmaps, int counts[], float4 *invmaps[], int *nimaps );
      
      void _gpuPrintInvMaps( int nmaps, int nparticles, int counts[], float4 *invmap[], FILE* logFile );
      
      /* 
       * We are still plagued by kernel call overheads. This is for a big fat
       * merged inverse gather kernel:
       * Since we have 32 bit floats, we have 23 bits of mantissa or the largest
       * integer we can represent is 2^23. So it should be quite safe to add 
       * 100000 * n to the index where n is the stream in which we should do the
       * lookup. This assumes that nints < 100000, preferably nints << 100000
       * which should always be true
       *
       **/
      int _gpuCalcInvMap_merged( int nints, int nparticles, int *particles, int nmaps, int counts[], float4 *invmaps[], int *nimaps );
      
      /* 
       * Repacks the invmap streams for more efficient access in the
       * merged inverse gather kernel
       *
       * buf should be nimaps * nparticles large.
       *
       **/
      int _gpuRepackInvMap_merged( int nparticles, int nmaps, int *counts, float4 *invmaps[], float4 *buf );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_BONDED_H_ */
