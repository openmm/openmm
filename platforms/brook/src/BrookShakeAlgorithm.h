#ifndef OPENMM_BROOK_SHAKE_ALGORITHM_H_
#define OPENMM_BROOK_SHAKE_ALGORITHM_H_

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
#include <set>
#include <sstream>

#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"
#include "openmm/OpenMMException.h"

namespace OpenMM {

struct ShakeCluster {

   int _centralID;
   int _peripheralID[3];
   int _size;
   float _distance;
   float _centralInvMass, _peripheralInvMass;

   ShakeCluster(){};
   ShakeCluster( int centralID, float invMass) : _centralID(centralID), _centralInvMass(invMass), _size(0) {
      _peripheralID[0] = _peripheralID[1] = _peripheralID[2] = -1;
   }; 

   void addAtom( int id, float dist, float invMass){
      if( _size == 3 ){
          std::stringstream message;
          message << "ShakeCluster::addAtom: " << "atom " << id << " has more than 3 constraints!." << std::endl; 
          throw OpenMMException( message.str() );
      }   
      if( _size > 0 && dist != _distance ){
          std::stringstream message;
          message << "ShakeCluster::addAtom: " << "atom " << id << " has different constraint distances: " <<  dist << " and " << _distance << std::endl; 
          throw OpenMMException( message.str() );
      }   
      if( _size > 0 && invMass != _peripheralInvMass ){
          std::stringstream message;
          message << "ShakeCluster::addAtom: " << " constrainted atoms associated w/ atom " << id << " have different masses: " <<  invMass << " and " << _peripheralInvMass << std::endl; 
          throw OpenMMException( message.str() );
      }   
      _peripheralID[_size++]  = id; 
      _distance              = dist;
      _peripheralInvMass     = invMass;
   }   
};

/**
 *
 * Encapsulates stochastic dynamics algorithm 
 *
 */

class BrookShakeAlgorithm : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookShakeAlgorithm( );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookShakeAlgorithm();
  
      /** 
       * Get number of constraints
       * 
       * @return   number of constraints
       *
       */
      
      int getNumberOfConstraints( void ) const;
      
      /** 
       * Get max iterations
       * 
       * @return   max iterations
       *
       */
      
      int getMaxIterations( void ) const;
      
      /** 
       * Set  max iterations
       * 
       * @param   max iterations
       *
       * @return DefaultReturnValue
       *
       */
      
      int setMaxIterations( int maxIterations );
      
      /** 
       * Get SHAKE tolerance
       * 
       * @return  SHAKE tolerance  
       *
       */
      
      BrookOpenMMFloat getShakeTolerance( void ) const;
      
      /** 
       * Set SHAKE tolerance
       * 
       * @param  SHAKE tolerance  
       *
       * @return DefaultReturnValue
       *
       */
      
      int setShakeTolerance( BrookOpenMMFloat tolerance );
      
      /**
       * Get Shake particle stream width
       *
       * @return particle stream width
       */

      int getShakeParticleStreamWidth( void ) const; 

      /**
       * Get Shake particle stream height
       *
       * @return particle stream height
       */

      int getShakeParticleStreamHeight( void ) const;

      /**
       * Get Shake particle stream size
       * 
       * @return particle stream size
       */

      int getShakeParticleStreamSize( void ) const; 

      /**
       * Get Shake constraint stream width
       *
       * @return constraint stream width
       */

      int getShakeConstraintStreamWidth( void ) const; 

      /**
       * Get Shake constraint stream height
       *
       * @return constraint stream height
       */

      int getShakeConstraintStreamHeight( void ) const;

      /**
       * Get Shake constraint stream size
       * 
       * @return constraint stream size
       */

      int getShakeConstraintStreamSize( void ) const; 

      /** 
       * Get array of Shake streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /*  
       * Setup of Shake parameters
       *
       * @param masses                masses
       * @param constraintIndices     constraint particle indices
       * @param constraintLengths     constraint lengths
       * @param platform              Brook platform
       *
       * @return ErrorReturnValue if error
       *
       */
              
      int setup( const std::vector<double>& masses, const std::vector< std::vector<int> >& constraintIndices,
                 const std::vector<double>& constraintLengths, const Platform& platform );
          
      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContentsString( int level = 0 ) const;

      /** 
       * Get Shake particle indices stream
       *
       * @return  Shake particle indices stream
       *
       */
      
      BrookFloatStreamInternal* getShakeParticleIndicesStream( void ) const;
      
      /** 
       * Get  Shake particle parameter stream
       *
       * @return   Shake particle parameter stream
       *
       */
      
      BrookFloatStreamInternal* getShakeParticleParameterStream( void ) const;
      
      /** 
       * Get XCons0 stream
       *
       * @return  XCons0 stream
       *
       */
      
      BrookFloatStreamInternal* getShakeXCons0Stream( void ) const;
      
      /** 
       * Get XCons1 stream
       *
       * @return  XCons1 stream
       *
       */
      
      BrookFloatStreamInternal* getShakeXCons1Stream( void ) const;
      
      /** 
       * Get XCons2 stream
       *
       * @return  XCons2 stream
       *
       */
      
      BrookFloatStreamInternal* getShakeXCons2Stream( void ) const;
      
      /** 
       * Get XCons3 stream
       *
       * @return  XCons3 stream
       *
       */
      
      BrookFloatStreamInternal* getShakeXCons3Stream( void ) const;
      
      /** 
       * Get Shake inverse map stream
       *
       * @return  Shake inverse map stream
       *
       */
      
      BrookFloatStreamInternal* getShakeInverseMapStream( void ) const;

      /*  
       * Check constraints
       *
       * @param positions             atom positions
       * @param outputString          output message
       * @param tolerance             tolerance to compare (if < 0, then use algorithm tolerance
       *
       * @return number of errors
       *
       */
              
      int checkConstraints( BrookStreamInternal* positions, std::string& outputString, float tolerance ) const;
              
   private:
   
      // streams indices

      enum BrookShakeAlgorithmStreams { 
              ShakeParticleIndicesStream,
              ShakeParticleParameterStream,
              ShakeXCons0Stream,
              ShakeXCons1Stream,
              ShakeXCons2Stream,
              ShakeXCons3Stream,
              ShakeInverseMapStream,
              LastStreamIndex
           };

      // number of constraints

      int _numberOfConstraints;

      // max iterations

      int _maxIterations;

      // particle stream dimensions

      int _shakeParticleStreamWidth;
      int _shakeParticleStreamHeight;
      int _shakeParticleStreamSize;

      // constraint stream dimensions

      int _shakeConstraintStreamSize;
      int _shakeConstraintStreamWidth;
      int _shakeConstraintStreamHeight;

      // SHAKE tolerance

      BrookOpenMMFloat _shakeTolerance;

      // inverse sqrt masses

      BrookOpenMMFloat* _inverseSqrtMasses;

      // internal streams

      BrookFloatStreamInternal* _shakeStreams[LastStreamIndex];

      // Shake cluster (used for debugging)

      std::map<int, ShakeCluster> _clusters;

      /* 
       * Setup of stream dimensions
       *
       * @param particleStreamSize        particle stream size
       * @param particleStreamWidth       particle stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValueValue
       *
       * */
      
      int _initializeStreamSizes( int particleStreamSize, int particleStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfParticles         number of particles
       * @param numberOfConstraints       number of constraints
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValueValue
       *
       */
      
      int _initializeStreamSizes( int numberOfParticles, int numberOfConstraints, const Platform& platform );
      
      /** 
       * Initialize stream dimensions and streams
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int _initializeStreams( const Platform& platform );

      /*  
       * Set Shake streams
       *
       * @param masses                masses
       * @param constraintIndices     constraint particle indices
       * @param constraintLengths     constraint lengths
       * @param platform              platform reference
       *
       * @return ErrorReturnValue if error
       *
       * @throw OpenMMException if constraintIndices.size() != constraintLengths.size()
       *
       */
           
      int _setShakeStreams( const std::vector<double>& masses, const std::vector< std::vector<int> >& constraintIndices,
                            const std::vector<double>& constraintLengths, const Platform& platform );
          
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_SHAKE_ALGORITHM_H_ */
