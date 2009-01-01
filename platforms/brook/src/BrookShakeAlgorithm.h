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
#include <set>

#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

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
