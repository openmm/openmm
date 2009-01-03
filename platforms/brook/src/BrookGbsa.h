#ifndef OPENMM_BROOK_GBSA_H_
#define OPENMM_BROOK_GBSA_H_

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
#include "BrookPlatform.h"
#include "BrookCommon.h"
#include "../../../platforms/reference/src/gbsa/CpuObc.h"

namespace OpenMM {

/**
 *
 * Used by BrookCalcGBSAOBCForceKernel kernel to execute OBC algorithm on GPU
 *
 */

class BrookGbsa : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookGbsa( );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookGbsa();
  
      /**
       * Return number of force streams
       * 
       * @return number of force streams
       *
       */

      int getNumberOfForceStreams( void ) const; 

      /** 
       * Get duplication factor
       * 
       * @return   duplication factor
       *
       */
      
      int getDuplicationFactor( void ) const;
      
      /** 
       * Get particle ceiling parameter
       * 
       * @return particle ceiling parameter
       *
       */
         
      int getParticleSizeCeiling( void ) const;
 
      /** 
       * Get outer loop unroll
       * 
       * @return   outer loop unroll (fixed value)
       *
       */
      
      int getOuterLoopUnroll( void ) const;
      
      /** 
       * Set outer loop unroll
       * 
       * @param  outer loop unroll (fixed value)
       *
       * @return updated outer loop unroll (fixed value)
       *
       */
      
      int setOuterLoopUnroll( int outerUnroll );
      
      /**
       * Return unrolling for inner loops
       * 
       * @return outer loop unrolling
       */

      int getInnerLoopUnroll( void ) const; 

      /**
       * Return true if ACE approximation is to be included
       * 
       * @return true if ACE approximation is to be included
       */

      int includeAce( void ) const; 

      /**
       * Get partial force stream width
       *
       * @return partial force stream width
       */

      int getPartialForceStreamWidth( void ) const; 

      /**
       * Get partial force stream height
       *
       * @return partial force stream height
       */

      int getPartialForceStreamHeight( void ) const;

      /**
       * Get partial force stream size
       * 
       * @return partial force stream size
       */

      int getPartialForceStreamSize( void ) const; 

      /**
       * Get Gbsa particle stream width
       *
       * @return particle stream width
       */

      int getGbsaParticleStreamWidth( void ) const; 

      /**
       * Get Gbsa particle stream height
       *
       * @return particle stream height
       */

      int getGbsaParticleStreamHeight( void ) const;

      /**
       * Get Gbsa particle stream size
       * 
       * @return particle stream size
       */

      int getGbsaParticleStreamSize( void ) const; 

      /**
       * Get solute dielectric
       * 
       * @return solute dielectric
       */

      float getSoluteDielectric( void ) const; 

      /**
       * Get solvent dielectric
       * 
       * @return solvent dielectric
       */

      float getSolventDielectric( void ) const; 

      /**
       * Get OBC dielectric offset
       * 
       * @return dielectric offset
       */

      float getDielectricOffset( void ) const; 

      /** 
       * Get particle radii
       *
       * @return   particle radii stream
       *
       */
      
      BrookFloatStreamInternal* getObcParticleRadii( void ) const;
      
      /** 
       * Get scaled particle radii
       *
       * @return   scaled particle radii stream
       *
       */
      
      BrookFloatStreamInternal* getObcScaledParticleRadii( void ) const;
      
      /** 
       * Get particle radii w/ dielectric offset
       *
       * @return   particle radii w/ dielectric offset stream
       *
       */
      
      BrookFloatStreamInternal* getObcParticleRadiiWithDielectricOffset( void ) const;
      
      /** 
       * Get Born radii stream 
       *
       * @return  Born radii stream
       *
       */
      
      BrookFloatStreamInternal* getObcBornRadii( void ) const;
      
      /** 
       * Get Born radii2 stream
       *
       * @return   Born radii2 stream
       *
       */
      
      BrookFloatStreamInternal* getObcBornRadii2( void ) const;
      
      /** 
       * Get Obc intermediate force stream 
       *
       * @return  Obc intermediate force stream
       *
       */
      
      BrookFloatStreamInternal* getObcIntermediateForce( void ) const;
      
      /** 
       * Get Obc chain stream
       *
       * @return   Obc chain  stream
       *
       */
      
      BrookFloatStreamInternal* getObcChain( void ) const;
      
      /** 
       * Get force streams 
       *
       * @return  force streams
       *
       */
      
      BrookFloatStreamInternal** getForceStreams( void );
      
      /** 
       * Get array of Gbsa streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /** 
       * Return true if force[index] stream is set 
       *
       * @return  true  if index is valid && force[index] stream is set; else false
       *
       */
      
      int isForceStreamSet( int index ) const;
      
      /** 
       * Return true if Born radii have been initialized
       *
       * @return  true  if Born radii have been initialized
       *
       */
      
      int haveBornRadiiBeenInitialized( void ) const;
      
      /** 
       * Calculate Born radii
       *
       * @return  calculate Born radii
       *
       */
      
      // int calculateBornRadii( const Stream& positions );
      
      /* 
       * Setup of Gbsa parameters
       *
       * @param particleParameters    vector of OBC parameters [particleI][0=charge]
       *                                                       [particleI][1=radius]
       *                                                       [particleI][2=scaling factor]
       * @param solventDielectric     solvent dielectric
       * @param soluteDielectric      solute dielectric
       * @param platform              Brook platform
       *
       * @return nonzero value if error
       *
       * */
      
      int setup( const std::vector<std::vector<double> >& particleParameters, 
                 double solventDielectric, double soluteDielectric, const Platform& platform  );
      
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
       * Compute forces
       * 
       */
      
      void computeForces( BrookStreamImpl& positionStream, BrookStreamImpl& forceStream );
      
   private:
   
      // fixed number of force streams

      static const int NumberOfForceStreams     = 4;

      // streams indices

      enum { 
              ObcParticleRadiiStream,
              ObcScaledParticleRadiiStream,
              ObcParticleRadiiWithDielectricOffsetStream,
              ObcBornRadiiStream,
              ObcBornRadii2Stream,
              ObcIntermediateForceStream,
              ObcChainStream,
              LastStreamIndex
           };

      // particle ceiling

      int _particleSizeCeiling;

      // unroll in i/j dimensions

      int _outerUnroll;
      int _innerUnroll;

      // duplication factor

      int _duplicationFactor;

      // include ACE approximation

      int _includeAce;

      // force stream width

      int _partialForceStreamWidth;
      int _partialForceStreamHeight;
      int _partialForceStreamSize;

      // Particle stream dimensions

      int _gbsaParticleStreamWidth;
      int _gbsaParticleStreamHeight;
      int _gbsaParticleStreamSize;

      // dielectrics

      double _solventDielectric;
      double _soluteDielectric;

      // dielectric offset

      double _dielectricOffset;

      // internal streams

      BrookFloatStreamInternal* _gbsaStreams[LastStreamIndex];
      BrookFloatStreamInternal* _gbsaForceStreams[NumberOfForceStreams];

      int _bornRadiiInitialized;

      // CpuObc reference -- was used to calculate initial Born radii
      // no longer used -- to be removed?

      // CpuObc* _cpuObc;

      /* 
       * Setup of stream dimensions
       *
       * @param particleStreamSize        particle stream size
       * @param particleStreamWidth       particle stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int _initializeStreamSizes( int particleStreamSize, int particleStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfParticles         number of particles
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       */
      
      int _initializeStreamSizes( int numberOfParticles, const Platform& platform );
      
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
       * Setup of stream dimensions for partial force streams
       *
       * @param particleStreamSize        particle stream size
       * @param particleStreamWidth       particle stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int _initializePartialForceStreamSize( int particleStreamSize, int particleStreamWidth );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_GBSA_H_ */
