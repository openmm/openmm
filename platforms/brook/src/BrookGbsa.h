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
#include <set>

#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class BrookGbsa : public BrookCommon {

   public:
  
      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      BrookGbsa( );
  
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
       * Get atom ceiling parameter
       * 
       * @return atom ceiling parameter
       *
       */
         
      int getAtomSizeCeiling( void ) const;
 
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
       * Get atomic radii
       *
       * @return   atomic radii stream
       *
       */
      
      BrookFloatStreamInternal* getObcAtomicRadii( void ) const;
      
      /** 
       * Get scaled atomic radii
       *
       * @return   scaled atomic radii stream
       *
       */
      
      BrookFloatStreamInternal* getObcScaledAtomicRadii( void ) const;
      
      /** 
       * Get atomic radii w/ dielectric offset
       *
       * @return   atomic radii w/ dielectric offset stream
       *
       */
      
      BrookFloatStreamInternal* getObcAtomicRadiiWithDielectricOffset( void ) const;
      
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
       * Return true if force[index] stream is set 
       *
       * @return  true  if index is valid && force[index] stream is set; else false
       *
       */
      
      int isForceStreamSet( int index ) const;
      
      /* 
       * Setup of Gbsa parameters
       *
       * @param atomParameters        vector of OBC parameters [atomI][0=index]
       *                                                       [atomI][1=charge]
       *                                                       [atomI][2=radius]
       *                                                       [atomI][2=scaling factor]
       * @param solventDielectric     solvent dielectric
       * @param soluteDielectric      solute dielectric
       *
       * @return nonzero value if error
       *
       * */
      
      int setup( const std::vector<std::vector<double> >& atomParameters, 
                 double solventDielectric, double soluteDielectric );
      
      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContentsString( int level = 0 ) const;

   private:
   
      // fixed number of force streams

      static const int NumberOfForceStreams     = 4;

      // streams indices

      enum { 
              ObcAtomicRadiiStream,
              ObcRadiiStreamStream,
              ObcScaledAtomicRadiiStream,
              ObcAtomicRadiiWithDielectricOffsetStream,
              ObcBornRadiiStream,
              ObcBornRadii2Stream,
              ObcIntermediateForceStream,
              ObcChainStream,
              LastStreamIndex
           };

      // atom ceiling

      int _atomSizeCeiling;

      // unroll in i/j dimensions

      int _outerUnroll;
      int _innerUnroll;

      // duplication factor

      int _duplicationFactor;

      // force stream width

      int _partialForceStreamWidth;
      int _partialForceStreamHeight;
      int _partialForceStreamSize;

      // Gbsa stream dimensions

      int _gbsaStreamWidth;
      int _gbsaStreamHeight;
      int _gbsaStreamSize;

      // dielectrics

      double _solventDielectric
      double _soluteDielectric

      // internal streams

      BrookFloatStreamInternal* _gbsaStreams[LastStreamIndex];
      BrookFloatStreamInternal* _gbsaForceStreams[NumberOfForceStreams];

      /* 
       * Setup of stream dimensions
       *
       * @param atomStreamSize        atom stream size
       * @param atomStreamWidth       atom stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int initializeStreamSizes( int atomStreamSize, int atomStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfAtoms             number of atoms
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       */
      
      int initializeStreamSizes( int numberOfAtoms, const Platform& platform );
      
      /** 
       * Initialize stream dimensions and streams
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int initializeStreams( const Platform& platform );
      
      /* 
       * Setup of stream dimensions for partial force streams
       *
       * @param atomStreamSize        atom stream size
       * @param atomStreamWidth       atom stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int initializePartialForceStreamSize( int atomStreamSize, int atomStreamWidth );
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_GBSA_H_ */
