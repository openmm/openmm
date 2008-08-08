#ifndef BrookNonBonded_H_
#define BrookNonBonded_H_

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

#include "BrookFloatStreamImpl.h"
#include "BrookIntStreamImpl.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 * This kernel is invoked by StandardMMForceField to calculate the forces acting on the system.
 */
class BrookNonBonded : public BrookCommon {

   public:
  
      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      BrookNonBonded( );
  
      ~BrookNonBonded();
  
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
       * Get j-stream width
       *
       * @param platform platform reference
       *
       * @return j-stream width
       */

      int getJStreamWidth( const Platform& platform ); 

      /**
       * Get j-stream width
       *
       * @return j-stream width
       */

      int getJStreamWidth( void ) const; 

      /**
       * Get j-stream height
       *
       * @param platform platform reference
       *
       * @return j-stream height
       */

      int getJStreamHeight( const Platform& platform ); 

      /**
       * Get j-stream height
       *
       * @return j-stream height
       */

      int getJStreamHeight( void ) const;

      /**
       * Get j-stream size
       * 
       * @param platform platform reference
       *
       * @return j-stream size
       */

      int getJStreamSize( const Platform& platform ); 

      /**
       * Get j-stream size
       * 
       * @return j-stream size
       */

      int getJStreamSize( void ) const; 

      /**
       * Get partial force stream width
       *
       * @param platform platform reference
       *
       * @return partial force stream width
       */

      int getPartialForceStreamWidth( const Platform& platform ); 

      /**
       * Get partial force stream width
       *
       * @return partial force stream width
       */

      int getPartialForceStreamWidth( void ) const; 

      /**
       * Get partial force stream height
       *
       * @param platform platform reference
       *
       * @return partial force stream height
       */

      int getPartialForceStreamHeight( const Platform& platform ); 

      /**
       * Get partial force stream height
       *
       * @return partial force stream height
       */

      int getPartialForceStreamHeight( void ) const;

      /**
       * Get partial force stream size
       * 
       * @param platform platform reference
       *
       * @return partial force stream size
       */

      int getPartialForceStreamSize( const Platform& platform ); 

      /**
       * Get partial force stream size
       * 
       * @return partial force stream size
       */

      int getPartialForceStreamSize( void ) const; 

      /**
       * Get partial force stream size
       * 
       * @return partial force stream size
       */

      /**
       * Get exclusion stream width
       * 
       * @return exclusion stream width
       */

      int getExclusionStreamWidth( void ) const; 

      /**
       * Get exclusion stream size
       * 
       * @return exclusion stream size
       */

      int getExclusionStreamSize( void ) const; 

      /** 
       * Get exclusion stream 
       *
       * @return  exclusion stream
       *
       */
      
      BrookFloatStreamImpl* getExclusionStream( void ) const;
      
      /** 
       * Get vdw stream 
       *
       * @return  vdw stream
       *
       */
      
      BrookFloatStreamImpl* getVdwStream( void ) const;
      
      /** 
       * Get charge stream 
       *
       * @return  charge stream
       *
       */
      
      BrookFloatStreamImpl* getChargeStream( void ) const;
      
      /** 
       * Get sigma-eps stream 
       *
       * @return  sigma-eps stream
       *
       */
      
      BrookFloatStreamImpl* getSigmaStream( void ) const;
      
      /** 
       * Get epsilon stream 
       *
       * @return  epsilon stream
       *
       */
      
      BrookFloatStreamImpl* getEpsilonStream( void ) const;
      
      /** 
       * Get force streams 
       *
       * @return  force streams
       *
       */
      
      BrookFloatStreamImpl** getForceStreams( void );
      
      /** 
       * Return true if force[index] stream is set 
       *
       * @return  true  if index is valid && force[index] stream is set; else false
       *
       */
      
      int isForceStreamSet( int index ) const;
      
      /* 
       * Setup of nonbonded ixns
       *
       * @param numberOfAtoms         number of atoms
       * @param nonbondedParameters   vector of nonbonded parameters [atomI][0=c6]
       *                                                             [atomI][1=c12]
       *                                                             [atomI][2=charge]
       * @param platform              Brook platform
       * @param log                   optional Log file reference
       *
       * @return nonzero value if error
       * */
      
      int setup( int numberOfAtoms, const std::vector<std::vector<double> >& nonbondedParameters,
                 const std::vector<std::set<int> >& exclusions,  const BrookPlatform& platform );
      
      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

   private:
   
      // fixed number of force streams

      static const int NumberOfForceStreams     = 4;

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

      // exclusions stream dimensions

      int _exclusionStreamWidth;
      int _exclusionStreamHeight;
      int _exclusionStreamSize;

      // j-stream dimensions

      int _jStreamWidth;
      int _jStreamHeight;
      int _jStreamSize;

      // streams

      BrookFloatStreamImpl* _exclusionStream;
      BrookFloatStreamImpl* _vdwStream;
      BrookFloatStreamImpl* _chargeStream;
      BrookFloatStreamImpl* _sigmaStream;
      BrookFloatStreamImpl* _epsilonStream;
      BrookFloatStreamImpl* _nonbondedForceStreams[NumberOfForceStreams];

      /** 
       * Initialize exclusion stream dimensions and stream
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int initializeExclusionStream( const Platform& platform );
            
      /** 
       * Set exclusion (4x4)
       * 
       * @param i                         atom i index
       * @param j                         atom j index
       * @param exclusionStreamWidth      exclusion stream width
       * @param exclusion                 array of packed exclusions
       *
       * @return nonzero value if error
       *
       */
      
      int setExclusion( int i, int j, int exclusionStreamWidth, BrookOpenMMFloat* exclusion );
      
      /** 
       * Initialize exclusions
       * 
       * @param exclusions                vector of sets containing exclusions (1 set entry for every atom)
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int initializeExclusions( const std::vector<std::set<int> >& exclusionsVector, const Platform& platform );

      /** 
       * Initialize stream dimensions and streams
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int initializeStreams( int numberOfAtoms, const Platform& platform );
      
      /** 
       * Set sigma & epsilon given c6 & c12 (geometric rule)
       * 
       * @param c6                        vdw c6
       * @param c12                       vdw c12
       * @param sigma                     massaged sigma
       * @param epsilon                   massaged epsilon
       *
       * @return nonzero value if error
       *
       */
      
      int setSigmaEpsilon( double c6, double c12, double* sigma , double* epsilon );

      /** 
       * Initialize vdw & charge
       * 
       * @param exclusions                vector of sets containing exclusions (1 set entry for every atom)
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int initializeVdwAndCharge( const std::vector<std::vector<double> >& nonbondedParameters, const Platform& platform );
      
};

} // namespace OpenMM

#endif /*OPENMM_BROOKKERNELS_H_*/
