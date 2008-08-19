#ifndef OPENMM_BROOK_STOCHASTIC_DYNAMCIS_H_
#define OPENMM_BROOK_STOCHASTIC_DYNAMCIS_H_

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

class BrookStochasticDynamics : public BrookCommon {

   public:
  
      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1;

      /** 
       * Constructor
       * 
       */
      
      BrookStochasticDynamics( );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookStochasticDynamics();
  
      /**
       * Get StochasticDynamics atom stream width
       *
       * @return atom stream width
       */

      int getStochasticDynamicsAtomStreamWidth( void ) const; 

      /**
       * Get StochasticDynamics atom stream height
       *
       * @return atom stream height
       */

      int getStochasticDynamicsAtomStreamHeight( void ) const;

      /**
       * Get StochasticDynamics atom stream size
       * 
       * @return atom stream size
       */

      int getStochasticDynamicsAtomStreamSize( void ) const; 

      /** 
       * Update parameters
       * 
       * @param  temperature     temperature
       * @param  friction        friction
       *
       * @return   solute dielectric
       *
       */
      
      int updateParameters( double temperature, double friction );
      
      /** 
       * Get array of StochasticDynamics streams 
       *
       * @return  array ofstreams
       *
       */
      
      BrookFloatStreamInternal** getStreams( void );
      
      /* 
       * Setup of StochasticDynamics parameters
       *
       * @param atomParameters        vector of OBC parameters [atomI][0=charge]
       *                                                       [atomI][1=radius]
       *                                                       [atomI][2=scaling factor]
       * @param solventDielectric     solvent dielectric
       * @param soluteDielectric      solute dielectric
       * @param platform              Brook platform
       *
       * @return nonzero value if error
       *
       * */
      
      int setup( const std::vector<std::vector<double> >& atomParameters, 
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

   private:
   
      enum FixedParameters   { GDT, EPH, EMH, EP, EM, B, C, D, V, X, Yv, Yx, MaxFixedParameters };
      enum TwoDArrayIndicies { X2D, V2D, OldV, xPrime2D, vPrime2D, Max2DArrays };
      enum OneDArrayIndicies { InverseMasses, Max1DArrays };

      RealOpenMM _fixedParameters[MaxFixedParameters];

      // streams indices

      enum { 
              SDPC1Stream,
              SDPC2Stream,
              SD2XStream,
              SD1VStream,
              LastStreamIndex
           };

      // Atom stream dimensions

      int _sdAtomStreamWidth;
      int _sdAtomStreamHeight;
      int _sdAtomStreamSize;

      // inverse masses

      vector<BrookOpenMMFloat>& _inverseMasses;

      // internal streams

      BrookFloatStreamInternal* _sdStreams[LastStreamIndex];

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
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_STOCHASTIC_DYNAMCIS_H_ */
