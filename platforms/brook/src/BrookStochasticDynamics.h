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
  
      /** 
       * Constructor
       * 
       * @param masses             atomic masses
       * @param randomNumberSeed   random number seed
       *
       */
      
      BrookStochasticDynamics( const std::vector<double>& masses, uint32_t randomNumberSeed = 1364 );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookStochasticDynamics();
  
      /**
       * Get tau
       *
       * @return tau
       */

      BrookOpenMMFloat getTau( void ) const; 

      /**
       * Get temperature
       *
       * @return temperature
       */

      BrookOpenMMFloat getTemperature( void ) const; 

      /**
       * Get step size
       *
       * @return step size
       */

      BrookOpenMMFloat getStepSize( void ) const; 

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

      /** 
       * Get SDPC1 stream 
       *
       * @return  SDPC1 stream
       *
       */
      
      BrookFloatStreamInternal* getSDPC1( void ) const;
      
      /** 
       * Get SDPC2 stream 
       *
       * @return  SDPC2 stream
       *
       */
      
      BrookFloatStreamInternal* getSDPC2( void ) const;
      
      /** 
       * Get SD2X stream 
       *
       * @return  SD2X stream
       *
       */
      
      BrookFloatStreamInternal* getSD2X( void ) const;
      
      /** 
       * Get SD1V stream 
       *
       * @return  SD1V stream
       *
       */
      
      BrookFloatStreamInternal* getSD1V( void ) const;

   private:
   
      enum DerivedParameters   { GDT, EPH, EMH, EP, EM, B, C, D, V, X, Yv, Yx,
                                 Sd1pc1, Sd1pc2, Sd1pc3, Sd2pc1, Sd2pc2, MaxDerivedParameters };

      BrookOpenMMFloat _derivedParameters[MaxDerivedParameters];

      // streams indices

      enum { 
              SDPC1Stream,
              SDPC2Stream,
              SD2XStream,
              SD1VStream,
              LastStreamIndex
           };

      // randomNumberSeed

      uint32_t _randomNumberSeed;

      BrookOpenMMFloat _tau;
      BrookOpenMMFloat _temperature;
      BrookOpenMMFloat _stepSize;

      // Atom stream dimensions

      int _sdAtomStreamWidth;
      int _sdAtomStreamHeight;
      int _sdAtomStreamSize;

      /** 
       * Get derived parameter string
       * 
       * @return  string
       *
       */

      std::string _getDerivedParametersString( BrookStochasticDynamics::DerivedParameters ) const;

      /** 
       * Update derived parameters
       * 
       * @return  DefaultReturn
       *
       */

      int _updateDerivedParameters( void );

      /*
       * Update streams
       * 
       * @return  DefaultReturn
       *
       */
      
      int _updateSdStreams( void );
      
      // inverse sqrt masses

      BrookOpenMMFloat* _inverseSqrtMasses;

      // internal streams

      BrookFloatStreamInternal* _sdStreams[LastStreamIndex];

      /** 
       * Set tau
       * 
       * @param tau   new tau value
       *
       * @return      DefaultReturn
       *
       */
      
      int _setTau( BrookOpenMMFloat tau );
      
      /** 
       * Set friction = 1/tau
       * 
       * @param friction   new friction value
       *
       * @return      DefaultReturn
       *
       */
      
      int _setFriction( BrookOpenMMFloat friction );
      
      /** 
       * Set temperature
       * 
       * @parameter   temperature
       *
       * @return      DefaultReturn
       *
       */
      
      int _setTemperature( BrookOpenMMFloat temperature );
      
      /** 
       * Set stepSize
       * 
       * @param   stepSize
       *
       * @return      DefaultReturn
       *
       */
      
      int _setStepSize( BrookOpenMMFloat stepSize );
      
      /* 
       * Setup of stream dimensions
       *
       * @param atomStreamSize        atom stream size
       * @param atomStreamWidth       atom stream width
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       * */
      
      int _initializeStreamSizes( int atomStreamSize, int atomStreamWidth );

      /** 
       * Initialize stream dimensions
       * 
       * @param numberOfAtoms             number of atoms
       * @param platform                  platform
       *
       * @return ErrorReturnValue if error, else DefaultReturnValue
       *
       */
      
      int _initializeStreamSizes( int numberOfAtoms, const Platform& platform );
      
      /** 
       * Initialize stream dimensions and streams
       * 
       * @param platform                  platform
       *
       * @return nonzero value if error
       *
       */
      
      int _initializeStreams( const Platform& platform );

      /** 
       * Set masses 
       * 
       * @param masses             atomic masses
       *
       */
      
      int _setInverseSqrtMasses( const std::vector<double>& masses );
      
      
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_STOCHASTIC_DYNAMCIS_H_ */
