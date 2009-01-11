#ifndef OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_
#define OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_

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

#include "BrookStreamImpl.h"
//#include "BrookFloatStreamInternal.h"
#include "BrookPlatform.h"
#include "BrookCommon.h"

namespace OpenMM {

/**
 *
 * Encapsulates removal of center of mass
 *
 */

class BrookVelocityCenterOfMassRemoval : public BrookCommon {

   public:
  
      /** 
       * Constructor
       * 
       */
      
      BrookVelocityCenterOfMassRemoval( );
  
      /** 
       * Destructor
       * 
       */
      
      ~BrookVelocityCenterOfMassRemoval( );
  
      /**
       * Get  particle stream width
       *
       * @return particle stream width
       */

      int getComParticleStreamWidth( void ) const; 

      /**
       * Get  particle stream height
       *
       * @return particle stream height
       */

      int getComParticleStreamHeight( void ) const;

      /**
       * Get  particle stream size
       * 
       * @return particle stream size
       */

      int getComParticleStreamSize( void ) const; 

      /** 
       * Remove velocity center-of-mass
       * 
       * @param  velocities                  particle velocities
       *
       * @return  DefaultReturnValue
       *
       */
      
      int removeVelocityCenterOfMass( BrookStreamImpl& velocityStream );

      /** 
       * Get velocity center-of-mass and kinetic energy (used for diagnostics)
       * 
       * @param  velocities                  particle velocities
       * @param  velocityCom                 output velocity com
       * @param  ke                          output kinetic energy
       *
       * @return  DefaultReturnValue
       *
       */
      
      int getVelocityCenterOfMass( BrookStreamImpl& vStream, BrookOpenMMFloat velocityCom[3], BrookOpenMMFloat* ke );

      /* 
       * Setup of  parameters
       *
       * @param masses                particle masses
       * @param platform              Brook platform
       *
       * @return ErrorReturnValue value if error, else DefaultReturnValue
       *
       * */
      
      int setup( const std::vector<double>& masses, const Platform& platform  );
      
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
       * Get inverse of total mass
       *
       * @return inverse of total mass
       *
       */
         
      BrookOpenMMFloat getTotalInverseMass( void ) const;
      
      /** 
       * Get inverse mass stream 
       *
       * @return  inverse mass stream
       *
       */
      
      BrookFloatStreamInternal* getMassStream( void ) const;
      
      /** 
       * Get work stream 
       *
       * @return  work stream
       *
       */
      
      BrookFloatStreamInternal* getWorkStream( void ) const;
      
      /** 
       * Get linear momentum stream 
       *
       * @return  linear momentum stream
       *
       */
      
      BrookFloatStreamInternal* getLinearMomentumStream( void ) const;
      
   private:
   
      // streams indices

      enum BrookVelocityCenterOfMassRemovalStreams { 
              MassStream,
              WorkStream,
              LinearMomentumStream,
              LastStreamIndex
           };

      BrookOpenMMFloat _totalInverseMass;

      // Particle stream dimensions

      int _particleStreamWidth;
      int _particleStreamHeight;
      int _particleStreamSize;

      // internal streams

      BrookFloatStreamInternal* _streams[LastStreamIndex];

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

      /** 
       * Set inverse masses 
       * 
       * @param masses             particle masses
       *
       */
      
      int _setMasses( const std::vector<double>& masses );
      
      
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_VELOCITY_CENTER_OF_MASS_REMOVAL_H_ */
