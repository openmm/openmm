#ifndef OPENMM_BROOKSTREAMFACTORY_H_
#define OPENMM_BROOKSTREAMFACTORY_H_

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

#include "openmm/StreamFactory.h"
#include "openmm/Platform.h"

namespace OpenMM {

/**
 * This StreamFactory creates all streams for BrookPlatform.
 */

class OPENMM_EXPORT BrookStreamFactory : public StreamFactory {

   public:

      /** 
       * BrookStreamFactory constructor
       * 
       */
      
      BrookStreamFactory( );

      /** 
       * BrookStreamFactory destructor
       * 
       */
      
	  ~BrookStreamFactory( );

      // 'external' streams 

      static const std::string ParticlePositions;
      static const std::string ParticleVelocities;
      static const std::string ParticleForces;

      /** 
       * Create StreamImpl
       *
       * @param name     stream name
       * @param size     stream size
       * @param type     data type (float, float2, ...)
       * @param platform platform reference
       * @param context  context (currently ignored)
       * 
       * @return StreamImpl
       */
      
      StreamImpl* createStreamImpl( std::string name, int size, Stream::DataType type, const Platform& platform, ContextImpl& context ) const;

      /** 
       * Create StreamImpl
       *
       * @param name     stream name
       * @param size     stream size
       * @param type     data type (float, float2, ...)
       * @param platform platform reference
       * 
       * @return StreamImpl
       */
      
      StreamImpl* createStreamImpl( std::string name, int size, Stream::DataType type, const Platform& platform ) const;

      /** 
       * Get particle stream width
       * 
       * @return particle stream width
       *
       *
       */

      int getDefaultParticleStreamWidth( void ) const;

      /** 
       * Set particle stream width
       * 
       * @param particleStreamWidth  particle stream width
       *
       * @return DefaultReturnValue
       *
       * @throw OpenMMException if particleStreamWidth < 1
       *
       */

      int setDefaultParticleStreamWidth( int particleStreamWidth );

      /** 
       * Get randomNumber stream width
       * 
       * @return randomNumber stream width
       *
       *
       */

      int getDefaultRandomNumberStreamWidth( void ) const;

      /** 
       * Set randomNumber stream width
       * 
       * @param randomNumberStreamWidth  randomNumber stream width
       *
       * @return DefaultReturnValue
       *
       * @throw OpenMMException if randomNumberStreamWidth < 1
       *
       */

      int setDefaultRandomNumberStreamWidth( int randomNumberStreamWidth );

      /** 
       * Get randomNumber stream size
       * 
       * @return randomNumber stream size
       *
       *
       */

      int getDefaultRandomNumberStreamSize( void ) const;

      /** 
       * Set randomNumber stream size
       * 
       * @param randomNumberStreamSize  randomNumber stream size
       *
       * @return DefaultReturnValue
       *
       * @throw OpenMMException if randomNumberStreamSize < 1
       *
       */

      int setDefaultRandomNumberStreamSize( int randomNumberStreamSize );

      /** 
       * Get default dangle value
       * 
       * @return default dangle value
       *
       */

      double getDefaultDangleValue( void ) const;

      /** 
       * Set default dangle value
       * 
       * @param DefaultDangleValue default dangle value
       *
       * @return DefaultReturnValue
       *
       */

      int setDefaultDangleValue( double defaultDangleValue );

   private:

      static const int DefaultStreamParticleWidth         = 32;

      static const int DefaultStreamRandomNumberWidth     = 1024;
      static const int DefaultStreamRandomNumberSize      = 1024*1024;

      static const double DefaultDangleValue;

      static const int DefaultReturnValue                 = 0;
      static const int ErrorReturnValue                   = -1; 

      int _defaultParticleStreamWidth;

      int _defaultStreamRandomNumberWidth;
      int _defaultStreamRandomNumberSize;

      double _defaultDangleValue;

};

} // namespace OpenMM

#endif /*OPENMM_BROOKSTREAMFACTORY_H_*/
