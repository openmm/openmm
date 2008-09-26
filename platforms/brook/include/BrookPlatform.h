#ifndef OPENMM_BROOKPLATFORM_H_
#define OPENMM_BROOKPLATFORM_H_

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

// default float size for Brook

#define BrookOpenMMFloat float

#include "Platform.h"
#include "BrookStreamFactory.h"

namespace OpenMM {

/**
 * This Platform subclass uses the Brook implementations of all the OpenMM kernels.
 */

class BrookPlatform : public Platform {

   public:

      // return values

      static const int DefaultReturnValue = 0;
      static const int DefaultErrorValue  = -1;
      /** 
       * BrookPlatform constructor
       *
       */
      
      BrookPlatform();

      /** 
       * BrookPlatform constructor
       *
       * @param defaultAtomStreamWidth  stream width
       * @param runtime                 Brook runtime (cal/cpu)
       * @param log                     log file reference
       *
       */
      
      BrookPlatform( int atomStreamWdith, const std::string& runtime, FILE* log = NULL );

      /** 
       * BrookPlatform destructor
       *
       */
      
	  ~BrookPlatform();

      /** 
       * Return platform name
       *
       * @return "Brook"
       */
      
      std::string getName() const;

      /** 
       * Return platform speed
       *
       * @return speed
       */

      double getSpeed( void ) const;

      /** 
       * Return true if BrookPlatform supports double precison
       *
       * @return true if BrookPlatform supports double precison
       */
      
      bool supportsDoublePrecision( void ) const;

      const StreamFactory& getDefaultStreamFactory( void ) const;
  
      /**
       * Get the runtime
       * 
       * @return runtime
       */

      std::string getRuntime( void );
  
      /** 
       * Get stream height and size given minimum number of elements in stream
       * 
       * @param size   input size of stream
       * @param width  width of stream
       * @param height output height of stream (may be NULL)
       *
       * @return stream size (=streamWidth*height)
       */

      static int getStreamSize( int size, int streamWidth, int* outputHeight );

      /** 
       * Get default stream width
       * 
       * @return default stream width
       */

      int getAtomStreamWidth( void ) const;

      /** 
       * Set log file reference
       * 
       * @param  log file reference
       *
       * @return DefaultReturnValue
       *
       */
      
      int setLog( FILE* log );

      /* 
       * Get contents of object
       *
       * @param level of dump
       *
       * @return string containing contents
       *
       * */
      
      std::string getContents( int level ) const;

      /** 
       * Get log file reference
       * 
       * @return  log file reference
       *
       */
      
      FILE* getLog( void ) const;
      
   private:

      // log file reference

      FILE* _log;

      BrookStreamFactory _defaultStreamFactory;

      // default stream width

      static const int DefaultAtomStreamWidth   = 32;

      // atom streamwidth

      int _atomStreamWidth;

      // Brook runtime

      std::string _runtime;

      /** 
       * Initialize kernel factory
       *
       */
      
      void _initializeKernelFactory( void );

      /** 
       * Set & validate runtime
       *
       * @param runtime    Brook runtime (cal/cpu)
       *
       * @throws exception if runtime is invalid
       */
      
      void _setBrookRuntime( const std::string& runtime );
      
};

} // namespace OpenMM

#endif /*OPENMM_BROOKPLATFORM_H_*/
