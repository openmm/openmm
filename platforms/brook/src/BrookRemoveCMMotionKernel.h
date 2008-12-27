#ifndef OPENMM_BROOK_REMOVE_CM_MOTION_KERNEL_H_
#define OPENMM_BROOK_REMOVE_CM_MOTION_KERNEL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright ( c ) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files ( the "Software" ), *
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

#include "kernels.h"
#include "OpenMMBrookInterface.h"
#include "BrookVelocityCenterOfMassRemoval.h"

namespace OpenMM {

/**
 * Brook class for removing center-of-mass motion
 */

class BrookRemoveCMMotionKernel : public RemoveCMMotionKernel {

   public:

      /**
       * BrookRemoveCMMotionKernel constructor
       * 
       * @param name                      name of the stream to create
       * @param platform                  platform
       * @param openMMBrookInterface      OpenMM-Brook interface
       * @param System                    System reference
       *
       */
  
      BrookRemoveCMMotionKernel( std::string name, const Platform& platform, OpenMMBrookInterface& openMMBrookInterface, System& system );

      /**
       * BrookRemoveCMMotionKernel destructor
       * 
       */
  
      ~BrookRemoveCMMotionKernel();

      /** 
       * Get COM removal frequency
       * 
       * @return frequency 
       *
       */
      
      int getFrequency( void ) const;
       
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

      /** 
       * Initialize the kernel
       * 
       * @param system   System reference
       * @param force    CMMotionRemover reference
       *
       */
      void initialize( const System& system, const CMMotionRemover& force );

      /** 
       * Execute the kernel.
       * 
       * @param context  OpenMMContextImpl reference
       *
       */

      void execute( OpenMMContextImpl& context );

   private:

      int _frequency;

      // log file reference

      FILE* _log;

      BrookVelocityCenterOfMassRemoval* _brookVelocityCenterOfMassRemoval;

      // interface

      OpenMMBrookInterface& _openMMBrookInterface;

      // System reference

      System& _system;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_REMOVE_CM_MOTION_KERNEL_H_ */
