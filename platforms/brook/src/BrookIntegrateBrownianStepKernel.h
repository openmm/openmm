#ifndef OPENMM_BROOK_INTEGRATE_BROWNIAN_STEP_KERNEL_H_
#define OPENMM_BROOK_INTEGRATE_BROWNIAN_STEP_KERNEL_H_

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
#include "BrookBrownianDynamics.h"
#include "BrookShakeAlgorithm.h"
#include "BrookRandomNumberGenerator.h"

namespace OpenMM {

/**
 * This is the base class of Float and Double streams in the Brook Platform.
 */

class BrookIntegrateBrownianStepKernel : public IntegrateBrownianStepKernel {

   public:

      // return values

      static const int DefaultReturnValue = 0;
      static const int ErrorReturnValue   = -1; 

      /**
       * BrookIntegrateBrownianStepKernel constructor
       * 
       * @param name        name of the stream to create
       * @param platform    platform
       *
       */
  
      BrookIntegrateBrownianStepKernel( std::string name, const Platform& platform );

      /**
       * BrookIntegrateBrownianStepKernel destructor
       * 
       */
  
      ~BrookIntegrateBrownianStepKernel();

      /** 
       * Initialize the kernel, setting up all parameters related to integrator.
       * 
       * @param masses             particle masses
       * @param constraintIndices  each element contains the indices of two particles whose distance should be constrained
       * @param constraintLengths  required distance between each pair of constrained particles
       */

      void initialize( const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
                       const std::vector<double>& constraintLengths );

      /** 
       * Execute kernel
       * 
       * @param positions          coordinates
       * @param velocities         velocities
       * @param forces             forces
       * @param temperature        heat bath temperature
       * @param friction           friction coefficient coupling the system to the heat bath
       * @param stepSize           step size
       *
       */

      void execute( Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize );

   protected:

      BrookBrownianDynamics*        _brookBrownianDynamics;
      BrookShakeAlgorithm*          _brookShakeAlgorithm;
      BrookRandomNumberGenerator*   _brookRandomNumberGenerator;

};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INTEGRATE_BROWNIAN_STEP_KERNEL_H_ */
