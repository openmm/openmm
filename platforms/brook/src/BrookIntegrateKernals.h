#ifndef OPENMM_BROOK_INTEGRATE_KERNELS_H_
#define OPENMM_BROOK_INTEGRATE_KERNELS_H_

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

namespace OpenMM {

/**
 * This is the base class of Float and Double streams in the Brook Platform.
 */

class BrookIntegrateVerletStepKernel : public IntegrateVerletStepKernel {

   public:

      /**
       * BrookIntegrateVerletStepKernel constructor
       * 
       * @param name        name of the stream to create
       * @param platform    platform
       *
       */
  
      BrookIntegrateVerletStepKernel( std::string name, const Platform& platform );

      /**
       * BrookIntegrateVerletStepKernel destructor
       * 
       */
  
      ~BrookIntegrateVerletStepKernel();

      /** 
       * Initialize the kernel, setting up all parameters related to integrator.
       * 
       * @param masses             the mass of each atom
       * @param constraintIndices  each element contains the indices of two atoms whose distance should be constrained
       * @param constraintLengths  the required distance between each pair of constrained atoms
       *
       */
      void initialize( const std::vector<double>& masses, const std::vector<std::vector<int> >& constraintIndices,
                       const std::vector<double>& constraintLengths );
      /** 
       * Execute the kernel.
       * 
       * @param positions          a Stream of type Double3 containing the position (x, y, z) of each atom
       * @param velocities         a Stream of type Double3 containing the velocity (x, y, z) of each atom
       * @param forces             a Stream of type Double3 containing the force (x, y, z) on each atom
       * @param stepSize           the integration step size
       *
       */

      void execute( Stream& positions, Stream& velocities, const Stream& forces, double stepSize );

   protected:

      int _width;
      int _streamWidth;
      int _streamHeight;
      int _streamSize;
  
};

} // namespace OpenMM

#endif /* OPENMM_BROOK_INTEGRATE_KERNELS_H_ */
