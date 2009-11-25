#ifndef INTEGRATE_NML_STEP_KERNEL_H
#define INTEGRATE_NML_STEP_KERNEL_H

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Chris Sweet                                                       *
 * Contributors: Christopher Bruns                                            *
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


#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include "NMLIntegrator.h"

namespace OpenMM {

/**
 * This kernel is invoked by NMLIntegrator to take one time step.
 */
class IntegrateNMLStepKernel : public KernelImpl {
    public:
        static std::string Name() {
            return "IntegrateNMLStep";
        }
        IntegrateNMLStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
        }
        /**
         * Initialize the kernel.
         *
         * @param system     the System this kernel will be applied to
         * @param integrator the NMLIntegrator this kernel will be used for
         */
        virtual void initialize(const System& system, const NMLIntegrator& integrator) = 0;
        /**
         * Execute the kernel.
         *
         * @param context    the context in which to execute this kernel
         * @param integrator the NMLIntegrator this kernel is being used for
         */
        virtual void execute(ContextImpl& context, const NMLIntegrator& integrator, const double currentPE, const int stepType) = 0;
};

} /* namespace OpenMM */

#endif /* INTEGRATE_NML_STEP_KERNEL_H */

