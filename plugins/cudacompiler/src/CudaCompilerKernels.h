#ifndef OPENMM_CUDACOMPILER_KERNELS_H_
#define OPENMM_CUDACOMPILER_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "internal/windowsExportCudaCompiler.h"
#include "CudaKernels.h"
#include <string>

namespace OpenMM {

/**
 * This kernel is a compiler for CUDA kernels based on the runtime compilation feature
 * introduced in CUDA 7.
 */
class OPENMM_EXPORT_CUDACOMPILER CudaRuntimeCompilerKernel : public CudaCompilerKernel {
public:
    CudaRuntimeCompilerKernel(const std::string& name, const Platform& platform);
    /**
     * Compile a kernel to PTX.
     *
     * @param source     the source code for the kernel
     * @param options    the flags to be passed to the compiler
     * @param cu         the CudaContext for which the kernel is being compiled
     */
    std::string createModule(const std::string& source, const std::string& flags, CudaContext& cu);
    /**
     * Get the maximum architecture version the compiler supports.
     */
    int getMaxSupportedArchitecture() const {
        return maxSupportedArchitecture;
    }
private:
    int maxSupportedArchitecture;
};

} // namespace OpenMM

#endif /*OPENMM_CUDACOMPILER_KERNELS_H_*/
