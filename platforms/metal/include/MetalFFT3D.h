/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform's VkFFT integration.                *
 * Source: platforms/opencl/src/OpenCLFFT3D.cpp                               *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALFFT3D_H_
#define OPENMM_METALFFT3D_H_

#include "openmm/common/FFT3D.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeKernel.h"
#include <memory>

namespace OpenMM {

class MetalContext;
class MetalVkFFT;

/**
 * @brief GPU-only, single-precision FFT using the bundled VkFFT Metal backend.
 *
 * This preserves OpenCL's unnormalized FFT3D layout and direction conventions.
 * The Objective-C++ runtime and public interface remain C++11; Metal-cpp is
 * isolated in the private VkFFT implementation. Calls on a plan must be
 * serialized, including GPU execution when switching between queues.
 */
class MetalFFT3D : public FFT3DImpl {
public:
    /**
     * @brief Prepare forward and inverse transforms of the specified dimensions.
     * @param context Borrowed Metal context, which must outlive this plan.
     * @param xsize First dimension; z is the fastest-varying index.
     * @param ysize Second dimension.
     * @param zsize Third dimension.
     * @param realToComplex Whether to use real input and a half-complex spectrum.
     * @throws OpenMMException If dimensions or GPU plan creation are invalid.
     */
    MetalFFT3D(MetalContext& context, int xsize, int ysize, int zsize, bool realToComplex=false);
    /** @brief Release the plan and its FFT resources. */
    ~MetalFFT3D();
    /**
     * @brief Append a transform to the context's currently selected queue.
     * @param in Source array, which may be overwritten as FFT workspace.
     * @param out Destination array, distinct from the source.
     * @param forward True for forward, false for inverse.
     * @note Both arrays must hold xsize*ysize*zsize complex float values, even
     *       for R2C. Its spectrum contains xsize*ysize*(zsize/2+1) complex values.
     *       Forward followed by inverse multiplies input by xsize*ysize*zsize.
     * @throws OpenMMException If buffers are incompatible or encoding fails.
     */
    void execFFT(ArrayInterface& in, ArrayInterface& out, bool forward=true) override;
    /** @return The smallest dimension >= minimum with no prime factor above 13. */
    static int findLegalDimension(int minimum);
private:
    MetalContext& context;
    size_t requiredBytes;
    bool realToComplex;
    bool packRealAsComplex;
    ComputeArray complexGrid;
    ComputeKernel packKernel, unpackKernel;
    std::unique_ptr<MetalVkFFT> plan;
};

} // namespace OpenMM

#endif // OPENMM_METALFFT3D_H_
