/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaProgram.h                               *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#ifndef OPENMM_METALPROGRAM_H_
#define OPENMM_METALPROGRAM_H_

#include "openmm/common/ComputeProgram.h"
#include <memory>

namespace OpenMM {

class MetalContext;

/**
 * @brief Wraps a compiled Metal library through the Common program interface.
 *
 * MetalContext::compileProgram() creates this object from native Metal 3.0
 * source. It does not translate existing Common kernel source into MSL.
 * The library is retained by this object; the context must outlive it and any
 * kernels created from it.
 */
class MetalProgram : public ComputeProgramImpl {
public:
    /**
     * @brief Retains a compiled library for kernel lookup and pipeline creation.
     * @param context The Metal context, which must outlive this program.
     * @param library A valid MTLLibrary from the context's device, bridged to void*.
     *                The caller's ownership is unchanged.
     */
    MetalProgram(MetalContext& context, void* library);
    /** @brief Releases this program's retained library without destroying its kernels. */
    ~MetalProgram();
    /**
     * @brief Creates an independently owned compute pipeline for a library entry point.
     * @param name The host-visible name of the Metal kernel function.
     * @return A shared-ownership ComputeKernel handle with no arguments bound.
     *         It may outlive this program, but not the context.
     * @throws OpenMMException If the function is missing or pipeline creation fails.
     */
    ComputeKernel createKernel(const std::string& name) override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
    MetalContext& context;
};

} // namespace OpenMM
#endif
