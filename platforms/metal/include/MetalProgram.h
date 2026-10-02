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
#include <map>
#include <memory>

namespace OpenMM {

class MetalContext;

/**
 * @brief Wraps a compiled Metal library through the Common program interface.
 *
 * MetalContext::compileProgram() creates this object from native Metal
 * source or Common source adapted to Metal's argument-buffer interface.
 * Common programs lazily retain fixed-point and floating-accumulator libraries.
 * Both libraries use the context's same runtime-selected language target.
 * Kernels share their library cache and may outlive this program; the context
 * must outlive both the program and every kernel created from it.
 */
class MetalProgram : public ComputeProgramImpl {
public:
    /**
     * @brief Retains a compiled library for kernel lookup and pipeline creation.
     * @param context The Metal context, which must outlive this program.
     * @param library A valid MTLLibrary from the context's device, bridged to void*.
     *                The caller's ownership is unchanged.
     * @param commonSource Whether entry points use translated Common bindings.
     * @param source Original Common source for a lazily compiled alternate ABI.
     * @param defines Final per-program definitions, including selected fast paths.
     * @param strictMath Preserve compensated arithmetic even in fixed-point mode.
     */
    MetalProgram(MetalContext& context, void* library, bool commonSource=false,
            const std::string& source="", const std::map<std::string, std::string>& defines={}, bool strictMath=false);
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
    std::shared_ptr<Impl> impl;
    MetalContext& context;
    bool commonSource;
};

} // namespace OpenMM
#endif
