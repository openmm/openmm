/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaKernel.h                                *
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

#ifndef OPENMM_METALKERNEL_H_
#define OPENMM_METALKERNEL_H_

#include "openmm/common/ComputeKernel.h"
#include "openmm/common/ComputeVectorTypes.h"
#include <memory>
#include <vector>

namespace OpenMM {

class MetalArray;
class MetalContext;

/**
 * @brief Executes a Metal compute pipeline through the Common kernel interface.
 *
 * Arguments use consecutive direct Metal buffer bindings starting at zero.
 * Arrays and copied primitive values share the same 31 slots (indices 0-30). The
 * caller must match the native MSL signature, byte layouts, and resource-access
 * rules; this class does not inspect or translate the shader's argument ABI.
 *
 * @note The context must outlive this kernel. Array arguments are non-owning
 *       references and must remain alive while bound; their current buffers are
 *       resolved at each launch so that resizing and rebinding are observed.
 */
class MetalKernel : public ComputeKernelImpl {
public:
    /**
     * @brief Retains a native compute pipeline for subsequent launches.
     * @param context The Metal context, which must outlive this kernel.
     * @param pipeline A valid MTLComputePipelineState from the context's device,
     *                 bridged to void*. The caller's ownership is unchanged.
     * @param name The shader entry-point name used for diagnostics.
     */
    MetalKernel(MetalContext& context, void* pipeline, const std::string& name);
    /** @brief Releases the retained pipeline and argument storage without waiting. */
    ~MetalKernel();
    /** @return The shader entry-point name. */
    std::string getName() const override { return name; }
    /** @return The pipeline's maximum threads per threadgroup. */
    int getMaxBlockSize() const override;
    /**
     * @brief Enqueues a one-dimensional launch on the context's current queue.
     * @param threads The nonnegative logical thread count; zero enqueues no work.
     * @param blockSize Threads per group, or -1 for ComputeContext::ThreadBlockSize.
     * @throws OpenMMException If the thread count or block size is invalid, a
     *         nonempty launch has more than 31 or unbound arguments, or submission fails.
     * @note Launches complete threadgroups, rounding up and capping the group
     *       count at context.getNumThreadBlocks(). The shader must handle bounds
     *       and use grid-stride iteration when the capped grid is smaller than
     *       its workload. This method does not wait for GPU completion.
     */
    void execute(int threads, int blockSize=-1) override;
protected:
    /**
     * @brief Appends a non-owning array argument at the next buffer binding.
     * @param value An initialized MetalArray, or ComputeArray wrapping one,
     *              from this kernel's context.
     * @throws OpenMMException If value is uninitialized or from another backend or context.
     */
    void addArrayArg(ArrayInterface& value) override;
    /**
     * @brief Appends a primitive argument by copying its host bytes.
     * @param value A non-null pointer to the value; it need not remain alive after this call.
     * @param size The byte count, from 1 through sizeof(mm_double4) (32 bytes).
     * @throws OpenMMException If value is null or size is unsupported.
     * @note Bytes are passed without type conversion; their layout must match MSL.
     */
    void addPrimitiveArg(const void* value, int size) override;
    /**
     * @brief Appends an unbound argument slot.
     * @note Bind it with setArg() before a nonempty launch. The 31-slot limit is
     *       checked at execution, not when the placeholder is appended.
     */
    void addEmptyArg() override;
    /**
     * @brief Replaces an existing slot with a non-owning array binding.
     * @param index The zero-based index of an already appended argument.
     * @param value An initialized MetalArray, or ComputeArray wrapping one,
     *              from this kernel's context.
     * @throws OpenMMException If the index is invalid or the array is incompatible.
     * @note The previous slot may hold a primitive value, an array, or a placeholder.
     */
    void setArrayArg(int index, ArrayInterface& value) override;
    /**
     * @brief Replaces an existing slot with a copy of a primitive value.
     * @param index The zero-based index of an already appended argument.
     * @param value A non-null pointer to the host bytes to copy.
     * @param size The byte count, from 1 through sizeof(mm_double4) (32 bytes).
     * @throws OpenMMException If the index, pointer, or size is invalid.
     * @note The slot may previously hold an array or a differently sized value.
     *       The new bytes and size must match the shader's binding contract.
     */
    void setPrimitiveArg(int index, const void* value, int size) override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
    MetalContext& context;
    std::string name;
    std::vector<mm_double4> primitiveArgs;
    std::vector<int> primitiveArgSizes;
    std::vector<MetalArray*> arrayArgs;
};

} // namespace OpenMM
#endif
