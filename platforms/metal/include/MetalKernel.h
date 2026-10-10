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
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

namespace OpenMM {

class MetalArray;
class MetalContext;

/**
 * @brief Executes a Metal compute pipeline through the Common kernel interface.
 *
 * Native MSL uses consecutive direct buffer bindings. Adapted Common kernels
 * use a reflected argument buffer, preserving Common argument order without
 * the 31-direct-buffer limit. Each dispatch snapshots its primitive arguments.
 *
 * @note The context must outlive this kernel. Array arguments are non-owning
 *       references and must remain alive while bound; their current buffers are
 *       resolved at each launch so that resizing and rebinding are observed.
 */
class MetalKernel : public ComputeKernelImpl {
public:
    /**
     * @brief Creates the current ABI pipeline and lazily caches its alternate.
     * @param context The Metal context, which must outlive this kernel.
     * @param name The shader entry-point name used for diagnostics.
     * @param commonSource Whether to use reflected Common argument bindings.
     * @param libraryLookup Returns a borrowed MTLLibrary for the requested ABI.
     *        The retained callback must keep each returned library alive.
     * @param pipelineMaximum Optional compiler threadgroup maximum, or zero to
     *        retain the default pipeline. MetalProgram supplies this only for
     *        opt-in tiled force kernels with a declared workgroup size.
     */
    MetalKernel(MetalContext& context, const std::string& name, bool commonSource,
            const std::function<void*(bool)>& libraryLookup, int pipelineMaximum=0);
    /** @brief Releases the retained pipeline and argument storage without waiting. */
    ~MetalKernel();
    /** @return The shader entry-point name. */
    std::string getName() const override { return name; }
    /**
     * @return The active ABI pipeline's maximum threads per threadgroup, capped
     *         by an enabled force-pipeline compiler hint for the three tiled
     *         Nonbonded/GBSA entry points.  Other pipelines retain their limit.
     */
    int getMaxBlockSize() const override;
    /**
     * @brief Enqueues a one-dimensional launch on the context's current queue.
     * @param threads The nonnegative logical thread count; zero enqueues no work.
     * @param blockSize Threads per group, or -1 for ComputeContext::ThreadBlockSize.
     * @throws OpenMMException If the thread count or block size is invalid,
     *         the block size differs from an active exact-size requirement,
     *         arguments do not match the binding ABI, total threadgroup storage
     *         exceeds device limits, or submission fails. Native MSL is limited
     *         to 31 direct buffer arguments; adapted Common kernels are not.
     * @note Launches complete threadgroups, rounding up and capping the group
     *       count at context.getNumThreadBlocks(). The shader must handle bounds
     *       and use grid-stride iteration when the capped grid is smaller than
     *       its workload. This method does not wait for GPU completion.
     */
    void execute(int threads, int blockSize=-1) override;
    /**
     * @brief Append dynamic threadgroup storage for a translated OpenCL kernel.
     * @param bytes Nonzero byte count, validated with all other local storage at launch.
     */
    void addLocalArg(size_t bytes);
    /**
     * @brief Bind dynamic threadgroup storage at an original Common argument index.
     * @param index An existing placeholder for a LOCAL_ARG or __local parameter.
     * @param bytes Required byte count; it must be nonzero and fit device limits.
     */
    void setLocalArg(int index, size_t bytes);
protected:
    /**
     * @brief Appends a non-owning array argument at the next buffer binding.
     * @param value An initialized MetalArray, or ComputeArray wrapping one,
     *              from the same Metal device.
     * @throws OpenMMException If value is uninitialized or from another backend or device.
     * @note Linked contexts share a queue; callers must otherwise establish any
     *       required ordering with the owning context before launching work.
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
     * @note Bind it with setArg() before a nonempty launch. Native MSL's 31-slot
     *       limit is checked at execution, not when the placeholder is appended.
     */
    void addEmptyArg() override;
    /**
     * @brief Replaces an existing slot with a non-owning array binding.
     * @param index The zero-based index of an already appended argument.
     * @param value An initialized MetalArray, or ComputeArray wrapping one,
     *              from the same Metal device.
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
    /** @return Cached pipeline and binding metadata for the context's current mode. */
    Impl& getActiveImpl() const;
    mutable std::unique_ptr<Impl> variants[2];
    mutable std::mutex variantMutex;
    MetalContext& context;
    std::string name;
    bool commonSource;
    int pipelineMaximum;
    std::function<void*(bool)> libraryLookup;
    std::vector<mm_double4> primitiveArgs;
    std::vector<int> primitiveArgSizes;
    std::vector<MetalArray*> arrayArgs;
    std::vector<size_t> localArgSizes;
};

} // namespace OpenMM
#endif
