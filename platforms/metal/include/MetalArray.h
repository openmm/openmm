/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaArray.h                                 *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2009-2022 Stanford University and the Authors.      *
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

#ifndef OPENMM_METALARRAY_H_
#define OPENMM_METALARRAY_H_

#include "openmm/common/ArrayInterface.h"
#include <memory>

namespace OpenMM {

class MetalContext;

/**
 * @brief Owns a private Metal buffer through the Common ArrayInterface.
 *
 * Transfers and copies use the context's current queue at the time of the call,
 * following CudaArray/HipArray. The context must outlive this array. Synchronize
 * explicitly when accessing the same data from different queues.
 *
 * @note Nonblocking host transfers only accept a range within this context's
 *       pinned buffer. Do not read or reuse that range until the transfer's
 *       queue, or an event recorded on that queue after the transfer, has completed.
 */
class MetalArray : public ArrayInterface {
public:
    /**
     * @brief Creates an uninitialized array with no device storage.
     * @note Call initialize() before using storage or transferring data.
     */
    MetalArray();
    /**
     * @brief Creates an array and allocates its private device storage.
     * @param context The Metal context, which must outlive this array.
     * @param size The number of elements, including zero for an empty array.
     * @param elementSize The size of one element in bytes.
     * @param name The array name used for diagnostics and the buffer label.
     * @throws OpenMMException If the size is invalid or allocation fails.
     * @see initialize(ComputeContext&, size_t, int, const std::string&)
     */
    MetalArray(MetalContext& context, size_t size, int elementSize, const std::string& name);
    /** @brief Releases this array's ownership of its Metal buffer without waiting. */
    ~MetalArray();
    /**
     * @brief Initializes this array once and allocates its device storage.
     * @param context A MetalContext that must outlive this array.
     * @param size The number of elements; it must fit in an int.
     * @param elementSize A positive element size in bytes.
     * @param name The array name used for diagnostics and the buffer label.
     * @throws OpenMMException If already initialized, the byte count overflows,
     *         either size is invalid, or allocation fails.
     * @throws std::bad_cast If context is not a MetalContext.
     * @note A zero-element array is initialized with a minimal backing buffer.
     *       Newly allocated contents are not initialized by this class.
     */
    void initialize(ComputeContext& context, size_t size, int elementSize, const std::string& name) override;
    /**
     * @brief Initializes storage with sizeof(T) bytes per element.
     * @tparam T The host type whose size defines the element stride.
     * @param context A MetalContext that must outlive this array.
     * @param size The number of elements.
     * @param name The array name used for diagnostics and the buffer label.
     * @note The caller must ensure T matches the shader's element layout.
     * @see initialize(ComputeContext&, size_t, int, const std::string&)
     */
    template <class T>
    void initialize(ComputeContext& context, size_t size, const std::string& name) {
        initialize(context, size, sizeof(T), name);
    }
    /**
     * @brief Replaces the buffer, retaining the context, element size, and name.
     * @param size The new number of elements.
     * @throws OpenMMException If uninitialized, the size is invalid, or allocation fails.
     * @warning Contents are discarded, even if the new size is unchanged.
     *          Previously borrowed buffer handles refer to the old allocation.
     */
    void resize(size_t size) override;
    /** @return Whether this array owns an allocated Metal buffer. */
    bool isInitialized() const override;
    /** @return The logical element count, or zero before initialization. */
    size_t getSize() const override { return size; }
    /** @return The element size in bytes, or zero before initialization. */
    int getElementSize() const override { return elementSize; }
    /** @return The diagnostic name, or an empty string before initialization. */
    const std::string& getName() const override { return name; }
    /**
     * @return A borrowed reference to the context to which this array belongs.
     * @throws OpenMMException If the array has not been initialized.
     */
    ComputeContext& getContext() override;
    /**
     * @brief Returns a borrowed native MTLBuffer handle, not a host data pointer.
     * @return The buffer object bridged to void*, with no ownership transfer.
     * @throws OpenMMException If the array has not been initialized.
     * @note The handle is borrowed until resize() or destruction. Do not release it.
     */
    void* getBuffer() const;
    using ArrayInterface::upload;
    using ArrayInterface::download;
    /**
     * @brief Copies all elements from host memory on the current queue.
     * @param data The source of getSize()*getElementSize() bytes.
     * @param blocking Whether to wait for transfer completion before returning.
     * @note With blocking=false, data must be within the context's pinned buffer.
     * @see uploadSubArray()
     */
    void upload(const void* data, bool blocking=true) override {
        uploadSubArray(data, 0, getSize(), blocking);
    }
    /**
     * @brief Copies a host range into this array on the current queue.
     * @param data The source of elements*getElementSize() bytes; may be null for zero elements.
     * @param offset The first destination element index, not a byte offset.
     * @param elements The number of elements to copy.
     * @param blocking Whether to stage the host data and wait for completion.
     * @throws OpenMMException If uninitialized, the destination range is invalid,
     *         a nonempty source is null, a nonblocking source is outside the
     *         context's pinned buffer, or a Metal operation fails.
     * @note With blocking=false, data may include an offset into the pinned buffer,
     *       but the entire source range must fit. Keep it unchanged until completion.
     */
    void uploadSubArray(const void* data, int offset, int elements, bool blocking=true) override;
    /**
     * @brief Copies all elements to host memory on the current queue.
     * @param data The destination for getSize()*getElementSize() bytes; may be null for an empty array.
     * @param blocking Whether to wait and copy staged results into data before returning.
     * @throws OpenMMException If uninitialized, a nonempty destination is null,
     *         a nonblocking destination is outside the context's pinned buffer,
     *         or a Metal operation fails.
     * @note With blocking=false, the entire destination range must fit within the
     *       context's pinned buffer. Do not read or reuse it until completion.
     */
    void download(void* data, bool blocking=true) const override;
    /**
     * @brief Enqueues a device-to-device copy on the current queue without waiting.
     * @param dest An initialized MetalArray, or ComputeArray wrapping one, from
     *             the same context with matching element count and size.
     * @throws OpenMMException If either array is uninitialized, the destination
     *         is incompatible, or a Metal operation fails.
     * @note Copying to this array itself, or copying an empty array, does no work.
     */
    void copyTo(ArrayInterface& dest) const override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
    MetalContext* context;
    size_t size;
    int elementSize;
    std::string name;
};

} // namespace OpenMM
#endif
