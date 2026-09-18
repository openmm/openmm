/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA Platform.                                      *
 * Source: platforms/cuda/include/CudaContext.h                               *
 *                                                                            *
 * Original CUDA Platform code:                                               *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
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

#ifndef OPENMM_METALCONTEXT_H_
#define OPENMM_METALCONTEXT_H_

#include "MetalArray.h"
#include "openmm/common/ComputeContext.h"
#include <memory>

namespace OpenMM {

class MetalQueue;

/**
 * @brief Single-device, single-precision ComputeContext for the Metal Platform.
 *
 * This runtime foundation follows CudaContext/HipContext. It provides arrays,
 * queues, events, and native MSL compilation, but is not yet registered as a
 * simulation Platform and does not translate Common Compute kernel sources.
 * Unsupported simulation interfaces throw OpenMMException.
 *
 * @warning The System must outlive this context. Externally owned arrays,
 * programs, kernels, and events created for this context must be destroyed before
 * it. Finish work using context-owned data before destroying the context. Host
 * access to the context and its shared transfer workspace must be serialized.
 */
class MetalContext : public ComputeContext {
public:
    /**
     * @brief Create the default Metal queue and initialize the standard buffers.
     * @param system System supplying particle masses and default box vectors;
     *               the context borrows this object.
     * @throws OpenMMException If no supported Apple silicon device is available,
     *         the particle count is too large, or resource initialization fails.
     * @note Atom storage is padded to a multiple of TileSize, with at least one
     *       tile even for an empty System. Initialization finishes before return.
     */
    explicit MetalContext(const System& system);
    /** @brief Release context-owned resources and its references to the queues. */
    ~MetalContext();
    /** @return One; this implementation uses a single device context. */
    int getNumContexts() const override { return 1; }
    /** @return Zero, the index of this single device context. */
    int getContextIndex() const override { return 0; }
    /** @return A vector containing a borrowed pointer to this context. */
    std::vector<ComputeContext*> getAllContexts() override { return {this}; }
    /**
     * @brief Access the associated simulation context (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; this runtime has no simulation ContextImpl.
     */
    ContextImpl* getContextImpl() override;
    /** @return The mutable host energy accumulator, initially zero. */
    double& getEnergyWorkspace() override { return energyWorkspace; }
    /**
     * @brief Create a queue on this context's device without selecting it.
     * @return A shared owner of the new MetalQueue.
     * @throws OpenMMException If queue creation fails.
     */
    ComputeQueue createQueue() override;
    /**
     * @return A borrowed reference to the currently selected MetalQueue.
     * @throws OpenMMException If the current queue is not a MetalQueue on this device.
     */
    MetalQueue& getCurrentMetalQueue();
    /** @return A borrowed native @c id<MTLDevice> handle; do not release it. */
    void* getDevice() const;
    /** @return A new, uninitialized MetalArray; the caller owns the returned object. */
    MetalArray* createArray() override;
    /**
     * @brief Resolve a MetalArray or a ComputeArray wrapper belonging to this context.
     * @param array Array to resolve; ownership is not transferred.
     * @return A borrowed reference to the underlying MetalArray.
     * @throws OpenMMException If the array is uninitialized, has a different backend,
     *         or belongs to another context.
     */
    MetalArray& unwrap(ArrayInterface& array) const;
    /** @return A shared owner of a new, initially unrecorded MetalEvent. */
    ComputeEvent createEvent() override;
    /**
     * @brief Create a sorting utility (not implemented).
     * @param trait Sort description; ownership is taken and it is deleted before throwing.
     * @param length Requested number of elements (unused).
     * @param uniform Requested distribution hint (unused).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal sorting is not implemented.
     */
    ComputeSort createSort(ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform=true) override;
    /**
     * @brief Compile native MSL 3.0 source with fast math enabled.
     * @param source Native Metal source, not untranslated Common Compute source.
     * @param defines Macro names and replacement text prepended as preprocessor definitions.
     * @return A shared owner of the compiled MetalProgram.
     * @throws OpenMMException If Metal compilation fails, including compiler diagnostics.
     */
    ComputeProgram compileProgram(const std::string source,
            const std::map<std::string, std::string>& defines={}) override;
    /**
     * @brief Choose a SIMD-aligned group size using the device's memory and thread limits.
     * @param memory Threadgroup-memory bytes required per thread; zero ignores memory usage.
     * @return The largest permitted multiple of getSIMDWidth().
     * @throws OpenMMException If memory is negative or nonfinite, or one SIMD group cannot fit.
     * @note The selected kernel may impose a smaller maximum group size.
     */
    int computeThreadBlockSize(double memory) const override;
    /**
     * @brief Submit a GPU blit that zeros all bytes in an array on the current queue.
     * @param array Initialized array belonging to this context.
     * @throws OpenMMException If array validation or command submission fails.
     * @note This does not wait for completion. A zero-length array requires no command.
     */
    void clearBuffer(ArrayInterface& array) override;
    /**
     * @brief Register an array for subsequent calls to clearAutoclearBuffers().
     * @param array Initialized array from this context, borrowed without deduplication.
     * @throws OpenMMException If array validation fails.
     * @warning There is no removal operation; the array or wrapper must remain valid
     *          for every later clearAutoclearBuffers() call.
     */
    void addAutoclearBuffer(ArrayInterface& array) override;
    /**
     * @brief Submit clears for all registered arrays on the current queue without waiting.
     * @note Force and energy buffers are registered by the constructor. This method
     *       is not yet connected to a simulation force-evaluation lifecycle.
     * @throws OpenMMException If any array validation or clear submission fails.
     */
    void clearAutoclearBuffers();
    /** @return False; this implementation targets the GPU. */
    bool getIsCPU() const override { return false; }
    /** @return The 32-lane SIMD width assumed by this Apple silicon implementation. */
    int getSIMDWidth() const override { return 32; }
    /** @return False; 64-bit integer storage does not provide 64-bit atomic accumulation. */
    bool getSupports64BitGlobalAtomics() const override { return false; }
    /** @return False; double-precision device arithmetic is not supported. */
    bool getSupportsDoublePrecision() const override { return false; }
    /** @return False; device data and calculations use single precision. */
    bool getUseDoublePrecision() const override { return false; }
    /** @return False; mixed-precision position storage is not implemented. */
    bool getUseMixedPrecision() const override { return false; }
    /** @return The padded atom count divided by TileSize. */
    int getNumAtomBlocks() const override { return paddedNumAtoms/TileSize; }
    /** @return The fixed launch-grid cap of 128 threadgroups, not a hardware core count. */
    int getNumThreadBlocks() const override { return 128; }
    /** @return The device limit for a one-dimensional group; individual kernels may allow fewer threads. */
    int getMaxThreadBlockSize() const override;
    /** @return The owned, padded @c mm_float4 position/charge array, initially zero. */
    ArrayInterface& getPosq() override { return posq; }
    /**
     * @brief Access mixed-precision position corrections (not supported).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; no correction buffer is allocated.
     */
    ArrayInterface& getPosqCorrection() override;
    /** @return The owned, padded @c mm_float4 array of velocity xyz and inverse mass w. */
    ArrayInterface& getVelm() override { return velm; }
    /**
     * @brief Access floating-point force buffers (not supported).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; only the long force buffer is allocated.
     */
    ArrayInterface& getForceBuffers() override;
    /**
     * @brief Access a floating-point force buffer (not supported).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; only the long force buffer is allocated.
     */
    ArrayInterface& getFloatForceBuffer() override;
    /**
     * @return The owned @c int64_t force storage in consecutive x, y, z component planes,
     *         each containing the padded atom count.
     * @note This allocation alone does not implement fixed-point accumulation or force evaluation.
     */
    ArrayInterface& getLongForceBuffer() override { return force; }
    /** @return The owned, initially zero @c float energy buffer, also registered for autoclear. */
    ArrayInterface& getEnergyBuffer() override { return energyBuffer; }
    /**
     * @brief Access the energy-parameter derivative buffer (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; no derivative buffer is allocated.
     */
    ArrayInterface& getEnergyParamDerivBuffer() override;
    /**
     * @return A borrowed host pointer to the context-owned shared transfer buffer.
     * @note Its capacity covers the largest standard array at context construction;
     *       resizing arrays does not grow this buffer.
     * @warning Nonblocking array transfers require a range within this buffer. Do not
     *          read a pending download or modify/reuse a pending transfer's range until
     *          its queue or a suitably ordered event has completed.
     */
    void* getPinnedBuffer() override;
    /** @return A borrowed reference to the context-owned host thread pool. */
    ThreadPool& getThreadPool() override;
    /** @return The owned, padded @c int atom-index array, initially the identity mapping. */
    ArrayInterface& getAtomIndexArray() override { return atomIndexArray; }
    /** @return True if any stored box vector has a nonzero off-diagonal component. */
    bool getBoxIsTriclinic() const override;
    /**
     * @brief Copy the stored periodic box vectors, in nanometers.
     * @param[out] a First box vector.
     * @param[out] b Second box vector.
     * @param[out] c Third box vector.
     */
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const override;
    /**
     * @brief Store periodic box vectors without validating them or changing coordinates.
     * @param a First box vector, in nanometers.
     * @param b Second box vector, in nanometers.
     * @param c Third box vector, in nanometers.
     */
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) override;
    /**
     * @brief Access integration utilities (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal integration utilities are not implemented.
     */
    IntegrationUtilities& getIntegrationUtilities() override;
    /**
     * @brief Access expression utilities (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal expression utilities are not implemented.
     */
    ExpressionUtilities& getExpressionUtilities() override;
    /**
     * @brief Access bonded utilities (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal bonded utilities are not implemented.
     */
    BondedUtilities& getBondedUtilities() override;
    /**
     * @brief Access nonbonded utilities (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal nonbonded utilities are not implemented.
     */
    NonbondedUtilities& getNonbondedUtilities() override;
    /**
     * @brief Create nonbonded utilities (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal nonbonded utilities are not implemented.
     */
    NonbondedUtilities* createNonbondedUtilities() override;
    /**
     * @brief Create an FFT utility (not implemented).
     * @param xsize Requested x dimension (unused).
     * @param ysize Requested y dimension (unused).
     * @param zsize Requested z dimension (unused).
     * @param realToComplex Requested transform type (unused).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; Metal FFT is not implemented.
     */
    FFT3D createFFT(int xsize, int ysize, int zsize, bool realToComplex=false) override;
    /** @brief No-op; runtime buffers are initialized by the constructor, not simulation setup. */
    void initializeContexts() override {}
    /**
     * @brief Set nonbonded charges (not implemented).
     * @param charges Requested charges in elementary-charge units (unused).
     * @throws OpenMMException Always; nonbonded charge handling is not implemented.
     */
    void setCharges(const std::vector<double>& charges) override;
    /**
     * @brief Request use of the position array's charge component (not implemented).
     * @return Never returns; this method always throws.
     * @throws OpenMMException Always; nonbonded charge handling is not implemented.
     */
    bool requestPosqCharges() override;
    /** @return The empty derivative-name list; derivative registration is not implemented. */
    const std::vector<std::string>& getEnergyParamDerivNames() const override { return energyParamDerivNames; }
    /** @return Mutable host derivative bookkeeping; no device derivative evaluation is implemented. */
    std::map<std::string, double>& getEnergyParamDerivWorkspace() override { return energyParamDerivWorkspace; }
    /**
     * @brief Register an energy-parameter derivative (not implemented).
     * @param param Requested parameter name (unused).
     * @throws OpenMMException Always; energy-parameter derivatives are not implemented.
     */
    void addEnergyParameterDerivative(const std::string& param) override;
    /**
     * @brief Wait for tracked work on the current queue, matching CUDA/HIP flush behavior.
     * @throws OpenMMException If the current queue is invalid or GPU execution fails.
     * @note This is a host completion wait, not just a request to submit pending commands.
     */
    void flushQueue() override;
private:
    friend class MetalArray;
    /** @return A borrowed native @c id<MTLBuffer> handle for the shared transfer workspace. */
    void* getPinnedBufferHandle() const;
    struct Impl;
    std::unique_ptr<Impl> impl;
    MetalArray posq, velm, force, energyBuffer, atomIndexArray;
    std::vector<ArrayInterface*> autoclearBuffers;
    Vec3 periodicBoxVectors[3];
    double energyWorkspace;
    std::vector<std::string> energyParamDerivNames;
    std::map<std::string, double> energyParamDerivWorkspace;
};

} // namespace OpenMM
#endif
