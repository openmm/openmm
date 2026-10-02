/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA and OpenCL Platforms.                          *
 * Sources: platforms/cuda/include/CudaContext.h                              *
 *          platforms/opencl/include/OpenCLContext.h                         *
 *                                                                            *
 * Original CUDA/OpenCL Platform code:                                       *
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
 * @brief Single-device, single-precision implementation of Common Compute.
 *
 * Common host algorithms and generated kernels are reused through a Metal-only
 * source/binding adapter. The System and optional simulation ContextImpl must
 * outlive this object. Host calls and access to the pinned workspace must be
 * serialized; synchronize explicitly between independent queues.
 */
class MetalContext : public ComputeContext {
public:
    /** @brief Allocate the default device and state; linked contexts share a queue.
     *  @param linked Optional parent context, used by Common's inner-context forces.
     *  @param floatingAccumulators Internal test mode, inherited by linked contexts.
     */
    explicit MetalContext(const System& system, ContextImpl* simulation=nullptr,
            MetalContext* linked=nullptr, bool floatingAccumulators=false);
    /** @return Whether this context and its linked contexts use floating accumulators. */
    bool getUseFloatingPointAccumulators() const;
    /**
     * @brief Switch the accumulator format for this context and all linked contexts.
     *
     * Finish pending work before changing interpretation, then clear accumulators
     * and invalidate forces. Allocations and initialized force parameters stay
     * unchanged. Host operations on the context group must be serialized.
     */
    void setUseFloatingPointAccumulators(bool enabled);
    /** @brief Restore the mode during exception unwinding without masking the original error. */
    void restoreFloatingPointAccumulators(bool enabled) noexcept;
    /** @brief Enable and reset the optional fixed-point minimizer's contribution-range diagnostic. */
    void setCheckFixedPointRange(bool enabled);
    /** @brief Throw if any linked context recorded an unrepresentable Q32.32 contribution. */
    void checkFixedPointRange();
    /** @brief Disable diagnostic recording on every linked context during exception unwinding. */
    void restoreCheckFixedPointRange() noexcept;
    /** @return Two diagnostic words: enabled, sticky overflow; not an accumulator buffer. */
    MetalArray& getFixedPointRangeBuffer() { return fixedPointRange; }
    /** @brief Finish queued work and release owned utilities and resources. */
    ~MetalContext();
    /** @brief Initialize force-dependent buffers after all Forces have initialized. */
    void initialize();
    /** @brief Initialize this single context once. */
    void initializeContexts() override { initialize(); }
    int getNumContexts() const override { return 1; }
    int getContextIndex() const override { return 0; }
    std::vector<ComputeContext*> getAllContexts() override { return {this}; }
    /** @return The borrowed simulation context; standalone runtime contexts throw. */
    ContextImpl* getContextImpl() override;
    double& getEnergyWorkspace() override { return energyWorkspace; }
    ComputeQueue createQueue() override;
    /** @return The selected queue, validated to belong to this device. */
    MetalQueue& getCurrentMetalQueue();
    /** @return A borrowed native id<MTLDevice>; do not release it. */
    void* getDevice() const;
    MetalArray* createArray() override;
    /** @brief Unwrap a Metal array on this device, including linked-context arrays. */
    MetalArray& unwrap(ArrayInterface& array) const;
    ComputeEvent createEvent() override;
    ComputeSort createSort(ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform=true) override;
    /**
     * @brief Compile Common/OpenCL compute sources or native MSL at target 3.0.
     * @param source Runtime-generated or static kernel source.
     * @param defines Per-program macros, overriding the context defaults.
     *        Presence of OPENMM_METAL_REQUIRE_SAFE_MATH requests precise math
     *        for both accumulator variants when non-finite semantics matter.
     * @throws OpenMMException Compilation diagnostics on failure.
     */
    ComputeProgram compileProgram(const std::string source,
            const std::map<std::string, std::string>& defines={}) override;
    int computeThreadBlockSize(double memory) const override;
    /** @brief Encode a zero-fill on the current queue without waiting. */
    void clearBuffer(ArrayInterface& array) override;
    /** @brief Register a borrowed array that must outlive future clear calls. */
    void addAutoclearBuffer(ArrayInterface& array) override;
    void clearAutoclearBuffers();
    bool getIsCPU() const override { return false; }
    int getSIMDWidth() const override { return 32; }
    /** @return False: emulated accumulation does not imply native 64-bit atomics. */
    bool getSupports64BitGlobalAtomics() const override { return false; }
    bool getSupportsDoublePrecision() const override { return false; }
    bool getUseDoublePrecision() const override { return false; }
    bool getUseMixedPrecision() const override { return false; }
    int getNumAtomBlocks() const override { return paddedNumAtoms/TileSize; }
    /** @return OpenCL's Apple launch budget, or the existing fallback if unknown. */
    int getNumThreadBlocks() const override { return numComputeUnits == 0 ? 128 : 12*numComputeUnits; }
    /** @return Best-effort GPU core count; zero means the registry did not report it. */
    int getNumComputeUnits() const { return numComputeUnits; }
    int getMaxThreadBlockSize() const override;
    /** @return The hardware threadgroup-memory capacity in bytes. */
    size_t getMaxThreadgroupMemory() const;
    ArrayInterface& getPosq() override { return posq; }
    /** @throws OpenMMException Mixed precision is unsupported. */
    ArrayInterface& getPosqCorrection() override;
    ArrayInterface& getVelm() override { return velm; }
    /** @throws OpenMMException Floating-point force mirrors are not used. */
    ArrayInterface& getForceBuffers() override;
    /** @throws OpenMMException Floating-point force mirrors are not used. */
    ArrayInterface& getFloatForceBuffer() override;
    /** @return Eight-byte-capacity force planes, interpreted as packed floats during minimization. */
    ArrayInterface& getLongForceBuffer() override { return force; }
    /** @brief Decode the currently selected fixed-point or floating accumulator format. */
    void downloadFixedPointBuffer(ArrayInterface& array, std::vector<double>& values) override;
    ArrayInterface& getEnergyBuffer() override { return energyBuffer; }
    ArrayInterface& getEnergyParamDerivBuffer() override;
    /**
     * @return Borrowed shared transfer workspace.
     * @warning Do not read or reuse a pending transfer's range before completion.
     */
    void* getPinnedBuffer() override;
    ThreadPool& getThreadPool() override;
    ArrayInterface& getAtomIndexArray() override { return atomIndexArray; }
    bool getBoxIsTriclinic() const override;
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const override;
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) override;
    mm_float4 getPeriodicBoxSize() const;
    mm_float4 getInvPeriodicBoxSize() const;
    mm_float4 getPeriodicBoxVecX() const;
    mm_float4 getPeriodicBoxVecY() const;
    mm_float4 getPeriodicBoxVecZ() const;
    IntegrationUtilities& getIntegrationUtilities() override;
    ExpressionUtilities& getExpressionUtilities() override;
    BondedUtilities& getBondedUtilities() override;
    NonbondedUtilities& getNonbondedUtilities() override;
    NonbondedUtilities* createNonbondedUtilities() override;
    FFT3D createFFT(int xsize, int ysize, int zsize, bool realToComplex=false) override;
    int findLegalFFTDimension(int minimum) override;
    void setCharges(const std::vector<double>& charges) override;
    bool requestPosqCharges() override;
    const std::vector<std::string>& getEnergyParamDerivNames() const override { return energyParamDerivNames; }
    std::map<std::string, double>& getEnergyParamDerivWorkspace() override { return energyParamDerivWorkspace; }
    void addEnergyParameterDerivative(const std::string& param) override;
    /** @brief No reduction is needed: all producers accumulate into the long buffer. */
    void reduceForces() {}
    /** @brief GPU reduction followed by the usual host sum of partial energies. */
    double reduceEnergy();
    /** @brief Submit and wait, preserving the CUDA/HIP Common completion contract. */
    void flushQueue() override;
private:
    friend class MetalArray;
    /** @return The native shared transfer buffer, borrowed by MetalArray. */
    void* getPinnedBufferHandle() const;
    /** @brief Grow transfer workspace after finishing pending uses. */
    void resizePinnedBuffer(size_t bytes);
    /** @brief Compile context-owned charge and energy kernels on first use. */
    void initializeUtilityKernels();
    struct Impl;
    struct AccumulatorState;
    std::unique_ptr<Impl> impl;
    ContextImpl* simulation;
    std::shared_ptr<AccumulatorState> accumulatorState;
    MetalArray posq, velm, force, energyBuffer, energySum, atomIndexArray, chargeBuffer, energyParamDerivBuffer, fixedPointRange;
    std::unique_ptr<IntegrationUtilities> integration;
    std::unique_ptr<ExpressionUtilities> expression;
    std::unique_ptr<BondedUtilities> bonded;
    std::unique_ptr<NonbondedUtilities> nonbonded;
    ComputeKernel reduceEnergyKernel, setChargesKernel;
    std::map<std::string, std::string> compilationDefines;
    std::vector<ArrayInterface*> autoclearBuffers;
    Vec3 periodicBoxVectors[3];
    bool initialized, hasAssignedPosqCharges, flexibleBox;
    int numComputeUnits;
    double energyWorkspace;
    std::vector<std::string> energyParamDerivNames;
    std::map<std::string, double> energyParamDerivWorkspace;
};

} // namespace OpenMM
#endif
