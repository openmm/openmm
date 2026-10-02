/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM CUDA and OpenCL Platforms.                          *
 * Sources: platforms/cuda/src/CudaContext.cpp                                *
 *          platforms/opencl/src/OpenCLContext.cpp                            *
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

#include "MetalContext.h"
#include "MetalEvent.h"
#include "MetalFFT3D.h"
#include "MetalIntegrationUtilities.h"
#include "MetalKernel.h"
#include "MetalKernelSources.h"
#include "MetalNonbondedUtilities.h"
#include "MetalOpenCLKernelSources.h"
#include "MetalProgram.h"
#include "MetalQueue.h"
#include "MetalSort.h"
#include "MetalSourceAdapter.h"
#include "CommonKernelSources.h"
#include "openmm/MonteCarloFlexibleBarostat.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#import <Metal/Metal.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>

using namespace OpenMM;
using namespace std;

struct MetalContext::Impl {
    id<MTLDevice> device = nil;
    id<MTLBuffer> pinnedBuffer = nil;
    ThreadPool threads;
    Impl() : threads(1) {
    }
};

/** One mode switch covers every linked Context that can share accumulator arrays. */
struct MetalContext::AccumulatorState {
    bool floating;
    vector<MetalContext*> contexts;
    explicit AccumulatorState(bool floating) : floating(floating) {
    }
};

MetalContext::MetalContext(const System& system, ContextImpl* simulation, MetalContext* linked, bool floatingAccumulators) : ComputeContext(system),
        simulation(simulation), accumulatorState(linked == nullptr ?
                make_shared<AccumulatorState>(floatingAccumulators) : linked->accumulatorState),
        initialized(false), hasAssignedPosqCharges(false), flexibleBox(false), energyWorkspace(0.0) {
    try {
        impl.reset(new Impl());
        @autoreleasepool {
            impl->device = (linked == nullptr ? MTLCreateSystemDefaultDevice() : (__bridge id<MTLDevice>) linked->getDevice());
            if (impl->device == nil)
                throw OpenMMException("No Metal device is available");
            if (![impl->device supportsFamily:MTLGPUFamilyApple7])
                throw OpenMMException("The Metal Platform requires Apple silicon");
            // Inner CustomCV/ATM contexts share the parent's stream, just as in
            // OpenCL, so state copies and force evaluation stay device-ordered.
            defaultQueue = (linked == nullptr ? createQueue() : linked->getCurrentQueue());
            restoreDefaultQueue();
            numAtoms = system.getNumParticles();
            if (numAtoms > numeric_limits<int>::max()-TileSize)
                throw OpenMMException("Too many atoms for the Metal Platform");
            paddedNumAtoms = max(TileSize, TileSize*((numAtoms+TileSize-1)/TileSize));
            posq.initialize<mm_float4>(*this, paddedNumAtoms, "posq");
            velm.initialize<mm_float4>(*this, paddedNumAtoms, "velm");
            // Keep allocations stable across scoped minimization. Float kernels
            // use the first N packed floats of each N-element fixed-point array.
            force.initialize<int64_t>(*this, 3*size_t(paddedNumAtoms), "force");
            fixedPointRange.initialize<uint32_t>(*this, 2, "fixedPointRange");
            fixedPointRange.upload(vector<uint32_t>{0, 0});
            energyBuffer.initialize<float>(*this, getNumThreadBlocks()*ThreadBlockSize, "energyBuffer");
            atomIndexArray.initialize<int>(*this, paddedNumAtoms, "atomIndex");
            atomIndex.resize(paddedNumAtoms);
            iota(atomIndex.begin(), atomIndex.end(), 0);
            atomIndexArray.upload(atomIndex);
            posCellOffsets.resize(paddedNumAtoms, mm_int4(0, 0, 0, 0));
            vector<mm_float4> velocities(paddedNumAtoms, mm_float4(0, 0, 0, 0));
            for (int i = 0; i < numAtoms; i++) {
                double mass = system.getParticleMass(i);
                velocities[i].w = (mass == 0.0 ? 0.0f : (float) (1.0/mass));
            }
            velm.upload(velocities);
            clearBuffer(posq);
            clearBuffer(force);
            clearBuffer(energyBuffer);
            addAutoclearBuffer(force);
            addAutoclearBuffer(energyBuffer);
            size_t pinnedBytes = max(force.getSize()*force.getElementSize(),
                    energyBuffer.getSize()*energyBuffer.getElementSize());
            pinnedBytes = max(pinnedBytes, posq.getSize()*posq.getElementSize());
            impl->pinnedBuffer = [impl->device newBufferWithLength:pinnedBytes options:MTLResourceStorageModeShared];
            if (impl->pinnedBuffer == nil)
                throw OpenMMException("Error creating Metal pinned transfer buffer");
            system.getDefaultPeriodicBoxVectors(periodicBoxVectors[0], periodicBoxVectors[1], periodicBoxVectors[2]);
            for (int i = 0; i < system.getNumForces(); i++)
                if (dynamic_cast<const MonteCarloFlexibleBarostat*>(&system.getForce(i)) != nullptr)
                    flexibleBox = true;
            getCurrentMetalQueue().finish();
            compilationDefines["OPENMM_METAL_NATIVE_FLOAT_ATOMICS"] = to_string(OPENMM_METAL_NATIVE_FLOAT_ATOMICS);
            compilationDefines["OPENMM_METAL_FAST_SPARSE_FORCE_AGGREGATION"] = to_string(OPENMM_METAL_FAST_SPARSE_FORCE_AGGREGATION);
            // MSL 3.0 exposes only uint64 min/max, not uint64 add/CAS/load.
            // These operations additionally require Apple8 or newer on macOS.
            compilationDefines["OPENMM_METAL_HAS_UINT64_MIN_MAX"] =
                [impl->device supportsFamily:MTLGPUFamilyApple8] ? "1" : "0";
            compilationDefines["OPENMM_METAL_CHECK_FIXED_POINT_RANGE"] = to_string(!OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS);
            // Preserve the OpenCL periodic-image convention for shared kernels.
            if (getBoxIsTriclinic()) {
                compilationDefines["APPLY_PERIODIC_TO_DELTA(delta)"] =
                    "{"
                    "real scale3 = floor(delta.z*invPeriodicBoxSize.z+0.5f); \\\n"
                    "delta.xyz -= scale3*periodicBoxVecZ.xyz; \\\n"
                    "real scale2 = floor(delta.y*invPeriodicBoxSize.y+0.5f); \\\n"
                    "delta.xy -= scale2*periodicBoxVecY.xy; \\\n"
                    "real scale1 = floor(delta.x*invPeriodicBoxSize.x+0.5f); \\\n"
                    "delta.x -= scale1*periodicBoxVecX.x;}";
                compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
                    "{"
                    "real scale3 = floor(pos.z*invPeriodicBoxSize.z); \\\n"
                    "pos.xyz -= scale3*periodicBoxVecZ.xyz; \\\n"
                    "real scale2 = floor(pos.y*invPeriodicBoxSize.y); \\\n"
                    "pos.xy -= scale2*periodicBoxVecY.xy; \\\n"
                    "real scale1 = floor(pos.x*invPeriodicBoxSize.x); \\\n"
                    "pos.x -= scale1*periodicBoxVecX.x;}";
                compilationDefines["APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)"] =
                    "{"
                    "real scale3 = floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f); \\\n"
                    "pos.x -= scale3*periodicBoxVecZ.x; \\\n"
                    "pos.y -= scale3*periodicBoxVecZ.y; \\\n"
                    "pos.z -= scale3*periodicBoxVecZ.z; \\\n"
                    "real scale2 = floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f); \\\n"
                    "pos.x -= scale2*periodicBoxVecY.x; \\\n"
                    "pos.y -= scale2*periodicBoxVecY.y; \\\n"
                    "real scale1 = floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f); \\\n"
                    "pos.x -= scale1*periodicBoxVecX.x;}";
            }
            else {
                compilationDefines["APPLY_PERIODIC_TO_DELTA(delta)"] =
                    "delta.xyz -= floor(delta.xyz*invPeriodicBoxSize.xyz+0.5f)*periodicBoxSize.xyz;";
                compilationDefines["APPLY_PERIODIC_TO_POS(pos)"] =
                    "pos.xyz -= floor(pos.xyz*invPeriodicBoxSize.xyz)*periodicBoxSize.xyz;";
                compilationDefines["APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)"] =
                    "{"
                    "pos.x -= floor((pos.x-center.x)*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x; \\\n"
                    "pos.y -= floor((pos.y-center.y)*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y; \\\n"
                    "pos.z -= floor((pos.z-center.z)*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;}";
            }
        }
        accumulatorState->contexts.push_back(this);
    }
    catch (...) {
        delete workThread;
        workThread = nullptr;
        throw;
    }
}

MetalContext::~MetalContext() {
    auto& contexts = accumulatorState->contexts;
    contexts.erase(remove(contexts.begin(), contexts.end(), this), contexts.end());
    // Drain before destroying utility-owned buffers; destruction must not throw.
    try {
        if (currentQueue)
            getCurrentMetalQueue().finish();
        if (defaultQueue && defaultQueue != currentQueue)
            static_cast<MetalQueue&>(*defaultQueue).finish();
    }
    catch (...) {
    }
    delete workThread;
    for (auto* value : forces)
        delete value;
    for (auto* value : reorderListeners)
        delete value;
    for (auto* value : preComputations)
        delete value;
    for (auto* value : postComputations)
        delete value;
    // Queue destructors finish pending work without throwing from destruction.
    currentQueue.reset();
    defaultQueue.reset();
}

bool MetalContext::getUseFloatingPointAccumulators() const {
    return accumulatorState->floating;
}

void MetalContext::setCheckFixedPointRange(bool enabled) {
    // Uploads are queue-ordered. Keep the sticky bit through every force
    // evaluation in one solve, including constraint and linked-context kernels.
    for (auto* context : accumulatorState->contexts)
        context->fixedPointRange.upload(vector<uint32_t>{enabled ? 1u : 0u, 0u});
}

void MetalContext::checkFixedPointRange() {
    bool overflow = false;
    for (auto* context : accumulatorState->contexts) {
        vector<uint32_t> words;
        context->fixedPointRange.download(words);
        overflow |= words[1] != 0;
    }
    if (overflow)
        throw OpenMMException("Minimization exceeded the GPU Q32.32 force range; CPU fallback is disabled for this platform");
}

void MetalContext::restoreCheckFixedPointRange() noexcept {
    for (auto* context : accumulatorState->contexts) {
        try {
            context->fixedPointRange.upload(vector<uint32_t>{0, 0});
        }
        catch (...) {
        }
    }
}

void MetalContext::setUseFloatingPointAccumulators(bool enabled) {
    if (enabled == accumulatorState->floating)
        return;
    // Drain every queue even when one reports an error. No queued dispatch may
    // still interpret a shared buffer in the old format when the mode changes.
    exception_ptr error;
    for (auto* context : accumulatorState->contexts) {
        for (const auto& queue : {context->currentQueue, context->defaultQueue}) {
            try {
                if (queue)
                    static_cast<MetalQueue&>(*queue).finish();
            }
            catch (...) {
                if (!error)
                    error = current_exception();
            }
        }
    }
    if (error)
        rethrow_exception(error);
    accumulatorState->floating = enabled;
    for (auto* context : accumulatorState->contexts)
        context->setForcesValid(false);
    for (auto* context : accumulatorState->contexts) {
        if (context->simulation != nullptr)
            context->simulation->systemChanged();
    }
    for (auto* context : accumulatorState->contexts)
        context->clearAutoclearBuffers();
}

void MetalContext::restoreFloatingPointAccumulators(bool enabled) noexcept {
    try {
        setUseFloatingPointAccumulators(enabled);
    }
    catch (...) {
        // Preserve the initiating exception, but never leave the host selecting
        // the temporary ABI. A later force evaluation clears its accumulators.
        accumulatorState->floating = enabled;
        for (auto* context : accumulatorState->contexts) {
            context->setForcesValid(false);
            try {
                if (context->simulation != nullptr)
                    context->simulation->systemChanged();
            }
            catch (...) {
            }
        }
    }
}

void* MetalContext::getDevice() const {
    return (__bridge void*) impl->device;
}

ContextImpl* MetalContext::getContextImpl() {
    if (simulation == nullptr)
        throw OpenMMException("The Metal Platform is not attached to a simulation Context");
    return simulation;
}

ComputeQueue MetalContext::createQueue() {
    return ComputeQueue(new MetalQueue(getDevice()));
}

MetalQueue& MetalContext::getCurrentMetalQueue() {
    MetalQueue* queue = dynamic_cast<MetalQueue*>(currentQueue.get());
    if (queue == nullptr)
        throw OpenMMException("The current queue is not a MetalQueue");
    id<MTLCommandQueue> native = (__bridge id<MTLCommandQueue>) queue->getQueue();
    if (native.device != impl->device)
        throw OpenMMException("The current Metal queue belongs to a different device");
    return *queue;
}

MetalArray* MetalContext::createArray() {
    return new MetalArray();
}

MetalArray& MetalContext::unwrap(ArrayInterface& array) const {
    ComputeArray* wrapper = dynamic_cast<ComputeArray*>(&array);
    MetalArray* metalArray = dynamic_cast<MetalArray*>(wrapper == nullptr ? &array : &wrapper->getArray());
    if (metalArray == nullptr || dynamic_cast<MetalContext&>(metalArray->getContext()).getDevice() != getDevice())
        throw OpenMMException("Array argument is not a MetalArray from this device");
    if (dynamic_cast<MetalContext&>(metalArray->getContext()).getUseFloatingPointAccumulators() != getUseFloatingPointAccumulators())
        throw OpenMMException("Cannot share accumulator arrays between ordinary and minimization Metal contexts");
    return *metalArray;
}

ComputeEvent MetalContext::createEvent() {
    return ComputeEvent(new MetalEvent(*this));
}

void MetalContext::downloadFixedPointBuffer(ArrayInterface& array, vector<double>& values) {
    unwrap(array);
    if (!getUseFloatingPointAccumulators()) {
        ComputeContext::downloadFixedPointBuffer(array, values);
        return;
    }
    if (array.getElementSize() != sizeof(float) && array.getElementSize() != sizeof(int64_t))
        throw OpenMMException("Floating accumulator decoding requires four- or eight-byte elements");
    values.resize(array.getSize());
    if (values.empty())
        return;
    // Keep fixed-point allocation capacity while decoding only the logical
    // element count from the temporary packed-float representation.
    vector<float> storage(array.getSize()*array.getElementSize()/sizeof(float));
    array.download(storage.data());
    for (size_t i = 0; i < values.size(); i++)
        values[i] = storage[i];
}

ComputeProgram MetalContext::compileProgram(const string source, const map<string, string>& defines) {
    @autoreleasepool {
        const bool floatingAccumulators = getUseFloatingPointAccumulators();
        string code;
        bool commonSource = MetalSourceAdapter::isCommonSource(source);
        map<string, string> allDefines = commonSource ? compilationDefines : map<string, string>();
#if OPENMM_METAL_FAST_MINIMIZE_SHUFFLE
        // Scope this switch to the exact minimizer program, not other solvers
        // that happen to use the same Common CUDA/HIP shuffle macro.
        if (source == CommonKernelSources::minimize)
            allDefines["WARP_SHUFFLE_DOWN(value, offset)"] = "simd_shuffle_down(value, offset)";
#endif
#if OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION
        if (source == CommonKernelSources::constantPotential)
            allDefines["WARP_SHUFFLE_DOWN(value, offset)"] = "simd_shuffle_down(value, offset)";
#endif
#if OPENMM_METAL_FAST_CONSTANT_POTENTIAL_CG_REDUCTION
        if (source == CommonKernelSources::constantPotentialCGSolver)
            allDefines["WARP_SHUFFLE_DOWN(value, offset)"] = "simd_shuffle_down(value, offset)";
#endif
#if OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_REDUCTION
        if (source == CommonKernelSources::constantPotentialMatrixSolver)
            allDefines["WARP_SHUFFLE_DOWN(value, offset)"] = "simd_shuffle_down(value, offset)";
#endif
#if OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_BROADCAST
        if (source == CommonKernelSources::constantPotentialMatrixSolver)
            allDefines["WARP_SHUFFLE(value, index)"] = "simd_shuffle(value, index)";
#endif
        for (auto& define : defines)
            allDefines[define.first] = define.second;
        if (commonSource && floatingAccumulators)
            allDefines["OPENMM_METAL_FLOAT_ACCUMULATORS"] = "1";
        for (auto& define : allDefines)
            code += "#define "+define.first+" "+define.second+"\n";
        code += commonSource ? MetalKernelSources::common+MetalSourceAdapter::translate(source, floatingAccumulators) : source;
        MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
        options.languageVersion = MTLLanguageVersion3_0;
        // CG needs compensated low terms; FP16 bounds require conservative Inf
        // overflow behavior. Floating minimization needs finite-value checks.
        const bool strictMath = source == CommonKernelSources::constantPotentialCGSolver ||
                allDefines.count("OPENMM_METAL_REQUIRE_SAFE_MATH") != 0;
        const bool fastMath = OPENMM_METAL_FAST_MATH && !floatingAccumulators && !strictMath;
        if (@available(macOS 15.0, *)) {
            options.mathMode = fastMath ? MTLMathModeFast : MTLMathModeSafe;
            options.mathFloatingPointFunctions = fastMath ? MTLMathFloatingPointFunctionsFast : MTLMathFloatingPointFunctionsPrecise;
        }
        else {
            // Compatibility only: these replacement properties require macOS 15.
            options.fastMathEnabled = fastMath;
        }
        NSError* error = nil;
        id<MTLLibrary> library = [impl->device newLibraryWithSource:[NSString stringWithUTF8String:code.c_str()]
                options:options error:&error];
        if (library == nil)
            throw OpenMMException("Error compiling Metal program: "+
                    (error == nil ? string("unknown error") : string(error.localizedDescription.UTF8String)));
        return ComputeProgram(new MetalProgram(*this, (__bridge void*) library, commonSource,
                source, allDefines, strictMath));
    }
}

int MetalContext::getMaxThreadBlockSize() const {
    return impl->device.maxThreadsPerThreadgroup.width;
}

size_t MetalContext::getMaxThreadgroupMemory() const {
    return impl->device.maxThreadgroupMemoryLength;
}

int MetalContext::computeThreadBlockSize(double memory) const {
    if (!isfinite(memory) || memory < 0)
        throw OpenMMException("Invalid shared memory requirement");
    int maximum = getMaxThreadBlockSize();
    if (memory > 0)
        maximum = min(double(maximum), floor(impl->device.maxThreadgroupMemoryLength/memory));
    if (maximum < getSIMDWidth())
        throw OpenMMException("Shared memory requirement exceeds Metal threadgroup capacity");
    return (maximum/getSIMDWidth())*getSIMDWidth();
}

void MetalContext::clearBuffer(ArrayInterface& array) {
    MetalArray& metalArray = unwrap(array);
    if (metalArray.getSize() == 0)
        return;
    @autoreleasepool {
        MetalQueue& queue = getCurrentMetalQueue();
        auto queueLock = queue.lock();
        id<MTLCommandBuffer> command = (__bridge id<MTLCommandBuffer>) queue.getCommandBuffer();
        id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
        if (encoder == nil)
            throw OpenMMException("Error creating Metal clear command");
        [encoder fillBuffer:(__bridge id<MTLBuffer>) metalArray.getBuffer()
                range:NSMakeRange(0, metalArray.getSize()*metalArray.getElementSize()) value:0];
        [encoder endEncoding];
        queue.submit((__bridge void*) command);
    }
}

void MetalContext::addAutoclearBuffer(ArrayInterface& array) {
    unwrap(array);
    autoclearBuffers.push_back(&array);
}

void MetalContext::clearAutoclearBuffers() {
    // OpenCL clears up to six buffers per dispatch. Preserve that submission
    // granularity in immediate mode too, without general command batching.
    for (size_t base = 0; base < autoclearBuffers.size(); base += 6) {
        MetalArray* arrays[6];
        int count = 0;
        for (size_t i = base; i < min(base+6, autoclearBuffers.size()); i++) {
            MetalArray& array = unwrap(*autoclearBuffers[i]);
            if (array.getSize() != 0)
                arrays[count++] = &array;
        }
        if (count == 0)
            continue;
        @autoreleasepool {
            MetalQueue& queue = getCurrentMetalQueue();
            auto queueLock = queue.lock();
            id<MTLCommandBuffer> command = (__bridge id<MTLCommandBuffer>) queue.getCommandBuffer();
            id<MTLBlitCommandEncoder> encoder = [command blitCommandEncoder];
            if (encoder == nil)
                throw OpenMMException("Error creating Metal autoclear command");
            for (int i = 0; i < count; i++)
                [encoder fillBuffer:(__bridge id<MTLBuffer>) arrays[i]->getBuffer()
                        range:NSMakeRange(0, arrays[i]->getSize()*arrays[i]->getElementSize()) value:0];
            [encoder endEncoding];
            queue.submit((__bridge void*) command);
        }
    }
}

void* MetalContext::getPinnedBuffer() {
    return impl->pinnedBuffer.contents;
}

void* MetalContext::getPinnedBufferHandle() const {
    return (__bridge void*) impl->pinnedBuffer;
}

ThreadPool& MetalContext::getThreadPool() {
    return impl->threads;
}

bool MetalContext::getBoxIsTriclinic() const {
    return flexibleBox || periodicBoxVectors[0][1] != 0.0 || periodicBoxVectors[0][2] != 0.0 ||
            periodicBoxVectors[1][0] != 0.0 || periodicBoxVectors[1][2] != 0.0 ||
            periodicBoxVectors[2][0] != 0.0 || periodicBoxVectors[2][1] != 0.0;
}

void MetalContext::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void MetalContext::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

void MetalContext::flushQueue() {
    // Match CudaContext/HipContext: flush also waits for the current stream.
    getCurrentMetalQueue().finish();
}

ArrayInterface& MetalContext::getPosqCorrection() {
    throw OpenMMException("Metal does not use mixed precision");
}

ArrayInterface& MetalContext::getForceBuffers() {
    throw OpenMMException("Metal does not use floating point force buffers");
}

ArrayInterface& MetalContext::getFloatForceBuffer() {
    throw OpenMMException("Metal does not use floating point force buffers");
}

ArrayInterface& MetalContext::getEnergyParamDerivBuffer() {
    return energyParamDerivBuffer;
}

void MetalContext::addEnergyParameterDerivative(const string& param) {
    if (find(energyParamDerivNames.begin(), energyParamDerivNames.end(), param) == energyParamDerivNames.end()) {
        if (initialized)
            throw OpenMMException("Energy parameter derivatives must be registered before initialization");
        energyParamDerivNames.push_back(param);
    }
}

ComputeSort MetalContext::createSort(ComputeSortImpl::SortTrait* trait, unsigned int length, bool uniform) {
    return ComputeSort(new MetalSort(*this, trait, length, uniform));
}

IntegrationUtilities& MetalContext::getIntegrationUtilities() {
    if (!integration)
        integration.reset(new MetalIntegrationUtilities(*this, system));
    return *integration;
}

ExpressionUtilities& MetalContext::getExpressionUtilities() {
    if (!expression)
        expression.reset(new ExpressionUtilities(*this));
    return *expression;
}

BondedUtilities& MetalContext::getBondedUtilities() {
    if (!bonded)
        bonded.reset(new BondedUtilities(*this));
    return *bonded;
}

NonbondedUtilities& MetalContext::getNonbondedUtilities() {
    if (!nonbonded)
        nonbonded.reset(createNonbondedUtilities());
    return *nonbonded;
}

NonbondedUtilities* MetalContext::createNonbondedUtilities() {
    return new MetalNonbondedUtilities(*this);
}

FFT3D MetalContext::createFFT(int xsize, int ysize, int zsize, bool realToComplex) {
    return FFT3D(new MetalFFT3D(*this, xsize, ysize, zsize, realToComplex));
}

int MetalContext::findLegalFFTDimension(int minimum) {
    return MetalFFT3D::findLegalDimension(minimum);
}

void MetalContext::setCharges(const vector<double>& charges) {
    if (charges.size() != numAtoms)
        throw OpenMMException("Wrong number of Metal particle charges");
    initializeUtilityKernels();
    if (!chargeBuffer.isInitialized())
        chargeBuffer.initialize<float>(*this, numAtoms, "chargeBuffer");
    vector<float> values(charges.begin(), charges.end());
    chargeBuffer.upload(values.data());
    setChargesKernel->setArg(0, chargeBuffer);
    setChargesKernel->setArg(1, posq);
    setChargesKernel->setArg(2, atomIndexArray);
    setChargesKernel->setArg(3, numAtoms);
    setChargesKernel->execute(numAtoms);
}

bool MetalContext::requestPosqCharges() {
    bool available = !hasAssignedPosqCharges;
    hasAssignedPosqCharges = true;
    return available;
}

void MetalContext::resizePinnedBuffer(size_t bytes) {
    if (impl->pinnedBuffer.length >= bytes)
        return;
    getCurrentMetalQueue().finish();
    id<MTLBuffer> buffer = [impl->device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    if (buffer == nil)
        throw OpenMMException("Error growing the Metal pinned transfer buffer");
    impl->pinnedBuffer = buffer;
}

void MetalContext::initialize() {
    if (initialized)
        return;
    getBondedUtilities().initialize(system);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(getNonbondedUtilities());
    int energySize = max(getNumThreadBlocks()*ThreadBlockSize, nb.getNumEnergyBuffers());
    if (energyBuffer.getSize() != energySize)
        energyBuffer.resize(energySize);
    if (!energyParamDerivNames.empty()) {
        energyParamDerivBuffer.initialize<float>(*this, energyParamDerivNames.size()*energySize, "energyParamDerivBuffer");
        addAutoclearBuffer(energyParamDerivBuffer);
    }
    resizePinnedBuffer(max(force.getSize()*force.getElementSize(), energyBuffer.getSize()*energyBuffer.getElementSize()));
    findMoleculeGroups();
    nb.initialize(system);
    initialized = true;
}

void MetalContext::initializeUtilityKernels() {
    if (reduceEnergyKernel)
        return;
    // Reuse only the required OpenCL utilities. Native accuracy calibration and
    // the kernel-to-kernel clear calls are not needed by the Metal blit path.
    const string& utilities = MetalOpenCLKernelSources::utilities;
    size_t begin = utilities.find("__kernel void reduceEnergy(");
    size_t end = utilities.find("/**", begin);
    string source = utilities.substr(begin, end-begin);
    source += utilities.substr(utilities.find("__kernel void setCharges("));
    ComputeProgram program = compileProgram(source);
    reduceEnergyKernel = program->createKernel("reduceEnergy");
    for (int i = 0; i < 4; i++)
        reduceEnergyKernel->addArg();
    static_cast<MetalKernel&>(*reduceEnergyKernel).addLocalArg(512*sizeof(float));
    setChargesKernel = program->createKernel("setCharges");
    for (int i = 0; i < 4; i++)
        setChargesKernel->addArg();
    energySum.initialize<float>(*this, getNumThreadBlocks(), "energySum");
}

double MetalContext::reduceEnergy() {
    initializeUtilityKernels();
    int blockSize = min(512, reduceEnergyKernel->getMaxBlockSize());
    reduceEnergyKernel->setArg(0, energyBuffer);
    reduceEnergyKernel->setArg(1, energySum);
    reduceEnergyKernel->setArg(2, (int) energyBuffer.getSize());
    reduceEnergyKernel->setArg(3, blockSize);
    static_cast<MetalKernel&>(*reduceEnergyKernel).setLocalArg(4, blockSize*sizeof(float));
    reduceEnergyKernel->execute(blockSize*energySum.getSize(), blockSize);
    vector<float> partials;
    energySum.download(partials);
    return accumulate(partials.begin(), partials.end(), 0.0);
}

mm_float4 MetalContext::getPeriodicBoxSize() const {
    return mm_float4(periodicBoxVectors[0][0], periodicBoxVectors[1][1], periodicBoxVectors[2][2], 0);
}

mm_float4 MetalContext::getInvPeriodicBoxSize() const {
    mm_float4 box = getPeriodicBoxSize();
    return mm_float4(1/box.x, 1/box.y, 1/box.z, 0);
}

mm_float4 MetalContext::getPeriodicBoxVecX() const {
    const Vec3& v = periodicBoxVectors[0];
    return mm_float4(v[0], v[1], v[2], 0);
}

mm_float4 MetalContext::getPeriodicBoxVecY() const {
    const Vec3& v = periodicBoxVectors[1];
    return mm_float4(v[0], v[1], v[2], 0);
}

mm_float4 MetalContext::getPeriodicBoxVecZ() const {
    const Vec3& v = periodicBoxVectors[2];
    return mm_float4(v[0], v[1], v[2], 0);
}
