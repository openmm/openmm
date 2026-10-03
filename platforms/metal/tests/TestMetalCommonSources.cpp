/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalKernel.h"
#include "MetalQueue.h"
#include "CommonKernelSources.h"
#include "openmm/System.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>

using namespace OpenMM;
using namespace std;

/** Execute an unchanged Common coordinate-copy kernel. */
void testCoordinates(MetalContext& context) {
    ComputeArray input, output;
    input.initialize<float>(context, 9, "coordinateInput");
    output.initialize<mm_float4>(context, 3, "coordinateOutput");
    vector<float> values{1,2,3,4,5,6,7,8,9};
    input.upload(values);
    context.clearBuffer(output);
    ComputeKernel kernel = context.compileProgram(CommonKernelSources::copyCoordinateBuffers)->createKernel("copyFloatBuffer");
    kernel->addArg(input);
    kernel->addArg(output);
    kernel->addArg(3);
    kernel->execute(3);
    vector<mm_float4> result;
    output.download(result);
    for (int i = 0; i < 3; i++) {
        ASSERT_EQUAL(values[3*i], result[i].x);
        ASSERT_EQUAL(values[3*i+1], result[i].y);
        ASSERT_EQUAL(values[3*i+2], result[i].z);
        ASSERT_EQUAL(0, result[i].w);
    }
}

/** Preserve preprocessor branches and Common's argument order beyond 31 slots. */
void testArgumentBuffer(MetalContext& context) {
    string source = "KERNEL void manyArgs(GLOBAL int* output\n";
    for (int i = 0; i < 40; i++) source += ", int a"+to_string(i)+"\n";
    source += "#ifdef ENABLE_EXTRA\n, int extra\n#endif\n, int tail) { int sum = tail;\n";
    for (int i = 0; i < 40; i++) source += "sum += a"+to_string(i)+";\n";
    source += "#ifdef ENABLE_EXTRA\nsum += extra;\n#endif\nif (GLOBAL_ID == 0) output[0] = sum; }";
    ComputeArray output;
    output.initialize<int>(context, 1, "manyArgumentsOutput");
    for (int extra = 0; extra < 2; extra++) {
        map<string,string> defines;
        if (extra) defines["ENABLE_EXTRA"] = "1";
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("manyArgs");
        kernel->addArg(output);
        for (int i = 0; i < 40; i++) kernel->addArg(i);
        if (extra) kernel->addArg(101);
        kernel->addArg(5);
        kernel->execute(1);
        vector<int> result;
        output.download(result);
        ASSERT_EQUAL(785+101*extra, result[0]);
    }
}

/** Reused argument storage must preserve each queued scalar and array binding. */
void testArgumentSnapshots(MetalContext& context) {
    const string source = "KERNEL void snapshot(GLOBAL int* output, int index, int value) { if (GLOBAL_ID == 0) output[index] = value; }";
    const int count = 130;
    ComputeArray first, second;
    first.initialize<int>(context, count, "firstArgumentSnapshots");
    second.initialize<int>(context, count, "secondArgumentSnapshots");
    ComputeKernel kernel = context.compileProgram(source)->createKernel("snapshot");
    kernel->addArg(first);
    kernel->addArg(0);
    kernel->addArg(0);
    for (int wave = 0; wave < 3; wave++) {
        for (int i = 0; i < count; i++) {
            kernel->setArg(1, i);
            kernel->setArg(0, first);
            kernel->setArg(2, i*7-19+wave);
            kernel->execute(1);
            kernel->setArg(0, second);
            kernel->setArg(2, 31-i*5-wave);
            kernel->execute(1);
        }
        vector<int> result;
        first.download(result);
        for (int i = 0; i < count; i++) ASSERT_EQUAL(i*7-19+wave, result[i]);
        second.download(result);
        for (int i = 0; i < count; i++) ASSERT_EQUAL(31-i*5-wave, result[i]);
    }

    // A rejected binding must not strand or corrupt a reusable snapshot.
    kernel->setArg(2, int64_t(5));
    bool rejected = false;
    try {
        kernel->execute(1);
    }
    catch (const OpenMMException& error) {
        rejected = string(error.what()).find("Wrong primitive size") != string::npos;
    }
    ASSERT(rejected);
    kernel->setArg(1, 0);
    kernel->setArg(2, 73);
    kernel->execute(1);
    // The submitted command must retain its argument buffer after kernel destruction.
    kernel.reset();
    vector<int> result;
    second.download(result);
    ASSERT_EQUAL(73, result[0]);
}

/** Completion on one queue cannot release argument snapshots used by another. */
void testArgumentSnapshotsAcrossQueues(MetalContext& context) {
    const string source = "KERNEL void queueSnapshot(GLOBAL int* output, int index, int value) { if (GLOBAL_ID == 0) output[index] = value; }";
    ComputeArray first, second;
    first.initialize<int>(context, 130, "firstQueueSnapshots");
    second.initialize<int>(context, 130, "secondQueueSnapshots");
    ComputeQueue original = context.getCurrentQueue(), secondary = context.createQueue();
    ComputeKernel kernel = context.compileProgram(source)->createKernel("queueSnapshot");
    kernel->addArg(first);
    kernel->addArg(0);
    kernel->addArg(0);
    for (int wave = 0; wave < 3; wave++) {
        for (int i = 0; i < 130; i++) {
            context.setCurrentQueue(original);
            kernel->setArg(0, first);
            kernel->setArg(1, i);
            kernel->setArg(2, i+wave*17);
            kernel->execute(1);
            context.setCurrentQueue(secondary);
            kernel->setArg(0, second);
            kernel->setArg(2, 23-i-wave);
            kernel->execute(1);
        }
        vector<int> result;
        second.download(result);
        for (int i = 0; i < 130; i++) ASSERT_EQUAL(23-i-wave, result[i]);
        context.setCurrentQueue(original);
        first.download(result);
        for (int i = 0; i < 130; i++) ASSERT_EQUAL(i+wave*17, result[i]);
    }
    // Finish both queues before the output arrays leave scope.
    static_cast<MetalQueue&>(*secondary).finish();
    context.setCurrentQueue(original);
}

/** Exercise helper built-ins, private pointers, vector constructors and LOCAL_ARG. */
void testHelpersAndLocalMemory(MetalContext& context) {
    string source = R"(
DEVICE int helper(int* value) { return *value+LOCAL_ID; }
DEVICE int nested(int value) { return helper(&value); }
__kernel void scratch(__global int* output, __local int* values) {
    const int thread = get_local_id(0);
    int2 pair = (int2) (nested(2), 3);
    values[thread] = pair.x+pair.y;
    barrier(CLK_LOCAL_MEM_FENCE);
    if (thread == 0) output[0] = values[0]+values[31];
}
)";
    ComputeArray output;
    output.initialize<int>(context, 1, "localMemoryOutput");
    ComputeKernel kernel = context.compileProgram(source)->createKernel("scratch");
    kernel->addArg(output);
    static_cast<MetalKernel&>(*kernel).addLocalArg(32*sizeof(int));
    kernel->execute(32, 32);
    vector<int> result;
    output.download(result);
    ASSERT_EQUAL(41, result[0]);
}

/** Reject the aggregate dynamic allocation before starting a command encoder. */
void testLocalMemoryLimit(MetalContext& context) {
    const string source = R"(
KERNEL void localLimit(GLOBAL int* output, LOCAL_ARG int* first, LOCAL_ARG int* second) {
    first[LOCAL_ID] = LOCAL_ID;
    second[LOCAL_ID] = LOCAL_ID+1;
    SYNC_THREADS;
    if (LOCAL_ID == 0) output[0] = first[1]+second[1];
}
)";
    ComputeArray output;
    output.initialize<int>(context, 1, "localLimitOutput");
    ComputeKernel kernel = context.compileProgram(source)->createKernel("localLimit");
    kernel->addArg(output);
    MetalKernel& metalKernel = static_cast<MetalKernel&>(*kernel);
    size_t allocation = context.getMaxThreadgroupMemory()/2+16;
    metalKernel.addLocalArg(allocation);
    metalKernel.addLocalArg(allocation);
    bool rejected = false;
    try {
        kernel->execute(32, 32);
    }
    catch (const OpenMMException& error) {
        rejected = string(error.what()).find("Total threadgroup memory") != string::npos;
    }
    ASSERT(rejected);
    // A validation failure must leave the queue usable for the next dispatch.
    metalKernel.setLocalArg(1, 32*sizeof(int));
    metalKernel.setLocalArg(2, 32*sizeof(int));
    kernel->execute(32, 32);
    vector<int> result;
    output.download(result);
    ASSERT_EQUAL(3, result[0]);
}

/** Preserve mutually exclusive braces, helper macro calls and null arguments. */
void testConditionalSource(MetalContext& context) {
    const string source = R"(
#define READ_VALUE(value) readValue(value)
typedef struct { int value; } Value;
DEVICE int readValue(Value* value) { return value->value+LOCAL_ID; }
DEVICE int matrixValue(int (*matrix)[2]) {
    int (*row)[2] = matrix;
    return row[0][1];
}
KERNEL void conditionalSource(GLOBAL int* output, GLOBAL int* unused) {
    Value value = {5};
    int matrix[1][2] = {{1, 2}};
#ifdef ALTERNATIVE
    if (GLOBAL_ID == 0) {
#else
    if (LOCAL_ID == 0) {
#endif
        output[0] = READ_VALUE(&value)+matrixValue(matrix);
    }
}
)";
    ComputeArray output;
    output.initialize<int>(context, 1, "conditionalSourceOutput");
    for (int branch = 0; branch < 2; branch++) {
        map<string, string> defines;
        if (branch) defines["ALTERNATIVE"] = "1";
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("conditionalSource");
        kernel->addArg(output);
        kernel->addArg(nullptr);
        kernel->execute(32, 32);
        vector<int> result;
        output.download(result);
        ASSERT_EQUAL(7, result[0]);
    }
}

/** Validate carry and signed contributions under concurrent 64-bit reductions. */
void testAtomicReduction(MetalContext& context) {
    const string source = R"(
KERNEL void accumulate(GLOBAL mm_ulong* sums, int count) {
    for (int i = GLOBAL_ID; i < count; i += GLOBAL_SIZE) {
        ATOMIC_ADD(sums, (mm_ulong) realToFixedPoint(1.25f));
        ATOMIC_ADD(sums+1, (mm_ulong) realToFixedPoint(-0.75f));
    }
}
)";
    ComputeArray sums;
    sums.initialize<int64_t>(context, 2, "atomicSums");
    context.clearBuffer(sums);
    ComputeKernel kernel = context.compileProgram(source)->createKernel("accumulate");
    kernel->addArg(sums);
    kernel->addArg(10001);
    kernel->execute(10001);
    vector<int64_t> result;
    sums.download(result);
    ASSERT_EQUAL(int64_t(10001)*5368709120LL, result[0]);
    ASSERT_EQUAL(int64_t(10001)*-3221225472LL, result[1]);
}

/** Execute unchanged Common Verlet kernels in single precision. */
void testVerlet(MetalContext& context) {
    ComputeArray dt, position, velocity, force, delta;
    dt.initialize<mm_float2>(context, 1, "verletDt");
    position.initialize<mm_float4>(context, 1, "verletPosition");
    velocity.initialize<mm_float4>(context, 1, "verletVelocity");
    force.initialize<int64_t>(context, 3, "verletForce");
    delta.initialize<mm_float4>(context, 1, "verletDelta");
    dt.upload(vector<mm_float2>{mm_float2(0.25f, 0.25f)});
    position.upload(vector<mm_float4>{mm_float4(1,2,3,0)});
    velocity.upload(vector<mm_float4>{mm_float4(1,0,0,1)});
    force.upload(vector<int64_t>{int64_t(1)<<32,0,0});
    ComputeProgram program = context.compileProgram(CommonKernelSources::verlet);
    ComputeKernel first = program->createKernel("integrateVerletPart1");
    first->addArg(1); first->addArg(1); first->addArg(dt); first->addArg(position);
    first->addArg(velocity); first->addArg(force); first->addArg(delta);
    first->execute(1);
    ComputeKernel second = program->createKernel("integrateVerletPart2");
    second->addArg(1); second->addArg(dt); second->addArg(position);
    second->addArg(velocity); second->addArg(delta);
    second->execute(1);
    vector<mm_float4> result;
    position.download(result);
    ASSERT_EQUAL(1.3125f, result[0].x);
    ASSERT_EQUAL(2, result[0].y);
    ASSERT_EQUAL(3, result[0].z);
    velocity.download(result);
    ASSERT_EQUAL(1.25f, result[0].x);
}

/** Compile every integration helper and execute its unchanged energy reduction. */
void testIntegrationUtilities(MetalContext& context) {
    map<string,string> defines;
    defines["NUM_ATOMS"] = "1";
    defines["PADDED_NUM_ATOMS"] = "32";
    defines["NUM_CCMA_ATOMS"] = "1";
    defines["NUM_CCMA_CONSTRAINTS"] = "1";
    defines["KE_WORK_GROUP_SIZE"] = "64";
    defines["NUM_2_AVERAGE"] = "1";
    defines["NUM_3_AVERAGE"] = "1";
    defines["NUM_OUT_OF_PLANE"] = "1";
    defines["NUM_LOCAL_COORDS"] = "1";
    defines["NUM_SYMMETRY"] = "1";
    ComputeProgram program = context.compileProgram(CommonKernelSources::integrationUtilities, defines);
    ComputeArray velocity, energy;
    velocity.initialize<mm_float4>(context, 1, "integrationVelocity");
    energy.initialize<float>(context, 1, "integrationEnergy");
    velocity.upload(vector<mm_float4>{mm_float4(1,2,3,0.5f)});
    ComputeKernel kernel = program->createKernel("computeKineticEnergy");
    kernel->addArg(velocity);
    kernel->addArg(energy);
    kernel->execute(64, 64);
    vector<float> result;
    energy.download(result);
    ASSERT_EQUAL(14, result[0]);
}

/** Common's double-float charge reduction must retain low-order cancellation. */
void testCompensatedChargeReduction(MetalContext& context) {
    map<string, string> defines;
    defines["NUM_ELECTRODE_PARTICLES"] = "4";
    defines["THREAD_BLOCK_SIZE"] = "32";
    defines["THREAD_BLOCK_COUNT"] = "1";
    defines["ERROR_TARGET"] = "1e-6f";
    defines["USE_CHARGE_CONSTRAINT"] = "1";
    ComputeProgram program = context.compileProgram(CommonKernelSources::constantPotentialCGSolver, defines);
    ComputeArray charges, previous;
    charges.initialize<float>(context, 4, "compensatedCharges");
    previous.initialize<float>(context, 4, "compensatedPreviousCharges");
    vector<float> values{16777216.0f, 1.0f, -16777216.0f, 0.0f};
    charges.upload(values);
    previous.upload(values);
    ComputeKernel kernel = program->createKernel("solveInitializeStep1");
    kernel->addArg(charges);
    kernel->addArg(previous);
    kernel->addArg(0.0f);
    kernel->execute(32, 32);
    vector<float> result;
    charges.download(result);
    ASSERT_EQUAL(0.75f, result[1]);
    ASSERT_EQUAL(-0.25f, result[3]);
}

/** The private ABI keeps large forces floating, but tile indices stay 64-bit. */
void testFloatingAccumulators() {
    System system;
    system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, true);
    ASSERT_EQUAL(sizeof(int64_t), context.getLongForceBuffer().getElementSize());
    const string source = R"(
KERNEL void accumulateLarge(GLOBAL mm_ulong* forceBuffers, mm_long numTiles) {
    if (GLOBAL_ID == 0) ATOMIC_ADD(forceBuffers, (mm_ulong) realToFixedPoint(1.0e20f));
    if (GLOBAL_ID == 1) ATOMIC_ADD(forceBuffers, (mm_ulong) realToFixedPoint(-1.0e19f));
    if (GLOBAL_ID == 2) forceBuffers[1] = (numTiles >> 32);
}
KERNEL void decodeLarge(GLOBAL const mm_long* longForces, GLOBAL float* output) {
    real scale = 1/(real) 0x100000000;
    if (GLOBAL_ID == 0) output[0] = scale*longForces[0];
    if (GLOBAL_ID == 1) output[1] = scale*longForces[1];
}
)";
    ComputeArray sums, output;
    sums.initialize<float>(context, 2, "largeForceSums");
    output.initialize<float>(context, 2, "largeForceOutput");
    context.clearBuffer(sums);
    ComputeProgram program = context.compileProgram(source);
    ComputeKernel accumulate = program->createKernel("accumulateLarge");
    accumulate->addArg(sums);
    accumulate->addArg(int64_t(7)<<32);
    accumulate->execute(32, 32);
    ComputeKernel decode = program->createKernel("decodeLarge");
    decode->addArg(sums);
    decode->addArg(output);
    decode->execute(32, 32);
    vector<float> result;
    output.download(result);
    ASSERT_EQUAL_TOL(9.0e19f, result[0], 1e-6);
    ASSERT_EQUAL(7.0f, result[1]);
    ComputeArray malformed;
    malformed.initialize<short>(context, 3, "invalidAccumulatorElementSize");
    bool threw = false;
    try {
        vector<double> invalid;
        context.downloadFixedPointBuffer(malformed, invalid);
    }
    catch (const OpenMMException&) {
        threw = true;
    }
    ASSERT(threw);
}

/** Cached kernels keep their bindings and original integer ABI across linked mode switches. */
void testAccumulatorModeSwitching(bool initiallyFloating) {
    System system;
    system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, initiallyFloating);
    MetalContext linked(system, nullptr, &context);
    const string source = R"(
KERNEL void accumulate(GLOBAL mm_ulong* forceBuffers, real value, mm_long tile) {
    if (GLOBAL_ID == 0) ATOMIC_ADD(forceBuffers, (mm_ulong) realToFixedPoint(value));
    if (GLOBAL_ID == 1) forceBuffers[1] = realToFixedPoint((real) (tile >> 32));
}
KERNEL void decode(GLOBAL const mm_long* forceBuffers, GLOBAL real* result) {
    if (GLOBAL_ID < 2) result[GLOBAL_ID] = forceBuffers[GLOBAL_ID]/(real) 0x100000000;
}
)";
    ComputeArray result;
    result.initialize<float>(context, 2, "switchedForceResult");
    // Kernels deliberately outlive the programs that created them.
    ComputeKernel accumulate = context.compileProgram(source)->createKernel("accumulate");
    ComputeKernel decode = linked.compileProgram(source)->createKernel("decode");
    accumulate->addArg(context.getLongForceBuffer());
    accumulate->addArg(0.0f);
    accumulate->addArg(int64_t(7)<<32);
    decode->addArg(context.getLongForceBuffer());
    decode->addArg(result);
    for (bool floating : {false, true, false, true, false}) {
        // Enter through either member: the group must switch as a whole.
        linked.setUseFloatingPointAccumulators(floating);
        ASSERT_EQUAL(floating, context.getUseFloatingPointAccumulators());
        ASSERT_EQUAL(floating, linked.getUseFloatingPointAccumulators());
        const float value = floating ? 1.0e20f : 1.25f;
        context.clearBuffer(context.getLongForceBuffer());
        accumulate->setArg(1, value);
        accumulate->execute(32, 32);
        decode->execute(32, 32);
        vector<float> values;
        result.download(values);
        ASSERT_EQUAL_TOL(value, values[0], 1e-6);
        ASSERT_EQUAL(7.0f, values[1]);
    }
}

/** Revalidate cached cross-context bindings and decoders after an independent mode change. */
void testIndependentAccumulatorModes(bool initiallyFloating) {
    System system;
    system.addParticle(1);
    MetalContext first(system, nullptr, nullptr, initiallyFloating);
    MetalContext second(system, nullptr, nullptr, initiallyFloating);
    const string source = R"(
KERNEL void accumulateIndependent(GLOBAL mm_ulong* forceBuffers, real value) {
    if (GLOBAL_ID == 0) ATOMIC_ADD(forceBuffers, (mm_ulong) realToFixedPoint(value));
}
)";
    ComputeKernel firstKernel = first.compileProgram(source)->createKernel("accumulateIndependent");
    ComputeKernel secondKernel = second.compileProgram(source)->createKernel("accumulateIndependent");
    // Bind while the independent contexts agree, then keep these bindings unchanged.
    firstKernel->addArg(second.getLongForceBuffer());
    firstKernel->addArg(1.25f);
    secondKernel->addArg(first.getLongForceBuffer());
    secondKernel->addArg(-2.5f);
    auto verifySharedBuffers = [&] {
        first.clearBuffer(first.getLongForceBuffer());
        second.clearBuffer(second.getLongForceBuffer());
        // Independent queues require explicit ordering before cross-context access.
        first.flushQueue();
        second.flushQueue();
        firstKernel->execute(32, 32);
        first.flushQueue();
        secondKernel->execute(32, 32);
        second.flushQueue();
        vector<double> values;
        first.downloadFixedPointBuffer(second.getLongForceBuffer(), values);
        ASSERT_EQUAL(1.25, values[0]);
        for (size_t i = 1; i < values.size(); i++) ASSERT_EQUAL(0.0, values[i]);
        second.downloadFixedPointBuffer(first.getLongForceBuffer(), values);
        ASSERT_EQUAL(-2.5, values[0]);
        for (size_t i = 1; i < values.size(); i++) ASSERT_EQUAL(0.0, values[i]);
    };
    auto expectModeMismatch = [](const function<void()>& operation) {
        bool rejected = false;
        try {
            operation();
        }
        catch (const OpenMMException& error) {
            rejected = string(error.what()).find("Cannot share accumulator arrays") != string::npos;
        }
        ASSERT(rejected);
    };
    verifySharedBuffers();
    first.setUseFloatingPointAccumulators(!initiallyFloating);
    ASSERT_EQUAL(!initiallyFloating, first.getUseFloatingPointAccumulators());
    ASSERT_EQUAL(initiallyFloating, second.getUseFloatingPointAccumulators());
    expectModeMismatch([&] { firstKernel->execute(32, 32); });
    expectModeMismatch([&] { secondKernel->execute(32, 32); });
    vector<double> values;
    expectModeMismatch([&] { first.downloadFixedPointBuffer(second.getLongForceBuffer(), values); });
    expectModeMismatch([&] { second.downloadFixedPointBuffer(first.getLongForceBuffer(), values); });
    first.setUseFloatingPointAccumulators(initiallyFloating);
    // Validation failures must leave both queues and the original bindings reusable.
    verifySharedBuffers();
}

/** Zero-weight ATM states must not propagate an inactive state's Inf/NaN force. */
void testFloatingInactiveATMState() {
    System system;
    system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, true);
    ComputeArray force, first, second, displacement, order;
    force.initialize<float>(context, 3, "atmHybridForce");
    first.initialize<float>(context, 3, "atmFirstForce");
    second.initialize<float>(context, 3, "atmSecondForce");
    displacement.initialize<float>(context, 3, "atmDisplacementForce");
    order.initialize<int>(context, 1, "atmAtomOrder");
    context.clearBuffer(displacement);
    order.upload(vector<int>{0});
    ComputeKernel kernel = context.compileProgram(CommonKernelSources::atmforce)->createKernel("hybridForce");
    kernel->addArg(1); kernel->addArg(1);
    kernel->addArg(force); kernel->addArg(first); kernel->addArg(second);
    kernel->addArg(displacement); kernel->addArg(displacement);
    kernel->addArg(order); kernel->addArg(order); kernel->addArg(order);
    kernel->addArg(1.0f); kernel->addArg(0.0f);
    vector<float> active{2.0f, -3.0f, 4.0f};
    vector<float> inactive{numeric_limits<float>::infinity(), -numeric_limits<float>::infinity(), numeric_limits<float>::quiet_NaN()};
    for (int side = 0; side < 2; side++) {
        first.upload(side == 0 ? active : inactive);
        second.upload(side == 0 ? inactive : active);
        kernel->setArg(10, side == 0 ? 1.0f : 0.0f);
        kernel->setArg(11, side == 0 ? 0.0f : 1.0f);
        context.clearBuffer(force);
        kernel->execute(32, 32);
        vector<float> result;
        force.download(result);
        ASSERT_EQUAL_CONTAINERS(active, result);
    }
}

/** The selected per-function native policy must still meet OpenCL's tolerance. */
void testNativeMathPolicy(MetalContext& context) {
    const string source = R"(
KERNEL void evaluateMath(GLOBAL const float* input, GLOBAL float* output) {
    int i = GLOBAL_ID;
    if (i >= NUM_VALUES) return;
    float v = input[i];
    output[5*i] = SQRT(v);
    output[5*i+1] = RSQRT(v);
    output[5*i+2] = RECIP(v);
    output[5*i+3] = EXP(v);
    output[5*i+4] = LOG(v);
}
)";
    const int count = 40;
    ComputeArray input, output;
    input.initialize<float>(context, count, "mathInput");
    output.initialize<float>(context, 5*count, "mathOutput");
    vector<float> values(count), result;
    float nextValue = 1e-4f;
    for (int i = 0; i < count/2; i++) {
        values[i] = 0.01f+0.1f*i;
        values[count/2+i] = nextValue;
        nextValue *= (float) M_PI;
    }
    input.upload(values);
    for (int strict = 0; strict < 2; strict++) {
        map<string, string> defines{{"NUM_VALUES", to_string(count)}};
        if (strict)
            defines["OPENMM_METAL_REQUIRE_SAFE_MATH"] = "1";
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("evaluateMath");
        kernel->addArg(input);
        kernel->addArg(output);
        kernel->execute(count);
        output.download(result);
        for (int i = 0; i < count; i++) {
            double v = values[i];
            double expected[] = {sqrt(v), 1/sqrt(v), 1/v, exp(v), log(v)};
            for (int j = 0; j < 5; j++) {
                if (expected[j] > numeric_limits<float>::max()) {
                    ASSERT_EQUAL(numeric_limits<float>::infinity(), result[5*i+j]);
                }
                else {
                    ASSERT(isfinite(result[5*i+j]));
                    ASSERT(fabs((result[5*i+j]-expected[j])/expected[j]) < 1e-6);
                }
            }
        }
    }
}

int main() {
    try {
        System system;
        system.addParticle(1);
        MetalContext context(system);
        testCoordinates(context);
        testNativeMathPolicy(context);
        testArgumentBuffer(context);
        testArgumentSnapshots(context);
        testArgumentSnapshotsAcrossQueues(context);
        testHelpersAndLocalMemory(context);
        testLocalMemoryLimit(context);
        testConditionalSource(context);
        testAtomicReduction(context);
        testVerlet(context);
        testIntegrationUtilities(context);
        testCompensatedChargeReduction(context);
        testFloatingAccumulators();
        testAccumulatorModeSwitching(false);
        testAccumulatorModeSwitching(true);
        testIndependentAccumulatorModes(false);
        testIndependentAccumulatorModes(true);
        testFloatingInactiveATMState();
    }
    catch (const exception& error) {
        if (string(error.what()).find("No Metal device") != string::npos) {
            cout << "No Metal device; Common source GPU tests skipped" << endl;
            return 77;
        }
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Metal Common source tests passed" << endl;
    return 0;
}
