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
#include "MetalSourceAdapter.h"
#include "CommonKernelSources.h"
#include "openmm/System.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

using namespace OpenMM;
using namespace std;

/** Verify each independent switch actually selects its exact template function. */
void testSelection() {
    for (bool floating : {false, true}) {
        string groups = MetalSourceAdapter::translate(CommonKernelSources::customNonbondedGroups, floating);
        ASSERT_EQUAL(bool(OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE), groups.find("simd_shuffle_xor(") != string::npos);
        for (const string& source : {groups}) {
            ASSERT(source.find("#define __CUDA_ARCH__") == string::npos);
            ASSERT(source.find("#define USE_HIP") == string::npos);
        }
    }
    // A user helper with the same name must not acquire an unrelated fast path.
    string unrelated = "DEVICE int reduceMax(int val, LOCAL_ARG int* temp) { return val; }\n"
        "KERNEL void probe(GLOBAL int* result) { result[GLOBAL_ID] = 0; }\n";
    ASSERT(MetalSourceAdapter::translate(unrelated).find("simd_shuffle_xor(") == string::npos);
    unrelated = "DEVICE void atomicAddMixed(GLOBAL mixed* target, mixed value) { *target = value; }\n"
        "KERNEL void probe(GLOBAL int* result) { result[GLOBAL_ID] = 0; }\n";
    ASSERT(MetalSourceAdapter::translate(unrelated).find("metalAtomicAdd(target, value)") == string::npos);
}

/** Stress both float atomic primitives and their fetch-add return contract. */
void testFloatAtomics() {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    const string& minimize = CommonKernelSources::minimize;
    size_t begin = minimize.find("DEVICE void atomicAddMixed(");
    size_t end = minimize.find("KERNEL void recordInitialPos(", begin);
    ASSERT(begin != string::npos && end != string::npos);
    const string source = minimize.substr(begin, end-begin)+R"(
KERNEL void accumulate(GLOBAL float* sums, GLOBAL float* previous, int count) {
    for (int i = GLOBAL_ID; i < count; i += GLOBAL_SIZE) {
        atomicAddMixed(sums, 1.0f);
        ATOMIC_ADD(sums+1, -1.0f);
        previous[i] = ATOMIC_ADD(sums+2, 1.0f);
    }
}

)";
    const int count = 10001;
    ComputeArray sums, previous, mode;
    sums.initialize<float>(context, 3, "floatAtomicSums");
    previous.initialize<float>(context, count, "floatAtomicPrevious");
    mode.initialize<int>(context, 1, "configuredFloatAtomicMode");
    ComputeKernel configured = context.compileProgram(
        "KERNEL void configuredMode(GLOBAL int* output) { if (GLOBAL_ID == 0) output[0] = OPENMM_METAL_NATIVE_FLOAT_ATOMICS; }")
        ->createKernel("configuredMode");
    configured->addArg(mode);
    configured->execute(32, 32);
    vector<int> selected;
    mode.download(selected);
    ASSERT_EQUAL(OPENMM_METAL_NATIVE_FLOAT_ATOMICS, selected[0]);
    for (int native = 0; native <= 1; native++) {
        // Per-program overrides verify both primitives even in the OFF build.
        map<string, string> defines;
        defines["OPENMM_METAL_NATIVE_FLOAT_ATOMICS"] = to_string(native);
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("accumulate");
        kernel->addArg(sums);
        kernel->addArg(previous);
        kernel->addArg(count);
        context.clearBuffer(sums);
        kernel->execute(count);
        vector<float> values;
        sums.download(values);
        ASSERT_EQUAL(float(count), values[0]);
        ASSERT_EQUAL(-float(count), values[1]);
        ASSERT_EQUAL(float(count), values[2]);
        previous.download(values);
        sort(values.begin(), values.end());
        for (int i = 0; i < count; i++) ASSERT_EQUAL(float(i), values[i]);
    }
}

/** Verify the checked conversion and its hidden Common buffer binding on the GPU. */
void testFixedPointRangeDiagnostic() {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    const string source = R"(
KERNEL void diagnosticState(GLOBAL uint* state, int action) {
    if (GLOBAL_ID != 0) return;
#if OPENMM_METAL_CHECK_FIXED_POINT_RANGE
    if (action == 1) {
        atomic_store_explicit(_metal.fixedPointRange, 1u, memory_order_relaxed);
        atomic_store_explicit(_metal.fixedPointRange+1, 0u, memory_order_relaxed);
    }
    if (action == 2)
        atomic_store_explicit(_metal.fixedPointRange, 0u, memory_order_relaxed);
    state[0] = 1;
    state[1] = atomic_load_explicit(_metal.fixedPointRange+1, memory_order_relaxed);
#else
    state[0] = 0;
    state[1] = 0;
#endif
}
KERNEL void convertChecked(GLOBAL const float* input, GLOBAL mm_long* output) {
    if (GLOBAL_ID < 8) output[GLOBAL_ID] = realToFixedPoint(input[GLOBAL_ID]);
}
)";
    ComputeArray state, input, output;
    state.initialize<unsigned int>(context, 2, "fixedPointDiagnosticState");
    input.initialize<float>(context, 8, "fixedPointDiagnosticInput");
    output.initialize<int64_t>(context, 8, "fixedPointDiagnosticOutput");
    ComputeProgram program = context.compileProgram(source);
    ComputeKernel inspect = program->createKernel("diagnosticState");
    inspect->addArg(state);
    inspect->addArg(1);
    inspect->execute(32, 32);
    vector<unsigned int> flags;
    state.download(flags);
    if (flags[0] == 0) return; // Float-minimization builds do not require this diagnostic ABI.
    ASSERT_EQUAL(0u, flags[1]);
    const float limit = 2147483648.0f;
    input.upload(vector<float>{1.25f, -3.5f, limit, -limit,
        numeric_limits<float>::infinity(), -numeric_limits<float>::infinity(),
        numeric_limits<float>::quiet_NaN(), nextafter(limit, 0.0f)});
    ComputeKernel convert = program->createKernel("convertChecked");
    convert->addArg(input);
    convert->addArg(output);
    convert->execute(32, 32);
    vector<int64_t> values;
    output.download(values);
    ASSERT_EQUAL(int64_t(5368709120), values[0]);
    ASSERT_EQUAL(int64_t(-15032385536), values[1]);
    for (int i = 2; i <= 6; i++) ASSERT_EQUAL(int64_t(0), values[i]);
    ASSERT_EQUAL(int64_t(2147483520)*int64_t(4294967296), values[7]);
    inspect->setArg(1, 0);
    inspect->execute(32, 32);
    state.download(flags);
    ASSERT_EQUAL(1u, flags[1]);
    inspect->setArg(1, 2);
    inspect->execute(32, 32);
}

/** Exercise the exact MSL3 xor-shuffle reduction across two SIMD groups. */
void testShuffle() {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    const string source = R"(
#include <metal_stdlib>
using namespace metal;
kernel void shuffleMaximum(device uint* output [[buffer(0)]], uint gid [[thread_position_in_grid]]) {
    uint maximum = gid;
    for (int mask = 16; mask > 0; mask /= 2)
        maximum = max(maximum, simd_shuffle_xor(maximum, mask));
    output[gid] = maximum;
}
)";
    ComputeArray output;
    output.initialize<unsigned int>(context, 64, "shuffleMaximumOutput");
    ComputeKernel kernel = context.compileProgram(source)->createKernel("shuffleMaximum");
    kernel->addArg(output);
    kernel->execute(64, 64);
    vector<unsigned int> result;
    output.download(result);
    for (int i = 0; i < 64; i++)
        ASSERT_EQUAL(32u*(i/32)+31u, result[i]);
}

int main() {
    try {
        testSelection();
        testFloatAtomics();
        testFixedPointRangeDiagnostic();
        testShuffle();
    }
    catch (const exception& error) {
        if (string(error.what()).find("No Metal device") != string::npos)
            return 77;
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Metal source adapter fast-path tests passed" << endl;
    return 0;
}
