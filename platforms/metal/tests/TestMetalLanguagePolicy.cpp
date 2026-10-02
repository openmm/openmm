/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform tests: copyright (c) 2026 Chun-Chi Hung.                      *
 * Author: Chun-Chi Hung.                                                      *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalKernelSources.h"
#include "../src/MetalLanguagePolicy.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace OpenMM;
using namespace std;

/** Verify capability limits independently of this machine's installed SDK and OS. */
static void testSelection() {
    const int versions[] = {300, 310, 320, 400, 410};
    for (int requested : versions)
        for (int sdk : versions)
            for (int os : versions)
                for (bool apple : {false, true}) {
                    const int actual = MetalLanguagePolicy::selectVersion(requested, sdk, os, apple);
                    ASSERT_EQUAL(min(min(requested, sdk), min(os, apple ? 410 : 310)), actual);
                }
    ASSERT_EQUAL(0, MetalLanguagePolicy::selectVersion(410, 299, 410, true));
    ASSERT_EQUAL(0, MetalLanguagePolicy::selectVersion(410, 410, 299, true));
    ASSERT_EQUAL(320, MetalLanguagePolicy::selectVersion(399, 410, 410, true));
    ASSERT_EQUAL(300, MetalLanguagePolicy::selectVersion(300, 410, 410, true));
    ASSERT_EQUAL(310, MetalLanguagePolicy::selectVersion(410, 410, 410, false));

    // Optimization switches must not influence language selection, whether
    // built all OFF, individually ON, or together ON. Only the independent
    // language override may restrict automatic capability-based selection.
#if OPENMM_METAL_TUNE_LANGUAGE_VERSION
    const int requested = OPENMM_METAL_LANGUAGE_VERSION;
#else
    const int requested = 410;
#endif
    ASSERT_EQUAL(requested, MetalLanguagePolicy::requestedVersion());
    for (int available : versions) {
        ASSERT_EQUAL(min(requested, available), MetalLanguagePolicy::selectVersion(
                MetalLanguagePolicy::requestedVersion(), available, 410, true));
        ASSERT_EQUAL(min(requested, available), MetalLanguagePolicy::selectVersion(
                MetalLanguagePolicy::requestedVersion(), 410, available, true));
    }
}

/** Verify version-specific shader choices and their fallbacks live in MSL source. */
static void testShaderSource() {
    const string& common = MetalKernelSources::common;
    ASSERT(common.find("#if OPENMM_METAL_TUNE_FORCE_REQUIRED_THREADS && defined(FORCE_WORK_GROUP_SIZE) && __METAL_VERSION__ >= 400") != string::npos);
    ASSERT(common.find("#if OPENMM_METAL_USE_TILED_ACQ_REL_BARRIERS && __METAL_VERSION__ >= 410") != string::npos);
    ASSERT(common.find("#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup, memory_order_acq_rel, thread_scope_simdgroup);") != string::npos);
    ASSERT(common.find("#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup);") != string::npos);
    ASSERT(common.find("#define SYNC_THREADS threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);") != string::npos);
    const string& halfBounds = MetalKernelSources::neighborHalfBounds;
    ASSERT(halfBounds.find("#if OPENMM_METAL_FAST_FP16_BOUNDS_NEXTAFTER && __METAL_VERSION__ >= 310") != string::npos);
    ASSERT(halfBounds.find("nextafter(result, half(INFINITY))") != string::npos);
    ASSERT(halfBounds.find("as_type<half>(ushort(bits+1))") != string::npos);
    const string& mathPolicy = MetalKernelSources::mathPolicy;
    ASSERT(mathPolicy.find("#if OPENMM_METAL_HAS_MATH_PRAGMAS") != string::npos);
    ASSERT(mathPolicy.find("#if OPENMM_METAL_USE_FAST_MATH") != string::npos);
    ASSERT(mathPolicy.find("#pragma METAL fp math_mode(fast)") != string::npos);
    ASSERT(mathPolicy.find("#pragma METAL fp math_mode(safe)") != string::npos);
}

/** Verify target and tile-local synchronization in both initial and lazy accumulator libraries. */
static void testShaderTarget(bool marked) {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    const int version = context.getMetalLanguageVersion();
    ASSERT(version >= 300 && version <= MetalLanguagePolicy::requestedVersion());
    ComputeArray output;
    output.initialize<int>(context, 64, "languagePolicyOutput");
    map<string, string> defines;
    if (marked)
        defines["OPENMM_METAL_TILED_FORCE_PROGRAM"] = "1";
    // A caller cannot force the private selection marker onto an unreviewed program.
    defines["OPENMM_METAL_USE_TILED_ACQ_REL_BARRIERS"] = "1";
    const string source = R"(
        KERNEL void checkLanguageTarget(GLOBAL int* output) {
            LOCAL int values[64];
            values[LOCAL_ID] = 100+LOCAL_ID;
            SYNC_WARPS;
            int selected = 0;
        #if OPENMM_METAL_USE_TILED_ACQ_REL_BARRIERS && __METAL_VERSION__ >= 410
            selected = 1;
        #endif
            const int next = (LOCAL_ID&~31)+((LOCAL_ID+1)&31);
            output[GLOBAL_ID] = values[next]+1000*__METAL_VERSION__+1000000*selected;
        }
    )";
    ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("checkLanguageTarget");
    kernel->addArg(output);
    int enabled = 0;
#if OPENMM_METAL_FAST_TILED_ACQ_REL_BARRIERS
    enabled = marked && version >= 410;
#endif
    for (bool floating : {false, true, false}) {
        context.setUseFloatingPointAccumulators(floating);
        kernel->execute(64, 64);
        vector<int> values;
        output.download(values);
        for (int i = 0; i < 64; i++)
            ASSERT_EQUAL(100+(i&~31)+((i+1)&31)+1000*version+1000000*enabled, values[i]);
    }
    cout << "MSL " << version << ", tiled barrier " << enabled << endl;
}

/** Verify source math pragmas preserve strict and floating-accumulator safety. */
static void testMathPolicy(bool common, bool strict) {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    const vector<float> inputValues = {0.0f, -0.0f, 1.0f, -1.0f,
            numeric_limits<float>::infinity(), -numeric_limits<float>::infinity(),
            numeric_limits<float>::quiet_NaN()};
    ComputeArray input, output;
    input.initialize<float>(context, inputValues.size(), "mathPolicyInput");
    output.initialize<int>(context, inputValues.size(), "mathPolicyOutput");
    input.upload(inputValues);

    // Native programs do not use the Common source adapter or alternate ABI.
    const string source = common ? R"(
        #if OPENMM_METAL_HAS_MATH_PRAGMAS != 0 && OPENMM_METAL_HAS_MATH_PRAGMAS != 1
        #error The host must replace caller-provided pragma capability macros.
        #endif
        KERNEL void checkMathPolicy(GLOBAL const float* input, GLOBAL int* output) {
            const int index = GLOBAL_ID;
            if (index >= 7) return;
            int classification = 0;
        #if !OPENMM_METAL_USE_FAST_MATH
            float value = input[index];
            classification = 2*int(isinf(value))+4*int(isnan(value))+8*int(isfinite(value));
        #endif
            output[index] = OPENMM_METAL_USE_FAST_MATH+classification;
        }
    )" : R"(
        #if OPENMM_METAL_HAS_MATH_PRAGMAS != 0 && OPENMM_METAL_HAS_MATH_PRAGMAS != 1
        #error The host must replace caller-provided pragma capability macros.
        #endif
        #include <metal_stdlib>
        using namespace metal;
        kernel void checkMathPolicy(device const float* input [[buffer(0)]],
                device int* output [[buffer(1)]], uint index [[thread_position_in_grid]]) {
            if (index >= 7) return;
            int classification = 0;
        #if !OPENMM_METAL_USE_FAST_MATH
            float value = input[index];
            classification = 2*int(isinf(value))+4*int(isnan(value))+8*int(isfinite(value));
        #endif
            output[index] = OPENMM_METAL_USE_FAST_MATH+classification;
        }
    )";
    map<string, string> defines;
    if (strict)
        defines["OPENMM_METAL_REQUIRE_SAFE_MATH"] = "1";
    // Callers cannot override the host's runtime capability or safety decision.
    defines["OPENMM_METAL_HAS_MATH_PRAGMAS"] = "2";
    defines["OPENMM_METAL_USE_FAST_MATH"] = "2";
    ComputeKernel kernel;
    for (bool floating : {false, true, false}) {
        context.setUseFloatingPointAccumulators(floating);
        if (!common || !kernel) {
            kernel = context.compileProgram(source, defines)->createKernel("checkMathPolicy");
            kernel->addArg(input);
            kernel->addArg(output);
        }
        kernel->execute(32, 32);
        vector<int> values;
        output.download(values);
        const bool fast = OPENMM_METAL_FAST_MATH && !floating && !strict;
        for (int i = 0; i < inputValues.size(); i++) {
            // Fast mode makes no promise about nonfinite arithmetic. Only the
            // safe branch classifies these runtime inputs; no constant folding
            // of an embedded Inf/NaN can hide an incorrectly selected math mode.
            const int classification = fast ? 0 :
                    2*int(isinf(inputValues[i]))+4*int(isnan(inputValues[i]))+8*int(isfinite(inputValues[i]));
            ASSERT_EQUAL(int(fast)+classification, values[i]);
        }
    }
}

int main(int argc, char** argv) {
    try {
        testSelection();
        testShaderSource();
        if (argc == 1 || string(argv[1]) != "--selection-only") {
            testShaderTarget(false);
            testShaderTarget(true);
            for (bool common : {false, true})
                for (bool strict : {false, true})
                    testMathPolicy(common, strict);
        }
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
