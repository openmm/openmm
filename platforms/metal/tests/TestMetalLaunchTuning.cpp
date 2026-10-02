/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform tests: copyright (c) 2026 Chun-Chi Hung.                      *
 * Author: Chun-Chi Hung.                                                       *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalNonbondedUtilities.h"
#include "openmm/internal/AssertionUtilities.h"
#include <AvailabilityMacros.h>
#include <algorithm>
#include <iostream>

using namespace OpenMM;
using namespace std;

/** Verify shared source/allocation geometry and the actual dispatched grid. */
void testLaunchGeometry(bool floating) {
    System system;
    system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, floating);
    auto& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    int expectedThreads = min(256, context.getMaxThreadBlockSize());
    expectedThreads -= expectedThreads%32;
    int expectedGroups = context.getNumComputeUnits() == 0 ? context.getNumThreadBlocks() : 6*context.getNumComputeUnits();
#if OPENMM_METAL_TUNE_FORCE_THREADGROUP_SIZE
    expectedThreads = OPENMM_METAL_FORCE_THREADGROUP_SIZE;
#endif
#if OPENMM_METAL_TUNE_FORCE_GROUPS_PER_COMPUTE_UNIT
    expectedGroups = OPENMM_METAL_FORCE_GROUPS_PER_COMPUTE_UNIT*context.getNumComputeUnits();
#endif
    ASSERT_EQUAL(expectedThreads, nb.getForceThreadBlockSize());
    ASSERT_EQUAL(expectedGroups, nb.getNumForceThreadBlocks());
    ASSERT(expectedGroups <= context.getNumThreadBlocks());
    ASSERT_EQUAL(expectedThreads*expectedGroups, nb.getNumEnergyBuffers());
    context.initialize();
    ASSERT(context.getEnergyBuffer().getSize() >= nb.getNumEnergyBuffers());

    // Declared production geometry opts in; the entry-point name alone does not.
    map<string, string> defines;
    defines["FORCE_WORK_GROUP_SIZE"] = to_string(expectedThreads);
    const string body = R"((GLOBAL int* output) {
        output[GLOBAL_ID] = LOCAL_SIZE+10000*NUM_GROUPS;
    })";
    ComputeArray output;
    const int count = expectedThreads*expectedGroups;
    output.initialize<int>(context, count+1, "launchGeometryOutput");
    bool exactSizeEnabled = false;
#if OPENMM_METAL_TUNE_FORCE_REQUIRED_THREADS && __MAC_OS_X_VERSION_MAX_ALLOWED >= 260000
    // The context never selects a 4.x language on an older runtime. This
    // separate SDK gate mirrors the host property's compile-time availability.
    exactSizeEnabled = context.getMetalLanguageVersion() >= 400;
#endif
    for (const string& name : {"computeNonbonded", "computeBornSum", "computeGBSAForce1"}) {
        ComputeKernel kernel = context.compileProgram("KERNEL void "+name+body, defines)->createKernel(name);
        kernel->addArg(output);
        // Reuse one kernel across both lazily created accumulator pipelines.
        for (bool mode : {floating, !floating, floating}) {
            context.setUseFloatingPointAccumulators(mode);
#if OPENMM_METAL_TUNE_FORCE_PIPELINE_MAX_THREADS
            ASSERT(kernel->getMaxBlockSize() <= OPENMM_METAL_FORCE_PIPELINE_MAX_THREADS);
#endif
            ASSERT(kernel->getMaxBlockSize() >= expectedThreads);
            output.upload(vector<int>(count+1, -1));
            kernel->execute(count, expectedThreads);
            vector<int> values;
            output.download(values);
            for (int i = 0; i < count; i++)
                ASSERT_EQUAL(expectedThreads+10000*expectedGroups, values[i]);
            ASSERT_EQUAL(-1, values[count]);

            vector<int> invalidSizes{kernel->getMaxBlockSize()+32};
            if (exactSizeEnabled) {
                invalidSizes.push_back(expectedThreads-32);
                invalidSizes.push_back(expectedThreads+32);
            }
            for (int invalidSize : invalidSizes) {
                bool rejected = false;
                try {
                    kernel->execute(count, invalidSize);
                }
                catch (const OpenMMException&) {
                    rejected = true;
                }
                ASSERT(rejected);
            }
            // Neither a newer target with the switch OFF nor an older target
            // with the switch ON may impose the optional exact-size contract.
            if (!exactSizeEnabled) {
                const int smallerSize = expectedThreads/2;
                output.upload(vector<int>(count+1, -1));
                kernel->execute(expectedThreads, smallerSize);
                output.download(values);
                for (int i = 0; i < expectedThreads; i++)
                    ASSERT_EQUAL(smallerSize+20000, values[i]);
                ASSERT_EQUAL(-1, values[expectedThreads]);
            }
        }
    }

    // Identical unmarked programs must retain identical unrestricted limits,
    // including an unrelated function sharing the force geometry define.
    ComputeKernel unmarked = context.compileProgram("KERNEL void computeNonbonded"+body)->createKernel("computeNonbonded");
    ComputeKernel other = context.compileProgram("KERNEL void unrelated"+body, defines)->createKernel("unrelated");
    ASSERT_EQUAL(unmarked->getMaxBlockSize(), other->getMaxBlockSize());
    // Both controls must still accept a smaller shape, even when the matching
    // force program has an exact size. Verify the actual dispatch, not just max.
    for (ComputeKernel kernel : {unmarked, other}) {
        kernel->addArg(output);
        for (bool mode : {floating, !floating, floating}) {
            context.setUseFloatingPointAccumulators(mode);
            output.upload(vector<int>(count+1, -1));
            kernel->execute(32, 32);
            vector<int> values;
            output.download(values);
            for (int i = 0; i < 32; i++)
                ASSERT_EQUAL(10032, values[i]);
            ASSERT_EQUAL(-1, values[32]);
        }
    }
#if OPENMM_METAL_TUNE_FORCE_PIPELINE_MAX_THREADS
    // A mismatched source contract must be rejected before an unsafe dispatch.
    if (OPENMM_METAL_FORCE_PIPELINE_MAX_THREADS < 256) {
        defines["FORCE_WORK_GROUP_SIZE"] = "256";
        bool rejected = false;
        try {
            context.compileProgram("KERNEL void computeBornSum"+body, defines)->createKernel("computeBornSum");
        }
        catch (const OpenMMException&) {
            rejected = true;
        }
        ASSERT(rejected);
    }
#endif
}

/** Validate shader-owned launch metadata without submitting any GPU work. */
void testRequiredThreadMetadata() {
    System system;
    system.addParticle(1);
    MetalContext context(system);
    if (context.getMetalLanguageVersion() < 400)
        return;
    auto compile = [&](const string& annotation) {
        string source = "#include <metal_stdlib>\nusing namespace metal;\n"
            "[[required_threads_per_threadgroup(64, 1, 1)]]\n"
            "[[user_annotation(\""+annotation+"\")]]\n"
            "kernel void metadataProbe() {}\n";
        return context.compileProgram(source)->createKernel("metadataProbe");
    };
    ComputeKernel kernel = compile("openmm_required_threads=64");
    kernel->execute(0, 64);
    bool rejected = false;
    try {
        kernel->execute(0, 32);
    }
    catch (const OpenMMException&) {
        rejected = true;
    }
    ASSERT(rejected);
    // Unrelated developer metadata must not impose an exact launch size.
    compile("unrelated_annotation=64")->execute(0, 32);
    for (const string& value : {"", "0", "-1", "+64", " 64", "64x", "128", "999999999999999999999999"}) {
        rejected = false;
        try {
            compile("openmm_required_threads="+value);
        }
        catch (const OpenMMException&) {
            rejected = true;
        }
        ASSERT(rejected);
    }
}

int main() {
    try {
        testLaunchGeometry(false);
        testLaunchGeometry(true);
        testRequiredThreadMetadata();
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
