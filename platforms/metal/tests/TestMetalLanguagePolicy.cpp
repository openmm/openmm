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
#include <iostream>

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
    ASSERT(common.find("#define SYNC_WARPS simdgroup_barrier(mem_flags::mem_threadgroup);") != string::npos);
    ASSERT(common.find("#define SYNC_THREADS threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);") != string::npos);

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

int main(int argc, char** argv) {
    try {
        testSelection();
        testShaderSource();
        if (argc == 1 || string(argv[1]) != "--selection-only") {
            testShaderTarget(false);
            testShaderTarget(true);
        }
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
