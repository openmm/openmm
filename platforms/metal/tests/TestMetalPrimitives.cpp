/* -------------------------------------------------------------------------- *
 * OpenMM — Metal Platform                                                    *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                  *
 * Authors: Chun-Chi Hung                                                      *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                    *
 * See https://openmm.org/development.                                         *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Lesser General Public License as published    *
 * by the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version.                                         *
 * This program is distributed in the hope that it will be useful,              *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
 * GNU Lesser General Public License for more details.                         *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.         *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "openmm/System.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <cstdint>
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

/** @brief Check lane mapping, integer/float reductions, scans, and active masks. */
void testCollectives(MetalContext& context) {
    const string source = R"(
KERNEL void collectives(GLOBAL uint* output) {
    uint lane = LOCAL_ID&31, base = GLOBAL_ID*19;
    output[base] = simdBallot(lane%3 == 0);
    output[base+1] = simdAny(lane == 17);
    output[base+2] = simdAll(lane < 32);
    output[base+3] = simdBroadcast(lane, 5);
    output[base+4] = simdShuffle(lane, 31-lane);
    output[base+5] = simdShuffleDown(lane, 1);
    output[base+6] = simdShuffleXor(lane, 7);
    output[base+7] = simdRotate(lane, lane, 7);
    output[base+8] = simdReduceAdd(lane);
    output[base+9] = simdReduceMin(lane);
    output[base+10] = simdReduceMax(lane);
    output[base+11] = simdPrefixInclusiveAdd(lane);
    output[base+12] = simdPrefixExclusiveAdd(lane);
    output[base+13] = popcount(lane);
    output[base+14] = clz(1u<<lane);
    output[base+15] = ctz(1u<<lane);
    output[base+16] = uint(simdReduceAdd(float(lane)));
    // Ballots report active lanes, not all lanes in the physical group.
    if (lane < 16) {
        output[base+17] = simdBallot(true);
        output[base+18] = simdReduceAdd(lane);
    }
    else {
        output[base+17] = simdBallot(true);
        output[base+18] = simdReduceAdd(lane);
    }
}
)";
    ComputeArray output;
    output.initialize<unsigned int>(context, 19*64, "collectiveOutput");
    ComputeKernel kernel = context.compileProgram(source)->createKernel("collectives");
    kernel->addArg(output);
    kernel->execute(64, 64);
    vector<unsigned int> actual;
    output.download(actual);
    for (unsigned int i = 0; i < 64; i++) {
        unsigned int lane = i%32, base = i*19, count = 0;
        for (unsigned int bits = lane; bits != 0; bits >>= 1)
            count += bits&1;
        const unsigned int expected[] = {0x49249249u, 1, 1, 5, 31-lane,
            lane+1, lane^7u, (lane+7)&31u, 496, 0, 31,
            lane*(lane+1)/2, lane*(lane-1)/2, count, 31-lane, lane, 496,
            lane < 16 ? 0xffffu : 0xffff0000u, lane < 16 ? 120u : 376u};
        for (int j = 0; j < 19; j++)
            if (j != 5 || lane < 31) // Shuffle past the active group is unspecified.
                ASSERT_EQUAL(expected[j], actual[base+j]);
    }
}

/** @brief Device/threadgroup int32 and device float operations keep their distinct APIs. */
void testAtomics(MetalContext& context) {
    const string source = R"(
KERNEL void atomics(GLOBAL int* integers, GLOBAL float* floats, GLOBAL mm_ulong* wide) {
    LOCAL int subtotal;
    if (LOCAL_ID == 0) subtotal = 0;
    SYNC_THREADS;
    atomicAddInt32(&subtotal, 1);
    for (int i = GLOBAL_ID; i < 1025; i += GLOBAL_SIZE) {
        atomicAddInt32(integers, 1);
        atomicAddFloat32(floats, 1.0f);
#if OPENMM_METAL_HAS_UINT64_MIN_MAX
        atomicMinUInt64(wide, (mm_ulong(1)<<40)+mm_ulong(i));
        atomicMaxUInt64(wide+1, (mm_ulong(1)<<40)+mm_ulong(i));
#endif
    }
    SYNC_THREADS;
    if (LOCAL_ID == 0) integers[1+GROUP_ID] = subtotal;
    if (GLOBAL_ID == 0) integers[3] = OPENMM_METAL_HAS_UINT64_MIN_MAX;
}
)";
    ComputeArray integers, floats, wide;
    integers.initialize<int>(context, 4, "atomicIntegers");
    floats.initialize<float>(context, 1, "atomicFloats");
    wide.initialize<uint64_t>(context, 2, "atomicWide");
    for (int native = 0; native < 2; native++) {
        context.clearBuffer(integers);
        context.clearBuffer(floats);
        wide.upload(vector<uint64_t>{~uint64_t(0), 0});
        map<string, string> defines;
        defines["OPENMM_METAL_NATIVE_FLOAT_ATOMICS"] = to_string(native);
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("atomics");
        kernel->addArg(integers);
        kernel->addArg(floats);
        kernel->addArg(wide);
        kernel->execute(128, 64);
        vector<int> actualInts;
        vector<float> actualFloats;
        vector<uint64_t> actualWide;
        integers.download(actualInts);
        floats.download(actualFloats);
        wide.download(actualWide);
        ASSERT_EQUAL(1025, actualInts[0]);
        ASSERT_EQUAL(64, actualInts[1]);
        ASSERT_EQUAL(64, actualInts[2]);
        ASSERT_EQUAL(1025.0f, actualFloats[0]);
        if (actualInts[3]) {
            ASSERT_EQUAL(uint64_t(1)<<40, actualWide[0]);
            ASSERT_EQUAL((uint64_t(1)<<40)+1024, actualWide[1]);
        }
        else {
            ASSERT_EQUAL(~uint64_t(0), actualWide[0]);
            ASSERT_EQUAL(uint64_t(0), actualWide[1]);
        }
    }
    bool rejected = false;
    try {
        context.compileProgram("KERNEL void unsupported(GLOBAL mm_ulong* v) { atomicAddUInt64(v, mm_ulong(1)); }");
    }
    catch (const OpenMMException& error) {
        const string message = error.what();
        rejected = message.find("atomicAddUInt64") != string::npos && message.find("deleted") != string::npos;
    }
    ASSERT(rejected);
}

/** @brief Grouped Q32.32 writes preserve exact unsigned sums across carries, signs, and inactive lanes. */
void testGroupedFixedPoint(MetalContext& context) {
    const string source = R"(
KERNEL void grouped(GLOBAL const mm_ulong* input, GLOBAL mm_ulong* output,
        unsigned int count, unsigned int bins, unsigned int mask) {
    for (unsigned int i = GLOBAL_ID; i < count; i += GLOBAL_SIZE) {
        if ((mask>>(LOCAL_ID&31))&1u)
            METAL_ACCUMULATE_SPARSE_FORCE(output, (i*7+i/32)%bins, input[i]);
    }
}

)";
    const uint64_t patterns[] = {0, 1, ~uint64_t(0), uint64_t(1)<<32,
        (uint64_t(1)<<32)-1, uint64_t(1)<<63, 0x0000ffffffffffffULL,
        0xfedcba9876543210ULL, 0x123456789abcdef0ULL};
    const int count = 4099;
    vector<uint64_t> input(count);
    for (int i = 0; i < count; i++) input[i] = patterns[i%9];
    ComputeArray values;
    values.initialize<uint64_t>(context, count, "groupedFixedPointValues");
    values.upload(input);
    for (int enabled = 0; enabled < 2; enabled++) {
        map<string, string> defines;
        defines["OPENMM_METAL_FAST_SPARSE_FORCE_AGGREGATION"] = to_string(enabled);
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("grouped");
        for (int i = 0; i < 5; i++) kernel->addArg();
        for (int bins : {1, 3, 17, 32, 67}) {
            ComputeArray output;
            output.initialize<uint64_t>(context, bins, "groupedFixedPointOutput");
            for (unsigned int mask : {0u, 1u, 0x80000000u, 0xaaaaaaaau, 0x55555555u, 0xffffffffu}) {
                context.clearBuffer(output);
                kernel->setArg(0, values);
                kernel->setArg(1, output);
                kernel->setArg(2, (unsigned int) count);
                kernel->setArg(3, (unsigned int) bins);
                kernel->setArg(4, mask);
                kernel->execute(count, 64);
                vector<uint64_t> actual, expected(bins, 0);
                output.download(actual);
                for (int i = 0; i < count; i++)
                    if ((mask>>(i&31))&1u) expected[(i*7+i/32)%bins] += input[i];
                ASSERT_EQUAL_CONTAINERS(expected, actual);
            }
        }
    }
}

/** @brief CUDA's 64-bit parameter shuffles transport two words without converting to float. */
void testWideShuffle(MetalContext& context) {
    ComputeArray input, output;
    input.initialize<uint64_t>(context, 64, "wideShuffleInput");
    output.initialize<uint64_t>(context, 128, "wideShuffleOutput");
    vector<uint64_t> values(64);
    for (int i = 0; i < 64; i++) values[i] = 0xfedcba9876543210ULL+uint64_t(i)*0x123456789ULL;
    values[0] = 0;
    values[1] = ~uint64_t(0);
    values[2] = uint64_t(1)<<63;
    values[3] = (uint64_t(1)<<32)-1;
    values[4] = uint64_t(1)<<32;
    input.upload(values);
    ComputeKernel kernel = context.compileProgram(R"(
KERNEL void wideShuffle(GLOBAL const mm_ulong* input, GLOBAL mm_ulong* output) {
    uint lane = LOCAL_ID&31;
    mm_ulong bits = input[GLOBAL_ID];
    output[GLOBAL_ID] = simdShuffle(bits, (lane+13)&31u);
    output[GLOBAL_ID+64] = as_type<mm_ulong>(simdShuffle(as_type<mm_long>(bits), (lane+17)&31u));
}
)")->createKernel("wideShuffle");
    kernel->addArg(input);
    kernel->addArg(output);
    kernel->execute(64, 64);
    vector<uint64_t> actual;
    output.download(actual);
    for (int i = 0; i < 64; i++) {
        ASSERT_EQUAL(values[(i/32)*32+(i+13)%32], actual[i]);
        ASSERT_EQUAL(values[(i/32)*32+(i+17)%32], actual[i+64]);
    }
}

int main() {
    try {
        System system;
        system.addParticle(1);
        MetalContext context(system);
        testCollectives(context);
        testAtomics(context);
        testGroupedFixedPoint(context);
        testWideShuffle(context);
    }
    catch (const exception& error) {
        if (string(error.what()).find("No Metal device") != string::npos)
            return 77;
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Metal SIMD and atomic primitives passed" << endl;
    return 0;
}
