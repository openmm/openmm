/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
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
#include "MetalQueue.h"
#include "openmm/System.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeVectorTypes.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace OpenMM;
using namespace std;

// These tests exercise the Common runtime interfaces with native MSL smoke
// kernels.  They do not establish Common kernel-source or simulation support.

void expectException(const string& description, const function<void()>& operation) {
    try {
        operation();
    }
    catch (const OpenMMException&) {
        return;
    }
    throw OpenMMException("Expected an exception: "+description);
}

string loadSource() {
    const string path = string(METAL_TEST_SOURCE_DIR)+"/runtime.metal";
    ifstream input(path.c_str());
    if (!input)
        throw OpenMMException("Cannot open Metal test source: "+path);
    return string(istreambuf_iterator<char>(input), istreambuf_iterator<char>());
}

vector<int> inputValues(int count) {
    vector<int> values(count);
    for (int i = 0; i < count; i++)
        values[i] = i-17;
    return values;
}

void checkTransform(const vector<int>& input, const vector<int>& output, int offset) {
    ASSERT_EQUAL(input.size(), output.size());
    for (size_t i = 0; i < input.size(); i++)
        ASSERT_EQUAL(3*input[i]+offset, output[i]);
}

ComputeKernel createTransform(ComputeProgram program, ArrayInterface& input,
        ArrayInterface& output, int count, int offset) {
    ComputeKernel kernel = program->createKernel("transform");
    kernel->addArg(input);
    kernel->addArg(output);
    kernel->addArg(count);
    kernel->addArg(offset);
    return kernel;
}

void testContext(ComputeContext& context) {
    const size_t paddedAtoms = context.getPaddedNumAtoms();
    ASSERT_EQUAL(65, context.getNumAtoms());
    ASSERT(context.getPaddedNumAtoms() >= context.getNumAtoms());
    ASSERT_EQUAL(1, context.getNumContexts());
    ASSERT(!context.getUseDoublePrecision());
    ASSERT(!context.getUseMixedPrecision());
    ASSERT(!context.getSupports64BitGlobalAtomics());
    ASSERT_EQUAL(paddedAtoms, context.getPosq().getSize());
    ASSERT_EQUAL(sizeof(mm_float4), context.getPosq().getElementSize());
    ASSERT_EQUAL(paddedAtoms, context.getVelm().getSize());
    ASSERT_EQUAL(sizeof(mm_float4), context.getVelm().getElementSize());
    ASSERT_EQUAL(3*paddedAtoms, context.getLongForceBuffer().getSize());
    ASSERT_EQUAL(sizeof(int64_t), context.getLongForceBuffer().getElementSize());
    ASSERT_EQUAL(sizeof(float), context.getEnergyBuffer().getElementSize());
    ASSERT_EQUAL(paddedAtoms, context.getAtomIndexArray().getSize());
}

void testArrays(ComputeContext& context) {
    ComputeArray source, destination;
    ASSERT(!source.isInitialized());
    source.initialize<int>(context, 9, "source");
    destination.initialize<int>(context, 9, "destination");
    ArrayInterface& array = source;
    ASSERT(array.isInitialized());
    ASSERT_EQUAL(string("source"), array.getName());
    ASSERT(&array.getContext() == &context);
    vector<int> expected = inputValues(9), result;
    array.upload(expected);
    const int replacement[] = {101, 102, 103};
    array.uploadSubArray(replacement, 2, 3);
    copy(replacement, replacement+3, expected.begin()+2);
    array.download(result);
    ASSERT_EQUAL_CONTAINERS(expected, result);
    array.copyTo(destination);
    destination.download(result);
    ASSERT_EQUAL_CONTAINERS(expected, result);
    destination.resize(17);
    expected = inputValues(17);
    destination.upload(expected);
    destination.download(result);
    ASSERT_EQUAL_CONTAINERS(expected, result);
    expectException("copy size mismatch", [&] { array.copyTo(destination); });
    expectException("upload size mismatch", [&] { array.upload(expected); });
    expectException("negative subarray offset", [&] { array.uploadSubArray(replacement, -1, 1); });
    expectException("subarray past end", [&] { array.uploadSubArray(replacement, 8, 2); });
    expectException("negative subarray length", [&] { array.uploadSubArray(replacement, 0, -1); });
    expectException("non-pinned asynchronous upload", [&] { array.upload(expected.data(), false); });
    expectException("non-pinned asynchronous download", [&] { array.download(result.data(), false); });
    const size_t pinnedBytes = max(context.getLongForceBuffer().getSize()*context.getLongForceBuffer().getElementSize(),
            context.getEnergyBuffer().getSize()*context.getEnergyBuffer().getElementSize());
    unsigned char* pinned = static_cast<unsigned char*>(context.getPinnedBuffer());
    expectException("pinned upload exceeds capacity", [&] { array.upload(pinned+pinnedBytes-sizeof(int), false); });
    expectException("pinned download exceeds capacity", [&] { array.download(pinned+pinnedBytes-sizeof(int), false); });
    expectException("pinned offset past end", [&] { array.download(pinned+pinnedBytes, false); });

    // Use raw-pointer overloads: Common's vector helpers assume nonempty vectors.
    ComputeArray empty, emptyCopy;
    empty.initialize<int>(context, 0, "empty");
    emptyCopy.initialize<int>(context, 0, "emptyCopy");
    ASSERT_EQUAL(size_t(0), empty.getSize());
    empty.upload(static_cast<const void*>(nullptr));
    empty.download(static_cast<void*>(nullptr));
    empty.copyTo(emptyCopy);
    context.clearBuffer(empty);
    empty.resize(1);
    empty.upload(vector<int>(1, 9));
    empty.download(result);
    ASSERT_EQUAL(9, result[0]);
    empty.resize(0);
    context.clearBuffer(empty);
}

void testLaunches(ComputeContext& context, ComputeProgram program) {
    const int sizes[] = {1, 63, 64, 65, 129, context.getNumThreadBlocks()*64+129};
    for (int count : sizes) {
        ComputeArray input, output, groups;
        input.initialize<int>(context, count, "launchInput");
        output.initialize<int>(context, count, "launchOutput");
        vector<int> values = inputValues(count), result;
        input.upload(values);
        ComputeKernel kernel = createTransform(program, input, output, count, 5);
        // Cover both explicit and default launch sizes.
        kernel->execute(count, 64);
        output.download(result);
        checkTransform(values, result, 5);
        kernel->setArg(3, -9);
        kernel->execute(count);
        output.download(result);
        checkTransform(values, result, -9);

        const int numGroups = min((count+63)/64, context.getNumThreadBlocks());
        groups.initialize<unsigned int>(context, numGroups, "groupWidths");
        context.clearBuffer(groups);
        ComputeKernel widths = program->createKernel("recordGroupWidth");
        widths->addArg(groups);
        widths->execute(count, 64);
        vector<unsigned int> actualWidths;
        groups.download(actualWidths);
        for (unsigned int width : actualWidths)
            ASSERT_EQUAL(64, width);
    }
}

void testArguments(ComputeContext& context, ComputeProgram program) {
    ComputeArray input, output, replacement;
    input.initialize<int>(context, 1, "argumentInput");
    output.initialize<int>(context, 1, "argumentOutput");
    replacement.initialize<int>(context, 1, "replacementOutput");
    vector<int> values = inputValues(1), result;
    input.upload(values);
    ComputeKernel kernel = program->createKernel("transform");
    kernel->addArg();
    kernel->addArg();
    kernel->addArg();
    kernel->addArg();
    expectException("unbound kernel arguments", [&] { kernel->execute(1); });
    kernel->setArg(0, input);
    kernel->setArg(1, output);
    kernel->setArg(2, 1);
    kernel->setArg(3, 7);
    kernel->execute(1);
    output.download(result);
    checkTransform(values, result, 7);
    kernel->setArg(1, replacement);
    kernel->setArg(3, -2);
    kernel->execute(1);
    replacement.download(result);
    checkTransform(values, result, -2);
    output.download(result);
    checkTransform(values, result, 7);

    // Bound ArrayInterface objects must resolve their current storage after resize.
    input.resize(129);
    replacement.resize(129);
    values = inputValues(129);
    input.upload(values);
    kernel->setArg(2, 129);
    kernel->execute(129);
    replacement.download(result);
    checkTransform(values, result, -2);
    expectException("invalid argument index", [&] { kernel->setArg(4, 1); });
    struct OversizedArgument { double values[5]; } oversized = {};
    expectException("oversized primitive argument", [&] { kernel->setArg(3, oversized); });
    expectException("zero block size", [&] { kernel->execute(129, 0); });
    expectException("excessive block size", [&] { kernel->execute(129, kernel->getMaxBlockSize()+1); });
}

void testClearing(MetalContext& metal) {
    ComputeContext& context = metal;
    ArrayInterface& force = context.getLongForceBuffer();
    vector<int64_t> forceValues(force.getSize(), (int64_t(1)<<40)+7), forceResult;
    force.upload(forceValues);
    force.download(forceResult);
    ASSERT_EQUAL_CONTAINERS(forceValues, forceResult);
    // Common uses the pinned workspace for the full force buffer.  flushQueue()
    // must make the asynchronous device-to-host transfer observable on the CPU.
    void* pinned = context.getPinnedBuffer();
    ASSERT(pinned != nullptr);
    force.download(pinned, false);
    context.flushQueue();
    ASSERT_EQUAL(0, memcmp(forceValues.data(), pinned, forceValues.size()*sizeof(int64_t)));
    context.clearBuffer(force);
    force.download(forceResult);
    for (int64_t value : forceResult)
        ASSERT_EQUAL(0, value);

    // The registered array is context-owned and outlives subsequent clear calls.
    ArrayInterface& energy = context.getEnergyBuffer();
    vector<float> energyValues(energy.getSize(), 7.0f), energyResult;
    context.addAutoclearBuffer(energy);
    for (int repeat = 0; repeat < 2; repeat++) {
        energy.upload(energyValues);
        metal.clearAutoclearBuffers();
        energy.download(energyResult);
        for (float value : energyResult)
            ASSERT_EQUAL(0.0f, value);
    }
}

void testQueues(MetalContext& metal, ComputeProgram program) {
    ComputeContext& context = metal;
    ComputeQueue original = context.getCurrentQueue(), secondary = context.createQueue();
    ComputeArray input, output, copied;
    const int count = 129;
    input.initialize<int>(context, count, "queueInput");
    output.initialize<int>(context, count, "queueOutput");
    copied.initialize<int>(context, count, "queueCopy");
    vector<int> values = inputValues(count), result(count, 0);
    int* pinned = static_cast<int*>(context.getPinnedBuffer());
    // Exercise an interior pinned pointer as well as the workspace base.
    int* offsetPinned = pinned+4;
    const size_t bytes = count*sizeof(int);
    memcpy(offsetPinned, values.data(), bytes);
    input.upload(offsetPinned, false);
    ComputeEvent uploaded = context.createEvent();
    uploaded->enqueue();
    uploaded->queueWait(secondary);
    context.setCurrentQueue(secondary);
    ComputeKernel kernel = createTransform(program, input, output, count, 23);
    kernel->execute(count);
    output.copyTo(copied);
    // The uploaded event orders this write after the input transfer consumed
    // the same pinned memory.  No CPU access occurs while either is in flight.
    copied.download(offsetPinned, false);
    ComputeEvent computed = context.createEvent();
    computed->enqueue();
    computed->queueWait(original);
    context.restoreDefaultQueue();
    ASSERT(context.getCurrentQueue() == original);

    // A GPU event dependency must make the source queue's readback transitively
    // observable after finishing the target queue, without a source host wait.
    context.flushQueue();
    memcpy(result.data(), offsetPinned, bytes);
    checkTransform(values, result, 23);
    kernel->setArg(3, -4);
    kernel->execute(count);
    output.download(pinned, false);
    ComputeEvent downloaded = context.createEvent();
    downloaded->enqueue();
    downloaded->wait();
    memcpy(result.data(), pinned, bytes);
    checkTransform(values, result, -4);
    kernel->setArg(3, 17);
    kernel->execute(count);
    output.download(offsetPinned, false);
    downloaded->enqueue();
    downloaded->wait();
    memcpy(result.data(), offsetPinned, bytes);
    checkTransform(values, result, 17);
    metal.getCurrentMetalQueue().finish();
}

void testErrors(ComputeContext& context, ComputeProgram program) {
    expectException("simulation context is unavailable", [&] { context.getContextImpl(); });
    expectException("integration utilities are unavailable", [&] { context.getIntegrationUtilities(); });
    expectException("nonbonded utilities are unavailable", [&] { context.getNonbondedUtilities(); });
    expectException("invalid shader source", [&] { context.compileProgram("This is not valid MSL."); });
    expectException("missing kernel", [&] { program->createKernel("missingKernel"); });
    System emptySystem;
    MetalContext other(emptySystem);
    ASSERT_EQUAL(0, other.getNumAtoms());
    ComputeArray local, foreign;
    local.initialize<int>(context, 1, "local");
    foreign.initialize<int>(other, 1, "foreign");
    expectException("cross-context array copy", [&] { local.copyTo(foreign); });
    expectException("cross-context array clear", [&] { context.clearBuffer(foreign); });
    ComputeKernel kernel = program->createKernel("transform");
    expectException("cross-context kernel argument", [&] { kernel->addArg(foreign); });
}

int main() {
    System system;
    for (int i = 0; i < 65; i++)
        system.addParticle(1.0);
    unique_ptr<MetalContext> metal;
    try {
        metal.reset(new MetalContext(system));
    }
    catch (const exception& error) {
        if (string(error.what()).find("No Metal device") != string::npos) {
            cout << "SKIPPED: No Metal device; GPU tests were not executed." << endl;
            return 77;
        }
        cerr << "Exception creating Metal context: " << error.what() << endl;
        return 1;
    }
    try {
        ComputeContext& context = *metal;
        map<string, string> defines;
        defines["VALUE_SCALE"] = "3";
        ComputeProgram program = context.compileProgram(loadSource(), defines);
        testContext(context);
        testArrays(context);
        testLaunches(context, program);
        testArguments(context, program);
        testClearing(*metal);
        testQueues(*metal, program);
        testErrors(context, program);
        metal->getCurrentMetalQueue().finish();
    }
    catch (const exception& error) {
        cerr << "Exception: " << error.what() << endl;
        return 1;
    }
    cout << "Done: Metal Platform Common interface tests executed successfully." << endl;
    return 0;
}
