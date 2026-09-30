/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code: Portions copyright (c) 2026 Chun-Chi Hung.             *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.*
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalKernel.h"
#include "MetalReductionOptimizations.h"
#include "CommonKernelSources.h"
#include "MetalOpenCLKernelSources.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <numeric>

using namespace OpenMM;
using namespace std;

/** @brief Extract a test kernel without compiling unrelated force templates. */
string extractFunction(const string& source, const string& signature) {
    size_t begin = source.find(signature);
    ASSERT(begin != string::npos);
    size_t body = source.find('{', begin);
    int depth = 1;
    size_t end = body+1;
    for (; end < source.size() && depth != 0; end++) {
        if (source[end] == '{') depth++;
        if (source[end] == '}') depth--;
    }
    ASSERT_EQUAL(0, depth);
    return source.substr(begin, end-begin);
}

/** @brief Verify independent source selection and immunity of unrelated same-name helpers. */
void testSelection() {
    using Settings = MetalReductionOptimizations::Settings;
    Settings none;
    for (const string* source : {&CommonKernelSources::rg, &CommonKernelSources::rmsd,
            &CommonKernelSources::orientationRestraintForce, &CommonKernelSources::customCentroidBond,
            &CommonKernelSources::lcpo, &CommonKernelSources::customManyParticle, &MetalOpenCLKernelSources::sort})
        ASSERT_EQUAL(*source, MetalReductionOptimizations::apply(*source, none));
    for (int option = 0; option < 4; option++) {
        Settings settings;
        settings.centroid = option == 0;
        settings.rg = option == 1;
        settings.rmsd = option == 2;
        settings.orientation = option == 3;
        int index = 0;
        for (const string* source : {&CommonKernelSources::customCentroidBond, &CommonKernelSources::rg,
                &CommonKernelSources::rmsd, &CommonKernelSources::orientationRestraintForce}) {
            string optimized = MetalReductionOptimizations::apply(*source, settings);
            ASSERT_EQUAL(option == index++, optimized.find("simd_shuffle_down") != string::npos);
            ASSERT_EQUAL(optimized, MetalReductionOptimizations::apply(optimized, settings));
        }
    }
    for (int option = 0; option < 2; option++) {
        Settings settings;
        settings.lcpoScan = option == 0;
        settings.manyParticleScan = option == 1;
        int index = 0;
        for (const string* source : {&CommonKernelSources::lcpo, &CommonKernelSources::customManyParticle,
                &MetalOpenCLKernelSources::sort}) {
            string optimized = MetalReductionOptimizations::apply(*source, settings);
            ASSERT_EQUAL(option == index++, optimized.find("metalBlockInclusiveScan") != string::npos);
            ASSERT_EQUAL(optimized, MetalReductionOptimizations::apply(optimized, settings));
        }
    }
    Settings all;
    all.centroid = true;
    all.rg = true;
    all.rmsd = true;
    all.orientation = true;
    all.lcpoScan = true;
    all.manyParticleScan = true;
    const string unrelated = "DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) { return value; }";
    ASSERT_EQUAL(unrelated, MetalReductionOptimizations::apply(unrelated, all));
    // A recognizable function name is insufficient: an unreviewed body must
    // retain its original implementation, even when the option is enabled.
    const string* originals[] = {&CommonKernelSources::customCentroidBond, &CommonKernelSources::rg,
        &CommonKernelSources::rmsd, &CommonKernelSources::orientationRestraintForce,
        &CommonKernelSources::lcpo, &CommonKernelSources::customManyParticle,
        &MetalOpenCLKernelSources::sort, &MetalOpenCLKernelSources::sort};
    const char* signatures[] = {"KERNEL void computeGroupCenters(", "KERNEL void computeCenterPosition(",
        "KERNEL void computeRMSDPart1(", "KERNEL void computeCorrelationMatrix(",
        "KERNEL void computeNeighborStartIndices(", "KERNEL void computeNeighborStartIndices(",
        "__kernel void computeBucketPositions(", "__kernel void sortShortList("};
    for (int option = 0; option < 6; option++) {
        Settings settings;
        settings.centroid = option == 0;
        settings.rg = option == 1;
        settings.rmsd = option == 2;
        settings.orientation = option == 3;
        settings.lcpoScan = option == 4;
        settings.manyParticleScan = option == 5;
        string modified = *originals[option];
        size_t body = modified.find('{', modified.find(signatures[option]));
        ASSERT(body != string::npos);
        modified.insert(body+1, "\n/* Unreviewed template body. */\n");
        ASSERT_EQUAL(modified, MetalReductionOptimizations::apply(modified, settings));
    }
}

/** @brief Compare every lane's reduction result, including partially occupied SIMD groups and scratch reuse. */
void testReductions(MetalContext& context) {
    MetalReductionOptimizations::Settings settings;
    settings.rg = true;
    const string optimized = MetalReductionOptimizations::apply(CommonKernelSources::rg, settings);
    const string probe = R"(
KERNEL void reduceProbe(GLOBAL const real* input, GLOBAL real* output, int count) {
    LOCAL volatile real temp[1024];
    real first = 0, second = 0;
    for (int i = LOCAL_ID; i < count; i += LOCAL_SIZE) {
        first += input[i];
        second += input[i]*input[i];
    }
    first = reduceValue(first, temp);
    second = reduceValue(second, temp);
    output[LOCAL_ID] = first;
    output[LOCAL_ID+LOCAL_SIZE] = second;
}
)";
    const int count = 4001;
    vector<float> values(count);
    double first = 0, second = 0;
    for (int i = 0; i < count; i++) {
        values[i] = i%17-8;
        first += values[i];
        second += values[i]*values[i];
    }
    ComputeArray input, output;
    input.initialize<float>(context, count, "reductionInput");
    output.initialize<float>(context, 2048, "reductionOutput");
    input.upload(values);
    for (const string* source : {&CommonKernelSources::rg, &optimized}) {
        ComputeKernel kernel = context.compileProgram(extractFunction(*source, "DEVICE real reduceValue(")+probe)->createKernel("reduceProbe");
        kernel->addArg(input);
        kernel->addArg(output);
        kernel->addArg(count);
        for (int width : {1, 7, 31, 32, 33, 63, 64, 65, 127, 128, 255, 256, 1024}) {
            if (width > kernel->getMaxBlockSize())
                continue;
            kernel->execute(width, width);
            vector<float> result;
            output.download(result);
            for (int lane = 0; lane < width; lane++) {
                ASSERT_EQUAL(first, result[lane]);
                ASSERT_EQUAL(second, result[lane+width]);
            }
        }
    }
}

/** @brief Check weighted centroid reductions across short groups and repeated group iterations. */
void testCentroids(MetalContext& context) {
    vector<int> offsets{0};
    vector<int> particles;
    vector<float> weights;
    vector<mm_float4> positions;
    vector<Vec3> expected;
    for (int count : {1, 7, 31, 32, 33, 64, 257}) {
        Vec3 sum;
        for (int i = 0; i < count; i++) {
            mm_float4 value(i%17-8, i%9-4, i%13-6, 0);
            float weight = 1.0f/count;
            particles.push_back(positions.size());
            positions.push_back(value);
            weights.push_back(weight);
            sum += weight*Vec3(value.x, value.y, value.z);
        }
        expected.push_back(sum);
        offsets.push_back(particles.size());
    }
    ComputeArray posq, atomIndices, groupWeights, groupOffsets, centers;
    posq.initialize<mm_float4>(context, positions.size(), "centroidPositions");
    atomIndices.initialize<int>(context, particles.size(), "centroidParticles");
    groupWeights.initialize<float>(context, weights.size(), "centroidWeights");
    groupOffsets.initialize<int>(context, offsets.size(), "centroidOffsets");
    centers.initialize<mm_float4>(context, expected.size(), "centroidOutput");
    posq.upload(positions);
    atomIndices.upload(particles);
    groupWeights.upload(weights);
    groupOffsets.upload(offsets);
    for (bool enabled : {false, true}) {
        MetalReductionOptimizations::Settings settings;
        settings.centroid = enabled;
        string transformed = MetalReductionOptimizations::apply(CommonKernelSources::customCentroidBond, settings);
        string source = extractFunction(transformed, "KERNEL void computeGroupCenters(");
        if (enabled)
            source = extractFunction(transformed, "DEVICE real metalReduceGroupValue(")+source;
        source = context.replaceStrings(source, {{"computeGroupCenters", "centroidProbe"}});
        ComputeKernel kernel = context.compileProgram(source)->createKernel("centroidProbe");
        kernel->addArg((int) expected.size());
        kernel->addArg(posq);
        kernel->addArg(atomIndices);
        kernel->addArg(groupWeights);
        kernel->addArg(groupOffsets);
        kernel->addArg(centers);
        for (int groups : {1, 2, 7}) {
            kernel->execute(64*groups, 64);
            vector<mm_float4> result;
            centers.download(result);
            for (int i = 0; i < (int) expected.size(); i++) {
                ASSERT_EQUAL_VEC(expected[i], Vec3(result[i].x, result[i].y, result[i].z), 1e-5);
                ASSERT_EQUAL(0, result[i].w);
            }
        }
    }
}

/** @brief Exercise inclusive scans, chunk carries, overflow exits, and neighbor-counter clearing. */
void testNeighborScans(MetalContext& context) {
    for (bool lcpo : {false, true}) {
        const string& original = lcpo ? CommonKernelSources::lcpo : CommonKernelSources::customManyParticle;
        for (bool enabled : {false, true}) {
            MetalReductionOptimizations::Settings settings;
            settings.lcpoScan = lcpo && enabled;
            settings.manyParticleScan = !lcpo && enabled;
            string transformed = MetalReductionOptimizations::apply(original, settings);
            string kernelSource = extractFunction(transformed, "KERNEL void computeNeighborStartIndices(");
            if (enabled)
                kernelSource = extractFunction(transformed, "DEVICE unsigned int metalBlockInclusiveScan(")+kernelSource;
            // Preserve explicit OFF in an all-options-ON build.
            kernelSource = context.replaceStrings(kernelSource, {{"computeNeighborStartIndices", "neighborScanProbe"}});
            for (int count : {0, 1, 31, 33, 257, 4099}) {
                vector<int> counts(max(1, count));
                for (int i = 0; i < count; i++)
                    counts[i] = i%7;
                int total = accumulate(counts.begin(), counts.end(), 0);
                ComputeArray countsBuffer, startsBuffer, pairCount;
                countsBuffer.initialize<int>(context, counts.size(), "scanCounts");
                startsBuffer.initialize<int>(context, count+1, "scanStarts");
                pairCount.initialize<int>(context, 1, "scanPairCount");
                map<string, string> defines{{"NUM_ACTIVE", to_string(count)}, {"NUM_ATOMS", to_string(count)}, {"THREAD_BLOCK_SIZE", "256"}};
                ComputeKernel kernel = context.compileProgram(kernelSource, defines)->createKernel("neighborScanProbe");
                if (lcpo) {
                    kernel->addArg(pairCount);
                    kernel->addArg(countsBuffer);
                    kernel->addArg(startsBuffer);
                }
                else {
                    kernel->addArg(countsBuffer);
                    kernel->addArg(startsBuffer);
                    kernel->addArg(pairCount);
                }
                kernel->addArg(total);
                for (int width : {7, 32, 33, 64, 256}) {
                    for (bool overflow : {false, true}) {
                        countsBuffer.upload(counts);
                        startsBuffer.upload(vector<int>(count+1, -7));
                        pairCount.upload(vector<int>{overflow ? total+1 : total});
                        kernel->execute(width, width);
                        vector<int> starts, remaining;
                        startsBuffer.download(starts);
                        countsBuffer.download(remaining);
                        int expected = 0;
                        for (int i = 0; i <= count; i++) {
                            ASSERT_EQUAL(overflow ? (lcpo ? -7 : 0) : expected, starts[i]);
                            if (i < count) {
                                expected += counts[i];
                                ASSERT_EQUAL(overflow ? counts[i] : 0, remaining[i]);
                            }
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        testSelection();
        if (argc == 2 && string(argv[1]) == "--source-only") {
            cout << "Metal reduction/scan source-selection tests passed" << endl;
            return 0;
        }
        System system;
        system.addParticle(1.0);
        MetalContext context(system);
        testReductions(context);
        testCentroids(context);
        testNeighborScans(context);
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Metal reduction/scan optimization tests passed" << endl;
    return 0;
}
