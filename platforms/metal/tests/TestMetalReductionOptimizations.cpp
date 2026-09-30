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
    for (int option = 0; option < 1; option++) {
        Settings settings;
        settings.centroid = option == 0;
        int index = 0;
        for (const string* source : {&CommonKernelSources::customCentroidBond, &CommonKernelSources::rg,
                &CommonKernelSources::rmsd, &CommonKernelSources::orientationRestraintForce}) {
            string optimized = MetalReductionOptimizations::apply(*source, settings);
            ASSERT_EQUAL(option == index++, optimized.find("simd_shuffle_down") != string::npos);
            ASSERT_EQUAL(optimized, MetalReductionOptimizations::apply(optimized, settings));
        }
    }
    Settings all;
    all.centroid = true;
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
    for (int option = 0; option < 1; option++) {
        Settings settings;
        settings.centroid = option == 0;
        string modified = *originals[option];
        size_t body = modified.find('{', modified.find(signatures[option]));
        ASSERT(body != string::npos);
        modified.insert(body+1, "\n/* Unreviewed template body. */\n");
        ASSERT_EQUAL(modified, MetalReductionOptimizations::apply(modified, settings));
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
        testCentroids(context);
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Metal reduction/scan optimization tests passed" << endl;
    return 0;
}
