/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Adapted from Common reduction kernels and OpenCL sort kernels.             *
 * Original OpenMM code:                                                      *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.       *
 * Authors: Peter Eastman, Evan Pretti                                         *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
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

#include "MetalReductionOptimizations.h"
#include "CommonKernelSources.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace std;

namespace {

/** @brief Extract a known template function; these audited templates contain no conditional braces. */
string functionText(const string& source, const string& signature) {
    size_t begin = source.find(signature);
    if (begin == string::npos)
        throw OpenMMException("Missing Metal optimization template: "+signature);
    size_t body = source.find('{', begin);
    int depth = 1;
    size_t end = body+1;
    for (; end < source.size() && depth != 0; end++) {
        if (source[end] == '{') depth++;
        if (source[end] == '}') depth--;
    }
    if (depth != 0)
        throw OpenMMException("Unbalanced Metal optimization template: "+signature);
    return source.substr(begin, end-begin);
}

/** @brief Replace only a complete, unchanged function fingerprint. */
bool replaceFunction(string& source, const string& original, const string& replacement) {
    size_t begin = source.find(original);
    if (begin == string::npos)
        return false;
    source.replace(begin, original.size(), replacement);
    return true;
}

/** @brief Reduce within SIMD groups, then their totals, with shared scratch reuse protected. */
string reductionFunction(const string& name) {
    return "DEVICE real "+name+R"((real value, LOCAL_ARG volatile real* temp) {
    const int lane = LOCAL_ID%32;
    const int warp = LOCAL_ID/32;
    const int active = min(32, (int) LOCAL_SIZE-32*warp);
    const int warps = ((int) LOCAL_SIZE+31)/32;
    for (int offset = 16; offset > 0; offset /= 2) {
        real other = simd_shuffle_down(value, offset);
        if (lane+offset < active)
            value += other;
    }
    SYNC_THREADS;
    if (lane == 0)
        temp[warp] = value;
    SYNC_THREADS;
    value = (LOCAL_ID < warps ? temp[LOCAL_ID] : (real) 0);
    if (warp == 0) {
        for (int offset = 16; offset > 0; offset /= 2) {
            real other = simd_shuffle_down(value, offset);
            if (lane+offset < min(32, (int) LOCAL_SIZE))
                value += other;
        }
        if (lane == 0)
            temp[0] = value;
    }
    SYNC_THREADS;
    real result = temp[0];
    SYNC_THREADS;
    return result;
}
)";
}

} // namespace

MetalReductionOptimizations::Settings MetalReductionOptimizations::getBuildSettings() {
    Settings settings;
#if OPENMM_METAL_FAST_CENTROID_REDUCTION
    settings.centroid = true;
#endif
    return settings;
}

string MetalReductionOptimizations::apply(const string& source) {
    return apply(source, getBuildSettings());
}

string MetalReductionOptimizations::apply(const string& source, const Settings& settings) {
    string result = source;
    if (settings.centroid) {
        string original = functionText(CommonKernelSources::customCentroidBond, "KERNEL void computeGroupCenters(");
        string replacement = original.substr(0, original.find("        // Sum the values."));
        size_t temp = replacement.find("LOCAL volatile real3 temp[64]");
        replacement.replace(temp, string("LOCAL volatile real3 temp[64]").size(), "LOCAL volatile real temp[64]");
        replacement += R"(
        center.x = metalReduceGroupValue(center.x, temp);
        center.y = metalReduceGroupValue(center.y, temp);
        center.z = metalReduceGroupValue(center.z, temp);
        if (LOCAL_ID == 0)
            centerPositions[group] = make_real4(center.x, center.y, center.z, 0);
    }
})";
        replaceFunction(result, original, reductionFunction("metalReduceGroupValue")+replacement);
    }
    return result;
}
