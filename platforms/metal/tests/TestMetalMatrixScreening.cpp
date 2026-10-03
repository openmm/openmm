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
 * This program is distributed WITHOUT ANY WARRANTY; without even the          *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   *
 * See the GNU Lesser General Public License for more details.                 *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.          *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalKernelSources.h"
#include "MetalQueue.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

using namespace OpenMM;
using namespace std;

/** @brief Mirror the shader ABI without depending on host float3 alignment. */
struct MatrixScreenParameters {
    mm_float4 boxX, boxY, boxZ;
    float cutoffSquared;
    uint32_t periodic, tileCount, reserved;
};
static_assert(sizeof(MatrixScreenParameters) == 64, "Matrix-screen parameter size");
static_assert(offsetof(MatrixScreenParameters, cutoffSquared) == 48, "Matrix-screen scalar offset");

/** @brief Count diagnostic mask disagreements; candidate masks are NOT authoritative. */
unsigned int bitCount(uint32_t value) {
    unsigned int count = 0;
    for (; value; value &= value-1) count++;
    return count;
}

/** @brief Keep compilation, uploads, and allocation out of resident-data timings. */
class MatrixScreenFixture {
public:
    MetalContext& context;
    ComputeArray first, second, exact, candidates, distances, parameters, shape;
    ComputeKernel scalar, matrix;
    MatrixScreenParameters params;

    MatrixScreenFixture(MetalContext& context, int tiles, bool diagnostics) : context(context) {
        params = {mm_float4(8, 0, 0, 0), mm_float4(0, 8, 0, 0), mm_float4(0, 0, 8, 0),
                1.0f, 0, uint32_t(tiles), 0};
        first.initialize<mm_float4>(context, 32*tiles, "matrixFirst");
        second.initialize<mm_float4>(context, 32*tiles, "matrixSecond");
        exact.initialize<uint32_t>(context, 32*tiles, "matrixExactMasks");
        candidates.initialize<uint32_t>(context, 32*tiles, "matrixCandidateMasks");
        distances.initialize<float>(context, diagnostics ? 1024*tiles : 1, "matrixDistances");
        parameters.initialize(context, 1, sizeof(params), "matrixParameters");
        shape.initialize<uint32_t>(context, 2, "matrixExecutionShape");
        map<string, string> defines;
        defines["OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN"] = "1";
        defines["METAL_MATRIX_DIAGNOSTICS"] = diagnostics ? "1" : "0";
        ComputeProgram program = context.compileProgram(MetalKernelSources::neighborMatrix, defines);
        scalar = program->createKernel("scalarScreen");
        matrix = program->createKernel("matrixScreen");
        for (ComputeKernel kernel : {scalar, matrix}) {
            ASSERT(kernel->getMaxBlockSize() >= 32);
            kernel->addArg(first);
            kernel->addArg(second);
            kernel->addArg(exact);
            kernel->addArg(candidates);
            kernel->addArg(distances);
            kernel->addArg(parameters);
            kernel->addArg(shape);
        }
    }

    void upload(const vector<mm_float4>& a, const vector<mm_float4>& b) {
        first.upload(a);
        second.upload(b);
        parameters.upload(&params);
    }

    void execute(ComputeKernel kernel) {
        kernel->execute(32*params.tileCount, 32);
    }

    /** @brief Assert exact final-mask agreement and report, not hide, candidate errors. */
    void verify(const string& label, bool exactMatrixDistances=false) {
        vector<uint32_t> reference, result, candidate;
        vector<float> scalarDistance, matrixDistance;
        execute(scalar);
        vector<uint32_t> dimensions;
        shape.download(dimensions);
        ASSERT_EQUAL(32, dimensions[0]);
        ASSERT_EQUAL(32, dimensions[1]);
        exact.download(reference);
        if (exactMatrixDistances) distances.download(scalarDistance);
        execute(matrix);
        shape.download(dimensions);
        ASSERT_EQUAL(32, dimensions[0]);
        ASSERT_EQUAL(32, dimensions[1]);
        exact.download(result);
        candidates.download(candidate);
        ASSERT_EQUAL_CONTAINERS(reference, result);
        if (exactMatrixDistances) {
            distances.download(matrixDistance);
            ASSERT_EQUAL_CONTAINERS(scalarDistance, matrixDistance);
        }
        unsigned int falseNegative = 0, falsePositive = 0, neighbors = 0;
        for (size_t i = 0; i < reference.size(); i++) {
            falseNegative += bitCount(reference[i]&~candidate[i]);
            falsePositive += bitCount(candidate[i]&~reference[i]);
            neighbors += bitCount(reference[i]);
        }
        cout << "{\"case\":\"" << label << "\",\"pairs\":" << 1024*params.tileCount
             << ",\"scalar_neighbors\":" << neighbors << ",\"candidate_false_negatives\":" << falseNegative
             << ",\"candidate_false_positives\":" << falsePositive << ",\"final_mask_mismatches\":0}" << endl;
    }
};

/** @brief Known integer geometry detects matrix orientation and norm-packing defects. */
void testIntegerGeometry(MatrixScreenFixture& fixture) {
    vector<mm_float4> a(fixture.first.getSize()), b(a.size());
    for (int i = 0; i < a.size(); i++) {
        a[i] = mm_float4(i%8, (i/8)%4, i%3, 0);
        b[i] = mm_float4((3*i)%7, i%5, (2*i)%3, 0);
    }
    fixture.upload(a, b);
    fixture.verify("integer_geometry", true);
}

/** @brief Exercise strict cutoff decisions immediately below, at, and above one. */
void testCutoffBoundary(MatrixScreenFixture& fixture) {
    vector<mm_float4> a(fixture.first.getSize(), mm_float4(0, 0, 0, 0)), b(a.size());
    float below = nextafter(1.0f, 0.0f), above = nextafter(1.0f, 2.0f);
    for (int i = 0; i < b.size(); i++)
        b[i] = mm_float4(i%3 == 0 ? below : i%3 == 1 ? 1.0f : above, 0, 0, 0);
    fixture.upload(a, b);
    fixture.verify("cutoff_boundary");
    vector<uint32_t> masks;
    fixture.exact.download(masks);
    for (int i = 0; i < masks.size(); i++) {
        uint32_t expected = 0;
        for (int j = 0; j < 32; j++)
            if ((32*(i/32)+j)%3 == 0) expected |= 1u<<j;
        ASSERT_EQUAL(expected, masks[i]);
    }
}

/** @brief Cover centered and cancellation-prone coordinates, plus periodic image changes. */
void testRandomAndPeriodic(MatrixScreenFixture& fixture) {
    mt19937 random(20260930);
    uniform_real_distribution<float> value(-4.0f, 4.0f);
    vector<mm_float4> a(fixture.first.getSize()), b(a.size());
    for (int scenario = 0; scenario < 5; scenario++) {
        fixture.params.periodic = scenario >= 3;
        fixture.params.boxY.x = scenario == 4 ? 2.0f : 0.0f;
        fixture.params.boxZ.x = scenario == 4 ? 1.0f : 0.0f;
        fixture.params.boxZ.y = scenario == 4 ? 2.0f : 0.0f;
        for (int i = 0; i < a.size(); i++) {
            float origin = scenario == 1 ? 1048576.0f : 0.0f;
            float scale = scenario == 2 ? 1024.0f : 1.0f;
            a[i] = mm_float4(origin+scale*value(random), origin+scale*value(random), origin+scale*value(random), 0);
            b[i] = mm_float4(origin+scale*value(random), origin+scale*value(random), origin+scale*value(random), 0);
            // The wide tile intentionally contains close pairs far from its origin.
            if (scenario == 2) b[i] = mm_float4(a[i].x+0.5f, a[i].y, a[i].z, 0);
        }
        fixture.upload(a, b);
        const char* labels[] = {"random", "large_origin", "wide_tile_cancellation", "orthorhombic", "triclinic"};
        fixture.verify(labels[scenario]);
    }
}

/** @brief Exact binary coordinates verify periodic image shifts and strict cutoffs independently. */
void testKnownPeriodicImages(MatrixScreenFixture& fixture) {
    vector<mm_float4> a(fixture.first.getSize(), mm_float4(3.5f, -3.5f, 3.5f, 0)), b(a.size());
    for (int triclinic = 0; triclinic < 2; triclinic++) {
        fixture.params.periodic = 1;
        fixture.params.boxY = mm_float4(triclinic ? 2 : 0, 8, 0, 0);
        fixture.params.boxZ = mm_float4(triclinic ? 1 : 0, triclinic ? 2 : 0, 8, 0);
        for (int i = 0; i < b.size(); i++) {
            int image = (i%32)/4-4;
            b[i] = mm_float4(3.5f+0.5f*(i%4)+image*(8+fixture.params.boxY.x+fixture.params.boxZ.x),
                    -3.5f+image*(8+fixture.params.boxZ.y), 3.5f+image*8, 0);
        }
        fixture.upload(a, b);
        fixture.verify(triclinic ? "known_triclinic_images" : "known_orthorhombic_images");
        vector<uint32_t> masks;
        fixture.exact.download(masks);
        for (uint32_t mask : masks) ASSERT_EQUAL(0x33333333u, mask);
    }
}

/**
 * @brief Time resident-data end-to-end launches; includes packing, all MMAs,
 * masks, exact rechecks, host encoding, submission, and completion waits.
 * Compilation, allocations and identical coordinate uploads are excluded.
 */
void benchmark(MetalContext& context) {
    const int iterations = 200;
    MatrixScreenFixture fixture(context, 1024, false);
    mt19937 random(271828);
    uniform_real_distribution<float> value(-4.0f, 4.0f);
    vector<mm_float4> a(fixture.first.getSize()), b(a.size());
    for (int i = 0; i < a.size(); i++) {
        a[i] = mm_float4(value(random), value(random), value(random), 0);
        b[i] = mm_float4(value(random), value(random), value(random), 0);
    }
    for (int periodic = 0; periodic < 2; periodic++) {
        fixture.params.periodic = periodic;
        fixture.params.boxY.x = periodic ? 2.0f : 0.0f;
        fixture.params.boxZ.y = periodic ? 2.0f : 0.0f;
        fixture.upload(a, b);
        fixture.verify(periodic ? "benchmark_triclinic" : "benchmark_nonperiodic");
        for (ComputeKernel kernel : {fixture.scalar, fixture.matrix}) {
            for (int i = 0; i < 3; i++) fixture.execute(kernel);
            context.getCurrentMetalQueue().finish();
        }
        for (int repeat = 0; repeat < 5; repeat++) {
            // Alternate ordering to reduce systematic cache/thermal ordering bias.
            vector<ComputeKernel> kernels{fixture.scalar, fixture.matrix};
            if (repeat%2) reverse(kernels.begin(), kernels.end());
            for (ComputeKernel kernel : kernels) {
                auto start = chrono::steady_clock::now();
                for (int i = 0; i < iterations; i++) fixture.execute(kernel);
                context.getCurrentMetalQueue().finish();
                double elapsed = chrono::duration<double>(chrono::steady_clock::now()-start).count();
                cout << setprecision(9) << "{\"benchmark\":\"" << kernel->getName()
                     << "\",\"periodic\":" << periodic << ",\"repeat\":" << repeat
                     << ",\"tiles\":1024,\"iterations\":" << iterations << ",\"wall_seconds\":" << elapsed << "}" << endl;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 2 || (argc == 2 && string(argv[1]) != "--benchmark"))
            throw OpenMMException("Usage: TestMetalMatrixScreening [--benchmark]");
        System system;
        system.addParticle(1);
        MetalContext context(system);
        MatrixScreenFixture fixture(context, 17, true);
        testIntegerGeometry(fixture);
        testCutoffBoundary(fixture);
        testRandomAndPeriodic(fixture);
        testKnownPeriodicImages(fixture);
        if (argc == 2) benchmark(context);
        cout << "Done (experimental only; final masks always use scalar rechecks)" << endl;
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}
