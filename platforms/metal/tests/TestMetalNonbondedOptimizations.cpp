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
#include "MetalNonbondedSources.h"
#include "MetalOpenCLKernelSources.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

/** @brief Compare all four stored components, including radius or charge extent. */
void assertBoundsEqual(const mm_float4& expected, const mm_float4& actual) {
    ASSERT_EQUAL_TOL(expected.x, actual.x, 1e-6);
    ASSERT_EQUAL_TOL(expected.y, actual.y, 1e-6);
    ASSERT_EQUAL_TOL(expected.z, actual.z, 1e-6);
    ASSERT_EQUAL_TOL(expected.w, actual.w, 1e-6);
}

/**
 * @brief Compare the exact production bounds source to the unmodified OpenCL kernel.
 * Cover periodic ordered images (including triclinic shear and atoms more than
 * half a cell apart), short tiles, idle lanes, and grid-stride iterations. Both
 * sources are tested in OFF builds too; the source-selection assertion verifies
 * the independent FAST_BLOCK_BOUNDS switch used by actual neighbor kernels.
 */
void testBlockBounds(int periodicMode) {
    const int maxAtoms = 4097, maxBlocks = (maxAtoms+31)/32;
    System system;
    for (int i = 0; i < maxAtoms; i++) system.addParticle(1);
    Vec3 a(16, 0, 0), b(periodicMode == 2 ? 4 : 0, 16, 0);
    Vec3 c(periodicMode == 2 ? -2 : 0, periodicMode == 2 ? 2 : 0, 16);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    MetalContext context(system);
    map<string, string> defines;
    defines["TILE_SIZE"] = "32";
    defines["OPENMM_METAL_REQUIRE_SAFE_MATH"] = "1";
    if (periodicMode != 0) defines["USE_PERIODIC"] = "1";
    const string& original = MetalOpenCLKernelSources::findInteractingBlocks;
    size_t first = original.find("__kernel void findBlockBounds(");
    size_t last = original.find("__kernel void computeSortKeys(", first);
    ASSERT(first != string::npos && last != string::npos);
    const string selected = MetalNonbondedSources::neighbors(original, false);
    ASSERT_EQUAL(bool(OPENMM_METAL_FAST_BLOCK_BOUNDS),
                 selected.find(MetalKernelSources::neighborBounds) != string::npos);
    ComputeKernel reference = context.compileProgram(original.substr(first, last-first), defines)->createKernel("findBlockBounds");
    ComputeKernel optimized = context.compileProgram(MetalKernelSources::neighborBounds, defines)->createKernel("findBlockBounds");
    vector<mm_float4> positions(maxAtoms);
    for (int i = 0; i < maxAtoms; i++) {
        // The first three atoms deliberately distinguish ordered expansion
        // from wrapping all atoms relative to the first atom of the block.
        const int lane = i%32;
        float x = lane == 0 ? 0.0625f : (lane == 1 ? 0.46875f : (lane == 2 ? 0.75f : float(lane%17)/16));
        float y = float((3*lane)%17)/16, z = float((7*lane)%17)/16;
        x += (i%5)-2;
        y += (i%3)-1;
        z += (i%7)-3;
        Vec3 position = a*x+b*y+c*z;
        positions[i] = mm_float4(position[0], position[1], position[2], float(lane-15)/16);
    }
    ComputeArray posq, centers, bounds, ranges, rebuild;
    posq.initialize<mm_float4>(context, maxAtoms, "boundsPositions");
    centers.initialize<mm_float4>(context, maxBlocks+1, "boundsCenters");
    bounds.initialize<mm_float4>(context, maxBlocks+1, "boundsWidths");
    ranges.initialize<mm_float2>(context, context.getNumThreadBlocks()+1, "boundsSizeRanges");
    rebuild.initialize<int>(context, 1, "boundsRebuild");
    posq.upload(positions);
    for (ComputeKernel kernel : {reference, optimized}) {
        kernel->addArg(maxAtoms);
        kernel->addArg(mm_float4(16, 16, 16, 0));
        kernel->addArg(mm_float4(1.0f/16, 1.0f/16, 1.0f/16, 0));
        kernel->addArg(mm_float4(a[0], a[1], a[2], 0));
        kernel->addArg(mm_float4(b[0], b[1], b[2], 0));
        kernel->addArg(mm_float4(c[0], c[1], c[2], 0));
        kernel->addArg(posq);
        kernel->addArg(centers);
        kernel->addArg(bounds);
        kernel->addArg(rebuild);
        kernel->addArg(ranges);
    }
    const mm_float4 sentinel(-12345, -12345, -12345, -12345);
    for (int atoms : {1, 3, 31, 32, 33, 63, 64, 65, 129, maxAtoms}) {
        const int blocks = (atoms+31)/32;
        for (int groupLimit : {1, context.getNumThreadBlocks()}) {
            vector<mm_float4> referenceCenters, referenceBounds;
            for (bool fast : {false, true}) {
                const int blocksPerGroup = fast && periodicMode == 0 ? 2 : 64;
                const int groups = min((blocks+blocksPerGroup-1)/blocksPerGroup, groupLimit);
                ComputeKernel kernel = fast ? optimized : reference;
                centers.upload(vector<mm_float4>(maxBlocks+1, sentinel));
                bounds.upload(vector<mm_float4>(maxBlocks+1, sentinel));
                ranges.upload(vector<mm_float2>(ranges.getSize(), mm_float2(-12345, -12345)));
                rebuild.upload(vector<int>(1, 1));
                kernel->setArg(0, atoms);
                kernel->execute(groups*64, 64);
                vector<mm_float4> actualCenters, actualBounds;
                vector<mm_float2> actualRanges;
                vector<int> actualRebuild;
                centers.download(actualCenters);
                bounds.download(actualBounds);
                ranges.download(actualRanges);
                rebuild.download(actualRebuild);
                ASSERT_EQUAL(0, actualRebuild[0]);
                assertBoundsEqual(sentinel, actualCenters[blocks]);
                assertBoundsEqual(sentinel, actualBounds[blocks]);
                ASSERT_EQUAL(-12345, actualRanges[groups].x);
                ASSERT_EQUAL(-12345, actualRanges[groups].y);
                if (periodicMode == 1 && atoms == 3) {
                    ASSERT_EQUAL_TOL(6.5, actualCenters[0].x, 1e-6);
                    ASSERT_EQUAL_TOL(5.5, actualBounds[0].x, 1e-6);
                }
                for (int group = 0; group < groups; group++) {
                    float minSize = 1e38f, maxSize = 0;
                    for (int base = group*blocksPerGroup; base < blocks; base += groups*blocksPerGroup) {
                        for (int i = base; i < min(base+blocksPerGroup, blocks); i++) {
                            float size = actualBounds[i].x+actualBounds[i].y+actualBounds[i].z;
                            minSize = min(minSize, size);
                            maxSize = max(maxSize, size);
                        }
                    }
                    ASSERT_EQUAL_TOL(minSize, actualRanges[group].x, 1e-6);
                    ASSERT_EQUAL_TOL(maxSize, actualRanges[group].y, 1e-6);
                }
                if (!fast) {
                    referenceCenters = actualCenters;
                    referenceBounds = actualBounds;
                }
                else {
                    for (int i = 0; i < blocks; i++) {
                        assertBoundsEqual(referenceCenters[i], actualCenters[i]);
                        assertBoundsEqual(referenceBounds[i], actualBounds[i]);
                    }
                }
            }
        }
    }
}

/** @brief Conservative rounding includes values next to half spacing and overflow. */
void testHalfBounds(bool initialFloating) {
    System system;
    system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, initialFloating);
    vector<float> input{0, 1.0e-8f, 0.00006099f, 0.10001f, 1.0001f, 100.01f, 65504, 65505, 1.0e10f};
    ComputeArray in, out, infinite;
    in.initialize<float>(context, input.size(), "halfBoundsInput");
    out.initialize<float>(context, input.size(), "halfBoundsOutput");
    infinite.initialize<int>(context, input.size(), "halfBoundsInfinite");
    in.upload(input);
    string source = MetalNonbondedSources::halfBoundsSource()+R"(
KERNEL void roundBounds(GLOBAL const float* input, GLOBAL float* output, GLOBAL int* infinite, int count) {
    for (int i = GLOBAL_ID; i < count; i += GLOBAL_SIZE) {
        output[i] = float(metalBoundsHalf(make_real4(input[i], 0, 0, 0)).x);
        infinite[i] = isinf(output[i]) ? 1 : 0;
    }
})";
    map<string, string> defines;
    defines["OPENMM_METAL_REQUIRE_SAFE_MATH"] = "1";
    ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("roundBounds");
    kernel->addArg(in);
    kernel->addArg(out);
    kernel->addArg(infinite);
    kernel->addArg((int) input.size());
    // Reuse the same kernel across lazy pipeline variants and then return to
    // the initial ABI. REQUIRE_SAFE_MATH must survive both compilation orders.
    for (bool floating : {initialFloating, !initialFloating, initialFloating}) {
        context.setUseFloatingPointAccumulators(floating);
        kernel->execute(input.size());
        vector<float> result;
        vector<int> flags;
        out.download(result);
        infinite.download(flags);
        for (int i = 0; i < input.size(); i++) {
            ASSERT(result[i] >= input[i]);
            ASSERT_EQUAL(input[i] > 65504 ? 1 : 0, flags[i]);
            if (input[i] <= 65504) {
                ASSERT(result[i]-input[i] <= max(6.0e-8f, input[i]*0.001f));
            }
            else {
                ASSERT(isinf(result[i]));
            }
        }
    }
}

/**
 * @brief Two long blocks contain one neighbor per atom, plus a one-atom tail.
 * Verify bounds, positive sparse-list use when enabled, both accumulator modes,
 * and repeated reuse/rebuild without changing the tiled shuffle option.
 */
void testSparsePairs(bool floating) {
    System system;
    for (int i = 0; i < 65; i++) system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, floating);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    vector<vector<int> > exclusions(65);
    for (int i = 0; i < 65; i++) exclusions[i].push_back(i);
    nb.addInteraction(true, false, true, 1.0, exclusions,
        "if (!isExcluded && r2 < CUTOFF_SQUARED) { tempEnergy = 0.5f*r2; dEdR = -1; }", 0, true, true);
    context.initialize();
    vector<mm_float4> positions(context.getPaddedNumAtoms(), mm_float4(0,0,0,0));
    for (int i = 0; i < 32; i++) {
        positions[i] = mm_float4(10*i, 0, 0, 0);
        positions[i+32] = mm_float4(10*i, 0.2f, 0, 0);
    }
    positions[64] = mm_float4(1000, 0, 0, 0);
    for (int iteration = 0; iteration < 3; iteration++) {
        float distance = iteration == 2 ? 0.4f : 0.2f;
        for (int i = 32; i < 64; i++) positions[i].y = distance;
        context.getPosq().upload(positions);
        context.clearBuffer(context.getLongForceBuffer());
        context.clearBuffer(context.getEnergyBuffer());
        nb.prepareInteractions(1);
        nb.computeInteractions(1, true, true);
        vector<unsigned int> count;
        nb.getInteractionCount().download(count);
#if OPENMM_METAL_FAST_SPARSE_PAIRS
        ASSERT_EQUAL(2, count.size());
        ASSERT_EQUAL(32, count[1]);
        ASSERT_EQUAL(0, count[0]);
#else
        ASSERT_EQUAL(1, count.size());
        ASSERT(count[0] > 0);
#endif
        ASSERT_EQUAL_TOL(16*distance*distance, context.reduceEnergy(), 1e-5);
        vector<double> force;
        context.downloadFixedPointBuffer(context.getLongForceBuffer(), force);
        for (int i = 0; i < 65; i++) {
            ASSERT_EQUAL_TOL(0, force[i], 1e-6);
            ASSERT_EQUAL_TOL(i < 32 ? distance : (i < 64 ? -distance : 0), force[i+context.getPaddedNumAtoms()], 1e-6);
            ASSERT_EQUAL_TOL(0, force[i+2*context.getPaddedNumAtoms()], 1e-6);
        }
        vector<mm_float4> centers, bounds;
        nb.getBlockCenters().download(centers);
        nb.getBlockBoundingBoxes().download(bounds);
        ASSERT_EQUAL_TOL(155, centers[0].x, 1e-6);
        ASSERT_EQUAL_TOL(155, bounds[0].x, 1e-6);
        ASSERT_EQUAL_TOL(1000, centers[2].x, 1e-6);
        ASSERT_EQUAL_TOL(0, bounds[2].x, 1e-6);
    }
}

/** @brief Overflow the initial sparse allocation and retry with cached force kernels. */
void testSparseResize(bool floating) {
    const int blocks = 16, atoms = 32*blocks;
    System system;
    for (int i = 0; i < atoms; i++) system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, floating);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    vector<vector<int> > exclusions(atoms);
    for (int i = 0; i < atoms; i++) exclusions[i].push_back(i);
    nb.addInteraction(true, false, true, 1.0, exclusions,
        "if (!isExcluded && r2 < CUTOFF_SQUARED) { tempEnergy = 0.5f*r2; dEdR = -1; }", 0, true, true);
    context.initialize();
    vector<mm_float4> positions(atoms);
    for (int block = 0; block < blocks; block++)
        for (int lane = 0; lane < 32; lane++)
            positions[32*block+lane] = mm_float4(10*lane, 0.02f*block, 0, 0);
    context.getPosq().upload(positions);
    int retries = 0;
    do {
        ASSERT(retries < 3);
        context.setForcesValid(true);
        context.clearBuffer(context.getLongForceBuffer());
        context.clearBuffer(context.getEnergyBuffer());
        nb.prepareInteractions(1);
        nb.computeInteractions(1, true, true);
        retries++;
    } while (!context.getForcesValid());
#if OPENMM_METAL_FAST_SPARSE_PAIRS
    ASSERT_EQUAL(2, retries);
    vector<unsigned int> count;
    nb.getInteractionCount().download(count);
    ASSERT_EQUAL(32*blocks*(blocks-1)/2, count[1]);
    ASSERT_EQUAL(0, count[0]);
#else
    ASSERT_EQUAL(1, retries);
#endif
    vector<double> forces;
    context.downloadFixedPointBuffer(context.getLongForceBuffer(), forces);
    double energy = 0;
    for (int block = 0; block < blocks; block++) {
        double fy = 0;
        for (int other = 0; other < blocks; other++) {
            double delta = (double) positions[32*other].y-positions[32*block].y;
            fy += delta;
            energy += 32*0.25*delta*delta;
        }
        for (int lane = 0; lane < 32; lane++) {
            ASSERT_EQUAL_TOL(0, forces[32*block+lane], 1e-6);
            ASSERT_EQUAL_TOL(fy, forces[atoms+32*block+lane], 1e-5);
            ASSERT_EQUAL_TOL(0, forces[2*atoms+32*block+lane], 1e-6);
        }
    }
    ASSERT_EQUAL_TOL(energy, context.reduceEnergy(), 1e-5);
}

/** @brief Host storage for the Common Compute packed, three-component ABI. */
struct PackedParameter3 {
    float x, y, z;
};

/** @brief Count source markers without depending on generated argument names. */
int countSourceOccurrences(const string& source, const string& marker) {
    int count = 0;
    for (size_t offset = source.find(marker); offset != string::npos; offset = source.find(marker, offset+marker.size()))
        count++;
    return count;
}

/** @brief Hybrid adaptation changes read-only exchange, not force writes or barriers. */
void testHybridSourceContract() {
    const string& original = MetalOpenCLKernelSources::nonbonded;
    const string hybrid = MetalNonbondedSources::openclHybridSource(original);
    ASSERT(hybrid.find("real4 hybridPosq = make_real4(0);") != string::npos);
    ASSERT_EQUAL(4, countSourceOccurrences(hybrid, "simdShuffle(hybridPosq.x, uint(atom2-tbx))"));
    ASSERT_EQUAL(4, countSourceOccurrences(hybrid, "LOAD_ATOM2_PARAMETERS"));
    ASSERT(hybrid.find("localData[atom2].x") == string::npos);
    ASSERT(hybrid.find("APPLY_PERIODIC_TO_POS_WITH_CENTER(hybridPosq, blockCenterX)") != string::npos);
    for (const string& marker : {"SYNC_WARPS;", "ATOMIC_ADD(&forceBuffers[", "realToFixedPoint(",
                                 "localData[tbx+tj].fx += delta.x;", "atomIndices[tbx+tj]", "skipTiles["})
        ASSERT_EQUAL(countSourceOccurrences(original, marker), countSourceOccurrences(hybrid, marker));
}

/**
 * @brief Check generated tiled forces against an independent all-pairs calculation.
 *
 * Every parameter component and both halves of the fixed-ABI 64-bit parameter
 * affect the answer.  This catches lane/atom confusion for diagonal, excluded, ordinary,
 * and noncontiguous neighbor tiles.  Parameters and positions have nonzero
 * padded entries, so missing tail masks cannot silently use zero-valued data.
 * The active build selects the baseline, rotating, or hybrid production path;
 * no test-only force kernel replaces that selection.
 *
 * @param atoms Particle count, including a deliberately incomplete final tile.
 * @param geometry 0: no cutoff; 1: nonperiodic cutoff; 2: compact orthorhombic
 *                 images; 3: triclinic positions spanning the periodic cell.
 * @param floating Whether to exercise the floating accumulator ABI as well.
 */
void testTiledParameterGather(int atoms, int geometry, bool floating) {
    cout << "Tiled parameter gather: atoms=" << atoms << ", geometry=" << geometry
         << ", floating=" << floating << endl;
    static_assert(sizeof(PackedParameter3) == 3*sizeof(float), "Packed parameter stride must be 12 bytes");
    const bool cutoff = geometry != 0, periodic = geometry >= 2;
    const double cutoffDistance = 2.5;
    const Vec3 a(8, 0, 0), b(geometry == 3 ? 2 : 0, 8, 0), c(geometry == 3 ? -1 : 0, geometry == 3 ? 1 : 0, 8);
    System system;
    for (int i = 0; i < atoms; i++) system.addParticle(1);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    MetalContext context(system, nullptr, nullptr, floating);
    context.setPeriodicBoxVectors(a, b, c);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    const int padded = context.getPaddedNumAtoms();
    ComputeArray scalar, pair, triple, quad, stamp;
    scalar.initialize<float>(context, padded, "gatherScalar");
    pair.initialize<mm_float2>(context, padded, "gatherPair");
    triple.initialize<PackedParameter3>(context, padded, "gatherTriple");
    quad.initialize<mm_float4>(context, padded, "gatherQuad");
    // The existing floating-accumulator adapter rewrites every long/ulong
    // pointer as a float pointer: arbitrary integer buffers are not supported
    // by that ABI. Test exact 64-bit transport in the fixed ABI, and explicitly
    // supply its equivalent float weight in the floating ABI. This is not a
    // relaxation of the oracle or a skipped interaction/parameter component.
    stamp.initialize(context, padded, floating ? sizeof(float) : sizeof(uint64_t), "gatherStamp");
    nb.addParameter(ComputeParameterInfo(scalar, "gatherScalar", "float", 1));
    nb.addParameter(ComputeParameterInfo(pair, "gatherPair", "float", 2));
    nb.addParameter(ComputeParameterInfo(triple, "gatherTriple", "float", 3));
    nb.addParameter(ComputeParameterInfo(quad, "gatherQuad", "float", 4));
    nb.addParameter(ComputeParameterInfo(stamp, "gatherStamp", floating ? "float" : "ulong", 1));
    vector<vector<int> > exclusions(atoms);
    for (int i = 0; i < atoms; i++) exclusions[i].push_back(i);
    auto exclude = [&](int i, int j) {
        exclusions[i].push_back(j);
        exclusions[j].push_back(i);
    };
    exclude(2, 17); // Diagonal exclusion.
    if (atoms > 34) exclude(5, 34); // Off-diagonal exclusion tile.
    nb.addInteraction(cutoff, periodic, true, cutoffDistance, exclusions, R"(
if (!isExcluded
#ifdef USE_CUTOFF
        && r2 < CUTOFF_SQUARED
#endif
        ) {
#ifdef OPENMM_METAL_FLOAT_ACCUMULATORS
    real stampWeight1 = gatherStamp1;
    real stampWeight2 = gatherStamp2;
#else
    real stampWeight1 = (real) ((uint) (gatherStamp1 >> 32))*(1.0f/256)
        +(real) ((uint) gatherStamp1)*(1.0f/8192);
    real stampWeight2 = (real) ((uint) (gatherStamp2 >> 32))*(1.0f/256)
        +(real) ((uint) gatherStamp2)*(1.0f/8192);
#endif
    real weight1 = gatherScalar1+posq1.w+gatherPair1.x+2*gatherPair1.y
        +gatherTriple1.x+2*gatherTriple1.y+3*gatherTriple1.z
        +gatherQuad1.x+2*gatherQuad1.y+3*gatherQuad1.z+4*gatherQuad1.w
        +stampWeight1;
    real weight2 = gatherScalar2+posq2.w+gatherPair2.x+2*gatherPair2.y
        +gatherTriple2.x+2*gatherTriple2.y+3*gatherTriple2.z
        +gatherQuad2.x+2*gatherQuad2.y+3*gatherQuad2.z+4*gatherQuad2.w
        +stampWeight2;
    real k = 0.125f*(weight1+weight2);
    tempEnergy = 0.5f*k*r2;
    dEdR = -k;
})", 0, true, false); // Keep these checks on tiled, rather than sparse-pair, kernels.
    context.initialize();
    vector<float> scalars(padded, 99);
    vector<mm_float2> pairs(padded, mm_float2(99, 99));
    vector<PackedParameter3> triples(padded, PackedParameter3{99, 99, 99});
    vector<mm_float4> quads(padded, mm_float4(99, 99, 99, 99));
    vector<uint64_t> stamps(padded, uint64_t(99)<<32 | 99);
    vector<float> floatingStamps(padded, 99.0f/256+99.0f/8192);
    vector<mm_float4> positions(padded, mm_float4(0.125f, 0.25f, 0.375f, 99));
    for (int iteration = 0; iteration < 2; iteration++) {
        vector<double> weights(atoms);
        for (int i = 0; i < atoms; i++) {
            scalars[i] = float((i+3*iteration)%17+1)/64;
            pairs[i] = mm_float2(float((3*i)%19+1)/128, float((5*i)%23+1)/128);
            triples[i] = PackedParameter3{float(i%11+1)/256, float(i%13+1)/256, float(i%7+1)/256};
            quads[i] = mm_float4(float(i%5+1)/256, float(i%7+1)/256, float(i%11+1)/256, float(i%13+1)/256);
            stamps[i] = uint64_t(i+1+iteration)<<32 | uint64_t((13*i+5)%251);
            // Both terms are binary fractions exactly representable at these
            // magnitudes, so this preconversion preserves the reference weight.
            floatingStamps[i] = float(double(stamps[i]>>32)/256+double(uint32_t(stamps[i]))/8192);
            Vec3 p;
            if (geometry == 3)
                p = a*(double((17*i+3*iteration)%101)/101)+b*(double((29*i)%103)/103)+c*(double((43*i)%107)/107);
            else {
                p = Vec3(0.19*(i%7), 0.21*((i/7)%5), 0.23*(i/35));
                p[0] += 0.015*iteration*(i%3);
                if (geometry == 1)
                    p *= 3; // Mix accepted and cutoff-rejected pairs within tiles.
                if (geometry == 2)
                    p += a*(i%3-1)+b*(i%5-2)+c*(i%7-3);
            }
            positions[i] = mm_float4(p[0], p[1], p[2], float(i%9+1)/64);
            weights[i] = scalars[i]+double(positions[i].w)+pairs[i].x+2*double(pairs[i].y)
                +triples[i].x+2*double(triples[i].y)+3*double(triples[i].z)
                +quads[i].x+2*double(quads[i].y)+3*double(quads[i].z)+4*double(quads[i].w)
                +double(stamps[i]>>32)/256+double(uint32_t(stamps[i]))/8192;
        }
        scalar.upload(scalars);
        pair.upload(pairs);
        triple.upload(triples);
        quad.upload(quads);
        if (floating)
            stamp.upload(floatingStamps);
        else
            stamp.upload(stamps);
        context.getPosq().upload(positions);
        vector<Vec3> expectedForces(padded);
        double expectedEnergy = 0;
        for (int i = 0; i < atoms; i++) {
            for (int j = i+1; j < atoms; j++) {
                if (find(exclusions[i].begin(), exclusions[i].end(), j) != exclusions[i].end()) continue;
                Vec3 delta(double(positions[j].x)-positions[i].x, double(positions[j].y)-positions[i].y,
                           double(positions[j].z)-positions[i].z);
                if (periodic) {
                    delta -= c*floor(delta[2]/c[2]+0.5);
                    delta -= b*floor(delta[1]/b[1]+0.5);
                    delta -= a*floor(delta[0]/a[0]+0.5);
                }
                double r2 = delta.dot(delta);
                if (cutoff && r2 >= cutoffDistance*cutoffDistance) continue;
                double k = 0.125*(weights[i]+weights[j]);
                expectedEnergy += 0.5*k*r2;
                expectedForces[i] += delta*k;
                expectedForces[j] -= delta*k;
            }
        }
        nb.prepareInteractions(1);
        if (cutoff) {
            vector<unsigned int> count;
            vector<int> neighborAtoms;
            nb.getInteractionCount().download(count);
            nb.getInteractingAtoms().download(neighborAtoms);
            ASSERT_EQUAL(1, count.size());
            ASSERT(count[0] > 0);
            ASSERT(count[0]*32 <= neighborAtoms.size());
            // A bijection preserves each tile's interaction set while forcing
            // lane IDs to differ from both particle IDs and contiguous offsets.
            vector<int> original = neighborAtoms;
            for (unsigned int tile = 0; tile < count[0]; tile++)
                for (int lane = 0; lane < 32; lane++)
                    neighborAtoms[32*tile+lane] = original[32*tile+(13*lane+7)%32];
            nb.getInteractingAtoms().upload(neighborAtoms);
        }
        for (int outputs = 0; outputs < 3; outputs++) {
            bool includeForces = outputs != 1, includeEnergy = outputs != 2;
            stringstream caseDescription;
            caseDescription << "Tiled parameter gather: atoms=" << atoms << ", geometry=" << geometry
                << ", floating=" << floating << ", iteration=" << iteration << ", outputs=" << outputs;
            context.setForcesValid(true);
            context.clearBuffer(context.getLongForceBuffer());
            context.clearBuffer(context.getEnergyBuffer());
            nb.computeInteractions(1, includeForces, includeEnergy);
            ASSERT(context.getForcesValid());
            try {
                ASSERT_EQUAL_TOL(includeEnergy ? expectedEnergy : 0, context.reduceEnergy(), 2e-5);
            }
            catch (const exception& error) {
                throw OpenMMException(caseDescription.str()+", energy: "+error.what());
            }
            vector<double> forces;
            context.downloadFixedPointBuffer(context.getLongForceBuffer(), forces);
            for (int atom = 0; atom < padded; atom++) {
                for (int axis = 0; axis < 3; axis++) {
                    try {
                        ASSERT_EQUAL_TOL(includeForces ? expectedForces[atom][axis] : 0, forces[atom+axis*padded], 3e-5);
                    }
                    catch (const exception& error) {
                        stringstream description;
                        description << caseDescription.str() << ", atom=" << atom << ", axis=" << axis << ": " << error.what();
                        throw OpenMMException(description.str());
                    }
                }
            }
        }
    }
}

/** @brief Exercise the separate large-block argument layout and FP16 consumers. */
void testLargeBlocks() {
    const int atoms = 100001; // The same threshold used by the OpenCL utility.
    System system;
    for (int i = 0; i < atoms; i++) system.addParticle(1);
    MetalContext context(system);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    vector<vector<int> > exclusions;
    nb.addInteraction(true, false, false, 1.0, exclusions, "", 0, true, true);
    context.initialize();
    vector<mm_float4> positions(context.getPaddedNumAtoms(), mm_float4(0,0,0,0));
    for (int i = 0; i < atoms; i++) positions[i].x = 2*i;
    context.getPosq().upload(positions);
    nb.prepareInteractions(1);
    vector<unsigned int> counts;
    nb.getInteractionCount().download(counts);
    ASSERT_EQUAL(0, counts[0]);
    if (counts.size() == 2) { ASSERT_EQUAL(0, counts[1]); }
    vector<mm_float4> centers, bounds;
    nb.getBlockCenters().download(centers);
    nb.getBlockBoundingBoxes().download(bounds);
    ASSERT_EQUAL_TOL(200000, centers.back().x, 1e-6);
    ASSERT_EQUAL_TOL(0, bounds.back().x, 1e-6);
}

int main() {
    try {
        testBlockBounds(0);
        testBlockBounds(1);
        testBlockBounds(2);
        testHalfBounds(false);
        testHalfBounds(true);
        testSparsePairs(false);
        testSparsePairs(true);
        testSparseResize(false);
        testSparseResize(true);
        testHybridSourceContract();
        for (bool floating : {false, true}) {
            testTiledParameterGather(31, 0, floating);
            testTiledParameterGather(97, 0, floating);
            for (int geometry : {1, 2, 3})
                testTiledParameterGather(97, geometry, floating);
        }
        testLargeBlocks();
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
