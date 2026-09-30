/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Uses the original OpenMM Common pairwise kernel templates.                 *
 * Original OpenMM code:                                                      *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.       *
 * Authors: Peter Eastman                                                     *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalPairwiseOptimizations.h"
#include "MetalSourceAdapter.h"
#include "CommonKernelSources.h"
#include "openmm/System.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

using namespace OpenMM;
using namespace std;

typedef MetalPairwiseOptimizations::Settings Settings;

/** Build switches are independent and do not recognize a kernel by name alone. */
void testSelection() {
    const vector<string> templates = {CommonKernelSources::customGBValueN2, CommonKernelSources::customGBEnergyN2,
        CommonKernelSources::gbsaObc, CommonKernelSources::dpd, CommonKernelSources::customHbondForce};
    for (const string& source : templates)
        ASSERT_EQUAL(source, MetalPairwiseOptimizations::apply(source, Settings()));
    vector<Settings> options(3);
    options[0].customGBValue = true;
    options[1].customGBEnergy = true;
    options[2].dpdParticles = true;
    const int selected[] = {0,1,3};
    for (int i = 0; i < options.size(); i++) {
        for (int j = 0; j < templates.size(); j++) {
            string source = MetalPairwiseOptimizations::apply(templates[j], options[i]);
            ASSERT_EQUAL(j == selected[i], source != templates[j]);
            if (source != templates[j]) {
                ASSERT(source.find("simdShuffle(") != string::npos);
                ASSERT(source.find("#define USE_HIP") == string::npos);
                ASSERT(source.find("#define __CUDA_ARCH__") == string::npos);
            }
        }
        const string unrelated = "KERNEL void computeN2Value(GLOBAL float* output) { output[GLOBAL_ID] = 0; }\n";
        ASSERT_EQUAL(unrelated, MetalPairwiseOptimizations::apply(unrelated, options[i]));
        string changed = templates[selected[i]];
        changed.insert(changed.find("KERNEL void"), "// modified template skeleton\n");
        ASSERT_EQUAL(changed, MetalPairwiseOptimizations::apply(changed, options[i]));
    }
}

/** Exact Common value template with asymmetric parameters and a parameter derivative. */
void testCustomGBValue(int count, bool floating) {
    System system;
    for (int i = 0; i < count; i++) system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, floating);
    const int padded = 32*((count+31)/32), blocks = padded/32;
    map<string,string> substitutions;
    substitutions["PARAMETER_ARGUMENTS"] = ", GLOBAL const float4* global_params1, GLOBAL mm_ulong* global_dValue0dParam1";
    substitutions["ATOM_PARAMETER_DATA"] = "LOCAL float4 local_params1[LOCAL_BUFFER_SIZE];\nLOCAL real local_dValue0dParam1[LOCAL_BUFFER_SIZE];\n";
    substitutions["LOAD_ATOM1_PARAMETERS"] = "float4 params11 = global_params1[atom1];\nreal dValue0dParam1 = 0;\n";
    substitutions["LOAD_LOCAL_PARAMETERS_FROM_1"] = "local_params1[localAtomIndex] = params11;\n";
    substitutions["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = "local_params1[localAtomIndex] = global_params1[j];\nlocal_dValue0dParam1[localAtomIndex] = 0;\n";
    substitutions["LOAD_ATOM2_PARAMETERS"] = "float4 params12 = local_params1[atom2];\nreal temp_dValue0dParam1_1 = 0, temp_dValue0dParam1_2 = 0;\n";
    substitutions["COMPUTE_VALUE"] = "tempValue1 = params11.x+2*params12.x+0.25f*r;\n"
            "tempValue2 = params12.x+2*params11.x+0.25f*r;\n"
            "temp_dValue0dParam1_1 = params12.y;\ntemp_dValue0dParam1_2 = params11.y;\n";
    substitutions["ADD_TEMP_DERIVS1"] = "dValue0dParam1 += temp_dValue0dParam1_1;\n";
    substitutions["ADD_TEMP_DERIVS2"] = "local_dValue0dParam1[tbx+tj] += temp_dValue0dParam1_2;\n";
    substitutions["STORE_PARAM_DERIVS1"] = "ATOMIC_ADD(&global_dValue0dParam1[offset1], (mm_ulong) realToFixedPoint(dValue0dParam1));\n";
    substitutions["STORE_PARAM_DERIVS2"] = "ATOMIC_ADD(&global_dValue0dParam1[offset2], (mm_ulong) realToFixedPoint(local_dValue0dParam1[LOCAL_ID]));\n";
    string source = context.replaceStrings(CommonKernelSources::customGBValueN2, substitutions);
    Settings settings;
    settings.customGBValue = true;
    const string fast = MetalPairwiseOptimizations::apply(source, settings);
    ASSERT(fast != source);
    vector<mm_int2> tiles;
    for (int i = 0; i < blocks; i++) tiles.push_back(mm_int2(i,i));
    if (blocks > 1) tiles.push_back(mm_int2(1,0));
    sort(tiles.begin(), tiles.end(), [blocks](const mm_int2& a, const mm_int2& b) {
        return a.x+a.y*blocks-a.y*(a.y+1)/2 < b.x+b.y*blocks-b.y*(b.y+1)/2;
    });
    vector<unsigned int> exclusions(tiles.size()*32, ~0u);
    for (int t = 0; t < tiles.size(); t++)
        for (int i = 0; i < 32; i++)
            for (int j = 0; j < 32; j++) {
                const int a = 32*tiles[t].x+i, b = 32*tiles[t].y+j;
                if (a == b || (a == 0 && b == 31) || (a == 31 && b == 0)) exclusions[t*32+i] &= ~(1u<<j);
            }
    map<string,string> defines;
    defines["NUM_ATOMS"] = to_string(count);
    defines["PADDED_NUM_ATOMS"] = to_string(padded);
    defines["NUM_BLOCKS"] = to_string(blocks);
    defines["TILE_SIZE"] = "32";
    defines["LOCAL_BUFFER_SIZE"] = "64";
    defines["FIRST_EXCLUSION_TILE"] = "0";
    defines["LAST_EXCLUSION_TILE"] = to_string(tiles.size());
    defines["NUM_TILES_WITH_EXCLUSIONS"] = to_string(tiles.size());
    defines["USE_EXCLUSIONS"] = "1";
    ComputeArray positions, params, exclusionBits, exclusionTiles, output, derivative;
    positions.initialize<mm_float4>(context, padded, "pairwisePositions");
    params.initialize<mm_float4>(context, padded, "pairwiseParams");
    exclusionBits.initialize<unsigned int>(context, exclusions.size(), "pairwiseExclusions");
    exclusionTiles.initialize<mm_int2>(context, tiles.size(), "pairwiseExclusionTiles");
    output.initialize<int64_t>(context, padded, "pairwiseValues");
    derivative.initialize<int64_t>(context, padded, "pairwiseDerivatives");
    vector<mm_float4> pos(padded), par(padded);
    for (int i = 0; i < padded; i++) {
        pos[i] = mm_float4(0.125f*i, 0, 0, 0);
        par[i] = mm_float4(0.25f*(i+1), 0.5f*(i+1), 0, 0);
    }
    positions.upload(pos);
    params.upload(par);
    exclusionBits.upload(exclusions);
    exclusionTiles.upload(tiles);
    for (const string& program : {"// unoptimized reference\n"+source, fast}) {
        context.clearBuffer(output);
        context.clearBuffer(derivative);
        ComputeKernel kernel = context.compileProgram(program, defines)->createKernel("computeN2Value");
        kernel->addArg(positions);
        kernel->addArg(exclusionBits);
        kernel->addArg(exclusionTiles);
        kernel->addArg(output);
        kernel->addArg(blocks*(blocks+1)/2);
        kernel->addArg(params);
        kernel->addArg(derivative);
        kernel->execute(64,64);
        vector<double> values, derivs;
        context.downloadFixedPointBuffer(output, values);
        context.downloadFixedPointBuffer(derivative, derivs);
        for (int i = 0; i < count; i++) {
            double expected = 0, expectedDeriv = 0;
            for (int j = 0; j < count; j++)
                if (i != j && !(i == 0 && j == 31) && !(i == 31 && j == 0)) {
                    expected += par[i].x+2*par[j].x+0.25*abs(pos[i].x-pos[j].x);
                    expectedDeriv += par[j].y;
                }
            ASSERT_EQUAL_TOL(expected, values[i], 1e-6);
            ASSERT_EQUAL_TOL(expectedDeriv, derivs[i], 1e-6);
        }
    }
}

/** Preserve DPD's per-lane random stream and masked pair order in a full tile. */
void testDPD(int count, bool floating, bool periodic) {
    System system;
    for (int i = 0; i < count; i++) system.addParticle(1);
    MetalContext context(system, nullptr, nullptr, floating);
    const int padded = 32*((count+31)/32), blocks = padded/32;
    vector<mm_int2> fixedTiles;
    vector<int> tiles, neighbors;
    for (int x = 0; x < blocks; x++) {
        fixedTiles.push_back(mm_int2(x,x));
        for (int y = 0; y < x; y++) {
            tiles.push_back(x);
            for (int lane = 0; lane < 32; lane++) neighbors.push_back(32*y+lane);
        }
    }
    const int numTiles = tiles.size();
    if (tiles.empty()) { tiles.push_back(0); neighbors.resize(32); }
    ComputeArray positions, velocity, output, dt, types, params, exclusions, tileArray, interactionCount, centers, sizes, atoms, counter;
    positions.initialize<mm_float4>(context, padded, "dpdPositions");
    velocity.initialize<mm_float4>(context, padded, "dpdVelocities");
    output.initialize<int64_t>(context, padded*3, "dpdVelocityDelta");
    dt.initialize<mm_float2>(context, 1, "dpdDt");
    types.initialize<int>(context, padded, "dpdTypes");
    params.initialize<mm_float2>(context, 1, "dpdParams");
    exclusions.initialize<mm_int2>(context, fixedTiles.size(), "dpdFixedTiles");
    tileArray.initialize<int>(context, tiles.size(), "dpdTiles");
    interactionCount.initialize<unsigned int>(context, 1, "dpdCount");
    centers.initialize<mm_float4>(context, blocks, "dpdCenters");
    sizes.initialize<mm_float4>(context, blocks, "dpdSizes");
    atoms.initialize<int>(context, neighbors.size(), "dpdNeighbors");
    counter.initialize<int>(context, 1, "dpdCounter");
    vector<mm_float4> pos(padded), vel(padded), blockSizes(blocks);
    for (int i = 0; i < padded; i++) {
        pos[i] = mm_float4(0.13f*(i%5), 0.14f*((i/5)%5), 0.15f*(i/25), 0);
        vel[i] = mm_float4(0.01f*(i%3), -0.02f*(i%5), 0.01f*(i%7), 1);
    }
    // Alternating blocks exercise both periodic-copy branches.
    for (int i = 0; i < blocks; i++) blockSizes[i] = mm_float4(i%2 ? 4 : 0.5f, 0.5f, 0.5f, 0);
    positions.upload(pos);
    velocity.upload(vel);
    dt.upload(vector<mm_float2>{mm_float2(0.002f,0.002f)});
    types.upload(vector<int>(padded,0));
    params.upload(vector<mm_float2>{mm_float2(2,2.5f)});
    exclusions.upload(fixedTiles);
    tileArray.upload(tiles);
    interactionCount.upload(vector<unsigned int>{static_cast<unsigned int>(numTiles)});
    centers.upload(vector<mm_float4>(blocks, mm_float4(0,0,0,0)));
    sizes.upload(blockSizes);
    atoms.upload(neighbors);
    map<string,string> defines;
    defines["M_PI"] = "3.14159265358979323846f";
    defines["MAX_CUTOFF"] = "2.5f";
    defines["TILE_SIZE"] = "32";
    defines["WORK_GROUP_SIZE"] = "32";
    if (periodic) defines["USE_PERIODIC"] = "1";
    vector<double> reference;
    for (int mode = 0; mode < 2; mode++) {
        Settings settings;
        settings.dpdParticles = (mode&1) != 0;
        string source = MetalPairwiseOptimizations::apply(CommonKernelSources::dpd, settings);
        if (mode == 0) source = "// unchanged DPD baseline\n"+source;
        context.clearBuffer(output);
        context.clearBuffer(counter);
        ComputeKernel kernel = context.compileProgram(source, defines)->createKernel("integrateDPDPart2");
        kernel->addArg(count);
        kernel->addArg(padded);
        kernel->addArg(positions);
        kernel->addArg(velocity);
        kernel->addArg(output);
        kernel->addArg(dt);
        kernel->addArg(types);
        kernel->addArg(1);
        kernel->addArg(params);
        kernel->addArg(int64_t(77));
        kernel->addArg(1.0f);
        kernel->addArg(mm_float4(10,10,10,0));
        kernel->addArg(mm_float4(0.1f,0.1f,0.1f,0));
        kernel->addArg(mm_float4(10,0,0,0));
        kernel->addArg(mm_float4(0,10,0,0));
        kernel->addArg(mm_float4(0,0,10,0));
        kernel->addArg(exclusions);
        kernel->addArg(int(fixedTiles.size()));
        kernel->addArg(tileArray);
        kernel->addArg(interactionCount);
        kernel->addArg(centers);
        kernel->addArg(sizes);
        kernel->addArg(atoms);
        kernel->addArg(counter);
        // A single SIMD group makes tile assignment and random streams deterministic.
        kernel->execute(32,32);
        vector<double> result;
        context.downloadFixedPointBuffer(output, result);
        if (mode == 0) reference = result;
        else
            for (int i = 0; i < result.size(); i++) ASSERT_EQUAL_TOL(reference[i], result[i], 1e-6);
    }
}

int main(int argc, char** argv) {
    try {
        testSelection();
        if (argc == 2 && string(argv[1]) == "--selection-only") {
            cout << "Metal pairwise template selection tests passed" << endl;
            return 0;
        }
        for (bool floating : {false, true}) {
            for (int count : {1,31,32,33,63,65}) testCustomGBValue(count, floating);
            for (int count : {31,33,65}) {
                testDPD(count, floating, false);
                testDPD(count, floating, true);
            }
        }
    }
    catch (const exception& error) {
        if (string(error.what()).find("No Metal device") != string::npos) return 77;
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Metal pairwise optimization tests passed" << endl;
    return 0;
}
