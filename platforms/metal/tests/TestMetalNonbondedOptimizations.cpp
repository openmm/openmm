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
#include "openmm/internal/AssertionUtilities.h"
#include <cmath>
#include <iostream>

using namespace OpenMM;
using namespace std;

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

/** @brief Exercise the separate large-block argument layout. */
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
        testSparsePairs(false);
        testSparsePairs(true);
        testSparseResize(false);
        testSparseResize(true);
        testLargeBlocks();
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
