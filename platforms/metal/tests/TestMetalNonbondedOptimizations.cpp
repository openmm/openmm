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
#include "openmm/internal/AssertionUtilities.h"
#include <cmath>
#include <iostream>

using namespace OpenMM;
using namespace std;

/** @brief Exercise the separate large-block argument layout. */
void testLargeBlocks() {
    const int atoms = 100001; // The same threshold used by the OpenCL utility.
    System system;
    for (int i = 0; i < atoms; i++) system.addParticle(1);
    MetalContext context(system);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    vector<vector<int> > exclusions;
    nb.addInteraction(true, false, false, 1.0, exclusions, "", 0);
    context.initialize();
    vector<mm_float4> positions(context.getPaddedNumAtoms(), mm_float4(0,0,0,0));
    for (int i = 0; i < atoms; i++) positions[i].x = 2*i;
    context.getPosq().upload(positions);
    nb.prepareInteractions(1);
    vector<unsigned int> counts;
    nb.getInteractionCount().download(counts);
    ASSERT_EQUAL(0, counts[0]);
    ASSERT_EQUAL(1, counts.size());
    vector<mm_float4> centers, bounds;
    nb.getBlockCenters().download(centers);
    nb.getBlockBoundingBoxes().download(bounds);
    ASSERT_EQUAL_TOL(200000, centers.back().x, 1e-6);
    ASSERT_EQUAL_TOL(0, bounds.back().x, 1e-6);
}

int main() {
    try {
        testLargeBlocks();
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
