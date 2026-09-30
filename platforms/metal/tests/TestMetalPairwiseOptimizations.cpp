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
}

int main(int argc, char** argv) {
    try {
        testSelection();
        if (argc == 2 && string(argv[1]) == "--selection-only") {
            cout << "Metal pairwise template selection tests passed" << endl;
            return 0;
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
