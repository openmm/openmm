/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Adapted from Common pairwise kernels used by CUDA, HIP, and OpenCL.        *
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

#include "MetalPairwiseOptimizations.h"
#include "CommonKernelSources.h"
#include "openmm/OpenMMException.h"
#include <cctype>
#include <regex>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace std;

MetalPairwiseOptimizations::Settings MetalPairwiseOptimizations::getBuildSettings() {
    Settings settings;
    return settings;
}

string MetalPairwiseOptimizations::apply(const string& source) {
    return apply(source, getBuildSettings());
}

string MetalPairwiseOptimizations::apply(const string& source, const Settings& settings) {
    return source;
}
