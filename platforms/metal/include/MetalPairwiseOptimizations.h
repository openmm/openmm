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

#ifndef OPENMM_METALPAIRWISEOPTIMIZATIONS_H_
#define OPENMM_METALPAIRWISEOPTIMIZATIONS_H_

#include <set>
#include <string>

namespace OpenMM {

/**
 * @brief Optional transformations for audited Common pairwise templates.
 *
 * Register paths change lane communication, not fixed-point representation,
 * exclusions, or pair enumeration. GBSA selects readable shader-side paths
 * with defines; other pairwise kernels retain audited template adaptation.
 * Register paths require 32 active lanes per SIMD group, including padded
 * lanes of the last tile. A complete identifier-template match identifies
 * read-only GBSA Born-force parameters without changing the expression.
 */
class MetalPairwiseOptimizations {
public:
    /** Independently selectable transformations; defaults are the Common path. */
    struct Settings {
        bool customGBValue = false;
        bool customGBEnergy = false;
        bool gbsaBorn = false;
        bool gbsaForce = false;
        bool dpdParticles = false;
        bool dpdTile = false;
        bool customHbond = false;
    };

    /** @return Settings selected by the independent build-time switches. */
    static Settings getBuildSettings();
    /** @brief Apply build-time settings before Metal's language/ABI adaptation. */
    static std::string apply(const std::string& source);
    /** @brief Explicit settings support focused tests of individual paths. */
    static std::string apply(const std::string& source, const Settings& settings);
    /**
     * @brief Identify read-only Born-force parameters in the original GBSA snippet.
     *
     * Only a complete Common gbsaObc2 template with consistent identifier-only
     * substitutions is recognized. Call before cutoff replacement. Unknown or modified snippets return an empty set.
     * The returned names include the particle suffix (for example, bornForce1).
     */
    static std::set<std::string> getGBSAChainRuleBornForceParameters(const std::string& source);

};

} // namespace OpenMM
#endif
