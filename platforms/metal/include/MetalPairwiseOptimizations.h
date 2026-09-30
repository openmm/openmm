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

#include <string>

namespace OpenMM {

/**
 * @brief Optional register communication for audited Common pairwise templates.
 *
 * These transformations change lane communication, not the pair expressions,
 * fixed-point representation, exclusions, or pair enumeration. Each SIMD group
 * must contain 32 active lanes, including the padded lanes of the last tile.
 * A full template-skeleton match limits generated-expression substitutions to
 * the holes already provided by Common. Unrecognized source is left unchanged.
 */
class MetalPairwiseOptimizations {
public:
    /** Independently selectable transformations; defaults are the Common path. */
    struct Settings {
        bool customGBValue = false;
        bool customGBEnergy = false;
        bool dpdParticles = false;
        bool dpdTile = false;
    };

    /** @return Settings selected by the independent build-time switches. */
    static Settings getBuildSettings();
    /** @brief Apply build-time settings before Metal's language/ABI adaptation. */
    static std::string apply(const std::string& source);
    /** @brief Explicit settings support focused tests of individual paths. */
    static std::string apply(const std::string& source, const Settings& settings);
};

} // namespace OpenMM
#endif
