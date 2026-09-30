/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
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

#ifndef OPENMM_METALSOURCEADAPTER_H_
#define OPENMM_METALSOURCEADAPTER_H_

#include <string>

namespace OpenMM {

/**
 * @brief Adapts the Common Compute kernel dialect to Metal 3.0.
 *
 * This is a deliberately small source adapter, not a general OpenCL compiler.
 * It preserves preprocessor branches and mathematical kernel bodies. Metal's
 * own preprocessor selects branches and assigns dense argument-buffer IDs.
 * Entry-point arguments and explicit execution built-ins use the Metal ABI;
 * the host program supplies the matching argument encoder. Scoped minimization
 * additionally adapts audited fixed-point accumulator types and conversions to
 * float, preserving ordinary 64-bit integer indices and random-number state.
 * Independent build switches can select existing CUDA branches within exact
 * Common template functions, without changing other kernels' platform defines.
 */
class MetalSourceAdapter {
public:
    /** @return Whether source contains a Common/OpenCL entry-point declaration. */
    static bool isCommonSource(const std::string& source);
    /**
     * @brief Translate Common/OpenCL source without invoking an external compiler.
     * @param source Source after Common's textual substitutions, before preprocessing.
     * @param floatingAccumulators Use the scoped minimization accumulator ABI.
     * @return MSL source to follow the Metal Common preamble and compile definitions.
     * @throws OpenMMException If a signature is outside the supported kernel dialect.
     */
    static std::string translate(const std::string& source, bool floatingAccumulators=false);
};

} // namespace OpenMM
#endif
