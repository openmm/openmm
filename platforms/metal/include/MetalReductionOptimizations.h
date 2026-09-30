/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.*
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALREDUCTIONOPTIMIZATIONS_H_
#define OPENMM_METALREDUCTIONOPTIMIZATIONS_H_

#include <string>

namespace OpenMM {

/** @brief Select narrow SIMD alternatives without changing Common force mathematics. */
class MetalReductionOptimizations {
public:
    /** @brief Independent choices; an explicit Settings starts with every choice disabled. */
    struct Settings {
        bool centroid = false;
        bool rg = false;
        bool rmsd = false;
        bool orientation = false;
        bool lcpoScan = false;
        bool manyParticleScan = false;
        bool sortBucketScan = false;
        bool sortRegisterBitonic = false;
    };
    /** @return The independent build-time choices, before a caller's type-specific eligibility checks. */
    static Settings getBuildSettings();
    /** @brief Apply the independent build-time switches to recognized template functions. */
    static std::string apply(const std::string& source);
    /**
     * @brief Apply explicit choices, also used by tests to compare ON and OFF in one process.
     *
     * A replacement requires a byte-for-byte match to the original function.
     * Unrelated helpers and already-transformed functions remain unchanged.
     * Call before Common-to-MSL adaptation and before replacing sort type tokens.
     */
    static std::string apply(const std::string& source, const Settings& settings);
};

} // namespace OpenMM
#endif
