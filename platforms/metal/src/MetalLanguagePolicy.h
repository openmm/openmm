/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform code: copyright (c) 2026 Chun-Chi Hung.                       *
 * Author: Chun-Chi Hung.                                                      *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALLANGUAGEPOLICY_H_
#define OPENMM_METALLANGUAGEPOLICY_H_

#include "openmm/OpenMMException.h"
#include <algorithm>
#include <string>
#ifdef __OBJC__
#import <Metal/Metal.h>
#endif

namespace OpenMM {

/** @brief Select an explicit MSL target without confusing SDK, OS, and GPU support. */
class MetalLanguagePolicy {
public:
    /** @return The language-version ceiling, before runtime capability limits are applied.
     *
     * Automatic selection considers every target known to this backend,
     * independently of optimization switches. An explicit language override
     * can lower the ceiling to test shader fallbacks. Each optional shader path
     * still requires both its switch and support in the selected language.
     */
    static int requestedVersion() {
#if OPENMM_METAL_TUNE_LANGUAGE_VERSION
        return OPENMM_METAL_LANGUAGE_VERSION;
#else
        return 410;
#endif
    }

    /** @return Highest listed version allowed by all limits, or zero below 3.0.
     *  @param appleGPU New language features from 3.2 require Apple silicon.
     */
    static int selectVersion(int requested, int sdkMaximum, int osMaximum, bool appleGPU) {
        const int maximum = std::min(std::min(requested, sdkMaximum),
                std::min(osMaximum, appleGPU ? 410 : 310));
        for (int version : {410, 400, 320, 310, 300})
            if (version <= maximum)
                return version;
        return 0;
    }

#ifdef __OBJC__
    /** @return Highest target whose host enum is present in the SDK used to build OpenMM. */
    static int sdkMaximum() {
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 270000
        return 410;
#elif __MAC_OS_X_VERSION_MAX_ALLOWED >= 260000
        return 400;
#elif __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
        return 320;
#elif __MAC_OS_X_VERSION_MAX_ALLOWED >= 140000
        return 310;
#elif __MAC_OS_X_VERSION_MAX_ALLOWED >= 130000
        return 300;
#else
        return 0;
#endif
    }

    /** @return Highest MSL target accepted by the running macOS version. */
    static int runtimeMaximum() {
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 270000
        if (@available(macOS 27.0, *)) return 410;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 260000
        if (@available(macOS 26.0, *)) return 400;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
        if (@available(macOS 15.0, *)) return 320;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 140000
        if (@available(macOS 14.0, *)) return 310;
#endif
        if (@available(macOS 13.0, *)) return 300;
        return 0;
    }

    /** @return SDK enum for a previously capability-checked target; never silently changes it. */
    static MTLLanguageVersion languageVersion(int version) {
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 270000
        if (@available(macOS 27.0, *))
            if (version == 410) return MTLLanguageVersion4_1;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 260000
        if (@available(macOS 26.0, *))
            if (version == 400) return MTLLanguageVersion4_0;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
        if (@available(macOS 15.0, *))
            if (version == 320) return MTLLanguageVersion3_2;
#endif
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 140000
        if (@available(macOS 14.0, *))
            if (version == 310) return MTLLanguageVersion3_1;
#endif
        if (@available(macOS 13.0, *))
            if (version == 300) return MTLLanguageVersion3_0;
        throw OpenMMException("Unsupported Metal language target");
    }
#endif
};

} // namespace OpenMM
#endif
