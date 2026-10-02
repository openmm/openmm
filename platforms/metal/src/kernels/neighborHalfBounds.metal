/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform code: copyright (c) 2026 Chun-Chi Hung.                       *
 * Author: Chun-Chi Hung.                                                      *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

/**
 * @brief Round nonnegative bounds upward so FP16 boxes never shrink.
 * Only a rounded-down value advances; exact values and infinity stay unchanged.
 * Older language targets retain the equivalent positive-half bit increment.
 */
DEVICE half metalBoundHalf(real value) {
    half result = half(value);
#if OPENMM_METAL_FAST_FP16_BOUNDS_NEXTAFTER && __METAL_VERSION__ >= 310
    if (float(result) < value) result = nextafter(result, half(INFINITY));
#else
    ushort bits = as_type<ushort>(result);
    if (float(result) < value) result = as_type<half>(ushort(bits+1));
#endif
    return result;
}

/** @brief Compress the three bounding-box half widths; the unused lane is zero. */
DEVICE half4 metalBoundsHalf(real4 value) {
    return half4(metalBoundHalf(value.x), metalBoundHalf(value.y), metalBoundHalf(value.z), half(0));
}
