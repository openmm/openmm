/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform code: copyright (c) 2026 Chun-Chi Hung.                       *
 * Author: Chun-Chi Hung.                                                      *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

/**
 * @brief Select arithmetic assumptions before any shader functions are defined.
 *
 * The host gates pragma support by the running compiler's OS availability and
 * disables fast math for floating accumulators and explicitly safe programs.
 * Older runtimes use the equivalent legacy compilation option instead.
 * FP32 library-function selection is a separate host compiler option; it is
 * not changed by this arithmetic-mode pragma.
 */
#if OPENMM_METAL_HAS_MATH_PRAGMAS
#if OPENMM_METAL_USE_FAST_MATH
#pragma METAL fp math_mode(fast)
#else
#pragma METAL fp math_mode(safe)
#endif
#endif
