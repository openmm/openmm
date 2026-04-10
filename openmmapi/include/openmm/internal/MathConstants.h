#ifndef OPENMM_MATHCONSTANTS_H_
#define OPENMM_MATHCONSTANTS_H_

namespace OpenMM {

/**
 * Portable pi as double. MSVC's standard library does not define M_PI in
 * <cmath> unless _USE_MATH_DEFINES is defined before the include; public API
 * headers must not rely on M_PI. Matches the usual M_PI value where provided.
 */
constexpr double OpenMM_Pi = 3.14159265358979323846264338327950288;

} // namespace OpenMM

#endif
