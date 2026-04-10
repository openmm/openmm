#ifndef OPENMM_BUSSIRESCALE_H_
#define OPENMM_BUSSIRESCALE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Shared Bussi stochastic velocity rescaling algebra (Bussi et al. 2007,     *
 * J. Chem. Phys. 126, 014101; sign Bussi et al. 2009 Eq. A8).                *
 * Used by Reference and Common kernels and by unit tests with fixed draws. *
 * -------------------------------------------------------------------------- */

#include "SimTKOpenMMRealType.h"
#include <cmath>

namespace OpenMM {
namespace BussiRescale {

/**
 * Result of applying the Bussi rescaling formula for pre-drawn random variates.
 * rGamma must be the chi^2(dof-1) draw used in the kernel (0 when dof==1).
 */
struct Result {
    double alphaSquared;
    double alpha;
    /** Reservoir bookkeeping: K * (1 - alpha^2). */
    double deltaE;
};

/**
 * Compute rescaling factor alpha and related quantities.
 * Preconditions: dof > 0, kineticEnergy > 0, tau > 0 (caller validates).
 */
inline Result computeRescale(double dof, double kineticEnergy, double temperature,
                             double tau, double dt, double R1, double rGamma) {
    double c = std::exp(-dt / tau);
    double targetKE = 0.5 * dof * BOLTZ * temperature;
    double ratio = targetKE / (dof * kineticEnergy);
    double alphaSquared = c + (1.0 - c) * ratio * (rGamma + R1 * R1)
                        + 2.0 * R1 * std::sqrt(c * (1.0 - c) * ratio);
    if (alphaSquared < 0.0)
        alphaSquared = 0.0;
    double alphaMagnitude = std::sqrt(alphaSquared);
    double signTerm = R1 + std::sqrt(c * dof * kineticEnergy / ((1.0 - c) * targetKE));
    double alpha = (signTerm >= 0.0) ? alphaMagnitude : -alphaMagnitude;
    double deltaE = kineticEnergy * (1.0 - alphaSquared);
    return Result{alphaSquared, alpha, deltaE};
}

} // namespace BussiRescale
} // namespace OpenMM

#endif
