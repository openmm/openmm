#ifndef OPENMM_MULTIMODECAVITYDISPLACER_H_
#define OPENMM_MULTIMODECAVITYDISPLACER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Muhammad Hasyim                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "internal/windowsExport.h"
#include <vector>
#include <cmath>

namespace OpenMM {

/**
 * Utility class for displacing multi-mode cavity particles to their equilibrium
 * positions. When the cavity coupling is switched on, each mode's photon particle
 * needs to be set to its equilibrium to avoid sudden energy injection.
 *
 * For mode n (1-indexed), the equilibrium position is:
 *   q_n_eq_xy = -(eps_n * f_n / K_n) * d_xy
 *   q_n_eq_z  = 0  (z coordinate is preserved)
 *
 * where:
 *   eps_n = lambda_n * omega_n  (effective coupling)
 *   K_n = photonMass * omega_n^2  (spring constant)
 *   f_n = sin(n * pi * z0 / L)  (spatial profile)
 *   d_xy = molecular dipole moment (x,y components)
 *
 * This simplifies to:
 *   q_n_eq_xy = -(lambda_n / (photonMass_au * omega_n)) * f_n * d_xy
 *             = -(lambda1 / (photonMass_au * n * omega1)) * f_n * d_xy
 *
 * Usage from Python:
 *   displacer = MultiModeCavityDisplacer(numModes, omega1, lambda1, cavityLength, moleculeZ, photonMass)
 *   equilibriumPositions = displacer.computeEquilibriumPositions(charges, positions, cavityIndices)
 *   # Then update context positions for each cavity particle
 */
class OPENMM_EXPORT MultiModeCavityDisplacer {
public:
    /**
     * Create a MultiModeCavityDisplacer.
     *
     * @param numModes       number of cavity modes
     * @param omega1         fundamental cavity frequency in atomic units
     * @param lambda1        dimensionless coupling for fundamental mode
     * @param cavityLength   cavity length in nm
     * @param moleculeZ      molecule z-position in nm
     * @param photonMass     photon mass in amu (default: 1/1822.888)
     */
    MultiModeCavityDisplacer(int numModes, double omega1, double lambda1,
                             double cavityLength, double moleculeZ,
                             double photonMass = 1.0/1822.888);
    /**
     * Compute the combined displacement factor for a specific mode.
     * This is getDisplacementFactor(modeIndex) * getSpatialProfile(modeIndex).
     * Multiply by dipoleX/Y to get the equilibrium x/y position.
     *
     * @param modeIndex  0-based mode index
     * @return the combined factor (units: 1/e)
     */
    double getCombinedFactor(int modeIndex) const;
    /**
     * Get the displacement factor for mode n: -(lambda_n / (photonMass_au * omega_n))
     * This is the factor that multiplies f_n * d_xy to give the equilibrium position.
     *
     * @param modeIndex  0-based mode index
     * @return the displacement factor (units: 1/e, converts d in e*nm to q_eq in nm)
     */
    double getDisplacementFactor(int modeIndex) const;
    /**
     * Get the spatial profile for a mode.
     */
    double getSpatialProfile(int modeIndex) const {
        return spatialProfiles[modeIndex];
    }
    int getNumModes() const { return numModes; }
private:
    int numModes;
    double omega1;
    double lambda1;
    double photonMass;
    double cavityLength;
    double moleculeZ;
    std::vector<double> spatialProfiles;
    std::vector<double> displacementFactors;
};

} // namespace OpenMM

#endif /*OPENMM_MULTIMODECAVITYDISPLACER_H_*/
