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

#include "openmm/MultiModeCavityDisplacer.h"
#include "openmm/OpenMMException.h"
#include <cmath>

using namespace OpenMM;

MultiModeCavityDisplacer::MultiModeCavityDisplacer(int numModes, double omega1, double lambda1,
                                                   double cavityLength, double moleculeZ,
                                                   double photonMass) :
        numModes(numModes), omega1(omega1), lambda1(lambda1),
        cavityLength(cavityLength), moleculeZ(moleculeZ), photonMass(photonMass) {
    if (numModes < 1)
        throw OpenMMException("MultiModeCavityDisplacer: numModes must be >= 1");
    if (omega1 <= 0)
        throw OpenMMException("MultiModeCavityDisplacer: omega1 must be positive");
    if (photonMass <= 0)
        throw OpenMMException("MultiModeCavityDisplacer: photonMass must be positive");
    if (cavityLength <= 0)
        throw OpenMMException("MultiModeCavityDisplacer: cavityLength must be positive");
    
    const double AMU_TO_AU = 1822.888;
    double photonMass_au = photonMass * AMU_TO_AU;
    
    spatialProfiles.resize(numModes);
    displacementFactors.resize(numModes);
    for (int i = 0; i < numModes; i++) {
        int n = i + 1;
        double omega_n = n * omega1;
        double lambda_n = lambda1;
        
        // Spatial profile: f_n = sin(n * pi * z0 / L)
        spatialProfiles[i] = std::sin(n * M_PI * moleculeZ / cavityLength);
        
        // Displacement factor: -(lambda_n / (photonMass_au * omega_n))
        // This has units that, when multiplied by dipole (e*nm), give position (nm)
        // via the unit conversion chain in the kernel
        displacementFactors[i] = -lambda_n / (photonMass_au * omega_n);
    }
}

double MultiModeCavityDisplacer::getDisplacementFactor(int modeIndex) const {
    if (modeIndex < 0 || modeIndex >= numModes)
        throw OpenMMException("MultiModeCavityDisplacer: mode index out of range");
    return displacementFactors[modeIndex];
}

double MultiModeCavityDisplacer::getCombinedFactor(int modeIndex) const {
    if (modeIndex < 0 || modeIndex >= numModes)
        throw OpenMMException("MultiModeCavityDisplacer: mode index out of range");
    return displacementFactors[modeIndex] * spatialProfiles[modeIndex];
}
