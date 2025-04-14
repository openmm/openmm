#ifndef OPENMM_DPDINTEGRATORUTILITIES_H_
#define OPENMM_DPDINTEGRATORUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/DPDIntegrator.h"
#include <vector>

namespace OpenMM {

/**
 * This class defines a set of utility functions that are useful in implementing DPDIntegrator.
 */

class OPENMM_EXPORT DPDIntegratorUtilities {
public:
    /**
     * Create tables of the parameters to use for each pair of types.
     * 
     * @param integrator         the integrator to analyze
     * @param numParticles       the number of particles in the System
     * @param[out] numTypes      the number of unique particle types
     * @param[out] particleType  the type of each particle, represented as an index between 0 and numTypes
     * @param[out] friction      element [i][j] is the friction to use between particles of type i and j
     * @param[out] cutoff        element [i][j] is the cutoff to use between particles of type i and j
     * @param[out] maxCutoff     the maximum cutoff distance between any pair of particles
     */
    static void createTypeTables(const DPDIntegrator& integrator, int numParticles, int& numTypes, std::vector<int>& particleType,
            std::vector<std::vector<double> >& friction, std::vector<std::vector<double> >& cutoff, double& maxCutoff);
};

} // namespace OpenMM

#endif /*OPENMM_DPDINTEGRATORUTILITIES_H_*/
