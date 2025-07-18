#ifndef OPENMM_QTBINTEGRATORUTILITIES_H_
#define OPENMM_QTBINTEGRATORUTILITIES_H_

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

#include "openmm/internal/ThreadPool.h"
#include <vector>

namespace OpenMM {

/**
 * This class defines a set of utility functions that are useful in implementing QTBIntegrator.
 */

class OPENMM_EXPORT QTBIntegratorUtilities {
public:
    /**
     * Calculate the target noise spectrum.
     * 
     * @param temperature  the target temperature, in Kelvin
     * @param friction     the friction coefficient in 1/ps
     * @param dt           the integration step size, in ps
     * @param numFreq      the number of frequency bins in the spectrum
     * @param theta        the standard version of the spectrum is stored into this
     * @param thetad       the deconvolved version of the spectrum is stored into this
     * @param threads      a ThreadPool to use for parallelization
     */
    static void calculateSpectrum(double temperature, double friction, double dt, int numFreq, std::vector<double>& theta, std::vector<double>& thetad, ThreadPool& threads);
};

} // namespace OpenMM

#endif /*OPENMM_QTBINTEGRATORUTILITIES_H_*/
