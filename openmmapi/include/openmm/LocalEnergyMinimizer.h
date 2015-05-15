#ifndef OPENMM_LOCALENERGYMINIMIZER_H_
#define OPENMM_LOCALENERGYMINIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "Context.h"

namespace OpenMM {

/**
 * Given a Context, this class searches for a new set of particle positions that represent
 * a local minimum of the potential energy.  The search is performed with the L-BFGS algorithm.
 * Distance constraints are enforced during minimization by adding a harmonic restraining
 * force to the potential function.  The strength of the restraining force is steadily increased
 * until the minimum energy configuration satisfies all constraints to within the tolerance
 * specified by the Context's Integrator.
 */

class OPENMM_EXPORT LocalEnergyMinimizer {
public:
    /**
     * Search for a new set of particle positions that represent a local potential energy minimum.
     * On exit, the Context will have been updated with the new positions.
     *
     * @param context        a Context specifying the System to minimize and the initial particle positions
     * @param tolerance      this specifies how precisely the energy minimum must be located.  Minimization
     *                       will be halted once the root-mean-square value of all force components reaches
     *                       this tolerance.  The default value is 10.
     * @param maxIterations  the maximum number of iterations to perform.  If this is 0, minimation is continued
     *                       until the results converge without regard to how many iterations it takes.  The
     *                       default value is 0.
     */
    static void minimize(Context& context, double tolerance = 10, int maxIterations = 0);
};

} // namespace OpenMM

#endif /*OPENMM_LOCALENERGYMINIMIZER_H_*/
