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
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * A MinimizationReporter can be passed to LocalEnergyMinimizer::minimize() to provide
 * periodic information on the progress of minimization, and to give you the chance to
 * stop minimization early.  Define a subclass that overrides report() and implement it
 * to take whatever action you want.
 * 
 * To correctly interpret the information passed to the reporter, you need to know a bit
 * about how the minimizer works.  The L-BFGS algorithm used by the minimizer does not
 * support constraints.  The minimizer therefore replaces all constraints with harmonic
 * restraints, then performs unconstrained minimization of a combined objective function
 * that is the sum of the system's potential energy and the restraint energy.  Once
 * minimization completes, it checks whether all constraints are satisfied to an acceptable
 * tolerance.  It not, it increases the strength of the harmonic restraints and performs
 * additional minimization.  If the error in constrained distances is especially large,
 * it may choose to throw out all work that has been done so far and start over with
 * stronger restraints.  This has several important consequences.
 * 
 * <ul>
 * <li>The objective function being minimized not actually the same as the potential energy.</li>
 * <li>The objective function and the potential energy can both increase between iterations.</li>
 * <li>The total number of iterations performed could be larger than the number specified
 *     by the maxIterations argument, if that many iterations leaves unacceptable constraint errors.</li>
 * <li>All work is provisional.  It is possible for the minimizer to throw it out and start over.</li>
 * </ul>
 */
class OPENMM_EXPORT MinimizationReporter {
public:
    MinimizationReporter() {
    }
    virtual ~MinimizationReporter() {
    }
    /**
     * This is called after each iteration to provide information about the current status
     * of minimization.  It receives the current particle coordinates, the gradient of the
     * objective function with respect to them, and a set of useful statistics.  In particular,
     * args contains these values:
     * 
     * "system energy": the current potential energy of the system
     * 
     * "restraint energy": the energy of the harmonic restraints
     * 
     * "restraint strength": the force constant of the restraints (in kJ/mol/nm^2)
     * 
     * "max constraint error": the maximum relative error in the length of any constraint
     * 
     * If this function returns true, it will cause the L-BFGS optimizer to immediately
     * exit.  If all constrained distances are sufficiently close to their target values,
     * minimize() will return.  If any constraint error is unacceptably large, it will instead
     * cause the minimizer to immediately increase the strength of the harmonic restraints and
     * perform additional optimization.
     * 
     * @param iteration  the index of the current iteration.  This refers to the current call
     *                   to the L-BFGS optimizer.  Each time the minimizer increases the restraint
     *                   strength, the iteration index is reset to 0.
     * @param x          the current particle positions in flattened order: the three coordinates
     *                   of the first particle, then the three coordinates of the second particle, etc.
     * @param grad       the current gradient of the objective function (potential energy plus
     *                   restraint energy) with respect to the particle coordinates, in flattened
     *                   order
     * @param args       additional statistics described above about the current state of minimization
     * @return whether to immediately stop minimization
     */
    virtual bool report(int iteration, const std::vector<double>& x, const std::vector<double>& grad, std::map<std::string, double>& args) {
        return false;
    }
};

/**
 * Given a Context, this class searches for a new set of particle positions that represent
 * a local minimum of the potential energy.  The search is performed with the L-BFGS algorithm.
 * Distance constraints are enforced during minimization by adding a harmonic restraining
 * force to the potential function.  The strength of the restraining force is steadily increased
 * until the minimum energy configuration satisfies all constraints to within the tolerance
 * specified by the Context's Integrator.
 * 
 * Energy minimization is done using the force groups defined by the Integrator.
 * If you have called setIntegrationForceGroups() on it to restrict the set of forces
 * used for integration, only the energy of the included forces will be minimized.
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
     *                       this tolerance (in kJ/mol/nm).  The default value is 10.
     * @param maxIterations  the maximum number of iterations to perform.  If this is 0, minimation is continued
     *                       until the results converge without regard to how many iterations it takes.  The
     *                       default value is 0.
     * @param reporter       an optional MinimizationReporter to invoke after each iteration.  This can be used
     *                       to monitor the progress of minimization or to stop minimization early.
     */
    static void minimize(Context& context, double tolerance = 10, int maxIterations = 0, MinimizationReporter* reporter = NULL);
};

} // namespace OpenMM

#endif /*OPENMM_LOCALENERGYMINIMIZER_H_*/
