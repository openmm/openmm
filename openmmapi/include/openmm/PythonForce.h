#ifndef OPENMM_PYTHONFORCE_H_
#define OPENMM_PYTHONFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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

#include "Force.h"
#include "State.h"
#include <map>
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This abstract class represents an interface for performing a computation.  It is not intended to
 * be used or subclassed directly by users.  The Python wrapper contains a subclass that implements
 * the interface using a Python function.
 * @private
 */
class OPENMM_EXPORT PythonForceComputation {
public:
    PythonForceComputation() {
    }
    virtual ~PythonForceComputation() {
    }
    /**
     * Compute forces and energy.  The State contains particle positions, parameters, and
     * optionally periodic box vectors.  Implementations should store the potential energy
     * and particle forces into the energy and forces arguments.
     */
    virtual void compute(const State& state, double& energy, std::vector<Vec3>& forces) const = 0;
};

/**
 * This class provides a mechanism for computing forces and energy with Python code.  To use it,
 * define a Python function that takes a State object as its only argument.  The State contains
 * particle positions and global parameters.  Based on it, the function should compute the
 * potential energy and forces, returning them as its two return values.  The forces should be
 * represented as a NumPy array of shape (# particles, 3).  For example,
 * 
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: python
 * 
 *    def compute(state):
 *        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
 *        k = state.getParameters()['k']
 *        energy = k*np.sum(pos*pos)
 *        force = -0.5*k*pos
 *        return energy*kilojoules_per_mole, force*kilojoules_per_mole/nanometer
 * 
 * \endverbatim
 * 
 * Attaching units to the return values is optional.  If units are omitted, the values are assumed
 * to be in the default units (energy in kJ/mol, forces in kJ/mol/nm).
 * 
 * Now create a Python force, passing the function to the constructor.  If you want the force
 * to depend on global parameters, pass a dict as the second parameter with the names and default
 * values of the parameters.
 * 
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: python
 * 
 *    force = PythonForce(compute, {'k':2.5})
 * 
 * \endverbatim
 * 
 * The default value of a parameter is its value in newly created Contexts.  After a Context is
 * created, you can change the values of parameters by calling setParameter() on it.
 * 
 * The PythonForce cannot tell whether the function you provide makes use of periodic boundary
 * conditions, so you must tell it.  To make the force periodic, call
 * setUsesPeriodicBoundaryConditions(True).  This will cause usesPeriodicBoundaryConditions()
 * to return True, and the State passed to the computation function will contain periodic
 * box vectors.
 * 
 * When using XmlSerializer to save a PythonForce, it uses the Python pickle module to save
 * the computation function.  If it cannot be pickled, you will not be able to serialize the
 * PythonForce.  Functions defined at the top level of a module can usually be pickled, but local
 * functions defined inside another function cannot.
 * 
 * Compared to other types of forces, computing a force with Python code is slow and has high
 * overhead.  When possible, using a different force class is usually preferred.  For example,
 * the Python force shown in the example code above (a harmonic force attracting every particle
 * to the origin) could be implemented just as easily with a CustomExternalForce, and would
 * execute much faster if done that way.
 */
class OPENMM_EXPORT PythonForce : public Force {
public:
    /**
     * Create a PythonForce.  This constructor is used internally, and is not intended for use
     * by users.  The Python wrapper defines an alternate constructor that takes a Python
     * function instead of a PythonForceComputation.
     *
     * @param computation        an object defining how the forces and energy should be computed
     * @param globalParameters   any global parameters used by the force.  Keys are the parameter
     *                           names, and the corresponding values are their default values.
     * @private
     */
    explicit PythonForce(PythonForceComputation* computation, const std::map<std::string, double>& globalParameters);
    ~PythonForce();
    /**
     * Get the PythonForceComputation that defines the computation.
     * @private
     */
    const PythonForceComputation& getComputation() const;
    /**
     * Get all global parameters defined by this force.  Keys are the parameter names, and the
     * corresponding values are their default values.
     */
    const std::map<std::string, double>& getGlobalParameters() const;
    /**
     * Get the pickled representation of the computation function.  If it cannot be pickled,
     * this will be an empty vector.
     */
    const std::vector<char>& getPickledFunction() const;
    /**
     * Set the pickled representation of the computation function.  This is called automatically
     * by the Python constructor.
     * @private
     */
    void setPickledFunction(char* function, int length);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
    /**
     * Set whether or not this force makes use of periodic boundary conditions.
     * If this is set to true, periodic box vectors can be retrieved from the
     * State passed to the computation function.
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
protected:
    ForceImpl* createImpl() const;
private:
    PythonForceComputation* computation;
    std::map<std::string, double> globalParameters;
    bool usePeriodic;
    std::vector<char> pickled;
};

} // namespace OpenMM

#endif /*OPENMM_PYTHONFORCE_H_*/
