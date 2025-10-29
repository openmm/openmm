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

class OPENMM_EXPORT PythonForceComputation {
public:
    PythonForceComputation() {
    }
    virtual ~PythonForceComputation() {
    }
    /**
     */
    virtual void compute(const State& state, double& energy, std::vector<Vec3>& forces) const = 0;
};

/**
 */

class OPENMM_EXPORT PythonForce : public Force {
public:
    /**
     * Create a PythonForce.
     *
     * @param computation
     * @param globalParameters
     */
    explicit PythonForce(PythonForceComputation* computation, const std::map<std::string, double>& globalParameters);
    ~PythonForce();
    const PythonForceComputation& getComputation() const;
    const std::map<std::string, double>& getGlobalParameters() const;
    const std::vector<char>& getPickledFunction() const;
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
    bool setUsesPeriodicBoundaryConditions(bool periodic);
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
