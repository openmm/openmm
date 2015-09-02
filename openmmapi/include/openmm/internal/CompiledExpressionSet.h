#ifndef OPENMM_COMPILEDEXPRESSIONSET_H_
#define OPENMM_COMPILEDEXPRESSIONSET_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2015 Stanford University and the Authors.      *
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

#include "lepton/CompiledExpression.h"
#include "windowsExport.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class simplifies the management of a set of related CompiledExpressions that share variables.
 */
class OPENMM_EXPORT CompiledExpressionSet {
public:
    CompiledExpressionSet();
    /**
     * Add a CompiledExpression to the set.
     */
    void registerExpression(Lepton::CompiledExpression& expression);
    /**
     * Get the index of a particular variable.
     */
    int getVariableIndex(const std::string& name);
    /**
     * Set the value of a variable on every CompiledExpression.
     * 
     * @param index    the index of the variable, as returned by getVariableIndex()
     * @param value    the value to set it to
     */
    void setVariable(int index, double value);
    /**
     * Get the total number of variables for which indices have been allocated.
     */
    int getNumVariables() const;
private:
    std::vector<Lepton::CompiledExpression*> expressions;
    std::vector<std::string> variables;
    std::vector<std::vector<double*> > variableReferences;
};

} // namespace OpenMM

#endif /*OPENMM_COMPILEDEXPRESSIONSET_H_*/
