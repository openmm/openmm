#ifndef OPENMM_VECTOR_EXPRESSION_H_
#define OPENMM_VECTOR_EXPRESSION_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2018 Stanford University and the Authors.           *
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

#include "openmm/Vec3.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/ParsedExpression.h"
#include "windowsExport.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class wraps a Lepton expression to support parsing and evaluating expressions
 * that operate on Vec3s instead of scalars.  In addition to the standard functions
 * supported by Lepton, vector expressions can use the functions dot(), cross(),
 * vector(), _x(), _y(), and _z().
 */
class OPENMM_EXPORT VectorExpression {
public:
    /**
     * Create a VectorExpression by parsing a mathematical expression.
     *
     * @param expression        the mathematical expression to parse
     * @param customFunctions   a map specifying user defined functions that may appear in the expression.
     *                          The keys are function names, and the values are corresponding CustomFunction objects.
     */
    VectorExpression(const std::string& expression, const std::map<std::string, Lepton::CustomFunction*>& customFunctions);
    /**
     * Create a VectorExpression for evaluating a ParsedExpression.
     *
     * @param expression        the mathematical expression to evaluate
     */
    VectorExpression(const Lepton::ParsedExpression& expression);
    VectorExpression(const VectorExpression& expression);
    ~VectorExpression();
    /**
     * Evaluate the expression.
     *
     * @param variables    a map specifying the values of all variables that appear in the expression.  If any
     *                     variable appears in the expression but is not included in this map, an exception
     *                     will be thrown.
     */
    Vec3 evaluate(const std::map<std::string, Vec3>& variables) const;
    VectorExpression& operator=(const VectorExpression& exp);
private:
    class VectorOperation;
    class OpWrapper1;
    class OpWrapper2;
    class OpWrapper3;
    class Variable;
    class Constant;
    class Dot;
    class Cross;
    class Component;
    class Vector;
    Lepton::ParsedExpression parsed;
    Lepton::ExpressionProgram program;
    mutable std::vector<Vec3> stack;
    std::vector<VectorOperation*> operations;
    void analyzeExpression(const Lepton::ParsedExpression& expression);
};

} // namespace OpenMM

#endif /*OPENMM_VECTOR_EXPRESSION_H_*/
