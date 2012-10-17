#ifndef OPENMM_OPENCLEXPRESSIONUTILITIES_H_
#define OPENMM_OPENCLEXPRESSIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2011 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "OpenCLContext.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <sstream>
#include <string>
#include <utility>

namespace OpenMM {

/**
 * This class is used by various classes to generate OpenCL source code implementing
 * user defined mathematical expressions.
 */

class OPENMM_EXPORT OpenCLExpressionUtilities {
public:
    OpenCLExpressionUtilities(OpenCLContext& context) : context(context) {
    }
    /**
     * Generate the source code for calculating a set of expressions.
     *
     * @param expressions    the expressions to generate code for (keys are the variables to store the output values in)
     * @param variables      defines the source code to generate for each variable that may appear in the expressions.  Keys are
     *                       variable names, and the values are the code to generate for them.
     * @param functions      defines the variable name for each tabulated function that may appear in the expressions
     * @param prefix         a prefix to put in front of temporary variables
     * @param functionParams the variable name containing the parameters for each tabulated function
     * @param tempType       the type of value to use for temporary variables (defaults to "real")
     */
    std::string createExpressions(const std::map<std::string, Lepton::ParsedExpression>& expressions, const std::map<std::string, std::string>& variables,
            const std::vector<std::pair<std::string, std::string> >& functions, const std::string& prefix, const std::string& functionParams, const std::string& tempType="real");
    /**
     * Generate the source code for calculating a set of expressions.
     *
     * @param expressions    the expressions to generate code for (keys are the variables to store the output values in)
     * @param variables      defines the source code to generate for each variable or precomputed sub-expression that may appear in the expressions.
     *                       Each entry is an ExpressionTreeNode, and the code to generate wherever an identical node appears.
     * @param functions      defines the variable name for each tabulated function that may appear in the expressions
     * @param prefix         a prefix to put in front of temporary variables
     * @param functionParams the variable name containing the parameters for each tabulated function
     * @param tempType       the type of value to use for temporary variables (defaults to "float")
     */
    std::string createExpressions(const std::map<std::string, Lepton::ParsedExpression>& expressions, const std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& variables,
            const std::vector<std::pair<std::string, std::string> >& functions, const std::string& prefix, const std::string& functionParams, const std::string& tempType="float");
    /**
     * Calculate the spline coefficients for a tabulated function that appears in expressions.
     *
     * @param values         the tabulated values of the function
     * @param min            the value of the independent variable corresponding to the first element of values
     * @param max            the value of the independent variable corresponding to the last element of values
     * @return the spline coefficients
     */
    std::vector<mm_float4> computeFunctionCoefficients(const std::vector<double>& values, double min, double max);
    class FunctionPlaceholder;
private:
    void processExpression(std::stringstream& out, const Lepton::ExpressionTreeNode& node,
            std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& temps,
            const std::vector<std::pair<std::string, std::string> >& functions, const std::string& prefix, const std::string& functionParams,
            const std::vector<Lepton::ParsedExpression>& allExpressions, const std::string& tempType);
    std::string getTempName(const Lepton::ExpressionTreeNode& node, const std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& temps);
    void findRelatedTabulatedFunctions(const Lepton::ExpressionTreeNode& node, const Lepton::ExpressionTreeNode& searchNode,
            const Lepton::ExpressionTreeNode*& valueNode, const Lepton::ExpressionTreeNode*& derivNode);
    void findRelatedPowers(const Lepton::ExpressionTreeNode& node, const Lepton::ExpressionTreeNode& searchNode,
            std::map<int, const Lepton::ExpressionTreeNode*>& powers);
    OpenCLContext& context;
};

/**
 * This class serves as a placeholder for custom functions in expressions.
 */

class OpenCLExpressionUtilities::FunctionPlaceholder : public Lepton::CustomFunction {
public:
    int getNumArguments() const {
        return 1;
    }
    double evaluate(const double* arguments) const {
        return 0.0;
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return 0.0;
    }
    CustomFunction* clone() const {
        return new FunctionPlaceholder();
    }
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLEXPRESSIONUTILITIES_H_*/
