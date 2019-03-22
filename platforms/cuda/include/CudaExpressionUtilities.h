#ifndef OPENMM_CUDAEXPRESSIONUTILITIES_H_
#define OPENMM_CUDAEXPRESSIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
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

#include "CudaContext.h"
#include "openmm/TabulatedFunction.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <sstream>
#include <string>
#include <utility>

namespace OpenMM {

/**
 * This class is used by various classes to generate CUDA source code implementing
 * user defined mathematical expressions.
 */

class OPENMM_EXPORT_CUDA CudaExpressionUtilities {
public:
    CudaExpressionUtilities(CudaContext& context);
    /**
     * Generate the source code for calculating a set of expressions.
     *
     * @param expressions    the expressions to generate code for (keys are the variables to store the output values in)
     * @param variables      defines the source code to generate for each variable that may appear in the expressions.  Keys are
     *                       variable names, and the values are the code to generate for them.
     * @param functions      the tabulated functions that may appear in the expressions
     * @param functionNames  defines the variable name for each tabulated function that may appear in the expressions
     * @param prefix         a prefix to put in front of temporary variables
     * @param tempType       the type of value to use for temporary variables (defaults to "real")
     */
    std::string createExpressions(const std::map<std::string, Lepton::ParsedExpression>& expressions, const std::map<std::string, std::string>& variables,
            const std::vector<const TabulatedFunction*>& functions, const std::vector<std::pair<std::string, std::string> >& functionNames,
            const std::string& prefix, const std::string& tempType="real");
    /**
     * Generate the source code for calculating a set of expressions.
     *
     * @param expressions    the expressions to generate code for (keys are the variables to store the output values in)
     * @param variables      defines the source code to generate for each variable or precomputed sub-expression that may appear in the expressions.
     *                       Each entry is an ExpressionTreeNode, and the code to generate wherever an identical node appears.
     * @param functions      the tabulated functions that may appear in the expressions
     * @param functionNames  defines the variable name for each tabulated function that may appear in the expressions
     * @param prefix         a prefix to put in front of temporary variables
     * @param tempType       the type of value to use for temporary variables (defaults to "real")
     */
    std::string createExpressions(const std::map<std::string, Lepton::ParsedExpression>& expressions, const std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& variables,
            const std::vector<const TabulatedFunction*>& functions, const std::vector<std::pair<std::string, std::string> >& functionNames,
            const std::string& prefix, const std::string& tempType="real");
    /**
     * Calculate the spline coefficients for a tabulated function that appears in expressions.
     *
     * @param function   the function for which to compute coefficients
     * @param width      on output, the number of floats used for each value
     * @return the spline coefficients
     */
    std::vector<float> computeFunctionCoefficients(const TabulatedFunction& function, int& width);
    /**
     * Get a Lepton::CustomFunction that can be used to represent a TabulatedFunction when parsing expressions.
     * 
     * @param function   the function for which to get a placeholder
     */
    Lepton::CustomFunction* getFunctionPlaceholder(const TabulatedFunction& function);
    /**
     * Get a Lepton::CustomFunction that can be used to represent the periodicdistance() function when parsing expressions.
     */
    Lepton::CustomFunction* getPeriodicDistancePlaceholder();
private:
    class FunctionPlaceholder : public Lepton::CustomFunction {
        public:
            FunctionPlaceholder(int numArgs) : numArgs(numArgs) {
            }
            int getNumArguments() const {
                return numArgs;
            }
            double evaluate(const double* arguments) const {
                return 0.0;
            }
            double evaluateDerivative(const double* arguments, const int* derivOrder) const {
                return 0.0;
            }
            CustomFunction* clone() const {
                return new FunctionPlaceholder(numArgs);
            }
        private:
            int numArgs;
    };
    void processExpression(std::stringstream& out, const Lepton::ExpressionTreeNode& node,
            std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& temps,
            const std::vector<const TabulatedFunction*>& functions, const std::vector<std::pair<std::string, std::string> >& functionNames,
            const std::string& prefix, const std::vector<std::vector<double> >& functionParams, const std::vector<Lepton::ParsedExpression>& allExpressions, const std::string& tempType);
    std::string getTempName(const Lepton::ExpressionTreeNode& node, const std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& temps);
    void findRelatedCustomFunctions(const Lepton::ExpressionTreeNode& node, const Lepton::ExpressionTreeNode& searchNode,
            std::vector<const Lepton::ExpressionTreeNode*>& nodes);
    void findRelatedPowers(const Lepton::ExpressionTreeNode& node, const Lepton::ExpressionTreeNode& searchNode,
            std::map<int, const Lepton::ExpressionTreeNode*>& powers);
    void callFunction(std::stringstream& out, std::string singleFn, std::string doubleFn, const std::string& arg, const std::string& tempType);
    void callFunction2(std::stringstream& out, std::string singleFn, std::string doubleFn, const std::string& arg1, const std::string& arg2, const std::string& tempType);
    std::vector<std::vector<double> > computeFunctionParameters(const std::vector<const TabulatedFunction*>& functions);
    CudaContext& context;
    FunctionPlaceholder fp1, fp2, fp3, periodicDistance;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAEXPRESSIONUTILITIES_H_*/
