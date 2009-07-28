/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "lepton/ParsedExpression.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/Operation.h"
#include <limits>
#include <vector>

using namespace Lepton;
using namespace std;

ParsedExpression::ParsedExpression(const ExpressionTreeNode& rootNode) : rootNode(rootNode) {
}

const ExpressionTreeNode& ParsedExpression::getRootNode() const {
    return rootNode;
}

double ParsedExpression::evaluate() const {
    return evaluate(getRootNode(), map<string, double>());
}

double ParsedExpression::evaluate(const std::map<std::string, double>& variables) const {
    return evaluate(getRootNode(), variables);
}

double ParsedExpression::evaluate(const ExpressionTreeNode& node, const map<string, double>& variables) {
    vector<double> args(node.getChildren().size());
    for (int i = 0; i < args.size(); i++)
        args[i] = evaluate(node.getChildren()[i], variables);
    return node.getOperation().evaluate(&args[0], variables);
}

ParsedExpression ParsedExpression::optimize() const {
    ExpressionTreeNode result = precalculateConstantSubexpressions(getRootNode());
    result = substituteSimplerExpression(result);
    return result;
}

ParsedExpression ParsedExpression::optimize(const map<string, double>& variables) const {
    ExpressionTreeNode result = preevaluateVariables(getRootNode(), variables);
    result = precalculateConstantSubexpressions(result);
    result = substituteSimplerExpression(result);
    return ParsedExpression(result);
}

ExpressionTreeNode ParsedExpression::preevaluateVariables(const ExpressionTreeNode& node, const map<string, double>& variables) {
    if (node.getOperation().getId() == Operation::VARIABLE) {
        const Operation::Variable& var = dynamic_cast<const Operation::Variable&>(node.getOperation());
        map<string, double>::const_iterator iter = variables.find(var.getName());
        if (iter == variables.end())
            return node;
        return ExpressionTreeNode(new Operation::Constant(iter->second));
    }
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < children.size(); i++)
        children[i] = preevaluateVariables(node.getChildren()[i], variables);
    return ExpressionTreeNode(node.getOperation().clone(), children);
}

ExpressionTreeNode ParsedExpression::precalculateConstantSubexpressions(const ExpressionTreeNode& node) {
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < children.size(); i++)
        children[i] = precalculateConstantSubexpressions(node.getChildren()[i]);
    ExpressionTreeNode result = ExpressionTreeNode(node.getOperation().clone(), children);
    if (node.getOperation().getId() == Operation::VARIABLE)
        return result;
    for (int i = 0; i < children.size(); i++)
        if (children[i].getOperation().getId() != Operation::CONSTANT)
            return result;
    return ExpressionTreeNode(new Operation::Constant(evaluate(result, map<string, double>())));
}

ExpressionTreeNode ParsedExpression::substituteSimplerExpression(const ExpressionTreeNode& node) {
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < children.size(); i++)
        children[i] = substituteSimplerExpression(node.getChildren()[i]);
    switch (node.getOperation().getId()) {
        case Operation::ADD:
        {
            double first = getConstantValue(children[0]);
            double second = getConstantValue(children[1]);
            if (first == 0.0) // Add 0
                return children[1];
            if (first == 1.0) // Add 1
                return ExpressionTreeNode(new Operation::Increment(), children[1]);
            if (second == 0.0) // Add 0
                return children[0];
            if (second == 1.0) // Add 1
                return ExpressionTreeNode(new Operation::Increment(), children[0]);
            break;
        }
        case Operation::SUBTRACT:
        {
            double first = getConstantValue(children[0]);
            if (first == 0.0) // Subtract from 0
                return ExpressionTreeNode(new Operation::Negate(), children[1]);
            double second = getConstantValue(children[1]);
            if (second == 0.0) // Subtract 0
                return children[0];
            if (second == 1.0) // Subtract 1
                return ExpressionTreeNode(new Operation::Decrement(), children[0]);
            break;
        }
        case Operation::MULTIPLY:
        {
            double first = getConstantValue(children[0]);
            double second = getConstantValue(children[1]);
            if (first == 0.0 || second == 0.0) // Multiply by 0
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (first == 1.0) // Multiply by 1
                return children[1];
            if (second == 1.0) // Multiply by 1
                return children[0];
            if (children[0].getOperation().getId() == Operation::CONSTANT) {
                if (children[1].getOperation().getId() == Operation::MULTIPLY) {
                    if (children[1].getChildren()[0].getOperation().getId() == Operation::CONSTANT) // Combine two multiplies into a single one
                        return ExpressionTreeNode(new Operation::Multiply(), children[1].getChildren()[1], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[1].getChildren()[0])*first)));
                    if (children[1].getChildren()[1].getOperation().getId() == Operation::CONSTANT) // Combine two multiplies into a single one
                        return ExpressionTreeNode(new Operation::Multiply(), children[1].getChildren()[0], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[1].getChildren()[1])*first)));
                }
            }
            if (children[1].getOperation().getId() == Operation::CONSTANT) {
                if (children[0].getOperation().getId() == Operation::MULTIPLY) {
                    if (children[0].getChildren()[0].getOperation().getId() == Operation::CONSTANT) // Combine two multiplies into a single one
                        return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[1], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[0].getChildren()[0])*second)));
                    if (children[0].getChildren()[1].getOperation().getId() == Operation::CONSTANT) // Combine two multiplies into a single one
                        return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[0], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[0].getChildren()[1])*second)));
                }
            }
            break;
        }
        case Operation::DIVIDE:
        {
            double numerator = getConstantValue(children[0]);
            if (numerator = 0.0) // 0 divided by something
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (numerator == 1.0) // 1 divided by something
                return ExpressionTreeNode(new Operation::Reciprocal(), children[1]);
            double denominator = getConstantValue(children[1]);
            if (denominator == 1.0) // Divide by 1
                return children[0];
            if (children[1].getOperation().getId() == Operation::CONSTANT) {
                if (children[0].getOperation().getId() == Operation::MULTIPLY) {
                    if (children[0].getChildren()[0].getOperation().getId() == Operation::CONSTANT) // Combine a multiply and a divide into one multiply
                        return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[1], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[0].getChildren()[0])/denominator)));
                    if (children[0].getChildren()[1].getOperation().getId() == Operation::CONSTANT) // Combine a multiply and a divide into one multiply
                        return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[0], ExpressionTreeNode(new Operation::Constant(getConstantValue(children[0].getChildren()[1])/denominator)));
                }
                return ExpressionTreeNode(new Operation::Multiply(), children[0], ExpressionTreeNode(new Operation::Constant(1.0/denominator))); // Replace a divide with a multiply
            }
            break;
        }
        case Operation::POWER:
        {
            double base = getConstantValue(children[0]);
            if (base == 0.0) // 0 to any power is 0
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (base == 1.0) // 1 to any power is 1
                return ExpressionTreeNode(new Operation::Constant(1.0));
            double exponent = getConstantValue(children[1]);
            if (exponent == 1.0) // x^1 = x
                return children[0];
            if (exponent == -1.0) // x^-1 = recip(x)
                return ExpressionTreeNode(new Operation::Reciprocal(), children[0]);
            if (exponent == 2.0) // x^2 = square(x)
                return ExpressionTreeNode(new Operation::Square(), children[0]);
            if (exponent == 3.0) // x^3 = cube(x)
                return ExpressionTreeNode(new Operation::Cube(), children[0]);
            if (exponent == 0.5) // x^0.5 = sqrt(x)
                return ExpressionTreeNode(new Operation::Sqrt(), children[0]);
            break;
        }
    }
    return ExpressionTreeNode(node.getOperation().clone(), children);
}

ParsedExpression ParsedExpression::differentiate(const std::string& variable) const {
    return differentiate(getRootNode(), variable);
}

ExpressionTreeNode ParsedExpression::differentiate(const ExpressionTreeNode& node, const std::string& variable) {
    vector<ExpressionTreeNode> childDerivs(node.getChildren().size());
    for (int i = 0; i < childDerivs.size(); i++)
        childDerivs[i] = differentiate(node.getChildren()[i], variable);
    return node.getOperation().differentiate(node.getChildren(),childDerivs, variable);
}

double ParsedExpression::getConstantValue(const ExpressionTreeNode& node) {
    if (node.getOperation().getId() == Operation::CONSTANT)
        return dynamic_cast<const Operation::Constant&>(node.getOperation()).getValue();
    return numeric_limits<double>::quiet_NaN();
}

ExpressionProgram ParsedExpression::createProgram() const {
    return ExpressionProgram(*this);
}

ostream& Lepton::operator<<(ostream& out, const ExpressionTreeNode& node) {
    out << node.getOperation().getName();
    if (node.getChildren().size() > 0) {
        out << "(";
        for (int i = 0; i < node.getChildren().size(); i++) {
            if (i > 0)
                out << ", ";
            out << node.getChildren()[i];
        }
        out << ")";
    }
    return out;
}

ostream& Lepton::operator<<(ostream& out, const ParsedExpression& exp) {
    out << exp.getRootNode();
    return out;
}
