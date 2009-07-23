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

#include "ParsedExpression.h"
#include "Operation.h"
#include <vector>

using namespace Lepton;
using namespace std;

ParsedExpression::ParsedExpression(ExpressionTreeNode rootNode) : rootNode(rootNode) {
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
    ParsedExpression result = precalculateConstantSubexpressions(getRootNode());
    result = substituteSimplerExpression(result.getRootNode());
    return result;
}

ParsedExpression ParsedExpression::optimize(const map<string, double>& variables) const {
    ParsedExpression result = preevaluateVariables(getRootNode(), variables);
    result = precalculateConstantSubexpressions(result.getRootNode());
    result = substituteSimplerExpression(result.getRootNode());
    return result;
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
        case Operation::DIVIDE:
            if (children[0].getOperation().getId() == Operation::CONSTANT) {
                if (dynamic_cast<const Operation::Constant&>(children[0].getOperation()).getValue() == 1.0)
                    return ExpressionTreeNode(new Operation::Reciprocal(), children[1]);
            }
            break;
        case Operation::POWER:
            if (children[1].getOperation().getId() == Operation::CONSTANT) {
                double exponent = dynamic_cast<const Operation::Constant&>(children[1].getOperation()).getValue();
                if (exponent == 1.0)
                    return children[0];
                if (exponent == -1.0)
                    return ExpressionTreeNode(new Operation::Reciprocal(), children[0]);
                if (exponent == 2.0)
                    return ExpressionTreeNode(new Operation::Square(), children[0]);
                if (exponent == 3.0)
                    return ExpressionTreeNode(new Operation::Cube(), children[0]);
                if (exponent == 0.5)
                    return ExpressionTreeNode(new Operation::Sqrt(), children[0]);
            }
            break;
    }
    return ExpressionTreeNode(node.getOperation().clone(), children);
}
