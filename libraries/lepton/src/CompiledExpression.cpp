/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include <utility>

using namespace Lepton;
using namespace std;

CompiledExpression::CompiledExpression() {
}

CompiledExpression::CompiledExpression(const ParsedExpression& expression) {
    ParsedExpression expr = expression.optimize(); // Just in case it wasn't already optimized.
    vector<pair<ExpressionTreeNode, int> > temps;
    compileExpression(expr.getRootNode(), temps);
}

CompiledExpression::~CompiledExpression() {
    for (int i = 0; i < (int) operation.size(); i++)
        if (operation[i] != NULL)
            delete operation[i];
}

CompiledExpression::CompiledExpression(const CompiledExpression& expression) {
    *this = expression;
}

CompiledExpression& CompiledExpression::operator=(const CompiledExpression& expression) {
    arguments = expression.arguments;
    target = expression.target;
    variableIndices = expression.variableIndices;
    variableNames = expression.variableNames;
    workspace.resize(expression.workspace.size());
    argValues.resize(expression.argValues.size());
    operation.resize(expression.operation.size());
    for (int i = 0; i < (int) operation.size(); i++)
        operation[i] = expression.operation[i]->clone();
    return *this;
}

void CompiledExpression::compileExpression(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps) {
    if (findTempIndex(node, temps) != -1)
        return; // We have already processed a node identical to this one.
    
    // Process the child nodes.
    
    vector<int> args;
    for (int i = 0; i < node.getChildren().size(); i++) {
        compileExpression(node.getChildren()[i], temps);
        args.push_back(findTempIndex(node.getChildren()[i], temps));
    }
    
    // Process this node.
    
    if (node.getOperation().getId() == Operation::VARIABLE) {
        variableIndices[node.getOperation().getName()] = (int) workspace.size();
        variableNames.insert(node.getOperation().getName());
    }
    else {
        int stepIndex = (int) arguments.size();
        arguments.push_back(vector<int>());
        target.push_back((int) workspace.size());
        operation.push_back(node.getOperation().clone());
        if (args.size() == 0)
            arguments[stepIndex].push_back(0); // The value won't actually be used.  We just need something there.
        else {
            // If the arguments are sequential, we can just pass a pointer to the first one.
            
            bool sequential = true;
            for (int i = 1; i < args.size(); i++)
                if (args[i] != args[i-1]+1)
                    sequential = false;
            if (sequential)
                arguments[stepIndex].push_back(args[0]);
            else {
                arguments[stepIndex] = args;
                if (args.size() > argValues.size())
                    argValues.resize(args.size(), 0.0);
            }
        }
    }
    temps.push_back(make_pair(node, workspace.size()));
    workspace.push_back(0.0);
}

int CompiledExpression::findTempIndex(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return i;
    return -1;
}

const set<string>& CompiledExpression::getVariables() const {
    return variableNames;
}

double& CompiledExpression::getVariableReference(const string& name) {
    map<string, int>::iterator index = variableIndices.find(name);
    if (index == variableIndices.end())
        throw Exception("getVariableReference: Unknown variable '"+name+"'");
    return workspace[index->second];
}

double CompiledExpression::evaluate() const {
    // Loop over the operations and evaluate each one.
    
    for (int step = 0; step < operation.size(); step++) {
        const vector<int>& args = arguments[step];
        if (args.size() == 1)
            workspace[target[step]] = operation[step]->evaluate(&workspace[args[0]], dummyVariables);
        else {
            for (int i = 0; i < args.size(); i++)
                argValues[i] = workspace[args[i]];
            workspace[target[step]] = operation[step]->evaluate(&argValues[0], dummyVariables);
        }
    }
    return workspace[workspace.size()-1];
}
