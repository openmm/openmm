/* Portions copyright (c) 2014-2015 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "openmm/internal/CompiledExpressionSet.h"

using namespace OpenMM;
using namespace Lepton;
using namespace std;

CompiledExpressionSet::CompiledExpressionSet() {
}

void CompiledExpressionSet::registerExpression(Lepton::CompiledExpression& expression) {
    expressions.push_back(&expression);
    for (int i = 0; i < (int) variables.size(); i++)
        if (expression.getVariables().find(variables[i]) != expression.getVariables().end())
            variableReferences[i].push_back(&expression.getVariableReference(variables[i]));
}

int CompiledExpressionSet::getVariableIndex(const std::string& name) {
    for (int i = 0; i < (int) variables.size(); i++)
        if (variables[i] == name)
            return i;
    int index = variables.size();
    variables.push_back(name);
    variableReferences.push_back(vector<double*>());
    for (auto expression : expressions)
        if (expression->getVariables().find(name) != expression->getVariables().end())
            variableReferences[index].push_back(&expression->getVariableReference(name));
    return index;
}

void CompiledExpressionSet::setVariable(int index, double value) {
    for (auto ref : variableReferences[index])
        *ref = value;
}

int CompiledExpressionSet::getNumVariables() const {
    return variables.size();
}
