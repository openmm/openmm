/* Portions copyright (c) 2018 Stanford University and Simbios.
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

#include "openmm/internal/VectorExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace Lepton;
using namespace std;

static map<string, double> noVariables;

/**
 * Similar to Lepton::Operation, but it operates on Vec3s instead of doubles.
 */
class VectorExpression::VectorOperation {
public:
    virtual ~VectorOperation() {
    }
    virtual Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const = 0;
};

class VectorExpression::OpWrapper1 : public VectorExpression::VectorOperation {
public:
    OpWrapper1(const Operation& op) : op(op) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        Vec3 result;
        double arg = args[0][0];
        result[0] = op.evaluate(&arg, noVariables);
        arg = args[0][1];
        result[1] = op.evaluate(&arg, noVariables);
        arg = args[0][2];
        result[2] = op.evaluate(&arg, noVariables);
        return result;
    }
private:
    const Operation& op;
};

class VectorExpression::OpWrapper2 : public VectorExpression::VectorOperation {
public:
    OpWrapper2(const Operation& op) : op(op) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        Vec3 result;
        double args0[] = {args[0][0], args[1][0]};
        result[0] = op.evaluate(args0, noVariables);
        double args1[] = {args[0][1], args[1][1]};
        result[1] = op.evaluate(args1, noVariables);
        double args2[] = {args[0][2], args[1][2]};
        result[2] = op.evaluate(args2, noVariables);
        return result;
    }
private:
    const Operation& op;
};

class VectorExpression::OpWrapper3 : public VectorExpression::VectorOperation {
public:
    OpWrapper3(const Operation& op) : op(op) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        Vec3 result;
        double args0[] = {args[0][0], args[1][0], args[2][0]};
        result[0] = op.evaluate(args0, noVariables);
        double args1[] = {args[0][1], args[1][1], args[2][1]};
        result[1] = op.evaluate(args1, noVariables);
        double args2[] = {args[0][2], args[1][2], args[2][2]};
        result[2] = op.evaluate(args2, noVariables);
        return result;
    }
private:
    const Operation& op;
};

class VectorExpression::Variable : public VectorExpression::VectorOperation {
public:
    Variable(const string& name) : name(name) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        map<string, Vec3>::const_iterator iter = variables.find(name);
        if (iter == variables.end())
            throw Exception("No value specified for variable "+name);
        return iter->second;
    }
private:
    string name;
};

class VectorExpression::Constant : public VectorExpression::VectorOperation {
public:
    Constant(double value) : value(Vec3(value, value, value)) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        return value;
    }
private:
    Vec3 value;
};

class VectorExpression::Dot : public VectorExpression::VectorOperation {
public:
    Dot() {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        double result = args[0].dot(args[1]);
        return Vec3(result, result, result);
    }
};

class VectorExpression::Cross : public VectorExpression::VectorOperation {
public:
    Cross() {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        return args[0].cross(args[1]);
    }
};

class VectorExpression::Component : public VectorExpression::VectorOperation {
public:
    Component(int index) : index(index) {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        double value = args[0][index];
        return Vec3(value, value, value);
    }
private:
    int index;
};

class VectorExpression::Vector : public VectorExpression::VectorOperation {
public:
    Vector() {
    }
    Vec3 evaluate(const Vec3* args, const map<string, Vec3>& variables) const {
        return Vec3(args[0][0], args[1][1], args[2][2]);
    }
};

VectorExpression::VectorExpression(const string& expression, const map<string, CustomFunction*>& customFunctions) {
    PlaceholderFunction fn1(1), fn2(2), fn3(3);
    PlaceholderFunction crossFunction(2);
    PlaceholderFunction xFunction(1);
    PlaceholderFunction yFunction(1);
    PlaceholderFunction zFunction(1);
    PlaceholderFunction vectorFunction(3);
    map<string, CustomFunction*> functions = customFunctions;
    functions["dot"] = &fn2;
    functions["cross"] = &fn2;
    functions["_x"] = &fn1;
    functions["_y"] = &fn1;
    functions["_z"] = &fn1;
    functions["vector"] = &fn3;
    analyzeExpression(Parser::parse(expression, functions));
}

VectorExpression::VectorExpression(const ParsedExpression& expression) {
    analyzeExpression(expression);
}

VectorExpression::VectorExpression(const VectorExpression& expression) {
    *this = expression;
}

void VectorExpression::analyzeExpression(const ParsedExpression& expression) {
    parsed = expression.optimize();
    program = parsed.createProgram();
    stack.resize(program.getStackSize()+1);
    for (int i = 0; i < program.getNumOperations(); i++) {
        const Operation& op = program.getOperation(i);
        if (op.getId() == Operation::VARIABLE)
            operations.push_back(new Variable(op.getName()));
        else if (op.getId() == Operation::CONSTANT)
            operations.push_back(new Constant(dynamic_cast<const Operation::Constant&>(op).getValue()));
        else if (op.getName() == "dot")
            operations.push_back(new Dot());
        else if (op.getName() == "cross")
            operations.push_back(new Cross());
        else if (op.getName() == "_x")
            operations.push_back(new Component(0));
        else if (op.getName() == "_y")
            operations.push_back(new Component(1));
        else if (op.getName() == "_z")
            operations.push_back(new Component(2));
        else if (op.getName() == "vector")
            operations.push_back(new Vector());
        else if (op.getNumArguments() == 1)
            operations.push_back(new OpWrapper1(op));
        else if (op.getNumArguments() == 2)
            operations.push_back(new OpWrapper2(op));
        else if (op.getNumArguments() == 3)
            operations.push_back(new OpWrapper3(op));
        else
            throw OpenMMException("Unsupported operator in vector expression: "+op.getName());
    }
}

VectorExpression::~VectorExpression() {
    for (int i = 0; i < operations.size(); i++)
        delete operations[i];
}

Vec3 VectorExpression::evaluate(const map<string, Vec3>& variables) const {
    int stackSize = program.getStackSize();
    int stackPointer = stackSize;
    for (int i = 0; i < (int) operations.size(); i++) {
        int numArgs = program.getOperation(i).getNumArguments();
        Vec3 result = operations[i]->evaluate(&stack[stackPointer], variables);
        stackPointer += numArgs-1;
        stack[stackPointer] = result;
    }
    return stack[stackSize-1];
}

VectorExpression& VectorExpression::operator=(const VectorExpression& exp) {
    analyzeExpression(exp.parsed);
    return *this;
}
