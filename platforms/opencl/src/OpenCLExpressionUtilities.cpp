/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "OpenCLExpressionUtilities.h"
#include "openmm/OpenMMException.h"
#include "lepton/Operation.h"
#include <sstream>

using namespace OpenMM;
using namespace Lepton;
using namespace std;

static string doubleToString(double value) {
    stringstream s;
    s.precision(8);
    s << scientific << value << "f";
    return s.str();
}

static string intToString(int value) {
    stringstream s;
    s << value;
    return s.str();
}

string OpenCLExpressionUtilities::createExpression(const ParsedExpression& expression, const map<string, string>& variables) {
    return processExpression(expression.getRootNode(), variables);
}

string OpenCLExpressionUtilities::processExpression(const ExpressionTreeNode& node, const map<string, string>& variables) {
    switch (node.getOperation().getId()) {
        case Operation::CONSTANT:
            return doubleToString(dynamic_cast<const Operation::Constant*>(&node.getOperation())->getValue());
        case Operation::VARIABLE:
        {
            map<string, string>::const_iterator iter = variables.find(node.getOperation().getName());
            if (iter == variables.end())
                throw OpenMMException("Unknown variable in expression: "+node.getOperation().getName());
            return iter->second;
        }
        case Operation::ADD:
             return "("+processExpression(node.getChildren()[0], variables)+")+("+processExpression(node.getChildren()[1], variables)+")";
        case Operation::SUBTRACT:
            return "("+processExpression(node.getChildren()[0], variables)+")-("+processExpression(node.getChildren()[1], variables)+")";
        case Operation::MULTIPLY:
            return "("+processExpression(node.getChildren()[0], variables)+")*("+processExpression(node.getChildren()[1], variables)+")";
        case Operation::DIVIDE:
            return "("+processExpression(node.getChildren()[0], variables)+")/("+processExpression(node.getChildren()[1], variables)+")";
        case Operation::POWER:
            return "pow(("+processExpression(node.getChildren()[0], variables)+"), ("+processExpression(node.getChildren()[1], variables)+"))";
        case Operation::NEGATE:
            return "-("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::SQRT:
            return "sqrt("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::EXP:
            return "exp("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::LOG:
            return "log("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::SIN:
            return "sin("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::COS:
            return "cos("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::SEC:
            return "1.0f/cos("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::CSC:
            return "1.0f/sin("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::TAN:
            return "tan("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::COT:
            return "1.0f/tan("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::ASIN:
            return "asin("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::ACOS:
            return "acos("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::ATAN:
            return "atan("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::SQUARE:
            return "pow(("+processExpression(node.getChildren()[0], variables)+"), 2.0f)";
        case Operation::CUBE:
            return "pow(("+processExpression(node.getChildren()[0], variables)+"), 3.0f)";
        case Operation::RECIPROCAL:
            return "1.0f/("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::ADD_CONSTANT:
            return doubleToString(dynamic_cast<const Operation::AddConstant*>(&node.getOperation())->getValue())+"+("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::MULTIPLY_CONSTANT:
            return doubleToString(dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue())+"*("+processExpression(node.getChildren()[0], variables)+")";
        case Operation::POWER_CONSTANT:
            return "pow(("+processExpression(node.getChildren()[0], variables)+"), "+doubleToString(dynamic_cast<const Operation::PowerConstant*>(&node.getOperation())->getValue())+")";
    }
    throw OpenMMException("Internal error: Unknown operation in user-defined expression: "+node.getOperation().getName());
}
