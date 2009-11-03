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

string OpenCLExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const map<string, string>& variables,
        const vector<pair<string, string> >& functions, const string& prefix, const string& functionParams) {
    stringstream out;
    vector<pair<ExpressionTreeNode, string> > temps;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter) {
        processExpression(out, iter->second.getRootNode(), temps, variables, functions, prefix, functionParams);
        out << iter->first << getTempName(iter->second.getRootNode(), temps) << ";\n";
    }
    return out.str();
}

void OpenCLExpressionUtilities::processExpression(stringstream& out, const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& temps,
        const map<string, string>& variables, const vector<pair<string, string> >& functions, const string& prefix, const string& functionParams) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        processExpression(out, node.getChildren()[i], temps, variables, functions, prefix, functionParams);
    string name = prefix+intToString(temps.size());
    
    out << "float " << name << " = ";
    switch (node.getOperation().getId()) {
        case Operation::CONSTANT:
            out << doubleToString(dynamic_cast<const Operation::Constant*>(&node.getOperation())->getValue());
            break;
        case Operation::VARIABLE:
        {
            map<string, string>::const_iterator iter = variables.find(node.getOperation().getName());
            if (iter == variables.end())
                throw OpenMMException("Unknown variable in expression: "+node.getOperation().getName());
            out << iter->second;
            break;
        }
        case Operation::CUSTOM:
        {
            int i;
            for (i = 0; i < (int) functions.size() && functions[i].first != node.getOperation().getName(); i++)
                ;
            if (i == functions.size())
                throw OpenMMException("Unknown function in expression: "+node.getOperation().getName());
            out << "0.0f;\n";
            out << "{\n";
            out << "float4 params = " << functionParams << "[" << i << "];\n";
            out << "float x = " << getTempName(node.getChildren()[0], temps) << ";\n";
            out << "if (x >= params.x && x <= params.y) {\n";
            out << "int index = (int) (floor((x-params.x)*params.z));\n";
            out << "float4 coeff = " << functions[i].second << "[index];\n";
            out << "x = (x-params.x)*params.z-index;\n";
            if (dynamic_cast<const Operation::Custom*>(&node.getOperation())->getDerivOrder()[0] == 0)
                out << name << " = coeff.x+x*(coeff.y+x*(coeff.z+x*coeff.w));\n";
            else
                out << name << " = (coeff.y+x*(2.0f*coeff.z+x*3.0f*coeff.w))*params.z;\n";
            out << "}\n";
            out << "}";
            break;
        }
        case Operation::ADD:
            out << getTempName(node.getChildren()[0], temps) << "+" << getTempName(node.getChildren()[1], temps);
            break;
        case Operation::SUBTRACT:
            out << getTempName(node.getChildren()[0], temps) << "-" << getTempName(node.getChildren()[1], temps);
            break;
        case Operation::MULTIPLY:
            out << getTempName(node.getChildren()[0], temps) << "*" << getTempName(node.getChildren()[1], temps);
            break;
        case Operation::DIVIDE:
            out << getTempName(node.getChildren()[0], temps) << "/" << getTempName(node.getChildren()[1], temps);
            break;
        case Operation::POWER:
            out << "pow(" << getTempName(node.getChildren()[0], temps) << ", " << getTempName(node.getChildren()[1], temps) << ")";
            break;
        case Operation::NEGATE:
            out << "-" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::SQRT:
            out << "sqrt(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::EXP:
            out << "exp(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::LOG:
            out << "log(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::SIN:
            out << "sin(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::COS:
            out << "cos(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::SEC:
            out << "1.0f/cos(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::CSC:
            out << "1.0f/sin(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::TAN:
            out << "tan(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::COT:
            out << "1.0f/tan(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::ASIN:
            out << "asin(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::ACOS:
            out << "acos(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::ATAN:
            out << "atan(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::SQUARE:
        {
            string arg = getTempName(node.getChildren()[0], temps);
            out << arg << "*" << arg;
            break;
        }
        case Operation::CUBE:
        {
            string arg = getTempName(node.getChildren()[0], temps);
            out << arg << "*" << arg << "*" << arg;
            break;
        }
        case Operation::RECIPROCAL:
            out << "1.0f/" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::ADD_CONSTANT:
            out << doubleToString(dynamic_cast<const Operation::AddConstant*>(&node.getOperation())->getValue()) << "+" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::MULTIPLY_CONSTANT:
            out << doubleToString(dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue()) << "*" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::POWER_CONSTANT:
        {
            double exponent = dynamic_cast<const Operation::PowerConstant*>(&node.getOperation())->getValue();
            if (exponent == 0.0)
                out << "1.0f";
            else if (exponent == (int) exponent) {
                out << "0.0f;\n";
                out << "{\n";
                out << "float multiplier = " << (exponent < 0.0 ? "1.0f/" : "") << getTempName(node.getChildren()[0], temps) << ";\n";
                int exp = (int) fabs(exponent);
                bool hasAssigned = false;
                while (exp != 0) {
                    if (exp%2 == 1) {
                        if (!hasAssigned)
                            out << name << " = multiplier;\n";
                        else
                            out << name << " *= multiplier;\n";
                        hasAssigned = true;
                    }
                    exp >>= 1;
                    if (exp != 0)
                        out << "multiplier *= multiplier;\n";
                }
                out << "}";
            }
            else
                out << "pow(" << getTempName(node.getChildren()[0], temps) << ", " << doubleToString(exponent) << ")";
            break;
        }
        default:
            throw OpenMMException("Internal error: Unknown operation in user-defined expression: "+node.getOperation().getName());
    }
    out << ";\n";
    temps.push_back(make_pair(node, name));
}

string OpenCLExpressionUtilities::getTempName(const ExpressionTreeNode& node, const vector<pair<ExpressionTreeNode, string> >& temps) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return temps[i].second;
    stringstream out;
    out << "Internal error: No temporary variable for expression node: " << node;
    throw OpenMMException(out.str());
}
