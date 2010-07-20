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

string OpenCLExpressionUtilities::doubleToString(double value) {
    stringstream s;
    s.precision(8);
    s << scientific << value << "f";
    return s.str();
}

string OpenCLExpressionUtilities::intToString(int value) {
    stringstream s;
    s << value;
    return s.str();
}

string OpenCLExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const map<string, string>& variables,
        const vector<pair<string, string> >& functions, const string& prefix, const string& functionParams) {
    stringstream out;
    vector<ParsedExpression> allExpressions;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter)
        allExpressions.push_back(iter->second);
    vector<pair<ExpressionTreeNode, string> > temps;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter) {
        processExpression(out, iter->second.getRootNode(), temps, variables, functions, prefix, functionParams, allExpressions);
        out << iter->first << getTempName(iter->second.getRootNode(), temps) << ";\n";
    }
    return out.str();
}

void OpenCLExpressionUtilities::processExpression(stringstream& out, const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& temps,
        const map<string, string>& variables, const vector<pair<string, string> >& functions, const string& prefix, const string& functionParams,
        const vector<ParsedExpression>& allExpressions) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        processExpression(out, node.getChildren()[i], temps, variables, functions, prefix, functionParams, allExpressions);
    string name = prefix+intToString(temps.size());
    bool hasRecordedNode = false;
    
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
            bool isDeriv = (dynamic_cast<const Operation::Custom*>(&node.getOperation())->getDerivOrder()[0] == 1);
            out << "0.0f;\n";
            temps.push_back(make_pair(node, name));
            hasRecordedNode = true;

            // If both the value and derivative of the function are needed, it's faster to calculate them both
            // at once, so check to see if both are needed.

            const ExpressionTreeNode* valueNode = NULL;
            const ExpressionTreeNode* derivNode = NULL;
            for (int j = 0; j < (int) allExpressions.size(); j++)
                findRelatedTabulatedFunctions(node, allExpressions[j].getRootNode(), valueNode, derivNode);
            string valueName = name;
            string derivName = name;
            if (valueNode != NULL && derivNode != NULL) {
                string name2 = prefix+intToString(temps.size());
                out << "float " << name2 << " = 0.0f;\n";
                if (isDeriv) {
                    valueName = name2;
                    temps.push_back(make_pair(*valueNode, name2));
                }
                else {
                    derivName = name2;
                    temps.push_back(make_pair(*derivNode, name2));
                }
            }
            out << "{\n";
            out << "float4 params = " << functionParams << "[" << i << "];\n";
            out << "float x = " << getTempName(node.getChildren()[0], temps) << ";\n";
            out << "if (x >= params.x && x <= params.y) {\n";
            out << "int index = (int) (floor((x-params.x)*params.z));\n";
            out << "float4 coeff = " << functions[i].second << "[index];\n";
            out << "x = (x-params.x)*params.z-index;\n";
            if (valueNode != NULL)
                out << valueName << " = coeff.x+x*(coeff.y+x*(coeff.z+x*coeff.w));\n";
            if (derivNode != NULL)
                out << derivName << " = (coeff.y+x*(2.0f*coeff.z+x*3.0f*coeff.w))*params.z;\n";
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
            out << "EXP(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::LOG:
            out << "LOG(" << getTempName(node.getChildren()[0], temps) << ")";
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
        case Operation::SINH:
            out << "sinh(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::COSH:
            out << "cosh(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::TANH:
            out << "tanh(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::ERF:
            out << "erf(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::ERFC:
            out << "erfc(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        case Operation::STEP:
            out << getTempName(node.getChildren()[0], temps) << " >= 0.0f ? 1.0f : 0.0f";
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
            out << "RECIP(" << getTempName(node.getChildren()[0], temps) << ")";
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
                temps.push_back(make_pair(node, name));
                hasRecordedNode = true;

                // If multiple integral powers of the same base are needed, it's faster to calculate all of them
                // at once, so check to see if others are also needed.

                map<int, const ExpressionTreeNode*> powers;
                powers[(int) exponent] = &node;
                for (int j = 0; j < (int) allExpressions.size(); j++)
                    findRelatedPowers(node, allExpressions[j].getRootNode(), powers);
                vector<int> exponents;
                vector<string> names;
                vector<bool> hasAssigned(powers.size(), false);
                exponents.push_back((int) fabs(exponent));
                names.push_back(name);
                for (map<int, const ExpressionTreeNode*>::const_iterator iter = powers.begin(); iter != powers.end(); ++iter) {
                    if (iter->first != exponent) {
                        exponents.push_back(iter->first >= 0 ? iter->first : -iter->first);
                        string name2 = prefix+intToString(temps.size());
                        names.push_back(name2);
                        temps.push_back(make_pair(*iter->second, name2));
                        out << "float " << name2 << " = 0.0f;\n";
                    }
                }
                out << "{\n";
                out << "float multiplier = " << (exponent < 0.0 ? "1.0f/" : "") << getTempName(node.getChildren()[0], temps) << ";\n";
                bool done = false;
                while (!done) {
                    done = true;
                    for (int i = 0; i < (int) exponents.size(); i++) {
                        if (exponents[i]%2 == 1) {
                            if (!hasAssigned[i])
                                out << names[i] << " = multiplier;\n";
                            else
                                out << names[i] << " *= multiplier;\n";
                            hasAssigned[i] = true;
                        }
                        exponents[i] >>= 1;
                        if (exponents[i] != 0)
                            done = false;
                    }
                    if (!done)
                        out << "multiplier *= multiplier;\n";
                }
                out << "}";
            }
            else
                out << "pow(" << getTempName(node.getChildren()[0], temps) << ", " << doubleToString(exponent) << ")";
            break;
        }
        case Operation::MIN:
            out << "min(" << getTempName(node.getChildren()[0], temps) << ", " << getTempName(node.getChildren()[1], temps) << ")";
            break;
        case Operation::MAX:
            out << "max(" << getTempName(node.getChildren()[0], temps) << ", " << getTempName(node.getChildren()[1], temps) << ")";
            break;
        case Operation::ABS:
            out << "fabs(" << getTempName(node.getChildren()[0], temps) << ")";
            break;
        default:
            throw OpenMMException("Internal error: Unknown operation in user-defined expression: "+node.getOperation().getName());
    }
    out << ";\n";
    if (!hasRecordedNode)
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

void OpenCLExpressionUtilities::findRelatedTabulatedFunctions(const ExpressionTreeNode& node, const ExpressionTreeNode& searchNode,
            const ExpressionTreeNode*& valueNode, const ExpressionTreeNode*& derivNode) {
    if (searchNode.getOperation().getId() == Operation::CUSTOM && node.getChildren()[0] == searchNode.getChildren()[0]) {
        if (dynamic_cast<const Operation::Custom*>(&searchNode.getOperation())->getDerivOrder()[0] == 0)
            valueNode = &searchNode;
        else
            derivNode = &searchNode;
    }
    else
        for (int i = 0; i < (int) searchNode.getChildren().size(); i++)
            findRelatedTabulatedFunctions(node, searchNode.getChildren()[i], valueNode, derivNode);
}

void OpenCLExpressionUtilities::findRelatedPowers(const ExpressionTreeNode& node, const ExpressionTreeNode& searchNode, map<int, const ExpressionTreeNode*>& powers) {
    if (searchNode.getOperation().getId() == Operation::POWER_CONSTANT && node.getChildren()[0] == searchNode.getChildren()[0]) {
        double realPower = dynamic_cast<const Operation::PowerConstant*>(&searchNode.getOperation())->getValue();
        int power = (int) realPower;
        if (power != realPower)
            return; // We are only interested in integer powers.
        if (powers.find(power) != powers.end())
            return; // This power is already in the map.
        if (powers.begin()->first*power < 0)
            return; // All powers must have the same sign.
        powers[power] = &searchNode;
    }
    else
        for (int i = 0; i < (int) searchNode.getChildren().size(); i++)
            findRelatedPowers(node, searchNode.getChildren()[i], powers);
}

vector<mm_float4> OpenCLExpressionUtilities::computeFunctionCoefficients(const vector<double>& values, bool interpolating) {
    // First create a padded set of function values.

    vector<double> padded(values.size()+2);
    padded[0] = 2*values[0]-values[1];
    for (int i = 0; i < (int) values.size(); i++)
        padded[i+1] = values[i];
    padded[padded.size()-1] = 2*values[values.size()-1]-values[values.size()-2];

    // Now compute the spline coefficients.

    vector<mm_float4> f(values.size()-1);
    for (int i = 0; i < (int) values.size()-1; i++) {
        if (interpolating)
            f[i] = mm_float4((cl_float) padded[i+1],
                                (cl_float) (0.5*(-padded[i]+padded[i+2])),
                                (cl_float) (0.5*(2.0*padded[i]-5.0*padded[i+1]+4.0*padded[i+2]-padded[i+3])),
                                (cl_float) (0.5*(-padded[i]+3.0*padded[i+1]-3.0*padded[i+2]+padded[i+3])));
        else
            f[i] = mm_float4((cl_float) ((padded[i]+4.0*padded[i+1]+padded[i+2])/6.0),
                                (cl_float) ((-3.0*padded[i]+3.0*padded[i+2])/6.0),
                                (cl_float) ((3.0*padded[i]-6.0*padded[i+1]+3.0*padded[i+2])/6.0),
                                (cl_float) ((-padded[i]+3.0*padded[i+1]-3.0*padded[i+2]+padded[i+3])/6.0));
    }
    return f;
}
