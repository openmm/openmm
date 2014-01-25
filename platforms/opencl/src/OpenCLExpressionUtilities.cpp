/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2014 Stanford University and the Authors.      *
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
#include "openmm/internal/SplineFitter.h"
#include "lepton/Operation.h"

using namespace OpenMM;
using namespace Lepton;
using namespace std;

OpenCLExpressionUtilities::OpenCLExpressionUtilities(OpenCLContext& context) : context(context), fp1(1), fp2(2), fp3(3) {
}

string OpenCLExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const map<string, string>& variables,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix,
        const string& functionParams, const string& tempType) {
    vector<pair<ExpressionTreeNode, string> > variableNodes;
    for (map<string, string>::const_iterator iter = variables.begin(); iter != variables.end(); ++iter)
        variableNodes.push_back(make_pair(ExpressionTreeNode(new Operation::Variable(iter->first)), iter->second));
    return createExpressions(expressions, variableNodes, functions, functionNames, prefix, functionParams, tempType);
}

string OpenCLExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const vector<pair<ExpressionTreeNode, string> >& variables,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix,
        const string& functionParams, const string& tempType) {
    stringstream out;
    vector<ParsedExpression> allExpressions;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter)
        allExpressions.push_back(iter->second);
    vector<pair<ExpressionTreeNode, string> > temps = variables;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter) {
        processExpression(out, iter->second.getRootNode(), temps, functions, functionNames, prefix, functionParams, allExpressions, tempType);
        out << iter->first << getTempName(iter->second.getRootNode(), temps) << ";\n";
    }
    return out.str();
}

void OpenCLExpressionUtilities::processExpression(stringstream& out, const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& temps,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix, const string& functionParams,
        const vector<ParsedExpression>& allExpressions, const string& tempType) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        processExpression(out, node.getChildren()[i], temps, functions, functionNames, prefix, functionParams, allExpressions, tempType);
    string name = prefix+context.intToString(temps.size());
    bool hasRecordedNode = false;
    
    out << tempType << " " << name << " = ";
    switch (node.getOperation().getId()) {
        case Operation::CONSTANT:
            out << context.doubleToString(dynamic_cast<const Operation::Constant*>(&node.getOperation())->getValue());
            break;
        case Operation::VARIABLE:
            throw OpenMMException("Unknown variable in expression: "+node.getOperation().getName());
        case Operation::CUSTOM:
        {
            int i;
            for (i = 0; i < (int) functionNames.size() && functionNames[i].first != node.getOperation().getName(); i++)
                ;
            if (i == functionNames.size())
                throw OpenMMException("Unknown function in expression: "+node.getOperation().getName());
            out << "0.0f;\n";
            temps.push_back(make_pair(node, name));
            hasRecordedNode = true;

            // If both the value and derivative of the function are needed, it's faster to calculate them both
            // at once, so check to see if both are needed.

            vector<const ExpressionTreeNode*> nodes;
            for (int j = 0; j < (int) allExpressions.size(); j++)
                findRelatedTabulatedFunctions(node, allExpressions[j].getRootNode(), nodes);
            vector<string> nodeNames;
            nodeNames.push_back(name);
            for (int j = 1; j < (int) nodes.size(); j++) {
                string name2 = prefix+context.intToString(temps.size());
                out << tempType << " " << name2 << " = 0.0f;\n";
                nodeNames.push_back(name2);
                temps.push_back(make_pair(*nodes[j], name2));
            }
            out << "{\n";
            if (dynamic_cast<const Continuous1DFunction*>(functions[i]) != NULL) {
                out << "float4 params = " << functionParams << "[" << i << "];\n";
                out << "real x = " << getTempName(node.getChildren()[0], temps) << ";\n";
                out << "if (x >= params.x && x <= params.y) {\n";
                out << "x = (x-params.x)*params.z;\n";
                out << "int index = (int) (floor(x));\n";
                out << "index = min(index, (int) params.w);\n";
                out << "float4 coeff = " << functionNames[i].second << "[index];\n";
                out << "real b = x-index;\n";
                out << "real a = 1.0f-b;\n";
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0)
                        out << nodeNames[j] << " = a*coeff.x+b*coeff.y+((a*a*a-a)*coeff.z+(b*b*b-b)*coeff.w)/(params.z*params.z);\n";
                    else
                        out << nodeNames[j] << " = (coeff.y-coeff.x)*params.z+((1.0f-3.0f*a*a)*coeff.z+(3.0f*b*b-1.0f)*coeff.w)/params.z;\n";
                }
                out << "}\n";
            }
            else if (dynamic_cast<const Discrete1DFunction*>(functions[i]) != NULL) {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0) {
                        out << "float4 params = " << functionParams << "[" << i << "];\n";
                        out << "real x = " << getTempName(node.getChildren()[0], temps) << ";\n";
                        out << "if (x >= 0 && x < params.x) {\n";
                        out << "int index = (int) round(x);\n";
                        out << nodeNames[j] << " = " << functionNames[i].second << "[index];\n";
                        out << "}\n";
                    }
                }
            }
            else if (dynamic_cast<const Discrete2DFunction*>(functions[i]) != NULL) {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0 && derivOrder[1] == 0) {
                        out << "float4 params = " << functionParams << "[" << i << "];\n";
                        out << "int x = (int) round(" << getTempName(node.getChildren()[0], temps) << ");\n";
                        out << "int y = (int) round(" << getTempName(node.getChildren()[1], temps) << ");\n";
                        out << "int xsize = (int) params.x;\n";
                        out << "int ysize = (int) params.y;\n";
                        out << "int index = x+y*xsize;\n";
                        out << "if (index >= 0 && index < xsize*ysize)\n";
                        out << nodeNames[j] << " = " << functionNames[i].second << "[index];\n";
                    }
                }
            }
            else if (dynamic_cast<const Discrete3DFunction*>(functions[i]) != NULL) {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 0) {
                        out << "float4 params = " << functionParams << "[" << i << "];\n";
                        out << "int x = (int) round(" << getTempName(node.getChildren()[0], temps) << ");\n";
                        out << "int y = (int) round(" << getTempName(node.getChildren()[1], temps) << ");\n";
                        out << "int z = (int) round(" << getTempName(node.getChildren()[2], temps) << ");\n";
                        out << "int xsize = (int) params.x;\n";
                        out << "int ysize = (int) params.y;\n";
                        out << "int zsize = (int) params.z;\n";
                        out << "int index = x+(y+z*ysize)*xsize;\n";
                        out << "if (index >= 0 && index < xsize*ysize*zsize)\n";
                        out << nodeNames[j] << " = " << functionNames[i].second << "[index];\n";
                    }
                }
            }
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
        {
            bool haveReciprocal = false;
            for (int i = 0; i < (int) temps.size(); i++)
                if (temps[i].first.getOperation().getId() == Operation::RECIPROCAL && temps[i].first.getChildren()[0] == node.getChildren()[1]) {
                    haveReciprocal = true;
                    out << getTempName(node.getChildren()[0], temps) << "*" << temps[i].second;
                }
            if (!haveReciprocal)
                out << getTempName(node.getChildren()[0], temps) << "/" << getTempName(node.getChildren()[1], temps);
            break;
        }
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
        case Operation::DELTA:
            out << getTempName(node.getChildren()[0], temps) << " == 0.0f ? 1.0f : 0.0f";
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
            out << context.doubleToString(dynamic_cast<const Operation::AddConstant*>(&node.getOperation())->getValue()) << "+" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::MULTIPLY_CONSTANT:
            out << context.doubleToString(dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue()) << "*" << getTempName(node.getChildren()[0], temps);
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
                        string name2 = prefix+context.intToString(temps.size());
                        names.push_back(name2);
                        temps.push_back(make_pair(*iter->second, name2));
                        out << tempType << " " << name2 << " = 0.0f;\n";
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
                out << "pow(" << getTempName(node.getChildren()[0], temps) << ", " << context.doubleToString(exponent) << ")";
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
            vector<const Lepton::ExpressionTreeNode*>& nodes) {
    if (searchNode.getOperation().getId() == Operation::CUSTOM && node.getChildren()[0] == searchNode.getChildren()[0])
        nodes.push_back(&searchNode);
    else
        for (int i = 0; i < (int) searchNode.getChildren().size(); i++)
            findRelatedTabulatedFunctions(node, searchNode.getChildren()[i], nodes);
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

vector<float> OpenCLExpressionUtilities::computeFunctionCoefficients(const TabulatedFunction& function, int& width) {
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL) {
        // Compute the spline coefficients.

        const Continuous1DFunction& fn = dynamic_cast<const Continuous1DFunction&>(function);
        vector<double> values;
        double min, max;
        fn.getFunctionParameters(values, min, max);
        int numValues = values.size();
        vector<double> x(numValues), derivs;
        for (int i = 0; i < numValues; i++)
            x[i] = min+i*(max-min)/(numValues-1);
        SplineFitter::createNaturalSpline(x, values, derivs);
        vector<float> f(4*(numValues-1));
        for (int i = 0; i < (int) values.size()-1; i++) {
            f[4*i] = (float) values[i];
            f[4*i+1] = (float) values[i+1];
            f[4*i+2] = (float) (derivs[i]/6.0);
            f[4*i+3] = (float) (derivs[i+1]/6.0);
        }
        width = 4;
        return f;
    }
    if (dynamic_cast<const Discrete1DFunction*>(&function) != NULL) {
        // Record the tabulated values.
        
        const Discrete1DFunction& fn = dynamic_cast<const Discrete1DFunction&>(function);
        vector<double> values;
        fn.getFunctionParameters(values);
        int numValues = values.size();
        vector<float> f(numValues);
        for (int i = 0; i < numValues; i++)
            f[i] = (float) values[i];
        width = 1;
        return f;
    }
    if (dynamic_cast<const Discrete2DFunction*>(&function) != NULL) {
        // Record the tabulated values.
        
        const Discrete2DFunction& fn = dynamic_cast<const Discrete2DFunction&>(function);
        int xsize, ysize;
        vector<double> values;
        fn.getFunctionParameters(xsize, ysize, values);
        int numValues = values.size();
        vector<float> f(numValues);
        for (int i = 0; i < numValues; i++)
            f[i] = (float) values[i];
        width = 1;
        return f;
    }
    if (dynamic_cast<const Discrete3DFunction*>(&function) != NULL) {
        // Record the tabulated values.
        
        const Discrete3DFunction& fn = dynamic_cast<const Discrete3DFunction&>(function);
        int xsize, ysize, zsize;
        vector<double> values;
        fn.getFunctionParameters(xsize, ysize, zsize, values);
        int numValues = values.size();
        vector<float> f(numValues);
        for (int i = 0; i < numValues; i++)
            f[i] = (float) values[i];
        width = 1;
        return f;
    }
    throw OpenMMException("computeFunctionCoefficients: Unknown function type");
}

vector<mm_float4> OpenCLExpressionUtilities::computeFunctionParameters(const vector<const TabulatedFunction*>& functions) {
    vector<mm_float4> params(functions.size());
    for (int i = 0; i < (int) functions.size(); i++) {
        if (dynamic_cast<const Continuous1DFunction*>(functions[i]) != NULL) {
            const Continuous1DFunction& fn = dynamic_cast<const Continuous1DFunction&>(*functions[i]);
            vector<double> values;
            double min, max;
            fn.getFunctionParameters(values, min, max);
            params[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), (float) values.size()-2);
        }
        else if (dynamic_cast<const Discrete1DFunction*>(functions[i]) != NULL) {
            const Discrete1DFunction& fn = dynamic_cast<const Discrete1DFunction&>(*functions[i]);
            vector<double> values;
            fn.getFunctionParameters(values);
            params[i] = mm_float4((float) values.size(), 0.0f, 0.0f, 0.0f);
        }
        else if (dynamic_cast<const Discrete2DFunction*>(functions[i]) != NULL) {
            const Discrete2DFunction& fn = dynamic_cast<const Discrete2DFunction&>(*functions[i]);
            int xsize, ysize;
            vector<double> values;
            fn.getFunctionParameters(xsize, ysize, values);
            params[i] = mm_float4(xsize, ysize, 0.0f, 0.0f);
        }
        else if (dynamic_cast<const Discrete3DFunction*>(functions[i]) != NULL) {
            const Discrete3DFunction& fn = dynamic_cast<const Discrete3DFunction&>(*functions[i]);
            int xsize, ysize, zsize;
            vector<double> values;
            fn.getFunctionParameters(xsize, ysize, zsize, values);
            params[i] = mm_float4(xsize, ysize, zsize, 0.0f);
        }
        else
            throw OpenMMException("computeFunctionParameters: Unknown function type");
    }
    return params;
}

Lepton::CustomFunction* OpenCLExpressionUtilities::getFunctionPlaceholder(const TabulatedFunction& function) {
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL)
        return &fp1;
    if (dynamic_cast<const Discrete1DFunction*>(&function) != NULL)
        return &fp1;
    if (dynamic_cast<const Discrete2DFunction*>(&function) != NULL)
        return &fp2;
    if (dynamic_cast<const Discrete3DFunction*>(&function) != NULL)
        return &fp3;
    throw OpenMMException("getFunctionPlaceholder: Unknown function type");
}
