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

#include "CudaExpressionUtilities.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/SplineFitter.h"
#include "lepton/Operation.h"

using namespace OpenMM;
using namespace Lepton;
using namespace std;

CudaExpressionUtilities::CudaExpressionUtilities(CudaContext& context) : context(context), fp1(1), fp2(2), fp3(3), periodicDistance(6) {
}

string CudaExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const map<string, string>& variables,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix, const string& tempType) {
    vector<pair<ExpressionTreeNode, string> > variableNodes;
    for (map<string, string>::const_iterator iter = variables.begin(); iter != variables.end(); ++iter)
        variableNodes.push_back(make_pair(ExpressionTreeNode(new Operation::Variable(iter->first)), iter->second));
    return createExpressions(expressions, variableNodes, functions, functionNames, prefix, tempType);
}

string CudaExpressionUtilities::createExpressions(const map<string, ParsedExpression>& expressions, const vector<pair<ExpressionTreeNode, string> >& variables,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix, const string& tempType) {
    stringstream out;
    vector<ParsedExpression> allExpressions;
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter)
        allExpressions.push_back(iter->second);
    vector<pair<ExpressionTreeNode, string> > temps = variables;
    vector<vector<double> > functionParams = computeFunctionParameters(functions);
    for (map<string, ParsedExpression>::const_iterator iter = expressions.begin(); iter != expressions.end(); ++iter) {
        processExpression(out, iter->second.getRootNode(), temps, functions, functionNames, prefix, functionParams, allExpressions, tempType);
        out << iter->first << getTempName(iter->second.getRootNode(), temps) << ";\n";
    }
    return out.str();
}

void CudaExpressionUtilities::processExpression(stringstream& out, const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& temps,
        const vector<const TabulatedFunction*>& functions, const vector<pair<string, string> >& functionNames, const string& prefix, const vector<vector<double> >& functionParams,
        const vector<ParsedExpression>& allExpressions, const string& tempType) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        processExpression(out, node.getChildren()[i], temps, functions, functionNames, prefix, functionParams, allExpressions, tempType);
    string name = prefix+context.intToString(temps.size());
    bool hasRecordedNode = false;
    bool isVecType = (tempType[tempType.size()-1] == '3');
    
    out << tempType << " " << name << " = ";
    switch (node.getOperation().getId()) {
        case Operation::CONSTANT:
        {
            string value = context.doubleToString(dynamic_cast<const Operation::Constant*>(&node.getOperation())->getValue());
            if (isVecType)
                out << "make_" << tempType << "(" << value << ")";
            else
                out << value;
            break;
        }
        case Operation::VARIABLE:
            throw OpenMMException("Unknown variable in expression: "+node.getOperation().getName());
        case Operation::CUSTOM:
        {
            if (isVecType)
                out << "make_" << tempType << "(0);\n";
            else
                out << "0;\n";
            temps.push_back(make_pair(node, name));
            hasRecordedNode = true;

            // If both the value and derivative of the function are needed, it's faster to calculate them both
            // at once, so check to see if both are needed.

            vector<const ExpressionTreeNode*> nodes;
            nodes.push_back(&node);
            for (int j = 0; j < (int) allExpressions.size(); j++)
                findRelatedCustomFunctions(node, allExpressions[j].getRootNode(), nodes);
            vector<string> nodeNames;
            nodeNames.push_back(name);
            for (int j = 1; j < (int) nodes.size(); j++) {
                string name2 = prefix+context.intToString(temps.size());
                out << tempType << " " << name2 << " = 0.0f;\n";
                nodeNames.push_back(name2);
                temps.push_back(make_pair(*nodes[j], name2));
            }
            out << "{\n";
            if (node.getOperation().getName() == "periodicdistance") {
                // This is the periodicdistance() function.

                out << tempType << "3 periodicDistance_delta = make_real3(";
                for (int i = 0; i < 3; i++) {
                    if (i > 0)
                        out << ", ";
                    out << getTempName(node.getChildren()[i], temps) << "-" << getTempName(node.getChildren()[i+3], temps);
                }
                out << ");\n";
                out << "APPLY_PERIODIC_TO_DELTA(periodicDistance_delta)\n";
                out << tempType << " periodicDistance_r2 = periodicDistance_delta.x*periodicDistance_delta.x + periodicDistance_delta.y*periodicDistance_delta.y + periodicDistance_delta.z*periodicDistance_delta.z;\n";
                out << tempType << " periodicDistance_rinv = RSQRT(periodicDistance_r2);\n";
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    int argIndex = -1;
                    for (int k = 0; k < 6; k++) {
                        if (derivOrder[k] > 0) {
                            if (derivOrder[k] > 1 || argIndex != -1)
                                throw OpenMMException("Unsupported derivative of periodicdistance"); // Should be impossible for this to happen.
                            argIndex = k;
                        }
                    }
                    if (argIndex == -1)
                        out << nodeNames[j] << " = RECIP(periodicDistance_rinv);\n";
                    else if (argIndex == 0)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? periodicDistance_delta.x*periodicDistance_rinv : 0);\n";
                    else if (argIndex == 1)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? periodicDistance_delta.y*periodicDistance_rinv : 0);\n";
                    else if (argIndex == 2)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? periodicDistance_delta.z*periodicDistance_rinv : 0);\n";
                    else if (argIndex == 3)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? -periodicDistance_delta.x*periodicDistance_rinv : 0);\n";
                    else if (argIndex == 4)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? -periodicDistance_delta.y*periodicDistance_rinv : 0);\n";
                    else if (argIndex == 5)
                        out << nodeNames[j] << " = (periodicDistance_r2 > 0 ? -periodicDistance_delta.z*periodicDistance_rinv : 0);\n";
                }
            }
            else if (node.getOperation().getName() == "dot") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    string child1 = getTempName(node.getChildren()[0], temps);
                    string child2 = getTempName(node.getChildren()[1], temps);
                    if (derivOrder[0] == 0 && derivOrder[1] == 0)
                        out << nodeNames[j] << " = make_" << tempType << "(dot(" << child1 << ", " << child2 << "));\n";
                    else
                        throw OpenMMException("Unsupported derivative order for dot()");
                }
            }
            else if (node.getOperation().getName() == "cross") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    string child1 = getTempName(node.getChildren()[0], temps);
                    string child2 = getTempName(node.getChildren()[1], temps);
                    if (derivOrder[0] == 0 && derivOrder[1] == 0)
                        out << nodeNames[j] << " = cross(" << child1 << ", " << child2 << ");\n";
                    else
                        throw OpenMMException("Unsupported derivative order for cross()");
                }
            }
            else if (node.getOperation().getName() == "vector") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 0) {
                        out << nodeNames[j] << ".x = " << getTempName(node.getChildren()[0], temps) << ".x;\n";
                        out << nodeNames[j] << ".y = " << getTempName(node.getChildren()[1], temps) << ".y;\n";
                        out << nodeNames[j] << ".z = " << getTempName(node.getChildren()[2], temps) << ".z;\n";
                    }
                    else if (derivOrder[0] == 1 && derivOrder[1] == 0 && derivOrder[2] == 0)
                        out << nodeNames[j] << ".x = 1;\n";
                    else if (derivOrder[0] == 0 && derivOrder[1] == 1 && derivOrder[2] == 0)
                        out << nodeNames[j] << ".y = 1;\n";
                    else if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 1)
                        out << nodeNames[j] << ".z = 1;\n";
                }
            }
            else if (node.getOperation().getName() == "_x") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0)
                        out << nodeNames[j] << " = make_" << tempType << "(" << getTempName(node.getChildren()[0], temps) << ".x);\n";
                    else
                        throw OpenMMException("Unsupported derivative order for _x()");
                }
            }
            else if (node.getOperation().getName() == "_y") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0)
                        out << nodeNames[j] << " = make_" << tempType << "(" << getTempName(node.getChildren()[0], temps) << ".y);\n";
                    else
                        throw OpenMMException("Unsupported derivative order for _y()");
                }
            }
            else if (node.getOperation().getName() == "_z") {
                for (int j = 0; j < nodes.size(); j++) {
                    const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                    if (derivOrder[0] == 0)
                        out << nodeNames[j] << " = make_" << tempType << "(" << getTempName(node.getChildren()[0], temps) << ".z);\n";
                    else
                        throw OpenMMException("Unsupported derivative order for _z()");
                }
            }
            else {
                // This is a tabulated function.
                
                int i;
                for (i = 0; i < (int) functionNames.size() && functionNames[i].first != node.getOperation().getName(); i++)
                    ;
                if (i == functionNames.size())
                    throw OpenMMException("Unknown function in expression: "+node.getOperation().getName());
                vector<string> paramsFloat, paramsInt;
                for (int j = 0; j < (int) functionParams[i].size(); j++) {
                    paramsFloat.push_back(context.doubleToString(functionParams[i][j]));
                    paramsInt.push_back(context.intToString((int) functionParams[i][j]));
                }
                vector<string> suffixes;
                if (isVecType) {
                    suffixes.push_back(".x");
                    suffixes.push_back(".y");
                    suffixes.push_back(".z");
                }
                else {
                    suffixes.push_back("");
                }
                for (auto& suffix : suffixes) {
                    out << "{\n";
                    if (dynamic_cast<const Continuous1DFunction*>(functions[i]) != NULL) {
                        out << "real x = " << getTempName(node.getChildren()[0], temps) << suffix << ";\n";
                        out << "if (x >= " << paramsFloat[0] << " && x <= " << paramsFloat[1] << ") {\n";
                        out << "x = (x - " << paramsFloat[0] << ")*" << paramsFloat[2] << ";\n";
                        out << "int index = (int) (floor(x));\n";
                        out << "index = min(index, (int) " << paramsInt[3] << ");\n";
                        out << "float4 coeff = " << functionNames[i].second << "[index];\n";
                        out << "real b = x-index;\n";
                        out << "real a = 1.0f-b;\n";
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0)
                                out << nodeNames[j] << suffix << " = a*coeff.x+b*coeff.y+((a*a*a-a)*coeff.z+(b*b*b-b)*coeff.w)/(" << paramsFloat[2] << "*" << paramsFloat[2] << ");\n";
                            else
                                out << nodeNames[j] << suffix << " = (coeff.y-coeff.x)*" << paramsFloat[2] << "+((1.0f-3.0f*a*a)*coeff.z+(3.0f*b*b-1.0f)*coeff.w)/" << paramsFloat[2] << ";\n";
                        }
                        out << "}\n";
                    }
                    else if (dynamic_cast<const Continuous2DFunction*>(functions[i]) != NULL) {
                        out << "real x = " << getTempName(node.getChildren()[0], temps) << suffix << ";\n";
                        out << "real y = " << getTempName(node.getChildren()[1], temps) << suffix << ";\n";
                        out << "if (x >= " << paramsFloat[2] << " && x <= " << paramsFloat[3] << " && y >= " << paramsFloat[4] << " && y <= " << paramsFloat[5] << ") {\n";
                        out << "x = (x - " << paramsFloat[2] << ")*" << paramsFloat[6] << ";\n";
                        out << "y = (y - " << paramsFloat[4] << ")*" << paramsFloat[7] << ";\n";
                        out << "int s = min((int) floor(x), " << paramsInt[0] << "-1);\n";
                        out << "int t = min((int) floor(y), " << paramsInt[1] << "-1);\n";
                        out << "int coeffIndex = 4*(s+" << paramsInt[0] << "*t);\n";
                        out << "float4 c[4];\n";
                        for (int j = 0; j < 4; j++)
                            out << "c[" << j << "] = " << functionNames[i].second << "[coeffIndex+" << j << "];\n";
                        out << "real da = x-s;\n";
                        out << "real db = y-t;\n";
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0 && derivOrder[1] == 0) {
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + ((c[3].w*db + c[3].z)*db + c[3].y)*db + c[3].x;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + ((c[2].w*db + c[2].z)*db + c[2].y)*db + c[2].x;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + ((c[1].w*db + c[1].z)*db + c[1].y)*db + c[1].x;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + ((c[0].w*db + c[0].z)*db + c[0].y)*db + c[0].x;\n";
                            }
                            else if (derivOrder[0] == 1 && derivOrder[1] == 0) {
                                out << nodeNames[j] << suffix << " = db*" << nodeNames[j] << suffix << " + (3.0f*c[3].w*da + 2.0f*c[2].w)*da + c[1].w;\n";
                                out << nodeNames[j] << suffix << " = db*" << nodeNames[j] << suffix << " + (3.0f*c[3].z*da + 2.0f*c[2].z)*da + c[1].z;\n";
                                out << nodeNames[j] << suffix << " = db*" << nodeNames[j] << suffix << " + (3.0f*c[3].y*da + 2.0f*c[2].y)*da + c[1].y;\n";
                                out << nodeNames[j] << suffix << " = db*" << nodeNames[j] << suffix << " + (3.0f*c[3].x*da + 2.0f*c[2].x)*da + c[1].x;\n";
                                out << nodeNames[j] << suffix << " *= " << paramsFloat[6] << ";\n";
                            }
                            else if (derivOrder[0] == 0 && derivOrder[1] == 1) {
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + (3.0f*c[3].w*db + 2.0f*c[3].z)*db + c[3].y;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + (3.0f*c[2].w*db + 2.0f*c[2].z)*db + c[2].y;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + (3.0f*c[1].w*db + 2.0f*c[1].z)*db + c[1].y;\n";
                                out << nodeNames[j] << suffix << " = da*" << nodeNames[j] << suffix << " + (3.0f*c[0].w*db + 2.0f*c[0].z)*db + c[0].y;\n";
                                out << nodeNames[j] << suffix << " *= " << paramsFloat[7] << ";\n";
                            }
                            else
                                throw OpenMMException("Unsupported derivative order for Continuous2DFunction");
                        }
                        out << "}\n";
                    }
                    else if (dynamic_cast<const Continuous3DFunction*>(functions[i]) != NULL) {
                        out << "real x = " << getTempName(node.getChildren()[0], temps) << suffix << ";\n";
                        out << "real y = " << getTempName(node.getChildren()[1], temps) << suffix << ";\n";
                        out << "real z = " << getTempName(node.getChildren()[2], temps) << suffix << ";\n";
                        out << "if (x >= " << paramsFloat[3] << " && x <= " << paramsFloat[4] << " && y >= " << paramsFloat[5] << " && y <= " << paramsFloat[6] << " && z >= " << paramsFloat[7] << " && z <= " << paramsFloat[8] << ") {\n";
                        out << "x = (x - " << paramsFloat[3] << ")*" << paramsFloat[9] << ";\n";
                        out << "y = (y - " << paramsFloat[5] << ")*" << paramsFloat[10] << ";\n";
                        out << "z = (z - " << paramsFloat[7] << ")*" << paramsFloat[11] << ";\n";
                        out << "int s = min((int) floor(x), " << paramsInt[0] << "-1);\n";
                        out << "int t = min((int) floor(y), " << paramsInt[1] << "-1);\n";
                        out << "int u = min((int) floor(z), " << paramsInt[2] << "-1);\n";
                        out << "int coeffIndex = 16*(s+" << paramsInt[0] << "*(t+" << paramsInt[1] << "*u));\n";
                        out << "float4 c[16];\n";
                        for (int j = 0; j < 16; j++)
                            out << "c[" << j << "] = " << functionNames[i].second << "[coeffIndex+" << j << "];\n";
                        out << "real da = x-s;\n";
                        out << "real db = y-t;\n";
                        out << "real dc = z-u;\n";
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 0) {
                                out << "real value[4] = {0, 0, 0, 0};\n";
                                for (int k = 3; k >= 0; k--)
                                    for (int m = 0; m < 4; m++) {
                                        int base = k + 4*m;
                                        out << "value[" << m << "] = db*value[" << m << "] + ((c[" << base << "].w*da + c[" << base << "].z)*da + c[" << base << "].y)*da + c[" << base << "].x;\n";
                                    }
                                out << nodeNames[j] << suffix << " = value[0] + dc*(value[1] + dc*(value[2] + dc*value[3]));\n";
                            }
                            else if (derivOrder[0] == 1 && derivOrder[1] == 0 && derivOrder[2] == 0) {
                                out << "real derivx[4] = {0, 0, 0, 0};\n";
                                for (int k = 3; k >= 0; k--)
                                    for (int m = 0; m < 4; m++) {
                                        int base = k + 4*m;
                                        out << "derivx[" << m << "] = db*derivx[" << m << "] + (3*c[" << base << "].w*da + 2*c[" << base << "].z)*da + c[" << base << "].y;\n";
                                    }
                                out << nodeNames[j] << suffix << " = derivx[0] + dc*(derivx[1] + dc*(derivx[2] + dc*derivx[3]));\n";
                                out << nodeNames[j] << suffix << " *= " << paramsFloat[9] << ";\n";
                            }
                            else if (derivOrder[0] == 0 && derivOrder[1] == 1 && derivOrder[2] == 0) {
                                const string suffixes[] = {".x", ".y", ".z", ".w"};
                                out << "real derivy[4] = {0, 0, 0, 0};\n";
                                for (int k = 3; k >= 0; k--)
                                    for (int m = 0; m < 4; m++) {
                                        int base = 4*m;
                                        string suffix = suffixes[k];
                                        out << "derivy[" << m << "] = da*derivy[" << m << "] + (3*c[" << (base+3) << "]" << suffix << "*db + 2*c[" << (base+2) << "]" << suffix << ")*db + c[" << (base+1) << "]" << suffix << ";\n";
                                    }
                                out << nodeNames[j] << suffix << " = derivy[0] + dc*(derivy[1] + dc*(derivy[2] + dc*derivy[3]));\n";
                                out << nodeNames[j] << suffix << " *= " << paramsFloat[10] << ";\n";
                            }
                            else if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 1) {
                                out << "real derivz[4] = {0, 0, 0, 0};\n";
                                for (int k = 3; k >= 0; k--)
                                    for (int m = 0; m < 4; m++) {
                                        int base = k + 4*m;
                                        out << "derivz[" << m << "] = db*derivz[" << m << "] + ((c[" << base << "].w*da + c[" << base << "].z)*da + c[" << base << "].y)*da + c[" << base << "].x;\n";
                                    }
                                out << nodeNames[j] << suffix << " = derivz[1] + dc*(2*derivz[2] + dc*3*derivz[3]);\n";
                                out << nodeNames[j] << suffix << " *= " << paramsFloat[11] << ";\n";
                            }
                            else
                                throw OpenMMException("Unsupported derivative order for Continuous3DFunction");
                        }
                        out << "}\n";
                    }
                    else if (dynamic_cast<const Discrete1DFunction*>(functions[i]) != NULL) {
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0) {
                                out << "real x = " << getTempName(node.getChildren()[0], temps) << suffix << ";\n";
                                out << "if (x >= 0 && x < " << paramsInt[0] << ") {\n";
                                out << "int index = (int) floor(x+0.5f);\n";
                                out << nodeNames[j] << suffix << " = " << functionNames[i].second << "[index];\n";
                                out << "}\n";
                            }
                        }
                    }
                    else if (dynamic_cast<const Discrete2DFunction*>(functions[i]) != NULL) {
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0 && derivOrder[1] == 0) {
                                out << "int x = (int) floor(" << getTempName(node.getChildren()[0], temps) << suffix << "+0.5f);\n";
                                out << "int y = (int) floor(" << getTempName(node.getChildren()[1], temps) << suffix << "+0.5f);\n";
                                out << "int xsize = (int) " << paramsInt[0] << ";\n";
                                out << "int ysize = (int) " << paramsInt[1] << ";\n";
                                out << "int index = x+y*xsize;\n";
                                out << "if (index >= 0 && index < xsize*ysize)\n";
                                out << nodeNames[j] << suffix << " = " << functionNames[i].second << "[index];\n";
                            }
                        }
                    }
                    else if (dynamic_cast<const Discrete3DFunction*>(functions[i]) != NULL) {
                        for (int j = 0; j < nodes.size(); j++) {
                            const vector<int>& derivOrder = dynamic_cast<const Operation::Custom*>(&nodes[j]->getOperation())->getDerivOrder();
                            if (derivOrder[0] == 0 && derivOrder[1] == 0 && derivOrder[2] == 0) {
                                out << "int x = (int) floor(" << getTempName(node.getChildren()[0], temps) << suffix << "+0.5f);\n";
                                out << "int y = (int) floor(" << getTempName(node.getChildren()[1], temps) << suffix << "+0.5f);\n";
                                out << "int z = (int) floor(" << getTempName(node.getChildren()[2], temps) << suffix << "+0.5f);\n";
                                out << "int xsize = (int) " << paramsInt[0] << ";\n";
                                out << "int ysize = (int) " << paramsInt[1] << ";\n";
                                out << "int zsize = (int) " << paramsInt[2] << ";\n";
                                out << "int index = x+(y+z*ysize)*xsize;\n";
                                out << "if (index >= 0 && index < xsize*ysize*zsize)\n";
                                out << nodeNames[j] << suffix << " = " << functionNames[i].second << "[index];\n";
                            }
                        }
                    }
                    out << "}\n";
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
            if (node.getChildren()[1].getOperation().getId() == Operation::RECIPROCAL) {
                for (int i = 0; i < (int) temps.size(); i++)
                    if (temps[i].first == node.getChildren()[1].getChildren()[1]) {
                        haveReciprocal = true;
                        out << getTempName(node.getChildren()[0], temps) << "*" << temps[i].second;
                    }
            }
            if (!haveReciprocal)
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
            out << "pow((" << tempType << ") " << getTempName(node.getChildren()[0], temps) << ", (" << tempType << ") " << getTempName(node.getChildren()[1], temps) << ")";
            break;
        case Operation::NEGATE:
            out << "-" << getTempName(node.getChildren()[0], temps);
            break;
        case Operation::SQRT:
            callFunction(out, "sqrtf", "sqrt", getTempName(node.getChildren()[0], temps), tempType); 
            break;
        case Operation::EXP:
            callFunction(out, "expf", "exp", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::LOG:
            callFunction(out, "logf", "log", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::SIN:
            callFunction(out, "sinf", "sin", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::COS:
            callFunction(out, "cosf", "cos", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::SEC:
            out << "1/";
            callFunction(out, "cosf", "cos", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::CSC:
            out << "1/";
            callFunction(out, "sinf", "sin", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::TAN:
            callFunction(out, "tanf", "tan", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::COT:
            out << "1/";
            callFunction(out, "tanf", "tan", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ASIN:
            callFunction(out, "asinf", "asin", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ACOS:
            callFunction(out, "acosf", "acos", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ATAN:
            callFunction(out, "atanf", "atan", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ATAN2:
            callFunction2(out, "atan2f", "atan2", getTempName(node.getChildren()[0], temps), getTempName(node.getChildren()[1], temps), tempType);
            break;
        case Operation::SINH:
            callFunction(out, "sinh", "sinh", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::COSH:
            callFunction(out, "cosh", "cosh", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::TANH:
            callFunction(out, "tanh", "tanh", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ERF:
            callFunction(out, "erf", "erf", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::ERFC:
            callFunction(out, "erfc", "erfc", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::STEP:
        {
            string compareVal = getTempName(node.getChildren()[0], temps);
            if (isVecType) {
                out << "make_" << tempType << "(0);\n";
                out << "{\n";
                out << tempType<<" tempCompareValue = " << compareVal << ";\n";
                out << name << ".x = (tempCompareValue.x >= 0 ? 1 : 0);\n";
                out << name << ".y = (tempCompareValue.y >= 0 ? 1 : 0);\n";
                out << name << ".z = (tempCompareValue.z >= 0 ? 1 : 0);\n";
                out << "}\n";
            }
            else
                out << compareVal << " >= 0 ? 1 : 0";
            break;
        }
        case Operation::DELTA:
        {
            string compareVal = getTempName(node.getChildren()[0], temps);
            if (isVecType) {
                out << "make_" << tempType << "(0);\n";
                out << "{\n";
                out << tempType<<" tempCompareValue = " << compareVal << ";\n";
                out << name << ".x = (tempCompareValue.x == 0 ? 1 : 0);\n";
                out << name << ".y = (tempCompareValue.y == 0 ? 1 : 0);\n";
                out << name << ".z = (tempCompareValue.z == 0 ? 1 : 0);\n";
                out << "}\n";
            }
            else
                out << compareVal << " == 0 ? 1 : 0";
            break;
        }
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
            if (isVecType) {
                string val = context.doubleToString(dynamic_cast<const Operation::AddConstant*>(&node.getOperation())->getValue());
                string arg = getTempName(node.getChildren()[0], temps);
                out << "make_" << tempType << "(";
                out << val << "+" << arg << ".x, ";
                out << val << "+" << arg << ".y, ";
                out << val << "+" << arg << ".z)";
            }
            else
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
                for (auto& power : powers) {
                    if (power.first != exponent) {
                        exponents.push_back(power.first >= 0 ? power.first : -power.first);
                        string name2 = prefix+context.intToString(temps.size());
                        names.push_back(name2);
                        temps.push_back(make_pair(*power.second, name2));
                        out << tempType << " " << name2 << " = 0.0f;\n";
                    }
                }
                out << "{\n";
                out << "real multiplier = " << (exponent < 0.0 ? "RECIP(" : "(") << getTempName(node.getChildren()[0], temps) << ");\n";
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
                out << "pow((" << tempType << ") " << getTempName(node.getChildren()[0], temps) << ", (" << tempType << ") " << context.doubleToString(exponent) << ")";
            break;
        }
        case Operation::MIN:
            callFunction2(out, "min", "min", getTempName(node.getChildren()[0], temps), getTempName(node.getChildren()[1], temps), tempType);
            break;
        case Operation::MAX:
            callFunction2(out, "max", "max", getTempName(node.getChildren()[0], temps), getTempName(node.getChildren()[1], temps), tempType);
            break;
        case Operation::ABS:
            callFunction(out, "fabs", "fabs", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::FLOOR:
            callFunction(out, "floor", "floor", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::CEIL:
            callFunction(out, "ceil", "ceil", getTempName(node.getChildren()[0], temps), tempType);
            break;
        case Operation::SELECT:
        {
            string compareVal = getTempName(node.getChildren()[0], temps);
            string val1 = getTempName(node.getChildren()[1], temps);
            string val2 = getTempName(node.getChildren()[2], temps);
            if (isVecType) {
                out << "make_" << tempType << "(0);\n";
                out << "{\n";
                out << tempType<<" tempCompareValue = " << compareVal << ";\n";
                out << name << ".x = (tempCompareValue.x != 0 ? " << val1 << ".x : " << val2 << ".x);\n";
                out << name << ".y = (tempCompareValue.y != 0 ? " << val1 << ".y : " << val2 << ".y);\n";
                out << name << ".z = (tempCompareValue.z != 0 ? " << val1 << ".z : " << val2 << ".z);\n";
                out << "}\n";
            }
            else
                out << "(" << compareVal << " != 0 ? " << val1 << " : " << val2 << ")";
            break;
        }
        default:
            throw OpenMMException("Internal error: Unknown operation in user-defined expression: "+node.getOperation().getName());
    }
    out << ";\n";
    if (!hasRecordedNode)
        temps.push_back(make_pair(node, name));
}

string CudaExpressionUtilities::getTempName(const ExpressionTreeNode& node, const vector<pair<ExpressionTreeNode, string> >& temps) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return temps[i].second;
    stringstream out;
    out << "Internal error: No temporary variable for expression node: " << node;
    throw OpenMMException(out.str());
}

void CudaExpressionUtilities::findRelatedCustomFunctions(const ExpressionTreeNode& node, const ExpressionTreeNode& searchNode,
            vector<const Lepton::ExpressionTreeNode*>& nodes) {
    if (searchNode.getOperation().getId() == Operation::CUSTOM && node.getOperation().getName() == searchNode.getOperation().getName()) {
        // Make sure the arguments are identical.
        
        for (int i = 0; i < (int) node.getChildren().size(); i++)
            if (node.getChildren()[i] != searchNode.getChildren()[i])
                return;
        
        // See if we already have an identical node.
        
        for (int i = 0; i < (int) nodes.size(); i++)
            if (*nodes[i] == searchNode)
                return;
        
        // Add the node.
        
        nodes.push_back(&searchNode);
    }
    else
        for (int i = 0; i < (int) searchNode.getChildren().size(); i++)
            findRelatedCustomFunctions(node, searchNode.getChildren()[i], nodes);
}

void CudaExpressionUtilities::findRelatedPowers(const ExpressionTreeNode& node, const ExpressionTreeNode& searchNode, map<int, const ExpressionTreeNode*>& powers) {
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

vector<float> CudaExpressionUtilities::computeFunctionCoefficients(const TabulatedFunction& function, int& width) {
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
    if (dynamic_cast<const Continuous2DFunction*>(&function) != NULL) {
        // Compute the spline coefficients.

        const Continuous2DFunction& fn = dynamic_cast<const Continuous2DFunction&>(function);
        vector<double> values;
        int xsize, ysize;
        double xmin, xmax, ymin, ymax;
        fn.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
        vector<double> x(xsize), y(ysize);
        for (int i = 0; i < xsize; i++)
            x[i] = xmin+i*(xmax-xmin)/(xsize-1);
        for (int i = 0; i < ysize; i++)
            y[i] = ymin+i*(ymax-ymin)/(ysize-1);
        vector<vector<double> > c;
        SplineFitter::create2DNaturalSpline(x, y, values, c);
        vector<float> f(16*c.size());
        for (int i = 0; i < (int) c.size(); i++) {
            for (int j = 0; j < 16; j++)
                f[16*i+j] = (float) c[i][j];
        }
        width = 4;
        return f;
    }
    if (dynamic_cast<const Continuous3DFunction*>(&function) != NULL) {
        // Compute the spline coefficients.

        const Continuous3DFunction& fn = dynamic_cast<const Continuous3DFunction&>(function);
        vector<double> values;
        int xsize, ysize, zsize;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        fn.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
        vector<double> x(xsize), y(ysize), z(zsize);
        for (int i = 0; i < xsize; i++)
            x[i] = xmin+i*(xmax-xmin)/(xsize-1);
        for (int i = 0; i < ysize; i++)
            y[i] = ymin+i*(ymax-ymin)/(ysize-1);
        for (int i = 0; i < zsize; i++)
            z[i] = zmin+i*(zmax-zmin)/(zsize-1);
        vector<vector<double> > c;
        SplineFitter::create3DNaturalSpline(x, y, z, values, c);
        vector<float> f(64*c.size());
        for (int i = 0; i < (int) c.size(); i++) {
            for (int j = 0; j < 64; j++)
                f[64*i+j] = (float) c[i][j];
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

vector<vector<double> > CudaExpressionUtilities::computeFunctionParameters(const vector<const TabulatedFunction*>& functions) {
    vector<vector<double> > params(functions.size());
    for (int i = 0; i < (int) functions.size(); i++) {
        if (dynamic_cast<const Continuous1DFunction*>(functions[i]) != NULL) {
            const Continuous1DFunction& fn = dynamic_cast<const Continuous1DFunction&>(*functions[i]);
            vector<double> values;
            double min, max;
            fn.getFunctionParameters(values, min, max);
            params[i].push_back(min);
            params[i].push_back(max);
            params[i].push_back((values.size()-1)/(max-min));
            params[i].push_back(values.size()-2);
        }
        else if (dynamic_cast<const Continuous2DFunction*>(functions[i]) != NULL) {
            const Continuous2DFunction& fn = dynamic_cast<const Continuous2DFunction&>(*functions[i]);
            vector<double> values;
            int xsize, ysize;
            double xmin, xmax, ymin, ymax;
            fn.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
            params[i].push_back(xsize-1);
            params[i].push_back(ysize-1);
            params[i].push_back(xmin);
            params[i].push_back(xmax);
            params[i].push_back(ymin);
            params[i].push_back(ymax);
            params[i].push_back((xsize-1)/(xmax-xmin));
            params[i].push_back((ysize-1)/(ymax-ymin));
        }
        else if (dynamic_cast<const Continuous3DFunction*>(functions[i]) != NULL) {
            const Continuous3DFunction& fn = dynamic_cast<const Continuous3DFunction&>(*functions[i]);
            vector<double> values;
            int xsize, ysize, zsize;
            double xmin, xmax, ymin, ymax, zmin, zmax;
            fn.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
            params[i].push_back(xsize-1);
            params[i].push_back(ysize-1);
            params[i].push_back(zsize-1);
            params[i].push_back(xmin);
            params[i].push_back(xmax);
            params[i].push_back(ymin);
            params[i].push_back(ymax);
            params[i].push_back(zmin);
            params[i].push_back(zmax);
            params[i].push_back((xsize-1)/(xmax-xmin));
            params[i].push_back((ysize-1)/(ymax-ymin));
            params[i].push_back((zsize-1)/(zmax-zmin));
        }
        else if (dynamic_cast<const Discrete1DFunction*>(functions[i]) != NULL) {
            const Discrete1DFunction& fn = dynamic_cast<const Discrete1DFunction&>(*functions[i]);
            vector<double> values;
            fn.getFunctionParameters(values);
            params[i].push_back(values.size());
        }
        else if (dynamic_cast<const Discrete2DFunction*>(functions[i]) != NULL) {
            const Discrete2DFunction& fn = dynamic_cast<const Discrete2DFunction&>(*functions[i]);
            int xsize, ysize;
            vector<double> values;
            fn.getFunctionParameters(xsize, ysize, values);
            params[i].push_back(xsize);
            params[i].push_back(ysize);
        }
        else if (dynamic_cast<const Discrete3DFunction*>(functions[i]) != NULL) {
            const Discrete3DFunction& fn = dynamic_cast<const Discrete3DFunction&>(*functions[i]);
            int xsize, ysize, zsize;
            vector<double> values;
            fn.getFunctionParameters(xsize, ysize, zsize, values);
            params[i].push_back(xsize);
            params[i].push_back(ysize);
            params[i].push_back(zsize);
        }
        else
            throw OpenMMException("computeFunctionParameters: Unknown function type");
    }
    return params;
}

Lepton::CustomFunction* CudaExpressionUtilities::getFunctionPlaceholder(const TabulatedFunction& function) {
    if (dynamic_cast<const Continuous1DFunction*>(&function) != NULL)
        return &fp1;
    if (dynamic_cast<const Continuous2DFunction*>(&function) != NULL)
        return &fp2;
    if (dynamic_cast<const Continuous3DFunction*>(&function) != NULL)
        return &fp3;
    if (dynamic_cast<const Discrete1DFunction*>(&function) != NULL)
        return &fp1;
    if (dynamic_cast<const Discrete2DFunction*>(&function) != NULL)
        return &fp2;
    if (dynamic_cast<const Discrete3DFunction*>(&function) != NULL)
        return &fp3;
    throw OpenMMException("getFunctionPlaceholder: Unknown function type");
}

Lepton::CustomFunction* CudaExpressionUtilities::getPeriodicDistancePlaceholder() {
    return &periodicDistance;
}

void CudaExpressionUtilities::callFunction(stringstream& out, string singleFn, string doubleFn, const string& arg, const string& tempType) {
    bool isDouble = (tempType[0] == 'd');
    bool isVector = (tempType[tempType.size()-1] == '3');
    string fn = (isDouble ? doubleFn : singleFn);
    if (isVector)
        out<<"make_"<<tempType<<"("<<fn<<"("<<arg<<".x), "<<fn<<"("<<arg<<".y), "<<fn<<"("<<arg<<".z))";
    else
        out<<fn<<"("<<arg<<")";
}

void CudaExpressionUtilities::callFunction2(stringstream& out, string singleFn, string doubleFn, const string& arg1, const string& arg2, const string& tempType) {
    bool isDouble = (tempType[0] == 'd');
    bool isVector = (tempType[tempType.size()-1] == '3');
    string fn = (isDouble ? doubleFn : singleFn);
    if (isVector) {
        out<<"make_"<<tempType<<"(";
        out<<fn<<"("<<arg1<<".x, "<<arg2<<".x), ";
        out<<fn<<"("<<arg1<<".y, "<<arg2<<".y), ";
        out<<fn<<"("<<arg1<<".z, "<<arg2<<".z))";
    }
    else
        out<<fn<<"(("<<tempType<<") "<<arg1<<", ("<<tempType<<") "<<arg2<<")";
}
