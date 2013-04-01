/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Yutong Zhao                                        *
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

#include "openmm/serialization/CustomIntegratorProxy.h"
#include <OpenMM.h>

using namespace std;
using namespace OpenMM;

CustomIntegratorProxy::CustomIntegratorProxy() : SerializationProxy("CustomIntegrator") {

}

void CustomIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomIntegrator& integrator = *reinterpret_cast<const CustomIntegrator*>(object);
    SerializationNode& globalVariablesNode = node.createChildNode("GlobalVariables");
    for(int i=0; i<integrator.getNumGlobalVariables(); i++) {
        globalVariablesNode.setDoubleProperty(integrator.getGlobalVariableName(i), integrator.getGlobalVariable(i));
    }
    SerializationNode& perDofVariablesNode = node.createChildNode("PerDofVariables");
    for(int i=0; i<integrator.getNumPerDofVariables(); i++) {
        SerializationNode& perDofValuesNode = perDofVariablesNode.createChildNode(integrator.getPerDofVariableName(i));
        vector<Vec3> perDofValues; integrator.getPerDofVariable(i, perDofValues);
        for(int j=0; j<perDofValues.size(); j++) {
            perDofValuesNode.createChildNode("Value").setDoubleProperty("x",perDofValues[j][0]).setDoubleProperty("y",perDofValues[j][1]).setDoubleProperty("z",perDofValues[j][2]);
        }
    }
    SerializationNode& computationsNode = node.createChildNode("Computations");
    for(int i=0; i<integrator.getNumComputations(); i++) {
        CustomIntegrator::ComputationType computationType;
        string computationVariable;
        string computationExpression;
        integrator.getComputationStep(i, computationType, computationVariable, computationExpression);
        computationsNode.createChildNode("Computation").setIntProperty("computationType",static_cast<int>(computationType))
            .setStringProperty("computationVariable",computationVariable).setStringProperty("computationExpression",computationExpression);
    }
    node.setStringProperty("kineticEnergyExpression",integrator.getKineticEnergyExpression());
    node.setIntProperty("randomSeed",integrator.getRandomNumberSeed());
    node.setDoubleProperty("stepSize",integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance",integrator.getConstraintTolerance());
}

void* CustomIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomIntegrator* integrator = new CustomIntegrator(node.getDoubleProperty("stepSize"));
    const SerializationNode& globalVariablesNode = node.getChildNode("GlobalVariables");
    const map<string, string> &globalVariableProp = globalVariablesNode.getProperties();
    for(map<string, string>::const_iterator cit = globalVariableProp.begin(); cit != globalVariableProp.end(); cit++) {
        integrator->addGlobalVariable(cit->first, globalVariablesNode.getDoubleProperty(cit->first));
    }
    const SerializationNode& perDofVariablesNode = node.getChildNode("PerDofVariables");
    const vector<SerializationNode>& perDofVariableList = perDofVariablesNode.getChildren();
    int count = 0;
    for(vector<SerializationNode>::const_iterator cit=perDofVariableList.begin(); cit != perDofVariableList.end(); cit++, count++) {
        const vector<SerializationNode>& perDofVariableVector = cit->getChildren();
        integrator->addPerDofVariable(cit->getName(),0);
        vector<Vec3> perDofValues;
        for(vector<SerializationNode>::const_iterator dit=perDofVariableVector.begin(); dit!=perDofVariableVector.end(); dit++) {
            perDofValues.push_back(Vec3(dit->getDoubleProperty("x"),dit->getDoubleProperty("y"),dit->getDoubleProperty("z")));
        }
        integrator->setPerDofVariable(count, perDofValues);
    }
    const SerializationNode& computationsNode = node.getChildNode("Computations");
    const vector<SerializationNode>& computationsList = computationsNode.getChildren();
    for(vector<SerializationNode>::const_iterator cit = computationsList.begin(); cit != computationsList.end(); cit++) {
        CustomIntegrator::ComputationType computationType = static_cast<CustomIntegrator::ComputationType>(cit->getIntProperty("computationType"));
        // make sure that the int casts to a valid enum
        if(computationType == CustomIntegrator::ComputeGlobal) {
            integrator->addComputeGlobal(cit->getStringProperty("computationVariable"), cit->getStringProperty("computationExpression"));
        } else if(computationType == CustomIntegrator::ComputePerDof) {
            integrator->addComputePerDof(cit->getStringProperty("computationVariable"), cit->getStringProperty("computationExpression"));
        } else if(computationType == CustomIntegrator::ComputeSum) {
            integrator->addComputeSum(cit->getStringProperty("computationVariable"), cit->getStringProperty("computationExpression"));
        } else if(computationType == CustomIntegrator::ConstrainPositions) {
            integrator->addConstrainPositions();
        } else if(computationType == CustomIntegrator::ConstrainVelocities) {
            integrator->addConstrainVelocities();
        } else if(computationType == CustomIntegrator::UpdateContextState) {
            integrator->addUpdateContextState();
        } else {
            throw(OpenMMException("Custom Integrator Deserialization: Unknown computation type"));
        }
    }
    integrator->setKineticEnergyExpression(node.getStringProperty("kineticEnergyExpression"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    return integrator;
}