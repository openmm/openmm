/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2017 Stanford University and the Authors.      *
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

#include "openmm/serialization/CustomCVForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomCVForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomCVForceProxy::CustomCVForceProxy() : SerializationProxy("CustomCVForce") {
}

void CustomCVForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 0);
    const CustomCVForce& force = *reinterpret_cast<const CustomCVForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("energy", force.getEnergyFunction());
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& energyDerivs = node.createChildNode("EnergyParameterDerivatives");
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        energyDerivs.createChildNode("Parameter").setStringProperty("name", force.getEnergyParameterDerivativeName(i));
    }
    SerializationNode& cvs = node.createChildNode("CollectiveVariables");
    for (int i = 0; i < force.getNumCollectiveVariables(); i++) {
        SerializationNode& node = cvs.createChildNode("CollectiveVariable").setStringProperty("name", force.getCollectiveVariableName(i));
        node.createChildNode("Force", &force.getCollectiveVariable(i));
    }
    SerializationNode& functions = node.createChildNode("Functions");
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions.createChildNode("Function", &force.getTabulatedFunction(i)).setStringProperty("name", force.getTabulatedFunctionName(i));
}

void* CustomCVForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 0)
        throw OpenMMException("Unsupported version number");
    CustomCVForce* force = NULL;
    try {
        CustomCVForce* force = new CustomCVForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        const SerializationNode& energyDerivs = node.getChildNode("EnergyParameterDerivatives");
        for (auto& parameter : energyDerivs.getChildren())
            force->addEnergyParameterDerivative(parameter.getStringProperty("name"));
        const SerializationNode& cvs = node.getChildNode("CollectiveVariables");
        for (auto& cv : cvs.getChildren()) {
            string name = cv.getStringProperty("name");
            force->addCollectiveVariable(name, cv.getChildren()[0].decodeObject<Force>());
        }
        const SerializationNode& functions = node.getChildNode("Functions");
        for (auto& function : functions.getChildren())
            force->addTabulatedFunction(function.getStringProperty("name"), function.decodeObject<TabulatedFunction>());
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}
