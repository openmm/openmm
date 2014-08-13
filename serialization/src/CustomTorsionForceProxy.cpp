/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
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

#include "openmm/serialization/CustomTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomTorsionForceProxy::CustomTorsionForceProxy() : SerializationProxy("CustomTorsionForce") {
}

void CustomTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomTorsionForce& force = *reinterpret_cast<const CustomTorsionForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("energy", force.getEnergyFunction());
    SerializationNode& perTorsionParams = node.createChildNode("PerTorsionParameters");
    for (int i = 0; i < force.getNumPerTorsionParameters(); i++) {
        perTorsionParams.createChildNode("Parameter").setStringProperty("name", force.getPerTorsionParameterName(i));
    }
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& torsions = node.createChildNode("Torsions");
    for (int i = 0; i < force.getNumTorsions(); i++) {
        int p1, p2, p3, p4;
        vector<double> params;
        force.getTorsionParameters(i, p1, p2, p3, p4, params);
        SerializationNode& node = torsions.createChildNode("Torsion").setIntProperty("p1", p1).setIntProperty("p2", p2).setIntProperty("p3", p3).setIntProperty("p4", p4);
        for (int j = 0; j < (int) params.size(); j++) {
            stringstream key;
            key << "param";
            key << j+1;
            node.setDoubleProperty(key.str(), params[j]);
        }
    }
}

void* CustomTorsionForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomTorsionForce* force = NULL;
    try {
        CustomTorsionForce* force = new CustomTorsionForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        const SerializationNode& perTorsionParams = node.getChildNode("PerTorsionParameters");
        for (int i = 0; i < (int) perTorsionParams.getChildren().size(); i++) {
            const SerializationNode& parameter = perTorsionParams.getChildren()[i];
            force->addPerTorsionParameter(parameter.getStringProperty("name"));
        }
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (int i = 0; i < (int) globalParams.getChildren().size(); i++) {
            const SerializationNode& parameter = globalParams.getChildren()[i];
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        }
        const SerializationNode& torsions = node.getChildNode("Torsions");
        vector<double> params(force->getNumPerTorsionParameters());
        for (int i = 0; i < (int) torsions.getChildren().size(); i++) {
            const SerializationNode& torsion = torsions.getChildren()[i];
            for (int j = 0; j < (int) params.size(); j++) {
                stringstream key;
                key << "param";
                key << j+1;
                params[j] = torsion.getDoubleProperty(key.str());
            }
            force->addTorsion(torsion.getIntProperty("p1"), torsion.getIntProperty("p2"), torsion.getIntProperty("p3"), torsion.getIntProperty("p4"), params);
        }
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}
