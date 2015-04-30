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

#include "openmm/serialization/CustomCompoundBondForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomCompoundBondForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomCompoundBondForceProxy::CustomCompoundBondForceProxy() : SerializationProxy("CustomCompoundBondForce") {
}

void CustomCompoundBondForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomCompoundBondForce& force = *reinterpret_cast<const CustomCompoundBondForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("particles", force.getNumParticlesPerBond());
    node.setStringProperty("energy", force.getEnergyFunction());
    SerializationNode& perBondParams = node.createChildNode("PerBondParameters");
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        perBondParams.createChildNode("Parameter").setStringProperty("name", force.getPerBondParameterName(i));
    }
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& bonds = node.createChildNode("Bonds");
    for (int i = 0; i < force.getNumBonds(); i++) {
        vector<int> particles;
        vector<double> params;
        force.getBondParameters(i, particles, params);
        SerializationNode& node = bonds.createChildNode("Bond");
        for (int j = 0; j < (int) particles.size(); j++) {
            stringstream key;
            key << "p";
            key << j+1;
            node.setIntProperty(key.str(), particles[j]);
        }
        for (int j = 0; j < (int) params.size(); j++) {
            stringstream key;
            key << "param";
            key << j+1;
            node.setDoubleProperty(key.str(), params[j]);
        }
    }
    SerializationNode& functions = node.createChildNode("Functions");
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions.createChildNode("Function", &force.getTabulatedFunction(i)).setStringProperty("name", force.getTabulatedFunctionName(i));
}

void* CustomCompoundBondForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomCompoundBondForce* force = NULL;
    try {
        CustomCompoundBondForce* force = new CustomCompoundBondForce(node.getIntProperty("particles"), node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        const SerializationNode& perBondParams = node.getChildNode("PerBondParameters");
        for (int i = 0; i < (int) perBondParams.getChildren().size(); i++) {
            const SerializationNode& parameter = perBondParams.getChildren()[i];
            force->addPerBondParameter(parameter.getStringProperty("name"));
        }
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (int i = 0; i < (int) globalParams.getChildren().size(); i++) {
            const SerializationNode& parameter = globalParams.getChildren()[i];
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        }
        const SerializationNode& bonds = node.getChildNode("Bonds");
        vector<int> particles(force->getNumParticlesPerBond());
        vector<double> params(force->getNumPerBondParameters());
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& bond = bonds.getChildren()[i];
            for (int j = 0; j < (int) particles.size(); j++) {
                stringstream key;
                key << "p";
                key << j+1;
                particles[j] = bond.getIntProperty(key.str());
            }
            for (int j = 0; j < (int) params.size(); j++) {
                stringstream key;
                key << "param";
                key << j+1;
                params[j] = bond.getDoubleProperty(key.str());
            }
            force->addBond(particles, params);
        }
        const SerializationNode& functions = node.getChildNode("Functions");
        for (int i = 0; i < (int) functions.getChildren().size(); i++) {
            const SerializationNode& function = functions.getChildren()[i];
            if (function.hasProperty("type")) {
                force->addTabulatedFunction(function.getStringProperty("name"), function.decodeObject<TabulatedFunction>());
            }
            else {
                // This is an old file created before TabulatedFunction existed.
                
                const SerializationNode& valuesNode = function.getChildNode("Values");
                vector<double> values;
                for (int j = 0; j < (int) valuesNode.getChildren().size(); j++)
                    values.push_back(valuesNode.getChildren()[j].getDoubleProperty("v"));
                force->addTabulatedFunction(function.getStringProperty("name"), new Continuous1DFunction(values, function.getDoubleProperty("min"), function.getDoubleProperty("max")));
            }
        }
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}
