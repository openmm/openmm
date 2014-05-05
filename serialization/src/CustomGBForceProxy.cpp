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

#include "openmm/serialization/CustomGBForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomGBForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomGBForceProxy::CustomGBForceProxy() : SerializationProxy("CustomGBForce") {
}

void CustomGBForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomGBForce& force = *reinterpret_cast<const CustomGBForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    SerializationNode& perParticleParams = node.createChildNode("PerParticleParameters");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        perParticleParams.createChildNode("Parameter").setStringProperty("name", force.getPerParticleParameterName(i));
    }
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& computedValues = node.createChildNode("ComputedValues");
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(i, name, expression, type);
        computedValues.createChildNode("Value").setStringProperty("name", name).setStringProperty("expression", expression).setIntProperty("type", (int) type);
    }
    SerializationNode& energyTerms = node.createChildNode("EnergyTerms");
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        energyTerms.createChildNode("Term").setStringProperty("expression", expression).setIntProperty("type", (int) type);
    }
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        vector<double> params;
        force.getParticleParameters(i, params);
        SerializationNode& node = particles.createChildNode("Particle");
        for (int j = 0; j < (int) params.size(); j++) {
            stringstream key;
            key << "param";
            key << j+1;
            node.setDoubleProperty(key.str(), params[j]);
        }
    }
    SerializationNode& exclusions = node.createChildNode("Exclusions");
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions.createChildNode("Exclusion").setIntProperty("p1", particle1).setIntProperty("p2", particle2);
    }
    SerializationNode& functions = node.createChildNode("Functions");
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions.createChildNode("Function", &force.getTabulatedFunction(i)).setStringProperty("name", force.getTabulatedFunctionName(i));
}

void* CustomGBForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomGBForce* force = NULL;
    try {
        CustomGBForce* force = new CustomGBForce();
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((CustomGBForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        const SerializationNode& perParticleParams = node.getChildNode("PerParticleParameters");
        for (int i = 0; i < (int) perParticleParams.getChildren().size(); i++) {
            const SerializationNode& parameter = perParticleParams.getChildren()[i];
            force->addPerParticleParameter(parameter.getStringProperty("name"));
        }
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (int i = 0; i < (int) globalParams.getChildren().size(); i++) {
            const SerializationNode& parameter = globalParams.getChildren()[i];
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        }
        const SerializationNode& computedValues = node.getChildNode("ComputedValues");
        for (int i = 0; i < (int) computedValues.getChildren().size(); i++) {
            const SerializationNode& value = computedValues.getChildren()[i];
            force->addComputedValue(value.getStringProperty("name"), value.getStringProperty("expression"), (CustomGBForce::ComputationType) value.getIntProperty("type"));
        }
        const SerializationNode& energyTerms = node.getChildNode("EnergyTerms");
        for (int i = 0; i < (int) energyTerms.getChildren().size(); i++) {
            const SerializationNode& term = energyTerms.getChildren()[i];
            force->addEnergyTerm(term.getStringProperty("expression"), (CustomGBForce::ComputationType) term.getIntProperty("type"));
        }
        const SerializationNode& particles = node.getChildNode("Particles");
        vector<double> params(force->getNumPerParticleParameters());
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            for (int j = 0; j < (int) params.size(); j++) {
                stringstream key;
                key << "param";
                key << j+1;
                params[j] = particle.getDoubleProperty(key.str());
            }
            force->addParticle(params);
        }
        const SerializationNode& exclusions = node.getChildNode("Exclusions");
        for (int i = 0; i < (int) exclusions.getChildren().size(); i++) {
            const SerializationNode& exclusion = exclusions.getChildren()[i];
            force->addExclusion(exclusion.getIntProperty("p1"), exclusion.getIntProperty("p2"));
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
