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

#include "openmm/serialization/CustomNonbondedForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomNonbondedForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomNonbondedForceProxy::CustomNonbondedForceProxy() : SerializationProxy("CustomNonbondedForce") {
}

void CustomNonbondedForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomNonbondedForce& force = *reinterpret_cast<const CustomNonbondedForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("energy", force.getEnergyFunction());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setBoolProperty("useSwitchingFunction", force.getUseSwitchingFunction());
    node.setDoubleProperty("switchingDistance", force.getSwitchingDistance());
    node.setBoolProperty("useLongRangeCorrection", force.getUseLongRangeCorrection());
    SerializationNode& perParticleParams = node.createChildNode("PerParticleParameters");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        perParticleParams.createChildNode("Parameter").setStringProperty("name", force.getPerParticleParameterName(i));
    }
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
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

    SerializationNode& interactionGroups = node.createChildNode("InteractionGroups");
    for (int i = 0; i < force.getNumInteractionGroups(); i++) {
        SerializationNode& interactionGroup = interactionGroups.createChildNode("InteractionGroup");
        std::set<int> set1;
        std::set<int> set2;
        force.getInteractionGroupParameters(i, set1, set2);
        SerializationNode& set1node = interactionGroup.createChildNode("Set1");
        for (std::set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
            set1node.createChildNode("Particle").setIntProperty("index", *it);
        SerializationNode& set2node = interactionGroup.createChildNode("Set2");
        for (std::set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
            set2node.createChildNode("Particle").setIntProperty("index", *it);
    }
}

void* CustomNonbondedForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomNonbondedForce* force = NULL;
    try {
        CustomNonbondedForce* force = new CustomNonbondedForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((CustomNonbondedForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setUseSwitchingFunction(node.getBoolProperty("useSwitchingFunction", false));
        force->setSwitchingDistance(node.getDoubleProperty("switchingDistance", -1.0));
        force->setUseLongRangeCorrection(node.getBoolProperty("useLongRangeCorrection", false));
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
        bool hasInteractionGroups = false; // Older files will be missing this block.
        for (int i = 0; i < (int) node.getChildren().size(); i++) {
            if (node.getChildren()[i].getName() == "InteractionGroups")
                hasInteractionGroups = true;
        }
        if (hasInteractionGroups) {
            const SerializationNode& interactionGroups = node.getChildNode("InteractionGroups");
            for (int i = 0; i < (int) interactionGroups.getChildren().size(); i++) {
                const SerializationNode& interactionGroup = interactionGroups.getChildren()[i];
                // Get set 1.
                const SerializationNode& set1node = interactionGroup.getChildNode("Set1");
                std::set<int> set1;
                for (int j = 0; j < (int) set1node.getChildren().size(); j++)
                    set1.insert(set1node.getChildren()[j].getIntProperty("index"));
                // Get set 2.
                const SerializationNode& set2node = interactionGroup.getChildNode("Set2");
                std::set<int> set2;
                for (int j = 0; j < (int) set2node.getChildren().size(); j++)
                    set2.insert(set2node.getChildren()[j].getIntProperty("index"));
                force->addInteractionGroup(set1, set2);
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
