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

#include "openmm/serialization/CustomManyParticleForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/CustomManyParticleForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

CustomManyParticleForceProxy::CustomManyParticleForceProxy() : SerializationProxy("CustomManyParticleForce") {
}

void CustomManyParticleForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const CustomManyParticleForce& force = *reinterpret_cast<const CustomManyParticleForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("particlesPerSet", force.getNumParticlesPerSet());
    node.setStringProperty("energy", force.getEnergyFunction());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setIntProperty("permutationMode", (int) force.getPermutationMode());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
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
        int type;
        force.getParticleParameters(i, params, type);
        SerializationNode& node = particles.createChildNode("Particle");
        node.setIntProperty("type", type);
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
    SerializationNode& filters = node.createChildNode("TypeFilters");
    for (int i = 0; i < force.getNumParticlesPerSet(); i++) {
        set<int> types;
        force.getTypeFilter(i, types);
        stringstream list;
        bool first = true;
        for (int type : types) {
            if (!first)
                list << ",";
            list << type;
            first = false;
        }
        filters.createChildNode("Filter").setIntProperty("index", i).setStringProperty("types", list.str());
    }
    SerializationNode& functions = node.createChildNode("Functions");
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions.createChildNode("Function", &force.getTabulatedFunction(i)).setStringProperty("name", force.getTabulatedFunctionName(i));
}

void* CustomManyParticleForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    CustomManyParticleForce* force = NULL;
    try {
        CustomManyParticleForce* force = new CustomManyParticleForce(node.getIntProperty("particlesPerSet"), node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((CustomManyParticleForce::NonbondedMethod) node.getIntProperty("method"));
        force->setPermutationMode((CustomManyParticleForce::PermutationMode) node.getIntProperty("permutationMode"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        const SerializationNode& perParticleParams = node.getChildNode("PerParticleParameters");
        for (auto& parameter : perParticleParams.getChildren())
            force->addPerParticleParameter(parameter.getStringProperty("name"));
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        const SerializationNode& particles = node.getChildNode("Particles");
        vector<double> params(force->getNumPerParticleParameters());
        for (auto& particle : particles.getChildren()) {
            for (int j = 0; j < (int) params.size(); j++) {
                stringstream key;
                key << "param";
                key << j+1;
                params[j] = particle.getDoubleProperty(key.str());
            }
            force->addParticle(params, particle.getIntProperty("type"));
        }
        const SerializationNode& exclusions = node.getChildNode("Exclusions");
        for (auto& exclusion : exclusions.getChildren())
            force->addExclusion(exclusion.getIntProperty("p1"), exclusion.getIntProperty("p2"));
        const SerializationNode& filters = node.getChildNode("TypeFilters");
        for (auto& filter : filters.getChildren()) {
            string typesString = filter.getStringProperty("types");
            vector<string> splitTypes;
            size_t searchPos = 0, nextPos;
            while ((nextPos = typesString.find_first_of(", ", searchPos)) != string::npos) {
                splitTypes.push_back(typesString.substr(searchPos, nextPos-searchPos));
                searchPos = nextPos+1;
            }
            splitTypes.push_back(typesString.substr(searchPos));
            set<int> types;
            for (auto& t : splitTypes) {
                if (t.size() > 0) {
                    int type;
                    stringstream(t) >> type;
                    types.insert(type);
                }
            }
            force->setTypeFilter(filter.getIntProperty("index"), types);
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
