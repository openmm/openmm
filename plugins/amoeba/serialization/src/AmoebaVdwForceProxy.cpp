/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2020 Stanford University and the Authors.      *
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

#include "openmm/serialization/AmoebaVdwForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaVdwForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaVdwForceProxy::AmoebaVdwForceProxy() : SerializationProxy("AmoebaVdwForce") {
}

void AmoebaVdwForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 4);
    const AmoebaVdwForce& force = *reinterpret_cast<const AmoebaVdwForce*>(object);
    bool useTypes = force.getUseParticleTypes();

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setStringProperty("SigmaCombiningRule", force.getSigmaCombiningRule());
    node.setStringProperty("EpsilonCombiningRule", force.getEpsilonCombiningRule());
    node.setDoubleProperty("VdwCutoff", force.getCutoffDistance());
    node.setIntProperty("method", (int) force.getNonbondedMethod());

    node.setDoubleProperty("n", force.getSoftcorePower());
    node.setDoubleProperty("alpha", force.getSoftcoreAlpha());
    node.setIntProperty("alchemicalMethod", (int) force.getAlchemicalMethod());
    node.setIntProperty("potentialFunction", (int) force.getPotentialFunction());
    node.setBoolProperty("useTypes", useTypes);

    SerializationNode& particles = node.createChildNode("VdwParticles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int ivIndex, typeIndex;
        double sigma, epsilon, reductionFactor;
        bool isAlchemical;
        force.getParticleParameters(i, ivIndex, sigma, epsilon, reductionFactor, isAlchemical, typeIndex);
        SerializationNode& particle = particles.createChildNode("Particle");
        if (useTypes)
            particle.setIntProperty("ivIndex", ivIndex).setIntProperty("type", typeIndex).setDoubleProperty("reductionFactor", reductionFactor).setBoolProperty("isAlchemical", isAlchemical);
        else
            particle.setIntProperty("ivIndex", ivIndex).setDoubleProperty("sigma", sigma).setDoubleProperty("epsilon", epsilon).setDoubleProperty("reductionFactor", reductionFactor).setBoolProperty("isAlchemical", isAlchemical); 

        std::vector< int > exclusions;
        force.getParticleExclusions(i,  exclusions);
        SerializationNode& particleExclusions = particle.createChildNode("ParticleExclusions");
        for (int j = 0; j < exclusions.size(); j++)
            particleExclusions.createChildNode("excl").setIntProperty("index", exclusions[j]);
    }
    if (useTypes) {
        SerializationNode& types = node.createChildNode("ParticleTypes");
        for (int i = 0; i < force.getNumParticleTypes(); i++) {
            double sigma, epsilon;
            force.getParticleTypeParameters(i, sigma, epsilon);
            SerializationNode& type = types.createChildNode("Type");
            type.setDoubleProperty("sigma", sigma).setDoubleProperty("epsilon", epsilon);
        }
        SerializationNode& pairs = node.createChildNode("TypePairs");
        for (int i = 0; i < force.getNumTypePairs(); i++) {
            int type1, type2;
            double sigma, epsilon;
            force.getTypePairParameters(i, type1, type2, sigma, epsilon);
            SerializationNode& pair = pairs.createChildNode("Pair");
            pair.setIntProperty("type1", type1).setIntProperty("type2", type2).setDoubleProperty("sigma", sigma).setDoubleProperty("epsilon", epsilon);
        }
    }
}

void* AmoebaVdwForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 4)
        throw OpenMMException("Unsupported version number");
    AmoebaVdwForce* force = new AmoebaVdwForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        force->setSigmaCombiningRule(node.getStringProperty("SigmaCombiningRule"));
        force->setEpsilonCombiningRule(node.getStringProperty("EpsilonCombiningRule"));
        force->setCutoffDistance(node.getDoubleProperty("VdwCutoff"));
        force->setNonbondedMethod((AmoebaVdwForce::NonbondedMethod) node.getIntProperty("method"));

        if (version > 2) {
           force->setAlchemicalMethod((AmoebaVdwForce::AlchemicalMethod) node.getIntProperty("alchemicalMethod"));
           force->setSoftcorePower(node.getDoubleProperty("n"));
           force->setSoftcoreAlpha(node.getDoubleProperty("alpha"));
        }
        
        bool useTypes = false;
        if (version > 3) {
            useTypes = node.getBoolProperty("useTypes");
           force->setPotentialFunction((AmoebaVdwForce::PotentialFunction) node.getIntProperty("potentialFunction"));
        }

        const SerializationNode& particles = node.getChildNode("VdwParticles");
        for (int i = 0; i < particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];

            if (version < 3) 
               force->addParticle(particle.getIntProperty("ivIndex"), particle.getDoubleProperty("sigma"), 
                                  particle.getDoubleProperty("epsilon"), particle.getDoubleProperty("reductionFactor"));
            else if (useTypes)
               force->addParticle(particle.getIntProperty("ivIndex"), particle.getIntProperty("type"), 
                                  particle.getDoubleProperty("reductionFactor"), particle.getBoolProperty("isAlchemical"));
            else
               force->addParticle(particle.getIntProperty("ivIndex"), particle.getDoubleProperty("sigma"), 
                                  particle.getDoubleProperty("epsilon"), particle.getDoubleProperty("reductionFactor"), 
                                  particle.getBoolProperty("isAlchemical"));

            // exclusions

            const SerializationNode& particleExclusions = particle.getChildNode("ParticleExclusions");
            std::vector<int> exclusions;
            for (int j = 0; j < particleExclusions.getChildren().size(); j++)
                exclusions.push_back(particleExclusions.getChildren()[j].getIntProperty("index"));
            force->setParticleExclusions(i, exclusions);
        }
        if (useTypes) {
            const SerializationNode& types = node.getChildNode("ParticleTypes");
            for (const auto& type : types.getChildren())
                force->addParticleType(type.getDoubleProperty("sigma"), type.getDoubleProperty("epsilon"));
            const SerializationNode& pairs = node.getChildNode("TypePairs");
            for (const auto& pair : pairs.getChildren())
                force->addTypePair(pair.getIntProperty("type1"), pair.getIntProperty("type2"), pair.getDoubleProperty("sigma"), pair.getDoubleProperty("epsilon"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
