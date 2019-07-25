/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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

#include "openmm/serialization/AmoebaCTForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaCTForce.h"
#include <sstream>
#include <cmath>

using namespace OpenMM;
using namespace std;

AmoebaCTForceProxy::AmoebaCTForceProxy()
    : SerializationProxy("AmoebaCTForce") {}

void AmoebaCTForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const AmoebaCTForce& force = *reinterpret_cast<const AmoebaCTForce*>(object);

    node.setIntProperty("numCTprTypes", force.getNumCTprTypes());
    node.setStringProperty("ApreCombiningRule", force.getApreCombiningRule());
    node.setStringProperty("BexpCombiningRule", force.getBexpCombiningRule());
    node.setDoubleProperty("CTCutoff", force.getCutoffDistance());

    node.setIntProperty("method", (int) force.getNonbondedMethod());

    SerializationNode& particles = node.createChildNode("CTParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumParticles()); ii++) {

        int CTType;
        double apre, bexp;
        force.getParticleParameters(ii, CTType, apre, bexp);

        SerializationNode& particle = particles.createChildNode("Particle");
        particle.setIntProperty("CTType", CTType).setDoubleProperty("apre", apre).setDoubleProperty("bexp", bexp);

        std::vector<int> exclusions;
        force.getParticleExclusions(ii, exclusions);

        SerializationNode& particleExclusions = particle.createChildNode("ParticleExclusions");
        for (unsigned int jj = 0; jj < exclusions.size(); jj++) {
            particleExclusions.createChildNode("excl").setIntProperty("index", exclusions[jj]);
        }
    }

    SerializationNode& oldTypes = node.createChildNode("oldTypes");
    for (int i = 0; i < force.getNumCTprTypes(); ++i) {
        SerializationNode& anOldType = oldTypes.createChildNode("anOldType");
        anOldType.setIntProperty("oldType", force.getOldCTprType(i));
    }

    SerializationNode& pairs = node.createChildNode("CTPairs");
    for (int i = 0; i < force.getNumCTprTypes(); ++i) {
        for (int j = 0; j < force.getNumCTprTypes(); ++j) {
            SerializationNode& onePair = pairs.createChildNode("onePair");
            double combinedApre, combinedBexp;
            force.getCTprParameters(i, j, combinedApre, combinedBexp);
            onePair.setIntProperty("typeI", i).setIntProperty("typeJ", j).setDoubleProperty("combinedApre", combinedApre).setDoubleProperty("combinedBexp", combinedBexp);
        }
    }
}

void* AmoebaCTForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    AmoebaCTForce* force = new AmoebaCTForce();
    try {
        force->setNumCTprTypes(node.getIntProperty("numCTprTypes"));
        force->setApreCombiningRule(node.getStringProperty("ApreCombiningRule"));
        force->setBexpCombiningRule(node.getStringProperty("BexpCombiningRule"));
        force->setCutoffDistance(node.getDoubleProperty("CTCutoff"));
        force->setNonbondedMethod((AmoebaCTForce::NonbondedMethod) node.getIntProperty("method"));

        const SerializationNode& particles = node.getChildNode("CTParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {
            const SerializationNode& particle = particles.getChildren()[ii];
            force->addParticle(particle.getIntProperty("CTType"),
                particle.getDoubleProperty("apre"), particle.getDoubleProperty("bexp"));

            // exclusions

            const SerializationNode& particleExclusions = particle.getChildNode("ParticleExclusions");
            std::vector<int> exclusions;
            for (unsigned int jj = 0; jj < particleExclusions.getChildren().size(); jj++) {
                exclusions.push_back(particleExclusions.getChildren()[jj].getIntProperty("index"));
            }
            force->setParticleExclusions(ii, exclusions);
        }

        const SerializationNode& oldTypes = node.getChildNode("oldTypes");
        int numCTprTypes = node.getIntProperty("numCTprTypes");
        force->resize(numCTprTypes);
        for (int i = 0; i < numCTprTypes; ++i) {
            force->setOldCTprType(i, oldTypes.getChildren()[i].getIntProperty("oldType"));
        }

        const SerializationNode& pairs = node.getChildNode("CTPairs");
        for (int i = 0; i < numCTprTypes; ++i) {
            for (int j = 0; j < numCTprTypes; ++j) {
                int k = i * numCTprTypes + j;
                const SerializationNode& onePair = pairs.getChildren()[k];
                force->addCTpr(onePair.getIntProperty("typeI"), onePair.getIntProperty("typeJ"),
                    onePair.getDoubleProperty("combinedApre"), onePair.getDoubleProperty("combinedBexp"));
            }
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
