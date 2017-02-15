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

#include "openmm/serialization/AmoebaVdwForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaVdwForce.h"
#include <sstream>
#include <cmath>

using namespace OpenMM;
using namespace std;

AmoebaVdwForceProxy::AmoebaVdwForceProxy()
    : SerializationProxy("AmoebaVdwForce") {}

void AmoebaVdwForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const AmoebaVdwForce& force = *reinterpret_cast<const AmoebaVdwForce*>(object);

    node.setIntProperty("method", (int)force.getNonbondedMethod());
    node.setDoubleProperty("VdwCutoff", force.getCutoff());
    node.setBoolProperty("useDispersionCorrection", force.getUseDispersionCorrection());
    node.setIntProperty("numVdwprTypes", force.getNumVdwprTypes());
    node.setStringProperty("SigmaCombiningRule", force.getSigmaCombiningRule());
    node.setStringProperty("EpsilonCombiningRule", force.getEpsilonCombiningRule());
    node.setStringProperty("FunctionalForm", force.getFunctionalForm());

    SerializationNode& particles = node.createChildNode("VdwParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumParticles()); ii++) {

        int ivIndex, vdwType;
        double sigma, epsilon, reductionFactor, lambda;
        force.getParticleParameters(ii, ivIndex, vdwType, sigma, epsilon, reductionFactor, lambda);

        SerializationNode& particle = particles.createChildNode("Particle");
        particle.setIntProperty("ivIndex", ivIndex).setIntProperty("vdwType", vdwType).setDoubleProperty("sigma", sigma).setDoubleProperty("epsilon", epsilon).setDoubleProperty("reductionFactor", reductionFactor).setDoubleProperty("lambda", lambda);

        std::vector<int> exclusions;
        force.getParticleExclusions(ii, exclusions);

        SerializationNode& particleExclusions = particle.createChildNode("ParticleExclusions");
        for (unsigned int jj = 0; jj < exclusions.size(); jj++) {
            particleExclusions.createChildNode("excl").setIntProperty("index", exclusions[jj]);
        }
    }

    SerializationNode& oldTypes = node.createChildNode("oldTypes");
    for (int i = 0; i < force.getNumVdwprTypes(); ++i) {
        SerializationNode& anOldType = oldTypes.createChildNode("anOldType");
        anOldType.setIntProperty("oldType", force.getOldVdwprType(i));
    }

    SerializationNode& pairs = node.createChildNode("vdwPairs");
    for (int i = 0; i < force.getNumVdwprTypes(); ++i) {
        for (int j = 0; j < force.getNumVdwprTypes(); ++j) {
            SerializationNode& onePair = pairs.createChildNode("onePair");
            double combinedSigma, combinedEpsilon;
            force.getVdwprParameters(i, j, combinedSigma, combinedEpsilon);
            onePair.setIntProperty("typeI", i).setIntProperty("typeJ", j).setDoubleProperty("combinedSigma", combinedSigma).setDoubleProperty("combinedEpsilon", combinedEpsilon);
        }
    }
}

void* AmoebaVdwForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    AmoebaVdwForce* force = new AmoebaVdwForce();
    try {
        force->setNonbondedMethod((AmoebaVdwForce::NonbondedMethod)node.getIntProperty("method"));
        force->setCutoff(node.getDoubleProperty("VdwCutoff"));
        force->setUseDispersionCorrection(node.getBoolProperty("useDispersionCorrection"));
        force->setNumVdwprTypes(node.getIntProperty("numVdwprTypes"));
        force->setSigmaCombiningRule(node.getStringProperty("SigmaCombiningRule"));
        force->setEpsilonCombiningRule(node.getStringProperty("EpsilonCombiningRule"));
        force->setFunctionalForm(node.getStringProperty("FunctionalForm"));
        const SerializationNode& particles = node.getChildNode("VdwParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {
            const SerializationNode& particle = particles.getChildren()[ii];
            force->addParticle(particle.getIntProperty("ivIndex"), particle.getIntProperty("vdwType"),
                particle.getDoubleProperty("sigma"), particle.getDoubleProperty("epsilon"),
                particle.getDoubleProperty("reductionFactor"), particle.getDoubleProperty("lambda"));

            // exclusions

            const SerializationNode& particleExclusions = particle.getChildNode("ParticleExclusions");
            std::vector<int> exclusions;
            for (unsigned int jj = 0; jj < particleExclusions.getChildren().size(); jj++) {
                exclusions.push_back(particleExclusions.getChildren()[jj].getIntProperty("index"));
            }
            force->setParticleExclusions(ii, exclusions);
        }

        const SerializationNode& oldTypes = node.getChildNode("oldTypes");
        int numVdwprTypes = node.getIntProperty("numVdwprTypes");
        force->resize(numVdwprTypes);
        for (int i = 0; i < numVdwprTypes; ++i) {
            force->setOldVdwprType(i, oldTypes.getChildren()[i].getIntProperty("oldType"));
        }

        const SerializationNode& pairs = node.getChildNode("vdwPairs");
        for (int i = 0; i < numVdwprTypes; ++i) {
            for (int j = 0; j < numVdwprTypes; ++j) {
                int k = i * numVdwprTypes + j;
                const SerializationNode& onePair = pairs.getChildren()[k];
                force->addVdwpr(onePair.getIntProperty("typeI"), onePair.getIntProperty("typeJ"),
                    onePair.getDoubleProperty("combinedSigma"), onePair.getDoubleProperty("combinedEpsilon"));
            }
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
