/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Andrew C. Simmonett, Andreas Krämer                               *
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

#include "openmm/serialization/NoseHooverIntegratorProxy.h"
#include <vector>
#include <OpenMM.h>

using namespace std;
using namespace OpenMM;

NoseHooverIntegratorProxy::NoseHooverIntegratorProxy() : SerializationProxy("NoseHooverIntegrator") {

}

void NoseHooverIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const NoseHooverIntegrator& integrator = *reinterpret_cast<const NoseHooverIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("maximumPairDistance", integrator.getMaximumPairDistance());
    node.setBoolProperty("hasSubsystemThermostats", integrator.hasSubsystemThermostats());
    if (integrator.hasSubsystemThermostats()) {
        // Serialize all thermostats separately
        for (int i = 0; i < integrator.getNumThermostats(); i++){
            const auto& chain = integrator.getThermostat(i);
            auto& chainNode = node.createChildNode("Thermostat");
            chainNode.setDoubleProperty("temperature", chain.getTemperature());
            chainNode.setDoubleProperty("collisionFrequency", chain.getCollisionFrequency());
            chainNode.setDoubleProperty("relativeTemperature", chain.getRelativeTemperature());
            chainNode.setDoubleProperty("relativeCollisionFrequency", chain.getRelativeCollisionFrequency());
            chainNode.setIntProperty("chainLength", chain.getChainLength());
            chainNode.setIntProperty("numMTS", chain.getNumMultiTimeSteps());
            chainNode.setIntProperty("numYS", chain.getNumYoshidaSuzukiTimeSteps());
            auto& particlesNode = chainNode.createChildNode("ThermostatedAtoms");
            for (int particle: chain.getThermostatedAtoms()){
                particlesNode.createChildNode("Particle").setIntProperty("index", particle);
            }
            auto& pairsNode = chainNode.createChildNode("ThermostatedPairs");
            for (auto& pair: chain.getThermostatedPairs()){
                auto& pairNode = pairsNode.createChildNode("Pair");
                pairNode.setIntProperty("index1", pair.first);
                pairNode.setIntProperty("index2", pair.second);
            }
        }
    } else { // Serialize standard thermostat
        node.setDoubleProperty("temperature", integrator.getTemperature());
        node.setDoubleProperty("collisionFrequency", integrator.getCollisionFrequency());
        node.setIntProperty("chainLength", integrator.getThermostat().getChainLength());
        node.setIntProperty("numMTS", integrator.getThermostat().getNumMultiTimeSteps());
        node.setIntProperty("numYS", integrator.getThermostat().getNumYoshidaSuzukiTimeSteps());
    }

}

void* NoseHooverIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    NoseHooverIntegrator* integrator;
    if (node.getBoolProperty("hasSubsystemThermostats")){
        // deserialize all chains
        integrator = new NoseHooverIntegrator(node.getDoubleProperty("stepSize"));
        for (auto& chainNode : node.getChildren()) {
            // particles
            const auto& particlesNode = chainNode.getChildNode("ThermostatedAtoms");
            vector<int> particles;
            for (auto& particleNode: particlesNode.getChildren()){
                particles.push_back(particleNode.getIntProperty("index"));
            }
            // pairs
            const auto& pairsNode = chainNode.getChildNode("ThermostatedPairs");
            vector<pair<int, int>> pairs;
            for (auto& pairNode: pairsNode.getChildren()){
                pairs.emplace_back(pairNode.getIntProperty("index1"), pairNode.getIntProperty("index2"));
            }
            integrator->addSubsystemThermostat(
                particles, pairs, 
                chainNode.getDoubleProperty("temperature"),
                chainNode.getDoubleProperty("collisionFrequency"),
                chainNode.getDoubleProperty("relativeTemperature"),
                chainNode.getDoubleProperty("relativeCollisionFrequency"),
                chainNode.getIntProperty("chainLength"),
                chainNode.getIntProperty("numMTS"),
                chainNode.getIntProperty("numYS")
            );
        }
    } else {
        integrator = new NoseHooverIntegrator(
            node.getDoubleProperty("temperature"),
            node.getDoubleProperty("collisionFrequency"),
            node.getDoubleProperty("stepSize"),
            node.getIntProperty("chainLength"),
            node.getIntProperty("numMTS"),
            node.getIntProperty("numYS")
        );
    }
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setMaximumPairDistance(node.getDoubleProperty("maximumPairDistance"));

    return integrator;
}
