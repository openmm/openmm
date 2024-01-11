/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2023 Stanford University and the Authors.      *
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

#include "openmm/serialization/ATMForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/ATMForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

ATMForceProxy::ATMForceProxy() : SerializationProxy("ATMForce") {
}

void ATMForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 0);
    const ATMForce& force = *reinterpret_cast<const ATMForce*> (object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setStringProperty("energy", force.getEnergyFunction());
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& energyDerivs = node.createChildNode("EnergyParameterDerivatives");
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        energyDerivs.createChildNode("Parameter").setStringProperty("name", force.getEnergyParameterDerivativeName(i));
    }
    SerializationNode& forces = node.createChildNode("Forces");
    for (int i = 0; i < force.getNumForces(); i++) {
        SerializationNode& f = forces.createChildNode("Force");
        f.createChildNode("Force", &force.getForce(i));
    }
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        Vec3 d1, d0;
        force.getParticleParameters(i, d1, d0);
        particles.createChildNode("Particle").setDoubleProperty("d1x", d1[0]).setDoubleProperty("d1y", d1[1]).setDoubleProperty("d1z", d1[2])
                .setDoubleProperty("d0x", d0[0]).setDoubleProperty("d0y", d0[1]).setDoubleProperty("d0z", d0[2]);
    }
}

void* ATMForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 0)
        throw OpenMMException("Unsupported version number");
    ATMForce* force = NULL;
    try {
        ATMForce* force = new ATMForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        const SerializationNode& energyDerivs = node.getChildNode("EnergyParameterDerivatives");
        for (auto& parameter : energyDerivs.getChildren())
            force->addEnergyParameterDerivative(parameter.getStringProperty("name"));
        const SerializationNode& forces = node.getChildNode("Forces");
        for (auto& f : forces.getChildren())
            force->addForce(f.getChildren()[0].decodeObject<Force>());
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& p : particles.getChildren())
            force->addParticle(Vec3(p.getDoubleProperty("d1x"), p.getDoubleProperty("d1y"), p.getDoubleProperty("d1z")),
                               Vec3(p.getDoubleProperty("d0x"), p.getDoubleProperty("d0y"), p.getDoubleProperty("d0z")));
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}
