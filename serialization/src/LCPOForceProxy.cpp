/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Evan Pretti                                        *
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

#include "openmm/serialization/LCPOForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/LCPOForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

LCPOForceProxy::LCPOForceProxy() : SerializationProxy("LCPOForce") {
}

void LCPOForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const LCPOForce& force = *reinterpret_cast<const LCPOForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
    node.setDoubleProperty("surfaceTension", force.getSurfaceTension());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double radius, p1, p2, p3, p4;
        force.getParticleParameters(i, radius, p1, p2, p3, p4);
        particles.createChildNode("Particle").setDoubleProperty("r", radius).setDoubleProperty("p1", p1).setDoubleProperty("p2", p2).setDoubleProperty("p3", p3).setDoubleProperty("p4", p4);
    }
}

void* LCPOForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 1) {
        throw OpenMMException("Unsupported version number");
    }
    LCPOForce* force = new LCPOForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
        force->setSurfaceTension(node.getDoubleProperty("surfaceTension"));
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren()) {
            force->addParticle(particle.getDoubleProperty("r"), particle.getDoubleProperty("p1"), particle.getDoubleProperty("p2"), particle.getDoubleProperty("p3"), particle.getDoubleProperty("p4"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
