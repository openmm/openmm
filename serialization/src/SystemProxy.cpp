/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/System.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

SystemProxy::SystemProxy() : SerializationProxy("System") {
}

void SystemProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const System& system = *reinterpret_cast<const System*>(object);
    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors(a, b, c);
    SerializationNode& box = node.createChildNode("PeriodicBoxVectors");
    box.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    box.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    box.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < system.getNumParticles(); i++)
        particles.createChildNode("Particle").setDoubleProperty("mass", system.getParticleMass(i));
    SerializationNode& constraints = node.createChildNode("Constraints");
    for (int i = 0; i < system.getNumConstraints(); i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraints.createChildNode("Constraint").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("d", distance);
    }
    SerializationNode& forces = node.createChildNode("Forces");
    for (int i = 0; i < system.getNumForces(); i++)
        forces.createChildNode("Force", &system.getForce(i));
}

void* SystemProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    System* system = new System();
    try {
        const SerializationNode& box = node.getChildNode("PeriodicBoxVectors");
        const SerializationNode& boxa = box.getChildNode("A");
        const SerializationNode& boxb = box.getChildNode("B");
        const SerializationNode& boxc = box.getChildNode("C");
        Vec3 a(boxa.getDoubleProperty("x"), boxa.getDoubleProperty("y"), boxa.getDoubleProperty("z"));
        Vec3 b(boxb.getDoubleProperty("x"), boxb.getDoubleProperty("y"), boxb.getDoubleProperty("z"));
        Vec3 c(boxc.getDoubleProperty("x"), boxc.getDoubleProperty("y"), boxc.getDoubleProperty("z"));
        system->setDefaultPeriodicBoxVectors(a, b, c);
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++)
            system->addParticle(particles.getChildren()[i].getDoubleProperty("mass"));
        const SerializationNode& constraints = node.getChildNode("Constraints");
        for (int i = 0; i < (int) constraints.getChildren().size(); i++) {
            const SerializationNode& constraint = constraints.getChildren()[i];
            system->addConstraint(constraint.getDoubleProperty("p1"), constraint.getDoubleProperty("p2"), constraint.getDoubleProperty("d"));
        }
        const SerializationNode& forces = node.getChildNode("Forces");
        for (int i = 0; i < (int) forces.getChildren().size(); i++) {
            system->addForce(forces.getChildren()[i].decodeObject<Force>());
        }
    }
    catch (...) {
        delete system;
        throw;
    }
    return system;
}