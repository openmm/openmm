/* -------------------------------------------------------------------------- *
 *                                OpenMMDrude                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "openmm/serialization/DrudeForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/DrudeForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

DrudeForceProxy::DrudeForceProxy() : SerializationProxy("DrudeForce") {
}

void DrudeForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const DrudeForce& force = *reinterpret_cast<const DrudeForce*>(object);
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p, p1, p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        force.getParticleParameters(i, p, p1, p2, p3, p4, charge, polarizability, aniso12, aniso34);
        particles.createChildNode("Particle").setIntProperty("p", p).setIntProperty("p1", p1).setIntProperty("p2", p2).setIntProperty("p3", p3).setIntProperty("p4", p4)
            .setDoubleProperty("charge", charge).setDoubleProperty("polarizability", polarizability).setDoubleProperty("a12", aniso12).setDoubleProperty("a34", aniso34);
    }
    SerializationNode& pairs = node.createChildNode("ScreenedPairs");
    for (int i = 0; i < force.getNumScreenedPairs(); i++) {
        int p1, p2;
        double thole;
        force.getScreenedPairParameters(i, p1, p2, thole);
        pairs.createChildNode("Pair").setIntProperty("p1", p1).setIntProperty("p2", p2).setDoubleProperty("thole", thole);
    }
}

void* DrudeForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    DrudeForce* force = new DrudeForce();
    try {
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren())
            force->addParticle(particle.getIntProperty("p"), particle.getIntProperty("p1"), particle.getIntProperty("p2"), particle.getIntProperty("p3"), particle.getIntProperty("p4"),
                    particle.getDoubleProperty("charge"), particle.getDoubleProperty("polarizability"), particle.getDoubleProperty("a12"), particle.getDoubleProperty("a34"));
        const SerializationNode& pairs = node.getChildNode("ScreenedPairs");
        for (auto& pair : pairs.getChildren())
            force->addScreenedPair(pair.getIntProperty("p1"), pair.getIntProperty("p2"), pair.getDoubleProperty("thole"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
