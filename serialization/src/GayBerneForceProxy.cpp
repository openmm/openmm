/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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

#include "openmm/serialization/GayBerneForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/GayBerneForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

GayBerneForceProxy::GayBerneForceProxy() : SerializationProxy("GayBerneForce") {
}

void GayBerneForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const GayBerneForce& force = *reinterpret_cast<const GayBerneForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setBoolProperty("useSwitchingFunction", force.getUseSwitchingFunction());
    node.setDoubleProperty("switchingDistance", force.getSwitchingDistance());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double sigma, epsilon, sx, sy, sz, ex, ey, ez;
        int xparticle, yparticle;
        force.getParticleParameters(i, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
        particles.createChildNode("Particle").setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon).setDoubleProperty("sx", sx)
                .setDoubleProperty("sy", sy).setDoubleProperty("sz", sz).setDoubleProperty("ex", ex).setDoubleProperty("ey", ey).setDoubleProperty("ez", ez)
                .setIntProperty("xparticle", xparticle).setIntProperty("yparticle", yparticle);
    }
    SerializationNode& exceptions = node.createChildNode("Exceptions");
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, sigma, epsilon);
        exceptions.createChildNode("Exception").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon);
    }
}

void* GayBerneForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    GayBerneForce* force = new GayBerneForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((GayBerneForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setUseSwitchingFunction(node.getBoolProperty("useSwitchingFunction", false));
        force->setSwitchingDistance(node.getDoubleProperty("switchingDistance", -1.0));
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren())
            force->addParticle(particle.getDoubleProperty("sig"), particle.getDoubleProperty("eps"), particle.getIntProperty("xparticle"),
                    particle.getIntProperty("yparticle"), particle.getDoubleProperty("sx"), particle.getDoubleProperty("sy"), particle.getDoubleProperty("sz"),
                    particle.getDoubleProperty("ex"), particle.getDoubleProperty("ey"), particle.getDoubleProperty("ez"));
        const SerializationNode& exceptions = node.getChildNode("Exceptions");
        for (auto& exception : exceptions.getChildren())
            force->addException(exception.getIntProperty("p1"), exception.getIntProperty("p2"), exception.getDoubleProperty("sig"), exception.getDoubleProperty("eps"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
