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

#include "openmm/serialization/AmoebaWcaDispersionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaWcaDispersionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaWcaDispersionForceProxy::AmoebaWcaDispersionForceProxy() : SerializationProxy("AmoebaWcaDispersionForce") {
}

void AmoebaWcaDispersionForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const AmoebaWcaDispersionForce& force = *reinterpret_cast<const AmoebaWcaDispersionForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setDoubleProperty("Epso",    force.getEpso());
    node.setDoubleProperty("Epsh",    force.getEpsh());
    node.setDoubleProperty("Rmino",   force.getRmino());
    node.setDoubleProperty("Rminh",   force.getRminh());
    node.setDoubleProperty("Awater",  force.getAwater());
    node.setDoubleProperty("Shctd",   force.getShctd());
    node.setDoubleProperty("Dispoff", force.getDispoff());
    node.setDoubleProperty("Slevy",   force.getSlevy());

    SerializationNode& particles = node.createChildNode("WcaDispersionParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumParticles()); ii++) {
        double radius, epsilon;
        force.getParticleParameters(ii,  radius, epsilon);
        particles.createChildNode("Particle").setDoubleProperty("radius", radius).setDoubleProperty("epsilon", epsilon);
    }

}

void* AmoebaWcaDispersionForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    AmoebaWcaDispersionForce* force = new AmoebaWcaDispersionForce();

    try {
        if (version > 1)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setEpso(   node.getDoubleProperty("Epso"));
        force->setEpsh(   node.getDoubleProperty("Epsh"));
        force->setRmino(  node.getDoubleProperty("Rmino"));
        force->setRminh(  node.getDoubleProperty("Rminh"));


        force->setAwater( node.getDoubleProperty("Awater"));
        force->setShctd(  node.getDoubleProperty("Shctd"));
        force->setDispoff(node.getDoubleProperty("Dispoff"));
        force->setSlevy(  node.getDoubleProperty("Slevy"));

        const SerializationNode& particles = node.getChildNode("WcaDispersionParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {
            const SerializationNode& particle = particles.getChildren()[ii];
            force->addParticle(particle.getDoubleProperty("radius"), particle.getDoubleProperty("epsilon"));
        }

    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
