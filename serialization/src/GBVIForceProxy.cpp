/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
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

#include "openmm/serialization/GBVIForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/GBVIForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

GBVIForceProxy::GBVIForceProxy() : SerializationProxy("GBVIForce") {
}

void GBVIForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const GBVIForce& force = *reinterpret_cast<const GBVIForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setIntProperty("scalingMethod", (int) force.getBornRadiusScalingMethod());
    node.setDoubleProperty("quinticLowerLimitFactor", force.getQuinticLowerLimitFactor());
    node.setDoubleProperty("quinticUpperBornRadiusLimit", force.getQuinticUpperBornRadiusLimit());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setDoubleProperty("soluteDielectric", force.getSoluteDielectric());
    node.setDoubleProperty("solventDielectric", force.getSolventDielectric());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, gamma;
        force.getParticleParameters(i, charge, radius, gamma);
        particles.createChildNode("Particle").setDoubleProperty("q", charge).setDoubleProperty("r", radius).setDoubleProperty("gamma", gamma);
    }
    SerializationNode& bonds = node.createChildNode("Bonds");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int particle1, particle2;
        double distance;
        force.getBondParameters(i, particle1, particle2, distance);
        bonds.createChildNode("Bond").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("d", distance);
    }
}

void* GBVIForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1 && node.getIntProperty("version") != 2)
        throw OpenMMException("Unsupported version number");
    GBVIForce* force = new GBVIForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((GBVIForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setSoluteDielectric(node.getDoubleProperty("soluteDielectric"));
        force->setSolventDielectric(node.getDoubleProperty("solventDielectric"));

        if( node.getIntProperty("version") >= 2 ){
            force->setBornRadiusScalingMethod( (GBVIForce::BornRadiusScalingMethod) node.getIntProperty( "scalingMethod"));
            force->setQuinticLowerLimitFactor(node.getDoubleProperty("quinticLowerLimitFactor"));
            force->setQuinticUpperBornRadiusLimit(node.getDoubleProperty("quinticUpperBornRadiusLimit"));
        }

        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getDoubleProperty("q"), particle.getDoubleProperty("r"), particle.getDoubleProperty("gamma"));
        }
        const SerializationNode& bonds = node.getChildNode("Bonds");
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& bond = bonds.getChildren()[i];
            force->addBond(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getDoubleProperty("d"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
