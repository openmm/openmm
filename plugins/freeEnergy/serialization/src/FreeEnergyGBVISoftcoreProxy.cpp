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

#include "openmm/serialization/FreeEnergyGBVISoftcoreForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/GBVISoftcoreForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

GBVISoftcoreForceProxy::GBVISoftcoreForceProxy() : SerializationProxy("GBVISoftcoreForce") {
}

void GBVISoftcoreForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const GBVISoftcoreForce& force = *reinterpret_cast<const GBVISoftcoreForce*>(object);
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setIntProperty("brScaleMethod", (int) force.getBornRadiusScalingMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setDoubleProperty("soluteDielectric", force.getSoluteDielectric());
    node.setDoubleProperty("solventDielectric", force.getSolventDielectric());
    node.setDoubleProperty("quinticLwrLmtFctr", force.getQuinticLowerLimitFactor());
    node.setDoubleProperty("quinticUpprBrLmt", force.getQuinticUpperBornRadiusLimit());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, gamma, bornRadiusScaleFactor;
        force.getParticleParameters(i, charge, radius, gamma, bornRadiusScaleFactor);
        particles.createChildNode("Particle").setDoubleProperty("q", charge).setDoubleProperty("r", radius).setDoubleProperty("gamma", gamma).setDoubleProperty("bRSclFctr", bornRadiusScaleFactor);
    }
    SerializationNode& bonds = node.createChildNode("Bonds");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int particle1, particle2;
        double distance;
        force.getBondParameters(i, particle1, particle2, distance);
        bonds.createChildNode("Bond").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("d", distance);
    }
}

void* GBVISoftcoreForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    GBVISoftcoreForce* force = new GBVISoftcoreForce();
    try {
        force->setNonbondedMethod((GBVISoftcoreForce::NonbondedSoftcoreMethod) node.getIntProperty("method"));
        force->setBornRadiusScalingMethod((GBVISoftcoreForce::BornRadiusScalingSoftcoreMethod) node.getIntProperty("brScaleMethod"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setSoluteDielectric(node.getDoubleProperty("soluteDielectric"));
        force->setSolventDielectric(node.getDoubleProperty("solventDielectric"));
        force->setQuinticLowerLimitFactor(node.getDoubleProperty("quinticLwrLmtFctr"));
        force->setQuinticUpperBornRadiusLimit(node.getDoubleProperty("quinticUpprBrLmt"));
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getDoubleProperty("q"), particle.getDoubleProperty("r"), particle.getDoubleProperty("gamma"), particle.getDoubleProperty("bRSclFctr"));
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
