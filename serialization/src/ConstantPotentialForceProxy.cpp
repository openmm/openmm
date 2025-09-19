/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
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

#include "openmm/serialization/ConstantPotentialForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/ConstantPotentialForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

ConstantPotentialForceProxy::ConstantPotentialForceProxy() : SerializationProxy("ConstantPotentialForce") {
}

void ConstantPotentialForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ConstantPotentialForce& force = *reinterpret_cast<const ConstantPotentialForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setDoubleProperty("ewaldTolerance", force.getEwaldErrorTolerance());
    double alpha;
    int nx, ny, nz;
    force.getPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("alpha", alpha);
    node.setIntProperty("nx", nx);
    node.setIntProperty("ny", ny);
    node.setIntProperty("nz", nz);
    node.setBoolProperty("exceptionsUsePeriodic", force.getExceptionsUsePeriodicBoundaryConditions());
    node.setIntProperty("constantPotentialMethod", (int) force.getConstantPotentialMethod());
    node.setDoubleProperty("cgTolerance", force.getCGErrorTolerance());
    node.setBoolProperty("usePreconditioner", force.getUsePreconditioner());
    node.setBoolProperty("useChargeConstraint", force.getUseChargeConstraint());
    node.setDoubleProperty("chargeConstraintTarget", force.getChargeConstraintTarget());
    Vec3 externalField;
    force.getExternalField(externalField);
    node.createChildNode("ExternalField").setDoubleProperty("x", externalField[0]).setDoubleProperty("y", externalField[1]).setDoubleProperty("z", externalField[2]);
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge;
        force.getParticleParameters(i, charge);
        particles.createChildNode("Particle").setDoubleProperty("q", charge);
    }
    SerializationNode& exceptions = node.createChildNode("Exceptions");
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd;
        force.getExceptionParameters(i, particle1, particle2, chargeProd);
        exceptions.createChildNode("Exception").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("q", chargeProd);
    }
    SerializationNode& electrodes = node.createChildNode("Electrodes");
    for(int i = 0; i < force.getNumElectrodes(); i++) {
        SerializationNode& electrode = electrodes.createChildNode("Electrode");
        std::set<int> electrodeParticles;
        double potential;
        double gaussianWidth;
        double thomasFermiScale;
        force.getElectrodeParameters(i, electrodeParticles, potential, gaussianWidth, thomasFermiScale);
        electrode.setDoubleProperty("potential", potential);
        electrode.setDoubleProperty("gaussianWidth", gaussianWidth);
        electrode.setDoubleProperty("thomasFermiScale", thomasFermiScale);
        SerializationNode& electrodeParticlesNode = electrode.createChildNode("Particles");
        for (int p : electrodeParticles) {
            electrodeParticlesNode.createChildNode("Particle").setIntProperty("index", p);
        }
    }
}

void* ConstantPotentialForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 1)
        throw OpenMMException("Unsupported version number");
    ConstantPotentialForce* force = new ConstantPotentialForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setEwaldErrorTolerance(node.getDoubleProperty("ewaldTolerance"));
        double alpha = node.getDoubleProperty("alpha", 0.0);
        int nx = node.getIntProperty("nx", 0);
        int ny = node.getIntProperty("ny", 0);
        int nz = node.getIntProperty("nz", 0);
        force->setPMEParameters(alpha, nx, ny, nz);
        force->setExceptionsUsePeriodicBoundaryConditions(node.getBoolProperty("exceptionsUsePeriodic"));
        force->setConstantPotentialMethod((ConstantPotentialForce::ConstantPotentialMethod) node.getIntProperty("constantPotentialMethod"));
        force->setCGErrorTolerance(node.getDoubleProperty("cgTolerance"));
        force->setUsePreconditioner(node.getBoolProperty("usePreconditioner"));
        force->setUseChargeConstraint(node.getBoolProperty("useChargeConstraint"));
        force->setChargeConstraintTarget(node.getDoubleProperty("chargeConstraintTarget"));
        const SerializationNode& externalFieldNode = node.getChildNode("ExternalField");
        Vec3 externalField(externalFieldNode.getDoubleProperty("x"), externalFieldNode.getDoubleProperty("y"), externalFieldNode.getDoubleProperty("z"));
        force->setExternalField(externalField);
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren())
            force->addParticle(particle.getDoubleProperty("q"));
        const SerializationNode& exceptions = node.getChildNode("Exceptions");
        for (auto& exception : exceptions.getChildren())
            force->addException(exception.getIntProperty("p1"), exception.getIntProperty("p2"), exception.getDoubleProperty("q"));
        const SerializationNode& electrodes = node.getChildNode("Electrodes");
        for (auto& electrode : electrodes.getChildren()) {
            std::set<int> electrodeParticles;
            const SerializationNode& electrodeParticlesNode = electrode.getChildNode("Particles");
            for (auto& child : electrodeParticlesNode.getChildren()) {
                electrodeParticles.insert(child.getIntProperty("index"));
            }
            force->addElectrode(electrodeParticles, electrode.getDoubleProperty("potential"), electrode.getDoubleProperty("gaussianWidth"), electrode.getDoubleProperty("thomasFermiScale"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
