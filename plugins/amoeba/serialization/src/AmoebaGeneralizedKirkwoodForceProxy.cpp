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

#include "openmm/serialization/AmoebaGeneralizedKirkwoodForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaGeneralizedKirkwoodForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaGeneralizedKirkwoodForceProxy::AmoebaGeneralizedKirkwoodForceProxy() : SerializationProxy("AmoebaGeneralizedKirkwoodForce") {
}

void AmoebaGeneralizedKirkwoodForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const AmoebaGeneralizedKirkwoodForce& force = *reinterpret_cast<const AmoebaGeneralizedKirkwoodForce*>(object);

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setDoubleProperty("GeneralizedKirkwoodSolventDielectric", force.getSolventDielectric());
    node.setDoubleProperty("GeneralizedKirkwoodSoluteDielectric",  force.getSoluteDielectric());
    //node.setDoubleProperty("GeneralizedKirkwoodDielectricOffset",  force.getDielectricOffset());
    node.setDoubleProperty("GeneralizedKirkwoodProbeRadius",       force.getProbeRadius());
    node.setDoubleProperty("GeneralizedKirkwoodSurfaceAreaFactor", force.getSurfaceAreaFactor());
    node.setIntProperty(  "GeneralizedKirkwoodIncludeCavityTerm", force.getIncludeCavityTerm());

    SerializationNode& particles = node.createChildNode("GeneralizedKirkwoodParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumParticles()); ii++) {
        double radius, charge, scalingFactor;
        force.getParticleParameters(ii, charge, radius, scalingFactor);
        particles.createChildNode("Particle").setDoubleProperty("charge", charge).setDoubleProperty("radius", radius).setDoubleProperty("scaleFactor", scalingFactor);
    }

}

void* AmoebaGeneralizedKirkwoodForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    AmoebaGeneralizedKirkwoodForce* force = new AmoebaGeneralizedKirkwoodForce();
    try {
        if (version > 1)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setSolventDielectric(  node.getDoubleProperty("GeneralizedKirkwoodSolventDielectric"));
        force->setSoluteDielectric(   node.getDoubleProperty("GeneralizedKirkwoodSoluteDielectric"));
        //force->setDielectricOffset(   node.getDoubleProperty("GeneralizedKirkwoodDielectricOffset"));
        force->setProbeRadius(        node.getDoubleProperty("GeneralizedKirkwoodProbeRadius"));
        force->setSurfaceAreaFactor(  node.getDoubleProperty("GeneralizedKirkwoodSurfaceAreaFactor"));
        force->setIncludeCavityTerm(  node.getIntProperty(   "GeneralizedKirkwoodIncludeCavityTerm"));

        const SerializationNode& particles = node.getChildNode("GeneralizedKirkwoodParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {
            const SerializationNode& particle = particles.getChildren()[ii];
            force->addParticle(particle.getDoubleProperty("charge"), particle.getDoubleProperty("radius"), particle.getDoubleProperty("scaleFactor"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}
