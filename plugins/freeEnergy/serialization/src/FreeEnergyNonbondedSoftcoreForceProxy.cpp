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

#include "openmm/serialization/FreeEnergyNonbondedSoftcoreForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/NonbondedSoftcoreForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

NonbondedSoftcoreForceProxy::NonbondedSoftcoreForceProxy() : SerializationProxy("NonbondedSoftcoreForce") {
}

void NonbondedSoftcoreForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const NonbondedSoftcoreForce& force = *reinterpret_cast<const NonbondedSoftcoreForce*>(object);
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setDoubleProperty("ewaldTolerance", force.getEwaldErrorTolerance());
    node.setDoubleProperty("rfDielectric", force.getReactionFieldDielectric());
    //node.setIntProperty("dispersionCorrection", force.getUseDispersionCorrection());
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, sigma, epsilon, softcoreLJLambda;
        force.getParticleParameters(i, charge, sigma, epsilon, softcoreLJLambda);
        particles.createChildNode("Particle").setDoubleProperty("q", charge).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon).setDoubleProperty("lambda", softcoreLJLambda);
    }
    SerializationNode& exceptions = node.createChildNode("Exceptions");
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon, softcoreLJLambda;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exceptions.createChildNode("Exception").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("q", chargeProd).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon).setDoubleProperty("lambda", softcoreLJLambda);
    }
}

void* NonbondedSoftcoreForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    NonbondedSoftcoreForce* force = new NonbondedSoftcoreForce();
    try {
        force->setNonbondedMethod((NonbondedSoftcoreForce::NonbondedSoftcoreMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setEwaldErrorTolerance(node.getDoubleProperty("ewaldTolerance"));
        force->setReactionFieldDielectric(node.getDoubleProperty("rfDielectric"));
        //force->setUseDispersionCorrection(node.getIntProperty("dispersionCorrection"));
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getDoubleProperty("q"), particle.getDoubleProperty("sig"), particle.getDoubleProperty("eps"), particle.getDoubleProperty("lambda"));
        }
        const SerializationNode& exceptions = node.getChildNode("Exceptions");
        for (int i = 0; i < (int) exceptions.getChildren().size(); i++) {
            const SerializationNode& exception = exceptions.getChildren()[i];
            force->addException(exception.getIntProperty("p1"), exception.getIntProperty("p2"), exception.getDoubleProperty("q"), exception.getDoubleProperty("sig"), exception.getDoubleProperty("eps"), exception.getDoubleProperty("lambda"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
