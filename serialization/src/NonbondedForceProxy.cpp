/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2018 Stanford University and the Authors.      *
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

#include "openmm/serialization/NonbondedForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/NonbondedForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

NonbondedForceProxy::NonbondedForceProxy() : SerializationProxy("NonbondedForce") {
}

void NonbondedForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 3);
    const NonbondedForce& force = *reinterpret_cast<const NonbondedForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setBoolProperty("useSwitchingFunction", force.getUseSwitchingFunction());
    node.setDoubleProperty("switchingDistance", force.getSwitchingDistance());
    node.setDoubleProperty("ewaldTolerance", force.getEwaldErrorTolerance());
    node.setDoubleProperty("rfDielectric", force.getReactionFieldDielectric());
    node.setIntProperty("dispersionCorrection", force.getUseDispersionCorrection());
    double alpha;
    int nx, ny, nz;
    force.getPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("alpha", alpha);
    node.setIntProperty("nx", nx);
    node.setIntProperty("ny", ny);
    node.setIntProperty("nz", nz);
    force.getLJPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("ljAlpha", alpha);
    node.setIntProperty("ljnx", nx);
    node.setIntProperty("ljny", ny);
    node.setIntProperty("ljnz", nz);
    node.setIntProperty("recipForceGroup", force.getReciprocalSpaceForceGroup());
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    SerializationNode& particleOffsets = node.createChildNode("ParticleOffsets");
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        int particle;
        double chargeScale, sigmaScale, epsilonScale;
        string parameter;
        force.getParticleParameterOffset(i, parameter, particle, chargeScale, sigmaScale, epsilonScale);
        particleOffsets.createChildNode("Offset").setStringProperty("parameter", parameter).setIntProperty("particle", particle).setDoubleProperty("q", chargeScale).setDoubleProperty("sig", sigmaScale).setDoubleProperty("eps", epsilonScale);
    }
    SerializationNode& exceptionOffsets = node.createChildNode("ExceptionOffsets");
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        int exception;
        double chargeProdScale, sigmaScale, epsilonScale;
        string parameter;
        force.getExceptionParameterOffset(i, parameter, exception, chargeProdScale, sigmaScale, epsilonScale);
        exceptionOffsets.createChildNode("Offset").setStringProperty("parameter", parameter).setIntProperty("exception", exception).setDoubleProperty("q", chargeProdScale).setDoubleProperty("sig", sigmaScale).setDoubleProperty("eps", epsilonScale);
    }
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        particles.createChildNode("Particle").setDoubleProperty("q", charge).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon);
    }
    SerializationNode& exceptions = node.createChildNode("Exceptions");
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exceptions.createChildNode("Exception").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setDoubleProperty("q", chargeProd).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon);
    }
}

void* NonbondedForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 3)
        throw OpenMMException("Unsupported version number");
    NonbondedForce* force = new NonbondedForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((NonbondedForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setUseSwitchingFunction(node.getBoolProperty("useSwitchingFunction", false));
        force->setSwitchingDistance(node.getDoubleProperty("switchingDistance", -1.0));
        force->setEwaldErrorTolerance(node.getDoubleProperty("ewaldTolerance"));
        force->setReactionFieldDielectric(node.getDoubleProperty("rfDielectric"));
        force->setUseDispersionCorrection(node.getIntProperty("dispersionCorrection"));
        double alpha = node.getDoubleProperty("alpha", 0.0);
        int nx = node.getIntProperty("nx", 0);
        int ny = node.getIntProperty("ny", 0);
        int nz = node.getIntProperty("nz", 0);
        force->setPMEParameters(alpha, nx, ny, nz);
        if (version >= 2) {
            alpha = node.getDoubleProperty("ljAlpha", 0.0);
            nx = node.getIntProperty("ljnx", 0);
            ny = node.getIntProperty("ljny", 0);
            nz = node.getIntProperty("ljnz", 0);
            force->setLJPMEParameters(alpha, nx, ny, nz);
        }
        force->setReciprocalSpaceForceGroup(node.getIntProperty("recipForceGroup", -1));
        if (version >= 3) {
            const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
            for (auto& parameter : globalParams.getChildren())
                force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
            const SerializationNode& particleOffsets = node.getChildNode("ParticleOffsets");
            for (auto& offset : particleOffsets.getChildren())
                force->addParticleParameterOffset(offset.getStringProperty("parameter"), offset.getIntProperty("particle"), offset.getDoubleProperty("q"), offset.getDoubleProperty("sig"), offset.getDoubleProperty("eps"));
            const SerializationNode& exceptionOffsets = node.getChildNode("ExceptionOffsets");
            for (auto& offset : exceptionOffsets.getChildren())
                force->addExceptionParameterOffset(offset.getStringProperty("parameter"), offset.getIntProperty("exception"), offset.getDoubleProperty("q"), offset.getDoubleProperty("sig"), offset.getDoubleProperty("eps"));
        }
        const SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren())
            force->addParticle(particle.getDoubleProperty("q"), particle.getDoubleProperty("sig"), particle.getDoubleProperty("eps"));
        const SerializationNode& exceptions = node.getChildNode("Exceptions");
        for (auto& exception : exceptions.getChildren())
            force->addException(exception.getIntProperty("p1"), exception.getIntProperty("p2"), exception.getDoubleProperty("q"), exception.getDoubleProperty("sig"), exception.getDoubleProperty("eps"));
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
