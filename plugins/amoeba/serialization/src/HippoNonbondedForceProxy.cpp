/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2019 Stanford University and the Authors.      *
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

#include "openmm/serialization/HippoNonbondedForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/HippoNonbondedForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

HippoNonbondedForceProxy::HippoNonbondedForceProxy() : SerializationProxy("HippoNonbondedForce") {
}

void HippoNonbondedForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 0);
    const HippoNonbondedForce& force = *reinterpret_cast<const HippoNonbondedForce*>(object);

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setIntProperty("nonbondedMethod", force.getNonbondedMethod());
    node.setDoubleProperty("cutoffDistance", force.getCutoffDistance());
    node.setDoubleProperty("switchingDistance", force.getSwitchingDistance());
    node.setDoubleProperty("ewaldErrorTolerance",force.getEwaldErrorTolerance());
    double alpha;
    int nx, ny, nz;
    force.getPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("pmeAlpha", alpha);
    node.setIntProperty("pmeGridX", nx);
    node.setIntProperty("pmeGridY", ny);
    node.setIntProperty("pmeGridZ", nz);
    force.getDPMEParameters(alpha, nx, ny, nz);
    node.setDoubleProperty("dpmeAlpha", alpha);
    node.setIntProperty("dpmeGridX", nx);
    node.setIntProperty("dpmeGridY", ny);
    node.setIntProperty("dpmeGridZ", nz);
    SerializationNode& coefficients = node.createChildNode("ExtrapolationCoefficients");
    vector<double> coeff = force.getExtrapolationCoefficients();
    for (int i = 0; i < coeff.size(); i++) {
        stringstream key;
        key << "c" << i;
        coefficients.setDoubleProperty(key.str(), coeff[i]);
    }
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getParticleParameters(i, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                    polarizability, axisType, atomZ, atomX, atomY);
        SerializationNode& particle = particles.createChildNode("Particle");
        particle.setDoubleProperty("charge", charge);
        particle.setDoubleProperty("coreCharge", coreCharge);
        particle.setDoubleProperty("alpha", alpha);
        particle.setDoubleProperty("epsilon", epsilon);
        particle.setDoubleProperty("damping", damping);
        particle.setDoubleProperty("c6", c6);
        particle.setDoubleProperty("pauliK", pauliK);
        particle.setDoubleProperty("pauliQ", pauliQ);
        particle.setDoubleProperty("pauliAlpha", pauliAlpha);
        particle.setDoubleProperty("polarizability", polarizability);
        particle.setIntProperty("axisType", axisType);
        particle.setIntProperty("atomZ", atomZ);
        particle.setIntProperty("atomX", atomX);
        particle.setIntProperty("atomY", atomY);
        SerializationNode& dipoleNode = particle.createChildNode("Dipole");
        dipoleNode.setDoubleProperty("d0", dipole[0]).setDoubleProperty("d1", dipole[1]).setDoubleProperty("d2", dipole[2]);
        SerializationNode& quadrupoleNode = particle.createChildNode("Quadrupole");
        quadrupoleNode.setDoubleProperty("q0", quadrupole[0]).setDoubleProperty("q1", quadrupole[1]).setDoubleProperty("q2", quadrupole[2]);
        quadrupoleNode.setDoubleProperty("q3", quadrupole[3]).setDoubleProperty("q4", quadrupole[4]).setDoubleProperty("q5", quadrupole[5]);
        quadrupoleNode.setDoubleProperty("q6", quadrupole[6]).setDoubleProperty("q7", quadrupole[7]).setDoubleProperty("q8", quadrupole[8]);
    }
    SerializationNode& exceptions = node.createChildNode("Exceptions");
    for (int i = 0; i < force.getNumExceptions(); i++) {
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        int p1, p2;
        force.getExceptionParameters(i, p1, p2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        SerializationNode& exception = exceptions.createChildNode("Exception");
        exception.setIntProperty("p1", p1);
        exception.setIntProperty("p2", p2);
        exception.setDoubleProperty("mmScale", multipoleMultipoleScale);
        exception.setDoubleProperty("dmScale", dipoleMultipoleScale);
        exception.setDoubleProperty("ddScale", dipoleDipoleScale);
        exception.setDoubleProperty("dispScale", dispersionScale);
        exception.setDoubleProperty("repScale", repulsionScale);
        exception.setDoubleProperty("ctScale", chargeTransferScale);
    }
}

void* HippoNonbondedForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 0 || version > 0)
        throw OpenMMException("Unsupported version number");
    HippoNonbondedForce* force = new HippoNonbondedForce();
    try {
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod(static_cast<HippoNonbondedForce::NonbondedMethod>(node.getIntProperty("nonbondedMethod")));
        force->setCutoffDistance(node.getDoubleProperty("cutoffDistance"));
        force->setSwitchingDistance(node.getDoubleProperty("switchingDistance"));
        force->setEwaldErrorTolerance(node.getDoubleProperty("ewaldErrorTolerance"));
        force->setPMEParameters(node.getDoubleProperty("pmeAlpha"), node.getIntProperty("pmeGridX"), node.getIntProperty("pmeGridY"), node.getIntProperty("pmeGridZ"));
        force->setDPMEParameters(node.getDoubleProperty("dpmeAlpha"), node.getIntProperty("dpmeGridX"), node.getIntProperty("dpmeGridY"), node.getIntProperty("dpmeGridZ"));
        const SerializationNode& coefficients = node.getChildNode("ExtrapolationCoefficients");
        vector<double> coeff;
        for (int i = 0; ; i++) {
            stringstream key;
            key << "c" << i;
            if (coefficients.getProperties().find(key.str()) == coefficients.getProperties().end())
                break;
            coeff.push_back(coefficients.getDoubleProperty(key.str()));
        }
        force->setExtrapolationCoefficients(coeff);
        const SerializationNode& particles = node.getChildNode("Particles");
        for (int i = 0; i < particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            std::vector<double> dipole;
            const SerializationNode& dipoleNode = particle.getChildNode("Dipole");
            dipole.push_back(dipoleNode.getDoubleProperty("d0"));
            dipole.push_back(dipoleNode.getDoubleProperty("d1"));
            dipole.push_back(dipoleNode.getDoubleProperty("d2"));
            std::vector<double> quadrupole;
            const SerializationNode& quadrupoleNode = particle.getChildNode("Quadrupole");
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q0"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q1"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q2"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q3"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q4"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q5"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q6"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q7"));
            quadrupole.push_back(quadrupoleNode.getDoubleProperty("q8"));
            force->addParticle(particle.getDoubleProperty("charge"), dipole, quadrupole, particle.getDoubleProperty("coreCharge"),
                    particle.getDoubleProperty("alpha"), particle.getDoubleProperty("epsilon"), particle.getDoubleProperty("damping"),
                    particle.getDoubleProperty("c6"), particle.getDoubleProperty("pauliK"), particle.getDoubleProperty("pauliQ"),
                    particle.getDoubleProperty("pauliAlpha"), particle.getDoubleProperty("polarizability"), particle.getIntProperty("axisType"),
                    particle.getIntProperty("atomZ"), particle.getIntProperty("atomX"), particle.getIntProperty("atomY"));
        }
        const SerializationNode& exceptions = node.getChildNode("Exceptions");
        for (int i = 0; i < exceptions.getChildren().size(); i++) {
            const SerializationNode& exception = exceptions.getChildren()[i];
            force->addException(exception.getIntProperty("p1"), exception.getIntProperty("p2"), exception.getDoubleProperty("mmScale"),
                    exception.getDoubleProperty("dmScale"), exception.getDoubleProperty("ddScale"),
                    exception.getDoubleProperty("dispScale"), exception.getDoubleProperty("repScale"),
                    exception.getDoubleProperty("ctScale"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }

    return force;
}
