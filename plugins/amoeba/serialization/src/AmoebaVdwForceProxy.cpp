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

#include "openmm/serialization/AmoebaVdwForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaVdwForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaVdwForceProxy::AmoebaVdwForceProxy() : SerializationProxy("AmoebaVdwForce") {
}

void AmoebaVdwForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 3);
    const AmoebaVdwForce& force = *reinterpret_cast<const AmoebaVdwForce*>(object);

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("SigmaCombiningRule", force.getSigmaCombiningRule());
    node.setStringProperty("EpsilonCombiningRule", force.getEpsilonCombiningRule());
    node.setDoubleProperty("VdwCutoff", force.getCutoffDistance());
    node.setIntProperty("method", (int) force.getNonbondedMethod());

    node.setDoubleProperty("n", force.getSoftcorePower());
    node.setDoubleProperty("alpha", force.getSoftcoreAlpha());
    node.setIntProperty("alchemicalMethod", (int) force.getAlchemicalMethod());

    SerializationNode& particles = node.createChildNode("VdwParticles");
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumParticles()); ii++) {

        int ivIndex;
        double sigma, epsilon, reductionFactor;
        bool isAlchemical;
        force.getParticleParameters(ii, ivIndex, sigma, epsilon, reductionFactor, isAlchemical);

        SerializationNode& particle = particles.createChildNode("Particle");
        particle.setIntProperty("ivIndex", ivIndex).setDoubleProperty("sigma", sigma).setDoubleProperty("epsilon", epsilon).setDoubleProperty("reductionFactor", reductionFactor).setBoolProperty("isAlchemical",isAlchemical);

        std::vector< int > exclusions;
        force.getParticleExclusions(ii,  exclusions);

        SerializationNode& particleExclusions = particle.createChildNode("ParticleExclusions");
        for (unsigned int jj = 0; jj < exclusions.size(); jj++) {
            particleExclusions.createChildNode("excl").setIntProperty("index", exclusions[jj]);
        }
    }

    node.setBoolProperty("useDispersionCorrection", force.getUseDispersionCorrection());
    node.setBoolProperty("usesVdwpr", force.getUsePairwiseVdw());
    node.setBoolProperty("usesLJ", force.getUseLennardJones());
    if (force.getUsePairwiseVdw()) {
        SerializationNode& vdwTypes = node.createChildNode("VdwTypes");
        for (int ii = 0; ii < force.getNumParticles(); ++ii) {
            int ityp;
            force.getCondensedType(ii, ityp);
            SerializationNode& vdwType = vdwTypes.createChildNode("VdwType");
            vdwType.setIntProperty("vdwtype", ityp);
        }

        int ntypes = force.getNumCondensedTypes();
        node.setIntProperty("NumVdwTypes", ntypes);
        SerializationNode& sigEpsTable = node.createChildNode("SigEpsTable");
        for (int ii = 0; ii < ntypes; ++ii) {
            for (int jj = 0; jj < ntypes; ++jj) {
                double sig, eps;
                force.getPairSigmaEpsilon(ii, jj, sig, eps);
                SerializationNode& sigEps = sigEpsTable.createChildNode("SigEps");
                sigEps.setDoubleProperty("sig", sig).setDoubleProperty("eps", eps);
            }
        }
    }
}

void* AmoebaVdwForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 3)
        throw OpenMMException("Unsupported version number");
    AmoebaVdwForce* force = new AmoebaVdwForce();
    try {
        if (version > 1)
            force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setSigmaCombiningRule(node.getStringProperty("SigmaCombiningRule"));
        force->setEpsilonCombiningRule(node.getStringProperty("EpsilonCombiningRule"));
        force->setCutoffDistance(node.getDoubleProperty("VdwCutoff"));
        force->setNonbondedMethod((AmoebaVdwForce::NonbondedMethod) node.getIntProperty("method"));

        const SerializationNode& particles = node.getChildNode("VdwParticles");
        for (unsigned int ii = 0; ii < particles.getChildren().size(); ii++) {
            const SerializationNode& particle = particles.getChildren()[ii];

            if (version < 3)
               force->addParticle(particle.getIntProperty("ivIndex"), particle.getDoubleProperty("sigma"),
                                  particle.getDoubleProperty("epsilon"), particle.getDoubleProperty("reductionFactor"));
            else
               force->addParticle(particle.getIntProperty("ivIndex"), particle.getDoubleProperty("sigma"),
                                  particle.getDoubleProperty("epsilon"), particle.getDoubleProperty("reductionFactor"),
                                  particle.getBoolProperty("isAlchemical"));

            // exclusions

            const SerializationNode& particleExclusions = particle.getChildNode("ParticleExclusions");
            std::vector< int > exclusions;
            for (unsigned int jj = 0; jj < particleExclusions.getChildren().size(); jj++) {
                exclusions.push_back(particleExclusions.getChildren()[jj].getIntProperty("index"));
            }
            force->setParticleExclusions(ii, exclusions);
        }

        if (version > 2) {
           force->setAlchemicalMethod((AmoebaVdwForce::AlchemicalMethod) node.getIntProperty("alchemicalMethod"));
           force->setSoftcorePower(node.getDoubleProperty("n"));
           force->setSoftcoreAlpha(node.getDoubleProperty("alpha"));

           force->setUseDispersionCorrection(node.getBoolProperty("useDispersionCorrection"));
           bool usesVdwpr = node.getBoolProperty("usesVdwpr");
           if (usesVdwpr) {
               bool usesLJ = node.getBoolProperty("usesLJ");
               force->setUseLennardJones();
               const SerializationNode& vdwTypes = node.getChildNode("VdwTypes");
               for (int ii = 0; ii < vdwTypes.getChildren().size(); ++ii) {
                   force->addCondensedType(vdwTypes.getChildren()[ii].getIntProperty("vdwtype"));
               }

               int ntypes = node.getIntProperty("NumVdwTypes");
               const SerializationNode& sigEpsTable = node.getChildNode("SigEpsTable");
               for (int ii = 0; ii < sigEpsTable.getChildren().size(); ++ii) {
                   int ityp = ii / ntypes;
                   int jtyp = ii % ntypes;
                   const SerializationNode& sigeps = sigEpsTable.getChildren()[ii];
                   force->setPairSigmaEpsilon(ityp, jtyp, sigeps.getDoubleProperty("sig"), sigeps.getDoubleProperty("eps"));
               }
           }
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
