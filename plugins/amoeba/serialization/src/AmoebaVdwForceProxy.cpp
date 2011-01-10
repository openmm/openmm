/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
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
    node.setIntProperty("version", 1);
    const AmoebaVdwForce& force = *reinterpret_cast<const AmoebaVdwForce*>(object);
    node.setStringProperty("SigmaCombiningRule", force.getSigmaCombiningRule());
    node.setStringProperty("EpsilonCombiningRule", force.getEpsilonCombiningRule());
    node.setDoubleProperty("VdwCutoff", force.getCutoff());
    node.setIntProperty("VdwUseNeighborList", force.getUseNeighborList());
    node.setIntProperty("VdwPBC", force.getPBC());
    SerializationNode& particles = node.createChildNode("VdwParticles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int ivIndex, classIndex;
        double sigma, epsilon, reductionFactor;
        force.getParticleParameters( i,  ivIndex, classIndex, sigma, epsilon, reductionFactor );
        particles.createChildNode("Particle").setIntProperty("index", i).setIntProperty("ivIndex", ivIndex).setIntProperty("classIndex", classIndex).setDoubleProperty("sig", sigma).setDoubleProperty("eps", epsilon).setDoubleProperty("red", reductionFactor);
    }

    SerializationNode& particleExclusions = node.createChildNode("VdwParticleExclusions");
    for (int i = 0; i < force.getNumParticles(); i++) {
        std::vector< int > exclusions;
        force.getParticleExclusions( i,  exclusions );
        SerializationNode& particle = particleExclusions.createChildNode("ParticleExclusion");
        particle.setIntProperty("size", exclusions.size() );
        particle.setIntProperty("index", i );
        for (unsigned int jj = 0; jj < exclusions.size(); jj++) {
            particle.createChildNode("ParticleExclusion").setIntProperty( "ex", exclusions[jj] );
        }
    }
}

void* AmoebaVdwForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    AmoebaVdwForce* force = new AmoebaVdwForce();
    try {

        force->setSigmaCombiningRule( node.getStringProperty( "SigmaCombiningRule" ) );
        force->setEpsilonCombiningRule( node.getStringProperty( "EpsilonCombiningRule" ) );
        force->setCutoff( node.getDoubleProperty( "VdwCutoff" ) );
        force->setUseNeighborList( node.getIntProperty( "VdwUseNeighborList" ) );
        force->setPBC( node.getIntProperty( "VdwPBC" ) );

        const SerializationNode& particles = node.getChildNode("VdwParticles");
        for (int i = 0; i < (int) particles.getChildren().size(); i++) {
            const SerializationNode& particle = particles.getChildren()[i];
            force->addParticle(particle.getIntProperty("ivIndex"), particle.getIntProperty("classIndex"), particle.getDoubleProperty("sig"), particle.getDoubleProperty("eps"), particle.getDoubleProperty("red"));
        }

        // exclusions

        const SerializationNode& particleExclusions = node.getChildNode("VdwParticleExclusions");
        for (int i = 0; i < (int) particleExclusions.getChildren().size(); i++) {
            const SerializationNode& particleExclusion = particleExclusions.getChildren()[i];
            std::vector< int > exclusions;
            for (unsigned int jj = 0; jj < particleExclusion.getChildren().size(); jj++) {
                exclusions.push_back( particleExclusion.getChildren()[jj].getIntProperty("ex") );
            }
            int index = particleExclusion.getIntProperty( "index" );
            force->setParticleExclusions( index, exclusions );
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
