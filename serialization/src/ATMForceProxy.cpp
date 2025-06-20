/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2023 Stanford University and the Authors.      *
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

#include "openmm/serialization/ATMForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/ATMForce.h"
#include "openmm/OpenMMException.h"
#include <vector>
#include <sstream>

using namespace OpenMM;
using namespace std;

ATMForceProxy::ATMForceProxy() : SerializationProxy("ATMForce") {
}

void ATMForceProxy::loadParams(int numParticles, const ATMForce& force, vector<int>& type, vector<Vec3>& d1, vector<Vec3>& d0, vector<int>& j1, vector<int>& i1, vector<int>& j0, vector<int>& i0) const {
      for (int p = 0; p < numParticles; p++) {
	const ATMTransformation* transformation = force.getParticleTransformation(p);
	if (transformation != NULL) {
	    string s = transformation->getName();
	    if (s == "FixedDisplacement") {
		type[p] = ATMTransformationType.at(s);
		d1[p] = dynamic_cast<const ATMFixedDisplacement*>(transformation)->getFixedDisplacement1();
		d0[p] = dynamic_cast<const ATMFixedDisplacement*>(transformation)->getFixedDisplacement0();
		j1[p] = i1[p] = j0[p] = i0[p] = -1;
	    }
	    else if (s == "VectordistanceDisplacement") {
		type[p] = ATMTransformationType.at(s);
		d1[p] = Vec3(0,0,0);
		d0[p] = Vec3(0,0,0);
		j1[p] = dynamic_cast<const ATMVectordistanceDisplacement*>(transformation)->getDestinationParticle1();
		i1[p] = dynamic_cast<const ATMVectordistanceDisplacement*>(transformation)->getOriginParticle1();
		j0[p] = dynamic_cast<const ATMVectordistanceDisplacement*>(transformation)->getDestinationParticle0();
		i0[p] = dynamic_cast<const ATMVectordistanceDisplacement*>(transformation)->getOriginParticle0();
	    }
	    else {
		throw OpenMMException("loadParams(): invalid particle Transformation");
	    }
	}
	else {
	    type[p] = ATMTransformationType.at("NoDisplacement");
	    d1[p] = Vec3(0,0,0);
	    d0[p] = Vec3(0,0,0);
	    j1[p] = i1[p] = j0[p] = i0[p] = -1;
	}
    }
}

void ATMForceProxy::storeParams(int numParticles, ATMForce& force, const SerializationNode& particles) const {
    for (auto& p : particles.getChildren()){

	//support older serialized ATMForces that did not store the transformation type
	//or the vector distance indexes
	int transformation_type;
        if (p.hasProperty("type")) {
	    transformation_type = p.getIntProperty("type");
	}
	else if (p.hasProperty("pj1")) {
	    if (p.getIntProperty("pj1") >= 0) {
		transformation_type = ATMTransformationType.at("VectordistanceDisplacement");
	    }
	    else {
		transformation_type = ATMTransformationType.at("FixedDisplacement");
	    }
	}
	else {
	    transformation_type = ATMTransformationType.at("FixedDisplacement");
	}

	//add particle based on transformation type
	if (transformation_type == ATMTransformationType.at("NoDisplacement")) {
	    force.addParticle();
	}
	else if (transformation_type == ATMTransformationType.at("FixedDisplacement")) {
	    force.addParticle(new ATMFixedDisplacement(Vec3(p.getDoubleProperty("d1x"), p.getDoubleProperty("d1y"), p.getDoubleProperty("d1z")), Vec3(p.getDoubleProperty("d0x"), p.getDoubleProperty("d0y"), p.getDoubleProperty("d0z"))));
	}
	else if (transformation_type == ATMTransformationType.at("VectordistanceDisplacement")) {
	    force.addParticle(new ATMVectordistanceDisplacement( p.getIntProperty("pj1"), p.getIntProperty("pi1"), p.getIntProperty("pj0"), p.getIntProperty("pi0")));
	}
	else {
	    throw OpenMMException("storeParams(): invalid particle Transformation");
	}

    }
}

void ATMForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 0);
    const ATMForce& force = *reinterpret_cast<const ATMForce*> (object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());
    node.setStringProperty("energy", force.getEnergyFunction());
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& energyDerivs = node.createChildNode("EnergyParameterDerivatives");
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        energyDerivs.createChildNode("Parameter").setStringProperty("name", force.getEnergyParameterDerivativeName(i));
    }
    SerializationNode& forces = node.createChildNode("Forces");
    for (int i = 0; i < force.getNumForces(); i++) {
        SerializationNode& f = forces.createChildNode("Force");
        f.createChildNode("Force", &force.getForce(i));
    }
    SerializationNode& particles = node.createChildNode("Particles");

    int numParticles = force.getNumParticles();
    vector<int> type(numParticles);
    vector<int> j1(numParticles);
    vector<int> i1(numParticles);
    vector<int> j0(numParticles);
    vector<int> i0(numParticles);
    vector<Vec3> d1(numParticles);
    vector<Vec3> d0(numParticles);
    loadParams(numParticles, force, type, d1, d0, j1, i1, j0, i0);

    for (int i = 0; i < force.getNumParticles(); i++) {
        particles.createChildNode("Particle").setIntProperty("type", type[i])
	  .setDoubleProperty("d1x", d1[i][0]).setDoubleProperty("d1y", d1[i][1]).setDoubleProperty("d1z", d1[i][2])
	  .setDoubleProperty("d0x", d0[i][0]).setDoubleProperty("d0y", d0[i][1]).setDoubleProperty("d0z", d0[i][2])
	  .setIntProperty("pj1", j1[i]).setIntProperty("pi1", i1[i]).setIntProperty("pj0", j0[i]).setIntProperty("pi0", i0[i]);
    }
}

void* ATMForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 0)
        throw OpenMMException("Unsupported version number");
    ATMForce* force = NULL;
    try {
        ATMForce* force = new ATMForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        const SerializationNode& energyDerivs = node.getChildNode("EnergyParameterDerivatives");
        for (auto& parameter : energyDerivs.getChildren())
            force->addEnergyParameterDerivative(parameter.getStringProperty("name"));
        const SerializationNode& forces = node.getChildNode("Forces");
        for (auto& f : forces.getChildren())
            force->addForce(f.getChildren()[0].decodeObject<Force>());
        const SerializationNode& particles = node.getChildNode("Particles");

	storeParams(force->getNumParticles(), *force, particles);

        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}
