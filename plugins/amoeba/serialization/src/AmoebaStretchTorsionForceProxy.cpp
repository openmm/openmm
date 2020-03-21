#include "openmm/serialization/AmoebaStretchTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaStretchTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaStretchTorsionForceProxy::AmoebaStretchTorsionForceProxy() : SerializationProxy("AmoebaStretchTorionForce") {
}

void AmoebaStretchTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
	node.setIntProperty("version", 1);
	const AmoebaStretchTorsionForce& force = *reinterpret_cast<const AmoebaStretchTorsionForce*>(object);
	node.setIntProperty("forceGroup", force.getForceGroup());
	SerializationNode& bonds = node.createChildNode("StretchTorsions");
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumStretchTorsions()); ii++) {
		int particle1, particle2, particle3, particle4;
		double lengthBA, lengthCB, lengthDC, k1, k2, k3, k4, k5, k6, k7, k8, k9;
		force.getStretchTorsionParameters(ii, particle1, particle2, particle3, particle4,
				lengthBA, lengthCB, lengthDC, k1, k2, k3, k4, k5, k6, k7, k8, k9);
		bonds.createChildNode("StretchTorsion").setIntProperty("p1", particle1).setIntProperty("p2", particle2).setIntProperty("p3", particle3).setIntProperty("p4", particle4).
			setDoubleProperty("lengthBA", lengthBA).setDoubleProperty("lengthCB", lengthCB).setDoubleProperty("lengthDC", lengthDC).
			setDoubleProperty("k1",k1).setDoubleProperty("k2",k2).setDoubleProperty("k3",k3).setDoubleProperty("k4",k4).setDoubleProperty("k5",k5).setDoubleProperty("k6",k6).
			setDoubleProperty("k7",k7).setDoubleProperty("k8",k8).setDoubleProperty("k9",k9);
	}

}

void* AmoebaStretchTorsionForceProxy::deserialize(const SerializationNode& node) const {
	int version = node.getIntProperty("version");
	if (version != 1)
		throw OpenMMException("Unsupported version number");
	AmoebaStretchTorsionForce* force = new AmoebaStretchTorsionForce();
	try {
		force->setForceGroup(node.getIntProperty("forceGroup", 0));
		const SerializationNode& bonds = node.getChildNode("StretchTorsions");
		for (unsigned int ii = 0; ii < (int) bonds.getChildren().size(); ii++) {
			const SerializationNode& bond = bonds.getChildren()[ii];
			force->addStretchTorsion(bond.getIntProperty("p1"), bond.getIntProperty("p2"), bond.getIntProperty("p3"), bond.getIntProperty("p4"),
					bond.getDoubleProperty("lengthBA"), bond.getDoubleProperty("lengthCB"), bond.getDoubleProperty("lengthDC"),
					bond.getDoubleProperty("k1"), bond.getDoubleProperty("k2"), bond.getDoubleProperty("k3"), 
					bond.getDoubleProperty("k4"), bond.getDoubleProperty("k5"), bond.getDoubleProperty("k6"),
					bond.getDoubleProperty("k7"), bond.getDoubleProperty("k8"), bond.getDoubleProperty("k9"));
		}
	}
	catch (...) {
		delete force;
		throw;
	}

	return force;
}
