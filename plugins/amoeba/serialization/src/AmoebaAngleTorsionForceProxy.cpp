#include "openmm/serialization/AmoebaAngleTorsionForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "openmm/AmoebaAngleTorsionForce.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

AmoebaAngleTorsionForceProxy::AmoebaAngleTorsionForceProxy() : SerializationProxy("AmoebaAngleTorionForce") {
}

void AmoebaAngleTorsionForceProxy::serialize(const void* object, SerializationNode& node) const {
	node.setIntProperty("version", 1);
	const AmoebaAngleTorsionForce& force = *reinterpret_cast<const AmoebaAngleTorsionForce*>(object);
	node.setIntProperty("forceGroup", force.getForceGroup());
	SerializationNode& bonds = node.createChildNode("AngleTorsions");
	for (unsigned int ii = 0; ii < static_cast<unsigned int>(force.getNumAngleTorsions()); ii++) {
		int particle1, particle2, particle3, particle4;
		double angleCBA, angleDCB, k1, k2, k3, k4, k5, k6;
		force.getAngleTorsionParameters(ii, particle1, particle2, particle3, particle4, 
				angleCBA, angleDCB, k1, k2, k3, k4, k5, k6);
		bonds.createChildNode("AngleTorsion").setIntProperty("p1",particle1).setIntProperty("p2",particle2).setIntProperty("p3",particle3).setIntProperty("p4",particle4).
			setDoubleProperty("angleCBA",angleCBA).setDoubleProperty("angleDCB",angleDCB).
			setDoubleProperty("k1",k1).setDoubleProperty("k2",k2).setDoubleProperty("k3",k3).setDoubleProperty("k4",k4).setDoubleProperty("k5",k5).
			setDoubleProperty("k6",k6);
	}

}

void* AmoebaAngleTorsionForceProxy::deserialize(const SerializationNode& node) const {
	int version = node.getIntProperty("version");
	if (version != 1)
		throw OpenMMException("Unsupported version number");
	AmoebaAngleTorsionForce* force = new AmoebaAngleTorsionForce();
	try {
		force->setForceGroup(node.getIntProperty("forceGroup", 0));
		const SerializationNode& bonds = node.getChildNode("AngleTorsions");
		for (unsigned int ii = 0; ii < (int) bonds.getChildren().size(); ii++) {
			const SerializationNode& bond = bonds.getChildren()[ii];
			force->addAngleTorsion(bond.getIntProperty("p1"),bond.getIntProperty("p2"), bond.getIntProperty("p3"), bond.getIntProperty("p4"),
					bond.getDoubleProperty("angleCBA"),bond.getDoubleProperty("angleDCB"),
					bond.getDoubleProperty("k1"), bond.getDoubleProperty("k2"), bond.getDoubleProperty("k3"),
					bond.getDoubleProperty("k4"), bond.getDoubleProperty("k5"), bond.getDoubleProperty("k6"));
		}
	}

	catch (...) {
		delete force;
		throw;
	}
	return force;
}
