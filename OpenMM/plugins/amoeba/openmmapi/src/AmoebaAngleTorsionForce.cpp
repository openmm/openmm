#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/AmoebaAngleTorsionForce.h"
#include "openmm/internal/AmoebaAngleTorsionForceImpl.h"

using namespace OpenMM;

AmoebaAngleTorsionForce::AmoebaAngleTorsionForce() {
}

int AmoebaAngleTorsionForce::addAngleTorsion(int particle1, int particle2, int particle3, int particle4, double angleCBA, double angleDCB,
		double k1, double k2, double k3, double k4, double k5, double k6) {
	angleTorsions.push_back(AngleTorsionInfo(particle1, particle2, particle3, particle4, angleCBA, angleDCB,
				k1, k2, k3, k4, k5, k6));
	return angleTorsions.size() - 1;
}

void AmoebaAngleTorsionForce::getAngleTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& angleCBA, double& angleDCB,
		double& k1, double& k2, double& k3, double& k4, double& k5, double& k6) const {
	particle1 = angleTorsions[index].particle1;
	particle2 = angleTorsions[index].particle2;
	particle3 = angleTorsions[index].particle3;
	particle4 = angleTorsions[index].particle4;
	angleCBA = angleTorsions[index].angleCBA;
	angleDCB = angleTorsions[index].angleDCB;
	k1 = angleTorsions[index].k1;
	k2 = angleTorsions[index].k2;
	k3 = angleTorsions[index].k3;
	k4 = angleTorsions[index].k4;
	k5 = angleTorsions[index].k5;
	k6 = angleTorsions[index].k6;
}

void AmoebaAngleTorsionForce::setAngleTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, double angleCBA, double angleDCB,
		double k1, double k2, double k3, double k4, double k5, double k6) {
	angleTorsions[index].particle1 = particle1;
	angleTorsions[index].particle2 = particle2;
	angleTorsions[index].particle3 = particle3;
	angleTorsions[index].particle4 = particle4;
	angleTorsions[index].angleCBA = angleCBA;
	angleTorsions[index].angleDCB = angleDCB;
	angleTorsions[index].k1 = k1;
	angleTorsions[index].k2 = k2;
	angleTorsions[index].k3 = k3;
	angleTorsions[index].k4 = k4;
	angleTorsions[index].k5 = k5;
	angleTorsions[index].k6 = k6;
}

void AmoebaAngleTorsionForce::updateParametersInContext(Context& context) {
	dynamic_cast<AmoebaAngleTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* AmoebaAngleTorsionForce::createImpl() const {
	return new AmoebaAngleTorsionForceImpl(*this);
}
