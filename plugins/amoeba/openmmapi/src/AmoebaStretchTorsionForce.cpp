#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/AmoebaStretchTorsionForce.h"
#include "openmm/internal/AmoebaStretchTorsionForceImpl.h"

using namespace OpenMM;

AmoebaStretchTorsionForce::AmoebaStretchTorsionForce() {
}

int AmoebaStretchTorsionForce::addStretchTorsion(int particle1, int particle2, int particle3, int particle4, double lengthBA, double lengthCB, double lengthDC,
		                                 double k1, double k2, double k3, double k4, double k5, double k6, double k7, double k8, double k9) {
    stretchTorsions.push_back(StretchTorsionInfo(particle1, particle2, particle3, particle4, lengthBA, lengthCB, lengthDC,
			                         k1, k2, k3, k4, k5, k6, k7, k8, k9));
    return stretchTorsions.size() - 1;   
}

void AmoebaStretchTorsionForce::getStretchTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& lengthBA, double& lengthCB, double& lengthDC,
		                                            double& k1, double& k2, double& k3, double& k4, double& k5, double& k6, double& k7, double& k8, double& k9) const {

	particle1 = stretchTorsions[index].particle1;
	particle2 = stretchTorsions[index].particle2;
	particle3 = stretchTorsions[index].particle3;
	particle4 = stretchTorsions[index].particle4;
	lengthBA = stretchTorsions[index].lengthBA;
	lengthCB = stretchTorsions[index].lengthCB;
	lengthDC = stretchTorsions[index].lengthDC;
	k1 = stretchTorsions[index].k1;
	k2 = stretchTorsions[index].k2;
	k3 = stretchTorsions[index].k3;
	k4 = stretchTorsions[index].k4;
	k5 = stretchTorsions[index].k5;
	k6 = stretchTorsions[index].k6;
	k7 = stretchTorsions[index].k7;
	k8 = stretchTorsions[index].k8;
	k9 = stretchTorsions[index].k9;
}

void AmoebaStretchTorsionForce::setStretchTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, double lengthBA, double lengthCB, double lengthDC,
			                                    double k1, double k2, double k3, double k4, double k5, double k6, double k7, double k8, double k9) {
        stretchTorsions[index].particle1 = particle1;
        stretchTorsions[index].particle2 = particle2;
        stretchTorsions[index].particle3 = particle3;
        stretchTorsions[index].particle4 = particle4;
        stretchTorsions[index].lengthBA = lengthBA;
        stretchTorsions[index].lengthCB = lengthCB;
        stretchTorsions[index].lengthDC = lengthDC;
        stretchTorsions[index].k1 = k1;
        stretchTorsions[index].k2 = k2;
        stretchTorsions[index].k3 = k3;
        stretchTorsions[index].k4 = k4;
        stretchTorsions[index].k5 = k5;
        stretchTorsions[index].k6 = k6;
        stretchTorsions[index].k7 = k7;
        stretchTorsions[index].k8 = k8;
        stretchTorsions[index].k9 = k9;

}

void AmoebaStretchTorsionForce::updateParametersInContext(Context& context) {
	dynamic_cast<AmoebaStretchTorsionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

ForceImpl* AmoebaStretchTorsionForce::createImpl() const {
	return new AmoebaStretchTorsionForceImpl(*this);
}
