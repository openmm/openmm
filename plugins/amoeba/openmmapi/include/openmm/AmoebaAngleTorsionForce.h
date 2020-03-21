#ifndef OPENMM_AMOEBA_ANGLE_TORSION_FORCE_H_
#define OPENMM_AMOEBA_ANGLE_TORSION_FORCE_H_

#include "openmm/Force.h"
#include "internal/windowsExportAmoeba.h"
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT_AMOEBA AmoebaAngleTorsionForce : public Force {

public:
	AmoebaAngleTorsionForce();

	int getNumAngleTorsions() const {
		return angleTorsions.size();
	}

	int addAngleTorsion(int particle1, int particle2, int particle3, int particle4, double angleCBA, double angleDCB,
			double k1, double k2, double k3, double k4, double k5, double k6);

	void getAngleTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& angleCBA, double& angleDCB,
			double& k1, double& k2, double& k3, double& k4, double& k5, double& k6) const;

	void setAngleTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, double angleCBA, double angleDCB,
			double k1, double k2, double k3, double k4, double k5, double k6);

	void updateParametersInContext(Context& context);

	bool usesPeriodicBoundaryConditions() const {
		return false;
	}

protected:
	ForceImpl* createImpl() const;

private:
	class AngleTorsionInfo;
	std::vector<AngleTorsionInfo> angleTorsions;
};

/**
 * This is an internal class used to record information about a angle-torsion.
 * @private
 */

class AmoebaAngleTorsionForce::AngleTorsionInfo {

public:
	int particle1, particle2, particle3, particle4;
	double angleCBA, angleDCB, k1, k2, k3, k4, k5, k6;
	AngleTorsionInfo() {
		particle1 = particle2 = particle3 = particle4 = -1;
		angleCBA = angleDCB = k1 = k2 = k3 = k4 = k5 = k6 = 0.0;
	}

	AngleTorsionInfo(int particle1, int particle2, int particle3, int particle4, double angleCBA, double angleDCB,
			double k1, double k2, double k3, double k4, double k5, double k6) : 
		particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4),
		angleCBA(angleCBA), angleDCB(angleDCB),
		k1(k1), k2(k2), k3(k3), k4(k4), k5(k5), k6(k6) {
	}

};

} // namespace OpenMM
#endif
