#ifndef OPENMM_AMOEBA_STRETCH_TORSION_FORCE_H_
#define OPENMM_AMOEBA_STRETCH_TORSION_FORCE_H_

#include "openmm/Force.h"
#include "internal/windowsExportAmoeba.h"
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT_AMOEBA AmoebaStretchTorsionForce : public Force {

public:

    AmoebaStretchTorsionForce();

    int getNumStretchTorsions() const {
	return stretchTorsions.size();
    }

    int addStretchTorsion(int particle1, int particle2, int particle3, int particle4, double lengthBA, double lengthCB, double lengthDC,
		          double k1, double k2, double k3, double k4, double k5, double k6, double k7, double k8, double k9);

    void getStretchTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4, double& lengthBA, double& lengthCB, double& lengthDC,
		                    double& k1, double& k2, double& k3, double& k4, double& k5, double& k6, double& k7, double& k8, double& k9) const;    

    void setStretchTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4, double lengthBA, double lengthCB, double lengthDC,
			             double k1, double k2, double k3, double k4, double k5, double k6, double k7, double k8, double k9);

    void updateParametersInContext(Context& context);

    bool usesPeriodicBoundaryConditions() const {
        return false;
    }

protected:
        ForceImpl* createImpl() const;
private:
    class StretchTorsionInfo;
    std::vector<StretchTorsionInfo> stretchTorsions;
};

/**
 * This is an internal class used to record information about a stretch-torsion.
 * @private
 */

class AmoebaStretchTorsionForce::StretchTorsionInfo {
public:
    int particle1, particle2, particle3, particle4;
    double lengthBA, lengthCB, lengthDC, k1, k2, k3, k4, k5, k6, k7, k8, k9;
    StretchTorsionInfo() {
	particle1 = particle2 = particle3 = particle4 = -1;  
        lengthBA = lengthCB = lengthDC = k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = k9 = 0.0;	
    }

    StretchTorsionInfo(int particle1, int particle2, int particle3, int particle4, double lengthBA, double lengthCB, double lengthDC,
		       double k1, double k2, double k3, double k4, double k5, double k6, double k7, double k8, double k9) :
	               particle1(particle1), particle2(particle2), particle3(particle3), particle4(particle4), 
		       lengthBA(lengthBA), lengthCB(lengthCB), lengthDC(lengthDC),
		       k1(k1), k2(k2), k3(k3), k4(k4), k5(k5), k6(k6), k7(k7), k8(k8), k9(k9) {
    }


};

} // namespace OpenMM

#endif /*OPENMM_AMOEBA_STRETCH_TORSION_FORCE_H_*/

