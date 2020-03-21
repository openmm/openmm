#ifndef OPENMM_AMOEBA_ANGLE_TORSION_FORCE_IMPL_H_
#define OPENMM_AMOEBA_ANGLE_TORSION_FORCE_IMPL_H_

#include "openmm/internal/ForceImpl.h"
#include "openmm/AmoebaAngleTorsionForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class AmoebaAngleTorsionForceImpl : public ForceImpl {

public:
	AmoebaAngleTorsionForceImpl(const AmoebaAngleTorsionForce& owner);
	~AmoebaAngleTorsionForceImpl();
	void initialize(ContextImpl& context);
	const AmoebaAngleTorsionForce& getOwner() const {
		return owner;
	}
	void updateContextState(ContextImpl& context) {
	}
	double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
	std::map<std::string, double> getDefaultParameters() {
		return std::map<std::string, double>();
	}
	std::vector<std::string> getKernelNames();
	void updateParametersInContext(ContextImpl& context);
protected:
	const AmoebaAngleTorsionForce& owner;
	Kernel kernel;
};

} // namespace OpenMM

#endif
