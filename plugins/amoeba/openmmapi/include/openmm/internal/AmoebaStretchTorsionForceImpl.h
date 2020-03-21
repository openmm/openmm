#ifndef OPENMM_AMOEBA_STRETCH_TORSION_FORCE_IMPL_H_
#define OPENMM_AMOEBA_STRETCH_TORSION_FORCE_IMPL_H_

#include "openmm/internal/ForceImpl.h"
#include "openmm/AmoebaStretchTorsionForce.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class AmoebaStretchTorsionForceImpl : public ForceImpl {
public:
	AmoebaStretchTorsionForceImpl(const AmoebaStretchTorsionForce& owner);
	~AmoebaStretchTorsionForceImpl();
	void initialize(ContextImpl& context);
	const AmoebaStretchTorsionForce& getOwner() const {
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
	const AmoebaStretchTorsionForce& owner;
	Kernel kernel;
};

} // namespace OpenMM

#endif
