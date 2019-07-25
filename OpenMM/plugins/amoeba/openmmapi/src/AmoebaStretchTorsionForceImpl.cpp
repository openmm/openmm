#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaStretchTorsionForceImpl.h"
#include "openmm/amoebaKernels.h"

using namespace OpenMM;

using std::pair;
using std::vector;
using std::set;

AmoebaStretchTorsionForceImpl::AmoebaStretchTorsionForceImpl(const AmoebaStretchTorsionForce& owner) : owner(owner) {
}

AmoebaStretchTorsionForceImpl::~AmoebaStretchTorsionForceImpl() {
}

void AmoebaStretchTorsionForceImpl::initialize(ContextImpl& context) {
	kernel = context.getPlatform().createKernel(CalcAmoebaStretchTorsionForceKernel::Name(), context);
	kernel.getAs<CalcAmoebaStretchTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaStretchTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
	if ((groups&(1<<owner.getForceGroup())) != 0)
		return kernel.getAs<CalcAmoebaStretchTorsionForceKernel>().execute(context, includeForces, includeEnergy);
	return 0.0;
}

std::vector<std::string> AmoebaStretchTorsionForceImpl::getKernelNames() {
	std::vector<std::string> names;
	names.push_back(CalcAmoebaStretchTorsionForceKernel::Name());
	return names;
}

void AmoebaStretchTorsionForceImpl::updateParametersInContext(ContextImpl& context) {
	kernel.getAs<CalcAmoebaStretchTorsionForceKernel>().copyParametersToContext(context, owner);
}
