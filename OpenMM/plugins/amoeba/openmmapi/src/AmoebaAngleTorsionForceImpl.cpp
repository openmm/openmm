#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaAngleTorsionForceImpl.h"
#include "openmm/amoebaKernels.h"

using namespace OpenMM;

using std::pair;
using std::vector;
using std::set;

AmoebaAngleTorsionForceImpl::AmoebaAngleTorsionForceImpl(const AmoebaAngleTorsionForce& owner) : owner(owner) {
}

AmoebaAngleTorsionForceImpl::~AmoebaAngleTorsionForceImpl() {
}

void AmoebaAngleTorsionForceImpl::initialize(ContextImpl& context) {
	kernel = context.getPlatform().createKernel(CalcAmoebaAngleTorsionForceKernel::Name(), context);
	kernel.getAs<CalcAmoebaAngleTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaAngleTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
	if ((groups&(1<<owner.getForceGroup())) != 0)
		return kernel.getAs<CalcAmoebaAngleTorsionForceKernel>().execute(context, includeForces, includeEnergy);
	return 0.0;
}

std::vector<std::string> AmoebaAngleTorsionForceImpl::getKernelNames() {
	std::vector<std::string> names;
	names.push_back(CalcAmoebaAngleTorsionForceKernel::Name());
	return names;
}

void AmoebaAngleTorsionForceImpl::updateParametersInContext(ContextImpl& context) {
	kernel.getAs<CalcAmoebaAngleTorsionForceKernel>().copyParametersToContext(context, owner);
}
