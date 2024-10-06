//
// Created by babaid on 05.10.24.
//

#include "openmm/internal/ExternalPuremdForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

ExternalPuremdForceImpl::ExternalPuremdForceImpl(const ExternalPuremdForce & owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

ExternalPuremdForceImpl::~ExternalPuremdForceImpl() {
}

void ExternalPuremdForceImpl::initialize(ContextImpl& context) {
    const System& system = context.getSystem();
    for (int i = 0; i < owner.getNumAtoms(); i++) {
        int particle;
        if (particle < 0 || particle >= system.getNumParticles())
        {
                stringstream msg;
                msg << "HarmonicBondForce: Illegal particle index for a bond: ";
                msg << particle;
                throw OpenMMException(msg.str());
        }
    }
    kernel = context.getPlatform().createKernel(CalcHarmonicBondForceKernel::Name(), context);
    kernel.getAs<CalcExternalPuremdForceKernel>().initialize(context.getSystem(), owner);
}

double ExternalPuremdForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcExternalPuremdForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> ExternalPuremdForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcHarmonicBondForceKernel::Name());
    return names;
}

vector<int> ExternalPuremdForceImpl::getSimulatedParticles() const {
    int numBonds = owner.getNumAtoms();
    vector<int> atoms(numBonds);
    vector<char> symbols(numBonds);
    vector<int> isQM(numBonds);
    for (int i = 0; i < numBonds; i++) {
        owner.getParticleParameters(i, atoms[i], symbols[i], isQM[i]);
    }
    return atoms;
}

void ExternalPuremdForceImpl::updateParametersInContext(ContextImpl& context, int firstBond, int lastBond) {
    kernel.getAs<CalcExternalPuremdForceKernel>().copyParametersToContext(context, owner, firstBond, lastBond);
    context.systemChanged();
}
