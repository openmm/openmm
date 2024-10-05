//
// Created by babaid on 05.10.24.
//

#include "openmm/ExternalPuremdForce.h"
#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ExternalPuremdForceImpl.h"

using namespace OpenMM;
using namespace std;

ExternalPuremdForce::ExternalPuremdForce() : usePeriodic(false), numContexts(0) {
}

int ExternalPuremdForce::addAtom(int particle) {
    atoms.push_back(AtomInfo(particle));
    return atoms.size()-1;
}

void ExternalPuremdForce::getAtom(int index, int &particle) const {
    ASSERT_VALID_INDEX(index, atoms)
    particle = atoms[index].particle;
}

ForceImpl*ExternalPuremdForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to atoms.
        firstChangedBond = atoms.size();
        lastChangedBond = -1;
    }
    numContexts++;
    return new ExternalPuremdForceImpl(*this);
}

void ExternalPuremdForce::updateParametersInContext(Context& context) {
    dynamic_cast<ExternalPuremdForceImpl &>(getImplInContext(context)).updateParametersInContext(getContextImpl(context), firstChangedBond, lastChangedBond);
    if (numContexts == 1) {
        // We just updated the only existing context for this force, so we can reset
        // the tracking of changed atoms.
        firstChangedBond = atoms.size();
        lastChangedBond = -1;
    }
}

