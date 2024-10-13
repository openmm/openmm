//
// Created by babaid on 05.10.24.
//

#include "openmm/ExternalPuremdForce.h"
#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ExternalPuremdForceImpl.h"

#include <cstring>
using namespace OpenMM;
using namespace std;

ExternalPuremdForce::ExternalPuremdForce() {}

ExternalPuremdForce::ExternalPuremdForce(const std::string& ffieldFile, const std::string& controlFile) :usePeriodic(false), numContexts(0), ffield_file(ffieldFile), control_file(controlFile) {
}

int ExternalPuremdForce::addAtom(int particle, char* symbol, bool isQM) {
    bool hasLen = std::strlen(symbol) > 1;
    allAtoms.push_back(particle);
    allIsQM.push_back(isQM);
    if (hasLen) {
      allSymbols.push_back(symbol[1]);
    }
    else
    {
      allSymbols.push_back('\0');
    }
    return allAtoms.size();
}

void ExternalPuremdForce::getParticleParameters(int index, int &particle, char* symbol, int &isQM)  const {
    particle = allAtoms[index];
    symbol[0] = allSymbols[index*2];
    symbol[1] = allSymbols[index*2 + 1];
    isQM = allIsQM[index];
}

ForceImpl*ExternalPuremdForce::createImpl() const {
    if (numContexts == 0) {
        // Begin tracking changes to atoms.
        firstChangedBond = allAtoms.size();
        lastChangedBond = -1;
    }
    numContexts++;
    return new ExternalPuremdForceImpl(*this);
}


