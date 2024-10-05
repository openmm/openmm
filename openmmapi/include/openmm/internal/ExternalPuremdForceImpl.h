//
// Created by babaid on 05.10.24.
//

#ifndef OPENMM_EXTERNALPUREMDFORCEIMPL_H
#define OPENMM_EXTERNALPUREMDFORCEIMPL_H
#include "ForceImpl.h"
#include "openmm/ExternalPuremdForce.h"
#include "openmm/Kernel.h"
#include <set>
#include <string>
#include <utility>

namespace OpenMM {

/**
 * This is the internal implementation of ExternalPuremdForce.
 */


    class ExternalPuremdForceImpl : public ForceImpl {
    public:
      ExternalPuremdForceImpl(const ExternalPuremdForce & owner);
        ~ExternalPuremdForceImpl();
        void initialize(ContextImpl& context);
        const ExternalPuremdForce & getOwner() const {
            return owner;
        }
        void updateContextState(ContextImpl& context, bool& forcesInvalid) {
            // This force field doesn't update the state directly.
        }
        double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
        std::map<std::string, double> getDefaultParameters() {
            return std::map<std::string, double>(); // This force field doesn't define any parameters.
        }
        std::vector<std::string> getKernelNames();
        std::vector<int> getSimulatedParticles() const;
        void updateParametersInContext(ContextImpl& context, int firstBond, int lastBond);
    private:
        const ExternalPuremdForce & owner;
        Kernel kernel;
    };

} // namespace OpenMM
#endif //OPENMM_EXTERNALPUREMDFORCE_H
