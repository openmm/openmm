#ifndef OPENMM_ATMFORCEFORCEIMPL_H_
#define OPENMM_ATMFORCEFORCEIMPL_H_


#include "openmm/ATMForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/windowsExport.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

//class System;

/**
 * This is the internal implementation of ATMForce.
 */

//class OPENMM_EXPORT_ATMFORCE ATMForceImpl : public OpenMM::ForceImpl {
class OPENMM_EXPORT ATMForceImpl : public OpenMM::ForceImpl {
public:
    ATMForceImpl(const ATMForce& owner);
    ~ATMForceImpl();
    void initialize(OpenMM::ContextImpl& context);
    const ATMForce& getOwner() const {
        return owner;
    }
    void updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    //{
    //    return std::map<std::string, double>(); // This force field doesn't define any parameters.
    //}
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(OpenMM::ContextImpl& context);
    double getPerturbationEnergy() const {
       return PerturbationEnergy;
    }

private:
    const ATMForce& owner;
    OpenMM::Kernel kernel;
    OpenMM::System innerSystem1, innerSystem2;
    OpenMM::VerletIntegrator innerIntegrator1, innerIntegrator2;
    OpenMM::Context *innerContext1, *innerContext2;
    double PerturbationEnergy;
    bool hasInitializedInnerContexts;
    int variable_force_groups_mask;
    void copysystem(const OpenMM::System& system, OpenMM::System& innerSystem);
};

} // namespace OpenMM

#endif /*OPENMM_ATMFORCEIMPL_H_*/
