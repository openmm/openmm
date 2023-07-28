#ifndef OPENMM_ATMFORCEFORCEIMPL_H_
#define OPENMM_ATMFORCEFORCEIMPL_H_


#include "openmm/ATMForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/windowsExport.h"
#include "lepton/CompiledExpression.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of ATMForce.
 */

class OPENMM_EXPORT ATMForceImpl : public ForceImpl {
public:
    ATMForceImpl(const ATMForce& owner);
    ~ATMForceImpl();
    void initialize(ContextImpl& context);
    const ATMForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    std::vector<std::pair<int, int> > getBondedParticles() const;
    void updateParametersInContext(ContextImpl& context);
    void getPerturbationEnergy(double& u0, double& u1, double& energy) const;

private:
    const ATMForce& owner;
    Kernel kernel;
    System innerSystem1, innerSystem2;
    VerletIntegrator innerIntegrator1, innerIntegrator2;
    Context *innerContext1, *innerContext2;
    Lepton::CompiledExpression energyExpression, u0DerivExpression, u1DerivExpression;
    double state1Energy, state2Energy, combinedEnergy;
    std::vector<std::string> globalParameterNames;
    std::vector<double> globalValues;
    void copySystem(ContextImpl& context, const System& system, System& innerSystem);
};

} // namespace OpenMM

#endif /*OPENMM_ATMFORCEIMPL_H_*/
