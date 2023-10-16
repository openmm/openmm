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
    void getPerturbationEnergy(ContextImpl& context, double& u1, double& u0, double& energy);

private:
    const ATMForce& owner;
    Kernel kernel;
    System innerSystem0, innerSystem1;
    VerletIntegrator innerIntegrator0, innerIntegrator1;
    Context *innerContext0, *innerContext1;
    Lepton::CompiledExpression energyExpression, u0DerivExpression, u1DerivExpression;
    double state0Energy, state1Energy, combinedEnergy;
    std::vector<std::string> globalParameterNames, paramDerivNames;
    std::vector<double> globalValues;
    std::vector<Lepton::CompiledExpression> paramDerivExpressions;
    void copySystem(ContextImpl& context, const System& system, System& innerSystem);
};

} // namespace OpenMM

#endif /*OPENMM_ATMFORCEIMPL_H_*/
