#ifndef OPENMM_RPMD_STOCHASTIC_CELL_RESCALING_BAROSTATIMPL_H_
#define OPENMM_RPMD_STOCHASTIC_CELL_RESCALING_BAROSTATIMPL_H_

#include "openmm/RPMDStochasticCellRescalingBarostat.h"
#include "openmm/RPMDUpdater.h"
#include "openmm/Kernel.h"
#include "openmm/Vec3.h"
#include "sfmt/SFMT.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

class RPMDStochasticCellRescalingBarostatImpl : public RPMDUpdater {
public:
    RPMDStochasticCellRescalingBarostatImpl(const RPMDStochasticCellRescalingBarostat& owner);
    void initialize(ContextImpl& context);
    const RPMDStochasticCellRescalingBarostat& getOwner() const {
        return owner;
    }
    void updateRPMDState(ContextImpl& context);
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
private:
    const RPMDStochasticCellRescalingBarostat& owner;
    int step;
    double engint;
    OpenMM_SFMT::SFMT random;
    Kernel kernel;
};

} // namespace OpenMM

#endif
