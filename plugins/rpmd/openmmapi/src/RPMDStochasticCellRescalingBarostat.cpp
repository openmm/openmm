#include "openmm/RPMDStochasticCellRescalingBarostat.h"
#include "openmm/internal/RPMDStochasticCellRescalingBarostatImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

RPMDStochasticCellRescalingBarostat::RPMDStochasticCellRescalingBarostat(
        double defaultPressure, double tauP, int frequency) {
    setDefaultPressure(defaultPressure);
    setTauP(tauP);
    setFrequency(frequency);
    setIsothermalCompressibility(4.5e-5);
    setScaleVelocities(true);
    setRandomNumberSeed(0);
}

void RPMDStochasticCellRescalingBarostat::setDefaultPressure(double pressure) {
    defaultPressure = pressure;
}

void RPMDStochasticCellRescalingBarostat::setTauP(double t) {
    if (t <= 0.0)
        throw OpenMMException("RPMDStochasticCellRescalingBarostat: tauP must be positive");
    tauP = t;
}

void RPMDStochasticCellRescalingBarostat::setIsothermalCompressibility(double betaT) {
    if (betaT <= 0.0)
        throw OpenMMException("RPMDStochasticCellRescalingBarostat: compressibility must be positive");
    isothermalCompressibility = betaT;
}

void RPMDStochasticCellRescalingBarostat::setFrequency(int freq) {
    if (freq <= 0)
        throw OpenMMException("RPMDStochasticCellRescalingBarostat: frequency must be positive");
    frequency = freq;
}

ForceImpl* RPMDStochasticCellRescalingBarostat::createImpl() const {
    return new RPMDStochasticCellRescalingBarostatImpl(*this);
}
