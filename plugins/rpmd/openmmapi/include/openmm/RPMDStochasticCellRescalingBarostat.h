#ifndef OPENMM_RPMD_STOCHASTIC_CELL_RESCALING_BAROSTAT_H_
#define OPENMM_RPMD_STOCHASTIC_CELL_RESCALING_BAROSTAT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Stochastic cell rescaling barostat for RPMD (Bernetti & Bussi, JCP 2020;
 * extends Bussi, Zykova-Timan & Parrinello, JCP 130, 074101, 2009).
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "openmm/internal/windowsExportRpmd.h"
#include <string>

namespace OpenMM {

/**
 * Stochastic cell rescaling (SCR) barostat for use with RPMDIntegrator.
 * Unlike RPMDMonteCarloBarostat, volume evolves by a stochastic differential
 * equation without Metropolis acceptance.
 */
class OPENMM_EXPORT_RPMD RPMDStochasticCellRescalingBarostat : public Force {
public:
    static const std::string& Pressure() {
        static const std::string key = "RPMDSCRPressure";
        return key;
    }
    /**
     * @param defaultPressure target pressure (bar)
     * @param tauP relaxation time (ps)
     * @param frequency apply barostat every this many integrator steps
     */
    RPMDStochasticCellRescalingBarostat(double defaultPressure, double tauP = 1.0, int frequency = 1);
    double getDefaultPressure() const {
        return defaultPressure;
    }
    void setDefaultPressure(double pressure);
    double getTauP() const {
        return tauP;
    }
    void setTauP(double t);
    /** Isothermal compressibility (1/bar). Default ~4.5e-5 for water. */
    double getIsothermalCompressibility() const {
        return isothermalCompressibility;
    }
    void setIsothermalCompressibility(double betaT);
    int getFrequency() const {
        return frequency;
    }
    void setFrequency(int freq);
    bool getScaleVelocities() const {
        return scaleVelocities;
    }
    void setScaleVelocities(bool scale) {
        scaleVelocities = scale;
    }
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
protected:
    ForceImpl* createImpl() const;
private:
    double defaultPressure;
    double tauP;
    double isothermalCompressibility;
    int frequency;
    bool scaleVelocities;
    int randomNumberSeed;
};

} // namespace OpenMM

#endif
