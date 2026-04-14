#include "openmm/internal/RPMDStochasticCellRescalingBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/kernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/RPMDIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include <cmath>
#include <limits>
#include <vector>

using namespace OpenMM;
using namespace OpenMM_SFMT;
using std::vector;

RPMDStochasticCellRescalingBarostatImpl::RPMDStochasticCellRescalingBarostatImpl(
        const RPMDStochasticCellRescalingBarostat& owner)
        : owner(owner), step(0), engint(0.0) {
}

void RPMDStochasticCellRescalingBarostatImpl::initialize(ContextImpl& context) {
    RPMDIntegrator* integrator = dynamic_cast<RPMDIntegrator*>(&context.getIntegrator());
    if (integrator == NULL)
        throw OpenMMException("RPMDStochasticCellRescalingBarostat must be used with RPMDIntegrator");
    if (!integrator->getApplyThermostat())
        throw OpenMMException("RPMDStochasticCellRescalingBarostat requires the integrator thermostat to be enabled");
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner, 1);
    int randSeed = owner.getRandomNumberSeed();
    if (randSeed == 0)
        randSeed = osrngseed();
    init_gen_rand(randSeed, random);
}

void RPMDStochasticCellRescalingBarostatImpl::updateRPMDState(ContextImpl& context) {
    if (++step < owner.getFrequency())
        return;
    step = 0;

    kernel.getAs<ApplyMonteCarloBarostatKernel>().synchronize(context);

    RPMDIntegrator& integrator = dynamic_cast<RPMDIntegrator&>(context.getIntegrator());
    const int numCopies = integrator.getNumCopies();
    const int numParticles = context.getSystem().getNumParticles();

    double kT = BOLTZ * integrator.getTemperature();
    double p0_bar = context.getParameter(RPMDStochasticCellRescalingBarostat::Pressure());
    double p0 = p0_bar * (AVOGADRO * 1e-25);
    // Isothermal compressibility in OpenMM internal units (mol*nm^3/kJ), consistent with kT (kJ/mol).
    double kappa_bar = owner.getIsothermalCompressibility();
    double betaT = kappa_bar / (AVOGADRO * 1e-25);
    double tauP = owner.getTauP();
    double dtBaro = integrator.getStepSize() * owner.getFrequency();

    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0] * box[1][1] * box[2][2];
    if (volume <= 0.0)
        return;

    double keAvg = 0.0;
    double virAvg = 0.0;
    vector<vector<Vec3> > savedPos(numCopies);

    for (int i = 0; i < numCopies; i++) {
        State state = integrator.getState(i, State::Positions | State::Forces | State::Energy);
        savedPos[i] = state.getPositions();
        keAvg += state.getKineticEnergy();
        const vector<Vec3>& pos = state.getPositions();
        const vector<Vec3>& forces = state.getForces();
        for (int j = 0; j < numParticles; j++) {
            virAvg += pos[j][0] * forces[j][0] + pos[j][1] * forces[j][1] + pos[j][2] * forces[j][2];
        }
    }
    keAvg /= numCopies;
    virAvg /= numCopies;

    double pint;
    if (owner.getScaleVelocities()) {
        pint = (2.0 * keAvg + virAvg) / (3.0 * volume);
    } else {
        int ndof = 3 * numParticles - context.getSystem().getNumConstraints();
        pint = (ndof * kT + virAvg) / (3.0 * volume);
    }

    double lambda = sqrt(volume);
    double lambdaD = 0.25 * kT * betaT / tauP;
    double lambdaF = -2.0 * lambda * (p0 - pint - kT / (2.0 * volume));

    // Box-Muller on seeded SFMT (reproducible with setRandomNumberSeed).
    double u1 = std::max(genrand_real2(random), std::numeric_limits<double>::min());
    double u2 = genrand_real2(random);
    const double twoPi = 2.0 * 3.14159265358979323846264338327950288;
    double xi = std::sqrt(-2.0 * std::log(u1)) * std::cos(twoPi * u2);
    double dLambda = (lambdaD / kT) * lambdaF * dtBaro + sqrt(2.0 * lambdaD * dtBaro) * xi;
    double depsilon = 2.0 * log(1.0 + dLambda / lambda);
    double scaling = exp(depsilon / 3.0);

    context.getOwner().setPeriodicBoxVectors(box[0] * scaling, box[1] * scaling, box[2] * scaling);

    vector<Vec3> centroid(numParticles, Vec3());
    for (int j = 0; j < numParticles; j++) {
        for (int i = 0; i < numCopies; i++)
            centroid[j] += savedPos[i][j];
        centroid[j] *= 1.0 / numCopies;
    }

    for (int j = 0; j < numParticles; j++) {
        Vec3 offset = centroid[j] * (scaling - 1.0);
        for (int i = 0; i < numCopies; i++)
            savedPos[i][j] += offset;
    }
    for (int i = 0; i < numCopies; i++)
        integrator.setPositions(i, savedPos[i]);

    if (owner.getScaleVelocities()) {
        double invScaling = 1.0 / scaling;
        for (int i = 0; i < numCopies; i++) {
            State state = integrator.getState(i, State::Velocities);
            vector<Vec3> vel = state.getVelocities();
            for (int j = 0; j < numParticles; j++)
                vel[j] *= invScaling;
            integrator.setVelocities(i, vel);
        }
    }

    double volumeNew = volume * scaling * scaling * scaling;
    engint += 0.5 * dLambda * lambdaF - (lambdaD / (4.0 * kT)) * dtBaro * lambdaF * lambdaF
            + 0.5 * kT * log(volumeNew / volume);
}

std::map<std::string, double> RPMDStochasticCellRescalingBarostatImpl::getDefaultParameters() {
    return {{RPMDStochasticCellRescalingBarostat::Pressure(), getOwner().getDefaultPressure()}};
}

std::vector<std::string> RPMDStochasticCellRescalingBarostatImpl::getKernelNames() {
    return {ApplyMonteCarloBarostatKernel::Name()};
}
