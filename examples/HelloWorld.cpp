// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
    #pragma warning(disable:4251)   // no dll interface for some classes
#endif

#include "OpenMM.h"

#include <iostream>
#include <string>
#include <limits>

using namespace OpenMM;

Vec3 operator*(const Vec3& v, double r) {
    return Vec3(v[0]*r, v[1]*r, v[2]*r);
}

Vec3 operator*(double r, const Vec3& v) {
    return Vec3(r*v[0], r*v[1], r*v[2]);
}

// This is not part of OpenMM; just a struct we can use to collect
// atom parameters for this example. Normally atom parameters would
// come from the force field's parameterization file.
// We're going to use data in Angstrom and Kilocalorie units and
// show how to safely convert to OpenMM's internal unit system
// which uses nanometers and kilojoules.
struct AtomInfo {
    char*  symbol;
    double mass, charge, vdwRadiusAng, vdwEnergyKcal;
    Vec3   startPosAng;
};

static AtomInfo atoms[] = {
    {"Na+", 22.99,  1, 1.8680, 0.00277, Vec3(-3,0,0)},
    {"Cl-", 35.45, -1, 2.4700, 0.1000,  Vec3( 3,0,0)},
    {""} // end of list
};

static const double Temperature         = 100;   // Kelvins
static const double FrictionInPerPs     = 91.;   // collisions per ps
static const double StepSizeFs          = 2;     // femtoseconds
static const double ReportIntervalFs    = 1000;
static const double SimulationTimePs    = 1000;  // total simulation time (ps)

static const double SigmaPerVdwRadius = 2*std::pow(2., -1./6.);

static double showState(const OpenMMContext&);

int main() {
try {
    // Load all available OpenMM plugins from their default location.
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

    // Create a System and a NonbondedForce object within the System. 
    System system;
    NonbondedForce* nonbond = new NonbondedForce();
    system.addForce(nonbond);

    int numAtoms = 0;
    for (; *atoms[numAtoms].symbol; ++numAtoms) {
        const AtomInfo& atom = atoms[numAtoms];
        system.addParticle(atom.mass);
        nonbond->addParticle(atom.charge,
                             atom.vdwRadiusAng  * NmPerAngstrom * SigmaPerVdwRadius,
                             atom.vdwEnergyKcal * KJPerKcal);
    }

    // Create an integrator object for advancing time.
    //LangevinIntegrator integrator(Temperature, Friction, StepSizeFs * PsPerFs);
    VerletIntegrator integrator(StepSizeFs * PsPerFs);

    // Create an OpenMM Context for execution; let it choose best platform.
    OpenMMContext context(system, integrator);

    std::cout << "Will use OpenMM platform " << context.getPlatform().getName() << std::endl;
    std::cout << "eps(float)=" << std::numeric_limits<float>::epsilon() << std::endl;

    // Fill in a vector of starting positions, one per atom.
    std::vector<Vec3> positions(numAtoms);
    for (int i=0; i < numAtoms; ++i)
        positions[i] = atoms[i].startPosAng * NmPerAngstrom;

    // Set the starting positions in the OpenMM context. Velocities will be zero.
    context.setPositions(positions);

    // Output the initial state.
    double time = showState(context);

    const int NumSilentSteps = (int)(ReportIntervalFs / StepSizeFs + 0.5);
    do {
        integrator.step(NumSilentSteps);
        time = showState(context);
    } while (time < SimulationTimePs);

    } catch(const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

static double
showState(const OpenMMContext& context) {
    // Caution: at the moment asking for energy requires use of slow reference calculation.
    const State                 state       = context.getState(  State::Positions | State::Velocities 
                                                               | State::Forces    | State::Energy);
    const double                energy      = state.getPotentialEnergy() + state.getKineticEnergy();
    const std::vector<Vec3>&    positions   = state.getPositions();
    const std::vector<Vec3>&    forces      = state.getForces();

    std::cout << "t=" << state.getTime() << " E=" << energy << std::endl;
    const double dt = StepSizeFs*PsPerFs;
    for (unsigned i=0; i < positions.size(); ++i) {
        const Vec3 dq = dt*dt/(2*atoms[i].mass) * forces[i];
        std::cout << i << " " << positions[i] * AngstromsPerNm << std::endl;
        std::cout << "f=" << forces[i] << " dq=" << dq << std::endl;
    }

    return state.getTime();
}

