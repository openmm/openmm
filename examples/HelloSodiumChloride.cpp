#include "OpenMM.h"

#include <iostream>
#include <iomanip>
#include <string>

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
    {"NA", 22.99,  1, 1.8680, 0.00277, Vec3(8,0,0)},
    {"CL", 35.45, -1, 2.4700, 0.1000,  Vec3(-8,0,0)},
    {"NA", 22.99,  1, 1.8680, 0.00277, Vec3(0,9,0)},
    {"CL", 35.45, -1, 2.4700, 0.1000,  Vec3(0,-9,0)},
    {"NA", 22.99,  1, 1.8680, 0.00277, Vec3(0,0,-10)},
    {"CL", 35.45, -1, 2.4700, 0.1000,  Vec3( 0,0,10)},
    {""} // end of list
};

static const double Temperature         = 300;    // Kelvins
static const double Friction            = 1./91.; // picoseconds between collisions
static const double StepSizeFs          = 2;      // femtoseconds
static const double ReportIntervalFs    = 10;
static const double SimulationTimePs    = 100;    // total simulation time (ps)

static const double SigmaPerVdwRadius = 2*std::pow(2., -1./6.);

static void writePDB(const OpenMMContext&);

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
    LangevinIntegrator integrator(Temperature, Friction, StepSizeFs * PsPerFs);
    //VerletIntegrator integrator(StepSizeFs * PsPerFs);

    // Create an OpenMM Context for execution; let it choose best platform.
    OpenMMContext context(system, integrator);

    const std::string platformName = context.getPlatform().getName();
    printf( "REMARK  Using OpenMM platform %s\n", context.getPlatform().getName().c_str() );

    // Fill in a vector of starting positions, one per atom.
    std::vector<Vec3> positions(numAtoms);
    for (int i=0; i < numAtoms; ++i)
        positions[i] = atoms[i].startPosAng * NmPerAngstrom;

    // Set the starting positions in the OpenMM context. Velocities will be zero.
    context.setPositions(positions);

    // Output the initial state.
    writePDB(context);

    const int NumSilentSteps = (int)(ReportIntervalFs / StepSizeFs + 0.5);
    do {
        integrator.step(NumSilentSteps);
        writePDB(context);
    } while (context.getTime() < SimulationTimePs);

    } catch(const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

static void
writePDB(const OpenMMContext& context) {
    // Caution: at the moment asking for energy requires use of slow reference calculation.
    const State                 state       = context.getState(State::Positions | State::Velocities | State::Energy);
    const double                energy      = state.getPotentialEnergy() + state.getKineticEnergy();
    const std::vector<Vec3>&    positions   = state.getPositions();
    static int modelFrameNumber = 0; // numbering for MODEL records in pdb output

    // write out in PDB format -- printf is so much more compact than formatted cout
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    printf("REMARK 250 time=%.3f picoseconds; Energy = %.3f kilojoules/mole\n", state.getTime(), energy);
    for (unsigned i=0; i < positions.size(); ++i) {
        const Vec3 pos = positions[i] * AngstromsPerNm;
        printf("ATOM    %3d %2s   SLT     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", 
            i+1, atoms[i].symbol, pos[0], pos[1], pos[2], atoms[i].symbol);
    }
    printf("ENDMDL\n");
}

