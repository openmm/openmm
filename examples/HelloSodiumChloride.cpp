/* -----------------------------------------------------------------------------
 *                OpenMM(tm) HelloSodiumChloride example (May 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated constant energy simulation of a very simple system with just 
 * nonbonded forces, consisting of several sodium (Na+) and chloride (Cl-) ions. 
 * A multi-frame PDB file is written to stdout which can be read by VMD or other 
 * visualization tool to produce an animation of the resulting trajectory.
 *
 * Pay particular attention to the handling of units in this example. Incorrect
 * handling of units is a very common error; this example shows how you can
 * continue to work with Amber-style units of Angstroms and kCals while correctly
 * communicating with OpenMM in nanometers and kJoules.
 * -------------------------------------------------------------------------- */

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
    #pragma warning(disable:4251)   // no dll interface for some classes
#endif

#include "OpenMM.h"

#include <cstdio>
#include <string>

using OpenMM::Vec3; // so we can just say "Vec3" below

// -----------------------------------------------------------------------------
//                   MODELING AND SIMULATION PARAMETERS
// -----------------------------------------------------------------------------
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 10;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 100;     // total simulation time (ps)

static void simulateNaCl();
static void writePDB(const OpenMM::OpenMMContext&); // PDB file writer; see below.

// -----------------------------------------------------------------------------
//                                MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        // Load all available OpenMM plugins from their default location.
        OpenMM::Platform::loadPluginsFromDirectory
           (OpenMM::Platform::getDefaultPluginsDirectory());

        simulateNaCl();

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}

// -----------------------------------------------------------------------------
//                          ATOM AND FORCE FIELD DATA
// -----------------------------------------------------------------------------
// This is not part of OpenMM; just a struct we can use to collect
// atom parameters for this example. Normally atom parameters would
// come from the force field's parameterization file.
// We're going to use data in Angstrom and Kilocalorie units and
// show how to safely convert to OpenMM's internal unit system
// which uses nanometers and kilojoules.
struct AtomInfo { 
    const char* pdb; 
    double      mass, charge, vdwRadiusInAng, vdwEnergyInKcal;
    Vec3        initPosInAngstroms;
} atoms[] = {
    // pdb   mass   charge   vdwRadius  vdwEnergy   initPos
    {" NA ", 22.99,   1,     1.8680,    0.00277,    Vec3(8,0,0)},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     Vec3(-8,0,0)},
    {" NA ", 22.99,   1,     1.8680,    0.00277,    Vec3(0,9,0)},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     Vec3(0,-9,0)},
    {" NA ", 22.99,   1,     1.8680,    0.00277,    Vec3(0,0,-10)},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     Vec3( 0,0,10)},
    {""} // end of list
};

// Add missing scalar product operators for OpenMM::Vec3.
Vec3 operator*(const Vec3& v, double r) {return Vec3(v[0]*r, v[1]*r, v[2]*r);}
Vec3 operator*(double r, const Vec3& v) {return Vec3(r*v[0], r*v[1], r*v[2]);}
// This is the conversion factor that takes you from a van der Waals radius
// (defined as 1/2 the minimum energy separation) to the related Lennard Jones 
// "sigma" parameter (defined as the zero crossing separation).
static const double SigmaPerVdwRadius = 2*std::pow(2., -1./6.);

// -----------------------------------------------------------------------------
//                              NaCl SIMULATION
// -----------------------------------------------------------------------------
static void simulateNaCl() {
    // -------------------------------------------------------------------------
    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: OpenMM owns
    // the objects and will take care of deleting them; don't do it yourself!
    // -------------------------------------------------------------------------
    OpenMM::System system;
    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
    system.addForce(nonbond);

    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond->setCutoffDistance(2);
    nonbond->setPeriodicBoxVectors(Vec3(5,0,0), Vec3(0,5,0), Vec3(0,0,5));

    // -------------------------------------------------------------------------
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    // -------------------------------------------------------------------------
    std::vector<Vec3> initialPositions;
    for (int n=0; *atoms[n].pdb; ++n) {
        const AtomInfo& atom = atoms[n];
        system.addParticle(atom.mass);
        nonbond->addParticle(atom.charge,
                             atom.vdwRadiusInAng  * OpenMM::NmPerAngstrom 
                                                  * SigmaPerVdwRadius,
                             atom.vdwEnergyInKcal * OpenMM::KJPerKcal);
        initialPositions.push_back(atoms[n].initPosInAngstroms 
                                                  * OpenMM::NmPerAngstrom);
    }

    // -------------------------------------------------------------------------
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    // -------------------------------------------------------------------------
    OpenMM::VerletIntegrator integrator(StepSizeInFs * OpenMM::PsPerFs);
    OpenMM::OpenMMContext    context(system, integrator);
    context.setPositions(initialPositions);

    // -------------------------------------------------------------------------
    // Run the simulation:
    //  (1) Write the first line of the PDB file and the initial configuration.
    //  (2) Run silently entirely within OpenMM between reporting intervals.
    //  (3) Write a PDB frame when the time comes.
    // -------------------------------------------------------------------------
    printf( "REMARK  Using OpenMM platform %s\n", context.getPlatform().getName().c_str() );
    writePDB(context);

    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    do {
        integrator.step(NumSilentSteps);
        writePDB(context);
    } while (context.getTime() < SimulationTimeInPs);
}

// -----------------------------------------------------------------------------
//                               PDB FILE WRITER
// -----------------------------------------------------------------------------
static void
writePDB(const OpenMM::OpenMMContext& context) {
    // Caution: at the moment asking for energy requires use of slow reference calculation.
    const OpenMM::State state = context.getState(  OpenMM::State::Positions 
                                                 | OpenMM::State::Velocities 
                                                 | OpenMM::State::Energy);
    const double energy = state.getPotentialEnergy() + state.getKineticEnergy();
    const std::vector<Vec3>& positions = state.getPositions();
    static int modelFrameNumber = 0; // numbering for MODEL records in pdb output

    // write out in PDB format -- printf is so much more compact than formatted cout
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    printf("REMARK 250 time=%.3f picoseconds; Energy = %.3f kilojoules/mole\n", 
           state.getTime(), energy);
    for (unsigned i=0; i < positions.size(); ++i) {
        const Vec3 pos = positions[i] * OpenMM::AngstromsPerNm;
        printf("ATOM  %5d %4s SLT     1    %8.3f%8.3f%8.3f  1.00  0.00            \n", 
            i+1, atoms[i].pdb, pos[0], pos[1], pos[2]);
    }
    printf("ENDMDL\n");
}

