// -----------------------------------------------------------------------------
//           OpenMM(tm) HelloSodiumChloride example in C++ (June 2009)
// -----------------------------------------------------------------------------
// This is a complete, self-contained "hello world" example demonstrating 
// GPU-accelerated constant temperature simulation of a very simple system with
// just nonbonded forces, consisting of several sodium (Na+) and chloride (Cl-) 
// ions in implicit solvent. A multi-frame PDB file is written to stdout which 
// can be read by VMD or other visualization tool to produce an animation of the 
// resulting trajectory.
//
// Pay particular attention to the handling of units in this example. Incorrect
// handling of units is a very common error; this example shows how you can
// continue to work with Amber-style units like Angstroms, kCals, and van der
// Waals radii while correctly communicating with OpenMM in nm, kJ, and sigma.
// -----------------------------------------------------------------------------

#include <exception>
#include <string>
#include <cstdio>
using std::printf;


// -----------------------------------------------------------------------------
//                   MODELING AND SIMULATION PARAMETERS
// -----------------------------------------------------------------------------
static const double Temperature         = 300;     // Kelvins
static const double FrictionInPerPs     = 91.;     // collisions per picosecond
static const double SolventDielectric   = 80.;     // typical for water
static const double SoluteDielectric    = 2.;      // typical for protein

static const double StepSizeInFs        = 2;       // integration step size (fs)
static const double ReportIntervalInFs  = 50;      // how often to issue PDB frame (fs)
static const double SimulationTimeInPs  = 100;     // total simulation time (ps)

static const bool   WantEnergy          = true;


// -----------------------------------------------------------------------------
//                          ATOM AND FORCE FIELD DATA
// -----------------------------------------------------------------------------
// This is not part of OpenMM; just a struct we can use to collect atom 
// parameters for this example. Normally atom parameters would come from the 
// force field's parameterization file. We're going to use data in Angstrom and 
// Kilocalorie units and show how to safely convert to OpenMM's internal unit 
// system which uses nanometers and kilojoules.
static struct MyAtomInfo { 
    const char* pdb; 
    double      mass, charge, vdwRadiusInAng, vdwEnergyInKcal,
                gbsaRadiusInAng, gbsaScaleFactor;
    double      initPosInAng[3];
    double      posInAng[3]; // leave room for runtime state info
} atoms[] = {
// pdb   mass  charge  vdwRad vdwEnergy   gbsaRad gbsaScale  initPos
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     8, 0,  0},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,    -8, 0,  0},
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 9,  0},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0,-9,  0},
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 0,-10},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0, 0, 10},
{""} // end of list
};


// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------
// These four functions and an opaque structure are used to interface our main
// program with OpenMM without the main program having any direct interaction
// with the OpenMM API. This is a clean approach for interfacing with any MD
// code, although the details of the interface routines will differ.
struct MyOpenMMData;
static MyOpenMMData* myInitializeOpenMM(const MyAtomInfo atoms[],
                                        double temperature,
                                        double frictionInPs,
                                        double solventDielectric,
                                        double soluteDielectric,
                                        double stepSizeInFs, 
                                        std::string& platformName);
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                                      double& time, double& energy, 
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);


// -----------------------------------------------------------------------------
//                               PDB FILE WRITER
// -----------------------------------------------------------------------------
// Given state data, output a single frame (pdb "model") of the trajectory.
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal, 
                const MyAtomInfo atoms[]) 
{
    // Write out in PDB format -- printf is so much more compact than formatted cout.
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n", 
           timeInPs, energyInKcal);
    for (int n=0; *atoms[n].pdb; ++n)
        printf("ATOM  %5d %4s SLT     1    %8.3f%8.3f%8.3f  1.00  0.00\n", 
            n+1, atoms[n].pdb, 
            atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    printf("ENDMDL\n");
}


// -----------------------------------------------------------------------------
//                                MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    const int NumReports     = (int)(SimulationTimeInPs*1000 / ReportIntervalInFs + 0.5);
    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);

    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        double        time, energy;
        std::string   platformName;

        // Set up OpenMM data structures; returns OpenMM Platform name.
        MyOpenMMData* omm = myInitializeOpenMM(atoms, Temperature, FrictionInPerPs,
                                               SolventDielectric, SoluteDielectric,
                                               StepSizeInFs, platformName);

        // Run the simulation:
        //  (1) Write the first line of the PDB file and the initial configuration.
        //  (2) Run silently entirely within OpenMM between reporting intervals.
        //  (3) Write a PDB frame when the time comes.
        printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
        myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
        myWritePDBFrame(1, time, energy, atoms);

        for (int frame=2; frame <= NumReports; ++frame) {
            myStepWithOpenMM(omm, NumSilentSteps);
            myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
            myWritePDBFrame(frame, time, energy, atoms);
        } 
 
        // Clean up OpenMM data structures.
        myTerminateOpenMM(omm);

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}


// -----------------------------------------------------------------------------
//                           OpenMM-USING CODE
// -----------------------------------------------------------------------------
// The OpenMM API is visible only at this point and below. Normally this would
// be in a separate compilation module; we're including it here for simplicity.

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
#endif

#include "OpenMM.h"
using OpenMM::Vec3; // so we can just say "Vec3" below

struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0) {}
    ~MyOpenMMData() {delete system; delete context; delete integrator;}
    OpenMM::System*         system;
    OpenMM::Context*        context;
    OpenMM::Integrator*     integrator;
};


// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
// (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     in a manner which is opaque to the caller.
// (3) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (4) Create an Integrator and a Context associating the Integrator with
//     the System.
// (5) Select the OpenMM platform to be used.
// (6) Return the MyOpenMMData struct and the name of the Platform in use.
//
// Note that this function must understand the calling MD code's molecule and
// force field data structures so will need to be customized for each MD code.
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPs,
                    double              solventDielectric,
                    double              soluteDielectric,
                    double              stepSizeInFs, 
                    std::string&        platformName) 
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
       (OpenMM::Platform::getDefaultPluginsDirectory());

    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the OpenMM
    // System takes ownership of the force objects; don't delete them yourself.
    omm->system = new OpenMM::System();
    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
    OpenMM::GBSAOBCForce*   gbsa    = new OpenMM::GBSAOBCForce();
    omm->system->addForce(nonbond);
    omm->system->addForce(gbsa);

    // Specify dielectrics for GBSA implicit solvation.
    gbsa->setSolventDielectric(solventDielectric);
    gbsa->setSoluteDielectric(soluteDielectric);

    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) GBSA needs charge, radius, and scale factor.
    //  (4) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; *atoms[n].pdb; ++n) {
        const MyAtomInfo& atom = atoms[n];

        omm->system->addParticle(atom.mass);

        nonbond->addParticle(atom.charge,
                             atom.vdwRadiusInAng  * OpenMM::NmPerAngstrom 
                                                  * OpenMM::SigmaPerVdwRadius,
                             atom.vdwEnergyInKcal * OpenMM::KJPerKcal);

        gbsa->addParticle(atom.charge, 
                          atom.gbsaRadiusInAng * OpenMM::NmPerAngstrom,
                          atom.gbsaScaleFactor);

        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atom.initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atom.initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atom.initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }

    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero but could
    // have been set here.
    omm->integrator = new OpenMM::LangevinIntegrator(temperature, frictionInPs, 
                                                     stepSizeInFs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initialPosInNm);

    platformName = omm->context->getPlatform().getName();
    return omm;
}


// -----------------------------------------------------------------------------
//                     COPY STATE BACK TO CPU FROM OPENMM
// -----------------------------------------------------------------------------
static void
myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy, 
                 double& timeInPs, double& energyInKcal,
                 MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).

    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;

    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKcal = 0;
    if (wantEnergy) 
        energyInKcal = (state.getPotentialEnergy() + state.getKineticEnergy())
                        * OpenMM::KcalPerKJ;
}


// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM 
// -----------------------------------------------------------------------------
static void 
myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    omm->integrator->step(numSteps);
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
static void 
myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}
