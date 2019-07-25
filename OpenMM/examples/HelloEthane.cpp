/* -----------------------------------------------------------------------------
 *               OpenMM(tm) HelloEthane example in C++ (June 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated simulation of a system with both bonded and nonbonded forces, 
 * using ethane (H3-C-C-H3) in a vacuum as an example. A multi-frame PDB file is 
 * written to stdout which can be read by VMD or other visualization tool to 
 * produce an animation of the resulting trajectory.
 *
 * Pay particular attention to the handling of units in this example. Incorrect
 * handling of units is a very common error; this example shows how you can
 * continue to work with Amber-style units of Angstroms and kCals while correctly
 * communicating with OpenMM in nanometers and kJoules.
 * -------------------------------------------------------------------------- */

#include <cstdio>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
//                                 MOCK MD CODE
// -----------------------------------------------------------------------------
// The code starting here and through main() below is meant to represent in 
// simplified form some pre-existing molecular dynamics code, which defines its 
// own data structures for force fields, the atoms in this simulation, and the 
// simulation parameters, and takes care of recording the trajectory. All this 
// has nothing to do with OpenMM; the OpenMM-dependent code comes later and is 
// clearly marked below.
// -----------------------------------------------------------------------------

//                   MODELING AND SIMULATION PARAMETERS
const bool   UseConstraints      = false;   // Should we constrain C-H bonds?
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 10;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 100;     // total simulation time (ps)
static const bool   WantEnergy   = true;

//                            FORCE FIELD DATA
// For this example we're using a tiny subset of the Amber99 force field.
// We want to keep the data in the original unit system to avoid conversion
// bugs; this requires conversion on the way in and out of OpenMM.

// Amber reduces nonbonded forces between 1-4 bonded atoms.
const double Coulomb14Scale      = 0.5;
const double LennardJones14Scale = 0.5;

struct AtomType {
    double mass, charge, vdwRadiusInAngstroms, vdwEnergyInKcal;
} atomType[] = {/*0 H*/ 1.008, 0.0605, 1.4870, 0.0157,
                /*1 C*/12.011, -.1815, 1.9080, 0.1094};
const int H = 0, C = 1;

struct BondType {
    double nominalLengthInAngstroms, stiffnessInKcalPerAngstrom2;
    bool   canConstrain;
} bondType[] = {/*0 CC*/1.526, 310., false,
                /*1 CH*/1.09 , 340., true};
const int CC = 0, CH = 1;

struct AngleType {
    double nominalAngleInDegrees, stiffnessInKcalPerRadian2;
} angleType[] = {/*0 HCC*/109.5, 50.,
                 /*1 HCH*/109.5, 35.};
const int HCC = 0, HCH = 1;

struct TorsionType {
    int    periodicity;
    double phaseInDegrees, amplitudeInKcal;
} torsionType[] = {/*0 HCCH*/3, 0., 0.150};
const int HCCH = 0;

//                                MOLECULE DATA
// Now describe an ethane molecule by listing its atoms, bonds, angles, and 
// torsions. We'll provide a default configuration which centers the molecule 
// at (0,0,0) with the C-C bond along the x axis.

// Use this as an "end of list" marker so that we do not have to count; let the
// computer do that!
const int EndOfList=-1;

struct MyAtomInfo
{   int type; const char* pdb; double initPosInAng[3]; double posInAng[3];} atoms[] = 
    {/*0*/C,       " C1 ",       { -.7605,  0,   0 },     {0,0,0},
     /*1*/C,       " C2 ",       {  .7605,  0,   0 },     {0,0,0},
     /*2*/H,       "1H1 ",       {-1.135, 1.03,  0 },     {0,0,0}, // bonded to C1
     /*3*/H,       "2H1 ",       {-1.135, -.51, .89},     {0,0,0},
     /*4*/H,       "3H1 ",       {-1.135, -.51,-.89},     {0,0,0},
     /*5*/H,       "1H2 ",       { 1.135, 1.03,  0 },     {0,0,0}, // bonded to C2
     /*6*/H,       "2H2 ",       { 1.135, -.51, .89},     {0,0,0},
     /*7*/H,       "3H2 ",       { 1.135, -.51,-.89},     {0,0,0},
     EndOfList};

static struct {int type; int atoms[2];} bonds[] = 
    {CC,0,1,
     CH,0,2,CH,0,3,CH,0,4,          // C1 methyl
     CH,1,5,CH,1,6,CH,1,7,          // C2 methyl     
     EndOfList};
static struct {int type; int atoms[3];} angles[] = 
    {HCC,2,0,1,HCC,3,0,1,HCC,4,0,1, // C1 methyl
     HCH,2,0,3,HCH,2,0,4,HCH,3,0,4,
     HCC,5,1,0,HCC,6,1,0,HCC,7,1,0, // C2 methyl
     HCH,5,1,6,HCH,5,1,7,HCH,6,1,7,             
     EndOfList};
static struct {int type; int atoms[4];} torsions[] = 
    {HCCH,2,0,1,5,HCCH,2,0,1,6,HCCH,2,0,1,7,
     HCCH,3,0,1,5,HCCH,3,0,1,6,HCCH,3,0,1,7,
     HCCH,4,0,1,5,HCCH,4,0,1,6,HCCH,4,0,1,7,    
     EndOfList};


//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal, 
                const MyAtomInfo atoms[]) 
{
    // Write out in PDB format -- printf is so much more compact than formatted cout.
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n", 
           timeInPs, energyInKcal);
    for (int n=0; atoms[n].type != EndOfList; ++n)
        printf("ATOM  %5d %4s ETH     1    %8.3f%8.3f%8.3f  1.00  0.00\n", 
            n+1, atoms[n].pdb, 
            atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    printf("ENDMDL\n");
}

// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------
// These four functions and an opaque structure are used to interface our main
// program with OpenMM without the main program having any direct interaction
// with the OpenMM API. This is a clean approach for interfacing with any MD
// code, although the details of the interface routines will differ. This is
// still just "locally written" code and is not required by OpenMM.
struct MyOpenMMData;
static MyOpenMMData* myInitializeOpenMM(const MyAtomInfo atoms[],
                                        double stepSizeInFs, 
                                        std::string& platformName);
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                                      double& time, double& energy, 
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);


// -----------------------------------------------------------------------------
//                           ETHANE MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        std::string   platformName;

        // Set up OpenMM data structures; returns OpenMM Platform name.
        MyOpenMMData* omm = myInitializeOpenMM(atoms, StepSizeInFs, platformName);

        // Run the simulation:
        //  (1) Write the first line of the PDB file and the initial configuration.
        //  (2) Run silently entirely within OpenMM between reporting intervals.
        //  (3) Write a PDB frame when the time comes.
        printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());

        const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
        for (int frame=1; ; ++frame) {
            double time, energy;
            myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
            myWritePDBFrame(frame, time, energy, atoms);

            if (time >= SimulationTimeInPs)
                break;

            myStepWithOpenMM(omm, NumSilentSteps);
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
// -----------------------------------------------------------------------------

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
#endif

#include "OpenMM.h"
using OpenMM::Vec3; // so we can just say "Vec3" below

// This is our opaque "handle" class containing all the OpenMM objects that
// must persist from call to call during a simulation. The main program gets 
// a pointer to one of these but sees it as essentially a void* since it 
// doesn't know the definition of this class.
struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0) {}
    ~MyOpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
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
                    double              stepSizeInFs, 
                    std::string&        platformName) 
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
       (OpenMM::Platform::getDefaultPluginsDirectory());

    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    OpenMM::NonbondedForce&         nonbond     = *new OpenMM::NonbondedForce();
    OpenMM::HarmonicBondForce&      bondStretch = *new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicAngleForce&     bondBend    = *new OpenMM::HarmonicAngleForce();
    OpenMM::PeriodicTorsionForce&   bondTorsion = *new OpenMM::PeriodicTorsionForce();
    system.addForce(&nonbond);
    system.addForce(&bondStretch);
    system.addForce(&bondBend);
    system.addForce(&bondTorsion);

    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atype.mass);
        nonbond.addParticle(atype.charge,
                            atype.vdwRadiusInAngstroms * OpenMM::NmPerAngstrom 
                                                       * OpenMM::SigmaPerVdwRadius,
                            atype.vdwEnergyInKcal      * OpenMM::KJPerKcal);
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }

    // Process the bonds:
    //  (1) If we're using constraints, tell System about constrainable bonds;
    //      otherwise, tell HarmonicBondForce the bond stretch parameters 
    //      (tricky units!).
    //  (2) Create a list of bonds for generating nonbond exclusions.
    std::vector< std::pair<int,int> > bondPairs;
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atom = bonds[i].atoms;
        const BondType& bond = bondType[bonds[i].type];

        if (UseConstraints && bond.canConstrain) {
            system.addConstraint(atom[0], atom[1],
                                 bond.nominalLengthInAngstroms * OpenMM::NmPerAngstrom);
        } else {
            // Note factor of 2 for stiffness below because Amber specifies the constant
            // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
            // it as used in the force term kx, with energy kx^2/2.
            bondStretch.addBond(atom[0], atom[1],
                                bond.nominalLengthInAngstroms    
                                    * OpenMM::NmPerAngstrom,
                                bond.stiffnessInKcalPerAngstrom2 
                                    * 2 * OpenMM::KJPerKcal 
                                    * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
        }

        bondPairs.push_back(std::make_pair(atom[0], atom[1]));
    }
    // Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    // Create the 1-2-3 bond angle harmonic terms.
    for (int i=0; angles[i].type != EndOfList; ++i) {
        const int*       atom  = angles[i].atoms;
        const AngleType& angle = angleType[angles[i].type];

        // See note under bond stretch above regarding the factor of 2 here.
        bondBend.addAngle(atom[0],atom[1],atom[2],
                          angle.nominalAngleInDegrees     * OpenMM::RadiansPerDegree,
                          angle.stiffnessInKcalPerRadian2 * 2 * OpenMM::KJPerKcal);
    }

    // Create the 1-2-3-4 bond torsion (dihedral) terms.
    for (int i=0; torsions[i].type != EndOfList; ++i) {
        const int*         atom = torsions[i].atoms;
        const TorsionType& torsion = torsionType[torsions[i].type];
        bondTorsion.addTorsion(atom[0],atom[1],atom[2],atom[3], 
            torsion.periodicity, 
            torsion.phaseInDegrees  * OpenMM::RadiansPerDegree,
            torsion.amplitudeInKcal * OpenMM::KJPerKcal);
    }

    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    omm->integrator = new OpenMM::VerletIntegrator(StepSizeInFs * OpenMM::PsPerFs);
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

