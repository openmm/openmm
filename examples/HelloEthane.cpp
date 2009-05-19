/* -----------------------------------------------------------------------------
 *                    OpenMM(tm) HelloEthane example (May 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated simulation of a system with both bonded and nonbonded forces, 
 * using ethane (H3-C-C-H3) as an example. A multi-frame PDB file is written
 * to stdout which can be read by VMD or other visualization tool to produce 
 * an animation of the resulting trajectory.
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

#include <iostream>
#include <string>
#include <vector>
#include <utility>

using namespace OpenMM;

// These are missing from the current version of OpenMM so we're adding them
// temporarily here.
Vec3 operator*(const Vec3& v, double r) {return Vec3(v[0]*r, v[1]*r, v[2]*r);}
Vec3 operator*(double r, const Vec3& v) {return Vec3(r*v[0], r*v[1], r*v[2]);}
static const double SigmaPerVdwRadius = 2*std::pow(2., -1./6.);

// -----------------------------------------------------------------------------
//              FORCE FIELD DATA
// -----------------------------------------------------------------------------
// These data structures are not part of OpenMM; they are a model of the kinds
// of data structures an MD code uses to hold a set of force field parameters.
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

// -----------------------------------------------------------------------------
//              MOLECULE DATA
// -----------------------------------------------------------------------------
// Now describe an ethane molecule by listing its atoms, bonds, angles, and 
// torsions. We'll provide a default configuration which centers the molecule 
// at (0,0,0) with the C-C bond along the x axis.

// Use this as an "end of list" marker so that we do not have to count; let the
// computer do that!
const int EndOfList=-1;

struct AtomInfo {
    int type; char* pdbSymbol; Vec3 initPosInAngstroms;
} atoms[] = {/*0*/C, "C1", Vec3( -.765,   0,   0 ),
             /*1*/C, "C2", Vec3(  .765,   0,   0 ),
             /*2*/H, "H1", Vec3(-1.135, 1.03,  0 ), // bonded to C1
             /*3*/H, "H2", Vec3(-1.135, -.51, .89),
             /*4*/H, "H3", Vec3(-1.135, -.51,-.89),
             /*5*/H, "H4", Vec3( 1.135, 1.03,  0 ), // bonded to C2
             /*6*/H, "H5", Vec3( 1.135, -.51, .89),
             /*7*/H, "H6", Vec3( 1.135, -.51,-.89),
             EndOfList};

struct {int type; int a[2];} bonds[]    = {CC,0,1,CH,0,2,CH,0,3,CH,0,4,
                                                  CH,1,5,CH,1,6,CH,1,7, 
                                           EndOfList};
struct {int type; int a[3];} angles[]   = {HCC,2,0,1,HCC,3,0,1,HCC,4,0,1,
                                           HCC,5,1,0,HCC,6,1,0,HCC,7,1,0,
                                           HCH,2,0,3,HCH,2,0,4,HCH,3,0,4,
                                           HCH,5,1,6,HCH,5,1,7,HCH,6,1,7, 
                                           EndOfList};
struct {int type; int a[4];} torsions[] = {HCCH,2,0,1,5,HCCH,2,0,1,6,HCCH,2,0,1,7,
                                           HCCH,3,0,1,5,HCCH,3,0,1,6,HCCH,3,0,1,7,
                                           HCCH,4,0,1,5,HCCH,4,0,1,6,HCCH,4,0,1,7, 
                                           EndOfList};

// -----------------------------------------------------------------------------
//              MODELING AND SIMULATION PARAMETERS
// -----------------------------------------------------------------------------
const bool   UseConstraints      = false;   // Should we constrain C-H bonds?
const double Temperature         = 300;     // bath temperature in Kelvins
const double FrictionInPs        = 1./91.;  // picoseconds between collisions

const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 10;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 100;     // total simulation time (ps)

// PDB file writer; see below.
static void writePDB(const OpenMMContext&);

// -----------------------------------------------------------------------------
//                                MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
  // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
  // usage and runtime errors are caught and reported.
  try {

    // -------------------------------------------------------------------------
    // Load all available OpenMM plugins from their default location.
    // -------------------------------------------------------------------------
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());

    // -------------------------------------------------------------------------
    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces.
    // -------------------------------------------------------------------------
    System system;
    NonbondedForce&         nonbond     = *new NonbondedForce();
    HarmonicBondForce&      bondStretch = *new HarmonicBondForce();
    HarmonicAngleForce&     bondBend    = *new HarmonicAngleForce();
    PeriodicTorsionForce&   bondTorsion = *new PeriodicTorsionForce();
    system.addForce(&nonbond);
    system.addForce(&bondStretch);
    system.addForce(&bondBend);
    system.addForce(&bondTorsion);

    // -------------------------------------------------------------------------
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    // -------------------------------------------------------------------------
    std::vector<Vec3> initialPositions;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atype.mass);
        nonbond.addParticle(atype.charge,
                            atype.vdwRadiusInAngstroms * NmPerAngstrom * SigmaPerVdwRadius,
                            atype.vdwEnergyInKcal      * KJPerKcal);
        initialPositions.push_back(atoms[n].initPosInAngstroms * NmPerAngstrom);
    }

    // -------------------------------------------------------------------------
    // Process the bonds:
    //  (1) HarmonicBondForce needs bond stretch parameters (in MD units!).
    //  (2) If we're using constraints, tell System about constrainable bonds.
    //  (3) Create a list of bonds for generating nonbond exclusions.
    // -------------------------------------------------------------------------
    std::vector< std::pair<int,int> > bondPairs;
    for (int i=0; bonds[i].type != EndOfList; ++i) {
        const int*      atoms = bonds[i].a;
        const BondType& bond  = bondType[bonds[i].type];

        // Note factor of 2 for stiffness below because Amber specifies the constant
        // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
        // it as used in the force term kx, with energy kx^2/2.
        bondStretch.addBond(atoms[0], atoms[1],
                            bond.nominalLengthInAngstroms    * NmPerAngstrom,
                            bond.stiffnessInKcalPerAngstrom2 * 2 * KJPerKcal * AngstromsPerNm * AngstromsPerNm);
        if (UseConstraints && bond.canConstrain)
            system.addConstraint(atoms[0], atoms[1],
                                 bond.nominalLengthInAngstroms * NmPerAngstrom);
        bondPairs.push_back(std::make_pair(atoms[0], atoms[1]));
    }
    // Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms.
    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    // -------------------------------------------------------------------------
    // Create the 1-2-3 bond angle harmonic terms.
    // -------------------------------------------------------------------------
    for (int i=0; angles[i].type != EndOfList; ++i) {
        const int*       atoms = angles[i].a;
        const AngleType& angle = angleType[angles[i].type];

        // See note under bond stretch above regarding the factor of 2 here.
        bondBend.addAngle(atoms[0],atoms[1],atoms[2],
                          angle.nominalAngleInDegrees     * RadiansPerDegree,
                          angle.stiffnessInKcalPerRadian2 * 2 * KJPerKcal);
    }

    // -------------------------------------------------------------------------
    // Create the 1-2-3-4 bond torsion (dihedral) terms.
    // -------------------------------------------------------------------------
    for (int i=0; torsions[i].type != EndOfList; ++i) {
        const int*         atoms = torsions[i].a;
        const TorsionType& torsion = torsionType[torsions[i].type];
        bondTorsion.addTorsion(atoms[0],atoms[1],atoms[2],atoms[3], 
            torsion.periodicity, 
            torsion.phaseInDegrees  * RadiansPerDegree,
            torsion.amplitudeInKcal * KJPerKcal);
    }

    // -------------------------------------------------------------------------
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    // -------------------------------------------------------------------------
    //LangevinIntegrator integrator(Temperature, FrictionInPs, StepSizeInFs * PsPerFs);
    VerletIntegrator integrator(StepSizeInFs * PsPerFs);
    OpenMMContext    context(system, integrator);
    context.setPositions(initialPositions);

    // -------------------------------------------------------------------------
    // Run the simulation:
    //  (1) Write the first line of the PDB file and the initial configuration.
    //  (2) Run silently entirely within OpenMM between reporting intervals.
    //  (3) Write a PDB frame when the time comes.
    // -------------------------------------------------------------------------
    printf("REMARK  Using OpenMM platform %s\n", context.getPlatform().getName().c_str() );
    writePDB(context);

    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    do {
        integrator.step(NumSilentSteps);
        writePDB(context);
    } while (context.getTime() < SimulationTimeInPs);

    // -------------------------------------------------------------------------
    // Normal return from main.
    // -------------------------------------------------------------------------
    return 0;

  // Catch and report usage and runtime errors detected by OpenMM and fail.
  } catch(const std::exception& e) {
    printf("EXCEPTION: %s\n", e.what());
    return 1;
  }
}


// -----------------------------------------------------------------------------
//                               PDB FILE WRITER
// -----------------------------------------------------------------------------
static void
writePDB(const OpenMMContext& context) {
    // Caution: at the moment asking for energy requires use of slow Reference 
    // platform calculation.
    const State                 state       = context.getState(State::Positions | State::Velocities | State::Energy);
    const double                energy      = state.getPotentialEnergy() + state.getKineticEnergy();
    const std::vector<Vec3>&    positions   = state.getPositions();
    static int modelFrameNumber = 0; // numbering for MODEL records in pdb output

    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    printf("REMARK 250 time=%.3f picoseconds; Energy = %.3f kilojoules/mole\n", state.getTime(), energy);
    for (unsigned i=0; i < positions.size(); ++i) {
        const Vec3 pos = positions[i] * AngstromsPerNm;
        printf("ATOM  %5d %2s   ETH     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", 
            i+1, atoms[i].pdbSymbol, pos[0], pos[1], pos[2], atoms[i].pdbSymbol);
    }
    printf("ENDMDL\n");
}

