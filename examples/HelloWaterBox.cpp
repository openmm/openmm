/* -----------------------------------------------------------------------------
 *                    OpenMM(tm) HelloWaterBox example (May 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated simulation of a system with both bonded and nonbonded forces, 
 * using water (H-O-H) in a periodic box as an example. A multi-frame PDB file is 
 * written to stdout which can be read by VMD or other visualization tool to 
 * produce an animation of the resulting trajectory.
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
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>

using OpenMM::Vec3; // so we can just say "Vec3" below

// -----------------------------------------------------------------------------
//                   MODELING AND SIMULATION PARAMETERS
// -----------------------------------------------------------------------------
const bool   UseConstraints      = true;   // Should we constrain O-H bonds?
const double StepSizeInFs        = 2;      // integration step size (fs)
const double ReportIntervalInFs  = 100;    // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 10;    // total simulation time (ps)

static void simulateWaterBox();
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
        simulateWaterBox();

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}


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

// This is the conversion factor that takes you from a van der Waals radius
// (defined as 1/2 the minimum energy separation) to the related Lennard Jones 
// "sigma" parameter (defined as the zero crossing separation).
static const double SigmaPerVdwRadius = 2*std::pow(2., -1./6.);

// -----------------------------------------------------------------------------
//                            WATER SIMULATION
// -----------------------------------------------------------------------------
static void simulateWaterBox() {
    // -------------------------------------------------------------------------
    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: OpenMM owns
    // the objects and will take care of deleting them; don't do it yourself!
    // -------------------------------------------------------------------------
    OpenMM::System system;

    OpenMM::NonbondedForce&         nonbond     = *new OpenMM::NonbondedForce();
    system.addForce(&nonbond);

    OpenMM::HarmonicBondForce&      bondStretch = *new OpenMM::HarmonicBondForce();
    system.addForce(&bondStretch);

    OpenMM::HarmonicAngleForce&     bondBend    = *new OpenMM::HarmonicAngleForce();
    system.addForce(&bondBend);

    OpenMM::AndersenThermostat&     thermostat  = *new OpenMM::AndersenThermostat(
            300, // temperature in kelvins
            91.0); // collision frequency in 1/picoseconds

    // -------------------------------------------------------------------------
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    // -------------------------------------------------------------------------
    
    // Volume of one water is 30 Angstroms cubed;
    // Thus length in one dimension is cube-root of 30,
    // or 3.107 Angstroms or 0.3107 nanometers
    double waterSize = 0.3107; // edge of cube containing one water, in nanometers

    // Place 1000 water molecules one at a time in a 10x10x10 rectilinear grid
    const int watersPerEdge = 10;
    double boxEdgeLength = waterSize * watersPerEdge;

    // Create periodic box
    nonbond.setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond.setCutoffDistance(2);
    nonbond.setPeriodicBoxVectors(Vec3(boxEdgeLength,0,0), Vec3(0,boxEdgeLength,0), Vec3(0,0,boxEdgeLength));

    // Create data structures to hold lists of initial positions and bonds
    std::vector<Vec3> initialPositions;
    std::vector< std::pair<int,int> > bondPairs;
    
    // Add water molecules one at a time in the 10x10x10 cubic lattice
    for (int latticeX = 0; latticeX < watersPerEdge; ++latticeX)
        for (int latticeY = 0; latticeY < watersPerEdge; ++latticeY)
            for (int latticeZ = 0; latticeZ < watersPerEdge; ++latticeZ)
            {
                // Add parameters for one water molecule
                
                // Add atom masses to system
                int  oIndex = system.addParticle(15.9994); // O
                int h1Index = system.addParticle(1.00794); // H1
                int h2Index = system.addParticle(1.00794); // H2
                
                // Add atom charge, sigma, and stiffness to nonbonded force
                nonbond.addParticle( // Oxygen
                        -0.834, // charge
                        1.7683 * OpenMM::NmPerAngstrom * SigmaPerVdwRadius,
                        0.1520 * OpenMM::KJPerKcal);
                nonbond.addParticle( // Hydrogen1
                        0.417, // charge
                        0.0001 * OpenMM::NmPerAngstrom * SigmaPerVdwRadius,
                        0.0000 * OpenMM::KJPerKcal);
                nonbond.addParticle( // Hydrogen2
                        0.417, // charge
                        0.0001 * OpenMM::NmPerAngstrom * SigmaPerVdwRadius,
                        0.0000 * OpenMM::KJPerKcal);
                
                // Add stretch parameters for two covalent bonds
                // Note factor of 2 for stiffness below because Amber specifies the constant
                // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
                // it as used in the force term kx, with energy kx^2/2.
                bondStretch.addBond(oIndex, h1Index,
                        0.9572 * OpenMM::NmPerAngstrom,
                        553.0 * 2 * OpenMM::KJPerKcal 
                            * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                bondStretch.addBond(oIndex, h2Index,
                        0.9572 * OpenMM::NmPerAngstrom,
                        553.0 * 2 * OpenMM::KJPerKcal 
                            * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);

                // constrain O-H bond lengths
                if (UseConstraints) {
                    system.addConstraint(oIndex, h1Index,
                                     0.9572 * OpenMM::NmPerAngstrom);
                    system.addConstraint(oIndex, h2Index,
                                     0.9572 * OpenMM::NmPerAngstrom);
                }

                // Store bonds for exclusion list
                bondPairs.push_back(std::make_pair(oIndex, h1Index));
                bondPairs.push_back(std::make_pair(oIndex, h2Index));
                            
                // Add bond bend parameters for one angle
                // See note under bond stretch above regarding the factor of 2 here.
                bondBend.addAngle(h1Index, oIndex, h2Index,
                        104.52 * OpenMM::RadiansPerDegree,
                        100.00 * 2 * OpenMM::KJPerKcal);
                       
                // Location of this molecule in the lattice
                Vec3 latticeVec(waterSize * latticeX, 
                                waterSize * latticeY, 
                                waterSize * latticeZ);

                // flip half the waters to prevent giant dipole
                int flip = (rand() % 100) > 49 ? 1 : -1;

                // place this water
                initialPositions.push_back(Vec3(0,0,0) + latticeVec); // O
                initialPositions.push_back(Vec3(0.09572*flip,0,0) + latticeVec); // H1
                initialPositions.push_back(Vec3(-0.02397*flip,0.09267*flip,0) + latticeVec); // H2
            }
    
    // Populate nonbonded exclusions
    nonbond.createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);    
            
    // -------------------------------------------------------------------------
    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    // -------------------------------------------------------------------------

    OpenMM::VerletIntegrator integrator(
            StepSizeInFs * OpenMM::PsPerFs);
    OpenMM::OpenMMContext    context(system, integrator);
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
}



// -----------------------------------------------------------------------------
//                               PDB FILE WRITER
// -----------------------------------------------------------------------------
static void
writePDB(const OpenMM::OpenMMContext& context) {
    // We don't calculate energy in this example because it would take most of the time
    const OpenMM::State state = context.getState(OpenMM::State::Positions);
    const std::vector<Vec3>& positions = state.getPositions();
    static int modelFrameNumber = 0; // numbering for MODEL records in pdb output

    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    printf("REMARK 250 time=%.3f picoseconds\n", state.getTime());
    int residueNumber = 1;
    char* atomName = " O  ";
    for (unsigned int atomIndex=0; atomIndex < positions.size(); ++atomIndex) 
    {
        // cycle through same three atom names
        if (0 == atomIndex%3) {
            atomName = " O  ";
            // increment residue number once every three atoms
            ++residueNumber;
        }
        else if (1 == atomIndex%3) atomName = " H1 ";
        else atomName = " H2 ";

        printf("HETATM%5d %4s HOH  %4d    ", // start of pdb ATOM line
            atomIndex, atomName, residueNumber);
        printf("%8.3f%8.3f%8.3f", // middle of pdb ATOM line
            positions[atomIndex][0] * OpenMM::AngstromsPerNm,  // x
            positions[atomIndex][1] * OpenMM::AngstromsPerNm,  // y
            positions[atomIndex][2] * OpenMM::AngstromsPerNm); // z
        printf("  1.00  0.00            \n"); // end of pdb ATOM line

    }
    printf("ENDMDL\n"); // end of trajectory frame
}
