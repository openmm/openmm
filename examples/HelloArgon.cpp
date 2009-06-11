// -----------------------------------------------------------------------------
//           OpenMM(tm) HelloArgon example in C++ (June 2009)
// -----------------------------------------------------------------------------
// This program demonstrates a simple molecular simulation using the OpenMM
// API for GPU-accelerated molecular dynamics simulation.  The system
// represents a small collection of argon atoms in a vaccuum.
// A multi-frame PDB file is written to stdout which  can be read by VMD or 
// other visualization tool to produce an animation of the resulting trajectory.
// -----------------------------------------------------------------------------

#include "OpenMM.h"
#include <cstdio>

// Forward declaration of writePdb() subroutine for printing atomic 
// coordinates, defined later in this source file.
void writePdb(const OpenMM::OpenMMContext& context);

void simulateArgon()
{
    // Load any shared libraries containing GPU implementations
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // Create a system with nonbonded forces
    OpenMM::System system;
    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce(); 
    system.addForce(nonbond);

    // Create three atoms
    std::vector<OpenMM::Vec3> initialPositions(3);
    for (int a = 0; a < 3; ++a) 
    {
        system.addParticle(39.95); // mass, grams per mole

        // charge, sigma, well depth
        nonbond->addParticle(0.0, 0.3350, 0.996); 

        initialPositions[a] = OpenMM::Vec3(0.5*a,0,0); // location, nanometers
    }

    OpenMM::VerletIntegrator integrator(0.004); // step size in picoseconds

    // Let OpenMM Context choose best platform.
    OpenMM::OpenMMContext context(system, integrator);
    printf( "REMARK  Using OpenMM platform %s\n", 
        context.getPlatform().getName().c_str() );

    // Set the starting positions of the atoms. Velocities will be zero.
    context.setPositions(initialPositions);

    // Simulate
    while(context.getTime() < 10.0) { // picoseconds
        writePdb(context); // output coordinates
        // Run 100 steps at a time, for efficient use of OpenMM
        integrator.step(100);
    }
    writePdb(context); // output final coordinates
}

int main() 
{
    try {
        simulateArgon();
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}

// writePdb() subroutine for quick-and-dirty trajectory output.
void writePdb(const OpenMM::OpenMMContext& context) 
{
    // Request atomic positions from OpenMM
    const OpenMM::State state = context.getState(OpenMM::State::Positions);
    const std::vector<OpenMM::Vec3>& pos = state.getPositions();

    // write out in PDB format

    // Use PDB MODEL cards to number trajectory frames
    static int modelFrameNumber = 0; 
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber); // start of frame
    for (int a = 0; a < context.getSystem().getNumParticles(); ++a)
    {
        printf("ATOM  %5d AR    AR     1    ", a+1); // atom number
        printf("%8.3f%8.3f%8.3f  1.00  0.00          AR\n", // coordinates
            // "*10" converts nanometers to Angstroms
            pos[a][0]*10, pos[a][1]*10, pos[a][2]*10);
    }
    printf("ENDMDL\n"); // end of frame
}
