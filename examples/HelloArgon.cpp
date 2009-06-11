#include "OpenMM.h"
#include <cstdio>

using namespace OpenMM;

// forward declaration of subroutine defined later in this file.
void writePdb(const OpenMMContext& context);

void simulateArgon() 
{
    // Load any shared libraries containing GPU implementations
    Platform::loadPluginsFromDirectory(
        Platform::getDefaultPluginsDirectory());

    System system;
    NonbondedForce* nonbond = new NonbondedForce(); 
    system.addForce(nonbond);

    // Create three atoms
    std::vector<Vec3> initialPositions(3);
    for (int a = 0; a < 3; ++a) 
    {
        system.addParticle(39.95); // mass
        nonbond->addParticle(0.0, 0.3350, 0.001603); // charge, sigma, well depth
        initialPositions[a] = Vec3(0.5*a,0,0); // location
    }

    VerletIntegrator integrator(0.020); // step size in picoseconds

    // Let OpenMM Context choose best platform.
    OpenMMContext context(system, integrator);
    printf( "REMARK  Using OpenMM platform %s\n", 
        context.getPlatform().getName().c_str() );

    // Set the starting positions of the atoms. Velocities will be zero.
    context.setPositions(initialPositions);

    // Simulate
    while(context.getTime() < 500.0) { // picoseconds
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
void writePdb(const OpenMMContext& context) 
{
    // Request atomic positions from OpenMM
    const State state = context.getState(State::Positions);
    const std::vector<Vec3>& pos = state.getPositions();

    // write out in PDB format

    // Use PDB MODEL cards to number trajectory frames
    static int modelFrameNumber = 0; 
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber); // start of frame
    for (int a = 0; a < context.getSystem().getNumParticles(); ++a)
    {
        printf("ATOM  %5d AR    AR     1    ", a+1); // atom number
        printf("%8.3f%8.3f%8.3f  1.00  0.00          AR\n",
            // notice "*10" converts nanometers to Angstroms
            pos[a][0]*10, pos[a][1]*10, pos[a][2]*10);
    }
    printf("ENDMDL\n"); // end of frame
}
