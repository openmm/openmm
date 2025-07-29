// -----------------------------------------------------------------------------
//           OpenMM(tm) HelloArgon example in C++ (June 2009)
// -----------------------------------------------------------------------------
// This program demonstrates a simple molecular simulation using the OpenMM
// API for GPU-accelerated molecular dynamics simulation. The primary goal is
// to make sure you can compile, link, and run with OpenMM and view the output.
// The example is available in C++, C, and Fortran 95.
//
// The system modeled here is a small number of argon atoms in a vacuum.
// A multi-frame PDB file is written to stdout which  can be read by VMD or 
// other visualization tool to produce an animation of the resulting trajectory.
// -----------------------------------------------------------------------------

#include "OpenMM.h"
#include <cstdio>

// Forward declaration of routine for printing one frame of the
// trajectory, defined later in this source file.
void writePdbFrame(int frameNum, const OpenMM::State&);

void simulateArgon()
{
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // Create a system with nonbonded forces.
    OpenMM::System system;
    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce(); 
    system.addForce(nonbond);

    // Create three atoms.
    std::vector<OpenMM::Vec3> initPosInNm(3);
    for (int a = 0; a < 3; ++a) 
    {
        initPosInNm[a] = OpenMM::Vec3(0.5*a,0,0); // location, nm

        system.addParticle(39.95); // mass of Ar, grams per mole

        // charge, L-J sigma (nm), well depth (kJ)
        nonbond->addParticle(0.0, 0.3350, 0.996); // vdWRad(Ar)=.188 nm
    }

    OpenMM::VerletIntegrator integrator(0.004); // step size in ps

    // Let OpenMM Context choose best platform.
    OpenMM::Context context(system, integrator);
    printf( "REMARK  Using OpenMM platform %s\n", 
        context.getPlatform().getName().c_str() );

    // Set starting positions of the atoms. Leave time and velocity zero.
    context.setPositions(initPosInNm);

    // Simulate.
    for (int frameNum=1; ;++frameNum) {
        // Output current state information.
        OpenMM::State state    = context.getState(OpenMM::State::Positions);
        const double  timeInPs = state.getTime();
        writePdbFrame(frameNum, state); // output coordinates

        if (timeInPs >= 10.)
            break;

        // Advance state many steps at a time, for efficient use of OpenMM.
        integrator.step(10); // (use a lot more than this normally)
    }
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

// Handy homebrew PDB writer for quick-and-dirty trajectory output.
void writePdbFrame(int frameNum, const OpenMM::State& state) 
{
    // Reference atomic positions in the OpenMM State.
    const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

    // Use PDB MODEL cards to number trajectory frames
    printf("MODEL     %d\n", frameNum); // start of frame
    for (int a = 0; a < (int)posInNm.size(); ++a)
    {
        printf("ATOM  %5d  AR   AR     1    ", a+1); // atom number
        printf("%8.3f%8.3f%8.3f  1.00  0.00\n",      // coordinates
            // "*10" converts nanometers to Angstroms
            posInNm[a][0]*10, posInNm[a][1]*10, posInNm[a][2]*10);
    }
    printf("ENDMDL\n"); // end of frame
}
