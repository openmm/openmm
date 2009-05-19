#include "OpenMM.h"
#include <cstdio>

using namespace OpenMM;

void writePdb(const OpenMMContext& context) {
    const State state = context.getState(State::Positions);
    const std::vector<Vec3>& pos = state.getPositions();
    static int modelFrameNumber = 0; // numbering for MODEL records in pdb output
    // write out in PDB format
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    for (int a = 0; a < context.getSystem().getNumParticles(); ++a)
        printf("ATOM  %5d AR    AR     1    %8.3f%8.3f%8.3f  1.00  0.00          AR\n",
            a+1, pos[a][0]*10, pos[a][1]*10, pos[a][2]*10);
    printf("ENDMDL\n");
}

int main() {
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
    System system;
    NonbondedForce* nonbond = new NonbondedForce(); 
    system.addForce(nonbond);

    // Create atoms
    int numAtoms = 3;
    for (int a = 0; a < numAtoms; ++a) {
        system.addParticle(39.95); // mass
        nonbond->addParticle(0.0, 0.3350, 0.001603); // charge, diameter, well depth
    }

    // Large step size may be required for stability with small forces
    VerletIntegrator integrator(0.020); // step size in picoseconds
    // Let OpenMM Context choose best platform.
    OpenMMContext context(system, integrator);
    printf( "REMARK  Using OpenMM platform %s\n", context.getPlatform().getName().c_str() );

    // Set the starting positions in the OpenMM context. Velocities will be zero.
    std::vector<Vec3> positions(numAtoms);
    for (int a = 0; a < numAtoms; ++a)
        positions[a] = Vec3(0.5*a,0,0);
    context.setPositions(positions);

    while(context.getTime() < 500.0) { // picoseconds
        writePdb(context);
        integrator.step(100); // number of steps between reports
    }
    writePdb(context);

    return 0;
}
