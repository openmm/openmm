// ((6.7929, 0.0, 0.0), (-2.264163559406279, 6.404455775962287, 0.0), (-2.264163559406279, -3.2019384603140684, 5.54658849047036))
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void testTrunactedOctahedron() {
    const int numMolecules = 5;
    const int numParticles = numMolecules*2;
    const float cutoff = 2.0;
    Vec3 a(6.7929, 0, 0);
    Vec3 b(-2.264163559406279, 6.404455775962287, 0);
    Vec3 c(-2.264163559406279, -3.2019384603140684, 5.54658849047036);
        
    System system;
    system.setDefaultPeriodicBoxVectors(a, b, c);
    NonbondedForce* force = new NonbondedForce();
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    
    force->setCutoffDistance(cutoff);
    force->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    
    for (int i = 0; i < numMolecules; i++) {
        system.addParticle(1.0);
        system.addParticle(1.0);
        force->addParticle(-1, 0.2, 0.2);
        force->addParticle(1, 0.2, 0.2);
        positions[2*i] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        positions[2*i+1] = positions[2*i] + Vec3(1.0, 0.0, 0.0);
        system.addConstraint(2*i, 2*i+1, 1.0);
    }
    system.addForce(force);
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, Platform::getPlatformByName("Reference"));
    context.setPositions(positions);
    State initialState = context.getState(State::Positions | State::Energy, true);
    double initialEnergy = initialState.getPotentialEnergy();

    context.setState(initialState);
    State finalState = context.getState(State::Positions | State::Energy, true);
    double finalEnergy = finalState.getPotentialEnergy();

    ASSERT_EQUAL_TOL(initialEnergy, finalEnergy, 1e-4);
}

int main(int argc, char* argv[]) {
    try {
        testTrunactedOctahedron();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
