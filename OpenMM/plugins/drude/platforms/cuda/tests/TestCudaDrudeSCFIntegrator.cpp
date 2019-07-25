/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the CUDA implementation of DrudeSCFIntegrator.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/DrudeForce.h"
#include "openmm/DrudeSCFIntegrator.h"
#include "SimTKOpenMMUtilities.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerDrudeCudaKernelFactories();

void testWater() {
    // Create a box of SWM4-NDP water molecules.  This involves constraints, virtual sites,
    // and Drude particles.
    
    const int gridSize = 4;
    const int numMolecules = gridSize*gridSize*gridSize;
    const double spacing = 0.6;
    const double boxSize = spacing*(gridSize+1);
    System system;
    NonbondedForce* nonbonded = new NonbondedForce();
    DrudeForce* drude = new DrudeForce();
    system.addForce(nonbonded);
    system.addForce(drude);
    system.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    for (int i = 0; i < numMolecules; i++) {
        int startIndex = system.getNumParticles();
        system.addParticle(15.6); // O
        system.addParticle(0.4);  // D
        system.addParticle(1.0);  // H1
        system.addParticle(1.0);  // H2
        system.addParticle(0.0);  // M
        nonbonded->addParticle(1.71636, 0.318395, 0.21094*4.184);
        nonbonded->addParticle(-1.71636, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(0.55733, 1, 0);
        nonbonded->addParticle(-1.11466, 1, 0);
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < j; k++)
                nonbonded->addException(startIndex+j, startIndex+k, 0, 1, 0);
        system.addConstraint(startIndex, startIndex+2, 0.09572);
        system.addConstraint(startIndex, startIndex+3, 0.09572);
        system.addConstraint(startIndex+2, startIndex+3, 0.15139);
        system.setVirtualSite(startIndex+4, new ThreeParticleAverageSite(startIndex, startIndex+2, startIndex+3, 0.786646558, 0.106676721, 0.106676721));
        drude->addParticle(startIndex+1, startIndex, -1, -1, -1, -1.71636, ONE_4PI_EPS0*1.71636*1.71636/(100000*4.184), 1, 1);
    }
    vector<Vec3> positions;
    for (int i = 0; i < gridSize; i++)
        for (int j = 0; j < gridSize; j++)
            for (int k = 0; k < gridSize; k++) {
                Vec3 pos(i*spacing, j*spacing, k*spacing);
                positions.push_back(pos);
                positions.push_back(pos);
                positions.push_back(pos+Vec3(0.09572, 0, 0));
                positions.push_back(pos+Vec3(-0.023999, 0.092663, 0));
                positions.push_back(pos);
            }
    
    // Simulate it and check energy conservation and the total force on the Drude particles.
    
    DrudeSCFIntegrator integ(0.0005);
    Platform& platform = Platform::getPlatformByName("CUDA");
    Context context(system, integ, platform);
    context.setPositions(positions);
    context.applyConstraints(1e-5);
    context.setVelocitiesToTemperature(300.0);
    State state = context.getState(State::Energy);
    double initialEnergy;
    int numSteps = 1000;
    double maxNorm = (platform.getPropertyValue(context, "Precision") == "double" ? 1.0 : 5.0);
    for (int i = 0; i < numSteps; i++) {
        integ.step(1);
        state = context.getState(State::Energy | State::Forces);
        if (i == 0)
            initialEnergy = state.getPotentialEnergy()+state.getKineticEnergy();
        else
            ASSERT_EQUAL_TOL(initialEnergy, state.getPotentialEnergy()+state.getKineticEnergy(), 0.01);
        const vector<Vec3>& force = state.getForces();
        double norm = 0.0;
        for (int j = 1; j < (int) force.size(); j += 5)
            norm += sqrt(force[j].dot(force[j]));
        norm = (norm/numMolecules);
        ASSERT(norm < maxNorm);
    }
}

int main(int argc, char* argv[]) {
    try {
        registerDrudeCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", string(argv[1]));
        testWater();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
