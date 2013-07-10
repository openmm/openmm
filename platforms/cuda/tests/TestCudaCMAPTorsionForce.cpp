/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2012 Stanford University and the Authors.      *
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
 * This tests the CUDA implementation of CMAPTorsionForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "CudaPlatform.h"
#include "openmm/CMAPTorsionForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

CudaPlatform platform;

void testCMAPTorsions() {
    const int mapSize = 36;

    // Create two systems: one with a pair of periodic torsions, and one with a CMAP torsion
    // that approximates the same force.

    System system1;
    for (int i = 0; i < 5; i++)
        system1.addParticle(1.0);
    PeriodicTorsionForce* periodic = new PeriodicTorsionForce();
    periodic->addTorsion(0, 1, 2, 3, 2, M_PI/4, 1.5);
    periodic->addTorsion(1, 2, 3, 4, 3, M_PI/3, 2.0);
    system1.addForce(periodic);
    System system2;
    for (int i = 0; i < 5; i++)
        system2.addParticle(1.0);
    CMAPTorsionForce* cmap = new CMAPTorsionForce();
    vector<double> mapEnergy(mapSize*mapSize);
    for (int i = 0; i < mapSize; i++) {
        double angle1 = i*2*M_PI/mapSize;
        double energy1 = 1.5*(1+cos(2*angle1-M_PI/4));
        for (int j = 0; j < mapSize; j++) {
            double angle2 = j*2*M_PI/mapSize;
            double energy2 = 2.0*(1+cos(3*angle2-M_PI/3));
            mapEnergy[i+j*mapSize] = energy1+energy2;
        }
    }
    cmap->addMap(mapSize, mapEnergy);
    cmap->addTorsion(0, 0, 1, 2, 3, 1, 2, 3, 4);
    system2.addForce(cmap);

    // Set the atoms in various positions, and verify that both systems give equal forces and energy.

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(5);
    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);
    Context c1(system1, integrator1, platform);
    Context c2(system2, integrator2, platform);
    for (int i = 0; i < 50; i++) {
        for (int j = 0; j < (int) positions.size(); j++)
            positions[j] = Vec3(5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt), 5.0*genrand_real2(sfmt));
        c1.setPositions(positions);
        c2.setPositions(positions);
        State s1 = c1.getState(State::Forces | State::Energy);
        State s2 = c2.getState(State::Forces | State::Energy);
        for (int i = 0; i < system1.getNumParticles(); i++)
            ASSERT_EQUAL_VEC(s1.getForces()[i], s2.getForces()[i], 0.05);
        ASSERT_EQUAL_TOL(s1.getPotentialEnergy(), s2.getPotentialEnergy(), 1e-3);
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc > 1)
            platform.setPropertyDefaultValue("CudaPrecision", string(argv[1]));
        testCMAPTorsions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
