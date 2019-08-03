/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

/**
 * Test the methods for manipulating System objects.
 */
void testCreateSystem() {
    int numParticles = 10;
    System system;
    
    // Test adding and modifying particles.
    
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0+0.1*i);
    system.setParticleMass(5, 100.0);
    ASSERT_EQUAL(numParticles, system.getNumParticles());
    for (int i = 0; i < numParticles; i++) {
        double mass = (i == 5 ? 100.0 : 1.0+0.1*i);
        ASSERT_EQUAL_TOL(mass, system.getParticleMass(i), 1e-15);
    }
    
    // Test adding, removing, and modifying constraints.
    
    for (int i = 0; i < numParticles-1; i++)
        system.addConstraint(i, i+1, 0.2*i);
    system.removeConstraint(5);
    system.setConstraintParameters(3, 0, 5, 99.0);
    ASSERT_EQUAL(numParticles-2, system.getNumConstraints());
    for (int i = 0; i < numParticles-2; i++) {
        int p1, p2;
        double dist;
        system.getConstraintParameters(i, p1, p2, dist);
        if (i == 3) {
            ASSERT_EQUAL(0, p1);
            ASSERT_EQUAL(5, p2);
            ASSERT_EQUAL(99.0, dist);
        }
        else {
            int j = (i < 5 ? i : i+1);
            ASSERT_EQUAL(j, p1);
            ASSERT_EQUAL(j+1, p2);
            ASSERT_EQUAL(0.2*j, dist);
        }
    }
    
    // Test adding and removing forces.
    
    HarmonicBondForce* bonds = new HarmonicBondForce();
    system.addForce(bonds);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    system.addForce(angles);
    ASSERT_EQUAL(2, system.getNumForces());
    ASSERT_EQUAL(bonds, &system.getForce(0));
    ASSERT_EQUAL(angles, &system.getForce(1));
    system.removeForce(0);
    ASSERT_EQUAL(1, system.getNumForces());
    ASSERT_EQUAL(angles, &system.getForce(0));
    
    // Test adding and removing virtual sites.
    
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL(false, system.isVirtualSite(i));
    }
    TwoParticleAverageSite* site = new TwoParticleAverageSite(2, 3, 0.4, 0.6);
    system.setVirtualSite(4, site);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL(i == 4, system.isVirtualSite(i));
    }
    ASSERT_EQUAL(site, &system.getVirtualSite(4));
    system.setVirtualSite(4, NULL);
    for (int i = 0; i < numParticles; i++) {
        ASSERT_EQUAL(false, system.isVirtualSite(i));
    }
}

int main() {
    try {
        testCreateSystem();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
