/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/Context.h"
#include "ReferencePlatform.h"
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/RBTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <limits>
#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-4;

/**
 * This tests a system with multiple forces, to make sure BondedUtilities is
 * processing them correctly.
 */
void testForces() {
    const int numParticles = 100;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    for (int i = 0; i < numParticles-1; i++)
        bonds->addBond(i, i+1, 1.0, 1.5);
    system.addForce(bonds);
    HarmonicAngleForce* angles = new HarmonicAngleForce();
    for (int i = 0; i < numParticles-2; i++)
        angles->addAngle(i, i+1, i+2, 2.0, 1.5);
    system.addForce(angles);
    PeriodicTorsionForce* periodic = new PeriodicTorsionForce();
    for (int i = 0; i < numParticles-3; i++)
        periodic->addTorsion(i, i+1, i+2, i+3, 2, PI_M/3, 1.1);
    system.addForce(periodic);
    RBTorsionForce* rb = new RBTorsionForce();
    for (int i = 0; i < numParticles-3; i += 3)
        rb->addTorsion(i, i+1, i+2, i+3, 1.0, 1.1, 1.2, 0.3, 0.4, 0.5);
    system.addForce(rb);
    ReferencePlatform ref;
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, ref);
    VerletIntegrator integrator2(0.01);
    Context context2(system, integrator2, platform);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(i, genrand_real2(sfmt), genrand_real2(sfmt));
    context1.setPositions(positions);
    context2.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    const vector<Vec3>& forces1 = state1.getForces();
    const vector<Vec3>& forces2 = state2.getForces();
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(forces1[i], forces2[i], TOL);
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy(), state2.getPotentialEnergy(), TOL);
}

void testComputeContextValidation() {
    // Verify molecule grouping when particle-group sizes differ (both p1 < p2 and p1 > p2)
    // and .y/.z NaN centroid detection during molecule reordering.
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(6.0, 0, 0), Vec3(0, 6.0, 0), Vec3(0, 0, 6.0));
    NonbondedForce* nb = new NonbondedForce();
    nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nb->setCutoffDistance(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    CustomCentroidBondForce* centroid = new CustomCentroidBondForce(2, "100.0*distance(g1, g2)^2");
    system.addForce(nb);
    system.addForce(bonds);
    system.addForce(centroid);
    vector<Vec3> positions;
    const double masses[5] = {12.0, 12.0, 12.0, 14.0, 14.0};
    const double baseX[5] = {4.5, 0.5, 5.2, 2.0, 3.5};
    const bool useThreeAtoms[5] = {true, false, false, false, true};
    for (int m = 0; m < 5; m++) {
        int base = system.getNumParticles();
        for (int a = 0; a < 3; a++) {
            system.addParticle(masses[m]);
            nb->addParticle(0.0, 0.25, 0.0);
            positions.push_back(Vec3(baseX[m], 1.0 + 0.2*(a > 0 ? 1 : 0), 1.0 + (a == 2 ? 0.5 + 0.05*m : 0.0)));
            if (a > 0) {
                bonds->addBond(base+a-1, base+a, 0.2, 10.0);
                nb->addException(base+a-1, base+a, 0.0, 1.0, 0.0);
            }
        }
        int g1 = centroid->addGroup({base});
        int g2 = useThreeAtoms[m] ? centroid->addGroup({base+1, base+2}) : centroid->addGroup({base+1});
        centroid->addBond({g1, g2});
    }
    ReferencePlatform ref;
    VerletIntegrator int1(0.001), int2(0.001);
    Context ctxRef(system, int1, ref), ctx(system, int2, platform);
    ctxRef.setPositions(positions);
    ctx.setPositions(positions);
    State stateRef = ctxRef.getState(State::Forces | State::Energy);
    State state = ctx.getState(State::Forces | State::Energy);
    ASSERT_EQUAL_TOL(stateRef.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
    for (int i = 0; i < system.getNumParticles(); i++)
        ASSERT_EQUAL_VEC(stateRef.getForces()[i], state.getForces()[i], TOL);

    for (int axis = 1; axis <= 2; axis++) {
        vector<Vec3> nanPos = positions;
        nanPos[4][axis] = numeric_limits<double>::quiet_NaN();
        VerletIntegrator nanInt(0.001);
        Context nanCtx(system, nanInt, platform);
        bool threw = false;
        try {
            nanCtx.setPositions(nanPos);
        }
        catch (const OpenMMException&) {
            threw = true;
        }
        ASSERT(threw);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testForces();
        testComputeContextValidation();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
