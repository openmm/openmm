/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2023 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Emilio Gallicchio                                  *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/ATMForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void test2Particles() {
  // A pair of particles tethered by an harmonic bond.
  // Displace the second one to test energy and forces at lambda=1/2

  System system;
  system.addParticle(1.0);
  system.addParticle(1.0);

  double lmbd = 0.5;
  double lambda1 = lmbd;
  double lambda2 = lmbd;
  double alpha = 0.0;
  double u0 = 0.0;
  double w0 = 0.0;
  double umax = 1.e6;
  double ubcore= 0.5*1.e6;
  double acore = 1./16.;
  double direction = 1.0;

  double displx = 1.0;
  double disply = 0.;
  double displz = 0.;

  vector<Vec3> positions(2);
  positions[0] = Vec3(0, 0, 0);
  positions[1] = Vec3(0, 0, 0);

  CustomBondForce* bond = new CustomBondForce("0.5*r^2");
  bond->addBond(0, 1);

  ATMForce* atm = new ATMForce(lambda1, lambda2, alpha, u0, w0, umax, ubcore, acore, direction);
  atm->addParticle(0, 0., 0., 0.);
  atm->addParticle(1, displx, disply, displz);
  atm->addForce(bond);
  system.addForce(atm);

  VerletIntegrator integrator(1.0);
  Context context(system, integrator, platform);
  context.setPositions(positions);

  State state = context.getState(State::Energy | State::Forces );
  double epot = state.getPotentialEnergy();
  double epert = atm->getPerturbationEnergy(context);

  ASSERT_EQUAL_TOL(lmbd, context.getParameter(atm->Lambda1()), 1e-6);
  ASSERT_EQUAL_TOL(lmbd, context.getParameter(atm->Lambda2()), 1e-6);
  ASSERT_EQUAL_TOL(lmbd*0.5*displx*displx, epot, 1e-6);
  ASSERT_EQUAL_TOL(0.5*displx*displx, epert, 1e-6);
  ASSERT_EQUAL_VEC(Vec3(-lmbd*displx, 0.0, 0.0), state.getForces()[1], 1e-6);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        test2Particles();
        // Insert tests here
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

