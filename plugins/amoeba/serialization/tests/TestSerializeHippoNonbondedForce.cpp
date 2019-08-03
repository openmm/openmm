/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/HippoNonbondedForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

void testSerialization() {
    // Create a Force.

    HippoNonbondedForce force1;
    force1.setForceGroup(3);
    force1.setNonbondedMethod(HippoNonbondedForce::PME);
    force1.setCutoffDistance(0.7);
    force1.setSwitchingDistance(0.6);
    force1.setEwaldErrorTolerance(1.0e-5); 
    force1.setPMEParameters(0.5, 20, 22, 24);
    force1.setDPMEParameters(0.4, 15, 16, 18);
    force1.setExtrapolationCoefficients({0.0, -0.1, 1.1});
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < 3; i++) {
        vector<double> dipole = {genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)};
        vector<double> quadrupole = {genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt),
                                     genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt),
                                     genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)};
        force1.addParticle(genrand_real2(sfmt), dipole, quadrupole, genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt),
                           genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt),
                           genrand_real2(sfmt), HippoNonbondedForce::Bisector, i+1, i+2, i+3);
    }
    for (int i = 0; i < 2; i++)
        force1.addException(i, i+1, genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt));

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<HippoNonbondedForce>(&force1, "Force", buffer);
    HippoNonbondedForce* copy = XmlSerializer::deserialize<HippoNonbondedForce>(buffer);

    // Compare the two forces to see if they are identical.  

    HippoNonbondedForce& force2 = *copy;
    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force1.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force1.getSwitchingDistance(), force2.getSwitchingDistance());
    ASSERT_EQUAL(force1.getEwaldErrorTolerance(), force2.getEwaldErrorTolerance());
    double alpha1, alpha2;
    int nx1, nx2, ny1, ny2, nz1, nz2;
    force1.getPMEParameters(alpha1, nx1, ny1, nz1);
    force2.getPMEParameters(alpha2, nx2, ny2, nz2);
    ASSERT_EQUAL(alpha1, alpha2);
    ASSERT_EQUAL(nx1, nx2);
    ASSERT_EQUAL(ny1, ny2);
    ASSERT_EQUAL(nz1, nz2);
    force1.getDPMEParameters(alpha1, nx1, ny1, nz1);
    force2.getDPMEParameters(alpha2, nx2, ny2, nz2);
    ASSERT_EQUAL(alpha1, alpha2);
    ASSERT_EQUAL(nx1, nx2);
    ASSERT_EQUAL(ny1, ny2);
    ASSERT_EQUAL(nz1, nz2);
    ASSERT_EQUAL_CONTAINERS(force1.getExtrapolationCoefficients(), force2.getExtrapolationCoefficients());
    ASSERT_EQUAL(force1.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force1.getNumParticles(); i++) {
        double charge1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1, polarizability1;
        double charge2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2, polarizability2;
        int axisType1, atomX1, atomY1, atomZ1;
        int axisType2, atomX2, atomY2, atomZ2;
        vector<double> dipole1, quadrupole1;
        vector<double> dipole2, quadrupole2;
        force1.getParticleParameters(i, charge1, dipole1, quadrupole1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1,
                                    polarizability1, axisType1, atomZ1, atomX1, atomY1);
        force2.getParticleParameters(i, charge2, dipole2, quadrupole2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2,
                                    polarizability2, axisType2, atomZ2, atomX2, atomY2);
        ASSERT_EQUAL(charge1, charge2);
        ASSERT_EQUAL(coreCharge1, coreCharge2);
        ASSERT_EQUAL(alpha1, alpha2);
        ASSERT_EQUAL(epsilon1, epsilon2);
        ASSERT_EQUAL(damping1, damping2);
        ASSERT_EQUAL(c61, c62);
        ASSERT_EQUAL(pauliK1, pauliK2);
        ASSERT_EQUAL(pauliQ1, pauliQ2);
        ASSERT_EQUAL(pauliAlpha1, pauliAlpha2);
        ASSERT_EQUAL(polarizability1, polarizability2);
        ASSERT_EQUAL(axisType1, axisType2);
        ASSERT_EQUAL(atomX1, atomX2);
        ASSERT_EQUAL(atomY1, atomY2);
        ASSERT_EQUAL(atomZ1, atomZ2);
        ASSERT_EQUAL_CONTAINERS(dipole1, dipole2);
        ASSERT_EQUAL_CONTAINERS(quadrupole1, quadrupole2);
    }
    ASSERT_EQUAL(force1.getNumExceptions(), force2.getNumExceptions());
    for (int i = 0; i < force1.getNumExceptions(); i++) {
        int p11, p21;
        int p12, p22;
        double mm1, dm1, dd1, disp1, rep1, ct1;
        double mm2, dm2, dd2, disp2, rep2, ct2;
        force1.getExceptionParameters(i, p11, p21, mm1, dm1, dd1, disp1, rep1, ct1);
        force2.getExceptionParameters(i, p12, p22, mm2, dm2, dd2, disp2, rep2, ct2);
        ASSERT_EQUAL(p11, p12);
        ASSERT_EQUAL(p21, p22);
        ASSERT_EQUAL(mm1, mm2);
        ASSERT_EQUAL(dm1, dm2);
        ASSERT_EQUAL(dd1, dd2);
        ASSERT_EQUAL(disp1, disp2);
        ASSERT_EQUAL(rep1, rep2);
        ASSERT_EQUAL(ct1, ct2);
    }
}

int main() {
    try {
        registerAmoebaSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

