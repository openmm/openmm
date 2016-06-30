/* -------------------------------------------------------------------------- *
 *                                OpenMMAmoeba                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2016 Stanford University and the Authors.      *
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
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaSerializationProxies();

static void getCovalentTypes(std::vector<std::string>& covalentTypes) {

    covalentTypes.push_back("Covalent12");
    covalentTypes.push_back("Covalent13");
    covalentTypes.push_back("Covalent14");
    covalentTypes.push_back("Covalent15");

    covalentTypes.push_back("PolarizationCovalent11");
    covalentTypes.push_back("PolarizationCovalent12");
    covalentTypes.push_back("PolarizationCovalent13");
    covalentTypes.push_back("PolarizationCovalent14");
}

void testSerialization() {
    // Create a Force.

    AmoebaMultipoleForce force1;
    force1.setForceGroup(3);
    force1.setNonbondedMethod(AmoebaMultipoleForce::NoCutoff);
    force1.setCutoffDistance(0.9);
    force1.setAEwald(0.544);
    //force1.setPmeBSplineOrder(4);

    std::vector<int> gridDimension;
    gridDimension.push_back(64);
    gridDimension.push_back(63);
    gridDimension.push_back(61);
    force1.setPmeGridDimensions(gridDimension); 
    //force1.setMutualInducedIterationMethod(AmoebaMultipoleForce::SOR); 
    force1.setMutualInducedMaxIterations(200); 
    force1.setMutualInducedTargetEpsilon(1.0e-05); 
    //force1.setElectricConstant(138.93); 
    force1.setEwaldErrorTolerance(1.0e-05); 
    
    vector<double> coeff;
    coeff.push_back(0.0);
    coeff.push_back(-0.1);
    coeff.push_back(1.1);
    force1.setExtrapolationCoefficients(coeff);

    std::vector<std::string> covalentTypes;
    getCovalentTypes(covalentTypes);

    for (unsigned int ii = 0; ii < 3; ii++) {
        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;
        molecularDipole.push_back(0.1); molecularDipole.push_back(rand()); molecularDipole.push_back(rand());
        for (unsigned int jj = 0; jj < 9; jj++) {
            molecularQuadrupole.push_back(static_cast<double>(rand()));
        }
        force1.addMultipole(static_cast<double>(ii+1), molecularDipole, molecularQuadrupole, AmoebaMultipoleForce::Bisector,
                            ii+1, ii+2, ii+3, static_cast<double>(rand()), static_cast<double>(rand()), static_cast<double>(rand()));

        for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
            std::vector< int > covalentMap;
            covalentMap.push_back(ii*jj); covalentMap.push_back(rand()); covalentMap.push_back(rand());
            force1.setCovalentMap(ii, static_cast<AmoebaMultipoleForce::CovalentType>(jj), covalentMap);
        }
    }

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<AmoebaMultipoleForce>(&force1, "Force", buffer);

    AmoebaMultipoleForce* copy = XmlSerializer::deserialize<AmoebaMultipoleForce>(buffer);

    // Compare the two forces to see if they are identical.  
    AmoebaMultipoleForce& force2 = *copy;

    ASSERT_EQUAL(force1.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force1.getCutoffDistance(),                force2.getCutoffDistance());
    ASSERT_EQUAL(force1.getNonbondedMethod(),               force2.getNonbondedMethod());
    ASSERT_EQUAL(force1.getAEwald(),                        force2.getAEwald());
    ASSERT_EQUAL(force1.getMutualInducedMaxIterations(),    force2.getMutualInducedMaxIterations());
    ASSERT_EQUAL(force1.getMutualInducedTargetEpsilon(),    force2.getMutualInducedTargetEpsilon());
    ASSERT_EQUAL(force1.getEwaldErrorTolerance(),           force2.getEwaldErrorTolerance());


    std::vector<int> gridDimension1;
    std::vector<int> gridDimension2;
    force1.getPmeGridDimensions(gridDimension1); 
    force2.getPmeGridDimensions(gridDimension2); 
    ASSERT_EQUAL(gridDimension1.size(),  gridDimension2.size());
    for (unsigned int jj = 0; jj < gridDimension1.size(); jj++) {
        ASSERT_EQUAL(gridDimension1[jj], gridDimension2[jj]);
    }
    
    ASSERT_EQUAL_CONTAINERS(force1.getExtrapolationCoefficients(), force2.getExtrapolationCoefficients());
    
    ASSERT_EQUAL(force1.getNumMultipoles(),  force2.getNumMultipoles());
    for (unsigned int ii = 0; ii < static_cast<unsigned int>(force1.getNumMultipoles()); ii++) {

        int axisType1, multipoleAtomZ1, multipoleAtomX1, multipoleAtomY1;
        int axisType2, multipoleAtomZ2, multipoleAtomX2, multipoleAtomY2;

        double charge1, thole1, dampingFactor1, polarity1;
        double charge2, thole2, dampingFactor2, polarity2;

        std::vector<double> molecularDipole1;
        std::vector<double> molecularQuadrupole1;

        std::vector<double> molecularDipole2;
        std::vector<double> molecularQuadrupole2;

        force1.getMultipoleParameters(ii, charge1, molecularDipole1, molecularQuadrupole1, axisType1, multipoleAtomZ1, multipoleAtomX1, multipoleAtomY1,
                                       thole1, dampingFactor1, polarity1);

        force2.getMultipoleParameters(ii, charge2, molecularDipole2, molecularQuadrupole2, axisType2, multipoleAtomZ2, multipoleAtomX2, multipoleAtomY2,
                                       thole2, dampingFactor2, polarity2);

        ASSERT_EQUAL(charge1,                        charge2);
        ASSERT_EQUAL(axisType1,                      axisType2);
        ASSERT_EQUAL(multipoleAtomZ1,                multipoleAtomZ2);
        ASSERT_EQUAL(multipoleAtomX1,                multipoleAtomX2);
        ASSERT_EQUAL(multipoleAtomY1,                multipoleAtomY2);
        ASSERT_EQUAL(thole1,                         thole2);
        ASSERT_EQUAL(dampingFactor1,                 dampingFactor2);
        ASSERT_EQUAL(polarity1,                      polarity2);

        ASSERT_EQUAL(molecularDipole1.size(),        molecularDipole2.size());
        ASSERT_EQUAL(molecularDipole1.size(),        3);
        for (unsigned int jj = 0; jj < molecularDipole1.size(); jj++) {
            ASSERT_EQUAL(molecularDipole1[jj], molecularDipole2[jj]);
        }
        ASSERT_EQUAL(molecularQuadrupole1.size(),        molecularQuadrupole2.size());
        ASSERT_EQUAL(molecularQuadrupole1.size(),        9);
        for (unsigned int jj = 0; jj < molecularQuadrupole1.size(); jj++) {
            ASSERT_EQUAL(molecularQuadrupole1[jj], molecularQuadrupole2[jj]);
        }

        for (unsigned int jj = 0; jj < covalentTypes.size(); jj++) {
            std::vector<int> covalentMap1;
            std::vector<int> covalentMap2;
            force1.getCovalentMap(ii, static_cast<AmoebaMultipoleForce::CovalentType>(jj), covalentMap1);
            force2.getCovalentMap(ii, static_cast<AmoebaMultipoleForce::CovalentType>(jj), covalentMap2);
            ASSERT_EQUAL(covalentMap1.size(),        covalentMap2.size());
            for (unsigned int kk = 0; kk < covalentMap1.size(); kk++) {
                ASSERT_EQUAL(covalentMap1[kk],        covalentMap2[kk]);
            }
        }
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

