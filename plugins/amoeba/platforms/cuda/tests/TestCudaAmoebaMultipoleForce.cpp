/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
 * This tests the CUDA implementation of AmoebaMultipoleForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol, testname) {double _norm_ = std::sqrt(expected.dot(expected)); double _scale_ = _norm_ > 1.0 ? _norm_ : 1.0; if ((std::abs((expected[0])-(found[0]))/_scale_ > (tol)) || (std::abs((expected[1])-(found[1]))/_scale_ > (tol)) || (std::abs((expected[2])-(found[2]))/_scale_ > (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};


using namespace OpenMM;
using namespace std;

const double TOL = 1e-4;

extern "C" void registerAmoebaCudaKernelFactories();

// setup for 2 ammonia molecules

static void setupMultipoleAmmonia(System& system, AmoebaMultipoleForce* amoebaMultipoleForce, AmoebaMultipoleForce::NonbondedMethod nonbondedMethod,
                                  AmoebaMultipoleForce::PolarizationType polarizationType,
                                  double cutoff, int inputPmeGridDimension) {

    // box

    double boxDimension                               = 0.6;
    Vec3 a(boxDimension, 0.0, 0.0);
    Vec3 b(0.0, boxDimension, 0.0);
    Vec3 c(0.0, 0.0, boxDimension);
    system.setDefaultPeriodicBoxVectors(a, b, c);

    int numberOfParticles                             = 8;

    amoebaMultipoleForce->setNonbondedMethod(nonbondedMethod);
    amoebaMultipoleForce->setPolarizationType(polarizationType);
    amoebaMultipoleForce->setCutoffDistance(cutoff);
    amoebaMultipoleForce->setMutualInducedTargetEpsilon(1.0e-06);
    amoebaMultipoleForce->setMutualInducedMaxIterations(500);
    amoebaMultipoleForce->setAEwald(1.4024714e+01);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-04);

    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    std::vector<double> nitrogenMolecularDipole(3);
    std::vector<double> nitrogenMolecularQuadrupole(9);

    nitrogenMolecularDipole[0]     =   8.3832254e-03;
    nitrogenMolecularDipole[1]     =   0.0000000e+00;
    nitrogenMolecularDipole[2]     =   3.4232474e-03;

    nitrogenMolecularQuadrupole[0] =  -4.0406249e-04;
    nitrogenMolecularQuadrupole[1] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[2] =  -2.6883671e-04;
    nitrogenMolecularQuadrupole[3] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[4] =   2.5463927e-04;
    nitrogenMolecularQuadrupole[5] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[6] =  -2.6883671e-04;
    nitrogenMolecularQuadrupole[7] =   0.0000000e+00;
    nitrogenMolecularQuadrupole[8] =   1.4942322e-04;

    // first N

    system.addParticle(1.4007000e+01);
    amoebaMultipoleForce->addMultipole( -5.7960000e-01, nitrogenMolecularDipole, nitrogenMolecularQuadrupole, 2, 1, 2, 3,  3.9000000e-01,  3.1996314e-01,  1.0730000e-03);

    // 3 H attached to first N

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -1.7388763e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -4.6837475e-03;

    hydrogenMolecularQuadrupole[0] =  -4.4253841e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =   1.5429571e-05;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =   4.1798924e-05;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =   1.5429571e-05;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   2.4549167e-06;

    system.addParticle(1.0080000e+00);
    system.addParticle(1.0080000e+00);
    system.addParticle(1.0080000e+00);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 2, 3, 3.9e-01,  2.8135002e-01,  4.96e-04);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 1, 3, 3.9e-01,  2.8135002e-01,  4.96e-04);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 0, 1, 2, 3.9e-01,  2.8135002e-01,  4.96e-04);

    // second N

    system.addParticle(  1.4007000e+01);
    amoebaMultipoleForce->addMultipole( -5.796e-01, nitrogenMolecularDipole, nitrogenMolecularQuadrupole, 2, 5, 6, 7,  3.9e-01,  3.1996314e-01,  1.073e-03);

    // 3 H attached to second N

    system.addParticle(  1.0080000e+00);
    system.addParticle(  1.0080000e+00);
    system.addParticle(  1.0080000e+00);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 6, 7, 3.9e-01,  2.8135002e-01,  4.96e-04);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 5, 7, 3.9e-01,  2.8135002e-01,  4.96e-04);
    amoebaMultipoleForce->addMultipole(  1.932e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 2, 4, 5, 6, 3.9e-01,  2.8135002e-01,  4.96e-04);

    // covalent maps

    std::vector< int > covalentMap;
    covalentMap.resize(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(0, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(0, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    amoebaMultipoleForce->setCovalentMap(1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    amoebaMultipoleForce->setCovalentMap(2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(1);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    amoebaMultipoleForce->setCovalentMap(3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    amoebaMultipoleForce->setCovalentMap(3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(0);
    covalentMap.push_back(1);
    covalentMap.push_back(2);
    covalentMap.push_back(3);
    amoebaMultipoleForce->setCovalentMap(3, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(4, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(4, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    amoebaMultipoleForce->setCovalentMap(5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(5, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    amoebaMultipoleForce->setCovalentMap(6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(5);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(6, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    amoebaMultipoleForce->setCovalentMap(7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    amoebaMultipoleForce->setCovalentMap(7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(4);
    covalentMap.push_back(5);
    covalentMap.push_back(6);
    covalentMap.push_back(7);
    amoebaMultipoleForce->setCovalentMap(7, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
    system.addForce(amoebaMultipoleForce);
}

static void getForcesEnergyMultipoleAmmonia(Context& context, std::vector<Vec3>& forces, double& energy) {
    std::vector<Vec3> positions(context.getSystem().getNumParticles());

    positions[0]              = Vec3(  1.5927280e-01,  1.7000000e-06,   1.6491000e-03);
    positions[1]              = Vec3(  2.0805540e-01, -8.1258800e-02,   3.7282500e-02);
    positions[2]              = Vec3(  2.0843610e-01,  8.0953200e-02,   3.7462200e-02);
    positions[3]              = Vec3(  1.7280780e-01,  2.0730000e-04,  -9.8741700e-02);
    positions[4]              = Vec3( -1.6743680e-01,  1.5900000e-05,  -6.6149000e-03);
    positions[5]              = Vec3( -2.0428260e-01,  8.1071500e-02,   4.1343900e-02);
    positions[6]              = Vec3( -6.7308300e-02,  1.2800000e-05,   1.0623300e-02);
    positions[7]              = Vec3( -2.0426290e-01, -8.1231400e-02,   4.1033500e-02);

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

// compare forces and energies 

static void compareForcesEnergy(std::string& testName, double expectedEnergy, double energy,
                                 const std::vector<Vec3>& expectedForces,
                                 const std::vector<Vec3>& forces, double tolerance) {
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        ASSERT_EQUAL_VEC_MOD(expectedForces[ii], forces[ii], tolerance, testName);
    }
    ASSERT_EQUAL_TOL_MOD(expectedEnergy, energy, tolerance, testName);
}

// compare relative differences in force norms and energies 

static void compareForceNormsEnergy(std::string& testName, double expectedEnergy, double energy,
                                     std::vector<Vec3>& expectedForces,
                                     const std::vector<Vec3>& forces, double tolerance) {
    for (unsigned int ii = 0; ii < forces.size(); ii++) {
        double expectedNorm = sqrt(expectedForces[ii][0]*expectedForces[ii][0] +
                                    expectedForces[ii][1]*expectedForces[ii][1] +
                                    expectedForces[ii][2]*expectedForces[ii][2]);

        double norm         = sqrt(forces[ii][0]*forces[ii][0] + forces[ii][1]*forces[ii][1] + forces[ii][2]*forces[ii][2]);
        double absDiff      = fabs(norm - expectedNorm);
        double relDiff      = 2.0*absDiff/(fabs(norm) + fabs(expectedNorm) + 1.0e-08);

        if (relDiff > tolerance && absDiff > 0.001) {
            std::stringstream details;
            details << testName << "Relative difference in norms " << relDiff << " larger than allowed tolerance at particle=" << ii;
            details << ": norms=" << norm << " expected norm=" << expectedNorm; 
            throwException(__FILE__, __LINE__, details.str());
        }
    }
    double energyAbsDiff = fabs(expectedEnergy - energy);   
    double energyRelDiff =  2.0*energyAbsDiff/(fabs(expectedEnergy) + fabs(energy) + 1.0e-08);   
    if (energyRelDiff > tolerance) {
        std::stringstream details;
        details << testName << "Relative difference in energies " << energyRelDiff << " larger than allowed tolerance.";
        details << "Energies=" << energy << " expected energy=" << expectedEnergy; 
        throwException(__FILE__, __LINE__, details.str());
    }
}

// test multipole direct polarization for system comprised of two ammonia molecules; no cutoff

static void testMultipoleAmmoniaDirectPolarization() {

    std::string testName      = "testMultipoleAmmoniaDirectPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Direct, 
                                             cutoff, inputPmeGridDimension);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    getForcesEnergyMultipoleAmmonia(context, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -1.7428832e+01;

    expectedForces[0]         = Vec3( -3.5574000e+02, -7.3919340e+00,  3.8989934e+01);
    expectedForces[1]         = Vec3(  3.0368045e+01, -8.7325694e+00,  6.9731151e+00);
    expectedForces[2]         = Vec3(  3.2358980e+01,  1.0234924e+01,  4.7203694e-01);
    expectedForces[3]         = Vec3(  2.1439022e+01,  5.8998414e+00, -3.8355239e+01);
    expectedForces[4]         = Vec3( -1.8052760e+02, -1.0618455e+00, -7.0030146e+01);
    expectedForces[5]         = Vec3(  4.2411304e+01, -1.6569222e+01,  1.9047581e+00);
    expectedForces[6]         = Vec3(  3.6823677e+02,  7.7839986e-01,  5.8404590e+01);
    expectedForces[7]         = Vec3(  4.1453480e+01,  1.6842405e+01,  1.6409513e+00);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test multipole mutual polarization for system comprised of two ammonia molecules; no cutoff

static void testMultipoleAmmoniaMutualPolarization() {

    std::string testName      = "testMultipoleAmmoniaMutualPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    getForcesEnergyMultipoleAmmonia(context, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = -1.7790449e+01;

    expectedForces[0]         = Vec3( -3.7523158e+02,  -7.9806295e+00,   3.7464051e+01);
    expectedForces[1]         = Vec3(  3.1352410e+01,  -9.4055551e+00,   8.5230415e+00);
    expectedForces[2]         = Vec3(  3.3504923e+01,   1.1029935e+01,   1.5052263e+00);
    expectedForces[3]         = Vec3(  2.3295507e+01,   6.3698827e+00,  -4.0403553e+01);
    expectedForces[4]         = Vec3( -1.9379275e+02,  -1.0903937e+00,  -7.3461740e+01);
    expectedForces[5]         = Vec3(  4.3278067e+01,  -1.6906589e+01,   1.5721909e+00);
    expectedForces[6]         = Vec3(  3.9529983e+02,   7.9661172e-01,   6.3499055e+01);
    expectedForces[7]         = Vec3(  4.2293601e+01,   1.7186738e+01,   1.3017270e+00);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    
    // Try changing the particle parameters and make sure it's still correct.
    
    for (int i = 0; i < numberOfParticles; i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        std::vector<double> dipole, quadrupole;
        amoebaMultipoleForce->getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        dipole[0] *= 0.7;
        quadrupole[2] *= 1.5;
        quadrupole[6] *= 1.5;
        amoebaMultipoleForce->setMultipoleParameters(i, 1.1*charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, 1.3*thole, 1.4*damping, 1.5*polarity);
    }
    LangevinIntegrator integrator2(0.0, 0.1, 0.01);
    Context context2(system, integrator2, context.getPlatform());
    context2.setPositions(context.getState(State::Positions).getPositions());
    State state1 = context.getState(State::Forces | State::Energy);
    State state2 = context2.getState(State::Forces | State::Energy);
    bool exceptionThrown = false;
    try {
        // This should throw an exception.
        compareForcesEnergy(testName, state2.getPotentialEnergy(), state1.getPotentialEnergy(), state2.getForces(), state1.getForces(), tolerance);
        for (int i = 0; i < numberOfParticles; i++)
            ASSERT_EQUAL_VEC(state1.getForces()[i], state2.getForces()[i], tolerance);
    }
    catch (std::exception ex) {
        exceptionThrown = true;
    }
    ASSERT(exceptionThrown);
    amoebaMultipoleForce->updateParametersInContext(context);
    state1 = context.getState(State::Forces | State::Energy);
    compareForcesEnergy(testName, state2.getPotentialEnergy(), state1.getPotentialEnergy(), state2.getForces(), state1.getForces(), tolerance);
}

// setup for box of 4 water molecules -- used to test PME

static void setupAndGetForcesEnergyMultipoleWater(AmoebaMultipoleForce::NonbondedMethod nonbondedMethod,
                                                   AmoebaMultipoleForce::PolarizationType polarizationType,
                                                   double cutoff, int inputPmeGridDimension, std::vector<Vec3>& forces,
                                                   double& energy) {

    // beginning of Multipole setup

    System system;

    // box dimensions

    double boxDimension                               = 1.8643;
    Vec3 a(boxDimension, 0.0, 0.0);
    Vec3 b(0.0, boxDimension, 0.0);
    Vec3 c(0.0, 0.0, boxDimension);
    system.setDefaultPeriodicBoxVectors(a, b, c);

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 12;
    amoebaMultipoleForce->setNonbondedMethod(nonbondedMethod);
    amoebaMultipoleForce->setPolarizationType(polarizationType);
    amoebaMultipoleForce->setCutoffDistance(cutoff);
    amoebaMultipoleForce->setMutualInducedTargetEpsilon(1.0e-06);
    amoebaMultipoleForce->setMutualInducedMaxIterations(500);
    amoebaMultipoleForce->setAEwald(5.4459052e+00);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-04);

    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        system.addParticle(1.5995000e+01);
        system.addParticle(1.0080000e+00);
        system.addParticle(1.0080000e+00);
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   3.5403072e-04;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =  -3.9025708e-04;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   3.6226356e-05;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =  -3.4284825e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =  -1.0024088e-04;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   1.3452570e-04;

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaMultipoleForce->addMultipole(-5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
    }

    // CovalentMaps

    std::vector< int > covalentMap;
    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

        covalentMap.resize(0);
        covalentMap.push_back(jj);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
    } 
 
    // 1-2 bonds needed

    AmoebaBondForce* amoebaBondForce  = new AmoebaBondForce();

    // addBond: particle1, particle2, length, quadraticK

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaBondForce->addBond(jj, jj+1,   0.0000000e+00,   0.0000000e+00);
        amoebaBondForce->addBond(jj, jj+2,   0.0000000e+00,   0.0000000e+00);
    }

    amoebaBondForce->setAmoebaGlobalBondCubic(-2.5500000e+01); 
    amoebaBondForce->setAmoebaGlobalBondQuartic(3.7931250e+02); 
    system.addForce(amoebaBondForce);

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3( -8.7387270e-01,   5.3220410e-01,    7.4214000e-03);
    positions[1]              = Vec3( -9.6050090e-01,   5.1173410e-01,   -2.2202700e-02);
    positions[2]              = Vec3( -8.5985900e-01,   4.9658230e-01,    1.0283390e-01);
    positions[3]              = Vec3(  9.1767100e-02,  -7.8956650e-01,    4.3804200e-01);
    positions[4]              = Vec3(  1.2333420e-01,  -7.0267430e-01,    4.2611550e-01);
    positions[5]              = Vec3(  1.7267090e-01,  -8.2320810e-01,    4.8124750e-01);
    positions[6]              = Vec3(  8.6290110e-01,   6.2153500e-02,    4.1280850e-01);
    positions[7]              = Vec3(  8.6385200e-01,   1.2684730e-01,    3.3887060e-01);
    positions[8]              = Vec3(  9.5063550e-01,   5.3173300e-02,    4.4799160e-01);
    positions[9]              = Vec3(  5.0844930e-01,   2.8684740e-01,   -6.9293750e-01);
    positions[10]             = Vec3(  6.0459330e-01,   3.0620510e-01,   -7.0100130e-01);
    positions[11]             = Vec3(  5.0590640e-01,   1.8880920e-01,   -6.8813470e-01);

    system.addForce(amoebaMultipoleForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);
    State state                      = context.getState(State::Forces | State::Energy);
    forces                           = state.getForces();
    energy                           = state.getPotentialEnergy();
}

// test multipole direct polarization using PME for box of water

static void testMultipoleWaterPMEDirectPolarization() {

    std::string testName      = "testMultipoleWaterDirectPolarization";

    int numberOfParticles     = 12;
    int inputPmeGridDimension = 20;
    double cutoff             = 0.70;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Direct, 
                                            cutoff, inputPmeGridDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     = 6.4585115e-01;

    expectedForces[0]         = Vec3( -1.2396731e+00,  -2.4231698e+01,   8.3348523e+00);
    expectedForces[1]         = Vec3( -3.3737276e+00,   9.9304523e+00,  -6.3917827e+00);
    expectedForces[2]         = Vec3(  4.4062247e+00,   1.9518971e+01,  -4.6552873e+00);
    expectedForces[3]         = Vec3( -1.3128824e+00,  -1.2887339e+00,  -1.4473147e+00);
    expectedForces[4]         = Vec3(  2.1137034e+00,   3.9457973e-01,   2.9269129e-01);
    expectedForces[5]         = Vec3(  1.0271174e+00,   1.2039367e+00,   1.2112214e+00);
    expectedForces[6]         = Vec3( -3.2082903e+00,   1.4979371e+01,  -1.0274832e+00);
    expectedForces[7]         = Vec3( -1.1880320e+00,  -1.5177166e+01,   2.5525509e+00);
    expectedForces[8]         = Vec3(  4.3607105e+00,  -7.0253274e+00,   2.9522580e-01);
    expectedForces[9]         = Vec3( -3.0175134e+00,   1.3607102e+00,   6.6883370e+00);
    expectedForces[10]        = Vec3(  9.2036949e-01,  -1.4717629e+00,  -3.3362339e+00);
    expectedForces[11]        = Vec3(  1.2523841e+00,  -1.9794292e+00,  -3.4670129e+00);

    double tolerance          = 1.0e-03;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// test multipole mutual polarization using PME for box of water

static void testMultipoleWaterPMEMutualPolarization() {

    std::string testName      = "testMultipoleWaterMutualPolarization";

    int numberOfParticles     = 12;
    int inputPmeGridDimension = 20;
    double cutoff             = 0.70;
    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                            cutoff, inputPmeGridDimension, forces, energy);
    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  6.5029855e-01;

    expectedForces[0]         = Vec3( -1.2367386e+00,  -2.4197036e+01,   8.3256759e+00);
    expectedForces[1]         = Vec3( -3.3825187e+00,   9.9387618e+00,  -6.4200475e+00);
    expectedForces[2]         = Vec3(  4.4108644e+00,   1.9486127e+01,  -4.6530661e+00);
    expectedForces[3]         = Vec3( -1.3129168e+00,  -1.2947383e+00,  -1.4438198e+00);
    expectedForces[4]         = Vec3(  2.1144837e+00,   3.9590305e-01,   2.9040889e-01);
    expectedForces[5]         = Vec3(  1.0287222e+00,   1.2100201e+00,   1.2103068e+00);
    expectedForces[6]         = Vec3( -3.2017550e+00,   1.4995985e+01,  -1.1036504e+00);
    expectedForces[7]         = Vec3( -1.2065398e+00,  -1.5192899e+01,   2.6233368e+00);
    expectedForces[8]         = Vec3(  4.3698604e+00,  -7.0550315e+00,   3.4204565e-01);
    expectedForces[9]         = Vec3( -3.0082825e+00,   1.3575082e+00,   6.6901032e+00);
    expectedForces[10]        = Vec3(  9.1775539e-01,  -1.4651882e+00,  -3.3322516e+00);
    expectedForces[11]        = Vec3(  1.2467701e+00,  -1.9832979e+00,  -3.4684052e+00);

    double tolerance          = 1.0e-03;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
}

// check validation of traceless/symmetric quadrupole tensor

static void testQuadrupoleValidation() {

    std::string testName      = "checkQuadrupoleValidation";

    int numberOfParticles     = 12;
    int pmeGridDimension      = 20;
    double cutoff             = 0.70;

    // beginning of Multipole setup

    System system;

    double boxDimension                               = 1.8643;
    Vec3 a(boxDimension, 0.0, 0.0);
    Vec3 b(0.0, boxDimension, 0.0);
    Vec3 c(0.0, 0.0, boxDimension);
    system.setDefaultPeriodicBoxVectors(a, b, c);

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    std::vector<Vec3> expectedForces(numberOfParticles);
    amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::PME);
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Direct);
    amoebaMultipoleForce->setCutoffDistance(0.7);
    amoebaMultipoleForce->setMutualInducedTargetEpsilon(1.0e-06);
    amoebaMultipoleForce->setMutualInducedMaxIterations(500);
    amoebaMultipoleForce->setAEwald(5.4459052e+00);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-04);

    std::vector<int> pmeGridDimensions(3);
    pmeGridDimensions[0] = pmeGridDimensions[1] = pmeGridDimensions[2] = pmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimensions);

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        system.addParticle(1.5995000e+01);
        system.addParticle(1.0080000e+00);
        system.addParticle(1.0080000e+00);
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   3.5403072e-04;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =  -3.9025708e-04;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   3.6226356e-05;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =  -3.4284825e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =  -1.0024088e-04;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   1.3452570e-04;

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaMultipoleForce->addMultipole(-5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
    }

    // CovalentMaps
/*
    std::vector< int > covalentMap;
    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

        covalentMap.resize(0);
        covalentMap.push_back(jj);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
    } 
*/ 
    AmoebaBondForce* amoebaBondForce  = new AmoebaBondForce();

    // addBond: particle1, particle2, length, quadraticK

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaBondForce->addBond(jj, jj+1,   0.0000000e+00,   0.0000000e+00);
        amoebaBondForce->addBond(jj, jj+2,   0.0000000e+00,   0.0000000e+00);
    }

    amoebaBondForce->setAmoebaGlobalBondCubic(-2.5500000e+01); 
    amoebaBondForce->setAmoebaGlobalBondQuartic(3.7931250e+02); 
    system.addForce(amoebaBondForce);

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3( -8.7387270e-01,   5.3220410e-01,    7.4214000e-03);
    positions[1]              = Vec3( -9.6050090e-01,   5.1173410e-01,   -2.2202700e-02);
    positions[2]              = Vec3( -8.5985900e-01,   4.9658230e-01,    1.0283390e-01);
    positions[3]              = Vec3(  9.1767100e-02,  -7.8956650e-01,    4.3804200e-01);
    positions[4]              = Vec3(  1.2333420e-01,  -7.0267430e-01,    4.2611550e-01);
    positions[5]              = Vec3(  1.7267090e-01,  -8.2320810e-01,    4.8124750e-01);
    positions[6]              = Vec3(  8.6290110e-01,   6.2153500e-02,    4.1280850e-01);
    positions[7]              = Vec3(  8.6385200e-01,   1.2684730e-01,    3.3887060e-01);
    positions[8]              = Vec3(  9.5063550e-01,   5.3173300e-02,    4.4799160e-01);
    positions[9]              = Vec3(  5.0844930e-01,   2.8684740e-01,   -6.9293750e-01);
    positions[10]             = Vec3(  6.0459330e-01,   3.0620510e-01,   -7.0100130e-01);
    positions[11]             = Vec3(  5.0590640e-01,   1.8880920e-01,   -6.8813470e-01);

    system.addForce(amoebaMultipoleForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);

    // traceless quadrupole

    try {
        oxygenMolecularQuadrupole[4] += 0.1;
        amoebaMultipoleForce->setMultipoleParameters(0, -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, 1, 2, -1,
                                                       3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        State state                      = context.getState(State::Forces | State::Energy);
        std::stringstream buffer;        
        buffer << "Exception not thrown for quadrupole tensor w/ nonzero trace.";
        throw OpenMMException(buffer.str());
    } catch(const std::exception& e) {
    }
    oxygenMolecularQuadrupole[4] -= 0.1;

    // symmetric quadrupole

    // XY and YX components

    try {
        oxygenMolecularQuadrupole[1] += 0.1;
        amoebaMultipoleForce->setMultipoleParameters(0, -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, 1, 2, -1,
                                                       3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        State state                      = context.getState(State::Forces | State::Energy);
        std::stringstream buffer;        
        buffer << "Exception not thrown for quadrupole tensor w/ nonzero trace.";
        throw OpenMMException(buffer.str());
    } catch(const std::exception& e) {
    }
    oxygenMolecularQuadrupole[1] -= 0.1;

    // XZ and ZX components

    try {
        oxygenMolecularQuadrupole[2] += 0.1;
        amoebaMultipoleForce->setMultipoleParameters(0, -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, 1, 2, -1,
                                                       3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        State state                      = context.getState(State::Forces | State::Energy);
        std::stringstream buffer;        
        buffer << "Exception not thrown for quadrupole tensor w/ nonzero trace.";
        throw OpenMMException(buffer.str());
    } catch(const std::exception& e) {
    }
    oxygenMolecularQuadrupole[2] -= 0.1;

    // YZ and ZY components

    try {
        oxygenMolecularQuadrupole[5] += 0.1;
        amoebaMultipoleForce->setMultipoleParameters(0, -5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, 1, 2, -1,
                                                       3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        State state                      = context.getState(State::Forces | State::Energy);
        std::stringstream buffer;        
        buffer << "Exception not thrown for quadrupole tensor w/ nonzero trace.";
        throw OpenMMException(buffer.str());
    } catch(const std::exception& e) {
    }
    oxygenMolecularQuadrupole[5] -= 0.1;

}

// setup for box of 2 water molecules and 3 ions

// this method does too much; I tried passing the context ptr back to 
// the tests methods, but the tests would seg fault w/ a bad_alloc error

static void setupAndGetForcesEnergyMultipoleIonsAndWater(AmoebaMultipoleForce::NonbondedMethod nonbondedMethod,
                                                          AmoebaMultipoleForce::PolarizationType polarizationType,
                                                          double cutoff, int inputPmeGridDimension, std::string testName,
                                                          std::vector<Vec3>& forces, double& energy) {

    // beginning of Multipole setup

    System system;

    // box dimensions

    double boxDimensions[3]                           = { 6.7538, 7.2977, 7.4897 };
    Vec3 a(boxDimensions[0], 0.0, 0.0);
    Vec3 b(0.0, boxDimensions[1], 0.0);
    Vec3 c(0.0, 0.0, boxDimensions[2]);
    system.setDefaultPeriodicBoxVectors(a, b, c);


    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 8;
    int numberOfWaters                                = 2;
    int numberOfIons                                  = numberOfParticles - numberOfWaters*3;

    amoebaMultipoleForce->setNonbondedMethod(nonbondedMethod);
    amoebaMultipoleForce->setPolarizationType(polarizationType);
    amoebaMultipoleForce->setCutoffDistance(cutoff);
    amoebaMultipoleForce->setMutualInducedTargetEpsilon(1.0e-06);
    amoebaMultipoleForce->setMutualInducedMaxIterations(500);
    amoebaMultipoleForce->setAEwald(5.4459052e+00);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-05);

    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    // 2 ions

    system.addParticle(3.5453000e+01 );
    system.addParticle(2.2990000e+01);

    std::vector<double> ionDipole(3);
    std::vector<double> ionQuadrupole(9);

    ionDipole[0]     =   0.0000000e+00;
    ionDipole[1]     =   0.0000000e+00;
    ionDipole[2]     =   0.0000000e+00;

    ionQuadrupole[0] =   0.0000000e+00;
    ionQuadrupole[1] =   0.0000000e+00;
    ionQuadrupole[2] =   0.0000000e+00;
    ionQuadrupole[3] =   0.0000000e+00;
    ionQuadrupole[4] =   0.0000000e+00;
    ionQuadrupole[5] =   0.0000000e+00;
    ionQuadrupole[6] =   0.0000000e+00;
    ionQuadrupole[7] =   0.0000000e+00;
    ionQuadrupole[8] =   0.0000000e+00;
    amoebaMultipoleForce->addMultipole( -1.0000000e+00, ionDipole, ionQuadrupole, 5, -1, -1, -1,   3.9000000e-01,   3.9842202e-01,   4.0000000e-03);
    amoebaMultipoleForce->addMultipole(  1.0000000e+00, ionDipole, ionQuadrupole, 5, -1, -1, -1,   3.9000000e-01,   2.2209062e-01,   1.2000000e-04);

    // waters

    for (unsigned int jj = 2; jj < numberOfParticles; jj += 3) {
        system.addParticle(1.5999000e+01);
        system.addParticle(1.0080000e+00);
        system.addParticle(1.0080000e+00);
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   3.5403072e-04;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =  -3.9025708e-04;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   3.6226356e-05;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =  -3.4284825e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =  -1.0024088e-04;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   1.3452570e-04;

    for (unsigned int jj = 2; jj < numberOfParticles; jj += 3) {
        amoebaMultipoleForce->addMultipole(-5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
    }

    // CovalentMaps

    std::vector< int > covalentMap;
    covalentMap.resize(0);
    covalentMap.push_back(0);
    amoebaMultipoleForce->setCovalentMap(0, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    covalentMap.resize(0);
    covalentMap.push_back(1);
    amoebaMultipoleForce->setCovalentMap(1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);

    for (unsigned int jj = 2; jj < numberOfParticles; jj += 3) {
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

        covalentMap.resize(0);
        covalentMap.push_back(jj);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
    } 
 
    // 1-2 bonds needed

    AmoebaBondForce* amoebaBondForce  = new AmoebaBondForce();

    // addBond: particle1, particle2, length, quadraticK

    for (unsigned int jj = 2; jj < numberOfParticles; jj += 3) {
        amoebaBondForce->addBond(jj, jj+1,   0.0000000e+00,   0.0000000e+00);
        amoebaBondForce->addBond(jj, jj+2,   0.0000000e+00,   0.0000000e+00);
    }

    amoebaBondForce->setAmoebaGlobalBondCubic(-2.5500000e+01); 
    amoebaBondForce->setAmoebaGlobalBondQuartic(3.7931250e+02); 
    system.addForce(amoebaBondForce);

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]              = Vec3( -1.4364000e+00,  -1.2848000e+00,    5.1940000e-01);
    positions[1]              = Vec3( -3.2644000e+00,   2.3620000e+00,    1.3643000e+00);
    positions[2]              = Vec3( -2.3780000e+00,   1.8976000e+00,   -1.5921000e+00);
    positions[3]              = Vec3( -2.3485183e+00,   1.8296632e+00,   -1.5310146e+00);
    positions[4]              = Vec3( -2.3784362e+00,   1.8623910e+00,   -1.6814092e+00);
    positions[5]              = Vec3( -2.1821000e+00,  -1.0808000e+00,    2.9547000e+00);
    positions[6]              = Vec3( -2.1198155e+00,  -1.0925202e+00,    2.8825940e+00);
    positions[7]              = Vec3( -2.1537255e+00,  -1.0076218e+00,    3.0099797e+00);

    system.addForce(amoebaMultipoleForce);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context = Context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);

        State state               = context.getState(State::Forces | State::Energy);
        forces                    = state.getForces();
        energy                    = state.getPotentialEnergy();

    return;

}

// test multipole mutual polarization using PME for system comprised of 2 ions and 2 waters

static void testMultipoleIonsAndWaterPMEDirectPolarization() {

    std::string testName      = "testMultipoleIonsAndWaterDirectPolarization";

    int numberOfParticles     = 8;
    int inputPmeGridDimension = 64;
    double cutoff             = 0.70;

    std::vector<Vec3> forces;
    double energy;

    setupAndGetForcesEnergyMultipoleIonsAndWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Direct, 
                                                  cutoff, inputPmeGridDimension, testName, forces, energy);

    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy     =  -4.6859568e+01;

    expectedForces[0]         = Vec3( -9.1266563e+00,   1.5193632e+01,  -4.0047974e+00);
    expectedForces[1]         = Vec3( -1.0497973e+00,   1.4622548e+01,   1.1789324e+01);
    expectedForces[2]         = Vec3( -3.2564644e+00,   6.5325105e+00,  -2.9698616e+00);
    expectedForces[3]         = Vec3(  3.0687040e+00,  -8.4253665e-01,  -3.4081010e+00);
    expectedForces[4]         = Vec3(  1.1407201e+00,  -3.1491550e+00,  -1.1326031e+00);
    expectedForces[5]         = Vec3( -6.1046529e+00,   9.5686061e-01,   1.1506333e-01);
    expectedForces[6]         = Vec3(  1.9275403e+00,  -5.6007439e-01,  -4.8387346e+00);
    expectedForces[7]         = Vec3(  4.0644209e+00,  -3.3666305e+00,  -1.7022384e+00);

    double tolerance          = 5.0e-04;
    compareForceNormsEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);

}

// test multipole mutual polarization using PME for system comprised of 2 ions and 2 waters

static void testMultipoleIonsAndWaterPMEMutualPolarization() {

    std::string testName            = "testMultipoleIonsAndWaterMutualPolarization";

    int numberOfParticles           = 8;
    int inputPmeGridDimension       = 64;
    double cutoff                   = 0.70;

    std::vector<Vec3> forces;
    double energy;

    std::vector<Vec3> inputGrid;

    setupAndGetForcesEnergyMultipoleIonsAndWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                                  cutoff, inputPmeGridDimension, testName, forces, energy);

    std::vector<Vec3> expectedForces(numberOfParticles);

    double expectedEnergy           = -4.6859424e+01;

    expectedForces[0]               = Vec3( -9.1272358e+00,   1.5191516e+01,  -4.0058826e+00);
    expectedForces[1]               = Vec3( -1.0497156e+00,   1.4622425e+01,   1.1789420e+01);
    expectedForces[2]               = Vec3( -3.2560478e+00,   6.5289712e+00,  -2.9779483e+00);
    expectedForces[3]               = Vec3(  3.0672153e+00,  -8.4407797e-01,  -3.4094884e+00);
    expectedForces[4]               = Vec3(  1.1382586e+00,  -3.1512949e+00,  -1.1387028e+00);
    expectedForces[5]               = Vec3( -6.1050295e+00,   9.5345692e-01,   1.1488832e-01);
    expectedForces[6]               = Vec3(  1.9319945e+00,  -5.5747599e-01,  -4.8469044e+00);
    expectedForces[7]               = Vec3(  4.0622614e+00,  -3.3687594e+00,  -1.6986575e+00);

    //double tolerance                = 1.0e-03;
    //compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    double tolerance                = 5.0e-04;
    compareForceNormsEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);

}

// setup for box of 216 water molecules -- used to test PME

static void setupAndGetForcesEnergyMultipoleLargeWater(AmoebaMultipoleForce::NonbondedMethod nonbondedMethod,
                                                        AmoebaMultipoleForce::PolarizationType polarizationType,
                                                        double cutoff, int inputPmeGridDimension, std::string& testName, 
                                                        std::vector<Vec3>& forces, double& energy,
                                                        std::vector< double >& outputMultipoleMoments,
                                                        std::vector< Vec3 >& inputGrid,
                                                        std::vector< double >& outputGridPotential) {

    // beginning of Multipole setup

    System system;

    // box dimensions

    double boxDimension                               = 1.8643;
    Vec3 a(boxDimension, 0.0, 0.0);
    Vec3 b(0.0, boxDimension, 0.0);
    Vec3 c(0.0, 0.0, boxDimension);
    system.setDefaultPeriodicBoxVectors(a, b, c);

    AmoebaMultipoleForce* amoebaMultipoleForce        = new AmoebaMultipoleForce();;
    int numberOfParticles                             = 648;
    amoebaMultipoleForce->setNonbondedMethod(nonbondedMethod);
    amoebaMultipoleForce->setPolarizationType(polarizationType);
    amoebaMultipoleForce->setCutoffDistance(cutoff);
    amoebaMultipoleForce->setMutualInducedTargetEpsilon(1.0e-06);
    amoebaMultipoleForce->setMutualInducedMaxIterations(500);
    amoebaMultipoleForce->setAEwald(5.4459052e+00);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-04);

    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        system.addParticle(1.5995000e+01);
        system.addParticle(1.0080000e+00);
        system.addParticle(1.0080000e+00);
    }

    std::vector<double> oxygenMolecularDipole(3);
    std::vector<double> oxygenMolecularQuadrupole(9);

    oxygenMolecularDipole[0]     =   0.0000000e+00;
    oxygenMolecularDipole[1]     =   0.0000000e+00;
    oxygenMolecularDipole[2]     =   7.5561214e-03;

    oxygenMolecularQuadrupole[0] =   3.5403072e-04;
    oxygenMolecularQuadrupole[1] =   0.0000000e+00;
    oxygenMolecularQuadrupole[2] =   0.0000000e+00;
    oxygenMolecularQuadrupole[3] =   0.0000000e+00;
    oxygenMolecularQuadrupole[4] =  -3.9025708e-04;
    oxygenMolecularQuadrupole[5] =   0.0000000e+00;
    oxygenMolecularQuadrupole[6] =   0.0000000e+00;
    oxygenMolecularQuadrupole[7] =   0.0000000e+00;
    oxygenMolecularQuadrupole[8] =   3.6226356e-05;

    std::vector<double> hydrogenMolecularDipole(3);
    std::vector<double> hydrogenMolecularQuadrupole(9);
    hydrogenMolecularDipole[0]     =  -2.0420949e-03;
    hydrogenMolecularDipole[1]     =   0.0000000e+00;
    hydrogenMolecularDipole[2]     =  -3.0787530e-03;

    hydrogenMolecularQuadrupole[0] =  -3.4284825e-05;
    hydrogenMolecularQuadrupole[1] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[2] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[3] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[4] =  -1.0024088e-04;
    hydrogenMolecularQuadrupole[5] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[6] =  -1.8948597e-06;
    hydrogenMolecularQuadrupole[7] =   0.0000000e+00;
    hydrogenMolecularQuadrupole[8] =   1.3452570e-04;

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaMultipoleForce->addMultipole(-5.1966000e-01, oxygenMolecularDipole, oxygenMolecularQuadrupole, 1, jj+1, jj+2, -1,
                                            3.9000000e-01, 3.0698765e-01, 8.3700000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+2, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
        amoebaMultipoleForce->addMultipole( 2.5983000e-01, hydrogenMolecularDipole, hydrogenMolecularQuadrupole, 0, jj, jj+1, -1, 
                                            3.9000000e-01, 2.8135002e-01, 4.9600000e-04);
    }

    // CovalentMaps

    std::vector< int > covalentMap;
    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);

        covalentMap.resize(0);
        covalentMap.push_back(jj);
        covalentMap.push_back(jj+1);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj,   static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(4), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(0), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+2);
        amoebaMultipoleForce->setCovalentMap(jj+1, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
        covalentMap.resize(0);
        covalentMap.push_back(jj+1);
        amoebaMultipoleForce->setCovalentMap(jj+2, static_cast<OpenMM::AmoebaMultipoleForce::CovalentType>(1), covalentMap);
    
    } 
    system.addForce(amoebaMultipoleForce);
 
    // 1-2 bonds needed

    AmoebaBondForce* amoebaBondForce  = new AmoebaBondForce();

    // addBond: particle1, particle2, length, quadraticK

    for (unsigned int jj = 0; jj < numberOfParticles; jj += 3) {
        amoebaBondForce->addBond(jj, jj+1,   0.0000000e+00,   0.0000000e+00);
        amoebaBondForce->addBond(jj, jj+2,   0.0000000e+00,   0.0000000e+00);
    }

    amoebaBondForce->setAmoebaGlobalBondCubic(0.0); 
    amoebaBondForce->setAmoebaGlobalBondQuartic(0.0); 
    system.addForce(amoebaBondForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3(  8.0394300e-01,   5.8680350e-01,    4.9277700e-02);
    positions[1]             = Vec3(  7.5814940e-01,   5.0226660e-01,    4.0375900e-02);
    positions[2]             = Vec3(  8.2870560e-01,   6.0624400e-01,   -3.9707400e-02);
    positions[3]             = Vec3(  1.1484000e-02,  -8.8765990e-01,    6.4458520e-01);
    positions[4]             = Vec3(  9.5892500e-02,  -8.4464940e-01,    6.4052470e-01);
    positions[5]             = Vec3(  2.4723500e-02,  -9.7944710e-01,    6.1378930e-01);
    positions[6]             = Vec3( -6.5763670e-01,  -2.5260000e-02,    8.1046320e-01);
    positions[7]             = Vec3( -6.6454990e-01,   6.8992500e-02,    7.8963560e-01);
    positions[8]             = Vec3( -6.6845370e-01,  -4.0076000e-02,    9.1037470e-01);
    positions[9]             = Vec3(  6.5831270e-01,   8.5501500e-02,   -6.6685290e-01);
    positions[10]            = Vec3(  6.2600580e-01,   8.8732600e-02,   -5.7651320e-01);
    positions[11]            = Vec3(  6.1694860e-01,   5.3229000e-03,   -7.0543230e-01);
    positions[12]            = Vec3(  5.4954790e-01,   6.4357640e-01,    1.8420070e-01);
    positions[13]            = Vec3(  4.7740750e-01,   6.5609280e-01,    1.2079650e-01);
    positions[14]            = Vec3(  6.2544340e-01,   6.3485600e-01,    1.2346110e-01);
    positions[15]            = Vec3( -4.6646340e-01,  -8.5021310e-01,   -2.6526210e-01);
    positions[16]            = Vec3( -4.5053590e-01,  -8.3883300e-01,   -3.6069710e-01);
    positions[17]            = Vec3( -5.5653260e-01,  -8.7510810e-01,   -2.5955820e-01);
    positions[18]            = Vec3( -7.7550740e-01,  -4.6613180e-01,    4.9045930e-01);
    positions[19]            = Vec3( -7.3577510e-01,  -5.4400590e-01,    5.3107060e-01);
    positions[20]            = Vec3( -7.0755520e-01,  -4.1773140e-01,    4.4037930e-01);
    positions[21]            = Vec3( -2.8190600e-02,   7.4872450e-01,   -7.6855300e-01);
    positions[22]            = Vec3( -7.9443300e-02,   7.4463600e-01,   -6.8256160e-01);
    positions[23]            = Vec3(  1.7033100e-02,   8.3813000e-01,   -7.6365310e-01);
    positions[24]            = Vec3( -3.7112750e-01,  -2.2624390e-01,   -1.9170030e-01);
    positions[25]            = Vec3( -4.4236150e-01,  -2.4258640e-01,   -2.4723220e-01);
    positions[26]            = Vec3( -4.0233380e-01,  -2.2106530e-01,   -9.7227800e-02);
    positions[27]            = Vec3( -5.8120030e-01,  -5.6157220e-01,    8.3549400e-02);
    positions[28]            = Vec3( -6.6764500e-01,  -5.7119710e-01,    1.2970660e-01);
    positions[29]            = Vec3( -5.1434340e-01,  -5.5317060e-01,    1.5597670e-01);
    positions[30]            = Vec3(  8.5281410e-01,   4.9997870e-01,    3.4439320e-01);
    positions[31]            = Vec3(  8.8661040e-01,   4.7595500e-01,    4.3409810e-01);
    positions[32]            = Vec3(  7.6829200e-01,   4.5403270e-01,    3.3783460e-01);
    positions[33]            = Vec3(  6.2913000e-03,   3.9622090e-01,   -6.4448110e-01);
    positions[34]            = Vec3( -6.4546800e-02,   4.4539620e-01,   -6.0008300e-01);
    positions[35]            = Vec3(  7.0262000e-03,   3.1229330e-01,   -5.9892730e-01);
    positions[36]            = Vec3(  1.6883500e-02,   6.5824910e-01,    6.0982750e-01);
    positions[37]            = Vec3(  2.9114400e-02,   6.3714540e-01,    7.0403040e-01);
    positions[38]            = Vec3( -3.9569500e-02,   5.9419720e-01,    5.6714930e-01);
    positions[39]            = Vec3(  3.7393550e-01,   6.2909200e-01,    8.1318410e-01);
    positions[40]            = Vec3(  4.1500630e-01,   6.1010560e-01,    9.0110400e-01);
    positions[41]            = Vec3(  4.3953600e-01,   5.9208230e-01,    7.5268270e-01);
    positions[42]            = Vec3(  3.2500410e-01,   4.5615770e-01,   -2.5643980e-01);
    positions[43]            = Vec3(  3.7432790e-01,   4.5313140e-01,   -3.3754880e-01);
    positions[44]            = Vec3(  2.6987370e-01,   5.3785040e-01,   -2.4860760e-01);
    positions[45]            = Vec3(  5.6184630e-01,   5.2015900e-01,    6.3763990e-01);
    positions[46]            = Vec3(  5.6189080e-01,   5.6140190e-01,    5.5312940e-01);
    positions[47]            = Vec3(  5.4901540e-01,   4.2688810e-01,    6.2109450e-01);
    positions[48]            = Vec3( -8.7750980e-01,   6.9408570e-01,   -6.1784650e-01);
    positions[49]            = Vec3( -8.2179580e-01,   7.3187880e-01,   -5.4705510e-01);
    positions[50]            = Vec3( -9.0362240e-01,   7.7367480e-01,   -6.6488210e-01);
    positions[51]            = Vec3( -6.9406820e-01,   2.2491740e-01,    7.1940890e-01);
    positions[52]            = Vec3( -7.3674620e-01,   2.2091000e-01,    6.3486690e-01);
    positions[53]            = Vec3( -7.4149900e-01,   2.8970280e-01,    7.7200060e-01);
    positions[54]            = Vec3(  4.8285280e-01,  -1.8445100e-02,    3.1521130e-01);
    positions[55]            = Vec3(  5.5574910e-01,   2.4338500e-02,    2.7236750e-01);
    positions[56]            = Vec3(  4.1347360e-01,   5.0063500e-02,    3.2371450e-01);
    positions[57]            = Vec3( -2.2024800e-01,  -3.1071870e-01,    9.1706370e-01);
    positions[58]            = Vec3( -2.3195790e-01,  -4.0722320e-01,    9.2465160e-01);
    positions[59]            = Vec3( -2.8015290e-01,  -2.9349640e-01,    8.4209880e-01);
    positions[60]            = Vec3(  1.6893780e-01,   6.6734280e-01,   -2.4352040e-01);
    positions[61]            = Vec3(  1.9716270e-01,   7.5186390e-01,   -2.0536790e-01);
    positions[62]            = Vec3(  8.7430700e-02,   6.4225300e-01,   -1.9539020e-01);
    positions[63]            = Vec3( -9.0804840e-01,  -6.2437310e-01,   -8.8188300e-02);
    positions[64]            = Vec3( -8.6732940e-01,  -7.0428590e-01,   -4.8030200e-02);
    positions[65]            = Vec3( -8.3644480e-01,  -5.8139450e-01,   -1.3828190e-01);
    positions[66]            = Vec3( -8.6567760e-01,  -8.6537570e-01,    5.6295900e-02);
    positions[67]            = Vec3( -8.1778220e-01,  -9.4654890e-01,    8.4163600e-02);
    positions[68]            = Vec3( -9.4534460e-01,  -8.6858770e-01,    1.0560810e-01);
    positions[69]            = Vec3( -5.7716930e-01,  -2.6316670e-01,   -4.5880740e-01);
    positions[70]            = Vec3( -5.4569620e-01,  -3.1693230e-01,   -5.2720970e-01);
    positions[71]            = Vec3( -5.5496000e-01,  -1.7071220e-01,   -4.7392400e-01);
    positions[72]            = Vec3(  7.2367810e-01,  -8.4678300e-01,   -6.9502250e-01);
    positions[73]            = Vec3(  7.9899670e-01,  -8.9648580e-01,   -7.2759260e-01);
    positions[74]            = Vec3(  7.5075030e-01,  -8.1725850e-01,   -6.0600380e-01);
    positions[75]            = Vec3( -2.3769060e-01,  -6.2523350e-01,    1.2921080e-01);
    positions[76]            = Vec3( -1.8309420e-01,  -6.2163180e-01,    4.8693900e-02);
    positions[77]            = Vec3( -2.3929030e-01,  -5.3708810e-01,    1.6453540e-01);
    positions[78]            = Vec3(  8.3347800e-02,  -5.0189060e-01,    5.4317800e-01);
    positions[79]            = Vec3(  1.0917180e-01,  -5.7641330e-01,    4.8632230e-01);
    positions[80]            = Vec3(  1.4837200e-02,  -5.5084220e-01,    5.9546910e-01);
    positions[81]            = Vec3(  7.4250070e-01,  -2.7418580e-01,    8.3795900e-02);
    positions[82]            = Vec3(  6.8666720e-01,  -2.4554090e-01,    1.6206940e-01);
    positions[83]            = Vec3(  7.1516850e-01,  -3.6419530e-01,    7.2493400e-02);
    positions[84]            = Vec3( -2.5059100e-02,   8.6314620e-01,    2.2861410e-01);
    positions[85]            = Vec3(  9.6445000e-03,   9.0720400e-01,    1.4964290e-01);
    positions[86]            = Vec3(  4.5097900e-02,   8.7155360e-01,    2.9051950e-01);
    positions[87]            = Vec3(  4.7779490e-01,   9.0242640e-01,    8.2515620e-01);
    positions[88]            = Vec3(  4.3957480e-01,   8.0786830e-01,    8.2489220e-01);
    positions[89]            = Vec3(  4.6833310e-01,   9.2867710e-01,    9.1788160e-01);
    positions[90]            = Vec3(  8.2204140e-01,   9.0145630e-01,   -2.5081510e-01);
    positions[91]            = Vec3(  8.5191840e-01,   8.1397830e-01,   -2.2168590e-01);
    positions[92]            = Vec3(  7.6397810e-01,   9.2011290e-01,   -1.8137750e-01);
    positions[93]            = Vec3( -7.9443650e-01,   1.7601300e-01,    4.6436790e-01);
    positions[94]            = Vec3( -7.9212150e-01,   2.3533020e-01,    3.8657500e-01);
    positions[95]            = Vec3( -8.7057070e-01,   1.1288830e-01,    4.4595260e-01);
    positions[96]            = Vec3(  3.2425690e-01,   3.8214720e-01,   -8.2471120e-01);
    positions[97]            = Vec3(  2.8321830e-01,   4.2912450e-01,   -7.4875880e-01);
    positions[98]            = Vec3(  2.7681870e-01,   2.9837230e-01,   -8.2620080e-01);
    positions[99]            = Vec3(  7.5575820e-01,  -8.9620900e-01,    2.3680670e-01);
    positions[100]           = Vec3(  6.6600420e-01,  -8.7027760e-01,    2.7104280e-01);
    positions[101]           = Vec3(  8.1544110e-01,  -9.1190240e-01,    3.1149610e-01);
    positions[102]           = Vec3( -8.4248740e-01,   3.5007110e-01,   -4.4389740e-01);
    positions[103]           = Vec3( -7.5693800e-01,   3.9510690e-01,   -4.4710480e-01);
    positions[104]           = Vec3( -8.6984880e-01,   3.5457480e-01,   -5.3702920e-01);
    positions[105]           = Vec3(  3.8837250e-01,  -4.8496240e-01,    6.5322550e-01);
    positions[106]           = Vec3(  4.1237110e-01,  -4.0401080e-01,    7.0255980e-01);
    positions[107]           = Vec3(  3.0065040e-01,  -4.6399160e-01,    6.0513040e-01);
    positions[108]           = Vec3(  6.2063930e-01,  -5.0831230e-01,    4.9540430e-01);
    positions[109]           = Vec3(  6.8959700e-01,  -5.3506820e-01,    5.6328860e-01);
    positions[110]           = Vec3(  5.3663630e-01,  -5.1121830e-01,    5.4640900e-01);
    positions[111]           = Vec3(  7.0354670e-01,  -5.1748580e-01,   -7.3878700e-02);
    positions[112]           = Vec3(  7.8529450e-01,  -5.6535940e-01,   -9.5943500e-02);
    positions[113]           = Vec3(  6.7807440e-01,  -4.7921810e-01,   -1.6187590e-01);
    positions[114]           = Vec3( -4.4116790e-01,  -4.7749880e-01,    3.0876830e-01);
    positions[115]           = Vec3( -5.0645290e-01,  -4.1075220e-01,    3.1159470e-01);
    positions[116]           = Vec3( -4.6594720e-01,  -5.2568230e-01,    3.8755370e-01);
    positions[117]           = Vec3( -9.1937480e-01,  -5.8400000e-05,   -2.5359570e-01);
    positions[118]           = Vec3( -8.5894750e-01,  -7.0402500e-02,   -2.2230370e-01);
    positions[119]           = Vec3( -8.7441760e-01,   8.3170500e-02,   -2.3447490e-01);
    positions[120]           = Vec3(  5.0867290e-01,   2.3568780e-01,    5.5935510e-01);
    positions[121]           = Vec3(  4.1446460e-01,   2.6088930e-01,    5.8683440e-01);
    positions[122]           = Vec3(  5.1853820e-01,   1.4937830e-01,    5.8561390e-01);
    positions[123]           = Vec3( -4.6831090e-01,  -6.1465890e-01,   -1.6794620e-01);
    positions[124]           = Vec3( -4.8688540e-01,  -5.9611250e-01,   -7.4636500e-02);
    positions[125]           = Vec3( -4.9162010e-01,  -7.0497770e-01,   -1.8127910e-01);
    positions[126]           = Vec3( -3.1791800e-01,  -5.4450000e-03,   -3.6397680e-01);
    positions[127]           = Vec3( -2.2253910e-01,  -2.4457600e-02,   -3.5240990e-01);
    positions[128]           = Vec3( -3.6044390e-01,  -3.5065000e-02,   -2.8414310e-01);
    positions[129]           = Vec3(  1.0461140e-01,   2.6758700e-01,   -2.2684050e-01);
    positions[130]           = Vec3(  1.8426490e-01,   3.2453330e-01,   -2.3574350e-01);
    positions[131]           = Vec3(  1.0569370e-01,   2.3628020e-01,   -1.3834830e-01);
    positions[132]           = Vec3( -1.4119340e-01,   4.1653970e-01,   -2.7320250e-01);
    positions[133]           = Vec3( -5.2065100e-02,   3.6979030e-01,   -2.6662970e-01);
    positions[134]           = Vec3( -1.3834110e-01,   4.7690560e-01,   -1.9435870e-01);
    positions[135]           = Vec3( -7.6602450e-01,  -2.1216400e-01,   -1.9516640e-01);
    positions[136]           = Vec3( -8.0191290e-01,  -2.8391260e-01,   -1.3557910e-01);
    positions[137]           = Vec3( -7.4415500e-01,  -2.6044280e-01,   -2.8169590e-01);
    positions[138]           = Vec3( -1.3600310e-01,   1.9674000e-01,    2.0349610e-01);
    positions[139]           = Vec3( -1.6201050e-01,   2.8693750e-01,    2.3123820e-01);
    positions[140]           = Vec3( -2.1785650e-01,   1.4514420e-01,    1.9201990e-01);
    positions[141]           = Vec3(  6.2897820e-01,  -4.2302590e-01,   -7.6557210e-01);
    positions[142]           = Vec3(  6.2334100e-01,  -4.4471660e-01,   -6.7174140e-01);
    positions[143]           = Vec3(  6.3346670e-01,  -5.0696850e-01,   -8.0495300e-01);
    positions[144]           = Vec3(  9.1588260e-01,  -3.9845200e-02,    3.5189180e-01);
    positions[145]           = Vec3(  9.8891550e-01,  -4.4673900e-02,    2.9156120e-01);
    positions[146]           = Vec3(  8.4126090e-01,  -2.2841000e-03,    2.9707980e-01);
    positions[147]           = Vec3(  4.8470900e-01,  -8.2561400e-02,    6.0082980e-01);
    positions[148]           = Vec3(  3.9021850e-01,  -6.2932500e-02,    6.0195610e-01);
    positions[149]           = Vec3(  5.0563070e-01,  -7.9866200e-02,    5.0777230e-01);
    positions[150]           = Vec3( -7.2845180e-01,  -3.4650580e-01,    7.5973620e-01);
    positions[151]           = Vec3( -7.6073760e-01,  -3.6974690e-01,    6.7323450e-01);
    positions[152]           = Vec3( -7.1326740e-01,  -2.4916760e-01,    7.5651020e-01);
    positions[153]           = Vec3( -3.0896820e-01,  -3.8029640e-01,    6.5520670e-01);
    positions[154]           = Vec3( -3.5019560e-01,  -4.5571260e-01,    6.1040330e-01);
    positions[155]           = Vec3( -2.8479430e-01,  -3.2175460e-01,    5.7933340e-01);
    positions[156]           = Vec3( -6.2826700e-02,  -6.4315900e-02,   -6.8812300e-02);
    positions[157]           = Vec3( -7.4971500e-02,   1.9900000e-02,   -1.8191100e-02);
    positions[158]           = Vec3(  3.2478400e-02,  -8.8932300e-02,   -5.6413600e-02);
    positions[159]           = Vec3(  1.1667520e-01,  -6.6784990e-01,    1.1452860e-01);
    positions[160]           = Vec3(  6.5194200e-02,  -7.3080350e-01,    6.5294000e-02);
    positions[161]           = Vec3(  1.6133150e-01,  -6.1778770e-01,    4.7196600e-02);
    positions[162]           = Vec3(  8.8627400e-02,  -7.1850240e-01,    3.7581390e-01);
    positions[163]           = Vec3(  1.2356120e-01,  -8.0690930e-01,    3.7094210e-01);
    positions[164]           = Vec3(  9.2028600e-02,  -6.8313750e-01,    2.8412340e-01);
    positions[165]           = Vec3(  2.1347270e-01,   8.4107000e-03,    6.0413030e-01);
    positions[166]           = Vec3(  1.8845570e-01,   4.3251500e-02,    5.1600410e-01);
    positions[167]           = Vec3(  1.6789670e-01,  -7.6656800e-02,    6.0793520e-01);
    positions[168]           = Vec3(  1.8425700e-02,   3.0164400e-02,    8.4213210e-01);
    positions[169]           = Vec3( -7.1641800e-02,   4.1848500e-02,    8.7065260e-01);
    positions[170]           = Vec3(  4.4510400e-02,  -6.2982500e-02,    8.6373290e-01);
    positions[171]           = Vec3( -3.1486750e-01,  -1.9966860e-01,   -5.7954700e-01);
    positions[172]           = Vec3( -3.2321140e-01,  -1.4613590e-01,   -5.0133480e-01);
    positions[173]           = Vec3( -3.4769180e-01,  -1.4810900e-01,   -6.5567720e-01);
    positions[174]           = Vec3(  2.2013690e-01,  -4.8207100e-02,   -6.6169910e-01);
    positions[175]           = Vec3(  1.3676160e-01,  -9.4600100e-02,   -6.4525960e-01);
    positions[176]           = Vec3(  2.7051720e-01,  -1.2158460e-01,   -6.9535940e-01);
    positions[177]           = Vec3( -1.5721060e-01,  -2.0015580e-01,    4.8442010e-01);
    positions[178]           = Vec3( -7.4675400e-02,  -2.0952300e-01,    5.3560160e-01);
    positions[179]           = Vec3( -1.8522760e-01,  -1.0781560e-01,    5.0024110e-01);
    positions[180]           = Vec3(  5.4002730e-01,   6.3800500e-01,   -8.0040500e-01);
    positions[181]           = Vec3(  5.0366070e-01,   7.1545920e-01,   -7.5257350e-01);
    positions[182]           = Vec3(  5.1480770e-01,   5.5941670e-01,   -7.4903220e-01);
    positions[183]           = Vec3( -6.3383580e-01,   5.7282910e-01,   -1.7429980e-01);
    positions[184]           = Vec3( -6.0668100e-01,   4.7712900e-01,   -1.7677570e-01);
    positions[185]           = Vec3( -5.6638740e-01,   6.1288510e-01,   -2.2951390e-01);
    positions[186]           = Vec3( -2.0998170e-01,  -2.7747820e-01,    7.0579400e-02);
    positions[187]           = Vec3( -1.4055440e-01,  -3.0201380e-01,    1.3644740e-01);
    positions[188]           = Vec3( -1.6881700e-01,  -2.1818660e-01,    6.9733000e-03);
    positions[189]           = Vec3( -7.6400000e-04,   5.6326380e-01,    1.4175360e-01);
    positions[190]           = Vec3( -7.3688000e-02,   5.0031150e-01,    1.5514670e-01);
    positions[191]           = Vec3( -2.5553000e-02,   6.4733770e-01,    1.7711800e-01);
    positions[192]           = Vec3(  3.9595890e-01,  -1.9078420e-01,   -1.9708050e-01);
    positions[193]           = Vec3(  4.3887020e-01,  -1.5694200e-01,   -1.1582060e-01);
    positions[194]           = Vec3(  3.7635540e-01,  -1.1834040e-01,   -2.5323660e-01);
    positions[195]           = Vec3(  3.9638900e-02,  -2.4093090e-01,    8.9424300e-01);
    positions[196]           = Vec3( -4.9643600e-02,  -2.7156660e-01,    8.9962920e-01);
    positions[197]           = Vec3(  8.4318200e-02,  -2.7149840e-01,    9.7721820e-01);
    positions[198]           = Vec3( -5.9039370e-01,  -3.5975630e-01,   -7.1984370e-01);
    positions[199]           = Vec3( -5.3914870e-01,  -4.0214860e-01,   -7.8361060e-01);
    positions[200]           = Vec3( -6.8562580e-01,  -3.9051900e-01,   -7.4071320e-01);
    positions[201]           = Vec3(  4.7759800e-01,   3.2863960e-01,   -5.4274200e-02);
    positions[202]           = Vec3(  4.5034450e-01,   3.6680450e-01,   -1.4201230e-01);
    positions[203]           = Vec3(  4.3083410e-01,   3.8043410e-01,    1.5118500e-02);
    positions[204]           = Vec3(  1.8100450e-01,   1.6674000e-01,   -8.4907090e-01);
    positions[205]           = Vec3(  1.0479500e-01,   1.5721720e-01,   -9.0737790e-01);
    positions[206]           = Vec3(  1.7365410e-01,   9.7140100e-02,   -7.7842430e-01);
    positions[207]           = Vec3( -6.9841710e-01,   8.5211760e-01,    4.9956020e-01);
    positions[208]           = Vec3( -6.3194850e-01,   9.0336360e-01,    4.5467020e-01);
    positions[209]           = Vec3( -6.7863830e-01,   7.5666570e-01,    5.1268950e-01);
    positions[210]           = Vec3(  8.0356880e-01,  -7.6669620e-01,    5.6240980e-01);
    positions[211]           = Vec3(  8.9444390e-01,  -7.9421520e-01,    5.4379860e-01);
    positions[212]           = Vec3(  8.0061200e-01,  -7.1151420e-01,    6.3743510e-01);
    positions[213]           = Vec3( -2.3686380e-01,   4.4018650e-01,    2.7494630e-01);
    positions[214]           = Vec3( -2.1006750e-01,   4.1932880e-01,    3.6593160e-01);
    positions[215]           = Vec3( -3.2910900e-01,   4.6299420e-01,    2.7725190e-01);
    positions[216]           = Vec3(  7.3324180e-01,   9.1021100e-02,    8.6347740e-01);
    positions[217]           = Vec3(  6.4934460e-01,   5.3444800e-02,    8.7843600e-01);
    positions[218]           = Vec3(  7.1407590e-01,   1.8691830e-01,    8.6323690e-01);
    positions[219]           = Vec3(  3.6906600e-02,   1.4742360e-01,    4.0082880e-01);
    positions[220]           = Vec3( -1.0515300e-02,   1.4450010e-01,    4.8531790e-01);
    positions[221]           = Vec3( -3.6861400e-02,   1.5333190e-01,    3.3364650e-01);
    positions[222]           = Vec3(  5.7666790e-01,  -9.2075640e-01,    5.7305300e-01);
    positions[223]           = Vec3(  5.4452540e-01,  -9.3954290e-01,    6.5798160e-01);
    positions[224]           = Vec3(  6.7020160e-01,  -8.8052280e-01,    5.6852240e-01);
    positions[225]           = Vec3(  4.1616300e-01,  -2.3723450e-01,    7.8105700e-02);
    positions[226]           = Vec3(  4.4947640e-01,  -2.2465620e-01,    1.6469280e-01);
    positions[227]           = Vec3(  3.6093380e-01,  -3.1332780e-01,    7.1125100e-02);
    positions[228]           = Vec3( -1.9830990e-01,  -6.8678560e-01,   -7.6648560e-01);
    positions[229]           = Vec3( -1.1489950e-01,  -6.8356660e-01,   -8.2028210e-01);
    positions[230]           = Vec3( -2.0935090e-01,  -5.9618710e-01,   -7.3178710e-01);
    positions[231]           = Vec3( -4.3741650e-01,  -7.8889500e-01,    1.7785560e-01);
    positions[232]           = Vec3( -3.6424030e-01,  -7.2995610e-01,    1.5380490e-01);
    positions[233]           = Vec3( -5.0710310e-01,  -7.4066850e-01,    1.3917790e-01);
    positions[234]           = Vec3(  5.1605280e-01,   6.8521860e-01,    4.5545030e-01);
    positions[235]           = Vec3(  5.3920960e-01,   7.6750670e-01,    4.8965960e-01);
    positions[236]           = Vec3(  5.4441350e-01,   6.8153880e-01,    3.6305340e-01);
    positions[237]           = Vec3( -9.1377180e-01,   9.0412110e-01,   -8.0577110e-01);
    positions[238]           = Vec3( -8.6299150e-01,   9.8552780e-01,   -7.9463610e-01);
    positions[239]           = Vec3( -9.1270510e-01,   8.7715830e-01,   -9.0107170e-01);
    positions[240]           = Vec3( -5.6874630e-01,  -3.9330600e-02,    5.3540210e-01);
    positions[241]           = Vec3( -6.0667690e-01,   3.6619200e-02,    4.9922460e-01);
    positions[242]           = Vec3( -5.8307630e-01,  -4.4694300e-02,    6.3380260e-01);
    positions[243]           = Vec3(  1.0312020e-01,   2.2809180e-01,    5.7525600e-02);
    positions[244]           = Vec3(  1.8161800e-02,   2.2164820e-01,    1.0293620e-01);
    positions[245]           = Vec3(  1.4691520e-01,   3.0734480e-01,    9.4432600e-02);
    positions[246]           = Vec3( -5.3437690e-01,  -9.0689060e-01,   -7.7012560e-01);
    positions[247]           = Vec3( -6.0761130e-01,  -8.5593580e-01,   -8.0463440e-01);
    positions[248]           = Vec3( -5.5313680e-01,  -9.9745020e-01,   -8.0224750e-01);
    positions[249]           = Vec3(  1.7436730e-01,  -4.6935620e-01,   -7.7408150e-01);
    positions[250]           = Vec3(  1.3315640e-01,  -4.6856170e-01,   -6.8363440e-01);
    positions[251]           = Vec3(  2.3486700e-01,  -3.9970620e-01,   -7.7872930e-01);
    positions[252]           = Vec3(  5.0382310e-01,   8.6391330e-01,   -6.1751380e-01);
    positions[253]           = Vec3(  5.7851670e-01,   9.1774780e-01,   -6.3741940e-01);
    positions[254]           = Vec3(  5.2100060e-01,   8.2278060e-01,   -5.3449130e-01);
    positions[255]           = Vec3( -2.3461000e-03,   8.8439120e-01,   -3.5703750e-01);
    positions[256]           = Vec3( -4.5869800e-02,   9.2025060e-01,   -4.4264850e-01);
    positions[257]           = Vec3(  7.7568300e-02,   8.3812640e-01,   -3.7824790e-01);
    positions[258]           = Vec3( -1.6677150e-01,  -9.0353490e-01,   -5.6323410e-01);
    positions[259]           = Vec3( -1.5077930e-01,  -8.7448310e-01,   -6.5150250e-01);
    positions[260]           = Vec3( -2.5054260e-01,  -8.5746520e-01,   -5.4471400e-01);
    positions[261]           = Vec3( -1.0245710e-01,  -4.1390500e-01,    2.9240710e-01);
    positions[262]           = Vec3( -1.6375100e-01,  -3.5806090e-01,    3.3803800e-01);
    positions[263]           = Vec3( -3.4371600e-02,  -4.4188880e-01,    3.6032470e-01);
    positions[264]           = Vec3(  6.7721230e-01,  -9.2755000e-01,   -6.1695000e-03);
    positions[265]           = Vec3(  6.3209610e-01,  -8.4066740e-01,    1.5854000e-03);
    positions[266]           = Vec3(  7.2195780e-01,  -9.3506790e-01,    7.6821700e-02);
    positions[267]           = Vec3( -5.2597410e-01,   5.0741940e-01,    2.8142130e-01);
    positions[268]           = Vec3( -5.3172740e-01,   5.6506650e-01,    2.0013640e-01);
    positions[269]           = Vec3( -5.9533220e-01,   4.4193270e-01,    2.6673520e-01);
    positions[270]           = Vec3(  4.3852700e-02,  -7.1092730e-01,   -3.0056810e-01);
    positions[271]           = Vec3(  1.2232900e-02,  -6.7601300e-01,   -2.1679320e-01);
    positions[272]           = Vec3(  3.0039200e-02,  -8.0474130e-01,   -3.0050550e-01);
    positions[273]           = Vec3(  4.7537430e-01,   6.7956000e-03,   -8.8926760e-01);
    positions[274]           = Vec3(  4.4972180e-01,  -8.1937800e-02,   -8.5037740e-01);
    positions[275]           = Vec3(  3.9238110e-01,   5.4650000e-02,   -9.0978500e-01);
    positions[276]           = Vec3(  8.4526190e-01,  -3.2384610e-01,    4.4702430e-01);
    positions[277]           = Vec3(  8.5335920e-01,  -2.3860050e-01,    4.1507690e-01);
    positions[278]           = Vec3(  9.3799800e-01,  -3.6222940e-01,    4.5249690e-01);
    positions[279]           = Vec3( -8.5624140e-01,  -3.3540460e-01,   -5.3955060e-01);
    positions[280]           = Vec3( -8.9833150e-01,  -3.2177130e-01,   -6.2636700e-01);
    positions[281]           = Vec3( -7.6568080e-01,  -3.0076830e-01,   -5.3672910e-01);
    positions[282]           = Vec3( -4.0866080e-01,  -7.0070860e-01,    9.2586930e-01);
    positions[283]           = Vec3( -4.5043520e-01,  -7.7640050e-01,    9.7012510e-01);
    positions[284]           = Vec3( -3.2086210e-01,  -6.9414110e-01,    9.6526100e-01);
    positions[285]           = Vec3( -2.9612090e-01,   2.9021400e-01,   -4.6137730e-01);
    positions[286]           = Vec3( -3.0085180e-01,   1.9752840e-01,   -4.3159520e-01);
    positions[287]           = Vec3( -2.4502340e-01,   3.3756140e-01,   -3.9070450e-01);
    positions[288]           = Vec3( -8.4956240e-01,  -3.3051010e-01,    4.2215900e-02);
    positions[289]           = Vec3( -8.2077940e-01,  -3.9086690e-01,    1.1548590e-01);
    positions[290]           = Vec3( -9.3822180e-01,  -3.1618550e-01,    5.2894000e-02);
    positions[291]           = Vec3( -8.6464030e-01,   7.5345250e-01,    1.9545370e-01);
    positions[292]           = Vec3( -9.2073720e-01,   6.7584430e-01,    1.8998460e-01);
    positions[293]           = Vec3( -8.9310500e-01,   7.8515510e-01,    2.8077440e-01);
    positions[294]           = Vec3(  5.3248170e-01,   6.8435100e-02,   -1.1431070e-01);
    positions[295]           = Vec3(  6.1630600e-01,   5.7417300e-02,   -6.9794300e-02);
    positions[296]           = Vec3(  4.9275030e-01,   1.5234490e-01,   -7.9235100e-02);
    positions[297]           = Vec3( -3.0166400e-02,   3.6028840e-01,   -9.2023940e-01);
    positions[298]           = Vec3(  2.5390700e-02,   4.3355180e-01,   -9.4581010e-01);
    positions[299]           = Vec3( -1.2837900e-02,   3.5198820e-01,   -8.2331230e-01);
    positions[300]           = Vec3( -7.6094250e-01,  -7.4142570e-01,   -7.6415170e-01);
    positions[301]           = Vec3( -7.5826150e-01,  -6.8315050e-01,   -8.4024930e-01);
    positions[302]           = Vec3( -7.8169550e-01,  -6.8557300e-01,   -6.8728990e-01);
    positions[303]           = Vec3( -7.1618050e-01,  -8.6617600e-02,   -7.8297100e-01);
    positions[304]           = Vec3( -6.9164460e-01,  -1.6643810e-01,   -7.3660090e-01);
    positions[305]           = Vec3( -8.1169890e-01,  -8.6541300e-02,   -7.7380060e-01);
    positions[306]           = Vec3(  8.6280550e-01,  -2.8731190e-01,   -7.5013210e-01);
    positions[307]           = Vec3(  8.4297110e-01,  -2.0142080e-01,   -7.9688520e-01);
    positions[308]           = Vec3(  7.7553640e-01,  -3.3421630e-01,   -7.5754400e-01);
    positions[309]           = Vec3(  2.9607200e-02,  -6.7251560e-01,   -9.1368960e-01);
    positions[310]           = Vec3(  1.0909000e-03,  -6.2708430e-01,   -9.9528360e-01);
    positions[311]           = Vec3(  8.0161300e-02,  -5.9814710e-01,   -8.7106130e-01);
    positions[312]           = Vec3( -2.2829370e-01,   4.6661410e-01,    7.7985190e-01);
    positions[313]           = Vec3( -2.4730820e-01,   5.6404020e-01,    7.8763210e-01);
    positions[314]           = Vec3( -1.7899690e-01,   4.3324110e-01,    8.5622400e-01);
    positions[315]           = Vec3( -5.1323270e-01,  -2.6480150e-01,    7.2113100e-02);
    positions[316]           = Vec3( -4.4180310e-01,  -2.8480730e-01,    1.4166490e-01);
    positions[317]           = Vec3( -5.5826690e-01,  -3.4508980e-01,    5.6782300e-02);
    positions[318]           = Vec3(  3.5970320e-01,  -7.1101700e-01,   -8.5706800e-01);
    positions[319]           = Vec3(  3.5573750e-01,  -6.6123030e-01,   -7.7069560e-01);
    positions[320]           = Vec3(  2.9308100e-01,  -6.7738800e-01,   -9.1162920e-01);
    positions[321]           = Vec3( -8.6077820e-01,  -8.3187420e-01,    3.5264550e-01);
    positions[322]           = Vec3( -7.9919290e-01,  -8.9965630e-01,    3.8875110e-01);
    positions[323]           = Vec3( -8.4377450e-01,  -8.2428940e-01,    2.5657630e-01);
    positions[324]           = Vec3( -6.9407750e-01,   8.5240530e-01,   -4.8975260e-01);
    positions[325]           = Vec3( -6.0369970e-01,   8.2005830e-01,   -4.7948010e-01);
    positions[326]           = Vec3( -6.8257340e-01,   9.0158170e-01,   -5.7057020e-01);
    positions[327]           = Vec3( -8.6181560e-01,   2.1174420e-01,    3.2775000e-02);
    positions[328]           = Vec3( -9.5070390e-01,   2.5868190e-01,    2.6787700e-02);
    positions[329]           = Vec3( -8.8015990e-01,   1.3696510e-01,    9.1486900e-02);
    positions[330]           = Vec3( -6.7034530e-01,  -7.0959980e-01,    5.7197940e-01);
    positions[331]           = Vec3( -6.3447070e-01,  -7.7970770e-01,    5.1435410e-01);
    positions[332]           = Vec3( -7.1147280e-01,  -7.6230200e-01,    6.4084900e-01);
    positions[333]           = Vec3( -4.2433970e-01,   1.6353470e-01,   -7.5364040e-01);
    positions[334]           = Vec3( -3.3715920e-01,   1.3734360e-01,   -7.8660110e-01);
    positions[335]           = Vec3( -4.5203330e-01,   2.3873860e-01,   -8.1607320e-01);
    positions[336]           = Vec3( -4.2091960e-01,  -8.1633330e-01,   -5.3063920e-01);
    positions[337]           = Vec3( -4.2728590e-01,  -7.1806470e-01,   -5.4109270e-01);
    positions[338]           = Vec3( -4.5013260e-01,  -8.3810340e-01,   -6.1998700e-01);
    positions[339]           = Vec3(  6.0367930e-01,   3.3084920e-01,   -8.4465460e-01);
    positions[340]           = Vec3(  5.0455880e-01,   3.3698360e-01,   -8.4011240e-01);
    positions[341]           = Vec3(  6.2487550e-01,   2.4834360e-01,   -8.0607210e-01);
    positions[342]           = Vec3(  1.8546120e-01,  -6.3282200e-02,    5.1304500e-02);
    positions[343]           = Vec3(  2.8101390e-01,  -7.7771500e-02,    5.1163200e-02);
    positions[344]           = Vec3(  1.7127760e-01,   2.0996700e-02,    9.0574100e-02);
    positions[345]           = Vec3( -3.5029200e-02,  -7.9917400e-02,   -3.4468400e-01);
    positions[346]           = Vec3( -6.3903800e-02,  -6.0213300e-02,   -2.5206780e-01);
    positions[347]           = Vec3( -4.6785200e-02,  -1.7349570e-01,   -3.5772680e-01);
    positions[348]           = Vec3(  2.5567190e-01,   6.2355480e-01,    4.2852620e-01);
    positions[349]           = Vec3(  1.9093710e-01,   6.4505930e-01,    4.9102940e-01);
    positions[350]           = Vec3(  3.4540670e-01,   6.4937420e-01,    4.5902510e-01);
    positions[351]           = Vec3( -7.3742490e-01,  -8.7628820e-01,   -2.6411710e-01);
    positions[352]           = Vec3( -7.3220480e-01,  -9.1540050e-01,   -3.5104230e-01);
    positions[353]           = Vec3( -7.9968040e-01,  -9.2863850e-01,   -2.1682500e-01);
    positions[354]           = Vec3(  5.1017210e-01,  -2.7173980e-01,    7.9174500e-01);
    positions[355]           = Vec3(  5.1045830e-01,  -2.0746280e-01,    7.2138780e-01);
    positions[356]           = Vec3(  5.9967910e-01,  -3.0815350e-01,    7.9296320e-01);
    positions[357]           = Vec3(  6.1703300e-02,  -6.0490320e-01,   -5.4304490e-01);
    positions[358]           = Vec3(  6.5202000e-03,  -6.6388800e-01,   -5.9525970e-01);
    positions[359]           = Vec3(  6.2525700e-02,  -6.3466150e-01,   -4.5175130e-01);
    positions[360]           = Vec3( -5.0181950e-01,   6.8138390e-01,   -8.8794760e-01);
    positions[361]           = Vec3( -4.0469720e-01,   6.5541180e-01,   -8.8475300e-01);
    positions[362]           = Vec3( -5.4953810e-01,   6.3245150e-01,   -8.1669610e-01);
    positions[363]           = Vec3( -3.5708340e-01,   8.1787480e-01,    1.0372050e-01);
    positions[364]           = Vec3( -4.3575160e-01,   7.6657380e-01,    8.8357500e-02);
    positions[365]           = Vec3( -3.8126100e-01,   9.1312250e-01,    1.2894930e-01);
    positions[366]           = Vec3( -1.0889180e-01,   6.4289110e-01,   -1.1000150e-01);
    positions[367]           = Vec3( -9.5792300e-02,   6.5121590e-01,   -1.2915400e-02);
    positions[368]           = Vec3( -1.4253020e-01,   7.3532640e-01,   -1.2649680e-01);
    positions[369]           = Vec3( -8.0675190e-01,   3.8993580e-01,   -9.3061890e-01);
    positions[370]           = Vec3( -8.4285770e-01,   4.7693320e-01,   -9.5868770e-01);
    positions[371]           = Vec3( -7.4065520e-01,   4.1059110e-01,   -8.6270860e-01);
    positions[372]           = Vec3( -7.3221050e-01,  -8.3486000e-02,    1.8651540e-01);
    positions[373]           = Vec3( -6.6332990e-01,  -2.5838100e-02,    1.5155080e-01);
    positions[374]           = Vec3( -7.5939010e-01,  -1.4675440e-01,    1.1813700e-01);
    positions[375]           = Vec3(  6.1370510e-01,  -3.7510720e-01,   -2.9444790e-01);
    positions[376]           = Vec3(  5.3141590e-01,  -3.1971250e-01,   -2.8369080e-01);
    positions[377]           = Vec3(  6.7472620e-01,  -3.0544670e-01,   -3.2680390e-01);
    positions[378]           = Vec3(  2.8333090e-01,   7.0116700e-01,    6.3582400e-02);
    positions[379]           = Vec3(  2.3304950e-01,   7.8436370e-01,    8.8113000e-02);
    positions[380]           = Vec3(  2.1603670e-01,   6.3345680e-01,    4.3706900e-02);
    positions[381]           = Vec3(  3.4046290e-01,  -5.8425160e-01,   -5.8383960e-01);
    positions[382]           = Vec3(  4.2396660e-01,  -5.6867730e-01,   -5.4787780e-01);
    positions[383]           = Vec3(  2.7987870e-01,  -5.6273080e-01,   -5.1485370e-01);
    positions[384]           = Vec3(  4.8651200e-01,   3.9384650e-01,   -5.0852640e-01);
    positions[385]           = Vec3(  4.8954070e-01,   2.9830160e-01,   -5.1540010e-01);
    positions[386]           = Vec3(  5.7513360e-01,   4.2777280e-01,   -4.8094980e-01);
    positions[387]           = Vec3( -4.9931530e-01,  -8.6556710e-01,    4.1410020e-01);
    positions[388]           = Vec3( -4.0971070e-01,  -9.0364250e-01,    4.1539320e-01);
    positions[389]           = Vec3( -5.0187830e-01,  -8.1863570e-01,    3.2854240e-01);
    positions[390]           = Vec3( -9.2923250e-01,  -9.5140200e-02,    7.7175180e-01);
    positions[391]           = Vec3( -1.0068535e+00,  -4.9193300e-02,    8.1361050e-01);
    positions[392]           = Vec3( -8.5382270e-01,  -3.5167000e-02,    7.7988780e-01);
    positions[393]           = Vec3(  5.8200510e-01,  -2.7347380e-01,    3.2175080e-01);
    positions[394]           = Vec3(  5.9114530e-01,  -2.1232990e-01,    3.9188270e-01);
    positions[395]           = Vec3(  6.2697690e-01,  -3.5436570e-01,    3.5518080e-01);
    positions[396]           = Vec3( -4.3869270e-01,   7.1030180e-01,   -3.4435510e-01);
    positions[397]           = Vec3( -3.5798370e-01,   6.6801330e-01,   -3.8293170e-01);
    positions[398]           = Vec3( -3.9584820e-01,   7.8582280e-01,   -3.0015890e-01);
    positions[399]           = Vec3(  3.0315060e-01,   2.0553140e-01,    3.3518590e-01);
    positions[400]           = Vec3(  2.0466680e-01,   2.0029920e-01,    3.3800050e-01);
    positions[401]           = Vec3(  3.1784090e-01,   2.6138240e-01,    4.0966770e-01);
    positions[402]           = Vec3(  7.3144120e-01,   1.1861840e-01,    2.1872590e-01);
    positions[403]           = Vec3(  6.9245610e-01,   2.0755440e-01,    2.3848660e-01);
    positions[404]           = Vec3(  7.4250960e-01,   1.1063670e-01,    1.1673060e-01);
    positions[405]           = Vec3(  3.0774670e-01,  -6.7782260e-01,   -6.9330000e-02);
    positions[406]           = Vec3(  3.0161020e-01,  -7.5652530e-01,   -1.2627210e-01);
    positions[407]           = Vec3(  3.7612340e-01,  -6.9199170e-01,   -2.0688000e-03);
    positions[408]           = Vec3(  4.8241200e-02,   1.4991530e-01,   -4.8562930e-01);
    positions[409]           = Vec3(  7.0825700e-02,   1.7883510e-01,   -4.0076820e-01);
    positions[410]           = Vec3( -1.4581300e-02,   7.7868400e-02,   -4.8044320e-01);
    positions[411]           = Vec3(  2.6566210e-01,  -4.7972300e-02,   -3.9240060e-01);
    positions[412]           = Vec3(  2.5708940e-01,  -2.6958700e-02,   -4.8906580e-01);
    positions[413]           = Vec3(  1.8079360e-01,  -1.7099600e-02,   -3.5945650e-01);
    positions[414]           = Vec3(  7.3593670e-01,   3.2192010e-01,    6.3185000e-03);
    positions[415]           = Vec3(  7.5313070e-01,   3.1236830e-01,   -9.0780600e-02);
    positions[416]           = Vec3(  6.4125230e-01,   3.3242850e-01,    5.3072000e-03);
    positions[417]           = Vec3( -7.2074000e-03,  -2.1935180e-01,   -6.7044710e-01);
    positions[418]           = Vec3( -7.9916200e-02,  -2.2604130e-01,   -7.3330810e-01);
    positions[419]           = Vec3( -3.2871000e-03,  -2.9557560e-01,   -6.1702790e-01);
    positions[420]           = Vec3(  8.0182800e-01,   3.3340310e-01,   -2.5836160e-01);
    positions[421]           = Vec3(  8.9266890e-01,   3.1760310e-01,   -2.9990300e-01);
    positions[422]           = Vec3(  7.7135080e-01,   4.0881250e-01,   -3.1490320e-01);
    positions[423]           = Vec3( -3.1753700e-01,   3.7248900e-02,    5.0846140e-01);
    positions[424]           = Vec3( -3.3276340e-01,   1.2794660e-01,    5.4135580e-01);
    positions[425]           = Vec3( -4.0442920e-01,  -2.1535000e-03,    5.2164500e-01);
    positions[426]           = Vec3(  7.7089090e-01,  -1.7749490e-01,   -4.1090550e-01);
    positions[427]           = Vec3(  8.0919970e-01,  -9.9267700e-02,   -3.6080690e-01);
    positions[428]           = Vec3(  8.4794900e-01,  -2.2265030e-01,   -4.4286640e-01);
    positions[429]           = Vec3( -5.0985980e-01,   6.5271910e-01,    5.1660950e-01);
    positions[430]           = Vec3( -4.1891080e-01,   6.9500010e-01,    5.0933000e-01);
    positions[431]           = Vec3( -5.2072650e-01,   6.0609800e-01,    4.2889530e-01);
    positions[432]           = Vec3(  8.8931480e-01,  -1.5854900e-02,   -7.9057690e-01);
    positions[433]           = Vec3(  8.4049130e-01,   2.2454500e-02,   -7.1223150e-01);
    positions[434]           = Vec3(  8.6392620e-01,   4.6002000e-02,   -8.5696830e-01);
    positions[435]           = Vec3( -4.2632820e-01,  -5.4538160e-01,   -5.2698140e-01);
    positions[436]           = Vec3( -3.4047810e-01,  -5.2088280e-01,   -5.5637760e-01);
    positions[437]           = Vec3( -4.9107950e-01,  -5.2513960e-01,   -5.9520410e-01);
    positions[438]           = Vec3(  8.8830700e-01,   7.8506050e-01,    4.7420010e-01);
    positions[439]           = Vec3(  9.6737760e-01,   8.0796480e-01,    5.2210120e-01);
    positions[440]           = Vec3(  8.3449840e-01,   7.2694370e-01,    5.2968560e-01);
    positions[441]           = Vec3( -3.0889500e-02,  -5.4040860e-01,   -7.7446500e-02);
    positions[442]           = Vec3(  2.4910200e-02,  -4.7046460e-01,   -5.3187100e-02);
    positions[443]           = Vec3( -1.0937030e-01,  -5.1212170e-01,   -1.2642620e-01);
    positions[444]           = Vec3(  5.0722190e-01,  -8.0898340e-01,    3.3208510e-01);
    positions[445]           = Vec3(  5.1254280e-01,  -8.4333670e-01,    4.2962250e-01);
    positions[446]           = Vec3(  4.8459280e-01,  -7.1548850e-01,    3.3664280e-01);
    positions[447]           = Vec3(  7.0974400e-02,  -8.6268490e-01,   -7.2122900e-01);
    positions[448]           = Vec3(  8.8211100e-02,  -8.1266230e-01,   -7.9698760e-01);
    positions[449]           = Vec3(  1.4856180e-01,  -8.7440360e-01,   -6.6601020e-01);
    positions[450]           = Vec3( -2.7264270e-01,   8.2117820e-01,    4.0979220e-01);
    positions[451]           = Vec3( -1.8893860e-01,   7.8611730e-01,    4.4435560e-01);
    positions[452]           = Vec3( -2.7256440e-01,   8.1557060e-01,    3.0746650e-01);
    positions[453]           = Vec3( -2.3667600e-01,   7.0807760e-01,    9.0055470e-01);
    positions[454]           = Vec3( -1.7087350e-01,   7.0278860e-01,    9.7330650e-01);
    positions[455]           = Vec3( -2.2325560e-01,   8.0596230e-01,    8.7050690e-01);
    positions[456]           = Vec3(  6.0904540e-01,  -5.3471490e-01,   -5.1588800e-01);
    positions[457]           = Vec3(  6.6627390e-01,  -6.1177680e-01,   -4.9309950e-01);
    positions[458]           = Vec3(  6.1303950e-01,  -4.7414890e-01,   -4.3691960e-01);
    positions[459]           = Vec3( -6.9432470e-01,   5.5588670e-01,   -7.2750070e-01);
    positions[460]           = Vec3( -6.8524660e-01,   5.1427650e-01,   -6.4407660e-01);
    positions[461]           = Vec3( -7.7219850e-01,   6.0882800e-01,   -7.1352640e-01);
    positions[462]           = Vec3( -6.5544400e-01,   5.6801890e-01,    7.6654940e-01);
    positions[463]           = Vec3( -5.9853210e-01,   5.8150060e-01,    6.8630620e-01);
    positions[464]           = Vec3( -6.0728400e-01,   6.2604000e-01,    8.2970960e-01);
    positions[465]           = Vec3( -1.7725100e-01,  -7.5128040e-01,    4.8288320e-01);
    positions[466]           = Vec3( -1.1106490e-01,  -7.1604590e-01,    4.2681180e-01);
    positions[467]           = Vec3( -1.2808000e-01,  -8.2063050e-01,    5.2385060e-01);
    positions[468]           = Vec3(  5.0880810e-01,  -1.7782370e-01,   -5.5526690e-01);
    positions[469]           = Vec3(  4.7579150e-01,  -1.6757400e-01,   -4.6732050e-01);
    positions[470]           = Vec3(  6.0010540e-01,  -1.6566020e-01,   -5.4639700e-01);
    positions[471]           = Vec3(  7.9737120e-01,  -5.3326000e-03,   -2.2789800e-02);
    positions[472]           = Vec3(  7.5436910e-01,  -9.2537600e-02,   -2.7176000e-03);
    positions[473]           = Vec3(  8.4035540e-01,  -2.5845500e-02,   -1.0913300e-01);
    positions[474]           = Vec3(  2.4805290e-01,  -4.5182680e-01,   -2.5649240e-01);
    positions[475]           = Vec3(  2.6536400e-01,  -5.1313010e-01,   -1.8699050e-01);
    positions[476]           = Vec3(  2.8661880e-01,  -3.6531040e-01,   -2.2184290e-01);
    positions[477]           = Vec3(  8.9407190e-01,   6.4140150e-01,   -2.2838520e-01);
    positions[478]           = Vec3(  8.6394270e-01,   5.7649930e-01,   -2.9124340e-01);
    positions[479]           = Vec3(  9.8698980e-01,   6.3685520e-01,   -2.2087390e-01);
    positions[480]           = Vec3( -5.0297400e-01,   3.8595440e-01,   -9.1329410e-01);
    positions[481]           = Vec3( -4.2761710e-01,   4.4573350e-01,   -8.9961960e-01);
    positions[482]           = Vec3( -5.7109730e-01,   4.3620760e-01,   -9.6075240e-01);
    positions[483]           = Vec3( -5.7912630e-01,   3.1473530e-01,   -1.3174480e-01);
    positions[484]           = Vec3( -6.4359930e-01,   2.3775370e-01,   -1.3815260e-01);
    positions[485]           = Vec3( -5.1910350e-01,   2.9841740e-01,   -5.0386200e-02);
    positions[486]           = Vec3( -2.3287450e-01,  -4.5325250e-01,   -2.6295780e-01);
    positions[487]           = Vec3( -3.1705790e-01,  -5.0582880e-01,   -2.5755610e-01);
    positions[488]           = Vec3( -2.4988940e-01,  -3.6717760e-01,   -2.2456790e-01);
    positions[489]           = Vec3(  7.3902040e-01,   6.0596960e-01,    8.7531410e-01);
    positions[490]           = Vec3(  6.8510920e-01,   5.6584400e-01,    8.0394270e-01);
    positions[491]           = Vec3(  6.7885220e-01,   6.2492120e-01,    9.4854880e-01);
    positions[492]           = Vec3( -1.7342650e-01,  -4.4833620e-01,   -6.2689720e-01);
    positions[493]           = Vec3( -1.2120170e-01,  -4.7044370e-01,   -5.5178260e-01);
    positions[494]           = Vec3( -2.1756220e-01,  -3.7095040e-01,   -5.9046920e-01);
    positions[495]           = Vec3( -9.7909800e-02,   4.1047140e-01,    5.5154950e-01);
    positions[496]           = Vec3( -1.5085520e-01,   4.3078300e-01,    6.2923170e-01);
    positions[497]           = Vec3( -2.2121000e-02,   3.5854830e-01,    5.8628920e-01);
    positions[498]           = Vec3( -2.9249090e-01,   4.2502460e-01,   -7.0552100e-01);
    positions[499]           = Vec3( -2.3858950e-01,   3.8279980e-01,   -7.7129020e-01);
    positions[500]           = Vec3( -2.8537830e-01,   3.6194940e-01,   -6.3036240e-01);
    positions[501]           = Vec3( -4.1927200e-01,  -1.0765570e-01,   -8.1010100e-01);
    positions[502]           = Vec3( -4.5513170e-01,  -1.8389200e-01,   -8.5055730e-01);
    positions[503]           = Vec3( -4.6978720e-01,  -3.2915400e-02,   -8.4249770e-01);
    positions[504]           = Vec3( -8.3022800e-01,  -5.9366610e-01,   -5.2440890e-01);
    positions[505]           = Vec3( -8.3569020e-01,  -5.0053960e-01,   -5.4596070e-01);
    positions[506]           = Vec3( -7.7653500e-01,  -5.9680800e-01,   -4.4872510e-01);
    positions[507]           = Vec3(  4.7451900e-02,   2.4985900e-01,    7.1027380e-01);
    positions[508]           = Vec3(  5.2750000e-03,   2.6682820e-01,    8.0047760e-01);
    positions[509]           = Vec3(  9.2790500e-02,   1.6390540e-01,    7.2751450e-01);
    positions[510]           = Vec3(  9.8318300e-02,  -2.4834430e-01,    6.2217110e-01);
    positions[511]           = Vec3(  7.1376800e-02,  -2.3868900e-01,    7.1029050e-01);
    positions[512]           = Vec3(  1.0725160e-01,  -3.3946690e-01,    5.9525570e-01);
    positions[513]           = Vec3( -1.7389390e-01,   6.3857050e-01,   -4.3802350e-01);
    positions[514]           = Vec3( -1.0857550e-01,   7.0876020e-01,   -4.2023360e-01);
    positions[515]           = Vec3( -1.6180390e-01,   5.6775180e-01,   -3.7084920e-01);
    positions[516]           = Vec3( -8.3384410e-01,  -7.8320210e-01,    7.9714340e-01);
    positions[517]           = Vec3( -8.6597850e-01,  -8.7176550e-01,    7.7689410e-01);
    positions[518]           = Vec3( -9.1332720e-01,  -7.1912210e-01,    7.9807020e-01);
    positions[519]           = Vec3(  3.0122650e-01,   4.4099240e-01,    1.7747380e-01);
    positions[520]           = Vec3(  3.0879580e-01,   3.5962220e-01,    2.2668340e-01);
    positions[521]           = Vec3(  2.9198270e-01,   5.0655710e-01,    2.4655760e-01);
    positions[522]           = Vec3(  3.8346200e-01,  -2.8443150e-01,   -8.3961770e-01);
    positions[523]           = Vec3(  4.1227770e-01,  -2.9408340e-01,   -9.3409110e-01);
    positions[524]           = Vec3(  4.5498420e-01,  -3.3552520e-01,   -7.8643110e-01);
    positions[525]           = Vec3(  5.4535540e-01,   1.2249720e-01,   -4.0869350e-01);
    positions[526]           = Vec3(  6.0755050e-01,   1.6343320e-01,   -3.4805580e-01);
    positions[527]           = Vec3(  4.8362230e-01,   8.8573600e-02,   -3.4405000e-01);
    positions[528]           = Vec3(  1.3637990e-01,  -3.3186850e-01,    1.0338270e-01);
    positions[529]           = Vec3(  1.5761460e-01,  -2.5187340e-01,    1.5683210e-01);
    positions[530]           = Vec3(  7.8556700e-02,  -3.8461200e-01,    1.6118390e-01);
    positions[531]           = Vec3(  8.4245020e-01,   3.8084570e-01,   -6.9184990e-01);
    positions[532]           = Vec3(  9.0750590e-01,   3.9283710e-01,   -7.7288830e-01);
    positions[533]           = Vec3(  7.5053500e-01,   3.8878480e-01,   -7.2751780e-01);
    positions[534]           = Vec3(  2.7768360e-01,  -8.5899240e-01,   -5.3138620e-01);
    positions[535]           = Vec3(  2.8386750e-01,  -7.7018020e-01,   -5.6323660e-01);
    positions[536]           = Vec3(  3.4891330e-01,  -9.1242960e-01,   -5.6853820e-01);
    positions[537]           = Vec3(  2.6823810e-01,  -7.8504070e-01,    6.9926380e-01);
    positions[538]           = Vec3(  3.3824260e-01,  -8.3764610e-01,    7.3839250e-01);
    positions[539]           = Vec3(  3.0089590e-01,  -6.9098950e-01,    7.0290360e-01);
    positions[540]           = Vec3(  9.5946000e-02,   5.9757730e-01,    8.8417370e-01);
    positions[541]           = Vec3(  1.9084960e-01,   5.8892180e-01,    8.6811780e-01);
    positions[542]           = Vec3(  7.0090900e-02,   6.3001980e-01,    9.7622150e-01);
    positions[543]           = Vec3( -3.2687830e-01,  -9.5478000e-03,    2.1684540e-01);
    positions[544]           = Vec3( -3.2605730e-01,  -1.4225700e-02,    3.1463820e-01);
    positions[545]           = Vec3( -2.7582100e-01,  -8.4479800e-02,    1.8809180e-01);
    positions[546]           = Vec3( -8.3433230e-01,  -5.5202940e-01,    2.1864880e-01);
    positions[547]           = Vec3( -8.2396710e-01,  -5.4694370e-01,    3.1266070e-01);
    positions[548]           = Vec3( -9.3264700e-01,  -5.5452020e-01,    1.9359510e-01);
    positions[549]           = Vec3( -3.7479050e-01,   2.2505660e-01,    7.1205330e-01);
    positions[550]           = Vec3( -4.7509020e-01,   2.3675960e-01,    7.1906840e-01);
    positions[551]           = Vec3( -3.3344270e-01,   3.0911900e-01,    7.2096390e-01);
    positions[552]           = Vec3(  5.4909720e-01,  -6.8048160e-01,    7.2400200e-02);
    positions[553]           = Vec3(  6.0527360e-01,  -6.6696760e-01,    1.5170900e-01);
    positions[554]           = Vec3(  5.8614280e-01,  -6.1178520e-01,    1.7524700e-02);
    positions[555]           = Vec3( -2.3127640e-01,   9.0287820e-01,   -1.3411380e-01);
    positions[556]           = Vec3( -2.8615520e-01,   9.5668910e-01,   -1.9830460e-01);
    positions[557]           = Vec3( -2.9306830e-01,   8.7146310e-01,   -6.8234400e-02);
    positions[558]           = Vec3( -5.4794480e-01,   6.9927600e-02,    4.9211700e-02);
    positions[559]           = Vec3( -4.8467110e-01,   1.3673600e-02,    9.6662900e-02);
    positions[560]           = Vec3( -5.5944570e-01,   3.5041600e-02,   -4.0422400e-02);
    positions[561]           = Vec3( -4.0842490e-01,  -6.1610810e-01,    5.3013490e-01);
    positions[562]           = Vec3( -3.5055240e-01,  -6.7988460e-01,    4.9398580e-01);
    positions[563]           = Vec3( -4.6296070e-01,  -6.7880320e-01,    5.8633470e-01);
    positions[564]           = Vec3(  4.6585780e-01,   7.8746100e-01,   -1.2817710e-01);
    positions[565]           = Vec3(  5.3858490e-01,   8.3094890e-01,   -7.7410200e-02);
    positions[566]           = Vec3(  4.0552000e-01,   7.4979180e-01,   -6.1891900e-02);
    positions[567]           = Vec3( -1.6560700e-02,  -3.7062430e-01,   -3.6569060e-01);
    positions[568]           = Vec3( -9.0792700e-02,  -4.1378610e-01,   -3.2720710e-01);
    positions[569]           = Vec3(  5.1374900e-02,  -4.3774530e-01,   -3.4403280e-01);
    positions[570]           = Vec3( -5.9512760e-01,   1.7073000e-02,   -2.2772060e-01);
    positions[571]           = Vec3( -5.8225940e-01,   7.1421900e-02,   -3.0604790e-01);
    positions[572]           = Vec3( -6.2819960e-01,  -6.3276600e-02,   -2.6202260e-01);
    positions[573]           = Vec3( -4.7641750e-01,  -4.2323550e-01,    8.9604240e-01);
    positions[574]           = Vec3( -5.4796980e-01,  -4.1341290e-01,    8.3129580e-01);
    positions[575]           = Vec3( -4.6422920e-01,  -5.2061790e-01,    8.9834640e-01);
    positions[576]           = Vec3(  1.4489300e-02,  -8.9340740e-01,   -3.4831200e-02);
    positions[577]           = Vec3( -6.9252500e-02,  -9.1064710e-01,   -8.0576000e-02);
    positions[578]           = Vec3(  8.8146500e-02,  -9.1521790e-01,   -9.6184600e-02);
    positions[579]           = Vec3( -6.0237270e-01,   6.8170090e-01,    6.7672100e-02);
    positions[580]           = Vec3( -6.3353490e-01,   6.5944010e-01,   -1.4737000e-02);
    positions[581]           = Vec3( -6.7945440e-01,   7.0260310e-01,    1.2553510e-01);
    positions[582]           = Vec3( -7.9759390e-01,  -4.8566970e-01,   -8.8075620e-01);
    positions[583]           = Vec3( -7.6587590e-01,  -4.4277470e-01,   -9.5861500e-01);
    positions[584]           = Vec3( -8.8262650e-01,  -4.4354590e-01,   -8.6497650e-01);
    positions[585]           = Vec3(  6.0913180e-01,   7.5063640e-01,   -3.7944500e-01);
    positions[586]           = Vec3(  6.8958950e-01,   8.0236210e-01,   -3.5044320e-01);
    positions[587]           = Vec3(  5.5351750e-01,   7.5362410e-01,   -2.9669720e-01);
    positions[588]           = Vec3(  7.4485800e-01,   5.3041050e-01,   -4.4708420e-01);
    positions[589]           = Vec3(  7.0182180e-01,   6.1806940e-01,   -4.3652910e-01);
    positions[590]           = Vec3(  8.0156580e-01,   5.2857300e-01,   -5.2411300e-01);
    positions[591]           = Vec3( -6.9004280e-01,  -5.9012070e-01,   -2.9270410e-01);
    positions[592]           = Vec3( -7.1539690e-01,  -6.8384200e-01,   -2.8572180e-01);
    positions[593]           = Vec3( -5.9319910e-01,  -5.8219810e-01,   -2.7391860e-01);
    positions[594]           = Vec3( -2.0769030e-01,  -9.0263320e-01,    8.2559380e-01);
    positions[595]           = Vec3( -1.2326710e-01,  -9.0347650e-01,    7.7889800e-01);
    positions[596]           = Vec3( -2.4674410e-01,  -8.1114260e-01,    8.1400270e-01);
    positions[597]           = Vec3( -5.9770390e-01,  -2.5353030e-01,    3.7815410e-01);
    positions[598]           = Vec3( -5.7799760e-01,  -1.8503970e-01,    4.3781640e-01);
    positions[599]           = Vec3( -6.3056510e-01,  -1.9169960e-01,    3.0646360e-01);
    positions[600]           = Vec3(  2.5756560e-01,  -9.0983610e-01,   -2.2681580e-01);
    positions[601]           = Vec3(  3.3909840e-01,  -9.6122750e-01,   -2.1952540e-01);
    positions[602]           = Vec3(  2.5286730e-01,  -8.8095350e-01,   -3.1936560e-01);
    positions[603]           = Vec3( -5.7980030e-01,   4.5624440e-01,   -4.8053250e-01);
    positions[604]           = Vec3( -5.0283550e-01,   4.0235700e-01,   -4.7629430e-01);
    positions[605]           = Vec3( -5.5234760e-01,   5.5030960e-01,   -4.6445780e-01);
    positions[606]           = Vec3(  2.6417710e-01,   3.6149920e-01,    5.5726940e-01);
    positions[607]           = Vec3(  1.8510390e-01,   3.4649300e-01,    6.1450570e-01);
    positions[608]           = Vec3(  2.5328370e-01,   4.4988030e-01,    5.1753010e-01);
    positions[609]           = Vec3(  8.0108540e-01,  -7.3935090e-01,   -4.6186460e-01);
    positions[610]           = Vec3(  8.1693510e-01,  -7.9012220e-01,   -3.8198140e-01);
    positions[611]           = Vec3(  8.8304810e-01,  -6.9238110e-01,   -4.8097940e-01);
    positions[612]           = Vec3( -5.8628640e-01,   1.5133800e-02,   -5.2805090e-01);
    positions[613]           = Vec3( -5.0874980e-01,   5.8718200e-02,   -5.5884230e-01);
    positions[614]           = Vec3( -6.4503990e-01,   1.9133400e-02,   -6.0165090e-01);
    positions[615]           = Vec3(  7.6453220e-01,  -5.9994620e-01,    2.8797170e-01);
    positions[616]           = Vec3(  7.0859250e-01,  -5.4012040e-01,    3.3515640e-01);
    positions[617]           = Vec3(  7.9449730e-01,  -6.7260900e-01,    3.4844210e-01);
    positions[618]           = Vec3( -4.1271350e-01,   6.8162960e-01,   -6.2517570e-01);
    positions[619]           = Vec3( -3.4841290e-01,   6.1054470e-01,   -6.4194430e-01);
    positions[620]           = Vec3( -3.5808100e-01,   7.5958210e-01,   -6.0333290e-01);
    positions[621]           = Vec3( -2.3867290e-01,   5.9441400e-02,    9.1386800e-01);
    positions[622]           = Vec3( -2.9103650e-01,  -1.4337800e-02,    9.5259360e-01);
    positions[623]           = Vec3( -2.8602000e-01,   1.0405050e-01,    8.3648420e-01);
    positions[624]           = Vec3(  6.2908620e-01,  -6.6369160e-01,   -8.8313160e-01);
    positions[625]           = Vec3(  5.3309000e-01,  -6.6824080e-01,   -8.8386380e-01);
    positions[626]           = Vec3(  6.6687380e-01,  -7.2037270e-01,   -8.1674370e-01);
    positions[627]           = Vec3(  2.5101170e-01,  -8.8838680e-01,    2.2900940e-01);
    positions[628]           = Vec3(  2.4302200e-01,  -8.1686710e-01,    1.6969450e-01);
    positions[629]           = Vec3(  3.4457660e-01,  -8.9596990e-01,    2.4839760e-01);
    positions[630]           = Vec3( -9.1418940e-01,   8.0389630e-01,    7.8826000e-01);
    positions[631]           = Vec3( -8.3833600e-01,   7.4209380e-01,    7.8290720e-01);
    positions[632]           = Vec3( -9.9161100e-01,   7.5608500e-01,    8.2971860e-01);
    positions[633]           = Vec3(  7.9708930e-01,  -3.2882190e-01,    7.1789600e-01);
    positions[634]           = Vec3(  8.5609970e-01,  -2.5716920e-01,    7.4938090e-01);
    positions[635]           = Vec3(  7.9853320e-01,  -3.2248890e-01,    6.2155040e-01);
    positions[636]           = Vec3(  7.9743030e-01,  -6.0061740e-01,    7.6822330e-01);
    positions[637]           = Vec3(  8.2105340e-01,  -5.0895770e-01,    7.5902860e-01);
    positions[638]           = Vec3(  7.2970170e-01,  -6.0508550e-01,    8.3860140e-01);
    positions[639]           = Vec3( -1.1738970e-01,  -5.9305270e-01,    7.0381050e-01);
    positions[640]           = Vec3( -1.5290840e-01,  -6.5518590e-01,    6.3431800e-01);
    positions[641]           = Vec3( -1.6038250e-01,  -5.0776740e-01,    6.8496070e-01);
    positions[642]           = Vec3(  5.8567050e-01,   3.6131160e-01,    3.0656670e-01);
    positions[643]           = Vec3(  5.1450330e-01,   4.2381370e-01,    2.6162660e-01);
    positions[644]           = Vec3(  5.3597340e-01,   3.2574600e-01,    3.8272470e-01);
    positions[645]           = Vec3( -7.5114680e-01,   3.5944460e-01,    2.4369600e-01);
    positions[646]           = Vec3( -8.1938720e-01,   4.1907000e-01,    2.7265900e-01);
    positions[647]           = Vec3( -7.8315770e-01,   3.2541650e-01,    1.6165560e-01);

    std::string platformName;
    platformName = "CUDA";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName(platformName));

    context.setPositions(positions);

    if (testName == "testSystemMultipoleMoments") {
        amoebaMultipoleForce->getSystemMultipoleMoments(context, outputMultipoleMoments);
    } else if (testName == "testMultipoleGridPotential") {
        amoebaMultipoleForce->getElectrostaticPotential(inputGrid, context, outputGridPotential);
    } else {
        State state               = context.getState(State::Forces | State::Energy);
        forces                    = state.getForces();
        energy                    = state.getPotentialEnergy();
    }

    if (nonbondedMethod == AmoebaMultipoleForce::PME) {
        double actualAlpha;
        int actualSize[3];
        amoebaMultipoleForce->getPMEParametersInContext(context, actualAlpha, actualSize[0], actualSize[1], actualSize[2]);
        ASSERT_EQUAL_TOL(amoebaMultipoleForce->getAEwald(), actualAlpha, 1e-5);
        for (int i = 0; i < 3; i++) {
            ASSERT(actualSize[i] >= inputPmeGridDimension);
            ASSERT(actualSize[i] < inputPmeGridDimension+10);
        }
    }
}

// test multipole mutual polarization using PME for box of water

static void testPMEMutualPolarizationLargeWater() {

    std::string testName      = "testPMEMutualPolarizationLargeWater";

    int numberOfParticles     = 648;
    int inputPmeGridDimension = 24;
    double cutoff             = 0.70;
    std::vector<Vec3> forces;
    double energy;
    std::vector< double > outputMultipoleMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;

    setupAndGetForcesEnergyMultipoleLargeWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                               cutoff, inputPmeGridDimension, testName,
                                               forces, energy, outputMultipoleMoments, inputGrid, outputGridPotential);
    static std::vector<Vec3> expectedForces; // Static to work around bug in Visual Studio that makes compilation very very slow.
    expectedForces.resize(numberOfParticles);

    double expectedEnergy     =  -1.3268930e+04;

    expectedForces[0]         = Vec3( -5.2764003e+02,   6.5154502e+02,   7.8683284e+02);
    expectedForces[1]         = Vec3(  3.3391292e+02,  -1.1162521e+03,  -2.9173021e+02);
    expectedForces[2]         = Vec3(  9.4613310e+01,  -5.8988099e+01,  -6.2871010e+02);
    expectedForces[3]         = Vec3( -1.6179658e+03,   1.8222798e+02,  -1.0373083e+02);
    expectedForces[4]         = Vec3(  1.0177962e+03,   1.7118957e+02,   4.0084976e+02);
    expectedForces[5]         = Vec3(  3.0582540e+02,  -5.5215191e+02,   4.8287747e+01);
    expectedForces[6]         = Vec3( -5.4931641e+01,  -6.7689368e+02,  -1.2260977e+03);
    expectedForces[7]         = Vec3( -2.9961754e+02,   1.3602079e+03,  -2.2087612e+02);
    expectedForces[8]         = Vec3( -3.1788045e+02,  -1.3540025e+02,   1.1250958e+03);
    expectedForces[9]         = Vec3(  1.0275111e+03,   1.9279625e+02,  -6.0872252e+02);
    expectedForces[10]        = Vec3( -5.7512410e+02,  -5.6280482e+01,   7.5039261e+02);
    expectedForces[11]        = Vec3( -4.7291981e+02,  -3.5875411e+01,   5.2335135e+00);
    expectedForces[12]        = Vec3( -1.5619899e+02,  -1.3743370e+02,   1.3587339e+03);
    expectedForces[13]        = Vec3( -4.7423876e+02,   1.0032873e+02,  -2.0875847e+02);
    expectedForces[14]        = Vec3(  7.5621178e+02,   3.7979419e+01,  -4.1221012e+02);
    expectedForces[15]        = Vec3(  6.7318876e+02,   4.9937148e+02,   1.1665351e+03);
    expectedForces[16]        = Vec3( -7.0858091e+01,   8.4758233e+01,  -1.3126942e+03);
    expectedForces[17]        = Vec3( -1.0636296e+03,  -1.6254851e+02,  -1.0982809e+02);
    expectedForces[18]        = Vec3( -1.0249489e+03,   6.7073895e+02,  -2.0754278e+02);
    expectedForces[19]        = Vec3(  5.6396078e+02,  -1.1702137e+03,   1.5324295e+02);
    expectedForces[20]        = Vec3(  6.1635144e+02,   5.6122242e+02,  -7.9410893e+01);
    expectedForces[21]        = Vec3( -1.5785582e+02,  -1.0857764e+03,  -1.1010131e+03);
    expectedForces[22]        = Vec3(  2.3733181e+01,   2.0158373e+02,   3.2583783e+02);
    expectedForces[23]        = Vec3(  1.8855396e+02,   1.2408889e+03,   4.5759149e+02);
    expectedForces[24]        = Vec3(  9.3643408e+02,  -6.1331107e+01,  -3.8601431e+02);
    expectedForces[25]        = Vec3( -3.1972504e+02,  -1.5783340e+02,  -1.3838003e+02);
    expectedForces[26]        = Vec3( -5.1178039e+02,  -2.5998898e+02,   6.1021975e+02);
    expectedForces[27]        = Vec3(  9.1836250e+02,  -4.3288112e+02,  -9.5579501e+02);
    expectedForces[28]        = Vec3( -9.3325184e+02,   6.2859013e+01,   5.0743198e+02);
    expectedForces[29]        = Vec3(  3.1839963e+02,   5.4102155e+02,   1.1261345e+03);
    expectedForces[30]        = Vec3(  7.3842640e+02,   1.0907054e+02,  -2.1126612e+02);
    expectedForces[31]        = Vec3( -1.8470175e+02,  -4.0377182e+01,   1.2464073e+02);
    expectedForces[32]        = Vec3( -5.8931386e+02,  -3.5021489e+02,  -7.2769901e+01);
    expectedForces[33]        = Vec3(  1.4804452e+02,  -1.2138869e+02,  -1.2423820e+03);
    expectedForces[34]        = Vec3( -2.4850935e+02,   2.3113218e+01,   1.4041979e+02);
    expectedForces[35]        = Vec3(  1.3188525e+02,  -6.1243481e+02,   6.4886373e+02);
    expectedForces[36]        = Vec3(  5.0217683e+02,   9.2857480e+02,  -7.2500648e+02);
    expectedForces[37]        = Vec3(  9.9446320e+01,  -3.0337270e+02,   8.0439741e+02);
    expectedForces[38]        = Vec3( -1.7167491e+02,  -9.3911948e+02,   5.6412500e+01);
    expectedForces[39]        = Vec3( -6.6143649e+02,   1.0876882e+03,   3.3801466e+02);
    expectedForces[40]        = Vec3(  5.6748150e+02,  -2.2459572e+02,   3.8360728e+02);
    expectedForces[41]        = Vec3(  7.7688834e+02,  -6.7341171e+02,  -6.4292796e+02);
    expectedForces[42]        = Vec3(  1.0253263e+02,  -1.3212797e+03,   9.3279892e+02);
    expectedForces[43]        = Vec3(  1.5424247e+02,  -1.0197346e+01,  -8.0230728e+02);
    expectedForces[44]        = Vec3( -6.6499384e+02,   1.2534353e+03,  -1.6748341e+02);
    expectedForces[45]        = Vec3( -1.9282945e+02,   5.1764789e+02,   1.5108154e+03);
    expectedForces[46]        = Vec3( -4.9204245e+02,   8.8664089e+02,  -8.8381078e+02);
    expectedForces[47]        = Vec3(  1.3207641e+01,  -5.7765859e+02,  -3.5921690e+02);
    expectedForces[48]        = Vec3(  4.1691955e+02,  -1.1292452e+03,  -6.2339128e+02);
    expectedForces[49]        = Vec3(  4.8613761e+02,   6.8325845e+02,   2.3424652e+02);
    expectedForces[50]        = Vec3(  3.1827501e+01,   6.8442319e+02,  -5.8685438e+02);
    expectedForces[51]        = Vec3(  1.1834790e+03,  -1.3418299e+03,   5.6280882e+02);
    expectedForces[52]        = Vec3( -5.8629616e+02,  -6.9100489e+01,  -9.7784543e+02);
    expectedForces[53]        = Vec3( -4.0294748e+02,   5.0669032e+02,   5.3893512e+02);
    expectedForces[54]        = Vec3(  6.3555947e+01,  -9.0345643e+02,   7.2235587e+02);
    expectedForces[55]        = Vec3(  4.9353706e+02,   4.3952680e+02,  -1.5413911e+02);
    expectedForces[56]        = Vec3( -3.4663273e+02,   8.4708544e+02,  -8.8983178e+01);
    expectedForces[57]        = Vec3(  1.5746383e+03,   5.8209931e+02,   3.6645765e+02);
    expectedForces[58]        = Vec3( -4.0076588e+02,  -1.2311235e+02,  -1.6226002e+02);
    expectedForces[59]        = Vec3( -3.3837369e+02,  -5.3890252e+02,  -5.5997831e+02);
    expectedForces[60]        = Vec3(  8.6512912e+02,  -1.0775827e+03,  -3.6920406e+02);
    expectedForces[61]        = Vec3(  1.0588621e+01,   4.4531974e+02,   1.1838800e+02);
    expectedForces[62]        = Vec3( -7.0505255e+02,   1.5737739e+02,   3.0134165e+02);
    expectedForces[63]        = Vec3( -1.2833932e+03,   5.1206074e+02,  -9.3423828e+01);
    expectedForces[64]        = Vec3(  2.0863538e+02,  -6.7108922e+02,   3.5777952e+02);
    expectedForces[65]        = Vec3(  4.9710436e+02,  -6.8550896e+01,  -2.7502680e+02);
    expectedForces[66]        = Vec3(  4.0388289e+02,   1.0051571e+03,  -2.3823633e+02);
    expectedForces[67]        = Vec3( -3.3895858e+02,  -6.7577709e+02,   3.1521120e+02);
    expectedForces[68]        = Vec3( -4.8763670e+02,  -2.9427128e+02,   2.9500438e+02);
    expectedForces[69]        = Vec3( -8.4080783e+02,  -1.7777886e+02,   6.5087880e+02);
    expectedForces[70]        = Vec3(  2.1136356e+02,   6.6812142e+01,  -5.5437438e+02);
    expectedForces[71]        = Vec3( -1.3829189e+01,   7.0936421e+02,  -5.0847380e+02);
    expectedForces[72]        = Vec3( -1.5684247e+03,   1.6837339e+02,  -4.6044105e+02);
    expectedForces[73]        = Vec3(  9.9163605e+02,  -4.1397181e+02,  -3.9624142e+02);
    expectedForces[74]        = Vec3(  6.7355000e+02,   5.7652252e+02,   1.0565003e+03);
    expectedForces[75]        = Vec3( -1.3722105e+03,  -1.4035644e+03,   2.7882380e+02);
    expectedForces[76]        = Vec3(  5.2187588e+02,   3.9602986e+02,  -2.2564009e+02);
    expectedForces[77]        = Vec3(  2.6381283e+02,   4.2945748e+02,   1.8339201e+02);
    expectedForces[78]        = Vec3(  5.3099327e+02,   1.7312849e+03,   3.5385638e+01);
    expectedForces[79]        = Vec3( -2.1320977e+02,  -1.1795453e+03,  -5.1691859e+02);
    expectedForces[80]        = Vec3( -7.1445218e+02,  -5.8918750e+02,   7.6494816e+02);
    expectedForces[81]        = Vec3(  1.0849545e+03,   6.5058603e+02,  -5.6051049e+02);
    expectedForces[82]        = Vec3( -4.6986151e+02,  -3.3868080e+02,   7.9000645e+02);
    expectedForces[83]        = Vec3( -3.8404718e+02,  -4.4915037e+02,  -4.2069600e+02);
    expectedForces[84]        = Vec3( -3.0678477e+02,  -7.2558567e+02,   2.0276275e+02);
    expectedForces[85]        = Vec3(  6.2969625e+01,   3.4257826e+02,  -6.8348290e+02);
    expectedForces[86]        = Vec3(  3.1042613e+02,   9.6563445e+01,  -1.8788848e+02);
    expectedForces[87]        = Vec3(  2.4427535e+01,   5.8460648e+02,  -1.3720684e+03);
    expectedForces[88]        = Vec3( -4.9904469e+02,  -9.1309605e+02,   1.4350031e+02);
    expectedForces[89]        = Vec3( -4.0910847e+01,   5.2845046e+01,   2.7341591e+02);
    expectedForces[90]        = Vec3( -4.5190765e+02,   2.2820859e+02,  -1.1052181e+03);
    expectedForces[91]        = Vec3(  9.1626336e+01,  -1.0496400e+03,   1.0216011e+02);
    expectedForces[92]        = Vec3( -1.2592678e+02,  -4.4050520e+01,   1.0764584e+03);
    expectedForces[93]        = Vec3(  8.4216135e+02,   1.7603150e+02,   1.3698861e+03);
    expectedForces[94]        = Vec3(  1.2583734e+02,   3.4343971e+02,  -7.3070182e+02);
    expectedForces[95]        = Vec3( -4.5597370e+02,  -6.8170738e+02,  -5.2258658e+02);
    expectedForces[96]        = Vec3(  1.3295639e+03,   2.5031118e+02,  -3.5171095e+02);
    expectedForces[97]        = Vec3( -2.0136001e+02,  -2.2978012e+02,   2.1122610e+02);
    expectedForces[98]        = Vec3( -1.1279722e+03,  -9.9794730e+02,  -1.1321006e+02);
    expectedForces[99]        = Vec3(  2.9553179e+02,  -2.0727600e+02,  -1.8912317e+03);
    expectedForces[100]       = Vec3( -9.3777487e+02,   4.3417972e+02,   5.7436670e+02);
    expectedForces[101]       = Vec3(  3.1120364e+02,  -4.3241206e+00,   5.6666340e+02);
    expectedForces[102]       = Vec3( -6.6764840e+02,  -2.3563274e+02,   9.9017992e+02);
    expectedForces[103]       = Vec3(  8.2857598e+02,   1.9030059e+02,  -2.4295780e+02);
    expectedForces[104]       = Vec3( -2.6709334e+02,   7.3284851e+01,  -3.9157038e+02);
    expectedForces[105]       = Vec3(  6.6235629e+02,  -8.4502811e+02,  -5.3670716e+02);
    expectedForces[106]       = Vec3(  3.8449086e+02,   8.0775263e+02,   5.2003911e+02);
    expectedForces[107]       = Vec3( -3.4338589e+02,   1.7870086e+02,  -6.1047615e+01);
    expectedForces[108]       = Vec3(  4.3523975e+02,   2.8776092e+02,  -1.5333067e+03);
    expectedForces[109]       = Vec3(  1.9214022e+01,  -9.9449061e+01,   3.9153169e+02);
    expectedForces[110]       = Vec3( -7.7525496e+02,   7.7309387e+01,   7.8886699e+02);
    expectedForces[111]       = Vec3( -1.1510224e+03,  -2.8207337e+02,   1.4866455e+03);
    expectedForces[112]       = Vec3(  1.1281316e+03,  -3.0961130e+02,  -1.3894255e+02);
    expectedForces[113]       = Vec3( -2.4829062e+02,   5.9466066e+02,  -9.0749265e+02);
    expectedForces[114]       = Vec3( -1.6228531e+02,  -4.1281378e+02,  -1.4860946e+03);
    expectedForces[115]       = Vec3( -2.9071192e+02,   7.5335271e+02,   6.5205720e+02);
    expectedForces[116]       = Vec3(  3.2042337e+02,  -4.0237231e+02,   1.0068554e+03);
    expectedForces[117]       = Vec3( -1.4081272e+03,  -5.5227579e+02,   2.8376428e+02);
    expectedForces[118]       = Vec3(  9.6978276e+02,  -9.6232683e+02,   2.3399918e+02);
    expectedForces[119]       = Vec3(  2.0781312e+02,   9.0531947e+01,   8.3513783e+01);
    expectedForces[120]       = Vec3(  2.9200105e+02,   1.0746865e+03,  -6.0773727e+02);
    expectedForces[121]       = Vec3( -8.7546167e+02,   3.1969683e+02,   8.1033199e+01);
    expectedForces[122]       = Vec3( -5.3422658e+01,  -3.8456585e+02,   1.4595003e+02);
    expectedForces[123]       = Vec3( -9.5134908e+01,   8.6818098e+02,  -1.3111403e+03);
    expectedForces[124]       = Vec3( -7.2576830e+02,   6.3811172e+01,   7.0993488e+02);
    expectedForces[125]       = Vec3(  5.6200572e+02,  -1.0972338e+03,  -3.5278717e+02);
    expectedForces[126]       = Vec3( -1.6903718e+02,   2.6804112e+02,  -8.9811032e+02);
    expectedForces[127]       = Vec3(  7.0934165e+02,  -3.2328700e+02,   1.0037562e+02);
    expectedForces[128]       = Vec3( -2.5937738e+00,  -3.3758145e+02,   3.6354335e+02);
    expectedForces[129]       = Vec3( -1.1993601e+03,   1.1618922e+02,  -1.1329865e+03);
    expectedForces[130]       = Vec3(  8.3480845e+02,   5.2646427e+02,   9.1816763e+01);
    expectedForces[131]       = Vec3(  1.9407174e+02,   7.3877492e+00,   8.8648162e+02);
    expectedForces[132]       = Vec3( -1.0820540e+03,   2.7482390e+02,  -1.4876848e+03);
    expectedForces[133]       = Vec3(  9.5717027e+02,  -4.7758356e+02,   4.0692835e+02);
    expectedForces[134]       = Vec3(  2.4979896e+02,   7.1530552e+02,   6.5834243e+02);
    expectedForces[135]       = Vec3( -2.7549874e+02,   1.6613806e+03,  -5.0628865e+02);
    expectedForces[136]       = Vec3( -1.5960544e+02,  -4.5920421e+02,   9.4068762e+02);
    expectedForces[137]       = Vec3(  6.9764295e+01,  -4.4024430e+02,  -2.3545402e+02);
    expectedForces[138]       = Vec3(  2.0207569e+03,  -2.6846158e+02,   3.3977998e+02);
    expectedForces[139]       = Vec3( -9.9905342e+02,   9.5155462e+02,   8.3259854e+01);
    expectedForces[140]       = Vec3( -8.0836255e+02,  -6.3005196e+02,   1.2603525e+02);
    expectedForces[141]       = Vec3(  5.2279913e+02,   1.5356061e+03,  -1.4399009e+02);
    expectedForces[142]       = Vec3( -1.3088783e+02,  -9.0861449e+02,   9.7040422e+02);
    expectedForces[143]       = Vec3(  2.8006682e+01,  -1.1783245e+03,  -3.7285390e+02);
    expectedForces[144]       = Vec3(  6.1040927e+01,   2.1747286e+01,   9.0366011e+02);
    expectedForces[145]       = Vec3(  7.0396235e+02,  -2.3865001e+02,  -6.8794682e+02);
    expectedForces[146]       = Vec3( -5.1132933e+02,   7.7102556e+02,  -4.7805042e+02);
    expectedForces[147]       = Vec3(  6.2879945e+02,  -6.9758821e+02,   1.0983043e+03);
    expectedForces[148]       = Vec3( -9.2449123e+02,   5.0248042e+02,  -2.0065211e+02);
    expectedForces[149]       = Vec3( -3.9485154e+02,   3.4679511e+02,  -7.3889603e+02);
    expectedForces[150]       = Vec3(  4.9834751e+02,  -1.0117501e+03,   1.2510946e+03);
    expectedForces[151]       = Vec3( -1.4942699e+02,  -7.4362790e+01,  -6.8008427e+02);
    expectedForces[152]       = Vec3( -1.8353666e+02,   5.1697483e+02,  -8.4397292e+01);
    expectedForces[153]       = Vec3(  7.0253008e+02,  -1.5361338e+02,   8.6830836e+02);
    expectedForces[154]       = Vec3( -4.1108898e+02,  -6.8181350e+02,  -5.6516951e+02);
    expectedForces[155]       = Vec3(  3.9294374e+02,   3.7294948e+02,  -4.9795080e+02);
    expectedForces[156]       = Vec3( -7.3716758e+02,  -5.5465931e+02,  -9.3434525e+02);
    expectedForces[157]       = Vec3(  1.8240105e+02,   1.5494359e+02,   1.3626243e+02);
    expectedForces[158]       = Vec3(  9.3172059e+02,   1.4431275e+02,   6.6916065e+02);
    expectedForces[159]       = Vec3(  7.3399513e+01,  -3.6000187e+02,   1.9853052e+03);
    expectedForces[160]       = Vec3( -2.4842520e+02,  -4.0819514e+02,  -8.0624794e+02);
    expectedForces[161]       = Vec3(  4.3585911e+02,  -2.6521764e+02,  -6.7985842e+02);
    expectedForces[162]       = Vec3( -8.2255783e+02,   1.1430906e+03,   1.2991183e+03);
    expectedForces[163]       = Vec3(  5.6437196e+02,  -1.5391213e+02,  -5.0013503e+02);
    expectedForces[164]       = Vec3(  3.8485296e+02,  -1.1754088e+02,  -1.4315511e+03);
    expectedForces[165]       = Vec3(  9.8205241e+02,   1.3278205e+01,   3.0367584e+02);
    expectedForces[166]       = Vec3( -3.2909916e+02,   4.5169387e+01,  -3.2019488e+02);
    expectedForces[167]       = Vec3( -4.2382795e+02,  -8.1938351e+02,  -1.5360675e+02);
    expectedForces[168]       = Vec3(  9.1689997e+02,   1.4443113e+03,  -1.8886795e+02);
    expectedForces[169]       = Vec3( -1.2161000e+03,  -2.3826187e+02,   3.4004889e+02);
    expectedForces[170]       = Vec3( -2.2837927e+02,  -9.8334152e+02,   3.0313790e+02);
    expectedForces[171]       = Vec3(  5.2598040e+02,  -1.2766936e+03,   8.1247934e+01);
    expectedForces[172]       = Vec3( -6.7839687e+01,   6.9220160e+02,   3.1975320e+02);
    expectedForces[173]       = Vec3( -4.8519844e+02,   4.3008979e+02,  -7.6788762e+02);
    expectedForces[174]       = Vec3(  3.9550092e+01,   1.0900122e+03,   5.0955281e+02);
    expectedForces[175]       = Vec3( -6.3000735e+02,  -7.2335263e+02,  -2.2590661e+02);
    expectedForces[176]       = Vec3(  1.5048762e+02,  -4.0086126e+02,  -1.2159849e+02);
    expectedForces[177]       = Vec3( -7.1939527e+02,  -1.0640815e+03,  -3.1826393e+02);
    expectedForces[178]       = Vec3(  7.7492135e+02,   1.0620129e+02,   3.6733165e+02);
    expectedForces[179]       = Vec3( -3.5466457e+02,   7.0081887e+02,   2.2629152e+02);
    expectedForces[180]       = Vec3(  3.0338275e+02,  -1.1243396e+02,  -1.3408996e+03);
    expectedForces[181]       = Vec3(  9.0394088e+00,   4.8818819e+02,   7.3781403e+02);
    expectedForces[182]       = Vec3( -5.6111982e+01,  -2.5224925e+02,   3.0344460e+02);
    expectedForces[183]       = Vec3( -7.4693225e+02,   1.1370727e+03,   9.3958175e+02);
    expectedForces[184]       = Vec3(  3.1310327e+02,  -1.2133761e+03,   2.6768166e+02);
    expectedForces[185]       = Vec3(  6.0706899e+02,   9.8403362e+01,  -5.0781782e+02);
    expectedForces[186]       = Vec3( -6.2721609e+02,   1.7751567e+02,   2.6230643e+02);
    expectedForces[187]       = Vec3(  3.5425475e+02,  -3.3336896e+02,   2.9529930e+02);
    expectedForces[188]       = Vec3(  4.7499029e+02,   4.5552354e+02,  -2.4939310e+02);
    expectedForces[189]       = Vec3(  1.1556321e+02,   2.5995186e+02,  -8.4171376e+02);
    expectedForces[190]       = Vec3( -4.5095184e+02,  -1.6639362e+01,   4.2956175e+02);
    expectedForces[191]       = Vec3( -7.5904218e+01,   3.9416063e+02,   4.2573815e+02);
    expectedForces[192]       = Vec3(  2.3041013e+02,  -1.2883920e+03,  -5.1943189e+02);
    expectedForces[193]       = Vec3(  8.6086944e+01,   1.8626094e+02,   3.7641641e+02);
    expectedForces[194]       = Vec3( -5.9266805e+02,   4.4766229e+02,  -4.5835497e+02);
    expectedForces[195]       = Vec3(  3.9368476e+02,   1.1050329e+03,  -1.3906133e+03);
    expectedForces[196]       = Vec3( -1.2912212e+03,  -2.9660534e+02,   3.9659215e+02);
    expectedForces[197]       = Vec3( -2.4803663e+01,  -2.7945607e+02,   3.5315630e+02);
    expectedForces[198]       = Vec3(  1.4147479e+02,   3.2766723e+02,   1.1137899e+03);
    expectedForces[199]       = Vec3(  7.8383167e+01,   2.2385715e+00,  -1.0346820e+03);
    expectedForces[200]       = Vec3( -2.9231551e+02,  -5.4648720e+02,  -6.5500265e+02);
    expectedForces[201]       = Vec3(  1.3708284e+03,  -1.0716376e+03,   4.5781190e+02);
    expectedForces[202]       = Vec3( -7.5533654e+02,   6.5396686e+02,  -4.2683276e+02);
    expectedForces[203]       = Vec3( -6.5558825e+02,   3.1164131e+02,   2.6323555e+02);
    expectedForces[204]       = Vec3(  1.3173432e+03,   1.3775889e+03,   1.2933072e+02);
    expectedForces[205]       = Vec3( -4.0922505e+02,  -7.5709333e+02,  -3.4967434e+02);
    expectedForces[206]       = Vec3(  4.2270572e+01,  -7.4393359e+02,   5.2985443e+02);
    expectedForces[207]       = Vec3( -1.5327880e+03,   7.3186706e+02,  -2.3440058e+02);
    expectedForces[208]       = Vec3(  1.4650210e+03,   2.8661572e+02,  -2.4233049e+02);
    expectedForces[209]       = Vec3(  8.4124688e+02,  -3.4020518e+02,  -3.1773592e+01);
    expectedForces[210]       = Vec3( -1.1480436e+03,  -7.4089732e+02,  -4.8883682e+02);
    expectedForces[211]       = Vec3(  3.5682075e+02,   2.1918110e+01,  -1.2558000e+02);
    expectedForces[212]       = Vec3(  1.6505256e+02,   7.9867933e+02,   9.9760760e+02);
    expectedForces[213]       = Vec3(  1.2477725e+03,  -8.8687803e+02,  -8.2498529e+02);
    expectedForces[214]       = Vec3(  1.7357042e+02,   2.0185456e+02,   8.0965301e+02);
    expectedForces[215]       = Vec3( -9.0937500e+02,   5.5968713e+01,   1.7581913e+02);
    expectedForces[216]       = Vec3(  1.3030711e+03,  -8.2642443e+02,   1.6679926e+02);
    expectedForces[217]       = Vec3( -8.4295510e+02,   2.0237547e+01,   4.0163431e+02);
    expectedForces[218]       = Vec3( -3.2088719e+02,   3.0824476e+02,   9.9166532e+01);
    expectedForces[219]       = Vec3(  1.4585495e+03,  -4.1497503e+01,  -2.2539636e+02);
    expectedForces[220]       = Vec3( -2.6148433e+02,   1.7592075e+02,   2.0638296e+02);
    expectedForces[221]       = Vec3( -9.2100962e+02,   2.6090225e+02,  -8.4461867e+02);
    expectedForces[222]       = Vec3( -9.0597287e+02,  -4.8246138e+02,  -1.8012972e+03);
    expectedForces[223]       = Vec3( -1.1775667e+02,   1.7275816e+02,   1.1286513e+03);
    expectedForces[224]       = Vec3(  1.0682239e+03,   8.0853797e+02,   2.5601774e+02);
    expectedForces[225]       = Vec3( -1.6637330e+02,   7.6165466e+02,  -5.9394315e+02);
    expectedForces[226]       = Vec3(  2.6206832e+02,  -1.5011747e+02,   5.4532792e+02);
    expectedForces[227]       = Vec3( -2.6474701e+02,  -1.4218144e+02,   1.1099125e+02);
    expectedForces[228]       = Vec3( -9.7365994e+02,  -5.5282012e+02,  -6.3158992e+02);
    expectedForces[229]       = Vec3(  1.1482718e+03,   1.4958662e+02,  -5.6523355e+02);
    expectedForces[230]       = Vec3(  4.4084982e+02,   7.9662251e+02,   4.9732039e+02);
    expectedForces[231]       = Vec3( -4.8746806e+02,  -1.5522005e+03,   1.0133483e+03);
    expectedForces[232]       = Vec3(  1.1332903e+03,   1.2836176e+03,  -2.0897082e+02);
    expectedForces[233]       = Vec3( -7.9973372e+01,   6.3066961e+02,  -5.8656505e+02);
    expectedForces[234]       = Vec3( -8.3714021e+02,  -1.5385924e+03,   1.0179672e+03);
    expectedForces[235]       = Vec3(  3.5358323e+02,   9.3554369e+02,   1.6565605e+02);
    expectedForces[236]       = Vec3(  1.7532341e+02,   3.0381323e+01,  -1.2538424e+03);
    expectedForces[237]       = Vec3( -1.0437301e+03,  -4.4804918e+02,   1.5276435e+03);
    expectedForces[238]       = Vec3(  8.5470872e+02,   1.0725550e+03,  -9.7359443e+01);
    expectedForces[239]       = Vec3(  1.8192015e+02,  -7.6967211e+01,  -8.6054596e+02);
    expectedForces[240]       = Vec3(  1.5910180e+03,  -7.9601605e+02,  -8.4030326e+02);
    expectedForces[241]       = Vec3( -6.4638370e+02,   4.0538540e+02,   2.3396878e+01);
    expectedForces[242]       = Vec3( -5.2648918e+02,   1.3689201e+02,   9.3014716e+02);
    expectedForces[243]       = Vec3(  2.2248997e+02,  -6.3647370e+02,  -8.3913128e+02);
    expectedForces[244]       = Vec3( -7.5759832e+02,   4.7985429e+01,   5.8577833e+02);
    expectedForces[245]       = Vec3(  1.2983550e+02,   4.4401985e+02,   2.3556259e+02);
    expectedForces[246]       = Vec3(  7.8075637e+02,   9.9893821e+02,   7.3511307e+02);
    expectedForces[247]       = Vec3( -8.6411159e+02,   2.3812628e+02,   1.9007971e+02);
    expectedForces[248]       = Vec3(  2.0762263e+02,  -7.0123789e+02,  -4.8175733e+02);
    expectedForces[249]       = Vec3( -8.2749187e+02,  -7.8772346e+02,  -9.8536337e+02);
    expectedForces[250]       = Vec3(  8.3905893e+01,  -2.8770518e+02,   6.0690561e+02);
    expectedForces[251]       = Vec3(  7.2280803e+02,   4.9294103e+02,  -1.2915033e+02);
    expectedForces[252]       = Vec3( -1.3615809e+03,  -1.8564754e+02,  -4.2172817e+02);
    expectedForces[253]       = Vec3(  9.2102155e+02,   5.4063967e+02,  -2.7291760e+02);
    expectedForces[254]       = Vec3(  7.0204283e+02,  -3.0338605e+02,   6.6795953e+02);
    expectedForces[255]       = Vec3( -7.2261285e+00,   3.2570478e+02,   5.9010590e+02);
    expectedForces[256]       = Vec3( -8.3167819e+02,   3.0568475e+02,  -9.9230441e+02);
    expectedForces[257]       = Vec3(  2.1074297e+02,  -4.0099946e+01,  -1.2770634e+02);
    expectedForces[258]       = Vec3(  1.1218559e+03,  -8.8304672e+02,   1.0217945e+03);
    expectedForces[259]       = Vec3( -4.4730264e+01,   5.3616918e+02,  -3.2704110e+02);
    expectedForces[260]       = Vec3( -1.1106141e+03,   3.3541536e+02,  -1.0552986e+02);
    expectedForces[261]       = Vec3( -9.3728497e+01,   4.7397796e+01,  -1.1294332e+03);
    expectedForces[262]       = Vec3(  5.1424956e+01,   4.0854995e+02,   5.6832607e+02);
    expectedForces[263]       = Vec3(  1.4286074e+02,   2.9172921e+01,   5.1355833e+02);
    expectedForces[264]       = Vec3( -4.4318133e+02,  -9.5478293e+02,  -1.8067439e+03);
    expectedForces[265]       = Vec3( -2.8083245e+02,   8.1414586e+02,   6.6089404e+02);
    expectedForces[266]       = Vec3(  1.0434876e+02,   3.0664078e+02,   1.3658532e+03);
    expectedForces[267]       = Vec3(  1.2999089e+03,   5.6746563e+02,   1.2200612e+03);
    expectedForces[268]       = Vec3( -5.9546942e+02,   3.6897149e+02,  -8.3331437e+02);
    expectedForces[269]       = Vec3( -1.1094881e+03,  -5.1903335e+02,  -3.5729033e+02);
    expectedForces[270]       = Vec3(  4.1873464e+02,   6.9788983e+02,  -1.3191929e+03);
    expectedForces[271]       = Vec3( -1.0086080e+02,   3.3114447e+02,   6.3419860e+02);
    expectedForces[272]       = Vec3( -2.8783906e+02,  -8.8308011e+02,  -5.4816768e+01);
    expectedForces[273]       = Vec3(  1.0254638e+03,   2.8710026e+02,  -1.9779676e+02);
    expectedForces[274]       = Vec3( -3.5906548e+02,  -5.5531449e+02,  -6.0205492e+01);
    expectedForces[275]       = Vec3( -3.9578349e+02,   5.9217272e+01,   1.0735648e+02);
    expectedForces[276]       = Vec3( -8.7298871e+02,  -1.5584711e+02,   9.2571354e+02);
    expectedForces[277]       = Vec3(  4.1520558e+02,   5.6442126e+02,  -1.8142289e+02);
    expectedForces[278]       = Vec3(  6.5423421e+02,  -4.1631573e+02,   1.0996479e+02);
    expectedForces[279]       = Vec3( -8.0063899e+01,  -1.2880940e+03,   2.9304469e+02);
    expectedForces[280]       = Vec3( -9.5930841e+02,   5.8293485e+02,  -1.1632205e+03);
    expectedForces[281]       = Vec3(  6.3492338e+02,   2.5550931e+02,   1.7380517e+02);
    expectedForces[282]       = Vec3( -1.9188175e+02,   1.0257290e+03,  -8.6438702e+02);
    expectedForces[283]       = Vec3(  4.0063530e+01,  -5.1808902e+02,   4.9556015e+02);
    expectedForces[284]       = Vec3(  5.8250018e+02,  -2.5030700e+02,   1.1762860e+03);
    expectedForces[285]       = Vec3( -4.3195459e+02,   7.4733530e+02,  -1.3002210e+03);
    expectedForces[286]       = Vec3(  2.7567334e+01,  -5.1662654e+02,   3.3258050e+02);
    expectedForces[287]       = Vec3(  6.1700857e+02,   3.0016354e+02,   9.5519139e+02);
    expectedForces[288]       = Vec3(  7.7553094e+02,   8.6453551e+02,  -1.0122996e+03);
    expectedForces[289]       = Vec3( -2.6730796e+02,  -7.7560875e+02,   4.3195620e+02);
    expectedForces[290]       = Vec3( -8.5982146e+02,   1.6692057e+02,   4.3838169e+02);
    expectedForces[291]       = Vec3(  1.2945126e+03,   6.3231748e+02,  -8.3020592e+02);
    expectedForces[292]       = Vec3( -6.9510729e+02,  -3.1930013e+02,  -1.3919425e+01);
    expectedForces[293]       = Vec3( -4.1154200e+02,  -3.3562358e+02,   6.3292682e+02);
    expectedForces[294]       = Vec3( -4.0919783e+02,  -3.8282298e+02,  -4.9125465e+02);
    expectedForces[295]       = Vec3(  6.3932145e+02,  -1.8769713e+01,   9.9241332e+01);
    expectedForces[296]       = Vec3(  8.6847663e+01,   8.7234739e+02,   2.7124112e+02);
    expectedForces[297]       = Vec3( -1.0307576e+03,  -6.2447562e+02,  -1.5796976e+03);
    expectedForces[298]       = Vec3(  6.2464595e+02,   1.0608165e+03,   3.1422521e+01);
    expectedForces[299]       = Vec3(  3.7767184e+02,   4.1170186e+02,   1.0696495e+03);
    expectedForces[300]       = Vec3(  1.0578030e+02,  -1.9544726e+03,  -4.6243957e+02);
    expectedForces[301]       = Vec3( -3.1074896e+02,   8.6738333e+02,  -2.0241448e+02);
    expectedForces[302]       = Vec3( -2.7350519e+02,   6.9945273e+02,   7.8755130e+02);
    expectedForces[303]       = Vec3(  7.8083070e+02,   6.5199614e+02,  -6.6698950e+02);
    expectedForces[304]       = Vec3(  1.3647451e+02,  -4.2490411e+02,   1.3236061e+01);
    expectedForces[305]       = Vec3( -1.1048984e+03,   3.7582184e+02,  -6.4718844e+01);
    expectedForces[306]       = Vec3(  1.5976050e+03,  -8.9091819e+02,   9.7113419e+02);
    expectedForces[307]       = Vec3( -4.5545384e+01,   8.8683300e+02,  -6.5927739e+01);
    expectedForces[308]       = Vec3( -1.3883497e+03,  -4.6171498e+02,  -2.9117829e+02);
    expectedForces[309]       = Vec3( -6.6661140e+02,  -8.1394964e+02,   1.2397900e+03);
    expectedForces[310]       = Vec3( -2.5293546e+02,   1.8568554e+02,  -6.8919479e+02);
    expectedForces[311]       = Vec3(  5.2052057e+02,   1.0288555e+03,   3.1435700e+02);
    expectedForces[312]       = Vec3( -2.4245555e+02,  -1.1100993e+03,  -1.6937710e+03);
    expectedForces[313]       = Vec3(  1.6647470e+02,   9.5272347e+02,   9.2528380e+02);
    expectedForces[314]       = Vec3(  8.6946518e+02,  -1.6295251e+02,   5.7452409e+02);
    expectedForces[315]       = Vec3(  1.4059831e+02,   5.6780959e+02,  -4.6024149e+02);
    expectedForces[316]       = Vec3(  2.3327811e+02,  -2.2376697e+02,   1.3399582e+01);
    expectedForces[317]       = Vec3( -2.6583612e+01,  -4.5801841e+02,   2.9595361e+01);
    expectedForces[318]       = Vec3(  1.2361796e+03,  -1.9473934e+02,  -2.6179421e+02);
    expectedForces[319]       = Vec3( -4.2330105e+02,   5.0768290e+02,   6.8352494e+02);
    expectedForces[320]       = Vec3( -2.0826312e+02,   1.4720747e+02,  -9.8828425e-01);
    expectedForces[321]       = Vec3( -7.3226106e+02,  -1.5366771e+01,   2.7882968e+02);
    expectedForces[322]       = Vec3(  4.2634684e+02,  -6.2926647e+02,   3.6300784e+02);
    expectedForces[323]       = Vec3(  1.7826269e+02,   1.0038378e+02,  -4.2408556e+02);
    expectedForces[324]       = Vec3( -1.1690005e+03,   2.1777241e+02,   9.1980300e+02);
    expectedForces[325]       = Vec3(  5.6841370e+02,  -2.6614377e+02,  -6.4968364e+01);
    expectedForces[326]       = Vec3(  3.5710749e+02,   1.3228373e+02,  -3.6567433e+02);
    expectedForces[327]       = Vec3(  8.1957745e+02,   7.3903486e+02,   3.9487193e+02);
    expectedForces[328]       = Vec3( -8.5354740e+02,   8.9297958e+01,   9.0615539e+01);
    expectedForces[329]       = Vec3( -2.3935807e+02,  -2.2950021e+02,  -4.6193868e+01);
    expectedForces[330]       = Vec3(  8.6120406e+01,   1.4046499e+03,  -1.5899345e+02);
    expectedForces[331]       = Vec3(  4.6319634e+02,  -4.6309406e+02,  -2.2891416e+02);
    expectedForces[332]       = Vec3( -4.2592255e+02,  -2.6503000e+02,   6.1788141e+02);
    expectedForces[333]       = Vec3( -2.4468457e+02,  -7.7827760e+02,   4.2470013e+02);
    expectedForces[334]       = Vec3(  2.6006448e+02,   7.6289112e+01,  -4.3430411e+02);
    expectedForces[335]       = Vec3( -2.5268926e+02,   7.7381529e+02,  -5.0896414e+02);
    expectedForces[336]       = Vec3(  1.0414844e+03,  -6.1885512e+02,   1.4495539e+03);
    expectedForces[337]       = Vec3( -5.9522011e+01,   1.3607073e+03,  -7.3705640e+01);
    expectedForces[338]       = Vec3( -5.9857094e+02,  -2.7213045e+02,  -9.7516268e+02);
    expectedForces[339]       = Vec3(  9.2852178e+02,   3.0121600e+02,  -7.4982031e+00);
    expectedForces[340]       = Vec3( -1.0201472e+03,   1.6269359e+02,   1.9280960e+02);
    expectedForces[341]       = Vec3( -1.8744984e+02,  -4.9790658e+02,   4.2841303e+02);
    expectedForces[342]       = Vec3( -1.0893114e+03,  -4.6044565e+02,  -2.0537532e+02);
    expectedForces[343]       = Vec3(  5.1068148e+02,  -3.0889884e+02,   1.7703226e+01);
    expectedForces[344]       = Vec3(  2.0128423e+02,   5.0813056e+02,  -5.0941772e+01);
    expectedForces[345]       = Vec3( -6.0195519e+02,   1.1710803e+03,  -5.8271481e+02);
    expectedForces[346]       = Vec3(  1.3898239e+02,  -3.2598252e+02,   1.0877023e+03);
    expectedForces[347]       = Vec3(  3.1630812e+02,  -8.9974673e+02,   9.6573268e+01);
    expectedForces[348]       = Vec3( -5.5334313e+01,  -1.1529065e+03,  -2.1949997e+02);
    expectedForces[349]       = Vec3( -4.4904784e+02,   2.4036076e+02,   4.1328142e+02);
    expectedForces[350]       = Vec3(  1.0611611e+03,   4.1620143e+02,   9.1677657e+01);
    expectedForces[351]       = Vec3(  1.3489840e+03,   9.9500659e+02,  -4.5894902e+01);
    expectedForces[352]       = Vec3( -9.8051252e+01,  -9.0794586e+02,  -8.9918421e+02);
    expectedForces[353]       = Vec3( -3.5567408e+02,  -7.2914902e+01,   4.7977644e+01);
    expectedForces[354]       = Vec3( -1.5976501e+03,  -1.2202674e+03,   7.2159213e+02);
    expectedForces[355]       = Vec3(  1.7391266e+02,   1.0197773e+03,  -9.1284547e+02);
    expectedForces[356]       = Vec3(  7.8937287e+02,   1.1964969e+02,  -3.8520683e+02);
    expectedForces[357]       = Vec3(  2.8170878e+02,   1.0377979e+03,  -5.0609230e+02);
    expectedForces[358]       = Vec3( -4.0852118e+01,  -4.3087314e+02,   4.1855459e+01);
    expectedForces[359]       = Vec3( -3.2767902e+02,  -7.8083477e+02,   1.1111190e+03);
    expectedForces[360]       = Vec3( -1.0691030e+03,   3.1877408e+02,  -7.9684323e+02);
    expectedForces[361]       = Vec3(  7.7603235e+02,   3.6577850e+02,  -1.9093841e+02);
    expectedForces[362]       = Vec3( -4.7607360e+02,  -5.1710653e+02,   7.2740737e+02);
    expectedForces[363]       = Vec3(  9.3461900e+02,   7.9988609e+01,  -5.8055314e+02);
    expectedForces[364]       = Vec3( -9.8128918e+02,  -2.6706371e+02,  -3.5178135e+01);
    expectedForces[365]       = Vec3( -6.5196668e+02,   8.7618054e+02,   3.3040412e+02);
    expectedForces[366]       = Vec3(  5.5458970e+02,  -1.1281839e+03,  -8.2774754e+02);
    expectedForces[367]       = Vec3(  1.8491757e+02,  -5.1421593e+01,   4.7068191e+02);
    expectedForces[368]       = Vec3( -4.3945944e+02,   8.2740025e+02,  -2.1736033e+01);
    expectedForces[369]       = Vec3( -2.5609394e+02,  -6.7141305e+02,  -3.2964376e+02);
    expectedForces[370]       = Vec3(  1.4932730e+02,   2.2746635e+02,  -3.7606156e+01);
    expectedForces[371]       = Vec3(  3.0873674e+02,   5.9974800e+02,   4.2207331e+02);
    expectedForces[372]       = Vec3( -3.8828932e+02,  -2.1491002e+02,   1.5266506e+03);
    expectedForces[373]       = Vec3(  4.4495676e+02,   6.4708482e+02,  -9.2222368e+02);
    expectedForces[374]       = Vec3( -4.8420767e+01,  -4.3781484e+02,  -5.0107314e+02);
    expectedForces[375]       = Vec3(  5.7593719e+02,  -1.8140066e+03,  -4.0189721e+02);
    expectedForces[376]       = Vec3( -5.1276233e+02,   6.6981030e+02,   3.9050744e+02);
    expectedForces[377]       = Vec3(  4.0809391e+02,   1.1596412e+03,  -4.1325341e+02);
    expectedForces[378]       = Vec3(  1.4975560e+03,   1.7852942e+02,  -7.6514466e+02);
    expectedForces[379]       = Vec3( -2.2890988e+02,   2.6128742e+02,   3.8545036e+02);
    expectedForces[380]       = Vec3( -3.8899762e+02,  -2.5609958e+02,   2.0655882e+02);
    expectedForces[381]       = Vec3( -1.9500869e+02,  -1.0947633e+03,  -9.1786660e+02);
    expectedForces[382]       = Vec3(  8.6146884e+02,   2.5767444e+02,   3.6801549e+02);
    expectedForces[383]       = Vec3( -2.0008562e+02,   2.1549793e+02,   2.5175877e+02);
    expectedForces[384]       = Vec3( -5.6491749e+02,   5.4714989e+02,   3.1934114e+02);
    expectedForces[385]       = Vec3(  1.7665111e+02,  -4.5297277e+02,   2.2387580e+02);
    expectedForces[386]       = Vec3(  6.2697664e+02,   1.1271358e+02,   2.9498078e+00);
    expectedForces[387]       = Vec3( -2.1918090e+03,  -7.8914005e+01,   1.0632280e+03);
    expectedForces[388]       = Vec3(  7.5750278e+02,  -5.0217666e+02,  -1.3057335e+02);
    expectedForces[389]       = Vec3(  1.0621096e+03,   1.6022661e+01,  -1.1645539e+03);
    expectedForces[390]       = Vec3( -4.7808938e+02,  -1.2425496e+03,  -1.5543074e+02);
    expectedForces[391]       = Vec3( -3.4676860e+02,   8.5391303e+02,   3.5351618e+01);
    expectedForces[392]       = Vec3(  5.4017203e+02,   8.6467815e+01,   1.1898140e+02);
    expectedForces[393]       = Vec3( -1.8873828e+01,   2.2133074e+02,  -1.3378739e+03);
    expectedForces[394]       = Vec3(  7.5888012e+01,  -1.6517743e+01,   1.1817186e+02);
    expectedForces[395]       = Vec3(  1.1912944e+02,  -1.6083226e+02,   7.8733940e+02);
    expectedForces[396]       = Vec3( -1.2400005e+03,  -3.9434827e+02,   1.4071802e+02);
    expectedForces[397]       = Vec3(  9.6249284e+02,   1.1394516e+02,  -3.4977717e+02);
    expectedForces[398]       = Vec3(  2.4187486e+02,   1.3373651e+02,   6.5894105e+01);
    expectedForces[399]       = Vec3(  8.8628445e+02,   7.8083892e+01,  -1.0288922e+03);
    expectedForces[400]       = Vec3( -9.4613057e+02,  -1.9076053e+02,   6.7506442e+02);
    expectedForces[401]       = Vec3( -6.2746897e+02,   2.9376858e+02,   9.2767458e+02);
    expectedForces[402]       = Vec3(  5.5409175e+01,  -1.2583442e+03,   7.1728490e+02);
    expectedForces[403]       = Vec3( -2.9076856e+02,   7.5539656e+02,   7.0121763e+00);
    expectedForces[404]       = Vec3(  2.6220897e+02,  -9.5606102e+01,  -7.7725998e+02);
    expectedForces[405]       = Vec3( -8.6582014e+02,   9.5597761e+02,   1.5941783e+02);
    expectedForces[406]       = Vec3( -1.6910265e+02,  -7.2646192e+02,  -3.5476587e+02);
    expectedForces[407]       = Vec3(  9.7629076e+02,  -8.3969437e+01,   2.5523637e+02);
    expectedForces[408]       = Vec3( -4.9553396e+01,   6.3557270e+02,  -7.5908312e+02);
    expectedForces[409]       = Vec3(  3.3210227e+01,   1.5198340e+02,   7.0322192e+02);
    expectedForces[410]       = Vec3( -2.6532687e+01,  -4.1589199e+02,   4.2771258e+02);
    expectedForces[411]       = Vec3(  1.1412630e+03,  -2.7366656e+02,   9.8419548e+02);
    expectedForces[412]       = Vec3( -3.7239786e+02,   7.2244667e+01,  -9.9502150e+02);
    expectedForces[413]       = Vec3( -4.6476941e+02,  -6.1433607e+01,  -2.5459288e+02);
    expectedForces[414]       = Vec3(  9.7028883e+02,   5.5981848e+02,   8.4425262e+02);
    expectedForces[415]       = Vec3(  7.2811577e+01,   2.2872709e+02,  -1.3822700e+03);
    expectedForces[416]       = Vec3( -1.0105821e+03,  -2.8444698e+02,  -7.2506392e+02);
    expectedForces[417]       = Vec3(  1.0027865e+03,   6.0924816e+02,  -5.7818592e+01);
    expectedForces[418]       = Vec3( -3.4173955e+02,  -1.1932310e+02,  -2.9495449e+01);
    expectedForces[419]       = Vec3( -2.7281265e+02,  -1.8869212e+02,   1.9643932e+02);
    expectedForces[420]       = Vec3( -7.4036328e+02,  -4.8733524e+02,   1.5862094e+03);
    expectedForces[421]       = Vec3(  5.6387864e+02,   3.0991659e+02,  -8.6746725e+02);
    expectedForces[422]       = Vec3(  6.2954421e-01,   4.9665567e+02,  -1.0114913e+03);
    expectedForces[423]       = Vec3(  8.9668630e+02,  -1.0121499e+03,  -1.1520006e+03);
    expectedForces[424]       = Vec3( -2.5699741e+02,   3.3243520e+02,   7.6038873e+02);
    expectedForces[425]       = Vec3( -1.2755484e+03,  -2.7786159e+01,   3.0900583e+02);
    expectedForces[426]       = Vec3( -1.2587339e+03,  -8.6851333e+02,   1.6295957e+02);
    expectedForces[427]       = Vec3(  6.4923970e+02,   4.4600015e+02,   3.3033348e+02);
    expectedForces[428]       = Vec3(  5.8905637e+02,  -1.1589062e+02,  -1.8168184e+02);
    expectedForces[429]       = Vec3( -1.1223714e+03,   3.0322406e+01,   8.7272053e+02);
    expectedForces[430]       = Vec3(  3.0982625e+02,   3.1418220e+02,  -3.9301742e+02);
    expectedForces[431]       = Vec3(  2.4213933e+02,  -7.2382572e+02,  -1.0155346e+03);
    expectedForces[432]       = Vec3(  1.1479692e+03,  -1.6837721e+03,   1.0545407e+02);
    expectedForces[433]       = Vec3( -8.2496264e+02,   6.0540594e+02,   1.3979931e+02);
    expectedForces[434]       = Vec3( -6.8171514e+02,   4.0392791e+02,  -3.4712316e+02);
    expectedForces[435]       = Vec3( -1.5568889e+02,  -1.4652975e+03,   5.1518148e+01);
    expectedForces[436]       = Vec3(  5.6573856e+02,   3.5494119e+02,  -3.3796345e+02);
    expectedForces[437]       = Vec3( -1.3938679e+02,   4.2296152e+02,  -2.0539863e+02);
    expectedForces[438]       = Vec3(  8.9203493e+01,   2.8764597e+02,  -7.3273496e+02);
    expectedForces[439]       = Vec3(  4.3114292e+02,   9.4778082e+01,   1.8894761e+02);
    expectedForces[440]       = Vec3( -2.2798193e+01,  -4.5460891e+01,   9.8310963e+01);
    expectedForces[441]       = Vec3( -2.1793579e+02,  -1.0807542e+03,  -2.3470465e+01);
    expectedForces[442]       = Vec3(  1.3881360e+02,   3.9806903e+02,   4.1089741e+01);
    expectedForces[443]       = Vec3( -4.5604974e+02,   4.8515999e+02,  -6.6025174e+02);
    expectedForces[444]       = Vec3(  3.6437700e+02,  -8.1622511e+02,  -9.6454258e+02);
    expectedForces[445]       = Vec3(  3.1998690e+02,  -3.3342314e+02,   1.2441763e+03);
    expectedForces[446]       = Vec3(  3.0822409e+01,   2.0596612e+02,   2.0937122e+02);
    expectedForces[447]       = Vec3( -8.4938718e+02,  -1.1366483e+03,  -1.3638049e+02);
    expectedForces[448]       = Vec3( -5.5232451e+01,   4.7335097e+02,  -5.4433565e+02);
    expectedForces[449]       = Vec3(  6.5928852e+02,   1.7488804e+02,   6.7879378e+02);
    expectedForces[450]       = Vec3( -9.4572079e+02,   1.9162420e+02,   4.7935043e+02);
    expectedForces[451]       = Vec3(  2.5029547e+02,   3.6606767e+01,  -9.3423713e+01);
    expectedForces[452]       = Vec3(  1.6543841e+02,  -1.1892943e+02,  -5.8737050e+02);
    expectedForces[453]       = Vec3( -9.8495219e+02,  -1.6882428e+03,  -4.0576035e+02);
    expectedForces[454]       = Vec3(  7.9182557e+02,   5.3997699e+02,   3.3231603e+02);
    expectedForces[455]       = Vec3(  1.8277090e+01,   1.4623150e+03,  -9.4985721e+01);
    expectedForces[456]       = Vec3( -1.0471948e+03,   4.6746395e+02,  -1.5020000e+03);
    expectedForces[457]       = Vec3(  7.2055316e+02,  -5.4467216e+02,   4.0459421e+02);
    expectedForces[458]       = Vec3(  7.5475531e+01,   4.9532870e+02,   1.3002474e+03);
    expectedForces[459]       = Vec3(  6.3589773e+02,  -2.1402397e+02,  -1.4147509e+03);
    expectedForces[460]       = Vec3(  2.8391861e+02,  -8.0142885e+01,   5.6445225e+02);
    expectedForces[461]       = Vec3( -8.2442674e+02,   5.3520687e+02,   1.1406667e+03);
    expectedForces[462]       = Vec3( -5.7988771e+02,  -3.3512887e+02,   4.5461752e+02);
    expectedForces[463]       = Vec3(  4.5746059e+02,   3.1029926e+02,  -5.0060176e+02);
    expectedForces[464]       = Vec3(  7.7747136e+02,   4.2131732e+02,   5.5836239e+02);
    expectedForces[465]       = Vec3( -1.2315491e+03,   9.4066088e+02,   9.6145313e+02);
    expectedForces[466]       = Vec3(  7.9043841e+02,  -2.2172613e+02,  -2.6508587e+02);
    expectedForces[467]       = Vec3(  7.8197894e+02,  -3.1383822e+02,   2.5444013e+02);
    expectedForces[468]       = Vec3( -7.1139498e+01,  -7.1203189e+01,  -4.6845544e+02);
    expectedForces[469]       = Vec3( -8.9774100e+01,   1.4389287e+02,   9.8957451e+01);
    expectedForces[470]       = Vec3(  3.0397922e+02,  -2.8247850e+01,   4.4896034e+02);
    expectedForces[471]       = Vec3( -9.7808367e+02,   1.0170553e+03,   8.1594649e+02);
    expectedForces[472]       = Vec3(  1.8933676e+02,  -8.8053046e+02,   6.8784042e+01);
    expectedForces[473]       = Vec3(  5.1636668e+02,   1.2970985e+02,  -8.9380858e+02);
    expectedForces[474]       = Vec3( -5.8549608e+02,  -1.8351156e+02,  -5.8043066e+02);
    expectedForces[475]       = Vec3(  2.7877663e+02,  -3.6576715e+02,   3.5497095e+02);
    expectedForces[476]       = Vec3(  3.0642370e+02,   5.3438407e+02,   1.9409584e+02);
    expectedForces[477]       = Vec3( -4.8523805e+02,   1.2253320e+03,   9.6414379e+02);
    expectedForces[478]       = Vec3( -1.5213013e+02,  -2.0017205e+02,  -6.6602643e+02);
    expectedForces[479]       = Vec3(  4.6435225e+02,  -4.1945066e+02,  -8.6774270e+01);
    expectedForces[480]       = Vec3(  2.4946889e+02,  -1.0456281e+03,   4.7404434e+02);
    expectedForces[481]       = Vec3( -3.0813505e+01,   2.8598938e+02,   7.7374350e+01);
    expectedForces[482]       = Vec3( -2.4538851e+02,   3.8306987e+02,  -5.0235520e+02);
    expectedForces[483]       = Vec3(  4.0348805e+01,   1.3050360e+03,  -6.8725580e+02);
    expectedForces[484]       = Vec3( -1.1095514e+02,  -2.7711572e+02,   2.6929959e+02);
    expectedForces[485]       = Vec3( -6.8071081e+01,  -3.4398150e+02,   2.5209743e+02);
    expectedForces[486]       = Vec3(  2.1534786e+03,   1.3249493e+01,   7.4171165e+01);
    expectedForces[487]       = Vec3( -7.2618755e+02,  -3.2914219e+02,   2.5917332e+02);
    expectedForces[488]       = Vec3( -1.0231949e+03,   7.2426062e+02,   1.9111862e+02);
    expectedForces[489]       = Vec3(  1.1408866e+03,   6.7858353e+02,  -2.1410916e+02);
    expectedForces[490]       = Vec3( -4.5559047e+02,  -1.3597950e+02,  -3.2284091e+02);
    expectedForces[491]       = Vec3( -8.5558133e+02,  -7.1748324e+01,   5.3332261e+02);
    expectedForces[492]       = Vec3( -7.1393886e+02,  -1.1275222e+03,  -6.2147584e+02);
    expectedForces[493]       = Vec3(  4.4029614e+02,   1.0518224e+02,   4.6519788e+02);
    expectedForces[494]       = Vec3( -3.5770378e+02,   1.0311834e+03,   2.2141802e+02);
    expectedForces[495]       = Vec3( -2.6023787e+02,   1.0070248e+03,  -1.1113552e+03);
    expectedForces[496]       = Vec3( -3.4950242e+02,   2.8418846e+01,   1.2865161e+03);
    expectedForces[497]       = Vec3(  1.5030225e+02,  -7.7257505e+02,   8.7232628e+02);
    expectedForces[498]       = Vec3( -1.7416695e+02,   9.9945569e+02,  -1.6742483e+02);
    expectedForces[499]       = Vec3(  2.0837481e+02,  -1.9388709e+02,   9.7409815e+01);
    expectedForces[500]       = Vec3( -1.0129845e+02,  -6.8652976e+02,   8.4930186e+02);
    expectedForces[501]       = Vec3(  1.2276099e+03,   1.6531986e+02,   3.6342506e+02);
    expectedForces[502]       = Vec3( -1.9998077e+02,  -2.6164759e+02,  -1.6601933e+02);
    expectedForces[503]       = Vec3( -4.0356327e+02,   4.0549152e+02,   2.9710052e+02);
    expectedForces[504]       = Vec3( -7.3045150e+02,  -1.6677706e+03,  -9.8671765e+02);
    expectedForces[505]       = Vec3( -2.1441879e+02,   1.4962605e+03,   5.6350856e+02);
    expectedForces[506]       = Vec3(  7.1102593e+02,   5.9407971e+02,   9.5452239e+02);
    expectedForces[507]       = Vec3(  3.2291041e+02,   1.3943261e+03,  -1.8477432e+03);
    expectedForces[508]       = Vec3( -1.7778104e+02,   2.8034575e+02,   1.1385089e+03);
    expectedForces[509]       = Vec3( -2.4943297e+02,  -9.0965773e+02,   5.2737490e+02);
    expectedForces[510]       = Vec3( -3.8318508e+02,   1.3999092e+03,  -5.9960180e+02);
    expectedForces[511]       = Vec3(  2.6930065e+01,  -3.6603454e+02,   1.0460530e+03);
    expectedForces[512]       = Vec3( -1.8964668e+02,  -1.1663184e+03,  -2.8438291e+02);
    expectedForces[513]       = Vec3( -9.4679273e+02,   8.7408257e+01,  -2.5740702e+02);
    expectedForces[514]       = Vec3(  1.6558849e+02,   3.3561310e+02,   3.8532718e+02);
    expectedForces[515]       = Vec3(  4.0403504e+02,  -8.5545173e+02,   4.0647038e+02);
    expectedForces[516]       = Vec3(  8.3825526e+02,   5.0343638e+01,  -3.9171626e+02);
    expectedForces[517]       = Vec3( -4.0069871e+02,  -6.4140295e+02,   7.5218861e+01);
    expectedForces[518]       = Vec3( -6.5213072e+02,   4.1689581e+02,  -8.9280241e+01);
    expectedForces[519]       = Vec3(  4.0832996e+02,  -3.4272605e+02,  -9.2740030e+02);
    expectedForces[520]       = Vec3( -2.1128040e+01,  -7.3467005e+02,   6.5449252e+02);
    expectedForces[521]       = Vec3( -2.8042421e+02,   3.3193211e+02,   5.1272473e+02);
    expectedForces[522]       = Vec3( -1.2057322e+03,   3.3276609e+02,   6.8838073e+02);
    expectedForces[523]       = Vec3(  8.9270663e+02,   1.3914748e+02,  -8.7213182e+02);
    expectedForces[524]       = Vec3(  8.8063936e+02,  -2.4952564e+02,  -8.3942753e+01);
    expectedForces[525]       = Vec3(  3.6527919e+02,   3.3758036e+02,  -1.3449147e+03);
    expectedForces[526]       = Vec3(  7.1836993e+01,  -3.9380031e+01,   4.7169292e+02);
    expectedForces[527]       = Vec3( -7.2114939e+01,  -1.3679457e+02,   3.9080695e+02);
    expectedForces[528]       = Vec3(  1.1465645e+02,  -4.4994719e+02,  -3.4180224e+02);
    expectedForces[529]       = Vec3( -6.0950292e+01,   4.2825716e+02,  -1.4233246e+02);
    expectedForces[530]       = Vec3( -1.9112828e+02,  -2.8511065e+00,   2.6894050e+02);
    expectedForces[531]       = Vec3(  4.4955085e+02,   2.0459389e+02,   8.9486972e+02);
    expectedForces[532]       = Vec3(  6.6339624e+01,  -4.2734850e+01,  -4.8004048e+02);
    expectedForces[533]       = Vec3( -5.0521801e+02,  -1.6330265e+02,  -6.2374206e+02);
    expectedForces[534]       = Vec3( -1.1019306e+03,  -4.5760347e+02,   4.6134602e+01);
    expectedForces[535]       = Vec3(  4.9278700e+02,   8.0495847e+02,  -5.7470233e+01);
    expectedForces[536]       = Vec3(  9.9361557e+02,  -2.4453172e+02,  -3.7696369e+02);
    expectedForces[537]       = Vec3( -1.1168499e+03,  -2.1108925e+02,  -4.1233970e+02);
    expectedForces[538]       = Vec3(  5.6559058e+02,  -2.8153614e+02,   3.0539700e+02);
    expectedForces[539]       = Vec3(  4.0805303e+02,   2.8796083e+02,  -3.0009135e+01);
    expectedForces[540]       = Vec3( -8.5249075e+02,  -1.0560038e+03,  -6.5111795e+02);
    expectedForces[541]       = Vec3(  8.4653565e+02,   2.0615416e+02,  -1.5085906e+02);
    expectedForces[542]       = Vec3( -2.8200251e+02,   6.6950966e+02,   5.9661450e+02);
    expectedForces[543]       = Vec3( -9.3339994e+01,   8.7084190e+02,  -7.8375352e+02);
    expectedForces[544]       = Vec3(  1.4187243e+02,   1.2829809e+02,   8.0626841e+02);
    expectedForces[545]       = Vec3(  1.5403073e+02,  -4.5932125e+02,  -8.3051878e+00);
    expectedForces[546]       = Vec3(  1.3919197e+03,   5.4728737e+02,  -1.0548033e+03);
    expectedForces[547]       = Vec3( -2.7215065e+02,   1.5112718e+02,   6.7366542e+02);
    expectedForces[548]       = Vec3( -7.8677637e+02,  -2.6895175e+02,   4.3374996e+02);
    expectedForces[549]       = Vec3(  6.6014864e+02,  -1.1474704e+03,  -3.7199889e+02);
    expectedForces[550]       = Vec3( -5.8741315e+02,   2.6286109e+02,   1.9036264e+02);
    expectedForces[551]       = Vec3(  1.4325969e+02,   7.8928381e+02,   3.6304214e+02);
    expectedForces[552]       = Vec3( -9.6661866e+02,  -1.0462801e+03,  -6.3261994e+02);
    expectedForces[553]       = Vec3(  5.5171916e+02,   2.8672800e+02,   2.7240662e+02);
    expectedForces[554]       = Vec3(  1.1625052e+03,   7.7864110e+02,  -5.2588322e+02);
    expectedForces[555]       = Vec3(  1.7042269e+03,  -5.8973225e+02,   1.4814500e+02);
    expectedForces[556]       = Vec3( -7.7324022e+02,   2.9167426e+02,  -1.0262563e+02);
    expectedForces[557]       = Vec3( -5.4655014e+02,  -8.4339654e+01,   9.7375629e+02);
    expectedForces[558]       = Vec3( -6.9618127e+02,  -1.6878530e+02,   7.1078501e+02);
    expectedForces[559]       = Vec3(  6.9623995e+02,  -3.7038954e+01,   2.2671528e+02);
    expectedForces[560]       = Vec3( -7.1025694e+01,  -1.8643963e+02,  -8.0181026e+02);
    expectedForces[561]       = Vec3( -1.4598467e+02,   1.7170707e+03,  -4.3305456e+02);
    expectedForces[562]       = Vec3(  7.9582463e+02,  -5.9551245e+02,   8.3485625e+01);
    expectedForces[563]       = Vec3( -4.6812643e+02,  -3.4598874e+02,  -3.6753356e+01);
    expectedForces[564]       = Vec3( -1.8586268e+02,   3.1377143e+02,  -1.6616022e+03);
    expectedForces[565]       = Vec3(  7.4400496e+02,   5.5058789e+02,   7.5374599e+02);
    expectedForces[566]       = Vec3( -6.3347529e+02,  -3.1407006e+02,   8.2964691e+02);
    expectedForces[567]       = Vec3( -4.2499528e+02,   9.4575975e+02,  -6.6700103e+02);
    expectedForces[568]       = Vec3( -9.2637696e+02,  -6.6301622e+02,   3.9913705e+02);
    expectedForces[569]       = Vec3(  5.0588430e+02,  -8.1942697e+01,   4.2128365e+02);
    expectedForces[570]       = Vec3(  2.7165158e+02,   2.7387902e+02,   1.0850277e+03);
    expectedForces[571]       = Vec3( -6.8391471e+00,  -1.2886307e+02,  -5.3501450e+02);
    expectedForces[572]       = Vec3( -2.9936240e+02,  -2.2085475e+02,  -6.9670281e+01);
    expectedForces[573]       = Vec3(  1.6853422e+02,   4.8410887e+02,   1.0517406e+03);
    expectedForces[574]       = Vec3( -4.9287521e+02,   1.6557733e+01,  -3.4008916e+02);
    expectedForces[575]       = Vec3(  1.1132397e+02,  -9.8734894e+02,  -1.1154151e+02);
    expectedForces[576]       = Vec3(  3.8645582e+02,   4.3709361e+02,   1.4861132e+03);
    expectedForces[577]       = Vec3( -1.0169967e+03,  -3.5279110e+02,  -5.4860861e+02);
    expectedForces[578]       = Vec3(  2.8875246e+02,   5.3726711e+01,  -5.9395936e+02);
    expectedForces[579]       = Vec3(  1.5169994e+03,   8.7535209e+01,   1.0103159e+03);
    expectedForces[580]       = Vec3( -1.9439968e+02,  -7.9639200e+02,  -1.0439069e+03);
    expectedForces[581]       = Vec3( -8.7374025e+02,   2.0597360e+02,   1.0218319e+02);
    expectedForces[582]       = Vec3(  6.4275846e+02,  -6.2671529e+02,   8.9274592e+02);
    expectedForces[583]       = Vec3( -5.4449356e+01,   7.2615928e+02,  -9.1355057e+02);
    expectedForces[584]       = Vec3( -2.7242606e+02,   3.9376960e+02,  -4.7692581e+01);
    expectedForces[585]       = Vec3(  2.6578654e+02,  -6.7385441e+02,  -1.5306127e+03);
    expectedForces[586]       = Vec3(  5.1049761e+02,   7.3430855e+02,   6.1178168e+02);
    expectedForces[587]       = Vec3( -3.3593193e+02,   3.5079666e+01,   9.7444298e+02);
    expectedForces[588]       = Vec3( -7.8227654e+01,  -1.1050481e+03,   1.4161574e+03);
    expectedForces[589]       = Vec3( -7.6753407e+02,   9.7291633e+02,   1.9943116e+02);
    expectedForces[590]       = Vec3(  2.9274184e+02,  -2.5706985e+01,  -5.8944290e+02);
    expectedForces[591]       = Vec3( -1.2502829e+03,   2.2304820e+02,  -8.8573767e+02);
    expectedForces[592]       = Vec3(  1.2218682e+02,  -8.2110835e+02,   3.4276058e+02);
    expectedForces[593]       = Vec3(  8.9529566e+02,  -1.8417573e+02,   8.4717919e+02);
    expectedForces[594]       = Vec3( -4.6475772e+01,  -1.5493620e+03,   4.8004365e+02);
    expectedForces[595]       = Vec3(  6.7104033e+02,   4.7267915e+02,  -6.2630817e+02);
    expectedForces[596]       = Vec3( -2.2771561e+02,   2.7463551e+02,  -2.9704610e+00);
    expectedForces[597]       = Vec3(  1.7757815e+02,  -1.8906761e+03,  -1.6656886e+02);
    expectedForces[598]       = Vec3( -7.4745588e+01,   1.0353962e+03,   4.2924498e+02);
    expectedForces[599]       = Vec3( -5.5648197e+02,   7.7838572e+02,  -4.1677346e+02);
    expectedForces[600]       = Vec3( -5.2972489e+02,   1.1967117e+02,   1.0565388e+03);
    expectedForces[601]       = Vec3(  7.8379665e+02,  -4.1803420e+02,   2.0308387e+02);
    expectedForces[602]       = Vec3(  1.2705528e+02,   1.1663435e+02,  -7.6514830e+02);
    expectedForces[603]       = Vec3( -1.5203094e+03,  -8.5825144e+01,  -2.3039380e+02);
    expectedForces[604]       = Vec3(  5.2063739e+02,  -1.7259834e+02,  -2.9854160e+01);
    expectedForces[605]       = Vec3(  4.7448961e+02,   2.8282389e+02,   3.4257463e+01);
    expectedForces[606]       = Vec3(  1.4519531e+03,  -1.2031444e+03,  -6.7084095e+02);
    expectedForces[607]       = Vec3( -7.5308994e+02,  -3.2117501e+02,   5.7700017e+02);
    expectedForces[608]       = Vec3( -1.7643179e+02,   8.5655779e+02,  -2.6285482e+02);
    expectedForces[609]       = Vec3( -1.3868959e+03,  -1.6195505e+01,  -1.3783528e+03);
    expectedForces[610]       = Vec3(  2.8034476e+02,  -4.2360778e+02,   4.5904202e+02);
    expectedForces[611]       = Vec3(  1.0938532e+03,   6.2455629e+02,  -1.0108636e+02);
    expectedForces[612]       = Vec3(  9.1439974e+01,  -7.2813810e+02,   1.0197075e+03);
    expectedForces[613]       = Vec3(  1.4048294e+02,   2.5667744e+02,  -5.1914694e+02);
    expectedForces[614]       = Vec3( -2.5604629e+01,   2.9915218e+01,  -5.1582702e+02);
    expectedForces[615]       = Vec3(  2.0331673e+02,   1.5538650e+02,  -1.1836865e+03);
    expectedForces[616]       = Vec3( -2.3788957e+02,  -1.1872269e+02,   9.7052698e+02);
    expectedForces[617]       = Vec3( -3.2419551e+01,  -1.8828745e+02,   4.3037423e+02);
    expectedForces[618]       = Vec3( -4.8585762e+02,   1.4114027e+02,   1.6021114e+02);
    expectedForces[619]       = Vec3(  2.8744134e+02,  -5.7197623e+02,  -1.8323508e+02);
    expectedForces[620]       = Vec3(  2.3652930e+02,  -1.5448350e+01,   2.1721186e+01);
    expectedForces[621]       = Vec3(  1.1684581e+03,   2.2790484e+02,   2.7798348e+02);
    expectedForces[622]       = Vec3( -4.6823976e+02,  -6.1804544e+02,   2.7193441e+02);
    expectedForces[623]       = Vec3( -5.8173089e+02,   4.4573870e+02,  -6.0699430e+02);
    expectedForces[624]       = Vec3(  7.9759250e+02,   1.8480473e+03,  -5.1434277e+02);
    expectedForces[625]       = Vec3( -1.0945975e+03,  -6.0402685e+02,   1.5843017e+02);
    expectedForces[626]       = Vec3(  1.2334392e+02,  -9.5602838e+02,   7.8137824e+02);
    expectedForces[627]       = Vec3( -8.5863208e+02,  -5.7413454e+02,   1.7014217e+02);
    expectedForces[628]       = Vec3( -1.2232569e+02,   6.0505774e+02,  -1.6806693e+02);
    expectedForces[629]       = Vec3(  6.1488135e+02,   4.1603705e+02,   2.1121948e+02);
    expectedForces[630]       = Vec3(  2.1858899e+02,   1.2961013e+03,   4.8534041e+02);
    expectedForces[631]       = Vec3(  2.1292701e+02,  -4.1049165e+02,   3.9057142e+00);
    expectedForces[632]       = Vec3( -6.3342469e+02,  -8.7529004e+02,   6.1928550e+01);
    expectedForces[633]       = Vec3( -7.3207079e+02,  -1.3874807e+03,   8.6303032e+02);
    expectedForces[634]       = Vec3(  6.6660351e+02,   9.5813895e+02,  -1.4555416e+02);
    expectedForces[635]       = Vec3(  4.1800653e+02,   2.4677420e+02,  -1.1424272e+03);
    expectedForces[636]       = Vec3(  8.3570430e+02,  -1.5342969e+03,  -1.2510269e+03);
    expectedForces[637]       = Vec3( -3.0635229e+02,   1.2465088e+03,  -1.9679529e+01);
    expectedForces[638]       = Vec3( -6.5317178e+02,  -1.6432929e+02,   9.2427724e+02);
    expectedForces[639]       = Vec3(  1.5665680e+03,   9.8040004e+01,   1.0786497e+02);
    expectedForces[640]       = Vec3( -4.6590380e+02,  -3.7603795e+02,  -9.4955273e+02);
    expectedForces[641]       = Vec3( -8.7572108e+02,   5.1750570e+02,  -1.2268062e+02);
    expectedForces[642]       = Vec3(  1.1575067e+03,  -2.2857204e+02,  -2.3816900e+02);
    expectedForces[643]       = Vec3( -2.6229644e+02,   3.1069110e+02,  -1.2098536e+02);
    expectedForces[644]       = Vec3( -1.8195659e+02,  -3.8984698e+02,   6.4622752e+02);
    expectedForces[645]       = Vec3(  1.0756277e+03,   2.7459106e+01,   1.0850508e+03);
    expectedForces[646]       = Vec3( -6.4620310e+02,   2.5885783e+02,  -2.0567224e+02);
    expectedForces[647]       = Vec3( -4.7388806e+02,  -5.5561844e+02,  -8.5019295e+02);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);

}

// test querying particle induced dipoles

static void testParticleInducedDipoles() {
    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    getForcesEnergyMultipoleAmmonia(context, forces, energy);
    std::vector<Vec3> dipole;
    amoebaMultipoleForce->getInducedDipoles(context, dipole);
    
    // Compare to values calculated by TINKER.
    
    std::vector<Vec3> expectedDipole(numberOfParticles);
    expectedDipole[0] = Vec3(0.0031710288, 9.3687453e-7, -0.0006919963);
    expectedDipole[1] = Vec3(8.0279737504e-5, -0.000279376, 4.778060103e-5);
    expectedDipole[2] = Vec3(0.000079322, 0.0002789804, 4.8696656126e-5);
    expectedDipole[3] = Vec3(-0.0001407394, 1.540638116e-6, -0.0007077775);
    expectedDipole[4] = Vec3(0.0019564439, -1.0409717e-7, 0.0007332188);
    expectedDipole[5] = Vec3(0.0008213891, -0.0007749618, -0.0003883865);
    expectedDipole[6] = Vec3(0.0046133992, -7.2868019e-7, 0.0002500622);
    expectedDipole[7] = Vec3(0.0008204731, 0.0007772727, -0.0003856176);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(expectedDipole[i], dipole[i], 1e-4);
}

// test querying particle lab frame permanent dipoles

static void testParticleLabFramePermanentDipoles() {
    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    getForcesEnergyMultipoleAmmonia(context, forces, energy);
    std::vector<Vec3> dipole;
    amoebaMultipoleForce->getLabFramePermanentDipoles(context, dipole);
    
    // Compare to values calculated by TINKER.
    
    std::vector<Vec3> expectedDipole(numberOfParticles);
    expectedDipole[0] = Vec3(0.00876454250, -2.04310718E-06, -0.00227593519);
    expectedDipole[1] = Vec3(0.000780382180, -0.00432882849, 0.00236926761);
    expectedDipole[2] = Vec3(0.000801345883, 0.00431830946, 0.00238143437);
    expectedDipole[3] = Vec3(-0.00109746996, 1.16087953e-5, -0.00487407492);
    expectedDipole[4] = Vec3(0.00203814102, -2.26554196e-5, 0.00882284298);
    expectedDipole[5] = Vec3(-0.00239443187, 0.00432388648, 0.000729303209);
    expectedDipole[6] = Vec3(0.00491086743, 2.86430963e-6, -0.000918996348);
    expectedDipole[7] = Vec3(-0.00239301946, -0.00432743976, 0.000712674115);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(expectedDipole[i], dipole[i], 1e-4);
}

// test querying particle total dipoles (fixed + induced)

static void testParticleTotalDipoles() {
    int numberOfParticles     = 8;
    int inputPmeGridDimension = 0;
    double cutoff             = 9000000.0;
    std::vector<Vec3> forces;
    double energy;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Mutual, 
                                             cutoff, inputPmeGridDimension);
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    getForcesEnergyMultipoleAmmonia(context, forces, energy);
    std::vector<Vec3> dipole;
    amoebaMultipoleForce->getTotalDipoles(context, dipole);
    
    // Compare to values calculated by TINKER.
    
    std::vector<Vec3> expectedDipole(numberOfParticles);
    expectedDipole[0] = Vec3(0.0119356307, -1.11302433e-6, -0.00296793872);
    expectedDipole[1] = Vec3(8.60636211e-4, -0.00460821816, 0.00241705344);
    expectedDipole[2] = Vec3(8.80646403e-4, 0.00459728769, 0.00243013245);
    expectedDipole[3] = Vec3(-0.00123822377, 1.31555550e-5, -0.00558185336);
    expectedDipole[4] = Vec3(0.00399455556, -2.27511931e-5, 0.00955607952);
    expectedDipole[5] = Vec3(-0.00157302682, 0.00354892386, 3.40921137e-4);
    expectedDipole[6] = Vec3(0.00952428069, 2.14171505e-6, -6.68945865e-4);
    expectedDipole[7] = Vec3(-0.00157252460, -0.00355015528, 3.27055162e-4);
    for (int i = 0; i < numberOfParticles; i++)
        ASSERT_EQUAL_VEC(expectedDipole[i], dipole[i], 1e-4);
}



// test computation of system multipole moments

static void testSystemMultipoleMoments() {

    std::string testName      = "testSystemMultipoleMoments";
    
    int numberOfParticles     = 648;
    int inputPmeGridDimension = 24;
    double cutoff             = 0.70;

    std::vector<Vec3> forces;
    double energy;
    std::vector<double> outputMultipoleMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;

    setupAndGetForcesEnergyMultipoleLargeWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                               cutoff, inputPmeGridDimension, testName,
                                               forces, energy, outputMultipoleMoments, inputGrid, outputGridPotential);

    std::vector<double> tinkerMoments(4);

    tinkerMoments[0]  =   0.0000000e+00;
    tinkerMoments[1]  =  -9.1118361e+00;
    tinkerMoments[2]  =   2.8371876e+01;
    tinkerMoments[3]  =   5.1518898e+01;
//    tinkerMoments[4]  =  -1.0768808e-01; // Quadrupole moments are not uniquely defined when using periodic boundary conditions
//    tinkerMoments[5]  =  -9.0458124e-01;
//    tinkerMoments[6]  =   1.8460385e+00;
//    tinkerMoments[7]  =  -9.0458124e-01;
//    tinkerMoments[8]  =   6.4395591e-02;
//    tinkerMoments[9]  =   1.6692567e-01;
//    tinkerMoments[10] =   1.8460385e-00;
//    tinkerMoments[11] =   1.6692567e-01;
//    tinkerMoments[12] =   4.3292490e-02;

    double tolerance = 1.0e-04;
    for (unsigned int ii = 0; ii < tinkerMoments.size(); ii++) {
        double difference;
        if (fabs(tinkerMoments[ii]) > 0.0) {
            difference = fabs(outputMultipoleMoments[ii] - tinkerMoments[ii])/fabs(tinkerMoments[ii]);
        } else {
            difference = fabs(outputMultipoleMoments[ii] - tinkerMoments[ii]);
        }
        if (difference > tolerance) {
            std::stringstream details;
            details << testName << "Multipole moment " << ii << " does not agree w/ TINKER computed moments: OpenMM=" << outputMultipoleMoments[ii];
            details << " TINKER=" <<  tinkerMoments[ii]  << " difference=" << difference;
            throwException(__FILE__, __LINE__, details.str());
        }

    }

}

// test computation of multipole potential on a grid 

static void testMultipoleGridPotential() {

    std::string testName      = "testMultipoleGridPotential";
    
    int numberOfParticles     = 648;
    int inputPmeGridDimension = 24;
    double cutoff             = 0.70;

    std::vector<Vec3> forces;
    double energy;

    // initialize grid

    int gridSize  = 27;
    std::vector<Vec3> inputGrid(gridSize);

    inputGrid[0]  = Vec3(-1.0318535e+00,  -1.0224502e+00,  -1.0202836e+00);
    inputGrid[1]  = Vec3(-1.0318535e+00,  -1.0224502e+00,  -3.4032700e-02);
    inputGrid[2]  = Vec3(-1.0318535e+00,  -1.0224502e+00,   9.5221820e-01);
    inputGrid[3]  = Vec3(-1.0318535e+00,  -3.0961200e-02,  -1.0202836e+00);
    inputGrid[4]  = Vec3(-1.0318535e+00,  -3.0961200e-02,  -3.4032700e-02);
    inputGrid[5]  = Vec3(-1.0318535e+00,  -3.0961200e-02,   9.5221820e-01);
    inputGrid[6]  = Vec3(-1.0318535e+00,   9.6052780e-01,  -1.0202836e+00);
    inputGrid[7]  = Vec3(-1.0318535e+00,   9.6052780e-01,  -3.4032700e-02);
    inputGrid[8]  = Vec3(-1.0318535e+00,   9.6052780e-01,   9.5221820e-01);
    inputGrid[9]  = Vec3(-3.3969000e-02,  -1.0224502e+00,  -1.0202836e+00);
    inputGrid[10] = Vec3(-3.3969000e-02,  -1.0224502e+00,  -3.4032700e-02);
    inputGrid[11] = Vec3(-3.3969000e-02,  -1.0224502e+00,   9.5221820e-01);
    inputGrid[12] = Vec3(-3.3969000e-02,  -3.0961200e-02,  -1.0202836e+00);
    inputGrid[13] = Vec3(-3.3969000e-02,  -3.0961200e-02,  -3.4032700e-02);
    inputGrid[14] = Vec3(-3.3969000e-02,  -3.0961200e-02,   9.5221820e-01);
    inputGrid[15] = Vec3(-3.3969000e-02,   9.6052780e-01,  -1.0202836e+00);
    inputGrid[16] = Vec3(-3.3969000e-02,   9.6052780e-01,  -3.4032700e-02);
    inputGrid[17] = Vec3(-3.3969000e-02,   9.6052780e-01,   9.5221820e-01);
    inputGrid[18] = Vec3( 9.6391550e-01,  -1.0224502e+00,  -1.0202836e+00);
    inputGrid[19] = Vec3( 9.6391550e-01,  -1.0224502e+00,  -3.4032700e-02);
    inputGrid[20] = Vec3( 9.6391550e-01,  -1.0224502e+00,   9.5221820e-01);
    inputGrid[21] = Vec3( 9.6391550e-01,  -3.0961200e-02,  -1.0202836e+00);
    inputGrid[22] = Vec3( 9.6391550e-01,  -3.0961200e-02,  -3.4032700e-02);
    inputGrid[23] = Vec3( 9.6391550e-01,  -3.0961200e-02,   9.5221820e-01);
    inputGrid[24] = Vec3( 9.6391550e-01,   9.6052780e-01,  -1.0202836e+00);
    inputGrid[25] = Vec3( 9.6391550e-01,   9.6052780e-01,  -3.4032700e-02);
    inputGrid[26] = Vec3( 9.6391550e-01,   9.6052780e-01,   9.5221820e-01);

    std::vector<double> outputGridPotential;
    std::vector< double > outputMultipoleMoments; 

    setupAndGetForcesEnergyMultipoleLargeWater(AmoebaMultipoleForce::PME, AmoebaMultipoleForce::Mutual, 
                                                cutoff, inputPmeGridDimension, testName, forces, energy,
                                                outputMultipoleMoments, inputGrid, outputGridPotential);

    // TINKER computed grid values

    std::vector<double> tinkerGridPotential(gridSize);

    tinkerGridPotential[0]  =   7.5696473e+01;
    tinkerGridPotential[1]  =  -1.5382673e+01;
    tinkerGridPotential[2]  =   6.9487135e+01;
    tinkerGridPotential[3]  =   1.5459661e+03;
    tinkerGridPotential[4]  =  -4.5138366e+02;
    tinkerGridPotential[5]  =   4.6553857e+01;
    tinkerGridPotential[6]  =   1.8486978e+01;
    tinkerGridPotential[7]  =  -5.9235079e+01;
    tinkerGridPotential[8]  =   1.0448125e+01;
    tinkerGridPotential[9]  =   1.9529845e+00;
    tinkerGridPotential[10] =   5.3587575e+00;
    tinkerGridPotential[11] =  -3.3794408e+01;
    tinkerGridPotential[12] =   1.5197977e+01;
    tinkerGridPotential[13] =   4.3613381e+02;
    tinkerGridPotential[14] =   1.6620315e+01;
    tinkerGridPotential[15] =   5.3724507e+01;
    tinkerGridPotential[16] =   6.9026934e+02;
    tinkerGridPotential[17] =  -1.2708548e+01;
    tinkerGridPotential[18] =  -9.8826500e+02;
    tinkerGridPotential[19] =   3.2952143e+01;
    tinkerGridPotential[20] =   9.1165100e+02;
    tinkerGridPotential[21] =   4.0844598e+01;
    tinkerGridPotential[22] =  -1.3509483e+01;
    tinkerGridPotential[23] =   1.4499749e+00;
    tinkerGridPotential[24] =   1.5011664e+02;
    tinkerGridPotential[25] =  -2.8581208e+02;
    tinkerGridPotential[26] =   1.3960144e+02;

    double tolerance = 4.0e-04;
    for (unsigned int ii = 0; ii < gridSize; ii++) {
        double difference = fabs((outputGridPotential[ii] - tinkerGridPotential[ii])/tinkerGridPotential[ii]);
        if (difference > tolerance) {
            std::stringstream details;
            details << testName << " potential for grid point " << ii << " does not agree w/ TINKER computed value: OpenMM=" << outputGridPotential[ii];
            details << " TINKER=" <<  tinkerGridPotential[ii]  << " difference=" << difference;
            throwException(__FILE__, __LINE__, details.str());
        }
    }

}

/**
 * Test whether the alternate kernels with no quadrupoles work correctly.
 */
static void testNoQuadrupoles(bool usePme) {

    string testName      = "testNoQuadrupoles";

    int inputPmeGridDimension = 10;
    double cutoff             = 0.3;

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    setupMultipoleAmmonia(system, amoebaMultipoleForce, usePme ? AmoebaMultipoleForce::PME : AmoebaMultipoleForce::NoCutoff, AmoebaMultipoleForce::Direct, 
                                             cutoff, inputPmeGridDimension);
    
    // First compute forces with quadrupoles.  This ensures they will be included in the kernel.
    
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    vector<Vec3> forces1;
    double energy1;
    getForcesEnergyMultipoleAmmonia(context, forces1, energy1);
    
    // Now set all quadrupoles to zero and recalculate.
    
    for (int i = 0; i < system.getNumParticles(); i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        amoebaMultipoleForce->getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        for (auto& q : quadrupole)
            q = 0;
        amoebaMultipoleForce->setMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
    }
    amoebaMultipoleForce->updateParametersInContext(context);
    vector<Vec3> forces2;
    double energy2;
    getForcesEnergyMultipoleAmmonia(context, forces2, energy2);
    ASSERT(energy1 != energy2);
    
    // Create a new context with no quadrupoles from the beginning.
    
    LangevinIntegrator integrator2(0.0, 0.1, 0.01);
    Context context2(system, integrator2, Platform::getPlatformByName("CUDA"));
    vector<Vec3> forces3;
    double energy3;
    getForcesEnergyMultipoleAmmonia(context2, forces3, energy3);
    double tolerance = 1e-5;
    compareForcesEnergy(testName, energy2, energy3, forces2, forces3, tolerance);
    
    // If we try to set a quadrupole, that should produce and exception.
    
    double charge, thole, damping, polarity;
    int axisType, atomX, atomY, atomZ;
    vector<double> dipole, quadrupole;
    amoebaMultipoleForce->getMultipoleParameters(0, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
    quadrupole[0] = 1.0;
    amoebaMultipoleForce->setMultipoleParameters(0, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
    bool threwException = false;
    try {
        amoebaMultipoleForce->updateParametersInContext(context2);
    }
    catch (OpenMMException& ex) {
        threwException = true;
    }
    ASSERT(threwException);
}

void testTriclinic() {
    // Create a triclinic box containing eight water molecules.

    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(1.8643, 0, 0), Vec3(-0.16248445120445926, 1.8572057756524414, 0), Vec3(0.16248445120445906, -0.14832299817478897, 1.8512735025730875));
    for (int i = 0; i < 24; i++)
        system.addParticle(1.0);
    AmoebaMultipoleForce* force = new AmoebaMultipoleForce();
    system.addForce(force);
    force->setNonbondedMethod(AmoebaMultipoleForce::PME);
    force->setCutoffDistance(0.7);
    force->setMutualInducedTargetEpsilon(1e-6);
    vector<int> grid(3);
    grid[0] = grid[1] = grid[2] = 24;
    force->setPmeGridDimensions(grid);
    force->setAEwald(5.4459051633620055);
    double o_charge = -0.42616, h_charge = 0.21308;
    vector<double> o_dipole(3), h_dipole(3), o_quadrupole(9, 0.0), h_quadrupole(9, 0.0);
    o_dipole[0] = 0;
    o_dipole[1] = 0;
    o_dipole[2] = 0.0033078867454609203;
    h_dipole[0] = -0.0053536858428776405;
    h_dipole[1] = 0;
    h_dipole[2] = -0.014378273997907321;

    o_quadrupole[0] = 0.00016405937591036892;
    o_quadrupole[4] = -0.00021618201787005826;
    o_quadrupole[8] = 5.212264195968935e-05;
    h_quadrupole[0] = 0.00011465301060008312;
    h_quadrupole[4] = 8.354184196619263e-05;
    h_quadrupole[8] = -0.00019819485256627578;
    h_quadrupole[2] = h_quadrupole[6] = -6.523731100577879e-05;

    for (int i = 0; i < 8; i++) {
        int atom1 = 3*i, atom2 = 3*i+1, atom3 = 3*i+2;
        force->addMultipole(o_charge, o_dipole, o_quadrupole, 1, atom2, atom3, -1, 0.39, pow(0.001*0.92, 1.0/6.0), 0.001*0.92);
        force->addMultipole(h_charge, h_dipole, h_quadrupole, 0, atom1, atom3, -1, 0.39, pow(0.001*0.539, 1.0/6.0), 0.001*0.539);
        force->addMultipole(h_charge, h_dipole, h_quadrupole, 0, atom1, atom2, -1, 0.39, pow(0.001*0.539, 1.0/6.0), 0.001*0.539);
        vector<int> coval1_12(2);
        coval1_12[0] = atom2;
        coval1_12[1] = atom3;
        force->setCovalentMap(atom1, AmoebaMultipoleForce::Covalent12, coval1_12);
        vector<int> coval2_12(1);
        coval2_12[0] = atom1;
        force->setCovalentMap(atom2, AmoebaMultipoleForce::Covalent12, coval2_12);
        force->setCovalentMap(atom3, AmoebaMultipoleForce::Covalent12, coval2_12);
        vector<int> coval2_13(1);
        coval2_13[0] = atom3;
        force->setCovalentMap(atom2, AmoebaMultipoleForce::Covalent13, coval2_13);
        vector<int> coval3_13(1);
        coval3_13[0] = atom2;
        force->setCovalentMap(atom3, AmoebaMultipoleForce::Covalent13, coval3_13);
        vector<int> polar(3);
        polar[0] = atom1;
        polar[1] = atom2;
        polar[2] = atom3;
        force->setCovalentMap(atom1, AmoebaMultipoleForce::PolarizationCovalent11, polar);
        force->setCovalentMap(atom2, AmoebaMultipoleForce::PolarizationCovalent11, polar);
        force->setCovalentMap(atom3, AmoebaMultipoleForce::PolarizationCovalent11, polar);
    }
    vector<Vec3> positions(24);
    positions[0] = Vec3(0.867966, 0.708769, -0.0696862);
    positions[1] = Vec3(0.780946, 0.675579, -0.0382259);
    positions[2] = Vec3(0.872223, 0.681424, -0.161756);
    positions[3] = Vec3(-0.0117313, 0.824445, 0.683762);
    positions[4] = Vec3(0.0216892, 0.789544, 0.605003);
    positions[5] = Vec3(0.0444268, 0.782601, 0.75302);
    positions[6] = Vec3(0.837906, -0.0092611, 0.681463);
    positions[7] = Vec3(0.934042, 0.0098069, 0.673406);
    positions[8] = Vec3(0.793962, 0.0573676, 0.626984);
    positions[9] = Vec3(0.658995, 0.184432, -0.692317);
    positions[10] = Vec3(0.588543, 0.240231, -0.671793);
    positions[11] = Vec3(0.618153, 0.106275, -0.727368);
    positions[12] = Vec3(0.71466, 0.575358, 0.233152);
    positions[13] = Vec3(0.636812, 0.612604, 0.286268);
    positions[14] = Vec3(0.702502, 0.629465, 0.15182);
    positions[15] = Vec3(-0.242658, -0.850419, -0.250483);
    positions[16] = Vec3(-0.169206, -0.836825, -0.305829);
    positions[17] = Vec3(-0.279321, -0.760247, -0.24031);
    positions[18] = Vec3(-0.803838, -0.360559, 0.230369);
    positions[19] = Vec3(-0.811375, -0.424813, 0.301849);
    positions[20] = Vec3(-0.761939, -0.2863, 0.270962);
    positions[21] = Vec3(-0.148063, 0.824409, -0.827221);
    positions[22] = Vec3(-0.20902, 0.868798, -0.7677);
    positions[23] = Vec3(-0.0700878, 0.882333, -0.832221);

    // Compute the forces and energy.

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);

    // Compare them to values computed by Gromacs.

    double expectedEnergy = 6.8278696;
    vector<Vec3> expectedForce(24);
    expectedForce[0] = Vec3(-104.755, 14.0833, 34.359);
    expectedForce[1] = Vec3(97.324, -1.41419, -142.976);
    expectedForce[2] = Vec3(29.3968, -6.31784, -3.8702);
    expectedForce[3] = Vec3(39.1915, 66.852, -28.5767);
    expectedForce[4] = Vec3(-8.17554, -6.71532, 7.63162);
    expectedForce[5] = Vec3(-87.3745, -91.8639, 35.4761);
    expectedForce[6] = Vec3(-19.7568, -40.8693, 5.60238);
    expectedForce[7] = Vec3(0.840984, 26.878, 10.7822);
    expectedForce[8] = Vec3(14.3469, 12.0583, -12.075);
    expectedForce[9] = Vec3(13.757, 16.9954, 46.9403);
    expectedForce[10] = Vec3(-5.04172, -14.008, -15.3804);
    expectedForce[11] = Vec3(-2.1715, 3.7405, -37.2209);
    expectedForce[12] = Vec3(-70.2284, -9.10438, -40.1287);
    expectedForce[13] = Vec3(4.46014, 3.89949, 4.64842);
    expectedForce[14] = Vec3(43.045, -4.79905, 151.879);
    expectedForce[15] = Vec3(20.2129, -0.895376, -27.2086);
    expectedForce[16] = Vec3(-5.10448, 3.57732, 17.0498);
    expectedForce[17] = Vec3(-13.7695, -1.03345, 12.3093);
    expectedForce[18] = Vec3(2.94972, 0.338904, -10.9914);
    expectedForce[19] = Vec3(0.69036, 1.22591, 4.50198);
    expectedForce[20] = Vec3(-4.61495, -2.76981, 3.57732);
    expectedForce[21] = Vec3(73.1489, 16.167, -99.5834);
    expectedForce[22] = Vec3(-31.8235, 6.11282, -21.125);
    expectedForce[23] = Vec3(13.167, 7.42242, 103.102);
    for (int i = 0; i < 24; i++) {
        ASSERT_EQUAL_VEC(expectedForce[i], state.getForces()[i], 1e-2);
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-2);
    
    // Make sure that the induced dipoles match values that were computed with the reference platform.
    
    vector<Vec3> dipole, expectedDipole(24);
    expectedDipole[0] = Vec3(0.001139674, 6.100086e-05, -9.047545e-05);
    expectedDipole[1] = Vec3(0.001196926, 4.956102e-05, -0.001802801);
    expectedDipole[2] = Vec3(0.0004220728, -0.0001239729, -0.0005876792);
    expectedDipole[3] = Vec3(-8.251352e-05, -0.0004572975, 0.000231773);
    expectedDipole[4] = Vec3(-9.605857e-05, -1.712989e-05, 0.0001950891);
    expectedDipole[5] = Vec3(-0.000673353, -0.000748613, 0.0004138051);
    expectedDipole[6] = Vec3(0.0001657927, 0.0003224215, -8.527183e-05);
    expectedDipole[7] = Vec3(4.321082e-05, 0.0001842417, 0.0001044274);
    expectedDipole[8] = Vec3(4.656055e-05, 0.0001173723, -0.000111185);
    expectedDipole[9] = Vec3(-0.0001234674, -0.0002114879, -0.0003179444);
    expectedDipole[10] = Vec3(3.941253e-05, -0.0001807048, -0.0001442925);
    expectedDipole[11] = Vec3(-2.111804e-05, -0.0001092564, -0.0003528088);
    expectedDipole[12] = Vec3(0.0006598309, -9.950432e-05, 0.0002676833);
    expectedDipole[13] = Vec3(-0.0001790531, 6.770892e-05, 0.0003118193);
    expectedDipole[14] = Vec3(4.071723e-05, -7.97693e-05, 0.001696929);
    expectedDipole[15] = Vec3(-0.0002070313, 1.716257e-05, 0.0002235802);
    expectedDipole[16] = Vec3(-0.0001630692, 6.042542e-05, 0.0002146744);
    expectedDipole[17] = Vec3(-0.0001572073, 9.014017e-05, 0.0001155261);
    expectedDipole[18] = Vec3(-1.415344e-05, -2.113717e-05, 8.32503e-05);
    expectedDipole[19] = Vec3(-5.514418e-06, -1.676616e-05, 5.401816e-05);
    expectedDipole[20] = Vec3(-3.962031e-05, -2.67952e-05, 2.322013e-05);
    expectedDipole[21] = Vec3(-0.0007098042, -0.0003138399, 0.000827001);
    expectedDipole[22] = Vec3(-0.0004787119, 0.0001518272, 2.344478e-05);
    expectedDipole[23] = Vec3(-0.0001339792, -1.357387e-05, 0.0008341705);
    force->getInducedDipoles(context, dipole);
    for (int i = 0; i < 24; i++) {
        ASSERT_EQUAL_VEC(expectedDipole[i], dipole[i], 1e-4);
    }

    // Make sure that the electrostatic potential matches values that were computed with the reference platform.

    vector<Vec3> points;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                points.push_back(Vec3(0.5*i, 0.5*j, 0.5*k));
    vector<double> potential;
    force->getElectrostaticPotential(points, context, potential);
    double expectedPotential[] = {15.79581, 11.61804, 7.143405, 15.55374, 20.84724, -117.8448, 10.65319, -23.70798,
            30.26213, -0.008813862, 7.123444, 22.09409, 41.93532, 41.92756, 24.00818, 20.77822, 22.8576, 29.88276,
            -16.35481, 44.62139, -64.10894, -64.28645, -26.05954, -1.25711, -54.76624, -6.754965, 10.83449};
    for (int i = 0; i < 27; i++)
        ASSERT_EQUAL_TOL(expectedPotential[i], potential[i], 1e-4);
}

void testZBisect() {
    System system;
    for (int i = 0; i < 7; i++)
        system.addParticle(1.0);
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    AmoebaMultipoleForce* force = new AmoebaMultipoleForce();
    system.addForce(force);
    force->setNonbondedMethod(AmoebaMultipoleForce::PME);
    force->setCutoffDistance(1.2);
    double charge[] = {-1.01875, 0, 0, 0, -0.51966, 0.25983, 0.25983};
    double dipole[7][3] = {
        {0.06620218576365969, 0.056934176095985306, 0.06298584667720743},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0.007556121391156931},
        {-0.05495981592297553, 0, -0.0030787530116780605},
        {-0.05495981592297553, 0, -0.0030787530116780605}};
    double quadrupole[7][9] = {
        {-0.0004042865090302655, 0.0010450291005955025, -0.0010871640586155112, 0.0010450291005955025, 0.0002512789255424535, -0.0009504541350087216, -0.0010871640586155112, -0.0009504541350087216, 0.00015300758348781198},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0.00035403072392177476, 0.009334284009749387, 0, 0.009334284009749387, -0.00039025708016361216, 0, 0, 0, 3.622635624183737e-05},
        {-3.42848251678095e-05, 0.07374084367702016, -1.894859653979126e-06, 0.07374084367702016, -0.00010024087598069868, 0, -1.894859653979126e-06, 0, 0.00013452570114850818},
        {-3.42848251678095e-05, 0.07374084367702016, -1.894859653979126e-06, 0.07374084367702016, -0.00010024087598069868, 0, -1.894859653979126e-06, 0, 0.00013452570114850818}
    };
    int axis[7][4] = {
        {2, 2, 1, 3},
        {5, -1, -1, -1},
        {5, -1, -1, -1},
        {5, -1, -1, -1},
        {1, 5, 6, -1},
        {0, 4, 6, -1},
        {0, 4, 5, -1}
    };
    double thole = 0.39;
    double damping[] = {0.33178695365189015, 0.33178695365189015, 0.33178695365189015, 0.33178695365189015, 0.306987653777382, 0.2813500172269554, 0.2813500172269554};
    double polarity[] = {0.001334, 0.001334, 0.001334, 0.001334, 0.000837, 0.000496, 0.000496};
    for (int i = 0; i < 7; i++) {
        vector<double> d, q;
        for (int j = 0; j < 3; j++)
            d.push_back(dipole[i][j]);
        for (int j = 0; j < 9; j++)
            q.push_back(quadrupole[i][j]);
        force->addMultipole(charge[i], d, q, axis[i][0], axis[i][1], axis[i][2], axis[i][3], thole, damping[i], polarity[i]);
    }
    for (int i = 0; i < 4; i++) {
        vector<int> map;
        if (i != 0) map.push_back(0);
        force->setCovalentMap(i, AmoebaMultipoleForce::Covalent12, map);
        map.clear();
        if (i != 1) map.push_back(1);
        if (i != 2) map.push_back(2);
        if (i != 3) map.push_back(3);
        force->setCovalentMap(i, AmoebaMultipoleForce::Covalent13, map);
        map.clear();
        map.push_back(0);
        map.push_back(1);
        map.push_back(2);
        map.push_back(3);
        force->setCovalentMap(i, AmoebaMultipoleForce::PolarizationCovalent11, map);
    }
    for (int i = 4; i < 7; i++) {
        vector<int> map;
        if (i != 4) map.push_back(4);
        force->setCovalentMap(i, AmoebaMultipoleForce::Covalent12, map);
        map.clear();
        if (i != 5) map.push_back(5);
        if (i != 6) map.push_back(6);
        force->setCovalentMap(i, AmoebaMultipoleForce::Covalent13, map);
        map.clear();
        map.push_back(4);
        map.push_back(5);
        map.push_back(6);
        force->setCovalentMap(i, AmoebaMultipoleForce::PolarizationCovalent11, map);
    }
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    vector<Vec3> positions;
    positions.push_back(Vec3(-0.06317711175870899, -0.04905009196658128, 0.0767217));
    positions.push_back(Vec3(-0.049166918626451395, -0.20747614470348363, 0.03979849999999996));
    positions.push_back(Vec3(-0.19317150000000005, -0.05811762921948427, 0.1632788999999999));
    positions.push_back(Vec3(0.04465103038516016, -0.018345116763806235, 0.18531239999999993));
    positions.push_back(Vec3(0.005630299999999998, 0.40965770000000035, 0.5731495));
    positions.push_back(Vec3(0.036148100000000016, 0.3627041999999996, 0.49299430000000033));
    positions.push_back(Vec3(0.07781149999999992, 0.4178183000000004, 0.6355703000000004));
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(-84.1532, state.getPotentialEnergy(), 0.01);
}

void testZOnly() {
    int numParticles = 3;
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(1.0);
    AmoebaMultipoleForce* force = new AmoebaMultipoleForce();
    system.addForce(force);
    vector<double> d(3), q(9, 0.0);
    d[0] = 0.05;
    d[1] = -0.05;
    d[2] = 0.1;
    force->addMultipole(0.0, d, q, AmoebaMultipoleForce::ZOnly, 1, 0, 0, 0.39, 0.33, 0.001);
    force->addMultipole(0.0, d, q, AmoebaMultipoleForce::Bisector, 0, 2, 0, 0.39, 0.33, 0.001);
    force->addMultipole(0.0, d, q, AmoebaMultipoleForce::ZOnly, 1, 0, 0, 0.39, 0.33, 0.001);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0.2, 0, 0);
    positions[2] = Vec3(0.2, 0.2, -0.05);
    
    // Evaluate the forces.
    
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(positions);
    State state = context.getState(State::Forces);
    double norm = 0.0;
    for (Vec3 f : state.getForces())
        norm += f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
    norm = std::sqrt(norm);

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    const double delta = 1e-3;
    double step = 0.5*delta/norm;
    vector<Vec3> positions2(numParticles), positions3(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        Vec3 p = positions[i];
        Vec3 f = state.getForces()[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(state2.getPotentialEnergy(), state3.getPotentialEnergy()+norm*delta, 1e-3)
}


int main(int argc, char* argv[]) {
    try {
        std::cout << "TestCudaAmoebaMultipoleForce running test..." << std::endl;
        registerAmoebaCudaKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));

        // tests using two ammonia molecules

        // test direct polarization, no cutoff

        testMultipoleAmmoniaDirectPolarization();

        // test mutual polarization, no cutoff

        testMultipoleAmmoniaMutualPolarization();

        // test multipole direct & mutual polarization using PME

        testMultipoleWaterPMEDirectPolarization();
        testMultipoleWaterPMEMutualPolarization();

        // check validation of traceless/symmetric quadrupole tensor

        testQuadrupoleValidation();

        // system w/ 2 ions and 2 water molecules

        testMultipoleIonsAndWaterPMEMutualPolarization();
        testMultipoleIonsAndWaterPMEDirectPolarization();

        // test querying induced dipoles
        
        testParticleInducedDipoles();
        
        // test computation of system multipole moments

        testSystemMultipoleMoments();

        // test computation of grid potential

        testMultipoleGridPotential();

        // large box of water
        testPMEMutualPolarizationLargeWater();
        
        testNoQuadrupoles(false);
        testNoQuadrupoles(true);
        
        // triclinic box of water
        
        testTriclinic();
        
        // test the ZBisect axis type.
        
        testZBisect();
        
        // test the ZOnly axis type.
        
        testZOnly();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
