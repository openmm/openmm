/* -------------------------------------------------------------------------- *
 *                                   OpenMMAmoeba                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
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
 * This tests the CUDA implementation of the extrapolated polarization algorithms in AmoebaMultipoleForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMAmoeba.h"
#include "openmm/System.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Vec3.h"
#include "openmm/Units.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol, testname) {double _norm_ = std::sqrt(expected.dot(expected)); double _scale_ = _norm_ > 1.0 ? _norm_ : 1.0; if ((std::abs((expected[0])-(found[0]))/_scale_ > (tol)) || (std::abs((expected[1])-(found[1]))/_scale_ > (tol)) || (std::abs((expected[2])-(found[2]))/_scale_ > (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};


using namespace OpenMM;
using namespace std;

extern "C" void registerAmoebaCudaKernelFactories();

const double TOL = 1e-4;


// print the energy and forces out, in AKMA units, to allow comparison with TINKER

static void printEnergyAndForces(double energy, vector<Vec3> &forces){
    size_t natoms = forces.size();
    double sf = 1.0;
    std::cout << "Energy (SI):" << std::setw(20) << std::setprecision(10) << energy << std::endl;
    std::cout << "Forces (SI):" << std::endl;
    for(int i = 0; i < natoms; ++i){
        std::cout << i+1 << "\t" << std::setw(20) << std::setprecision(10) << forces[i][0]*sf <<
                                    std::setw(20) << std::setprecision(10) << forces[i][1]*sf <<
                                    std::setw(20) << std::setprecision(10) << forces[i][2]*sf << std::endl;
    }

    sf = -OpenMM::KcalPerKJ/10.0;
    std::cout << "Energy (AKMA):" << std::setw(20) << std::setprecision(10) << energy*OpenMM::KcalPerKJ << std::endl;
    std::cout << "Forces (AKMA):" << std::endl;
    for(int i = 0; i < natoms; ++i){
        std::cout << i+1 << "\t" << std::setw(20) << std::setprecision(10) << forces[i][0]*sf <<
                                    std::setw(20) << std::setprecision(10) << forces[i][1]*sf <<
                                    std::setw(20) << std::setprecision(10) << forces[i][2]*sf << std::endl;
    }
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


vector<Vec3> setupWaterDimer(System& system,  AmoebaMultipoleForce* amoebaMultipoleForce, bool use_pol_groups) {
    const int NATOMS = 6;
    const char* atom_types[NATOMS] = {"O", "H1", "H2", "O", "H1", "H2"};
    const double coords[NATOMS][3] = {
        {  2.000000, 2.000000, 2.000000},
        {  2.500000, 2.000000, 3.000000},
        {  1.500000, 2.000000, 3.000000},
        {  0.000000, 0.000000, 0.000000},
        {  0.500000, 0.000000, 1.000000},
        { -0.500000, 0.000000, 1.000000}
    };

    std::map < std::string, double > tholemap;
    std::map < std::string, double > polarmap;
    std::map < std::string, double > chargemap;
    std::map < std::string, std::vector<double> > dipolemap;
    std::map < std::string, std::vector<double> > quadrupolemap;
    std::map < std::string, AmoebaMultipoleForce::MultipoleAxisTypes > axesmap;
    std::map < std::string, std::vector<int> > anchormap;
    std::map < std::string, double > massmap;
    std::map < std::string, std::vector<int> > polgrpmap;
    std::map < std::string, std::vector<int> > cov12map;
    std::map < std::string, std::vector<int> > cov13map;

    axesmap["O"]  = AmoebaMultipoleForce::Bisector;
    axesmap["H1"] = AmoebaMultipoleForce::ZThenX;
    axesmap["H2"] = AmoebaMultipoleForce::ZThenX;

    chargemap["O"]  = -0.51966;
    chargemap["H1"] = 0.25983;
    chargemap["H2"] = 0.25983;

    int oanc[3] = {1, 2, 0};
    int h1anc[3] = {-1, 1, 0};
    int h2anc[3] = {-2, -1, 0};
    std::vector<int> oancv(&oanc[0], &oanc[3]);
    std::vector<int> h1ancv(&h1anc[0], &h1anc[3]);
    std::vector<int> h2ancv(&h2anc[0], &h2anc[3]);
    anchormap["O"]  = oancv;
    anchormap["H1"] = h1ancv;
    anchormap["H2"] = h2ancv;

    double od[3] = {0.0, 0.0, 0.00755612136146};
    double hd[3] = {-0.00204209484795, 0.0, -0.00307875299958};
    std::vector<double> odv(&od[0], &od[3]);
    std::vector<double> hdv(&hd[0], &hd[3]);
    dipolemap["O"]  = odv;
    dipolemap["H1"] = hdv;
    dipolemap["H2"] = hdv;

    double oq[9] = {0.000354030721139,  0.0, 0.0,
                    0.0, -0.000390257077096, 0.0,
                    0.0, 0.0,  3.62263559571e-05};
    double hq[9] = {-3.42848248983e-05, 0.0, -1.89485963908e-06,
                     0.0,          -0.000100240875193,      0.0,
                    -1.89485963908e-06, 0.0,  0.000134525700091};
    std::vector<double> oqv(&oq[0], &oq[9]);
    std::vector<double> hqv(&hq[0], &hq[9]);
    quadrupolemap["O"]  = oqv;
    quadrupolemap["H1"] = hqv;
    quadrupolemap["H2"] = hqv;

    polarmap["O"]  = 0.3069876538;
    polarmap["H1"] = 0.2813500172;
    polarmap["H2"] = 0.2813500172;

    polarmap["O"]  = 0.000837;
    polarmap["H1"] = 0.000496;
    polarmap["H2"] = 0.000496;

    tholemap["O"]  = 0.3900;
    tholemap["H1"] = 0.3900;
    tholemap["H2"] = 0.3900;

    massmap["O"]  = 15.999;
    massmap["H1"] = 1.0080000;
    massmap["H2"] = 1.0080000;

    int opg[3] = {0,1,2};
    int h1pg[3] = {-1,0,1};
    int h2pg[3] = {-2,-1,0};
    std::vector<int> opgv(&opg[0], &opg[3]);
    std::vector<int> h1pgv(&h1pg[0], &h1pg[3]);
    std::vector<int> h2pgv(&h2pg[0], &h2pg[3]);
    if(!use_pol_groups){
        opgv.clear();
        h1pgv.clear();
        h2pgv.clear();
    }
    polgrpmap["O"] = opgv;
    polgrpmap["H1"] = h1pgv;
    polgrpmap["H2"] = h2pgv;

    int cov12o[2] = {1,2};
    int cov12h1[1] = {-1};
    int cov12h2[1] = {-2};
    std::vector<int> cov12ov(&cov12o[0], &cov12o[2]);
    std::vector<int> cov12h1v(&cov12h1[0], &cov12h1[1]);
    std::vector<int> cov12h2v(&cov12h2[0], &cov12h2[1]);
    cov12map["O"] = cov12ov;
    cov12map["H1"] = cov12h1v;
    cov12map["H2"] = cov12h2v;

    int cov13h1[1] = {1};
    int cov13h2[1] = {-1};
    std::vector<int> cov13h1v(&cov13h1[0], &cov13h1[1]);
    std::vector<int> cov13h2v(&cov13h2[0], &cov13h2[1]);
    cov13map["O"] = std::vector<int>();
    cov13map["H1"] = cov13h1v;
    cov13map["H2"] = cov13h2v;


    std::vector<Vec3> positions(NATOMS);
    for(int atom = 0; atom < NATOMS; ++atom){
        const char* element = atom_types[atom];
        double damp = polarmap[element];
        double alpha = pow(damp, 1.0/6.0);
        int atomz = atom + anchormap[element][0];
        int atomx = atom + anchormap[element][1];
        int atomy = anchormap[element][2]==0 ? -1 : atom + anchormap[element][2];
        amoebaMultipoleForce->addMultipole(chargemap[element], dipolemap[element], quadrupolemap[element],
                                           axesmap[element], atomz, atomx, atomy, tholemap[element], alpha, damp);
        system.addParticle(massmap[element]);
        double offset =0.0;
        positions[atom] = Vec3(coords[atom][0]+offset, coords[atom][1]+offset, coords[atom][2]+offset)*OpenMM::NmPerAngstrom;
        // Polarization groups
        std::vector<int> tmppol;
        std::vector<int>& polgrps = polgrpmap[element];
        for(int i=0; i < polgrps.size(); ++i)
            tmppol.push_back(polgrps[i]+atom);
        if(!tmppol.empty())
           amoebaMultipoleForce->setCovalentMap(atom, AmoebaMultipoleForce::PolarizationCovalent11, tmppol);
        // 1-2 covalent groups
        std::vector<int> tmp12;
        std::vector<int>& cov12s = cov12map[element];
        for(int i=0; i < cov12s.size(); ++i)
            tmp12.push_back(cov12s[i]+atom);
        if(!tmp12.empty())
           amoebaMultipoleForce->setCovalentMap(atom, AmoebaMultipoleForce::Covalent12, tmp12);
        // 1-3 covalent groups
        std::vector<int> tmp13;
        std::vector<int>& cov13s = cov13map[element];
        for(int i=0; i < cov13s.size(); ++i)
            tmp13.push_back(cov13s[i]+atom);
        if(!tmp13.empty())
           amoebaMultipoleForce->setCovalentMap(atom, AmoebaMultipoleForce::Covalent13, tmp13);
    }

    system.addForce(amoebaMultipoleForce);

    return positions;
}


static void check_finite_differences(vector<Vec3> analytic_forces, Context &context, vector<Vec3> positions)
{
    double tol = 1e-5;
    // We allow more permissive testing for single precision.
    if(Platform::getPlatformByName("CUDA").getPropertyValue(context, "Precision") != "double") tol = 5e-4;

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.
    double norm = 0.0;
    for (auto& f : analytic_forces)
        norm += f.dot(f);
    norm = std::sqrt(norm);
    const double stepSize = 1e-3;
    double step = 0.5*stepSize/norm;
    vector<Vec3> positions2(analytic_forces.size()), positions3(analytic_forces.size());
    for (int i = 0; i < (int) positions.size(); ++i) {
        Vec3 p = positions[i];
        Vec3 f = analytic_forces[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, tol);
}


static void testWaterDimerTriclinicPME() {

    std::string testName      = "testWaterDimerTriclinicPME";

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    vector<Vec3> coords = setupWaterDimer(system, amoebaMultipoleForce, true);

    system.setDefaultPeriodicBoxVectors(Vec3(2.0, 0.0, 0.0),
                                        Vec3(0.2, 2.0, 0.0),
                                        Vec3(0.1, 0.5, 2.0));
    amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::PME);
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Extrapolated);
    std::vector<double> coefs;
    coefs.push_back(0.0);  // The mu_0 coefficient
    coefs.push_back(-0.3); // The mu_1 coefficient
    coefs.push_back(0.0);  // The mu_2 coefficient
    coefs.push_back(1.3);  // The mu_3 coefficient
    amoebaMultipoleForce->setExtrapolationCoefficients(coefs);
    amoebaMultipoleForce->setCutoffDistance(9.0*OpenMM::NmPerAngstrom);
    amoebaMultipoleForce->setAEwald(4);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-06);
    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = 64;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(coords);
    OpenMM::State state = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces = state.getForces();
    double energy = state.getPotentialEnergy();
//    printEnergyAndForces(energy, forces);

    double expectedEnergy     = -1.945797427;
    std::vector<Vec3> expectedForces(forces.size());
    expectedForces[0] = Vec3(  -131.1099603,   -187.2725558,    36.94657685);
    expectedForces[1] = Vec3(    38.6397841,    2.410997985,    8.008437937);
    expectedForces[2] = Vec3(   38.69034185,    117.5018257,    32.43097836);
    expectedForces[3] = Vec3(  -117.3212339,   -102.3366145,   -30.50621066);
    expectedForces[4] = Vec3(   124.8343077,    169.7729804,   -24.10742414);
    expectedForces[5] = Vec3(   46.26244074, -0.07194110719,   -22.77727325);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    check_finite_differences(forces, context, coords);
}


static void testWaterDimerTriclinicPMENoPolGroups() {

    std::string testName      = "testWaterDimerTriclinicPMENoPolGroups";

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    vector<Vec3> coords = setupWaterDimer(system, amoebaMultipoleForce, false);

    system.setDefaultPeriodicBoxVectors(Vec3(2.0, 0.0, 0.0),
                                        Vec3(0.2, 2.0, 0.0),
                                        Vec3(0.1, 0.5, 2.0));
    amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::PME);
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Extrapolated);
    std::vector<double> coefs;
    coefs.push_back(0.0);  // The mu_0 coefficient
    coefs.push_back(-0.3); // The mu_1 coefficient
    coefs.push_back(0.0);  // The mu_2 coefficient
    coefs.push_back(1.3);  // The mu_3 coefficient
    amoebaMultipoleForce->setExtrapolationCoefficients(coefs);
    amoebaMultipoleForce->setCutoffDistance(9.0*OpenMM::NmPerAngstrom);
    amoebaMultipoleForce->setAEwald(4);
    amoebaMultipoleForce->setEwaldErrorTolerance(1.0e-06);
    std::vector<int> pmeGridDimension(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = 64;
    amoebaMultipoleForce->setPmeGridDimensions(pmeGridDimension);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(coords);
    OpenMM::State state = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces = state.getForces();
    double energy = state.getPotentialEnergy();
//    printEnergyAndForces(energy, forces);

    double expectedEnergy     =  -1.840068409;
    std::vector<Vec3> expectedForces(forces.size());

    expectedForces[0] = Vec3(  -69.85154559,  -104.2092334,   3.586495334);
    expectedForces[1] = Vec3(   19.50350452,   -14.5844519,   9.400418341);
    expectedForces[2] = Vec3(   16.75641493,   75.15006506,   19.14553199);
    expectedForces[3] = Vec3(  -67.24268213,  -47.39994175,  -18.81277222);
    expectedForces[4] = Vec3(   75.15808251,   110.6109313,   4.355432435);
    expectedForces[5] = Vec3(   25.67255306,  -19.56378113,  -17.68217953);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    check_finite_differences(forces, context, coords);
}


static void testWaterDimerNoCutoff() {

    std::string testName      = "testWaterDimerNoCutoff";

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    vector<Vec3> coords = setupWaterDimer(system, amoebaMultipoleForce, true);

    amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::NoCutoff);
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Extrapolated);
    std::vector<double> coefs;
    coefs.push_back(0.0);  // The mu_0 coefficient
    coefs.push_back(-0.3); // The mu_1 coefficient
    coefs.push_back(0.0);  // The mu_2 coefficient
    coefs.push_back(1.3);  // The mu_3 coefficient
    amoebaMultipoleForce->setExtrapolationCoefficients(coefs);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(coords);
    OpenMM::State state = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces = state.getForces();
    double energy = state.getPotentialEnergy();
//    printEnergyAndForces(energy, forces);

    double expectedEnergy     = -1.399194432;
    std::vector<Vec3> expectedForces(forces.size());

    expectedForces[0] = Vec3( -130.7294487,   -186.3287444,    41.40628056);
    expectedForces[1] = Vec3(  38.90143386,    2.140957908,    5.564712102);
    expectedForces[2] = Vec3(  38.32881448,    117.0462626,    29.90093041);
    expectedForces[3] = Vec3( -117.1147396,   -101.6981494,   -25.55733439);
    expectedForces[4] = Vec3(  124.7421318,    169.1571359,   -26.38724373);
    expectedForces[5] = Vec3(  45.87180816,  -0.3174626947,   -24.92734495);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    check_finite_differences(forces, context, coords);
}


static void testWaterDimerNoCutoffNoPolGroups() {

    std::string testName      = "testWaterDimerNoCutoffNoPolGroups";

    System system;
    AmoebaMultipoleForce* amoebaMultipoleForce = new AmoebaMultipoleForce();;
    vector<Vec3> coords = setupWaterDimer(system, amoebaMultipoleForce, false);

    amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::NoCutoff);
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Extrapolated);
    std::vector<double> coefs;
    coefs.push_back(0.0);  // The mu_0 coefficient
    coefs.push_back(-0.3); // The mu_1 coefficient
    coefs.push_back(0.0);  // The mu_2 coefficient
    coefs.push_back(1.3);  // The mu_3 coefficient
    amoebaMultipoleForce->setExtrapolationCoefficients(coefs);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName("CUDA"));
    context.setPositions(coords);
    OpenMM::State state = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces = state.getForces();
    double energy = state.getPotentialEnergy();
//    printEnergyAndForces(energy, forces);

    double expectedEnergy = -1.56926564;
    std::vector<Vec3> expectedForces(forces.size());

    expectedForces[0] = Vec3(   -69.623843,  -103.7006124,   6.162774255);
    expectedForces[1] = Vec3(  19.54326912,  -14.69441322,   8.014369439);
    expectedForces[2] = Vec3(  16.65441143,   74.88100242,   17.70364405);
    expectedForces[3] = Vec3( -67.10049929,  -47.08900953,  -16.01092086);
    expectedForces[4] = Vec3(  74.98800293,   110.2649458,   3.020145768);
    expectedForces[5] = Vec3(  25.53865881,  -19.66191302,  -18.89001266);

    double tolerance          = 1.0e-04;
    compareForcesEnergy(testName, expectedEnergy, energy, expectedForces, forces, tolerance);
    check_finite_differences(forces, context, coords);
}



int main(int numberOfArguments, char* argv[]) {

    try {
        registerAmoebaCudaKernelFactories();
        if (numberOfArguments > 1)
            Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", std::string(argv[1]));

        /*
         * Water dimer energy / force tests under various conditions.
         */

        // PME, triclinic
        testWaterDimerTriclinicPME();

        // PME, triclinic, no polarization groups
        testWaterDimerTriclinicPMENoPolGroups();

        // No cutoff
        testWaterDimerNoCutoff();

        // No cutoff, no polarization groups
        testWaterDimerNoCutoffNoPolGroups();

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
