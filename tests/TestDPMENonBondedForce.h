/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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
#include "openmm/HarmonicBondForce.h"
#include "openmm/DPMENonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/Units.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <iomanip>

#include <vector>

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

void testExclusionsAnd14() {
    System system;
    DPMENonbondedForce* nonbonded = new DPMENonbondedForce();
    for (int i = 0; i < 5; ++i) {
        system.addParticle(1.0);
        nonbonded->addParticle(0, 1.5, 0);
    }
    vector<pair<int, int> > bonds;
    bonds.push_back(pair<int, int>(0, 1));
    bonds.push_back(pair<int, int>(1, 2));
    bonds.push_back(pair<int, int>(2, 3));
    bonds.push_back(pair<int, int>(3, 4));
    nonbonded->createExceptionsFromBonds(bonds, 0.0, 0.0);
    int first14, second14;
    for (int i = 0; i < nonbonded->getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        nonbonded->getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if ((particle1 == 0 && particle2 == 3) || (particle1 == 3 && particle2 == 0))
            first14 = i;
        if ((particle1 == 1 && particle2 == 4) || (particle1 == 4 && particle2 == 1))
            second14 = i;
    }
    system.addForce(nonbonded);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    for (int i = 1; i < 5; ++i) {

        // Test LJ forces

        vector<Vec3> positions(5);
        const double r = 1.0;
        for (int j = 0; j < 5; ++j) {
            nonbonded->setParticleParameters(j, 0, 1.5, 0);
            positions[j] = Vec3(0, j, 0);
        }
        nonbonded->setParticleParameters(0, 0, 1.5, 1);
        nonbonded->setParticleParameters(i, 0, 1.5, 1);
        nonbonded->setExceptionParameters(first14, 0, 3, 0, 1.5, i == 3 ? 0.5 : 0.0);
        nonbonded->setExceptionParameters(second14, 1, 4, 0, 1.5, 0.0);
        positions[i] = Vec3(r, 0, 0);
        context.reinitialize();
        context.setPositions(positions);
        State state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces = state.getForces();
        double x = 1.5/r;
        double eps = 1.0;
        double force = 4.0*eps*(12*std::pow(x, 12.0)-6*std::pow(x, 6.0))/r;
        double energy = 4.0*eps*(std::pow(x, 12.0)-std::pow(x, 6.0));
        if (i == 3) {
            force *= 0.5;
            energy *= 0.5;
        }
        if (i < 3) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);

        // Test Coulomb forces

        nonbonded->setParticleParameters(0, 2, 1.5, 0);
        nonbonded->setParticleParameters(i, 2, 1.5, 0);
        nonbonded->setExceptionParameters(first14, 0, 3, i == 3 ? 4/1.2 : 0, 1.5, 0);
        nonbonded->setExceptionParameters(second14, 1, 4, 0, 1.5, 0);
        context.reinitialize();
        context.setPositions(positions);
        state = context.getState(State::Forces | State::Energy);
        const vector<Vec3>& forces2 = state.getForces();
        force = ONE_4PI_EPS0*4/(r*r);
        energy = ONE_4PI_EPS0*4/r;
        if (i == 3) {
            force /= 1.2;
            energy /= 1.2;
        }
        if (i < 3) {
            force = 0;
            energy = 0;
        }
        ASSERT_EQUAL_VEC(Vec3(-force, 0, 0), forces2[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(force, 0, 0), forces2[i], TOL);
        ASSERT_EQUAL_TOL(energy, state.getPotentialEnergy(), TOL);
    }
}


void make_waterbox(int natoms, double boxEdgeLength, DPMENonbondedForce *forceField,  vector<Vec3> &positions, vector<double>& eps, vector<double>& sig,
                   vector<pair<int, int> >& bonds, System &system, bool do_electrostatics) {
    const int RESSIZE = 3;
    const double masses[RESSIZE]    = {     8.0,    1.0,    1.0 };
    const double charges[RESSIZE]   = {  -0.834,  0.417,  0.417 };
    // Values from the CHARMM force field, in AKMA units
    const double epsilons[RESSIZE]  = { -0.1521, -0.046, -0.046 };
    const double halfrmins[RESSIZE] = {  1.7682, 0.2245, 0.2245 };
    positions.clear();
    if(natoms == 6){
        const double coords[6][3] = {
            {  2.000000, 2.000000, 2.000000},
            {  2.500000, 2.000000, 3.000000},
            {  1.500000, 2.000000, 3.000000},
            {  0.000000, 0.000000, 0.000000},
            {  0.500000, 0.000000, 1.000000},
            { -0.500000, 0.000000, 1.000000}
        };
        for(int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }else if(natoms == 375){
        const double coords[375][3] = {
            { -6.227577, -6.257691, -6.242372 },
            { -5.323856, -6.038433, -6.004096 },
            { -6.756728, -5.564157, -5.841659 },
            { -3.045385, -6.233543, -6.190951 },
            { -3.526738, -5.554781, -5.712262 },
            { -3.599325, -6.431115, -6.949685 },
            {  0.021426, -6.231744, -6.149739 },
            { -0.875209, -5.976613, -6.378901 },
            {  0.535110, -6.037842, -6.937199 },
            {  3.100433, -6.203101, -6.272233 },
            {  3.874100, -6.357529, -5.725284 },
            {  2.378191, -6.111572, -5.646492 },
            {  6.186877, -6.143892, -6.207648 },
            {  6.463778, -6.660529, -5.447403 },
            {  6.262846, -6.748390, -6.949519 },
            { -6.215886, -3.157497, -6.242725 },
            { -6.237827, -3.077604, -5.286327 },
            { -6.028329, -2.268904, -6.553829 },
            { -3.143522, -3.072478, -6.164691 },
            { -3.384628, -3.637034, -6.902733 },
            { -2.183705, -3.058977, -6.176044 },
            { -0.001261, -3.162923, -6.234377 },
            { -0.032904, -2.309579, -6.672970 },
            {  0.052913, -2.950816, -5.299691 },
            {  3.087208, -3.110696, -6.143281 },
            {  2.657287, -2.555869, -6.798183 },
            {  3.804920, -3.533514, -6.620435 },
            {  6.169768, -3.141467, -6.167741 },
            {  7.044757, -3.328081, -6.515780 },
            {  5.953438, -2.272871, -6.514589 },
            { -6.207742, -0.041917, -6.155951 },
            { -5.438834,  0.329527, -6.594553 },
            { -6.956596,  0.335819, -6.622957 },
            { -3.102424, -0.064365, -6.193126 },
            { -3.759220,  0.421493, -6.697232 },
            { -2.461457,  0.600160, -5.930203 },
            {  0.057084, -0.011721, -6.176128 },
            { -0.107739,  0.025184, -7.121132 },
            { -0.798344,  0.160864, -5.776104 },
            {  3.038705,  0.009068, -6.198274 },
            {  3.541849,  0.080397, -7.012718 },
            {  3.690240, -0.224328, -5.533000 },
            {  6.176831,  0.057180, -6.193629 },
            {  5.789021, -0.734630, -6.573369 },
            {  7.097056, -0.172982, -6.046071 },
            { -6.202562,  3.152015, -6.250632 },
            { -6.598511,  3.181133, -5.376596 },
            { -5.879135,  2.252396, -6.338049 },
            { -3.097153,  3.045646, -6.175594 },
            { -3.888514,  3.580848, -6.269758 },
            { -2.415842,  3.541069, -6.635956 },
            {  0.004020,  3.065952, -6.260601 },
            { -0.714153,  3.647384, -6.000359 },
            {  0.650339,  3.152212, -5.556049 },
            {  3.148082,  3.065474, -6.236447 },
            {  3.111314,  3.313051, -5.309670 },
            {  2.384648,  3.494128, -6.630129 },
            {  6.195396,  3.144020, -6.257042 },
            {  6.825505,  3.254686, -5.541307 },
            {  5.765894,  2.305751, -6.071584 },
            { -6.224762,  6.269436, -6.199499 },
            { -6.227138,  5.744470, -7.003219 },
            { -5.898129,  5.671700, -5.523050 },
            { -3.049795,  6.244848, -6.205918 },
            { -3.087748,  5.286175, -6.173255 },
            { -3.968304,  6.520275, -6.251129 },
            { -0.053210,  6.218043, -6.169688 },
            {  0.825799,  6.589399, -6.064810 },
            {  0.018791,  5.642523, -6.934642 },
            {  3.109156,  6.250395, -6.157883 },
            {  3.646142,  5.470269, -6.314764 },
            {  2.467680,  6.248139, -6.872071 },
            {  6.221542,  6.201019, -6.270430 },
            {  5.372468,  6.425629, -5.882893 },
            {  6.803922,  6.076521, -5.517505 },
            { -6.192636, -6.150665, -3.134919 },
            { -6.377934, -7.011364, -3.517572 },
            { -6.257279, -6.290045, -2.187314 },
            { -3.109530, -6.270465, -3.117828 },
            { -2.298394, -5.772194, -2.993950 },
            { -3.809490, -5.627637, -2.982220 },
            { -0.038820, -6.182503, -3.150607 },
            { -0.071984, -7.054630, -2.750788 },
            {  0.688173, -5.741421, -2.705086 },
            {  3.102774, -6.145922, -3.077989 },
            {  2.352681, -6.724399, -3.233832 },
            {  3.862449, -6.652292, -3.374712 },
            {  6.221089, -6.201529, -3.167266 },
            {  6.825244, -6.363707, -2.439079 },
            {  5.358331, -6.130338, -2.752367 },
            { -6.262103, -3.132226, -3.123586 },
            { -6.166340, -2.275929, -2.700331 },
            { -5.366218, -3.471706, -3.184445 },
            { -3.111137, -3.051004, -3.142406 },
            { -3.311768, -3.968529, -3.341042 },
            { -2.770611, -3.068346, -2.245019 },
            {  0.005724, -3.130303, -3.163012 },
            {  0.482202, -2.374690, -2.811473 },
            { -0.573062, -3.403472, -2.447507 },
            {  3.095078, -3.099676, -3.168462 },
            {  2.417453, -3.190640, -2.494584 },
            {  3.919837, -3.073655, -2.677893 },
            {  6.196105, -3.040136, -3.088337 },
            {  5.640465, -3.619831, -3.614436 },
            {  6.939673, -3.589544, -2.829842 },
            { -6.182065, -0.003400, -3.042390 },
            { -6.005810, -0.591008, -3.780776 },
            { -6.797182,  0.644978, -3.392815 },
            { -3.055877, -0.038021, -3.078205 },
            { -2.957574,  0.804401, -3.527898 },
            { -4.001938, -0.200904, -3.077218 },
            { -0.034735,  0.036962, -3.061441 },
            { -0.339387, -0.376429, -3.872525 },
            {  0.890739, -0.210258, -2.998671 },
            {  3.137299, -0.056306, -3.103061 },
            {  3.449365,  0.817518, -3.349247 },
            {  2.217754,  0.076224, -2.861330 },
            {  6.202733, -0.057422, -3.135143 },
            {  6.891866,  0.607339, -3.204067 },
            {  5.583070,  0.304112, -2.497277 },
            { -6.235328,  3.091028, -3.161196 },
            { -5.620224,  3.796590, -2.948140 },
            { -6.337344,  2.604974, -2.339669 },
            { -3.100761,  3.085670, -3.040001 },
            { -3.841514,  3.471502, -3.513266 },
            { -2.405567,  3.015122, -3.698249 },
            {  0.017686,  3.040997, -3.114710 },
            { -0.562190,  3.593571, -3.643829 },
            {  0.281457,  3.602138, -2.381834 },
            {  3.040493,  3.118625, -3.091261 },
            {  3.497428,  2.301545, -2.878768 },
            {  3.706278,  3.661985, -3.519106 },
            {  6.153747,  3.141770, -3.110452 },
            {  6.524448,  2.527782, -3.748544 },
            {  6.728033,  3.068367, -2.344704 },
            { -6.228412,  6.153262, -3.134991 },
            { -5.497842,  6.211200, -2.514933 },
            { -6.569498,  7.048988, -3.188810 },
            { -3.114261,  6.247531, -3.052990 },
            { -3.763440,  5.838093, -3.629608 },
            { -2.269351,  5.925773, -3.375729 },
            {  0.033856,  6.255675, -3.075505 },
            {  0.340270,  5.631131, -3.737032 },
            { -0.877663,  6.003468, -2.910926 },
            {  3.072932,  6.151993, -3.082847 },
            {  3.294846,  6.929869, -2.565910 },
            {  3.393954,  6.350465, -3.965522 },
            {  6.229389,  6.147520, -3.121044 },
            {  5.797138,  6.388219, -2.298375 },
            {  6.254689,  6.963103, -3.626757 },
            { -6.218080, -6.208893, -0.064724 },
            { -5.793234, -6.874918,  0.480691 },
            { -6.438098, -5.502241,  0.546665 },
            { -3.161793, -6.216735, -0.028339 },
            { -2.502753, -6.872853,  0.209857 },
            { -2.775567, -5.379832,  0.239959 },
            { -0.002532, -6.140991, -0.003478 },
            {  0.689039, -6.720671, -0.330985 },
            { -0.648843, -6.734293,  0.386184 },
            {  3.039844, -6.205646, -0.017228 },
            {  3.777609, -6.569800, -0.511855 },
            {  3.436401, -5.858902,  0.785313 },
            {  6.257004, -6.165220, -0.008551 },
            {  5.368434, -6.091008, -0.364215 },
            {  6.245060, -6.979374,  0.499948 },
            { -6.246883, -3.055033, -0.019337 },
            { -6.350454, -3.640127,  0.734652 },
            { -5.423699, -3.332793, -0.427723 },
            { -3.094331, -3.063472,  0.051555 },
            { -2.443051, -3.620051, -0.381590 },
            { -3.906094, -3.218919, -0.436745 },
            {  0.059986, -3.108175,  0.026527 },
            { -0.319331, -2.356822, -0.435152 },
            { -0.632825, -3.772569,  0.014099 },
            {  3.056825, -3.094971, -0.045775 },
            {  3.286789, -3.903484,  0.417890 },
            {  3.657675, -2.435504,  0.308690 },
            {  6.204418, -3.046271, -0.032672 },
            {  5.666837, -3.316225,  0.715455 },
            {  6.781360, -3.795765, -0.196864 },
            { -6.187331,  0.042326, -0.046499 },
            { -6.739113, -0.736038, -0.152558 },
            { -5.980298,  0.064208,  0.890635 },
            { -3.111868, -0.045870, -0.046351 },
            { -3.369378, -0.088582,  0.877460 },
            { -2.701406,  0.816674, -0.141741 },
            { -0.028600, -0.027927, -0.052228 },
            { -0.450128,  0.284536,  0.751664 },
            {  0.904097,  0.158746,  0.077351 },
            {  3.043640,  0.022778, -0.015695 },
            {  3.263633, -0.820206,  0.387469 },
            {  3.890123,  0.458652, -0.138341 },
            {  6.198501,  0.054145, -0.036443 },
            {  5.524424, -0.562672,  0.258040 },
            {  7.017693, -0.296767,  0.320414 },
            { -6.142883,  3.087172,  0.001701 },
            { -6.832739,  2.822585,  0.614608 },
            { -6.592184,  3.640199, -0.641611 },
            { -3.054590,  3.097169, -0.043753 },
            { -3.795838,  2.503772,  0.097642 },
            { -3.184109,  3.800316,  0.596849 },
            {  0.025258,  3.140300,  0.048161 },
            { -0.898451,  3.040556, -0.193447 },
            {  0.497539,  2.578923, -0.571006 },
            {  3.142639,  3.155831,  0.004528 },
            {  3.281734,  2.283294,  0.379885 },
            {  2.300622,  3.089671, -0.451752 },
            {  6.271003,  3.089878, -0.000229 },
            {  5.555786,  2.545332, -0.337120 },
            {  5.835504,  3.874496,  0.340762 },
            { -6.186874,  6.153181, -0.032864 },
            { -6.458269,  6.211097,  0.886130 },
            { -6.268404,  7.050370, -0.364476 },
            { -3.063571,  6.194209, -0.050917 },
            { -2.845329,  6.642761,  0.769286 },
            { -3.992058,  5.967481,  0.038907 },
            { -0.005511,  6.203833,  0.065479 },
            { -0.675428,  5.996461, -0.590089 },
            {  0.762902,  6.461025, -0.449259 },
            {  3.107182,  6.260055, -0.038149 },
            {  3.570052,  6.091314,  0.785768 },
            {  2.575104,  5.473755, -0.180227 },
            {  6.264803,  6.186329,  0.022994 },
            {  5.538470,  5.646253, -0.296882 },
            {  5.951245,  7.089071, -0.068103 },
            { -6.267283, -6.211735,  3.078419 },
            { -5.985306, -6.384271,  3.979686 },
            { -5.465039, -5.947785,  2.622023 },
            { -3.107839, -6.241119,  3.047498 },
            { -2.699991, -6.510712,  3.873658 },
            { -3.434740, -5.354930,  3.218863 },
            { -0.033434, -6.169301,  3.060705 },
            {  0.838698, -6.005633,  3.426984 },
            { -0.308002, -6.999962,  3.455898 },
            {  3.150365, -6.250863,  3.115724 },
            {  2.771880, -5.609654,  3.721663 },
            {  2.687846, -6.101322,  2.287910 },
            {  6.204347, -6.213550,  3.168674 },
            {  5.750267, -6.738865,  2.505789 },
            {  6.699059, -5.564375,  2.663310 },
            { -6.172325, -3.107328,  3.047333 },
            { -6.828965, -2.449862,  3.288422 },
            { -6.128632, -3.692984,  3.806714 },
            { -3.080296, -3.042262,  3.111628 },
            { -3.599745, -3.563632,  3.727994 },
            { -2.972174, -3.611996,  2.346596 },
            {  0.018427, -3.040709,  3.111508 },
            { -0.861291, -3.412145,  3.209954 },
            {  0.568794, -3.788138,  2.866545 },
            {  3.076991, -3.074154,  3.156063 },
            {  3.816301, -3.685968,  3.130285 },
            {  2.808084, -2.983447,  2.238990 },
            {  6.202874, -3.044570,  3.132656 },
            {  5.483336, -3.645941,  2.927282 },
            {  6.989368, -3.493060,  2.813525 },
            { -6.185279, -0.059803,  3.120451 },
            { -6.419455,  0.669043,  3.699681 },
            { -6.332522,  0.280211,  2.234860 },
            { -3.050408,  0.039017,  3.105488 },
            { -3.468427, -0.425632,  3.834133 },
            { -3.577897, -0.193686,  2.337916 },
            {  0.036504, -0.020609,  3.153990 },
            {  0.239664, -0.086591,  2.218077 },
            { -0.819081,  0.413707,  3.184104 },
            {  3.091492,  0.008311,  3.039147 },
            {  2.481020, -0.295977,  3.714642 },
            {  3.913193,  0.164051,  3.510436 },
            {  6.191365, -0.063471,  3.111335 },
            {  6.055864,  0.479518,  2.331358 },
            {  6.599523,  0.527950,  3.747880 },
            { -6.205614,  3.054919,  3.058779 },
            { -6.876203,  3.731733,  3.176233 },
            { -5.553000,  3.242989,  3.737223 },
            { -3.118412,  3.069232,  3.157498 },
            { -3.646030,  3.741133,  2.719629 },
            { -2.320885,  3.006400,  2.626871 },
            {  0.024980,  3.055432,  3.065769 },
            { -0.877220,  3.144581,  3.381441 },
            {  0.480709,  3.822007,  3.421065 },
            {  3.070975,  3.105245,  3.160550 },
            {  3.950017,  3.440960,  2.970440 },
            {  2.769856,  2.734952,  2.327620 },
            {  6.198446,  3.079479,  3.164063 },
            {  7.022878,  3.306229,  2.727648 },
            {  5.520109,  3.278660,  2.514642 },
            { -6.193602,  6.242270,  3.157025 },
            { -5.567869,  5.886374,  2.521920 },
            { -7.052009,  5.960994,  2.832091 },
            { -3.103387,  6.145098,  3.080057 },
            { -2.343780,  6.698458,  3.275890 },
            { -3.861615,  6.691317,  3.299822 },
            { -0.045554,  6.242642,  3.134180 },
            {  0.637737,  6.548568,  2.533277 },
            {  0.085332,  5.292898,  3.183344 },
            {  3.128530,  6.249892,  3.145938 },
            {  3.572612,  5.823157,  2.409562 },
            {  2.233691,  5.903232,  3.120423 },
            {  6.255670,  6.195132,  3.061028 },
            {  5.552429,  5.596136,  3.322210 },
            {  6.082243,  6.999450,  3.555555 },
            { -6.208072, -6.163006,  6.157412 },
            { -6.293854, -5.991009,  7.097954 },
            { -6.096342, -7.114505,  6.096361 },
            { -3.093215, -6.194828,  6.270532 },
            { -2.561298, -5.900834,  5.527435 },
            { -3.805554, -6.699586,  5.871339 },
            {  0.028309, -6.255582,  6.240192 },
            { -0.704628, -5.700497,  6.516318 },
            {  0.255286, -5.935581,  5.364040 },
            {  3.112376, -6.185324,  6.145057 },
            {  3.768593, -6.549569,  6.743611 },
            {  2.294122, -6.201701,  6.646808 },
            {  6.221282, -6.173949,  6.151454 },
            {  6.616602, -6.986255,  6.476173 },
            {  5.563911, -5.945569,  6.812716 },
            { -6.214945, -3.105559,  6.141418 },
            { -6.767312, -2.662372,  6.789520 },
            { -5.513791, -3.508545,  6.658668 },
            { -3.134984, -3.050634,  6.183218 },
            { -2.194062, -3.145360,  6.348303 },
            { -3.509793, -3.897382,  6.436395 },
            {  0.019709, -3.066046,  6.156136 },
            { -0.062167, -2.811154,  7.078030 },
            { -0.250671, -3.986954,  6.136538 },
            {  3.047052, -3.093053,  6.174132 },
            {  3.846770, -3.511425,  5.847060 },
            {  3.252831, -2.858011,  7.081861 },
            {  6.264526, -3.131008,  6.192287 },
            {  6.014957, -2.208761,  6.098811 },
            {  5.479151, -3.558207,  6.541928 },
            { -6.202784,  0.012060,  6.271161 },
            { -5.794767, -0.706266,  5.782190 },
            { -6.679360,  0.514844,  5.606597 },
            { -3.132137,  0.011406,  6.148615 },
            { -3.531033, -0.354855,  6.941266 },
            { -2.218018,  0.173816,  6.392684 },
            { -0.048817, -0.044068,  6.207411 },
            {  0.264785,  0.472342,  5.461397 },
            {  0.510082,  0.227141,  6.939285 },
            {  3.104615, -0.056762,  6.239817 },
            {  2.334121,  0.441331,  5.957323 },
            {  3.851786,  0.459650,  5.928986 },
            {  6.191962, -0.011674,  6.269664 },
            {  7.056171,  0.162434,  5.889678 },
            {  5.589740,  0.022864,  5.522872 },
            { -6.227582,  3.048278,  6.177333 },
            { -5.456950,  3.575652,  5.954701 },
            { -6.623561,  3.504488,  6.923412 },
            { -3.095950,  3.167462,  6.219110 },
            { -3.718176,  2.753599,  5.616521 },
            { -2.605272,  2.434738,  6.598467 },
            { -0.028196,  3.102085,  6.266027 },
            {  0.892169,  3.276065,  6.055742 },
            { -0.444613,  2.949998,  5.414541 },
            {  3.122055,  3.044925,  6.230400 },
            {  2.318548,  3.530579,  6.430594 },
            {  3.590528,  3.602778,  5.605184 },
            {  6.235943,  3.056340,  6.242112 },
            {  5.925321,  3.913199,  6.543556 },
            {  6.022471,  3.038970,  5.306329 },
            { -6.152078,  6.212531,  6.242895 },
            { -6.277517,  6.463813,  5.324917 },
            { -7.001466,  5.855602,  6.512529 },
            { -3.075211,  6.151648,  6.229296 },
            { -3.984620,  6.273567,  5.947033 },
            { -2.668015,  7.012240,  6.106276 },
            {  0.047399,  6.200405,  6.251791 },
            { -0.389422,  5.502373,  5.758333 },
            { -0.362937,  7.009521,  5.937917 },
            {  3.121019,  6.156392,  6.245211 },
            {  3.669300,  6.880502,  5.934383 },
            {  2.256232,  6.330001,  5.866310 },
            {  6.204797,  6.271399,  6.195395 },
            {  5.469074,  5.654753,  6.199970 },
            {  6.973096,  5.730251,  6.391449 }
        };
        for(int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }else{
        throw exception();
    }

    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLength, 0, 0),
                                        Vec3(0, boxEdgeLength, 0),
                                        Vec3(0, 0, boxEdgeLength));

    sig.clear();
    eps.clear();
    bonds.clear();
    for(int atom = 0; atom < natoms; ++atom){
        system.addParticle(masses[atom%RESSIZE]);
        double sigma = 2.0*pow(2.0, -1.0/6.0)*halfrmins[atom%RESSIZE]*OpenMM::NmPerAngstrom;
        double epsilon = fabs(epsilons[atom%RESSIZE])*OpenMM::KJPerKcal;
        sig.push_back(0.5*sigma);
        eps.push_back(2.0*sqrt(epsilon));
        if(atom%RESSIZE == 0){
            bonds.push_back(pair<int, int>(atom, atom+1));
            bonds.push_back(pair<int, int>(atom, atom+2));
        }
        double charge = do_electrostatics ? charges[atom] : 0;
        forceField->addParticle(charge, sigma, epsilon);
    }
}

void print_forces(const vector<Vec3>& forces){
    // Print the forces in AKMA units, to compare against other codes.
    std::cout << "Forces:" << std::endl;
    for(int n = 0; n < forces.size(); ++ n){
        std::cout << setw(3)<< n+1
                  << setw(16) << setprecision(10) << -forces[n][0]*OpenMM::NmPerAngstrom*OpenMM::KcalPerKJ
                  << setw(16) << setprecision(10) << -forces[n][1]*OpenMM::NmPerAngstrom*OpenMM::KcalPerKJ
                  << setw(16) << setprecision(10) << -forces[n][2]*OpenMM::NmPerAngstrom*OpenMM::KcalPerKJ << std::endl;
    }
}

void test_water125_dpme_vs_long_cutoff_with_exclusions() {
    const double cutoff = 8.0*OpenMM::NmPerAngstrom;
    const double alpha = 0.45*OpenMM::AngstromsPerNm;
    const double dalpha = 0.45*OpenMM::AngstromsPerNm;
    const int grid = 32;
    DPMENonbondedForce* forceField = new DPMENonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int NATOMS = 375;
    double boxEdgeLength = 16.01*OpenMM::NmPerAngstrom;

    make_waterbox(NATOMS, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->createExceptionsFromBonds(bonds, 1.0, 1.0);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setDispersionPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
//    std::cout << "Energy " << energy*OpenMM::KcalPerKJ << std::endl;
    const vector<Vec3>& forces = state.getForces();
//    print_forces(forces);


    // Find the exclusion information
    vector<set<int> > exclusions(NATOMS);
    for (int i = 0; i < forceField->getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        forceField->getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    const double longcutoff = 30.0*OpenMM::NmPerAngstrom;
    const double longcutoff2 = longcutoff*longcutoff;
    const double cutoff2 = cutoff*cutoff;
    const double cutoff6inv = 1.0 / (cutoff2*cutoff2*cutoff2);
    const int nboxes = ceil(longcutoff/boxEdgeLength);

    double refenergy = 0.0;
    // Loop over home box first...
    for(int i = 0; i < NATOMS; ++ i){
        for(int j = i+1; j < NATOMS; ++j){
            Vec3 dR = positions[i] - positions[j];
            double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
            double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
            double sig6 = sig2*sig2*sig2;
            double eps = epsvals[i]*epsvals[j];
            refenergy += 2.0*eps*(sig6-1.0)*sig6;
            if(R2 < cutoff2){
                // Add a shift term for direct space parts withing t
                refenergy += 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
            }
        }
    }

    // ... back out exclusions
    for (int ii = 0; ii < NATOMS; ii++){
        for (set<int>::const_iterator iter = exclusions[ii].begin(); iter != exclusions[ii].end(); ++iter) {
            if (*iter > ii) {
                int i = ii;
                int j = *iter;
                Vec3 dR = positions[i] - positions[j];
                double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
                double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
                double sig6 = sig2*sig2*sig2;
                double eps = epsvals[i]*epsvals[j];
                refenergy -= 2.0*eps*(sig6-1.0)*sig6;
                if(R2 < cutoff2){
                    // Add a shift term for direct space parts withing t
                    refenergy -= 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                }
            }
        }
    }

    // ... and now add in the image box terms
    for(int bx = -nboxes; bx <= nboxes; ++bx){
        for(int by = -nboxes; by <= nboxes; ++by){
            for(int bz = -nboxes; bz <= nboxes; ++bz){
                if(bx==0 && by==0 && bz==0) continue;
                Vec3 offset(bx*boxEdgeLength, by*boxEdgeLength, bz*boxEdgeLength);
                for(int i = 0; i < NATOMS; ++ i){
                    for(int j = 0; j < NATOMS; ++j){
                        Vec3 dR = positions[i] - positions[j] + offset;
                        double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
                        if(R2 > longcutoff2) continue;
                        double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
                        double sig6 = sig2*sig2*sig2;
                        double eps = epsvals[i]*epsvals[j];
                        refenergy += eps*(sig6-1.0)*sig6;
                        if(R2 < cutoff2){
                            // Add a shift term for direct space parts withing teh
                            refenergy += eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                        }
                    }
                }
            }
        }
    }
    refenergy *= 0.5;

    // For this test the reference energy is -394.662 kJ/mol, while the difference between DPME and 30A cutoffs
    // is just 0.064 kJ/mol.  The difference is due to the fact that arithmetic mean combination rules are used
    // up to the cutoff, while the reciprocal space uses the geometric mean.  See DOI: 10.1021/acs.jctc.5b00726
    ASSERT_EQUAL_TOL(refenergy, energy, 2E-4);
}


void test_water125_dpme_vs_long_cutoff_no_exclusions() {
    const double cutoff = 8.0*OpenMM::NmPerAngstrom;
    const double alpha = 0.45*OpenMM::AngstromsPerNm;
    const double dalpha = 0.45*OpenMM::AngstromsPerNm;
    const int grid = 32;
    DPMENonbondedForce* forceField = new DPMENonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int NATOMS = 375;
    double boxEdgeLength = 16.01*OpenMM::NmPerAngstrom;

    make_waterbox(NATOMS, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setDispersionPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
//    std::cout << "Energy " << energy*OpenMM::KcalPerKJ << std::endl;
    const vector<Vec3>& forces = state.getForces();
//    print_forces(forces);


    const double longcutoff = 30.0*OpenMM::NmPerAngstrom;
    const double longcutoff2 = longcutoff*longcutoff;
    const double cutoff2 = cutoff*cutoff;
    const double cutoff6inv = 1.0 / (cutoff2*cutoff2*cutoff2);
    const int nboxes = ceil(longcutoff/boxEdgeLength);

    double refenergy = 0.0;
    // Loop over home box first...
    for(int i = 0; i < NATOMS; ++ i){
        for(int j = i+1; j < NATOMS; ++j){
            Vec3 dR = positions[i] - positions[j];
            double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
            double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
            double sig6 = sig2*sig2*sig2;
            double eps = epsvals[i]*epsvals[j];
            refenergy += 2.0*eps*(sig6-1.0)*sig6;
            if(R2 < cutoff2){
                // Add a shift term for direct space parts withing t
                refenergy += 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
            }
        }
    }

    // ... and now add in the image box terms
    for(int bx = -nboxes; bx <= nboxes; ++bx){
        for(int by = -nboxes; by <= nboxes; ++by){
            for(int bz = -nboxes; bz <= nboxes; ++bz){
                if(bx==0 && by==0 && bz==0) continue;
                Vec3 offset(bx*boxEdgeLength, by*boxEdgeLength, bz*boxEdgeLength);
                for(int i = 0; i < NATOMS; ++ i){
                    for(int j = 0; j < NATOMS; ++j){
                        Vec3 dR = positions[i] - positions[j] + offset;
                        double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
                        if(R2 > longcutoff2) continue;
                        double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
                        double sig6 = sig2*sig2*sig2;
                        double eps = epsvals[i]*epsvals[j];
                        refenergy += eps*(sig6-1.0)*sig6;
                        if(R2 < cutoff2){
                            // Add a shift term for direct space parts withing teh
                            refenergy += eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                        }
                    }
                }
            }
        }
    }
    refenergy *= 0.5;

    // For this test the reference energy is 545537 kJ/mol, while the difference between DPME and 30A cutoffs
    // is just 1.624 kJ/mol.  The difference is due to the fact that arithmetic mean combination rules are used
    // up to the cutoff, while the reciprocal space uses the geometric mean.  See DOI: 10.1021/acs.jctc.5b00726
    ASSERT_EQUAL_TOL(refenergy, energy, 5E-6);
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);

        // Tests of energy; long cutoff LJ vs. DPME with short cutoffs.
        test_water125_dpme_vs_long_cutoff_with_exclusions();
        test_water125_dpme_vs_long_cutoff_no_exclusions();

        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
