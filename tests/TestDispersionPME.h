/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2017 Stanford University and the Authors.      *
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
#ifdef _MSC_VER
    // Prevent Windows from defining macros that interfere with other code.
    #define NOMINMAX
#endif

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/Units.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace std;

void testConvergence() {
    // Create a cloud of random particles.

    const int numParticles = 1000;
    const double boxWidth = 6.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(0.0, 0.1+0.3*genrand_real2(sfmt), genrand_real2(sfmt));
        while (true) {
            Vec3 pos = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
            double minDist = boxWidth;
            for (int j = 0; j < i; j++) {
                Vec3 delta = pos-positions[j];
                minDist = std::min(minDist, sqrt(delta.dot(delta)));
            }
            if (minDist > 0.15) {
                positions[i] = pos;
                break;
            }
        }
    }
    
    // Compute the energy with short and long cutoffs, and compare them to LJPME.
    // The long cutoff should match the LJPME result better than the short cutoff.
    
    force->setNonbondedMethod(NonbondedForce::LJPME);
    force->setCutoffDistance(1.0);
    force->setEwaldErrorTolerance(1e-5);
    force->setUseDispersionCorrection(false);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    double pmeEnergy = context.getState(State::Energy, false, 1).getPotentialEnergy();
    force->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    context.reinitialize();
    context.setPositions(positions);
    double shortEnergy = context.getState(State::Energy, false, 1).getPotentialEnergy();
    force->setCutoffDistance(3.0);
    context.reinitialize();
    context.setPositions(positions);
    double longEnergy = context.getState(State::Energy, false, 1).getPotentialEnergy();
    ASSERT(fabs(longEnergy-pmeEnergy) < fabs(shortEnergy-pmeEnergy))
}

void testErrorTolerance() {
    // Create a cloud of random particles.

    const int numParticles = 200;
    const double boxWidth = 5.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(0.0, 0.1+0.2*genrand_real2(sfmt), genrand_real2(sfmt));
        while (true) {
            Vec3 pos = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
            double minDist = boxWidth;
            for (int j = 0; j < i; j++) {
                Vec3 delta = pos-positions[j];
                minDist = std::min(minDist, sqrt(delta.dot(delta)));
            }
            if (minDist > 0.1) {
                positions[i] = pos;
                break;
            }
        }
    }
    force->setNonbondedMethod(NonbondedForce::LJPME);

    // For various values of the cutoff and error tolerance, see if the actual error is reasonable.

    for (double cutoff = 0.8; cutoff < boxWidth/2; cutoff *= 1.3) {
        force->setCutoffDistance(cutoff);
        vector<Vec3> refForces;
        double norm = 0.0;
        for (double tol = 5e-5; tol < 1e-3; tol *= 2.0) {
            force->setEwaldErrorTolerance(tol);
            VerletIntegrator integrator(0.01);
            Context context(system, integrator, platform);
            context.setPositions(positions);
            State state = context.getState(State::Forces);
            if (refForces.size() == 0) {
                refForces = state.getForces();
                for (int i = 0; i < numParticles; i++)
                    norm += refForces[i].dot(refForces[i]);
                norm = sqrt(norm);
            }
            else {
                double diff = 0.0;
                for (int i = 0; i < numParticles; i++) {
                    Vec3 delta = refForces[i]-state.getForces()[i];
                    diff += delta.dot(delta);
                }
                diff = sqrt(diff)/norm;
                ASSERT(diff < 2*tol);
            }

            // See if the PME parameters were calculated correctly.

            double expectedAlpha, actualAlpha;
            int expectedSize[3], actualSize[3];
            NonbondedForceImpl::calcPMEParameters(system, *force, expectedAlpha, expectedSize[0], expectedSize[1], expectedSize[2], true);
            force->getLJPMEParametersInContext(context, actualAlpha, actualSize[0], actualSize[1], actualSize[2]);
            ASSERT_EQUAL_TOL(expectedAlpha, actualAlpha, 1e-5);
            for (int i = 0; i < 3; i++) {
                ASSERT(actualSize[i] >= expectedSize[i]);
                ASSERT(actualSize[i] < expectedSize[i]+10);
            }
        }
    }
}

void testPMEParameters() {
    // Create a cloud of random particles.

    const int numParticles = 51;
    const double boxWidth = 4.7;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(0.0, 0.1+0.2*genrand_real2(sfmt), genrand_real2(sfmt));
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::LJPME);
    force->setCutoffDistance(0.5);
    
    // Compute the energy with an error tolerance of 0.1.

    force->setEwaldErrorTolerance(0.1);
    VerletIntegrator integrator1(0.01);
    Context context1(system, integrator1, platform);
    context1.setPositions(positions);
    double energy1 = context1.getState(State::Energy).getPotentialEnergy();
    double alpha;
    int gridx, gridy, gridz;
    force->getLJPMEParametersInContext(context1, alpha, gridx, gridy, gridz);
    
    // Try again with an error tolerance of 1e-4.

    force->setEwaldErrorTolerance(1e-4);
    VerletIntegrator integrator2(0.01);
    Context context2(system, integrator2, platform);
    context2.setPositions(positions);
    double energy2 = context2.getState(State::Energy).getPotentialEnergy();
    
    // Now explicitly set the parameters.  These should match the values that were
    // used for tolerance 0.1.

    force->setLJPMEParameters(alpha, gridx, gridy, gridz);
    VerletIntegrator integrator3(0.01);
    Context context3(system, integrator3, platform);
    context3.setPositions(positions);
    double energy3 = context3.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(energy1, energy3, 1e-5);
    ASSERT(fabs((energy1-energy2)/energy1) > 1e-5);
    force->getLJPMEParametersInContext(context2, alpha, gridx, gridy, gridz);
}

void testCoulombAndLJ() {
    // Create a cloud of random particles.

    const int numParticles = 200;
    const double boxWidth = 5.0;
    System system;
    system.setDefaultPeriodicBoxVectors(Vec3(boxWidth, 0, 0), Vec3(0, boxWidth, 0), Vec3(0, 0, boxWidth));
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    vector<double> charge(numParticles), sigma(numParticles), epsilon(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        charge[i] = -1.0+i*2.0/(numParticles-1);
        sigma[i] = 0.1+0.2*genrand_real2(sfmt);
        epsilon[i] = genrand_real2(sfmt);
        force->addParticle(charge[i], 1.0, 0.0);
        while (true) {
            Vec3 pos = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
            double minDist = boxWidth;
            for (int j = 0; j < i; j++) {
                Vec3 delta = pos-positions[j];
                minDist = std::min(minDist, sqrt(delta.dot(delta)));
            }
            if (minDist > 0.1) {
                positions[i] = pos;
                break;
            }
        }
    }
    force->setNonbondedMethod(NonbondedForce::LJPME);
    
    // Compute forces and energy with only Coulomb interactions.
    
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state1 = context.getState(State::Forces | State::Energy);
    
    // Now repeat with only LJ interactions.
    
    for (int i = 0; i < numParticles; i++)
        force->setParticleParameters(i, 0.0, sigma[i], epsilon[i]);
    context.reinitialize();
    context.setPositions(positions);
    State state2 = context.getState(State::Forces | State::Energy);
    
    // Finally compute with both Coulomb and LJ.
    
    for (int i = 0; i < numParticles; i++)
        force->setParticleParameters(i, charge[i], sigma[i], epsilon[i]);
    context.reinitialize();
    context.setPositions(positions);
    State state3 = context.getState(State::Forces | State::Energy);
    
    // Make sure the results agree.
    
    ASSERT_EQUAL_TOL(state1.getPotentialEnergy()+state2.getPotentialEnergy(), state3.getPotentialEnergy(), 1e-5);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(state1.getForces()[i]+state2.getForces()[i], state3.getForces()[i], 1e-5);
}

void make_dmfbox(int natoms, double boxEdgeLength, NonbondedForce *forceField,  vector<Vec3> &positions, vector<double>& eps, vector<double>& sig,
                   vector<pair<int, int> >& bonds, System &system, bool do_electrostatics) {

    const int RESSIZE = 12;
    const double charges[RESSIZE] = {
        0.08, 0.43, -0.54, -0.33, -0.09, 0.09,
        0.09, 0.09, -0.09, 0.09, 0.09, 0.09
    };
    const double masses[RESSIZE] = {
        1.008, 12.011, 15.9994, 14.007, 12.011, 1.008,
        1.008, 1.008, 12.011, 1.008, 1.008, 1.008
    };
    const double sigmas[RESSIZE] = {
        0.160361769265, 0.356359487256, 0.302905564168, 0.329632525712, 0.365268474438, 0.238760856462,
        0.238760856462, 0.238760856462, 0.365268474438, 0.238760856462, 0.238760856462, 0.238760856462
    };
    const double epsilons[RESSIZE] = {
        0.19246, 0.46024, 0.50208, 0.83680, 0.32635, 0.10042,
        0.10042, 0.10042, 0.32635, 0.10042, 0.10042, 0.10042
    };
    positions.clear();
    if (natoms == 12) {
        const double coords[12][3] = {
            {  0.620,  0.945,  0.128 },
            {  0.560,  0.905,  0.045 },
            {  0.502,  0.987, -0.026 },
            {  0.556,  0.770,  0.028 },
            {  0.643,  0.685,  0.106 },
            {  0.680,  0.741,  0.193 },
            {  0.729,  0.653,  0.045 },
            {  0.587,  0.597,  0.140 },
            {  0.477,  0.708, -0.072 },
            {  0.377,  0.688, -0.032 },
            {  0.524,  0.614, -0.102 },
            {  0.468,  0.775, -0.158 }
        };
        for (int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2]));
    }
    else
        throw exception();

    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLength, 0, 0),
                                        Vec3(0, boxEdgeLength, 0),
                                        Vec3(0, 0, boxEdgeLength));
    sig.clear();
    eps.clear();
    bonds.clear();
    for (int atom = 0; atom < natoms; ++atom) {
        system.addParticle(masses[atom%RESSIZE]);
        double sigma = sigmas[atom%RESSIZE];
        double epsilon = epsilons[atom%RESSIZE];
        sig.push_back(sigma);
        eps.push_back(epsilon);
        if (atom%RESSIZE == 0) {
            int offset = atom-1;
            bonds.push_back(pair<int, int>(offset+1, offset+ 2));
            bonds.push_back(pair<int, int>(offset+2, offset+ 3));
            bonds.push_back(pair<int, int>(offset+2, offset+ 4));
            bonds.push_back(pair<int, int>(offset+4, offset+ 5));
            bonds.push_back(pair<int, int>(offset+4, offset+ 9));
            bonds.push_back(pair<int, int>(offset+5, offset+ 6));
            bonds.push_back(pair<int, int>(offset+5, offset+ 7));
            bonds.push_back(pair<int, int>(offset+5, offset+ 8));
            bonds.push_back(pair<int, int>(offset+9, offset+10));
            bonds.push_back(pair<int, int>(offset+9, offset+11));
            bonds.push_back(pair<int, int>(offset+9, offset+12));
        }
        double charge = do_electrostatics ? charges[atom] : 0;
        forceField->addParticle(charge, sigma, epsilon);
    }
    // Make exceptions and replacement 1-4 parameters
    forceField->createExceptionsFromBonds(bonds, 1.0, 1.0);
    int nres = natoms / RESSIZE;
    double newqq = do_electrostatics ? charges[8]*charges[8] : 0.0;
    double newsig = 0.293996576986;
    double neweps = 0.144938011577;
    for(int i = 0; i < forceField->getNumExceptions(); ++i){
        for(int res = 0; res < nres; ++res){
            int particle1, particle2;
            double chargeProd, sigma, epsilon;
            forceField->getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
            if((particle1 == res+2 && particle2 == res+8) ||
               (particle1 == res+8 && particle2 == res+2) ||
               (particle1 == res+2 && particle2 == res+4) ||
               (particle1 == res+4 && particle2 == res+2))
                forceField->setExceptionParameters(i, particle1, particle2, newqq, newsig, neweps);
        }
    }
}


void make_waterbox(int natoms, double boxEdgeLength, NonbondedForce *forceField,  vector<Vec3> &positions, vector<double>& eps, vector<double>& sig,
                   vector<pair<int, int> >& bonds, System &system, bool do_electrostatics) {
    const int RESSIZE = 3;
    const double masses[RESSIZE]    = {     8.0,    1.0,    1.0 };
    const double charges[RESSIZE]   = {  -0.834,  0.417,  0.417 };
    // Values from the CHARMM force field, in AKMA units
    const double epsilons[RESSIZE]  = { -0.1521, -0.046, -0.046 };
    const double halfrmins[RESSIZE] = {  1.7682, 0.2245, 0.2245 };
    positions.clear();
    if (natoms == 6) {
        const double coords[6][3] = {
            {  2.000000, 2.000000, 2.000000},
            {  2.500000, 2.000000, 3.000000},
            {  1.500000, 2.000000, 3.000000},
            {  0.000000, 0.000000, 0.000000},
            {  0.500000, 0.000000, 1.000000},
            { -0.500000, 0.000000, 1.000000}
        };
        for (int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }
    else if (natoms == 375) {
        const double coords[375][3] = {
            { -6.22, -6.25, -6.24 },
            { -5.32, -6.03, -6.00 },
            { -6.75, -5.56, -5.84 },
            { -3.04, -6.23, -6.19 },
            { -3.52, -5.55, -5.71 },
            { -3.59, -6.43, -6.94 },
            {  0.02, -6.23, -6.14 },
            { -0.87, -5.97, -6.37 },
            {  0.53, -6.03, -6.93 },
            {  3.10, -6.20, -6.27 },
            {  3.87, -6.35, -5.72 },
            {  2.37, -6.11, -5.64 },
            {  6.18, -6.14, -6.20 },
            {  6.46, -6.66, -5.44 },
            {  6.26, -6.74, -6.94 },
            { -6.21, -3.15, -6.24 },
            { -6.23, -3.07, -5.28 },
            { -6.02, -2.26, -6.55 },
            { -3.14, -3.07, -6.16 },
            { -3.38, -3.63, -6.90 },
            { -2.18, -3.05, -6.17 },
            { -0.00, -3.16, -6.23 },
            { -0.03, -2.30, -6.67 },
            {  0.05, -2.95, -5.29 },
            {  3.08, -3.11, -6.14 },
            {  2.65, -2.55, -6.79 },
            {  3.80, -3.53, -6.62 },
            {  6.16, -3.14, -6.16 },
            {  7.04, -3.32, -6.51 },
            {  5.95, -2.27, -6.51 },
            { -6.20, -0.04, -6.15 },
            { -5.43,  0.32, -6.59 },
            { -6.95,  0.33, -6.62 },
            { -3.10, -0.06, -6.19 },
            { -3.75,  0.42, -6.69 },
            { -2.46,  0.60, -5.93 },
            {  0.05, -0.01, -6.17 },
            { -0.10,  0.02, -7.12 },
            { -0.79,  0.16, -5.77 },
            {  3.03,  0.00, -6.19 },
            {  3.54,  0.08, -7.01 },
            {  3.69, -0.22, -5.53 },
            {  6.17,  0.05, -6.19 },
            {  5.78, -0.73, -6.57 },
            {  7.09, -0.17, -6.04 },
            { -6.20,  3.15, -6.25 },
            { -6.59,  3.18, -5.37 },
            { -5.87,  2.25, -6.33 },
            { -3.09,  3.04, -6.17 },
            { -3.88,  3.58, -6.26 },
            { -2.41,  3.54, -6.63 },
            {  0.00,  3.06, -6.26 },
            { -0.71,  3.64, -6.00 },
            {  0.65,  3.15, -5.55 },
            {  3.14,  3.06, -6.23 },
            {  3.11,  3.31, -5.30 },
            {  2.38,  3.49, -6.63 },
            {  6.19,  3.14, -6.25 },
            {  6.82,  3.25, -5.54 },
            {  5.76,  2.30, -6.07 },
            { -6.22,  6.26, -6.19 },
            { -6.22,  5.74, -7.00 },
            { -5.89,  5.67, -5.52 },
            { -3.04,  6.24, -6.20 },
            { -3.08,  5.28, -6.17 },
            { -3.96,  6.52, -6.25 },
            { -0.05,  6.21, -6.16 },
            {  0.82,  6.58, -6.06 },
            {  0.01,  5.64, -6.93 },
            {  3.10,  6.25, -6.15 },
            {  3.64,  5.47, -6.31 },
            {  2.46,  6.24, -6.87 },
            {  6.22,  6.20, -6.27 },
            {  5.37,  6.42, -5.88 },
            {  6.80,  6.07, -5.51 },
            { -6.19, -6.15, -3.13 },
            { -6.37, -7.01, -3.51 },
            { -6.25, -6.29, -2.18 },
            { -3.10, -6.27, -3.11 },
            { -2.29, -5.77, -2.99 },
            { -3.80, -5.62, -2.98 },
            { -0.03, -6.18, -3.15 },
            { -0.07, -7.05, -2.75 },
            {  0.68, -5.74, -2.70 },
            {  3.10, -6.14, -3.07 },
            {  2.35, -6.72, -3.23 },
            {  3.86, -6.65, -3.37 },
            {  6.22, -6.20, -3.16 },
            {  6.82, -6.36, -2.43 },
            {  5.35, -6.13, -2.75 },
            { -6.26, -3.13, -3.12 },
            { -6.16, -2.27, -2.70 },
            { -5.36, -3.47, -3.18 },
            { -3.11, -3.05, -3.14 },
            { -3.31, -3.96, -3.34 },
            { -2.77, -3.06, -2.24 },
            {  0.00, -3.13, -3.16 },
            {  0.48, -2.37, -2.81 },
            { -0.57, -3.40, -2.44 },
            {  3.09, -3.09, -3.16 },
            {  2.41, -3.19, -2.49 },
            {  3.91, -3.07, -2.67 },
            {  6.19, -3.04, -3.08 },
            {  5.64, -3.61, -3.61 },
            {  6.93, -3.58, -2.82 },
            { -6.18, -0.00, -3.04 },
            { -6.00, -0.59, -3.78 },
            { -6.79,  0.64, -3.39 },
            { -3.05, -0.03, -3.07 },
            { -2.95,  0.80, -3.52 },
            { -4.00, -0.20, -3.07 },
            { -0.03,  0.03, -3.06 },
            { -0.33, -0.37, -3.87 },
            {  0.89, -0.21, -2.99 },
            {  3.13, -0.05, -3.10 },
            {  3.44,  0.81, -3.34 },
            {  2.21,  0.07, -2.86 },
            {  6.20, -0.05, -3.13 },
            {  6.89,  0.60, -3.20 },
            {  5.58,  0.30, -2.49 },
            { -6.23,  3.09, -3.16 },
            { -5.62,  3.79, -2.94 },
            { -6.33,  2.60, -2.33 },
            { -3.10,  3.08, -3.04 },
            { -3.84,  3.47, -3.51 },
            { -2.40,  3.01, -3.69 },
            {  0.01,  3.04, -3.11 },
            { -0.56,  3.59, -3.64 },
            {  0.28,  3.60, -2.38 },
            {  3.04,  3.11, -3.09 },
            {  3.49,  2.30, -2.87 },
            {  3.70,  3.66, -3.51 },
            {  6.15,  3.14, -3.11 },
            {  6.52,  2.52, -3.74 },
            {  6.72,  3.06, -2.34 },
            { -6.22,  6.15, -3.13 },
            { -5.49,  6.21, -2.51 },
            { -6.56,  7.04, -3.18 },
            { -3.11,  6.24, -3.05 },
            { -3.76,  5.83, -3.62 },
            { -2.26,  5.92, -3.37 },
            {  0.03,  6.25, -3.07 },
            {  0.34,  5.63, -3.73 },
            { -0.87,  6.00, -2.91 },
            {  3.07,  6.15, -3.08 },
            {  3.29,  6.92, -2.56 },
            {  3.39,  6.35, -3.96 },
            {  6.22,  6.14, -3.12 },
            {  5.79,  6.38, -2.29 },
            {  6.25,  6.96, -3.62 },
            { -6.21, -6.20, -0.06 },
            { -5.79, -6.87,  0.48 },
            { -6.43, -5.50,  0.54 },
            { -3.16, -6.21, -0.02 },
            { -2.50, -6.87,  0.20 },
            { -2.77, -5.37,  0.23 },
            { -0.00, -6.14, -0.00 },
            {  0.68, -6.72, -0.33 },
            { -0.64, -6.73,  0.38 },
            {  3.03, -6.20, -0.01 },
            {  3.77, -6.56, -0.51 },
            {  3.43, -5.85,  0.78 },
            {  6.25, -6.16, -0.00 },
            {  5.36, -6.09, -0.36 },
            {  6.24, -6.97,  0.49 },
            { -6.24, -3.05, -0.01 },
            { -6.35, -3.64,  0.73 },
            { -5.42, -3.33, -0.42 },
            { -3.09, -3.06,  0.05 },
            { -2.44, -3.62, -0.38 },
            { -3.90, -3.21, -0.43 },
            {  0.05, -3.10,  0.02 },
            { -0.31, -2.35, -0.43 },
            { -0.63, -3.77,  0.01 },
            {  3.05, -3.09, -0.04 },
            {  3.28, -3.90,  0.41 },
            {  3.65, -2.43,  0.30 },
            {  6.20, -3.04, -0.03 },
            {  5.66, -3.31,  0.71 },
            {  6.78, -3.79, -0.19 },
            { -6.18,  0.04, -0.04 },
            { -6.73, -0.73, -0.15 },
            { -5.98,  0.06,  0.89 },
            { -3.11, -0.04, -0.04 },
            { -3.36, -0.08,  0.87 },
            { -2.70,  0.81, -0.14 },
            { -0.02, -0.02, -0.05 },
            { -0.45,  0.28,  0.75 },
            {  0.90,  0.15,  0.07 },
            {  3.04,  0.02, -0.01 },
            {  3.26, -0.82,  0.38 },
            {  3.89,  0.45, -0.13 },
            {  6.19,  0.05, -0.03 },
            {  5.52, -0.56,  0.25 },
            {  7.01, -0.29,  0.32 },
            { -6.14,  3.08,  0.00 },
            { -6.83,  2.82,  0.61 },
            { -6.59,  3.64, -0.64 },
            { -3.05,  3.09, -0.04 },
            { -3.79,  2.50,  0.09 },
            { -3.18,  3.80,  0.59 },
            {  0.02,  3.14,  0.04 },
            { -0.89,  3.04, -0.19 },
            {  0.49,  2.57, -0.57 },
            {  3.14,  3.15,  0.00 },
            {  3.28,  2.28,  0.37 },
            {  2.30,  3.08, -0.45 },
            {  6.27,  3.08, -0.00 },
            {  5.55,  2.54, -0.33 },
            {  5.83,  3.87,  0.34 },
            { -6.18,  6.15, -0.03 },
            { -6.45,  6.21,  0.88 },
            { -6.26,  7.05, -0.36 },
            { -3.06,  6.19, -0.05 },
            { -2.84,  6.64,  0.76 },
            { -3.99,  5.96,  0.03 },
            { -0.00,  6.20,  0.06 },
            { -0.67,  5.99, -0.59 },
            {  0.76,  6.46, -0.44 },
            {  3.10,  6.26, -0.03 },
            {  3.57,  6.09,  0.78 },
            {  2.57,  5.47, -0.18 },
            {  6.26,  6.18,  0.02 },
            {  5.53,  5.64, -0.29 },
            {  5.95,  7.08, -0.06 },
            { -6.26, -6.21,  3.07 },
            { -5.98, -6.38,  3.97 },
            { -5.46, -5.94,  2.62 },
            { -3.10, -6.24,  3.04 },
            { -2.69, -6.51,  3.87 },
            { -3.43, -5.35,  3.21 },
            { -0.03, -6.16,  3.06 },
            {  0.83, -6.00,  3.42 },
            { -0.30, -6.99,  3.45 },
            {  3.15, -6.25,  3.11 },
            {  2.77, -5.60,  3.72 },
            {  2.68, -6.10,  2.28 },
            {  6.20, -6.21,  3.16 },
            {  5.75, -6.73,  2.50 },
            {  6.69, -5.56,  2.66 },
            { -6.17, -3.10,  3.04 },
            { -6.82, -2.44,  3.28 },
            { -6.12, -3.69,  3.80 },
            { -3.08, -3.04,  3.11 },
            { -3.59, -3.56,  3.72 },
            { -2.97, -3.61,  2.34 },
            {  0.01, -3.04,  3.11 },
            { -0.86, -3.41,  3.20 },
            {  0.56, -3.78,  2.86 },
            {  3.07, -3.07,  3.15 },
            {  3.81, -3.68,  3.13 },
            {  2.80, -2.98,  2.23 },
            {  6.20, -3.04,  3.13 },
            {  5.48, -3.64,  2.92 },
            {  6.98, -3.49,  2.81 },
            { -6.18, -0.05,  3.12 },
            { -6.41,  0.66,  3.69 },
            { -6.33,  0.28,  2.23 },
            { -3.05,  0.03,  3.10 },
            { -3.46, -0.42,  3.83 },
            { -3.57, -0.19,  2.33 },
            {  0.03, -0.02,  3.15 },
            {  0.23, -0.08,  2.21 },
            { -0.81,  0.41,  3.18 },
            {  3.09,  0.00,  3.03 },
            {  2.48, -0.29,  3.71 },
            {  3.91,  0.16,  3.51 },
            {  6.19, -0.06,  3.11 },
            {  6.05,  0.47,  2.33 },
            {  6.59,  0.52,  3.74 },
            { -6.20,  3.05,  3.05 },
            { -6.87,  3.73,  3.17 },
            { -5.55,  3.24,  3.73 },
            { -3.11,  3.06,  3.15 },
            { -3.64,  3.74,  2.71 },
            { -2.32,  3.00,  2.62 },
            {  0.02,  3.05,  3.06 },
            { -0.87,  3.14,  3.38 },
            {  0.48,  3.82,  3.42 },
            {  3.07,  3.10,  3.16 },
            {  3.95,  3.44,  2.97 },
            {  2.76,  2.73,  2.32 },
            {  6.19,  3.07,  3.16 },
            {  7.02,  3.30,  2.72 },
            {  5.52,  3.27,  2.51 },
            { -6.19,  6.24,  3.15 },
            { -5.56,  5.88,  2.52 },
            { -7.05,  5.96,  2.83 },
            { -3.10,  6.14,  3.08 },
            { -2.34,  6.69,  3.27 },
            { -3.86,  6.69,  3.29 },
            { -0.04,  6.24,  3.13 },
            {  0.63,  6.54,  2.53 },
            {  0.08,  5.29,  3.18 },
            {  3.12,  6.24,  3.14 },
            {  3.57,  5.82,  2.40 },
            {  2.23,  5.90,  3.12 },
            {  6.25,  6.19,  3.06 },
            {  5.55,  5.59,  3.32 },
            {  6.08,  6.99,  3.55 },
            { -6.20, -6.16,  6.15 },
            { -6.29, -5.99,  7.09 },
            { -6.09, -7.11,  6.09 },
            { -3.09, -6.19,  6.27 },
            { -2.56, -5.90,  5.52 },
            { -3.80, -6.69,  5.87 },
            {  0.02, -6.25,  6.24 },
            { -0.70, -5.70,  6.51 },
            {  0.25, -5.93,  5.36 },
            {  3.11, -6.18,  6.14 },
            {  3.76, -6.54,  6.74 },
            {  2.29, -6.20,  6.64 },
            {  6.22, -6.17,  6.15 },
            {  6.61, -6.98,  6.47 },
            {  5.56, -5.94,  6.81 },
            { -6.21, -3.10,  6.14 },
            { -6.76, -2.66,  6.78 },
            { -5.51, -3.50,  6.65 },
            { -3.13, -3.05,  6.18 },
            { -2.19, -3.14,  6.34 },
            { -3.50, -3.89,  6.43 },
            {  0.01, -3.06,  6.15 },
            { -0.06, -2.81,  7.07 },
            { -0.25, -3.98,  6.13 },
            {  3.04, -3.09,  6.17 },
            {  3.84, -3.51,  5.84 },
            {  3.25, -2.85,  7.08 },
            {  6.26, -3.13,  6.19 },
            {  6.01, -2.20,  6.09 },
            {  5.47, -3.55,  6.54 },
            { -6.20,  0.01,  6.27 },
            { -5.79, -0.70,  5.78 },
            { -6.67,  0.51,  5.60 },
            { -3.13,  0.01,  6.14 },
            { -3.53, -0.35,  6.94 },
            { -2.21,  0.17,  6.39 },
            { -0.04, -0.04,  6.20 },
            {  0.26,  0.47,  5.46 },
            {  0.51,  0.22,  6.93 },
            {  3.10, -0.05,  6.23 },
            {  2.33,  0.44,  5.95 },
            {  3.85,  0.45,  5.92 },
            {  6.19, -0.01,  6.26 },
            {  7.05,  0.16,  5.88 },
            {  5.58,  0.02,  5.52 },
            { -6.22,  3.04,  6.17 },
            { -5.45,  3.57,  5.95 },
            { -6.62,  3.50,  6.92 },
            { -3.09,  3.16,  6.21 },
            { -3.71,  2.75,  5.61 },
            { -2.60,  2.43,  6.59 },
            { -0.02,  3.10,  6.26 },
            {  0.89,  3.27,  6.05 },
            { -0.44,  2.94,  5.41 },
            {  3.12,  3.04,  6.23 },
            {  2.31,  3.53,  6.43 },
            {  3.59,  3.60,  5.60 },
            {  6.23,  3.05,  6.24 },
            {  5.92,  3.91,  6.54 },
            {  6.02,  3.03,  5.30 },
            { -6.15,  6.21,  6.24 },
            { -6.27,  6.46,  5.32 },
            { -7.00,  5.85,  6.51 },
            { -3.07,  6.15,  6.22 },
            { -3.98,  6.27,  5.94 },
            { -2.66,  7.01,  6.10 },
            {  0.04,  6.20,  6.25 },
            { -0.38,  5.50,  5.75 },
            { -0.36,  7.00,  5.93 },
            {  3.12,  6.15,  6.24 },
            {  3.66,  6.88,  5.93 },
            {  2.25,  6.33,  5.86 },
            {  6.20,  6.27,  6.19 },
            {  5.46,  5.65,  6.19 },
            {  6.97,  5.73,  6.39 }
        };
        for (int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }
    else
        throw exception();

    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLength, 0, 0),
                                        Vec3(0, boxEdgeLength, 0),
                                        Vec3(0, 0, boxEdgeLength));

    sig.clear();
    eps.clear();
    bonds.clear();
    for (int atom = 0; atom < natoms; ++atom) {
        system.addParticle(masses[atom%RESSIZE]);
        double sigma = 2.0*pow(2.0, -1.0/6.0)*halfrmins[atom%RESSIZE]*OpenMM::NmPerAngstrom;
        double epsilon = fabs(epsilons[atom%RESSIZE])*OpenMM::KJPerKcal;
        sig.push_back(0.5*sigma);
        eps.push_back(2.0*sqrt(epsilon));
        if (atom%RESSIZE == 0) {
            bonds.push_back(pair<int, int>(atom, atom+1));
            bonds.push_back(pair<int, int>(atom, atom+2));
        }
        double charge = do_electrostatics ? charges[atom] : 0;
        forceField->addParticle(charge, sigma, epsilon);
    }
}

void testWater2DpmeEnergiesForcesNoExclusions() {
    const double cutoff = 7.0*OpenMM::NmPerAngstrom;
    const double alpha = 4.0124063605;
    const double dalpha = 4.0124063605;
    const int grid = 32;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 6;
    double boxEdgeLength = 25*OpenMM::NmPerAngstrom;

    make_waterbox(numAtoms, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);
    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);


//Gromacs reference values from the following inputs
//
//------------------------------------------------
// water2.gro
//------------------------------------------------
//MD of 2 waters, t= 0.0
//    6
//    1SOL    OW1    1   0.200   0.200   0.200
//    1SOL    HW1    2   0.250   0.200   0.300
//    1SOL    HW2    3   0.150   0.200   0.300
//    2SOL    OW1    4   0.000   0.000   0.000
//    2SOL    HW1    5   0.050   0.000   0.100
//    2SOL    HW2    6  -0.050   0.000   0.100
//    2.50000   2.50000   2.50000
//------------------------------------------------
//
//------------------------------------------------
// water2.gro
//------------------------------------------------
//[ defaults ]
//; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
//  1             2               no              1.0     1.0
//
//[ atomtypes ]
//; full atom descriptions are available in ffoplsaa.atp
//; name  bond_type    mass    charge   ptype          sigma      epsilon
//OW   OW 8     15.99940    -0.834       A    3.15057e-01  6.36386e-01
//HW   HW 1      1.00800     0.417       A    4.00014e-02  1.92464e-01
//
//
//#ifdef NOEXCLUDE
//
//[ moleculetype ]
//; molname   nrexcl
//SOL     0
//
//#else
//
//[ moleculetype ]
//; molname   nrexcl
//SOL     3
//
//#endif
//
//#ifdef NOCHARGE
//[ atoms ]
//; id    at type res nr  residu name at name     cg nr   charge
//1      OW  1       SOL             OW1             1        0.000
//2      HW  1       SOL             HW1             1        0.000
//3      HW  1       SOL             HW2             1        0.000
//#else
//[ atoms ]
//; id    at type res nr  residu name at name     cg nr   charge
//1      OW  1       SOL             OW1             1       -0.834
//2      HW  1       SOL             HW1             1        0.417
//3      HW  1       SOL             HW2             1        0.417
//#endif
//
//[ settles ]
//; i j   funct   length
//1   1   0.09572 0.15139
//
//#ifndef NOEXCLUDE
//[ exclusions ]
//1   2   3
//2   1   3
//3   1   2
//#endif
//------------------------------------------------
//
//------------------------------------------------
// water2.top
//------------------------------------------------
//#include "tip3p.tpi"
//
//[ system ]
//Water dimer
//
//[ molecules ]
//SOL 2
//------------------------------------------------
//
//------------------------------------------------
// water2.mdp
//------------------------------------------------
//define = -DNOCHARGE -DNOEXCLUDE
//;define = -DNOCHARGE
//rvdw             = 0.7
//vdw-modifier     = Potential-Shift
//;vdw-modifier     = None
//rlist            = 0.9
//rcoulomb         = 0.7
//nstfout          = 1
//fourierspacing   = 0.08
//pme-order        = 5
//ewald-rtol       = 1.5E-2
//ewald-rtol-lj    = 1.5E-2
//vdwtype          = PME
//lj-pme-comb-rule = geometric
//continuation     = yes
//------------------------------------------------
//
//run commands:-
//    gmx grompp -c water2.gro -p water2.top -f water2.mdp -o run.tpr
//    gmx mdrun -s run.tpr -o traj.trr -pforce 1E-16
//    gmx dump -f traj.trr

    double refenergy = 1.34804e+03;
    vector<Vec3> refforces(6);
    refforces[0] = Vec3( 3.15301e-01, -3.91114e-03, -6.68965e+04);
    refforces[1] = Vec3( 1.67241e+04, -3.79846e-02,  3.34487e+04);
    refforces[2] = Vec3(-1.67243e+04, -9.55931e-02,  3.34486e+04);
    refforces[3] = Vec3(-1.71643e+00, -1.76696e+00, -6.68993e+04);
    refforces[4] = Vec3( 1.67254e+04,  1.59080e+00,  3.34495e+04);
    refforces[5] = Vec3(-1.67239e+04,  3.13897e-01,  3.34491e+04);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_TOL(refenergy, energy, 1E-4);
    for (int n = 0; n < numAtoms; ++n)
        ASSERT_EQUAL_VEC(refforces[n], forces[n], 1E-4);
}


void testDMFDpmeEnergiesForcesWithExclusions() {
    const double cutoff = 7.0*OpenMM::NmPerAngstrom;
    const double alpha = 4.0124063605;
    const double dalpha = 4.0124063605;
    const int grid = 32;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 12;
    double boxEdgeLength = 20*OpenMM::NmPerAngstrom;

    make_dmfbox(numAtoms, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setUseSwitchingFunction(false);
    forceField->setUseDispersionCorrection(false);

    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    double refenergy = 1.29545e+01;
    vector<Vec3> refforces(12);
//Gromacs reference values from the following inputs
//------------------------------------------------
// dmf.mdp
//------------------------------------------------
//rvdw             = 0.7
//vdw-modifier     = Potential-Shift
//rlist            = 0.9
//rcoulomb         = 0.7
//nstfout          = 1
//fourierspacing   = 0.0625
//pme-order        = 5
//ewald-rtol       = 1.5E-2
//ewald-rtol-lj    = 1.5E-2
//vdwtype          = PME
//lj-pme-comb-rule = geometric
//continuation     = yes
//------------------------------------------------
// dmf.gro
//------------------------------------------------
//Guyana Rwanda Oman Macau Angola Cameroon Senegal
//   12
//  101DMF     HA    1   0.620   0.945   0.128
//  101DMF      C    2   0.560   0.905   0.045
//  101DMF      O    3   0.502   0.987  -0.026
//  101DMF      N    4   0.556   0.770   0.028
//  101DMF     CC    5   0.643   0.685   0.106
//  101DMF    HC1    6   0.680   0.741   0.193
//  101DMF    HC2    7   0.729   0.653   0.045
//  101DMF    HC3    8   0.587   0.597   0.140
//  101DMF     CT    9   0.477   0.708  -0.072
//  101DMF    HT1   10   0.377   0.688  -0.032
//  101DMF    HT2   11   0.524   0.614  -0.102
//  101DMF    HT3   12   0.468   0.775  -0.158
//   2.00000   2.00000   2.00000
//------------------------------------------------
// dmf.top
//------------------------------------------------
//[ defaults ]
//; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
//1	2	yes	1.0	1.0
//
//[ atomtypes ]
//HGR52     1     1.008000    0.000  A  0.160361769265  0.19246
//CG2O1     6    12.011000    0.000  A  0.356359487256  0.46024
//OG2D1     8    15.999400    0.000  A  0.302905564168  0.50208
//NG2S0     7    14.007000    0.000  A  0.329632525712  0.83680
//CG331     6    12.011000    0.000  A  0.365268474438  0.32635
// HGA3     1     1.008000    0.000  A  0.238760856462  0.10042
//
//[ pairtypes ]
//OG2D1 CG331  1  0.293996576986  0.144938011577
//
//[ bondtypes ]
//;      i        j  func           b0           kb
//   CG2O1    HGR52     1   0.11000000    0.0000000000
//   CG2O1    OG2D1     1   0.12300000    0.0000000000
//   CG2O1    NG2S0     1   0.13500000    0.0000000000
//   CG331    NG2S0     1   0.14340000    0.0000000000
//   CG331     HGA3     1   0.11110000    0.0000000000
//
//[ angletypes ]
//;      i        j        k  func       theta0       ktheta          ub0          kub
//    HGA3    CG331     HGA3     5   108.400000 0.0000000000000   0.18020000  0.0000000
//   NG2S0    CG331     HGA3     5   105.000000 0.0000000000000   0.00000000  0.0000000
//   CG331    NG2S0    CG331     5   121.000000 0.0000000000000   0.00000000  0.0000000
//   CG2O1    NG2S0    CG331     5   119.500000 0.0000000000000   0.00000000  0.0000000
//   NG2S0    CG2O1    OG2D1     5   124.000000 0.0000000000000   0.00000000  0.0000000
//   NG2S0    CG2O1    HGR52     5   115.000000 0.0000000000000   0.00000000  0.0000000
//   OG2D1    CG2O1    HGR52     5   122.000000 0.0000000000000   0.00000000  0.0000000
//
//[ dihedraltypes ]
//;      i        j        k        l  func         phi0         kphi  mult
//   OG2D1    CG2O1    NG2S0    CG331     9   180.000000   0.00000000     2
//   HGR52    CG2O1    NG2S0    CG331     9   180.000000   0.00000000     2
//    HGA3    CG331    NG2S0    CG331     9     0.000000   0.00000000     3
//    HGA3    CG331    NG2S0    CG2O1     9     0.000000   0.00000000     3
//
//[ dihedraltypes ]
//; 'improper' dihedrals 
//;      i        j        k        l  func         phi0         kphi
//   CG2O1    NG2S0    OG2D1    HGR52     2     0.000000 0.0000000
//
//[ moleculetype ]
//; Name            nrexcl
//DMF       3
//
//[ atoms ]
//;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
//; residue 101 DMF rtp DMF  q  0.0
//     1      HGR52    101    DMF     HA      1    0.00      1.008   ; qtot 0.00
//     2      CG2O1    101    DMF      C      2    0.00     12.011   ; qtot 0.00
//     3      OG2D1    101    DMF      O      3    0.00    15.9994   ; qtot 0.00
//     4      NG2S0    101    DMF      N      4    0.00     14.007   ; qtot 0.00
//     5      CG331    101    DMF     CC      5    0.00     12.011   ; qtot 0.00
//     6       HGA3    101    DMF    HC1      6    0.00      1.008   ; qtot 0.00
//     7       HGA3    101    DMF    HC2      7    0.00      1.008   ; qtot 0.00
//     8       HGA3    101    DMF    HC3      8    0.00      1.008   ; qtot 0.00
//     9      CG331    101    DMF     CT      9    0.00     12.011   ; qtot 0.00
//    10       HGA3    101    DMF    HT1     10    0.00      1.008   ; qtot 0.00
//    11       HGA3    101    DMF    HT2     11    0.00      1.008   ; qtot 0.00
//    12       HGA3    101    DMF    HT3     12    0.00      1.008   ; qtot 0
//
//[ bonds ]
//;  ai    aj funct            c0            c1            c2            c3
//    1     2     1 
//    2     3     1 
//    2     4     1 
//    4     5     1 
//    4     9     1 
//    5     6     1 
//    5     7     1 
//    5     8     1 
//    9    10     1 
//    9    11     1 
//    9    12     1 
//
//[ pairs ]
//;  ai    aj funct            c0            c1            c2            c3
//    1     5     1 
//    1     9     1 
//    2     6     1 
//    2     7     1 
//    2     8     1 
//    2    10     1 
//    2    11     1 
//    2    12     1 
//    3     5     1 
//    3     9     1 
//    5    10     1 
//    5    11     1 
//    5    12     1 
//    6     9     1 
//    7     9     1 
//    8     9     1 
//
//[ angles ]
//;  ai    aj    ak funct            c0            c1            c2            c3
//    1     2     3     5 
//    1     2     4     5 
//    3     2     4     5 
//    2     4     5     5 
//    2     4     9     5 
//    5     4     9     5 
//    4     5     6     5 
//    4     5     7     5 
//    4     5     8     5 
//    6     5     7     5 
//    6     5     8     5 
//    7     5     8     5 
//    4     9    10     5 
//    4     9    11     5 
//    4     9    12     5 
//   10     9    11     5 
//   10     9    12     5 
//   11     9    12     5 
//
//[ dihedrals ]
//;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
//    1     2     4     5     9 
//    1     2     4     9     9 
//    3     2     4     5     9 
//    3     2     4     9     9 
//    2     4     5     6     9 
//    2     4     5     7     9 
//    2     4     5     8     9 
//    9     4     5     6     9 
//    9     4     5     7     9 
//    9     4     5     8     9 
//    2     4     9    10     9 
//    2     4     9    11     9 
//    2     4     9    12     9 
//    5     4     9    10     9 
//    5     4     9    11     9 
//    5     4     9    12     9 
//
//[ dihedrals ]
//;  ai    aj    ak    al funct            c0            c1            c2            c3
//    1     3     4     2     2 
//
//[ system ]
//; Name
//DMF Monomer
//
//[ molecules ]
//; Compound        #mols
//DMF       1
//
//
//run commands:-
//gmx_d grompp -c dmf.gro -p dmf.top -f dmf.mdp -o run.tpr
//gmx_d mdrun -s run.tpr -o traj.trr -pforce 1E-16 -debug
//gmx_d dump -f traj.trr

    // Gromacs reference values.  See inputs above.
    refforces[ 0] = Vec3(-3.19942e+00,  2.19766e+01,  2.72861e-01);
    refforces[ 1] = Vec3(-5.22879e+01,  2.82006e+02, -6.86534e+00);
    refforces[ 2] = Vec3( 1.23583e+01,  7.31235e+01,  4.08375e+01);
    refforces[ 3] = Vec3(-1.28656e-03,  1.14097e-03,  4.07527e-04);
    refforces[ 4] = Vec3( 1.51942e+02,  5.65684e+01,  2.41818e+02);
    refforces[ 5] = Vec3( 1.20290e+02, -1.65176e+02,  1.48309e+02);
    refforces[ 6] = Vec3( 4.57022e+01, -1.65198e+01,  1.88824e+01);
    refforces[ 7] = Vec3( 5.57228e+01, -5.68770e+01,  1.09979e+02);
    refforces[ 8] = Vec3(-9.58422e+01,  4.41982e+01, -1.27887e+02);
    refforces[ 9] = Vec3(-2.59754e+01, -1.42402e+01, -1.23475e+01);
    refforces[10] = Vec3(-1.37259e+02, -7.97690e+01, -2.39814e+02);
    refforces[11] = Vec3(-7.14520e+01, -1.45287e+02, -1.73189e+02);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();

    ASSERT_EQUAL_TOL(refenergy, energy, 1E-4);
    for (int n = 0; n < numAtoms; ++n)
        // The forces on atom 3 are much smaller in magnitude to the others,
        // which causes problems for testing in single precision
        if(n != 3)
            ASSERT_EQUAL_VEC(refforces[n], forces[n], 5E-4);
}


void testWater2DpmeEnergiesForcesWithExclusions() {
    const double cutoff = 7.0*OpenMM::NmPerAngstrom;
    const double alpha = 4.0124063605;
    const double dalpha = 4.0124063605;
    const int grid = 32;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 6;
    double boxEdgeLength = 25*OpenMM::NmPerAngstrom;

    make_waterbox(numAtoms, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->createExceptionsFromBonds(bonds, 1.0, 1.0);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setUseSwitchingFunction(false);
    forceField->setUseDispersionCorrection(false);

    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    double refenergy = -7.78060e-01;
    vector<Vec3> refforces(6);

    // Gromacs reference values.  See comments in testWater2DpmeEnergiesForcesNoExclusions() for details.
    refforces[0] = Vec3( 3.15301e-01, -3.91114e-03,  9.48373e-01);
    refforces[1] = Vec3(-4.74745e-02, -3.79846e-02, -5.69616e-02);
    refforces[2] = Vec3(-7.16866e-02, -9.55931e-02, -1.43348e-01);
    refforces[3] = Vec3(-1.78136e+00, -1.76696e+00, -1.70022e+00);
    refforces[4] = Vec3( 1.19309e+00,  1.59080e+00,  7.95443e-01);
    refforces[5] = Vec3( 3.92366e-01,  3.13897e-01,  1.56962e-01);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();


    ASSERT_EQUAL_TOL(refenergy, energy, 5E-4);
    for (int n = 0; n < numAtoms; ++n)
        ASSERT_EQUAL_VEC(refforces[n], forces[n], 5E-4);
}


void testWater125DpmeVsLongCutoffNoExclusions() {
    const double cutoff = 8.5*OpenMM::NmPerAngstrom;
    const double alpha = 0.45*OpenMM::AngstromsPerNm;
    const double dalpha = 0.45*OpenMM::AngstromsPerNm;
    const int grid = 32;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 375;
    double boxEdgeLength = 17.01*OpenMM::NmPerAngstrom;

    make_waterbox(numAtoms, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();

    // Make another call to the get the energy to test that the call is
    // idempotent when Coulomb terms are absent.
    double secondEnergy = context.getState(State::Energy).getPotentialEnergy();
    ASSERT_EQUAL_TOL(secondEnergy, energy, 5E-5);

//Gromacs reference values.  See comments in testWater2DpmeEnergiesForcesNoExclusions() for details.
//Coordinates are from make_waterbox, and the .gro file looks like
//
//;define = -DNOCHARGE -DNOEXCLUDE
//define = -DNOCHARGE
//vdw-modifier     = Potential-Shift
//rvdw             = 1.15
//rlist            = 1.15
//rcoulomb         = 1.15
//fourierspacing   = 0.04
//nstfout          = 1
//pme-order        = 5
//ewald-rtol       = 1E-8
//ewald-rtol-lj    = 1E-8
//vdwtype          = PME
//lj-pme-comb-rule = geometric
//continuation     = yes

    double gromacs_energy = 5.63157e+05;
    const vector<Vec3>& forces = state.getForces();
    vector<Vec3> refforces(numAtoms);
    refforces[  0] = Vec3(-1.12446e+05, -2.71952e+05, -1.91471e+05);
    refforces[  1] = Vec3( 2.70479e+05,  6.61172e+04,  7.21277e+04);
    refforces[  2] = Vec3(-1.58058e+05,  2.05779e+05,  1.19292e+05);
    refforces[  3] = Vec3( 3.16438e+05, -1.27994e+05,  1.08811e+05);
    refforces[  4] = Vec3(-1.36520e+05,  1.93405e+05,  1.36521e+05);
    refforces[  5] = Vec3(-1.79957e+05, -6.54374e+04, -2.45393e+05);
    refforces[  6] = Vec3( 1.30626e+05, -1.36728e+05,  2.93833e+05);
    refforces[  7] = Vec3(-2.74564e+05,  8.02095e+04, -7.09524e+04);
    refforces[  8] = Vec3( 1.43952e+05,  5.64523e+04, -2.22980e+05);
    refforces[  9] = Vec3(-4.22852e+04,  2.14660e+04, -3.23254e+05);
    refforces[ 10] = Vec3( 2.28059e+05, -4.44251e+04,  1.62898e+05);
    refforces[ 11] = Vec3(-1.85776e+05,  2.29044e+04,  1.60326e+05);
    refforces[ 12] = Vec3(-1.02075e+05,  3.27346e+05,  1.47993e+04);
    refforces[ 13] = Vec3( 7.77195e+04, -1.44336e+05,  2.10959e+05);
    refforces[ 14] = Vec3( 2.44142e+04, -1.83111e+05, -2.25837e+05);
    refforces[ 15] = Vec3(-4.81843e+04, -2.72889e+05, -1.75068e+05);
    refforces[ 16] = Vec3(-5.46666e+03,  2.18711e+04,  2.62456e+05);
    refforces[ 17] = Vec3( 5.35890e+04,  2.51021e+05, -8.74311e+04);
    refforces[ 18] = Vec3(-2.04721e+05,  1.58918e+05,  2.20448e+05);
    refforces[ 19] = Vec3(-7.05932e+04, -1.64717e+05, -2.17659e+05);
    refforces[ 20] = Vec3( 2.75339e+05,  5.73611e+03, -2.86718e+03);
    refforces[ 21] = Vec3(-5.65480e+03, -2.81809e+05, -1.38358e+05);
    refforces[ 22] = Vec3(-7.85576e+03,  2.25204e+05, -1.15217e+05);
    refforces[ 23] = Vec3( 1.34846e+04,  5.66345e+04,  2.53512e+05);
    refforces[ 24] = Vec3(-7.73167e+04, -4.43114e+04,  3.22354e+05);
    refforces[ 25] = Vec3(-1.24372e+05,  1.61973e+05, -1.88001e+05);
    refforces[ 26] = Vec3( 2.01689e+05, -1.17651e+05, -1.34455e+05);
    refforces[ 27] = Vec3(-1.79300e+05, -1.97922e+05,  1.94294e+05);
    refforces[ 28] = Vec3( 2.38945e+05, -4.88761e+04, -9.50350e+04);
    refforces[ 29] = Vec3(-5.95911e+04,  2.46878e+05, -9.93159e+04);
    refforces[ 30] = Vec3(-1.31892e+04, -2.15675e+05,  2.68757e+05);
    refforces[ 31] = Vec3( 2.31232e+05,  1.08107e+05, -1.32129e+05);
    refforces[ 32] = Vec3(-2.18089e+05,  1.07592e+05, -1.36669e+05);
    refforces[ 33] = Vec3( 1.90632e+04, -3.62889e+05,  8.61608e+04);
    refforces[ 34] = Vec3(-2.16178e+05,  1.59639e+05, -1.66288e+05);
    refforces[ 35] = Vec3( 1.97134e+05,  2.03295e+05,  8.00862e+04);
    refforces[ 36] = Vec3( 3.40071e+05, -6.87734e+04,  1.22584e+05);
    refforces[ 37] = Vec3(-4.17943e+04,  8.35901e+03, -2.64693e+05);
    refforces[ 38] = Vec3(-2.98358e+05,  6.03813e+04,  1.42075e+05);
    refforces[ 39] = Vec3(-3.21693e+05,  4.40885e+04,  1.41154e+04);
    refforces[ 40] = Vec3( 1.28817e+05,  2.02067e+04, -2.07114e+05);
    refforces[ 41] = Vec3( 1.92948e+05, -6.43160e+04,  1.92949e+05);
    refforces[ 42] = Vec3(-1.46001e+05,  3.20840e+05,  7.97305e+04);
    refforces[ 43] = Vec3(-1.27705e+05, -2.55411e+05, -1.24428e+05);
    refforces[ 44] = Vec3( 2.73737e+05, -6.54594e+04,  4.46323e+04);
    refforces[ 45] = Vec3( 1.50212e+04,  2.43629e+05, -2.20084e+05);
    refforces[ 46] = Vec3(-1.07432e+05,  8.26408e+03,  2.42419e+05);
    refforces[ 47] = Vec3( 9.23685e+04, -2.51915e+05, -2.23909e+04);
    refforces[ 48] = Vec3( 3.14327e+04, -2.94199e+05,  1.55484e+05);
    refforces[ 49] = Vec3(-2.23665e+05,  1.52883e+05, -2.54793e+04);
    refforces[ 50] = Vec3( 1.92227e+05,  1.41342e+05, -1.30033e+05);
    refforces[ 51] = Vec3( 5.73542e+04, -2.08689e+05, -2.68170e+05);
    refforces[ 52] = Vec3(-2.26784e+05,  1.85259e+05,  8.30474e+04);
    refforces[ 53] = Vec3( 1.69446e+05,  2.34613e+04,  1.85087e+05);
    refforces[ 54] = Vec3( 2.25490e+05, -1.91311e+05, -1.40103e+05);
    refforces[ 55] = Vec3(-8.20802e+03,  6.83988e+04,  2.54448e+05);
    refforces[ 56] = Vec3(-2.17315e+05,  1.22953e+05, -1.14373e+05);
    refforces[ 57] = Vec3(-7.09510e+04,  2.05639e+05, -2.69540e+05);
    refforces[ 58] = Vec3( 1.93603e+05,  3.38041e+04,  2.18192e+05);
    refforces[ 59] = Vec3(-1.22579e+05, -2.39458e+05,  5.13124e+04);
    refforces[ 60] = Vec3(-1.07250e+05,  3.35982e+05,  6.89600e+03);
    refforces[ 61] = Vec3( 6.90672e-01, -1.44228e+05, -2.24659e+05);
    refforces[ 62] = Vec3( 1.07222e+05, -1.91701e+05,  2.17695e+05);
    refforces[ 63] = Vec3( 2.64846e+05,  1.94002e+05,  5.27262e+03);
    refforces[ 64] = Vec3(-1.12985e+04, -2.71176e+05,  8.47510e+03);
    refforces[ 65] = Vec3(-2.53630e+05,  7.71894e+04, -1.37831e+04);
    refforces[ 66] = Vec3(-3.04551e+05,  4.21918e+04,  1.88948e+05);
    refforces[ 67] = Vec3( 2.87326e+05,  1.22193e+05,  3.30263e+04);
    refforces[ 68] = Vec3( 1.73006e+04, -1.64358e+05, -2.22023e+05);
    refforces[ 69] = Vec3( 2.45556e+04,  2.20596e+05,  2.41913e+05);
    refforces[ 70] = Vec3( 1.50804e+05, -2.17829e+05, -4.46813e+04);
    refforces[ 71] = Vec3(-1.75369e+05, -2.74087e+03, -1.97287e+05);
    refforces[ 72] = Vec3( 8.65787e+04, -2.77188e+04, -3.15014e+05);
    refforces[ 73] = Vec3(-2.42127e+05,  6.26659e+04,  1.11092e+05);
    refforces[ 74] = Vec3( 1.55591e+05, -3.48752e+04,  2.03883e+05);
    refforces[ 75] = Vec3( 7.06224e+04,  2.96642e+05, -1.51267e+05);
    refforces[ 76] = Vec3(-5.39280e+04, -2.57657e+05, -1.13850e+05);
    refforces[ 77] = Vec3(-1.67422e+04, -3.90659e+04,  2.65102e+05);
    refforces[ 78] = Vec3(-4.52577e+04, -3.21552e+05, -7.01075e+04);
    refforces[ 79] = Vec3( 2.35182e+05,  1.45173e+05,  3.48411e+04);
    refforces[ 80] = Vec3(-1.89931e+05,  1.76363e+05,  3.52722e+04);
    refforces[ 81] = Vec3(-2.29347e+05,  1.06977e+05, -2.70703e+05);
    refforces[ 82] = Vec3(-1.17929e+04, -2.56493e+05,  1.17929e+05);
    refforces[ 83] = Vec3( 2.41161e+05,  1.49452e+05,  1.52847e+05);
    refforces[ 84] = Vec3( 2.32633e+03,  3.03443e+05,  1.27473e+05);
    refforces[ 85] = Vec3(-2.11215e+05, -1.63335e+05, -4.50583e+04);
    refforces[ 86] = Vec3( 2.08887e+05, -1.40170e+05, -8.24544e+04);
    refforces[ 87] = Vec3( 5.83172e+04,  2.82133e+04, -3.26001e+05);
    refforces[ 88] = Vec3( 1.76890e+05, -4.71698e+04,  2.15220e+05);
    refforces[ 89] = Vec3(-2.35168e+05,  1.89225e+04,  1.10825e+05);
    refforces[ 90] = Vec3(-2.72444e+05, -1.46998e+05, -1.00637e+05);
    refforces[ 91] = Vec3( 2.78426e+04,  2.39442e+05,  1.16936e+05);
    refforces[ 92] = Vec3( 2.44570e+05, -9.23920e+04, -1.63042e+04);
    refforces[ 93] = Vec3(-3.10100e+04,  2.93390e+05, -1.87196e+05);
    refforces[ 94] = Vec3(-6.38830e+04, -2.90671e+05, -6.38833e+04);
    refforces[ 95] = Vec3( 9.48775e+04, -2.79002e+03,  2.51149e+05);
    refforces[ 96] = Vec3( 4.18753e+04, -1.23438e+05, -3.10188e+05);
    refforces[ 97] = Vec3( 1.29157e+05,  2.04500e+05,  9.41770e+04);
    refforces[ 98] = Vec3(-1.71038e+05, -8.10180e+04,  2.16049e+05);
    refforces[ 99] = Vec3(-5.61490e+04,  2.26993e+04, -3.44090e+05);
    refforces[100] = Vec3(-1.96230e+05, -2.88572e+04,  1.93344e+05);
    refforces[101] = Vec3( 2.52383e+05,  6.15570e+03,  1.50813e+05);
    refforces[102] = Vec3(-6.33126e+04,  3.55941e+05,  8.51300e+04);
    refforces[103] = Vec3(-1.75405e+05, -1.81783e+05, -1.69027e+05);
    refforces[104] = Vec3( 2.38759e+05, -1.74232e+05,  8.38891e+04);
    refforces[105] = Vec3( 1.51480e+05, -4.90606e+04,  3.17956e+05);
    refforces[106] = Vec3( 4.93229e+04, -1.61668e+05, -2.02771e+05);
    refforces[107] = Vec3(-2.00831e+05,  2.10712e+05, -1.15233e+05);
    refforces[108] = Vec3( 2.20204e+05, -2.33820e+05,  1.51386e+05);
    refforces[109] = Vec3( 3.36498e+04,  2.79296e+05, -1.51424e+05);
    refforces[110] = Vec3(-2.53896e+05, -4.54334e+04,  2.10242e-01);
    refforces[111] = Vec3(-1.94665e+05,  2.05890e+05,  2.40510e+05);
    refforces[112] = Vec3(-9.73233e+04, -1.29764e+05, -2.62775e+05);
    refforces[113] = Vec3( 2.92049e+05, -7.61856e+04,  2.22208e+04);
    refforces[114] = Vec3( 1.60261e+05, -3.43716e+05,  1.52456e+04);
    refforces[115] = Vec3( 1.11151e+05,  3.08358e+05, -8.60527e+04);
    refforces[116] = Vec3(-2.71443e+05,  3.54056e+04,  7.08103e+04);
    refforces[117] = Vec3(-4.27331e+04, -3.19871e+05, -1.68415e+05);
    refforces[118] = Vec3( 2.28408e+05,  2.15171e+05, -2.31720e+04);
    refforces[119] = Vec3(-1.85616e+05,  1.04782e+05,  1.91602e+05);
    refforces[120] = Vec3(-1.66057e+05, -9.58191e+04, -2.78448e+05);
    refforces[121] = Vec3( 1.91258e+05,  2.19476e+05,  6.89779e+04);
    refforces[122] = Vec3(-2.52379e+04, -1.23674e+05,  2.09489e+05);
    refforces[123] = Vec3( 6.55172e+03, -9.23173e+04,  3.29554e+05);
    refforces[124] = Vec3(-2.14685e+05,  1.13143e+05, -1.36353e+05);
    refforces[125] = Vec3( 2.08124e+05, -2.08124e+04, -1.93257e+05);
    refforces[126] = Vec3( 1.02696e+05, -3.39299e+05, -4.47114e+04);
    refforces[127] = Vec3(-1.81786e+05,  1.75407e+05, -1.69029e+05);
    refforces[128] = Vec3( 7.90534e+04,  1.63962e+05,  2.13738e+05);
    refforces[129] = Vec3(-3.45588e+05,  9.36838e+04,  5.67939e+04);
    refforces[130] = Vec3( 1.44967e+05, -2.60943e+05,  7.08730e+04);
    refforces[131] = Vec3( 2.00649e+05,  1.67207e+05, -1.27685e+05);
    refforces[132] = Vec3(-2.70181e+05,  2.05707e+05, -3.11852e+04);
    refforces[133] = Vec3( 1.09331e+05, -1.83207e+05, -1.86162e+05);
    refforces[134] = Vec3( 1.60883e+05, -2.25804e+04,  2.17340e+05);
    refforces[135] = Vec3(-1.04496e+05, -2.97000e+05, -1.63733e+05);
    refforces[136] = Vec3( 2.11303e+05,  1.73661e+04,  1.79462e+05);
    refforces[137] = Vec3(-1.06849e+05,  2.79694e+05, -1.57134e+04);
    refforces[138] = Vec3(-3.82271e+04,  2.11942e+05,  2.60118e+05);
    refforces[139] = Vec3(-1.96098e+05, -1.23692e+05, -1.71962e+05);
    refforces[140] = Vec3( 2.34334e+05, -8.82195e+04, -8.82188e+04);
    refforces[141] = Vec3( 2.17634e+05,  2.72527e+05,  1.42967e+05);
    refforces[142] = Vec3( 9.30924e+04, -1.86185e+05, -1.98197e+05);
    refforces[143] = Vec3(-3.10770e+05, -8.63246e+04,  5.52480e+04);
    refforces[144] = Vec3(-1.63880e+05, -2.98869e+05,  1.01299e+05);
    refforces[145] = Vec3( 6.83424e+04,  2.39196e+05,  1.61538e+05);
    refforces[146] = Vec3( 9.55787e+04,  5.97354e+04, -2.62846e+05);
    refforces[147] = Vec3( 1.06430e+05, -2.97083e+05, -7.97261e+04);
    refforces[148] = Vec3(-1.14920e+05,  6.41394e+04,  2.21823e+05);
    refforces[149] = Vec3( 8.52531e+03,  2.33043e+05, -1.42101e+05);
    refforces[150] = Vec3(-4.96390e+04, -4.12076e+04, -3.67834e+05);
    refforces[151] = Vec3( 1.25352e+05, -1.99962e+05,  1.61166e+05);
    refforces[152] = Vec3(-7.57841e+04,  2.41138e+05,  2.06689e+05);
    refforces[153] = Vec3(-3.06397e+05, -5.15222e+04, -1.37080e+05);
    refforces[154] = Vec3( 1.92949e+05, -1.92945e+05,  6.43163e+04);
    refforces[155] = Vec3( 1.13490e+05,  2.44443e+05,  7.27504e+04);
    refforces[156] = Vec3(-3.74471e+03,  3.83223e+05, -2.14754e+04);
    refforces[157] = Vec3( 2.17881e+05, -1.85835e+05, -1.05735e+05);
    refforces[158] = Vec3(-2.14191e+05, -1.97454e+05,  1.27176e+05);
    refforces[159] = Vec3(-3.33357e+05, -1.38307e+04, -1.17320e+05);
    refforces[160] = Vec3( 2.04157e+05, -9.93173e+04, -1.37945e+05);
    refforces[161] = Vec3( 1.29262e+05,  1.13105e+05,  2.55294e+05);
    refforces[162] = Vec3( 2.50185e+05,  2.64217e+05, -7.18219e+04);
    refforces[163] = Vec3(-2.46668e+05,  1.94013e+04, -9.97745e+04);
    refforces[164] = Vec3(-3.50266e+03, -2.83658e+05,  1.71598e+05);
    refforces[165] = Vec3(-2.05824e+05,  2.71177e+05, -1.16442e+05);
    refforces[166] = Vec3(-3.52163e+04, -1.88897e+05,  2.36923e+05);
    refforces[167] = Vec3( 2.41012e+05, -8.22954e+04, -1.20505e+05);
    refforces[168] = Vec3( 6.89329e+04,  2.09489e+05,  2.76583e+05);
    refforces[169] = Vec3( 1.87999e+05, -1.61968e+05, -1.24368e+05);
    refforces[170] = Vec3(-2.56932e+05, -4.75793e+04, -1.52255e+05);
    refforces[171] = Vec3( 3.39431e+05, -5.75509e+04,  1.62795e+05);
    refforces[172] = Vec3(-1.27767e+05,  2.66184e+05, -1.59709e+05);
    refforces[173] = Vec3(-2.11726e+05, -2.08613e+05, -3.11349e+03);
    refforces[174] = Vec3(-2.58601e+05,  4.61975e+04, -2.46016e+05);
    refforces[175] = Vec3( 7.15587e+04, -2.52015e+05,  1.40007e+05);
    refforces[176] = Vec3( 1.87109e+05,  2.05820e+05,  1.06028e+05);
    refforces[177] = Vec3( 3.92176e+03,  2.94819e+05, -1.84065e+05);
    refforces[178] = Vec3(-1.67230e+05, -8.36141e+04,  2.29167e+05);
    refforces[179] = Vec3( 1.63334e+05, -2.11213e+05, -4.50586e+04);
    refforces[180] = Vec3( 1.11151e+05,  2.40551e+05, -2.68210e+05);
    refforces[181] = Vec3(-1.76495e+05, -2.47099e+05, -3.53002e+04);
    refforces[182] = Vec3( 6.52872e+04,  6.52855e+03,  3.03586e+05);
    refforces[183] = Vec3(-4.84049e+04, -2.73305e+05, -2.95224e+05);
    refforces[184] = Vec3(-9.04206e+04, -1.44672e+04,  3.29135e+05);
    refforces[185] = Vec3( 1.38830e+05,  2.87820e+05, -3.38610e+04);
    refforces[186] = Vec3(-2.09072e+05, -1.53610e+05, -2.86654e+05);
    refforces[187] = Vec3(-1.30322e+05,  9.09223e+04,  2.42461e+05);
    refforces[188] = Vec3( 3.39380e+05,  6.27113e+04,  4.42669e+04);
    refforces[189] = Vec3(-3.15692e+05,  1.48900e+05, -9.20420e+04);
    refforces[190] = Vec3( 7.13687e+04, -2.72503e+05,  1.26518e+05);
    refforces[191] = Vec3( 2.44351e+05,  1.23612e+05, -3.44964e+04);
    refforces[192] = Vec3(-2.80727e+04,  3.15077e+05, -2.05427e+05);
    refforces[193] = Vec3(-2.29004e+05, -2.08496e+05,  9.57028e+04);
    refforces[194] = Vec3( 2.57096e+05, -1.06603e+05,  1.09738e+05);
    refforces[195] = Vec3( 3.33213e+05, -7.80040e+04, -5.04749e+03);
    refforces[196] = Vec3(-2.07680e+05, -7.82574e+04,  1.83605e+05);
    refforces[197] = Vec3(-1.25574e+05,  1.56274e+05, -1.78599e+05);
    refforces[198] = Vec3( 2.66795e+05, -2.82879e+04, -2.26620e+05);
    refforces[199] = Vec3(-2.28291e+05, -1.82015e+05,  4.01046e+04);
    refforces[200] = Vec3(-3.85027e+04,  2.10288e+05,  1.86592e+05);
    refforces[201] = Vec3( 1.93076e+05,  2.05300e+05,  2.64591e+05);
    refforces[202] = Vec3(-3.32259e+05, -3.65117e+04, -8.39768e+04);
    refforces[203] = Vec3( 1.39204e+05, -1.68822e+05, -1.80669e+05);
    refforces[204] = Vec3( 2.15422e+05,  2.88264e+05,  2.49786e+04);
    refforces[205] = Vec3( 4.29235e+04, -2.66744e+05,  1.13441e+05);
    refforces[206] = Vec3(-2.58339e+05, -2.15280e+04, -1.38395e+05);
    refforces[207] = Vec3( 3.27580e+05, -4.93781e+04,  7.43625e+03);
    refforces[208] = Vec3(-2.11622e+05, -1.58716e+05, -9.69925e+04);
    refforces[209] = Vec3(-1.15917e+05,  2.08125e+05,  8.95715e+04);
    refforces[210] = Vec3( 1.10967e+05, -2.71544e+05, -2.06279e+05);
    refforces[211] = Vec3(-8.86155e+04,  1.96918e+04,  2.98676e+05);
    refforces[212] = Vec3(-2.23913e+04,  2.51912e+05, -9.23690e+04);
    refforces[213] = Vec3( 1.91611e+05, -7.99998e+04, -2.83463e+05);
    refforces[214] = Vec3( 7.08725e+04,  1.44963e+05,  2.60941e+05);
    refforces[215] = Vec3(-2.62503e+05, -6.49202e+04,  2.25808e+04);
    refforces[216] = Vec3(-6.63168e+04, -2.84265e+04,  3.72700e+05);
    refforces[217] = Vec3(-2.02131e+05, -6.33547e+04, -1.96097e+05);
    refforces[218] = Vec3( 2.68466e+05,  9.18414e+04, -1.76621e+05);
    refforces[219] = Vec3(-6.80052e+03,  2.72744e+05, -2.21849e+05);
    refforces[220] = Vec3( 1.52709e+05, -5.52356e+04,  2.63181e+05);
    refforces[221] = Vec3(-1.45890e+05, -2.17461e+05, -4.12893e+04);
    refforces[222] = Vec3( 3.07522e+05, -1.21146e+05,  1.14593e+05);
    refforces[223] = Vec3(-2.11787e+05, -1.56664e+05, -8.99363e+04);
    refforces[224] = Vec3(-9.57076e+04,  2.77857e+05, -2.46987e+04);
    refforces[225] = Vec3(-3.24882e+05, -3.09840e+04, -1.31947e+05);
    refforces[226] = Vec3( 8.33124e+04, -5.05809e+04,  2.67793e+05);
    refforces[227] = Vec3( 2.41536e+05,  8.15185e+04, -1.35863e+05);
    refforces[228] = Vec3(-2.16561e+04, -1.67607e+05, -2.70252e+05);
    refforces[229] = Vec3( 1.10825e+05, -7.29802e+04,  2.24354e+05);
    refforces[230] = Vec3(-8.91995e+04,  2.40572e+05,  4.59512e+04);
    refforces[231] = Vec3(-2.22246e+05,  1.96766e+05, -2.46638e+05);
    refforces[232] = Vec3( 3.04748e+05,  5.66971e+04,  1.27568e+05);
    refforces[233] = Vec3(-8.24632e+04, -2.53495e+05,  1.19113e+05);
    refforces[234] = Vec3( 2.20621e+05, -2.03896e+05,  6.63128e+04);
    refforces[235] = Vec3(-9.59082e+04,  1.64055e+05,  1.53960e+05);
    refforces[236] = Vec3(-1.24758e+05,  3.98170e+04, -2.20320e+05);
    refforces[237] = Vec3(-7.79583e+03, -3.49688e+04,  3.64330e+05);
    refforces[238] = Vec3(-1.43293e+05, -1.65580e+05, -2.10163e+05);
    refforces[239] = Vec3( 1.51159e+05,  2.00522e+05, -1.54247e+05);
    refforces[240] = Vec3( 1.82058e+05, -3.72909e+04, -2.80373e+05);
    refforces[241] = Vec3(-1.95792e+05,  1.98809e+05,  7.22936e+04);
    refforces[242] = Vec3( 1.36910e+04, -1.61546e+05,  2.08094e+05);
    refforces[243] = Vec3( 1.40288e+05,  3.27379e+05,  4.78569e+03);
    refforces[244] = Vec3(-1.70016e+05, -1.73349e+05,  2.03353e+05);
    refforces[245] = Vec3( 2.97331e+04, -1.54072e+05, -2.08134e+05);
    refforces[246] = Vec3( 1.21932e+05,  3.52267e+05,  4.69291e+04);
    refforces[247] = Vec3(-2.91621e+05, -1.24021e+05,  3.01674e+04);
    refforces[248] = Vec3( 1.69673e+05, -2.28288e+05, -7.71238e+04);
    refforces[249] = Vec3(-1.41105e+05,  1.52818e+05,  2.59197e+05);
    refforces[250] = Vec3( 2.15511e+05, -1.77650e+05, -5.82442e+03);
    refforces[251] = Vec3(-7.43779e+04,  2.47926e+04, -2.53438e+05);
    refforces[252] = Vec3(-3.34214e+04,  3.09554e+05,  1.58194e+05);
    refforces[253] = Vec3(-2.05877e+05, -1.71563e+05, -6.00467e+04);
    refforces[254] = Vec3( 2.39329e+05, -1.38076e+05, -9.81875e+04);
    refforces[255] = Vec3( 1.32793e+05, -3.72271e+05,  2.88491e+04);
    refforces[256] = Vec3(-9.02633e+04,  2.78646e+05,  2.23700e+05);
    refforces[257] = Vec3(-4.25632e+04,  9.36427e+04, -2.52553e+05);
    refforces[258] = Vec3( 2.97252e+05,  2.17294e+05, -2.50945e+03);
    refforces[259] = Vec3(-1.35724e+05, -1.48966e+05,  2.41658e+05);
    refforces[260] = Vec3(-1.61536e+05, -6.83423e+04, -2.39199e+05);
    refforces[261] = Vec3( 2.50550e+05, -1.39926e+05,  2.48375e+05);
    refforces[262] = Vec3( 5.51783e+04, -1.65530e+04, -2.59341e+05);
    refforces[263] = Vec3(-3.05734e+05,  1.56505e+05,  1.09189e+04);
    refforces[264] = Vec3(-4.44616e+04,  4.17040e+04, -3.31507e+05);
    refforces[265] = Vec3(-1.79702e+05, -8.54319e+04,  2.00323e+05);
    refforces[266] = Vec3( 2.24182e+05,  4.37421e+04,  1.31227e+05);
    refforces[267] = Vec3(-9.89367e+04, -3.76136e+05,  2.17157e+04);
    refforces[268] = Vec3(-4.44424e+04,  1.68244e+05, -2.47606e+05);
    refforces[269] = Vec3( 1.43421e+05,  2.07965e+05,  2.25893e+05);
    refforces[270] = Vec3(-1.09030e+03, -2.44686e+05, -2.30145e+05);
    refforces[271] = Vec3(-1.86963e+05,  1.89757e+05,  3.34863e+04);
    refforces[272] = Vec3( 1.88003e+05,  5.49541e+04,  1.96680e+05);
    refforces[273] = Vec3(-1.15452e+05, -1.55243e+05,  2.81412e+05);
    refforces[274] = Vec3(-1.35893e+05,  1.74355e+05, -1.12817e+05);
    refforces[275] = Vec3( 2.51361e+05, -1.90905e+04, -1.68633e+05);
    refforces[276] = Vec3( 1.76202e+05, -2.31601e+05, -2.00886e+05);
    refforces[277] = Vec3(-2.96697e+05,  3.00025e+04,  1.06676e+05);
    refforces[278] = Vec3( 1.20457e+05,  2.01635e+05,  9.42700e+04);
    refforces[279] = Vec3(-1.66320e+05, -9.08355e+02,  2.65478e+05);
    refforces[280] = Vec3( 2.44823e+05,  9.45895e+04, -5.28589e+04);
    refforces[281] = Vec3(-7.84765e+04, -9.36654e+04, -2.12648e+05);
    refforces[282] = Vec3(-6.58002e+03, -1.21918e+05,  3.16459e+05);
    refforces[283] = Vec3( 2.15225e+05,  5.96416e+04, -1.14097e+05);
    refforces[284] = Vec3(-2.08615e+05,  6.22725e+04, -2.02387e+05);
    refforces[285] = Vec3( 7.09186e+04,  1.83619e+05,  2.71862e+05);
    refforces[286] = Vec3( 1.78911e+05, -1.02234e+05, -1.78910e+05);
    refforces[287] = Vec3(-2.49882e+05, -8.13582e+04, -9.29803e+04);
    refforces[288] = Vec3(-1.35529e+04, -3.20226e+05, -1.16293e+05);
    refforces[289] = Vec3( 2.28053e+05,  1.65034e+05,  5.70125e+04);
    refforces[290] = Vec3(-2.14515e+05,  1.55237e+05,  5.92732e+04);
    refforces[291] = Vec3(-2.65006e+05,  1.75233e+05,  1.91268e+05);
    refforces[292] = Vec3( 2.29905e+05,  1.02940e+05, -2.05886e+05);
    refforces[293] = Vec3( 3.51345e+04, -2.78156e+05,  1.46393e+04);
    refforces[294] = Vec3( 1.59458e+05,  2.25128e+05,  2.11611e+05);
    refforces[295] = Vec3( 1.24813e+05, -1.16492e+05, -2.05249e+05);
    refforces[296] = Vec3(-2.84281e+05, -1.08601e+05, -6.38824e+03);
    refforces[297] = Vec3( 2.61767e+05, -7.55993e+04, -2.32576e+05);
    refforces[298] = Vec3(-2.07803e+05, -1.78116e+05,  7.71831e+04);
    refforces[299] = Vec3(-5.39242e+04,  2.53755e+05,  1.55427e+05);
    refforces[300] = Vec3(-6.44117e+03,  2.31328e+05, -2.54922e+05);
    refforces[301] = Vec3(-2.61097e+04,  4.93208e+04,  2.72709e+05);
    refforces[302] = Vec3( 3.25053e+04, -2.80718e+05, -1.77304e+04);
    refforces[303] = Vec3( 7.06280e+04,  7.26140e+04,  3.28446e+05);
    refforces[304] = Vec3( 1.45890e+05,  7.98271e+04, -2.06449e+05);
    refforces[305] = Vec3(-2.16516e+05, -1.52473e+05, -1.21980e+05);
    refforces[306] = Vec3( 1.94879e+05, -2.83084e+05,  1.41824e+05);
    refforces[307] = Vec3(-2.57149e+05,  1.96432e+05,  9.64284e+04);
    refforces[308] = Vec3( 6.22633e+04,  8.66280e+04, -2.38228e+05);
    refforces[309] = Vec3( 3.26456e+04,  1.17141e+05, -3.28374e+05);
    refforces[310] = Vec3( 2.01298e+05, -1.11486e+05,  1.85810e+05);
    refforces[311] = Vec3(-2.33937e+05, -5.70482e+03,  1.42641e+05);
    refforces[312] = Vec3( 6.42964e+04,  1.88719e+05, -2.86572e+05);
    refforces[313] = Vec3( 1.22184e+05, -2.53768e+05,  1.00254e+05);
    refforces[314] = Vec3(-1.86433e+05,  6.49694e+04,  1.86429e+05);
    refforces[315] = Vec3(-4.12317e+04, -1.73607e+04, -3.68620e+05);
    refforces[316] = Vec3(-1.78981e+05,  1.43186e+05,  2.08269e+05);
    refforces[317] = Vec3( 2.20155e+05, -1.25802e+05,  1.60395e+05);
    refforces[318] = Vec3(-1.58613e+05,  3.01594e+05, -1.29345e+05);
    refforces[319] = Vec3( 2.79698e+05, -2.67790e+04,  4.76063e+04);
    refforces[320] = Vec3(-1.21063e+05, -2.74847e+05,  8.17976e+04);
    refforces[321] = Vec3( 1.00451e+05,  2.03427e+05, -2.75044e+05);
    refforces[322] = Vec3(-2.13956e+04,  7.64144e+04,  2.81199e+05);
    refforces[323] = Vec3(-7.91039e+04, -2.79910e+05, -6.08590e+03);
    refforces[324] = Vec3(-2.80664e+05,  5.26142e+04, -1.53706e+05);
    refforces[325] = Vec3( 2.23923e+05, -1.17559e+05, -9.23685e+04);
    refforces[326] = Vec3( 5.68062e+04,  6.49218e+04,  2.46158e+05);
    refforces[327] = Vec3( 2.88911e+05, -1.17893e+05, -7.40828e+04);
    refforces[328] = Vec3(-6.38612e+04,  2.37565e+05, -2.55450e+04);
    refforces[329] = Vec3(-2.25033e+05, -1.19637e+05,  9.96952e+04);
    refforces[330] = Vec3( 1.03510e+04,  7.35559e+04,  3.47129e+05);
    refforces[331] = Vec3( 1.26778e+05, -2.19542e+05, -1.51515e+05);
    refforces[332] = Vec3(-1.37190e+05,  1.45950e+05, -1.95573e+05);
    refforces[333] = Vec3(-1.31816e+05,  5.57951e+04, -2.81935e+05);
    refforces[334] = Vec3(-1.08367e+05, -9.75306e+04,  2.16731e+05);
    refforces[335] = Vec3( 2.40191e+05,  4.17715e+04,  6.52667e+04);
    refforces[336] = Vec3(-2.86680e+05, -2.63007e+05,  1.37964e+04);
    refforces[337] = Vec3( 1.03914e+05,  1.76653e+05, -2.56324e+05);
    refforces[338] = Vec3( 1.82779e+05,  8.64040e+04,  2.42594e+05);
    refforces[339] = Vec3( 1.09800e+03, -3.11636e+05,  1.85833e+05);
    refforces[340] = Vec3(-2.39750e+05,  1.52568e+05, -8.71817e+04);
    refforces[341] = Vec3( 2.38635e+05,  1.59089e+05, -9.86351e+04);
    refforces[342] = Vec3(-8.76706e+04, -6.10520e+04,  3.31680e+05);
    refforces[343] = Vec3( 2.64692e+05,  5.23233e+04, -1.16959e+05);
    refforces[344] = Vec3(-1.76973e+05,  8.70353e+03, -2.14689e+05);
    refforces[345] = Vec3(-1.15976e+05, -2.72303e+05, -1.33307e+05);
    refforces[346] = Vec3( 2.20682e+05,  1.51896e+05, -6.30516e+04);
    refforces[347] = Vec3(-1.04744e+05,  1.20457e+05,  1.96395e+05);
    refforces[348] = Vec3( 4.57064e+04,  3.43550e+05,  7.23245e+04);
    refforces[349] = Vec3(-1.91414e+05, -1.26580e+05, -1.85241e+05);
    refforces[350] = Vec3( 1.45684e+05, -2.17040e+05,  1.12977e+05);
    refforces[351] = Vec3(-1.88623e+05, -1.22936e+04,  3.10223e+05);
    refforces[352] = Vec3( 3.06931e+05,  5.73380e+04, -7.08300e+04);
    refforces[353] = Vec3(-1.18279e+05, -4.50586e+04, -2.39376e+05);
    refforces[354] = Vec3( 8.31603e+04, -2.75958e+05,  1.16877e+05);
    refforces[355] = Vec3(-2.08784e+05,  1.26300e+05,  5.15494e+04);
    refforces[356] = Vec3( 1.25608e+05,  1.49661e+05, -1.68370e+05);
    refforces[357] = Vec3( 1.44100e+05, -2.34453e+05,  1.73916e+05);
    refforces[358] = Vec3(-8.65072e+04,  2.39988e+05,  8.37141e+04);
    refforces[359] = Vec3(-5.75432e+04, -5.48037e+03, -2.57574e+05);
    refforces[360] = Vec3( 2.72504e+05,  2.99781e+04,  1.85778e+05);
    refforces[361] = Vec3(-3.41026e+04,  7.10483e+04, -2.61466e+05);
    refforces[362] = Vec3(-2.38458e+05, -1.00996e+05,  7.57451e+04);
    refforces[363] = Vec3( 1.45826e+05, -2.81284e+05,  1.15502e+05);
    refforces[364] = Vec3(-2.63407e+05,  3.47332e+04, -8.10477e+04);
    refforces[365] = Vec3( 1.17594e+05,  2.46657e+05, -3.44181e+04);
    refforces[366] = Vec3( 2.59381e+05, -5.73297e+04,  2.56690e+05);
    refforces[367] = Vec3(-1.25834e+05, -2.09725e+05, -1.49803e+05);
    refforces[368] = Vec3(-1.33554e+05,  2.67105e+05, -1.06844e+05);
    refforces[369] = Vec3( 7.18473e+04, -2.59021e+05,  1.89797e+05);
    refforces[370] = Vec3( 1.56664e+05,  2.11783e+05, -8.99369e+04);
    refforces[371] = Vec3(-2.28511e+05,  4.72761e+04, -9.98083e+04);
    refforces[372] = Vec3(-1.99083e+04,  3.17042e+05, -5.62792e+04);
    refforces[373] = Vec3(-1.96875e+05, -1.64948e+05, -7.96758e-01);
    refforces[374] = Vec3( 2.16840e+05, -1.52072e+05,  5.63215e+04);


    const double longcutoff = 30.0*OpenMM::NmPerAngstrom;
    const double longcutoff2 = longcutoff*longcutoff;
    const double cutoff2 = cutoff*cutoff;
    const double cutoff6inv = 1.0 / (cutoff2*cutoff2*cutoff2);
    const int nboxes = ceil(longcutoff/boxEdgeLength);

    double refenergy = 0.0;
    // Loop over home box first...
    for (int i = 0; i < numAtoms; ++ i) {
        for (int j = i+1; j < numAtoms; ++j) {
            Vec3 dR = positions[i] - positions[j];
            double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
            double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
            double sig6 = sig2*sig2*sig2;
            double eps = epsvals[i]*epsvals[j];
            refenergy += 2.0*eps*(sig6-1.0)*sig6;
            if (R2 < cutoff2) {
                // Add a shift term for direct space parts withing t
                refenergy += 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
            }
        }
    }

    // ... and now add in the image box terms
    for (int bx = -nboxes; bx <= nboxes; ++bx) {
        for (int by = -nboxes; by <= nboxes; ++by) {
            for (int bz = -nboxes; bz <= nboxes; ++bz) {
                if (bx==0 && by==0 && bz==0) continue;
                Vec3 offset(bx*boxEdgeLength, by*boxEdgeLength, bz*boxEdgeLength);
                for (int i = 0; i < numAtoms; ++ i) {
                    for (int j = 0; j < numAtoms; ++j) {
                        Vec3 dR = positions[i] - positions[j] + offset;
                        double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
                        if (R2 > longcutoff2) continue;
                        double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
                        double sig6 = sig2*sig2*sig2;
                        double eps = epsvals[i]*epsvals[j];
                        refenergy += eps*(sig6-1.0)*sig6;
                        if (R2 < cutoff2) {
                            // Add a shift term for direct space parts withing teh
                            refenergy += eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                        }
                    }
                }
            }
        }
    }
    refenergy *= 0.5;

    // For this test the reference energy is 545636 kJ/mol, while the difference between DPME and 30A cutoffs
    // is just 1.062 kJ/mol.  The difference is due to the fact that arithmetic mean combination rules are used
    // up to the cutoff, while the reciprocal space uses the geometric mean.  See DOI: 10.1021/acs.jctc.5b00726

    ASSERT_EQUAL_TOL(refenergy, energy, 5E-5);
    ASSERT_EQUAL_TOL(gromacs_energy, energy, 5E-4);

    for (int n = 0; n < numAtoms; ++n)
        ASSERT_EQUAL_VEC(refforces[n], forces[n], 5E-4);
}



void testWater125DpmeVsLongCutoffWithExclusions() {
    const double cutoff = 11.5*OpenMM::NmPerAngstrom;
    const double alpha = 0.45*OpenMM::AngstromsPerNm;
    const double dalpha = 4.2760443169;
    const int grid = 60;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int numAtoms = 375;
    double boxEdgeLength = 23.01*OpenMM::NmPerAngstrom;

    make_waterbox(numAtoms, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);

    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->createExceptionsFromBonds(bonds, 1.0, 1.0);
    forceField->setPMEParameters(alpha, grid, grid, grid);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    State state = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy();
    const vector<Vec3>& forces = state.getForces();
    // Gromacs reference values.  See comments in testWater2DpmeEnergiesForcesNoExclusions() for details.
    double gromacs_energy = -2.17481e+02;
    vector<Vec3> refforces(numAtoms);
    refforces[  0] = Vec3(-3.10130e+01, -6.01219e+01, -5.31302e+01);
    refforces[  1] = Vec3( 4.21376e+00,  9.13811e-01,  9.61319e-01);
    refforces[  2] = Vec3( 1.42363e+00,  2.97219e+00,  1.27185e+00);
    refforces[  3] = Vec3(-3.70172e+01, -2.95228e+01, -6.27987e+01);
    refforces[  4] = Vec3(-6.57386e-01,  2.39279e+00,  1.53243e+00);
    refforces[  5] = Vec3(-1.61426e+00,  7.11728e-01,  9.99752e-01);
    refforces[  6] = Vec3( 1.60813e+01, -6.77979e+01, -1.01302e+02);
    refforces[  7] = Vec3(-4.23118e+00,  9.25550e-01,  1.28003e+00);
    refforces[  8] = Vec3( 1.81589e+00,  1.10083e+00,  1.28615e+00);
    refforces[  9] = Vec3(-2.24427e+00, -5.71071e+01, -3.18064e+01);
    refforces[ 10] = Vec3( 3.41715e+00,  1.32918e+00,  1.27908e+00);
    refforces[ 11] = Vec3(-3.09930e+00,  9.91868e-01,  1.69245e+00);
    refforces[ 12] = Vec3( 6.08033e+01, -1.02990e+02, -8.26733e+01);
    refforces[ 13] = Vec3(-1.01346e+00,  1.35475e+00,  4.13725e+00);
    refforces[ 14] = Vec3(-6.71645e-01,  4.69536e-01,  4.64774e-01);
    refforces[ 15] = Vec3(-6.39056e+01, -9.85756e-01, -4.97837e+01);
    refforces[ 16] = Vec3( 9.97113e-01,  1.13195e-01,  4.49387e+00);
    refforces[ 17] = Vec3( 9.48864e-01,  4.22217e+00,  1.56364e+00);
    refforces[ 18] = Vec3( 2.18124e+01, -6.14078e+01, -8.10129e+01);
    refforces[ 19] = Vec3(-5.59701e-01, -1.47528e+00,  1.40157e+00);
    refforces[ 20] = Vec3( 4.29539e+00,  4.27738e-02,  1.03009e+00);
    refforces[ 21] = Vec3(-2.62377e+01,  2.55532e+01, -6.89576e+01);
    refforces[ 22] = Vec3( 1.17222e-01,  3.81618e+00,  1.74502e+00);
    refforces[ 23] = Vec3( 1.46448e-01,  1.06214e-01,  4.27946e+00);
    refforces[ 24] = Vec3(-1.90819e+00,  9.86513e+00, -1.05364e+02);
    refforces[ 25] = Vec3(-1.18537e+00,  1.75021e+00,  1.54448e+00);
    refforces[ 26] = Vec3( 2.96182e+00, -9.48323e-01,  1.61299e+00);
    refforces[ 27] = Vec3( 5.61257e+01,  7.69427e+01, -5.84534e+01);
    refforces[ 28] = Vec3(-8.11821e-01, -7.15455e-01,  6.60049e-01);
    refforces[ 29] = Vec3(-9.20901e-01,  3.76217e+00,  1.35685e+00);
    refforces[ 30] = Vec3(-5.03791e+01,  2.24324e+01, -4.31606e+01);
    refforces[ 31] = Vec3( 3.72530e+00,  3.48275e-01,  1.40060e+00);
    refforces[ 32] = Vec3( 7.61265e-01,  8.56416e-01,  6.23359e-01);
    refforces[ 33] = Vec3( 1.95683e+01,  4.16950e+01, -4.33266e+01);
    refforces[ 34] = Vec3(-2.19677e+00,  1.06287e+00,  1.58260e+00);
    refforces[ 35] = Vec3( 1.70923e+00,  2.09511e+00,  9.84492e-01);
    refforces[ 36] = Vec3(-7.84273e+01, -3.34396e+01, -3.66486e+01);
    refforces[ 37] = Vec3(-1.53567e-01,  1.74934e-01,  1.15266e+00);
    refforces[ 38] = Vec3(-3.45548e+00,  1.68633e-01,  1.15923e+00);
    refforces[ 39] = Vec3( 6.90725e+01, -2.09426e+01, -5.31532e+01);
    refforces[ 40] = Vec3( 1.33274e+00,  3.09504e-01,  1.40065e+00);
    refforces[ 41] = Vec3( 1.65828e+00, -1.96378e-01,  2.68350e+00);
    refforces[ 42] = Vec3( 3.39954e+01, -2.65418e+01, -6.81184e+01);
    refforces[ 43] = Vec3(-1.17491e+00, -2.88793e+00,  1.40541e+00);
    refforces[ 44] = Vec3(-9.63132e-01, -3.78178e-01,  1.06756e+00);
    refforces[ 45] = Vec3(-4.44610e+01, -1.81736e+01, -6.19409e+01);
    refforces[ 46] = Vec3( 1.48529e+00,  8.95548e-06,  4.60732e+00);
    refforces[ 47] = Vec3( 1.00451e+00, -3.88830e+00,  1.27726e+00);
    refforces[ 48] = Vec3(-4.23022e+00,  2.52963e+01, -3.04756e+01);
    refforces[ 49] = Vec3(-3.45228e+00,  7.61191e-01,  9.40731e-01);
    refforces[ 50] = Vec3( 2.78622e+00,  8.32081e-01,  1.37588e+00);
    refforces[ 51] = Vec3( 1.63360e+01,  3.02230e+01, -3.91446e+01);
    refforces[ 52] = Vec3(-2.64761e+00,  1.18509e+00,  9.81393e-01);
    refforces[ 53] = Vec3( 1.74883e+00, -1.51958e-01,  2.56504e+00);
    refforces[ 54] = Vec3(-3.01111e+01,  4.08021e+01, -3.45058e+01);
    refforces[ 55] = Vec3(-1.17129e-01,  9.74282e-02,  4.56402e+00);
    refforces[ 56] = Vec3(-3.03696e+00,  4.44012e-01,  1.44269e+00);
    refforces[ 57] = Vec3( 7.57956e+01, -1.12043e+01, -3.92084e+01);
    refforces[ 58] = Vec3(-1.55617e+00,  2.26444e-01,  2.88234e+00);
    refforces[ 59] = Vec3(-1.37287e+00, -3.79045e+00,  8.55326e-01);
    refforces[ 60] = Vec3(-2.88831e+01,  5.72484e+01, -7.17030e+01);
    refforces[ 61] = Vec3( 6.88458e-01, -2.06343e+00,  1.03552e+00);
    refforces[ 62] = Vec3( 6.94194e-01, -1.83673e+00,  2.90979e+00);
    refforces[ 63] = Vec3(-7.87489e+01,  2.17729e+01, -3.72371e+01);
    refforces[ 64] = Vec3( 1.39267e-01, -4.55150e+00,  9.53893e-01);
    refforces[ 65] = Vec3(-4.27928e+00, -1.14139e+00,  9.25405e-01);
    refforces[ 66] = Vec3( 7.13947e+01,  3.02722e+01, -5.16551e+01);
    refforces[ 67] = Vec3( 4.07933e+00, -1.28199e+00,  9.34707e-01);
    refforces[ 68] = Vec3(-7.53668e-02, -2.14888e+00,  1.27713e+00);
    refforces[ 69] = Vec3(-9.36762e+00,  2.98321e+01, -5.70761e+01);
    refforces[ 70] = Vec3( 1.35989e+00, -3.09988e+00,  9.89945e-01);
    refforces[ 71] = Vec3(-2.20314e+00, -7.81892e-01,  1.21781e+00);
    refforces[ 72] = Vec3( 4.89194e+01,  7.40749e+01, -4.22017e+01);
    refforces[ 73] = Vec3(-4.25696e+00, -1.18954e+00,  1.04383e+00);
    refforces[ 74] = Vec3(-1.40510e+00, -9.87719e-01,  3.30559e+00);
    refforces[ 75] = Vec3(-4.95912e+01, -8.24838e+01, -1.82510e+01);
    refforces[ 76] = Vec3( 7.39933e-01,  9.29483e-01, -1.09883e+00);
    refforces[ 77] = Vec3( 9.70359e-01,  1.03012e+00,  4.01759e+00);
    refforces[ 78] = Vec3(-7.87499e+00, -1.87981e+01,  5.77017e+00);
    refforces[ 79] = Vec3( 3.82847e+00,  1.01191e+00, -7.30715e-02);
    refforces[ 80] = Vec3(-2.64831e+00,  1.68402e+00,  1.06474e-01);
    refforces[ 81] = Vec3( 1.81335e+01, -6.72766e+01,  7.21397e+01);
    refforces[ 82] = Vec3(-1.89188e-01,  1.27803e+00,  8.45103e-01);
    refforces[ 83] = Vec3( 2.47584e+00,  1.54044e+00,  6.55440e-01);
    refforces[ 84] = Vec3(-1.05824e+00, -6.50833e+01, -3.88413e+01);
    refforces[ 85] = Vec3(-3.21762e+00,  1.58809e+00, -1.32879e-01);
    refforces[ 86] = Vec3( 3.30285e+00,  1.52154e+00, -3.56471e-01);
    refforces[ 87] = Vec3( 4.41972e+01, -3.58815e+01,  4.01227e+01);
    refforces[ 88] = Vec3(-1.22923e+00,  8.91826e-01,  2.84505e+00);
    refforces[ 89] = Vec3(-4.31649e+00,  1.00895e+00,  4.43366e-01);
    refforces[ 90] = Vec3(-3.65176e+01,  4.86616e+01, -6.44163e+00);
    refforces[ 91] = Vec3( 1.03989e+00,  3.91910e+00,  7.93373e-01);
    refforces[ 92] = Vec3( 4.11089e+00, -4.60733e-01, -7.78461e-02);
    refforces[ 93] = Vec3(-1.58332e+01, -6.82045e+01,  6.62934e+01);
    refforces[ 94] = Vec3(-7.75645e-02, -3.66897e+00, -3.55044e-01);
    refforces[ 95] = Vec3( 4.08399e-01,  2.84950e-01,  3.80641e+00);
    refforces[ 96] = Vec3(-5.05975e+00,  4.20346e+01,  3.49518e+01);
    refforces[ 97] = Vec3( 9.82673e-01,  2.64451e+00,  4.02971e-01);
    refforces[ 98] = Vec3(-1.24992e+00, -3.34479e-01,  1.95012e+00);
    refforces[ 99] = Vec3( 2.31815e+00, -2.30056e+00,  6.49918e+01);
    refforces[100] = Vec3(-2.16986e+00, -1.05184e-01,  1.89395e+00);
    refforces[101] = Vec3( 3.51302e+00,  1.83706e-01,  7.08534e-01);
    refforces[102] = Vec3( 4.46929e+01, -7.11604e+01, -6.62201e+00);
    refforces[103] = Vec3(-1.87730e+00, -1.10206e+00, -1.47504e+00);
    refforces[104] = Vec3(-1.41639e+00, -1.36347e+00,  6.10811e-01);
    refforces[105] = Vec3(-3.12256e+01, -1.73923e+01, -4.47428e+01);
    refforces[106] = Vec3( 7.92004e-01, -1.39066e+00, -2.68309e+00);
    refforces[107] = Vec3( 1.54891e+00,  2.50018e+00, -6.21395e-01);
    refforces[108] = Vec3(-3.83585e+01,  3.90231e+01, -3.68876e+01);
    refforces[109] = Vec3( 1.86446e-01,  3.49120e+00, -4.00101e-01);
    refforces[110] = Vec3(-4.05592e+00, -2.21609e-01,  2.13257e-01);
    refforces[111] = Vec3( 5.70991e+01, -5.95502e+01, -4.08285e+01);
    refforces[112] = Vec3(-4.25039e-01, -3.96206e-01, -3.23550e+00);
    refforces[113] = Vec3( 4.14811e+00, -1.47859e-01,  1.11071e-01);
    refforces[114] = Vec3(-2.74979e+01,  4.40260e+01,  3.06766e+00);
    refforces[115] = Vec3( 2.22804e-01,  3.48934e+00, -1.96162e-01);
    refforces[116] = Vec3(-4.21754e+00,  1.18655e-01,  3.36974e-01);
    refforces[117] = Vec3( 6.32783e+01,  8.02719e+01,  1.38219e+01);
    refforces[118] = Vec3(-1.49058e+00,  1.87045e+00, -1.32164e-01);
    refforces[119] = Vec3(-2.30776e+00,  1.73240e-01,  1.92445e+00);
    refforces[120] = Vec3(-4.04272e+01, -1.77143e+01,  1.57461e+01);
    refforces[121] = Vec3( 1.82378e+00,  2.60376e+00,  2.33866e-01);
    refforces[122] = Vec3( 1.33379e+00, -1.16675e+00,  3.34265e+00);
    refforces[123] = Vec3(-8.78170e+00,  1.35313e+01, -5.40700e+01);
    refforces[124] = Vec3(-2.72817e+00,  3.27413e-01, -6.95012e-01);
    refforces[125] = Vec3( 2.41935e+00, -2.66049e-01, -1.62940e+00);
    refforces[126] = Vec3(-3.63273e+01,  6.91610e+01, -3.72663e+00);
    refforces[127] = Vec3(-1.40546e+00,  7.48503e-01, -9.73230e-01);
    refforces[128] = Vec3( 5.39530e-01,  8.56302e-01,  2.57617e+00);
    refforces[129] = Vec3( 2.54426e+01, -5.04169e+01, -1.76229e+01);
    refforces[130] = Vec3( 6.82595e-01, -3.17362e+00,  2.20610e-01);
    refforces[131] = Vec3( 1.94539e+00,  1.53325e+00, -4.24624e-01);
    refforces[132] = Vec3( 3.61678e+01, -7.89030e+01, -9.61808e+00);
    refforces[133] = Vec3(-1.29295e+00, -1.55223e+00, -1.88241e+00);
    refforces[134] = Vec3(-1.51438e+00,  3.14594e-02,  3.60043e+00);
    refforces[135] = Vec3(-4.60860e+01,  6.18971e+01,  1.38188e+01);
    refforces[136] = Vec3( 3.07506e+00, -9.60659e-01,  1.81485e+00);
    refforces[137] = Vec3( 6.46515e-01, -7.30203e-01, -2.17471e-01);
    refforces[138] = Vec3( 7.76539e+00,  3.18751e+01, -6.18739e+01);
    refforces[139] = Vec3(-2.12790e+00, -1.06577e+00, -1.03098e+00);
    refforces[140] = Vec3( 3.70469e+00, -7.70591e-01, -2.00619e-01);
    refforces[141] = Vec3(-4.02748e+01,  1.93804e+01,  1.92911e+01);
    refforces[142] = Vec3( 5.74429e-01, -1.58840e+00, -2.07083e+00);
    refforces[143] = Vec3(-4.21689e+00, -7.25280e-01,  1.61839e-01);
    refforces[144] = Vec3( 4.10890e+01,  6.50772e+01, -6.64661e+00);
    refforces[145] = Vec3( 2.02022e-01, -1.46124e+00,  1.99785e+00);
    refforces[146] = Vec3( 1.91733e-01, -1.14906e+00, -4.18291e+00);
    refforces[147] = Vec3( 3.77785e+01,  1.01191e+02, -6.07436e+00);
    refforces[148] = Vec3(-1.03030e+00, -1.26650e+00,  3.57188e+00);
    refforces[149] = Vec3(-7.42295e-01, -1.05437e+00, -1.47491e+00);
    refforces[150] = Vec3(-7.40150e+01, -3.53477e+01,  1.81946e+01);
    refforces[151] = Vec3( 1.65267e+00,  1.39873e+00,  1.46780e+00);
    refforces[152] = Vec3( 1.18632e+00,  2.66165e+00,  1.61236e+00);
    refforces[153] = Vec3( 3.89536e+01, -2.96071e+01, -1.38610e+01);
    refforces[154] = Vec3( 2.05061e+00,  1.54050e+00,  5.35476e-01);
    refforces[155] = Vec3( 4.98336e-01,  3.85036e+00,  4.21301e-01);
    refforces[156] = Vec3(-5.53456e+01, -6.87157e+01, -3.55693e+01);
    refforces[157] = Vec3( 3.17070e+00,  1.65800e+00, -3.99476e-01);
    refforces[158] = Vec3(-2.02395e+00,  1.48143e+00,  1.00879e+00);
    refforces[159] = Vec3( 5.94087e+01, -4.57694e+01,  2.77585e+01);
    refforces[160] = Vec3( 2.05570e+00,  1.46552e+00, -1.41447e+00);
    refforces[161] = Vec3( 2.97741e-01,  9.12876e-01,  3.27902e+00);
    refforces[162] = Vec3( 1.93442e+01, -4.16391e+01,  3.34803e-01);
    refforces[163] = Vec3(-3.93865e+00,  8.46322e-01, -4.47285e-01);
    refforces[164] = Vec3(-7.34898e-01,  9.82679e-01,  1.44017e+00);
    refforces[165] = Vec3(-3.34824e+01, -1.38834e+01, -2.62952e+01);
    refforces[166] = Vec3( 1.26693e+00, -1.10171e+00,  3.22806e+00);
    refforces[167] = Vec3( 3.49758e+00, -1.84534e-01, -4.62412e-01);
    refforces[168] = Vec3( 1.82122e+00, -5.69776e+01, -3.86605e+01);
    refforces[169] = Vec3( 1.92404e+00, -1.05303e+00, -4.52186e-01);
    refforces[170] = Vec3(-3.25172e+00,  9.44050e-02, -5.87318e-01);
    refforces[171] = Vec3(-5.98475e+01,  1.85708e+01, -2.67060e+01);
    refforces[172] = Vec3(-3.77230e-01,  3.29609e+00, -4.81590e-01);
    refforces[173] = Vec3(-1.66002e+00, -2.41098e+00,  1.03425e-01);
    refforces[174] = Vec3( 6.47852e+01,  4.38661e+00,  1.77482e+01);
    refforces[175] = Vec3(-1.29162e-01, -3.40031e+00,  3.91521e-01);
    refforces[176] = Vec3( 1.28253e+00,  2.14633e+00,  2.64792e-01);
    refforces[177] = Vec3( 2.95010e+01, -5.60376e+00,  4.16733e+01);
    refforces[178] = Vec3(-1.43378e+00, -1.96178e-01,  2.40594e+00);
    refforces[179] = Vec3(-1.56255e+00, -3.08967e+00, -3.14335e-01);
    refforces[180] = Vec3(-6.02384e+01, -1.59585e+01,  7.16798e+01);
    refforces[181] = Vec3( 1.68669e+00, -3.54228e+00, -2.98122e-01);
    refforces[182] = Vec3( 9.46412e-01, -8.03852e-02,  4.17305e+00);
    refforces[183] = Vec3( 3.87313e+00,  4.46632e+01,  4.59130e+01);
    refforces[184] = Vec3(-9.83299e-02, -8.11669e-02,  4.04515e+00);
    refforces[185] = Vec3( 6.64172e-01,  3.66814e+00, -2.58278e-01);
    refforces[186] = Vec3(-1.65450e+01,  2.35336e+01,  7.09485e+01);
    refforces[187] = Vec3(-6.98821e-01,  2.45516e-01,  2.69560e+00);
    refforces[188] = Vec3( 4.04145e+00,  5.46258e-02, -1.82210e-02);
    refforces[189] = Vec3( 2.39311e+01,  1.22785e+01, -2.13260e+01);
    refforces[190] = Vec3(-3.33802e-02, -3.60205e+00,  6.21747e-01);
    refforces[191] = Vec3( 3.50616e+00,  6.09623e-01, -7.40284e-02);
    refforces[192] = Vec3( 2.32693e+01, -1.88836e+01,  1.24193e+01);
    refforces[193] = Vec3(-2.25948e+00, -1.74718e+00,  2.65705e-01);
    refforces[194] = Vec3(-1.32949e+00, -8.67446e-01,  7.75386e-01);
    refforces[195] = Vec3(-4.40626e+01,  1.17032e+01, -4.18284e+01);
    refforces[196] = Vec3( 1.61723e+00, -6.18309e-01,  2.45312e+00);
    refforces[197] = Vec3( 1.42535e+00,  1.71801e+00, -1.86477e+00);
    refforces[198] = Vec3( 3.37157e+00, -1.65818e+01,  7.59773e+01);
    refforces[199] = Vec3(-2.81989e+00, -1.37877e+00, -9.60015e-02);
    refforces[200] = Vec3(-8.77602e-02,  2.44299e+00,  1.17362e+00);
    refforces[201] = Vec3( 2.44644e+01, -3.26831e+01, -5.34506e+01);
    refforces[202] = Vec3(-4.14853e+00, -9.48873e-02, -2.18126e-01);
    refforces[203] = Vec3( 8.13810e-01, -1.30429e+00, -1.62539e+00);
    refforces[204] = Vec3( 1.03414e+01, -4.14898e+00,  2.56841e+01);
    refforces[205] = Vec3(-1.50990e-01, -3.90514e+00,  2.98649e-01);
    refforces[206] = Vec3(-3.60035e+00, -7.83546e-02, -7.95300e-01);
    refforces[207] = Vec3( 4.49059e+01,  2.88959e+01,  1.53626e+01);
    refforces[208] = Vec3(-2.55855e+00, -1.55231e+00, -3.96734e-01);
    refforces[209] = Vec3(-1.14072e+00,  3.25255e+00,  2.39860e-01);
    refforces[210] = Vec3(-4.12850e+01,  6.11654e+01,  2.43179e+01);
    refforces[211] = Vec3( 1.20781e+00, -8.28196e-01,  4.26212e+00);
    refforces[212] = Vec3( 7.30387e-01, -9.26965e-01, -9.32605e-01);
    refforces[213] = Vec3(-1.64972e+01,  4.58667e+01,  5.58154e+01);
    refforces[214] = Vec3( 4.43484e-01, -1.67205e+00,  3.42164e+00);
    refforces[215] = Vec3(-4.14170e+00, -9.68400e-01, -9.21231e-02);
    refforces[216] = Vec3( 1.79252e+01,  6.24205e+01, -1.47777e+01);
    refforces[217] = Vec3(-2.53437e+00, -9.16373e-01, -1.63636e+00);
    refforces[218] = Vec3( 3.32247e+00, -1.29776e+00, -1.07427e+00);
    refforces[219] = Vec3( 1.79752e+01,  5.16430e+01,  3.90249e+01);
    refforces[220] = Vec3( 7.14180e-01, -1.00018e+00,  3.14118e+00);
    refforces[221] = Vec3(-1.07965e+00, -3.41285e+00, -2.24161e-01);
    refforces[222] = Vec3( 3.16317e+01,  4.98025e+01, -4.17518e+01);
    refforces[223] = Vec3(-2.55899e+00, -1.63433e+00, -3.76149e-01);
    refforces[224] = Vec3(-1.19443e+00, -1.02227e+00, -5.37349e-02);
    refforces[225] = Vec3(-3.77249e+01, -4.85845e+01, -2.11073e+01);
    refforces[226] = Vec3( 7.98995e-01,  1.26056e+00,  4.28442e+00);
    refforces[227] = Vec3( 3.35648e+00,  9.11653e-01, -8.07025e-01);
    refforces[228] = Vec3(-3.16004e+01, -2.04156e+01,  5.01542e+01);
    refforces[229] = Vec3( 8.97766e-01,  1.34997e+00,  2.86585e+00);
    refforces[230] = Vec3(-2.09603e-01,  4.05651e+00, -1.91043e-02);
    refforces[231] = Vec3( 3.61582e+01, -3.35593e+01,  4.19940e+01);
    refforces[232] = Vec3( 3.66075e+00,  7.71021e-01,  3.45598e-01);
    refforces[233] = Vec3(-7.50361e-01,  1.22549e+00,  7.98677e-01);
    refforces[234] = Vec3(-4.42511e+01, -2.71049e+01, -4.56576e+01);
    refforces[235] = Vec3(-3.30181e-01,  1.98824e+00,  2.21005e+00);
    refforces[236] = Vec3(-7.68719e-01,  9.09063e-01, -3.77633e+00);
    refforces[237] = Vec3( 7.33817e+01, -3.03886e+01, -7.59758e+01);
    refforces[238] = Vec3(-1.71701e+00,  1.41048e+00, -1.84182e+00);
    refforces[239] = Vec3(-1.35491e+00,  2.27816e+00, -1.08799e+00);
    refforces[240] = Vec3(-4.58218e+01, -2.91931e+01,  1.12369e+01);
    refforces[241] = Vec3( 1.68546e+00,  2.81258e+00,  3.93755e-01);
    refforces[242] = Vec3( 8.28917e-01, -1.32237e+00,  2.91182e+00);
    refforces[243] = Vec3( 5.44446e+00, -3.99466e+01,  6.77923e+00);
    refforces[244] = Vec3(-1.09231e+00, -5.21927e-01,  1.83581e+00);
    refforces[245] = Vec3( 5.24849e-02, -7.21832e-01, -3.35896e+00);
    refforces[246] = Vec3(-1.29863e+01, -3.94851e+01, -2.70347e+01);
    refforces[247] = Vec3(-3.94725e+00, -2.01688e-01,  5.93703e-02);
    refforces[248] = Vec3( 1.30447e+00, -2.35362e+00, -2.94278e-01);
    refforces[249] = Vec3( 2.56542e+01, -3.79434e+01, -6.19782e+01);
    refforces[250] = Vec3( 2.48025e+00, -1.19994e+00,  1.41640e-01);
    refforces[251] = Vec3(-5.01379e-01,  2.48607e-01, -3.97535e+00);
    refforces[252] = Vec3( 3.46461e+01, -8.37214e+01, -3.92225e+01);
    refforces[253] = Vec3(-2.80216e+00, -1.13524e+00,  1.85582e-02);
    refforces[254] = Vec3(-1.28421e+00, -9.50174e-01, -5.56828e-01);
    refforces[255] = Vec3(-3.56755e+01,  1.44176e+01, -7.63345e-01);
    refforces[256] = Vec3( 1.25547e+00,  2.64846e+00,  1.29226e+00);
    refforces[257] = Vec3( 1.22691e+00,  5.72301e-01, -4.09853e+00);
    refforces[258] = Vec3(-6.02997e+00, -1.31579e+01, -5.03304e+01);
    refforces[259] = Vec3(-4.59947e-01, -8.09031e-01,  3.08383e+00);
    refforces[260] = Vec3(-8.93872e-01, -3.70797e-01, -2.89773e+00);
    refforces[261] = Vec3(-1.92458e+00,  2.50787e+01, -4.30937e+01);
    refforces[262] = Vec3( 2.10381e-01, -1.24125e-01, -4.01554e+00);
    refforces[263] = Vec3(-3.85853e+00,  6.72311e-01,  4.06873e-02);
    refforces[264] = Vec3( 1.66952e+01,  1.47383e+01,  4.13458e+01);
    refforces[265] = Vec3(-2.12261e+00, -5.40916e-01,  1.54532e+00);
    refforces[266] = Vec3( 3.54041e+00, -4.73740e-02,  4.73489e-01);
    refforces[267] = Vec3( 4.43789e+01,  6.96271e+01,  4.06335e+00);
    refforces[268] = Vec3(-9.59563e-01,  9.94730e-01, -3.12464e+00);
    refforces[269] = Vec3(-1.42363e+00,  1.62107e+00,  1.87365e+00);
    refforces[270] = Vec3(-5.38380e+01,  2.27365e+01,  1.95666e+01);
    refforces[271] = Vec3( 1.47405e+00,  2.11523e+00,  1.97777e-01);
    refforces[272] = Vec3( 2.54377e+00, -4.04439e-02,  1.97927e+00);
    refforces[273] = Vec3( 1.33568e+01,  1.88750e+01, -3.59233e+01);
    refforces[274] = Vec3(-1.09293e+00,  2.25841e+00, -4.61222e-01);
    refforces[275] = Vec3( 3.27762e+00, -1.71208e-01, -8.67824e-01);
    refforces[276] = Vec3(-3.41864e+01,  3.38579e+01,  5.95219e+01);
    refforces[277] = Vec3(-4.18249e+00, -8.81154e-02,  1.97498e-01);
    refforces[278] = Vec3( 1.06177e+00,  2.34089e+00,  2.42037e-01);
    refforces[279] = Vec3( 2.35073e+01,  1.61625e+01, -2.58858e+01);
    refforces[280] = Vec3( 3.87903e+00,  1.67313e-01, -4.03473e-02);
    refforces[281] = Vec3(-5.12341e-01, -5.62857e-01, -3.19237e+00);
    refforces[282] = Vec3( 3.39154e+01, -5.09470e+00, -2.29146e+01);
    refforces[283] = Vec3(-1.24992e+00,  5.24960e-01, -1.16022e+00);
    refforces[284] = Vec3(-2.39537e+00,  1.82818e-01, -1.53421e+00);
    refforces[285] = Vec3(-5.56195e+01,  2.89266e+01, -2.61137e+01);
    refforces[286] = Vec3( 2.49596e+00, -1.03213e+00, -1.51258e+00);
    refforces[287] = Vec3( 1.02471e+00, -1.05342e+00, -6.79809e-01);
    refforces[288] = Vec3(-1.52895e+01,  4.83734e+01, -7.33233e+00);
    refforces[289] = Vec3( 3.77837e+00, -1.59224e+00,  2.02675e-01);
    refforces[290] = Vec3(-3.51854e+00, -1.54234e+00,  2.89348e-01);
    refforces[291] = Vec3( 3.22262e+01,  2.32542e+01,  2.37269e+01);
    refforces[292] = Vec3( 1.88520e+00, -1.38396e+00, -1.98750e+00);
    refforces[293] = Vec3( 4.99036e-02, -4.53129e+00, -2.13859e-01);
    refforces[294] = Vec3(-7.48274e+00,  3.71318e+01, -2.35667e+01);
    refforces[295] = Vec3( 7.81671e-01, -1.07371e+00, -2.37009e+00);
    refforces[296] = Vec3(-3.98896e+00, -9.36146e-01,  2.62843e-02);
    refforces[297] = Vec3( 4.33081e+01,  4.24806e+01,  3.21190e+01);
    refforces[298] = Vec3(-2.52531e+00, -1.90285e+00,  2.72091e-01);
    refforces[299] = Vec3(-9.70927e-01, -1.13219e+00,  1.47351e+00);
    refforces[300] = Vec3(-4.71070e+01, -7.19253e+01,  5.82986e+01);
    refforces[301] = Vec3( 6.82929e-01,  1.00679e+00, -6.94522e-01);
    refforces[302] = Vec3( 8.31359e-01,  7.09804e-01, -8.12652e-01);
    refforces[303] = Vec3( 3.62309e+00, -3.42931e+01,  1.98231e+01);
    refforces[304] = Vec3( 1.23138e+00,  8.49297e-01, -2.39283e+00);
    refforces[305] = Vec3(-2.87428e+00,  1.41050e+00, -1.03948e+00);
    refforces[306] = Vec3(-3.83149e+00, -2.61638e+01,  2.98550e+01);
    refforces[307] = Vec3(-2.90010e+00,  1.38356e+00, -1.22906e+00);
    refforces[308] = Vec3( 2.62337e-01,  9.72795e-01, -4.02660e+00);
    refforces[309] = Vec3( 7.84655e+00, -5.09531e+01,  7.96762e+01);
    refforces[310] = Vec3( 2.46967e+00,  9.48901e-01, -1.17317e+00);
    refforces[311] = Vec3(-4.11596e+00,  7.95633e-01, -1.39335e+00);
    refforces[312] = Vec3( 5.03705e+01, -8.09993e+01,  1.13172e+02);
    refforces[313] = Vec3(-4.48925e-01,  4.29309e-01, -5.50572e-01);
    refforces[314] = Vec3(-2.66201e+00,  1.05662e+00, -1.44247e+00);
    refforces[315] = Vec3(-6.12660e+01,  2.31218e+01,  4.57287e+01);
    refforces[316] = Vec3( 8.42017e-01,  1.41042e+00, -7.78795e-01);
    refforces[317] = Vec3( 3.11298e+00, -8.04273e-01, -1.51273e+00);
    refforces[318] = Vec3( 1.83823e+01, -2.83588e+01,  6.13597e+01);
    refforces[319] = Vec3( 4.37781e+00,  8.12091e-02, -1.33992e+00);
    refforces[320] = Vec3(-5.81396e-01, -3.61991e+00, -1.22489e+00);
    refforces[321] = Vec3(-4.85270e+01, -6.50314e+01,  7.16441e+01);
    refforces[322] = Vec3(-3.40093e-02,  9.01514e-01, -1.26573e+00);
    refforces[323] = Vec3(-8.67193e-02, -3.94828e+00, -8.97127e-01);
    refforces[324] = Vec3( 6.25833e+01, -2.31957e+01,  8.66089e+01);
    refforces[325] = Vec3( 2.49809e+00, -7.64197e-01, -1.27775e+00);
    refforces[326] = Vec3( 2.11019e-01,  7.81764e-01, -1.23664e+00);
    refforces[327] = Vec3( 2.04801e+01,  3.14002e+01,  6.96727e+01);
    refforces[328] = Vec3(-7.99932e-01,  4.35586e+00, -8.77957e-01);
    refforces[329] = Vec3(-2.87129e+00, -1.02321e+00, -1.42702e+00);
    refforces[330] = Vec3(-6.36347e+01, -3.58428e+01,  4.52682e+01);
    refforces[331] = Vec3( 1.32210e+00, -2.58500e+00, -1.43327e+00);
    refforces[332] = Vec3( 1.52218e+00,  1.63513e+00, -2.29784e+00);
    refforces[333] = Vec3( 4.36105e+00,  3.70776e+01,  6.46133e+01);
    refforces[334] = Vec3(-1.19452e+00, -1.01252e+00, -1.43172e+00);
    refforces[335] = Vec3( 4.34983e+00, -6.52342e-02, -1.33662e+00);
    refforces[336] = Vec3( 1.13050e+01,  4.99187e+01,  7.16119e+01);
    refforces[337] = Vec3( 2.29221e-01,  6.35843e-01, -3.59259e+00);
    refforces[338] = Vec3( 1.65239e+00,  4.20995e-01, -1.42878e+00);
    refforces[339] = Vec3(-1.79210e+01,  1.92423e+01,  1.77516e+01);
    refforces[340] = Vec3(-2.91278e+00,  9.22214e-01, -8.41281e-01);
    refforces[341] = Vec3( 2.96465e+00,  1.05962e+00, -6.85280e-01);
    refforces[342] = Vec3( 5.19378e+01, -2.55266e+01,  3.63518e+01);
    refforces[343] = Vec3(-1.13088e+00,  4.54600e-01, -1.30862e+00);
    refforces[344] = Vec3(-2.05369e+00, -9.69459e-03, -2.69351e+00);
    refforces[345] = Vec3(-4.22214e+01,  4.88282e+01,  3.76741e+01);
    refforces[346] = Vec3( 3.46996e+00,  9.57280e-01, -7.60128e-01);
    refforces[347] = Vec3( 7.62585e-01,  1.14548e+00, -8.06697e-01);
    refforces[348] = Vec3(-2.35024e+01, -6.73585e+01,  6.49050e+01);
    refforces[349] = Vec3(-1.58984e+00, -6.44551e-01, -2.37758e+00);
    refforces[350] = Vec3( 1.32604e+00, -2.37576e+00, -1.54939e+00);
    refforces[351] = Vec3( 2.57776e+01, -1.40743e+01,  2.06580e+01);
    refforces[352] = Vec3( 4.24916e+00,  6.10665e-02, -8.19360e-01);
    refforces[353] = Vec3(-7.94423e-01, -1.75764e-01, -3.51414e+00);
    refforces[354] = Vec3(-1.24260e+01,  6.96334e-01,  5.97705e+01);
    refforces[355] = Vec3(-3.38899e+00,  9.44919e-01, -1.25780e+00);
    refforces[356] = Vec3( 6.81664e-01,  1.16965e+00, -2.34256e+00);
    refforces[357] = Vec3( 5.15292e+01,  5.13086e+01,  6.16724e+01);
    refforces[358] = Vec3(-1.04216e+00,  3.34685e+00, -1.33787e+00);
    refforces[359] = Vec3(-9.57017e-01, -1.59610e-01, -4.31963e+00);
    refforces[360] = Vec3(-5.88173e+01,  3.30524e+01,  6.23031e+01);
    refforces[361] = Vec3( 9.35413e-01, -1.09829e+00, -4.65202e+00);
    refforces[362] = Vec3( 6.95250e-01, -1.24375e+00, -6.50832e-01);
    refforces[363] = Vec3( 1.64454e+01,  1.08071e+02,  3.76968e+01);
    refforces[364] = Vec3(-4.25249e+00, -1.10053e+00, -7.56357e-01);
    refforces[365] = Vec3( 1.17898e+00, -1.10215e+00, -8.21366e-01);
    refforces[366] = Vec3(-5.20527e+00,  5.49522e+01,  4.55303e+01);
    refforces[367] = Vec3(-7.98794e-01, -2.78246e+00, -1.27815e+00);
    refforces[368] = Vec3(-1.02109e+00, -1.21921e+00, -1.26706e+00);
    refforces[369] = Vec3( 3.05185e+00,  3.99846e+01,  5.39902e+01);
    refforces[370] = Vec3( 1.93675e+00, -1.28986e+00, -1.21655e+00);
    refforces[371] = Vec3(-4.14952e+00, -1.06509e+00, -9.52154e-01);
    refforces[372] = Vec3( 6.06000e+01,  2.50872e+01,  4.29623e+01);
    refforces[373] = Vec3(-3.50284e+00, -1.47848e+00, -7.98049e-01);
    refforces[374] = Vec3(-8.66087e-01, -1.73056e+00, -6.59554e-01);

    // Find the exclusion information
    vector<set<int> > exclusions(numAtoms);
    for (int i = 0; i < forceField->getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        forceField->getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    const double longcutoff = 35.0*OpenMM::NmPerAngstrom;
    const double longcutoff2 = longcutoff*longcutoff;
    const double cutoff2 = cutoff*cutoff;
    const double cutoff6inv = 1.0 / (cutoff2*cutoff2*cutoff2);
    const int nboxes = ceil(longcutoff/boxEdgeLength);

    double refenergy = 0.0;
    // Loop over home box first...
    for (int i = 0; i < numAtoms; ++ i) {
        for (int j = i+1; j < numAtoms; ++j) {
            Vec3 dR = positions[i] - positions[j];
            double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
            double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
            double sig6 = sig2*sig2*sig2;
            double eps = epsvals[i]*epsvals[j];
            refenergy += 2.0*eps*(sig6-1.0)*sig6;
            if (R2 < cutoff2) {
                // Add a shift term for direct space parts withing t
                refenergy += 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
            }
        }
    }

    // ... back out exclusions
    for (int ii = 0; ii < numAtoms; ii++) {
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
                if (R2 < cutoff2) {
                    // Add a shift term for direct space parts withing t
                    refenergy -= 2.0*eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                }
            }
        }
    }

    // ... and now add in the image box terms
    for (int bx = -nboxes; bx <= nboxes; ++bx) {
        for (int by = -nboxes; by <= nboxes; ++by) {
            for (int bz = -nboxes; bz <= nboxes; ++bz) {
                if (bx==0 && by==0 && bz==0) continue;
                Vec3 offset(bx*boxEdgeLength, by*boxEdgeLength, bz*boxEdgeLength);
                for (int i = 0; i < numAtoms; ++ i) {
                    for (int j = 0; j < numAtoms; ++j) {
                        Vec3 dR = positions[i] - positions[j] + offset;
                        double R2 = dR[0]*dR[0] + dR[1]*dR[1] + dR[2]*dR[2];
                        if (R2 > longcutoff2) continue;
                        double sig2 = (sigvals[i] + sigvals[j])*(sigvals[i] + sigvals[j]) / R2;
                        double sig6 = sig2*sig2*sig2;
                        double eps = epsvals[i]*epsvals[j];
                        refenergy += eps*(sig6-1.0)*sig6;
                        if (R2 < cutoff2) {
                            // Add a shift term for direct space parts withing teh
                            refenergy += eps*(pow(sigvals[i]+sigvals[j], 6) - 64.0*pow(sigvals[i]*sigvals[j], 3))*cutoff6inv;
                        }
                    }
                }
            }
        }
    }
    refenergy *= 0.5;
    // For this test the reference energy is -294.078 kJ/mol, while the difference between DPME and 30A cutoffs
    // is just 0.064 kJ/mol.  The difference is due to the fact that arithmetic mean combination rules are used
    // up to the cutoff, while the reciprocal space uses the geometric mean.  See DOI: 10.1021/acs.jctc.5b00726
    ASSERT_EQUAL_TOL(refenergy, energy, 5E-4);
    ASSERT_EQUAL_TOL(gromacs_energy, energy, 5E-5);

    // Forces accumulated in single precision are tested to a more permissive criterion; the double
    // precision platform can match to 5E-5.
    for (int n = 0; n < numAtoms; ++n)
        ASSERT_EQUAL_VEC(refforces[n], forces[n], 1E-4);

}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testDMFDpmeEnergiesForcesWithExclusions();
        testConvergence();
        testErrorTolerance();
        testPMEParameters();
        testCoulombAndLJ();
        testWater2DpmeEnergiesForcesNoExclusions();
        testWater2DpmeEnergiesForcesWithExclusions();
        testWater125DpmeVsLongCutoffNoExclusions();
        testWater125DpmeVsLongCutoffWithExclusions();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
