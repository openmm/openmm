/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2020 Stanford University and the Authors.      *
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

#include "sfmt/SFMT.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/Integrator.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"

#include <cmath>

using namespace OpenMM;

Integrator::Integrator() : owner(NULL), context(NULL), forceGroups(0xFFFFFFFF) {
}

Integrator::~Integrator() {
    if (context != NULL) {
        // The Integrator is being deleted before the Context, so do cleanup now,
        // then notify the ContextImpl so its own destructor won't try to clean up
        // the (no longer existing) Integrator.

        cleanup();
        context->integratorDeleted();
    }
}

double Integrator::getStepSize() const {
    return stepSize;
}

void Integrator::setStepSize(double size) {
    stepSize = size;
}

double Integrator::getConstraintTolerance() const {
    return constraintTol;
}

void Integrator::setConstraintTolerance(double tol) {
    constraintTol = tol;
}

int Integrator::getIntegrationForceGroups() const {
    return forceGroups;
}

void Integrator::setIntegrationForceGroups(int groups) {
    forceGroups = groups;
}

std::vector<Vec3> Integrator::getVelocitiesForTemperature(const System &system, double temperature, int randomSeed) const {
    // Generate the list of Gaussian random numbers.
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(randomSeed, sfmt);
    std::vector<double> randoms;
    while (randoms.size() < system.getNumParticles()*3) {
        double x, y, r2;
        do {
            x = 2.0*genrand_real2(sfmt)-1.0;
            y = 2.0*genrand_real2(sfmt)-1.0;
            r2 = x*x + y*y;
        } while (r2 >= 1.0 || r2 == 0.0);
        double multiplier = sqrt((-2.0*std::log(r2))/r2);
        randoms.push_back(x*multiplier);
        randoms.push_back(y*multiplier);
    }

    // Assign the velocities.
    std::vector<Vec3> velocities(system.getNumParticles(), Vec3());
    int nextRandom = 0;
    for (int i = 0; i < system.getNumParticles(); i++) {
        double mass = system.getParticleMass(i);
        if (mass != 0) {
            double velocityScale = sqrt(BOLTZ*temperature/mass);
            velocities[i] = Vec3(randoms[nextRandom++], randoms[nextRandom++], randoms[nextRandom++])*velocityScale;
        }
    }
    return velocities;
}

