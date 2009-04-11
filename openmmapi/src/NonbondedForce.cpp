/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "Force.h"
#include "OpenMMException.h"
#include "NonbondedForce.h"
#include "internal/NonbondedForceImpl.h"

using namespace OpenMM;

NonbondedForce::NonbondedForce(int numParticles, int numNonbonded14) : particles(numParticles), nb14s(numNonbonded14),
        nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
    periodicBoxVectors[0] = Vec3(2, 0, 0);
    periodicBoxVectors[1] = Vec3(0, 2, 0);
    periodicBoxVectors[2] = Vec3(0, 0, 2);
}

NonbondedForce::NonbondedMethod NonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void NonbondedForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double NonbondedForce::getCutoffDistance() const {
    return cutoffDistance;
}

void NonbondedForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

void NonbondedForce::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void NonbondedForce::setPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c) {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw OpenMMException("First periodic box vector must be parallel to x.");
    if (b[0] != 0.0 || b[2] != 0.0)
        throw OpenMMException("Second periodic box vector must be parallel to y.");
    if (c[0] != 0.0 || c[1] != 0.0)
        throw OpenMMException("Third periodic box vector must be parallel to z.");
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

void NonbondedForce::getParticleParameters(int index, double& charge, double& sigma, double& epsilon) const {
    charge = particles[index].charge;
    sigma = particles[index].sigma;
    epsilon = particles[index].epsilon;
}

void NonbondedForce::setParticleParameters(int index, double charge, double sigma, double epsilon) {
    particles[index].charge = charge;
    particles[index].sigma = sigma;
    particles[index].epsilon = epsilon;
}

void NonbondedForce::getNonbonded14Parameters(int index, int& particle1, int& particle2, double& chargeProd, double& sigma, double& epsilon) const {
    particle1 = nb14s[index].particle1;
    particle2 = nb14s[index].particle2;
    chargeProd = nb14s[index].chargeProd;
    sigma = nb14s[index].sigma;
    epsilon = nb14s[index].epsilon;
}

void NonbondedForce::setNonbonded14Parameters(int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon) {
    nb14s[index].particle1 = particle1;
    nb14s[index].particle2 = particle2;
    nb14s[index].chargeProd = chargeProd;
    nb14s[index].sigma = sigma;
    nb14s[index].epsilon = epsilon;
}

ForceImpl* NonbondedForce::createImpl() {
    return new NonbondedForceImpl(*this);
}
