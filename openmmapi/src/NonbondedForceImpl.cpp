/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace OpenMM;
using namespace std;

NonbondedForceImpl::NonbondedForceImpl(NonbondedForce& owner) : owner(owner) {
}

NonbondedForceImpl::~NonbondedForceImpl() {
}

void NonbondedForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcNonbondedForceKernel::Name(), context);

    // Check for errors in the specification of exceptions.

    System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("NonbondedForce must have exactly as many particles as the System it belongs to.");
    vector<set<int> > exceptions(owner.getNumParticles());
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        owner.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "NonbondedForce: Illegal particle index for an exception: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "NonbondedForce: Illegal particle index for an exception: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exceptions[particle1].count(particle2) > 0 || exceptions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "NonbondedForce: Multiple exceptions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[particle1].insert(particle2);
        exceptions[particle2].insert(particle1);
    }
    dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner);
}

void NonbondedForceImpl::calcForces(ContextImpl& context) {
    dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).executeForces(context);
}

double NonbondedForceImpl::calcEnergy(ContextImpl& context) {
    return dynamic_cast<CalcNonbondedForceKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> NonbondedForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcNonbondedForceKernel::Name());
    return names;
}

class NonbondedForceImpl::ErrorFunction {
public:
    virtual double getValue(int arg) const = 0;
};

class NonbondedForceImpl::EwaldErrorFunction : public ErrorFunction {
public:
    EwaldErrorFunction(double width, double alpha, double target) : width(width), alpha(alpha), target(target) {
    }
    double getValue(int arg) const {
        double temp = arg*M_PI/(width*alpha);
        return target-0.05*sqrt(width*alpha)*arg*exp(-temp*temp);
    }
private:
    double width, alpha, target;
};

void NonbondedForceImpl::calcEwaldParameters(const System& system, const NonbondedForce& force, double& alpha, int& kmaxx, int& kmaxy, int& kmaxz) {
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    double tol = force.getEwaldErrorTolerance();
    alpha = (1.0/force.getCutoffDistance())*std::sqrt(-log(2.0*tol));
    kmaxx = findZero(EwaldErrorFunction(boxVectors[0][0], alpha, tol), 10);
    kmaxy = findZero(EwaldErrorFunction(boxVectors[1][1], alpha, tol), 10);
    kmaxz = findZero(EwaldErrorFunction(boxVectors[2][2], alpha, tol), 10);
    if (kmaxx%2 == 0)
        kmaxx++;
    if (kmaxy%2 == 0)
        kmaxy++;
    if (kmaxz%2 == 0)
        kmaxz++;
}

void NonbondedForceImpl::calcPMEParameters(const System& system, const NonbondedForce& force, double& alpha, int& xsize, int& ysize, int& zsize) {
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    double tol = force.getEwaldErrorTolerance();
    alpha = (1.0/force.getCutoffDistance())*std::sqrt(-log(2.0*tol));
    xsize = (int) ceil(alpha*boxVectors[0][0]/pow(0.5*tol, 0.2));
    ysize = (int) ceil(alpha*boxVectors[1][1]/pow(0.5*tol, 0.2));
    zsize = (int) ceil(alpha*boxVectors[2][2]/pow(0.5*tol, 0.2));
}

double NonbondedForceImpl::findZero(const NonbondedForceImpl::ErrorFunction& f, int initialGuess) {
    int arg = initialGuess;
    double value = f.getValue(arg);
    if (value > 0.0) {
        while (value > 0.0 && arg > 0)
            value = f.getValue(--arg);
        return arg+1;
    }
    while (value < 0.0)
        value = f.getValue(++arg);
    return arg;
}
