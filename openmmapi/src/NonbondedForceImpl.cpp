/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/Messages.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <map>
#include <sstream>
#include <algorithm>

using namespace OpenMM;
using namespace std;

NonbondedForceImpl::NonbondedForceImpl(const NonbondedForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
    recipForceGroup = owner.getReciprocalSpaceForceGroup();
    if (recipForceGroup < 0)
        recipForceGroup = owner.getForceGroup();
    includeDirectSpace = owner.getIncludeDirectSpace();
}

NonbondedForceImpl::~NonbondedForceImpl() {
}

void NonbondedForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcNonbondedForceKernel::Name(), context);

    // Check for errors in the specification of exceptions.

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("NonbondedForce must have exactly as many particles as the System it belongs to.");
    if (owner.getUseSwitchingFunction()) {
        if (owner.getSwitchingDistance() < 0 || owner.getSwitchingDistance() >= owner.getCutoffDistance())
            throw OpenMMException("NonbondedForce: Switching distance must satisfy 0 <= r_switch < r_cutoff");
    }
    for (int i = 0; i < owner.getNumParticles(); i++) {
        double charge, sigma, epsilon;
        owner.getParticleParameters(i, charge, sigma, epsilon);
        if (sigma < 0)
            throw OpenMMException("NonbondedForce: sigma for a particle cannot be negative");
        if (epsilon < 0)
            throw OpenMMException("NonbondedForce: epsilon for a particle cannot be negative");
    }
    vector<set<int> > exceptions(owner.getNumParticles());
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle[2];
        double chargeProd, sigma, epsilon;
        owner.getExceptionParameters(i, particle[0], particle[1], chargeProd, sigma, epsilon);
        int minp = min(particle[0], particle[1]);
        int maxp = max(particle[0], particle[1]);
        for (int j = 0; j < 2; j++) {
            if (particle[j] < 0 || particle[j] >= owner.getNumParticles()) {
                stringstream msg;
                msg << "NonbondedForce: Illegal particle index for an exception: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (exceptions[minp].count(maxp) > 0) {
            stringstream msg;
            msg << "NonbondedForce: Multiple exceptions are specified for particles ";
            msg << particle[0];
            msg << " and ";
            msg << particle[1];
            throw OpenMMException(msg.str());
        }
        exceptions[minp].insert(maxp);
        if (sigma < 0)
            throw OpenMMException("NonbondedForce: sigma for an exception cannot be negative");
        if (epsilon < 0)
            throw OpenMMException("NonbondedForce: epsilon for an exception cannot be negative");
    }
    for (int i = 0; i < owner.getNumParticleParameterOffsets(); i++) {
        string parameter;
        int particleIndex;
        double chargeScale, sigmaScale, epsilonScale;
        owner.getParticleParameterOffset(i, parameter, particleIndex, chargeScale, sigmaScale, epsilonScale);
        if (particleIndex < 0 || particleIndex >= owner.getNumParticles()) {
            stringstream msg;
            msg << "NonbondedForce: Illegal particle index for a particle parameter offset: ";
            msg << particleIndex;
            throw OpenMMException(msg.str());
        }
    }
    for (int i = 0; i < owner.getNumExceptionParameterOffsets(); i++) {
        string parameter;
        int exceptionIndex;
        double chargeScale, sigmaScale, epsilonScale;
        owner.getExceptionParameterOffset(i, parameter, exceptionIndex, chargeScale, sigmaScale, epsilonScale);
        if (exceptionIndex < 0 || exceptionIndex >= owner.getNumExceptions()) {
            stringstream msg;
            msg << "NonbondedForce: Illegal exception index for an exception parameter offset: ";
            msg << exceptionIndex;
            throw OpenMMException(msg.str());
        }
    }
    if (owner.getNonbondedMethod() != NonbondedForce::NoCutoff && owner.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("NonbondedForce: "+Messages::cutoffTooLarge);
        if (owner.getNonbondedMethod() == NonbondedForce::Ewald && (boxVectors[1][0] != 0.0 || boxVectors[2][0] != 0.0 || boxVectors[2][1] != 0))
            throw OpenMMException("NonbondedForce: Ewald is not supported with non-rectangular boxes.  Use PME instead.");
    }
    kernel.getAs<CalcNonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double NonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    bool includeDirect = (includeDirectSpace && (groups&(1<<forceGroup)) != 0);
    bool includeReciprocal = ((groups&(1<<recipForceGroup)) != 0);
    return kernel.getAs<CalcNonbondedForceKernel>().execute(context, includeForces, includeEnergy, includeDirect, includeReciprocal);
}

map<string, double> NonbondedForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
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
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
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

void NonbondedForceImpl::calcPMEParameters(const System& system, const NonbondedForce& force, double& alpha, int& xsize, int& ysize, int& zsize, bool lj) {
    if (lj)
        force.getLJPMEParameters(alpha, xsize, ysize, zsize);
    else
        force.getPMEParameters(alpha, xsize, ysize, zsize);
    if (alpha == 0.0) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double tol = force.getEwaldErrorTolerance();
        alpha = (1.0/force.getCutoffDistance())*std::sqrt(-log(2.0*tol));
        if (lj) {
            xsize = (int) ceil(alpha*boxVectors[0][0]/(3*pow(tol, 0.2)));
            ysize = (int) ceil(alpha*boxVectors[1][1]/(3*pow(tol, 0.2)));
            zsize = (int) ceil(alpha*boxVectors[2][2]/(3*pow(tol, 0.2)));
        }
        else {
            xsize = (int) ceil(2*alpha*boxVectors[0][0]/(3*pow(tol, 0.2)));
            ysize = (int) ceil(2*alpha*boxVectors[1][1]/(3*pow(tol, 0.2)));
            zsize = (int) ceil(2*alpha*boxVectors[2][2]/(3*pow(tol, 0.2)));
        }
        xsize = max(xsize, 6);
        ysize = max(ysize, 6);
        zsize = max(zsize, 6);
    }
}

int NonbondedForceImpl::findZero(const NonbondedForceImpl::ErrorFunction& f, int initialGuess) {
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

double NonbondedForceImpl::evalIntegral(double r, double rs, double rc, double sigma) {
    // Compute the indefinite integral of the LJ interaction multiplied by the switching function.
    // This is a large and somewhat horrifying expression, though it does grow on you if you look
    // at it long enough.  Perhaps it could be simplified further, but I got tired of working on it.
    
    double A = 1/(rc-rs);
    double A2 = A*A;
    double A3 = A2*A;
    double sig2 = sigma*sigma;
    double sig6 = sig2*sig2*sig2;
    double rs2 = rs*rs;
    double rs3 = rs*rs2;
    double r2 = r*r;
    double r3 = r*r2;
    double r4 = r*r3;
    double r5 = r*r4;
    double r6 = r*r5;
    double r9 = r3*r6;
    return sig6*A3*((
        sig6*(
            + rs3*28*(6*rs2*A2 + 15*rs*A + 10)
            - r*rs2*945*(rs2*A2 + 2*rs*A + 1)
            + r2*rs*1080*(2*rs2*A2 + 3*rs*A + 1)
            - r3*420*(6*rs2*A2 + 6*rs*A + 1)
            + r4*756*(2*rs*A2 + A)
            - r5*378*A2)
        -r6*(
            + rs3*84*(6*rs2*A2 + 15*rs*A + 10)
            - r*rs2*3780*(rs2*A2 + 2*rs*A + 1)
            + r2*rs*7560*(2*rs2*A2 + 3*rs*A + 1))
        )/(252*r9)
     - log(r)*10*(6*rs2*A2 + 6*rs*A + 1)
     + r*15*(2*rs*A2 + A)
     - r2*3*A2
    );
}

double NonbondedForceImpl::calcDispersionCorrection(const System& system, const NonbondedForce& force) {
    if (force.getNonbondedMethod() == NonbondedForce::NoCutoff || force.getNonbondedMethod() == NonbondedForce::CutoffNonPeriodic)
        return 0.0;

    // Record sigma and epsilon for every particle, including the default value
    // for every offset parameter.

    vector<double> sigma(force.getNumParticles()), epsilon(force.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge;
        force.getParticleParameters(i, charge, sigma[i], epsilon[i]);
    }
    map<string, double> param;
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        param[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string parameter;
        int index;
        double chargeScale, sigmaScale, epsilonScale;
        force.getParticleParameterOffset(i, parameter, index, chargeScale, sigmaScale, epsilonScale);
        sigma[index] += param[parameter]*sigmaScale;
        epsilon[index] += param[parameter]*epsilonScale;
    }

    // Identify all particle classes (defined by sigma and epsilon), and count the number of
    // particles in each class.

    map<pair<double, double>, int> classCounts;
    for (int i = 0; i < force.getNumParticles(); i++) {
        pair<double, double> key = make_pair(sigma[i], epsilon[i]);
        map<pair<double, double>, int>::iterator entry = classCounts.find(key);
        if (entry == classCounts.end())
            classCounts[key] = 1;
        else
            entry->second++;
    }

    // Loop over all pairs of classes to compute the coefficient.

    double sum1 = 0, sum2 = 0, sum3 = 0;
    bool useSwitch = force.getUseSwitchingFunction();
    double cutoff = force.getCutoffDistance();
    double switchDist = force.getSwitchingDistance();
    for (map<pair<double, double>, int>::const_iterator entry = classCounts.begin(); entry != classCounts.end(); ++entry) {
        double sigma = entry->first.first;
        double epsilon = entry->first.second;
        double count = (double) entry->second;
        count *= (count + 1) / 2;
        double sigma2 = sigma*sigma;
        double sigma6 = sigma2*sigma2*sigma2;
        sum1 += count*epsilon*sigma6*sigma6;
        sum2 += count*epsilon*sigma6;
        if (useSwitch)
            sum3 += count*epsilon*(evalIntegral(cutoff, switchDist, cutoff, sigma)-evalIntegral(switchDist, switchDist, cutoff, sigma));
    }
    for (map<pair<double, double>, int>::const_iterator class1 = classCounts.begin(); class1 != classCounts.end(); ++class1)
        for (map<pair<double, double>, int>::const_iterator class2 = classCounts.begin(); class2 != class1; ++class2) {
            double sigma = 0.5*(class1->first.first+class2->first.first);
            double epsilon = sqrt(class1->first.second*class2->first.second);
            double count = (double) class1->second;
            count *= (double) class2->second;
            double sigma2 = sigma*sigma;
            double sigma6 = sigma2*sigma2*sigma2;
            sum1 += count*epsilon*sigma6*sigma6;
            sum2 += count*epsilon*sigma6;
            if (useSwitch)
                sum3 += count*epsilon*(evalIntegral(cutoff, switchDist, cutoff, sigma)-evalIntegral(switchDist, switchDist, cutoff, sigma));
        }
    double numParticles = (double) system.getNumParticles();
    double numInteractions = (numParticles*(numParticles+1))/2;
    sum1 /= numInteractions;
    sum2 /= numInteractions;
    sum3 /= numInteractions;
    return 8*numParticles*numParticles*M_PI*(sum1/(9*pow(cutoff, 9))-sum2/(3*pow(cutoff, 3))+sum3);
}

void NonbondedForceImpl::updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle, int firstException, int lastException) {
    kernel.getAs<CalcNonbondedForceKernel>().copyParametersToContext(context, owner, firstParticle, lastParticle, firstException, lastException);
    context.systemChanged();
}

void NonbondedForceImpl::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcNonbondedForceKernel>().getPMEParameters(alpha, nx, ny, nz);
}

void NonbondedForceImpl::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcNonbondedForceKernel>().getLJPMEParameters(alpha, nx, ny, nz);
}
