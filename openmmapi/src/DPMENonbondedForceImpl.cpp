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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/DPMENonbondedForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <map>
#include <sstream>
#include <algorithm>

using namespace OpenMM;
using namespace std;

DPMENonbondedForceImpl::DPMENonbondedForceImpl(const DPMENonbondedForce& owner) : owner(owner) {
}

DPMENonbondedForceImpl::~DPMENonbondedForceImpl() {
}

void DPMENonbondedForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcDPMENonbondedForceKernel::Name(), context);

    // Check for errors in the specification of exceptions.

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("DPMENonbondedForce must have exactly as many particles as the System it belongs to.");
    vector<set<int> > exceptions(owner.getNumParticles());
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        owner.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "DPMENonbondedForce: Illegal particle index for an exception: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "DPMENonbondedForce: Illegal particle index for an exception: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exceptions[particle1].count(particle2) > 0 || exceptions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "DPMENonbondedForce: Multiple exceptions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exceptions[particle1].insert(particle2);
        exceptions[particle2].insert(particle1);
    }
    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    double cutoff = owner.getCutoffDistance();
    if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
        throw OpenMMException("DPMENonbondedForce: The cutoff distance cannot be greater than half the periodic box size.");
    kernel.getAs<CalcDPMENonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double DPMENonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    bool includeDirect = ((groups&(1<<owner.getForceGroup())) != 0);
    bool includeReciprocal = includeDirect;
    if (owner.getReciprocalSpaceForceGroup() >= 0)
        includeReciprocal = ((groups&(1<<owner.getReciprocalSpaceForceGroup())) != 0);
    return kernel.getAs<CalcDPMENonbondedForceKernel>().execute(context, includeForces, includeEnergy, includeDirect, includeReciprocal);
}

std::vector<std::string> DPMENonbondedForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcDPMENonbondedForceKernel::Name());
    return names;
}

class DPMENonbondedForceImpl::ErrorFunction {
public:
    virtual double getValue(int arg) const = 0;
};

class DPMENonbondedForceImpl::EwaldErrorFunction : public ErrorFunction {
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


void DPMENonbondedForceImpl::calcPMEParameters(const System& system, const DPMENonbondedForce& force, double& alpha, int& xsize, int& ysize, int& zsize, bool dispersion) {
    if(dispersion) {
        force.getDispersionPMEParameters(alpha, xsize, ysize, zsize);
    } else {
        force.getPMEParameters(alpha, xsize, ysize, zsize);
    }
    // The electrostatic grid dimensions
    if (alpha == 0.0) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double tol = force.getEwaldErrorTolerance();
        alpha = (1.0/force.getCutoffDistance())*std::sqrt(-log(2.0*tol));
        xsize = (int) ceil(2*alpha*boxVectors[0][0]/(3*pow(tol, 0.2)));
        ysize = (int) ceil(2*alpha*boxVectors[1][1]/(3*pow(tol, 0.2)));
        zsize = (int) ceil(2*alpha*boxVectors[2][2]/(3*pow(tol, 0.2)));
        xsize = max(xsize, 5);
        ysize = max(ysize, 5);
        zsize = max(zsize, 5);
    }
}


void DPMENonbondedForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcDPMENonbondedForceKernel>().copyParametersToContext(context, owner);
}

void DPMENonbondedForceImpl::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcDPMENonbondedForceKernel>().getPMEParameters(alpha, nx, ny, nz);
}
void DPMENonbondedForceImpl::getDispersionPMEParameters(double& dalpha, int& dnx, int& dny, int& dnz) const {
    kernel.getAs<CalcDPMENonbondedForceKernel>().getDispersionPMEParameters(dalpha, dnx, dny, dnz);
}
