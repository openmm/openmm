/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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
#include "openmm/internal/ConstantPotentialForceImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <map>
#include <sstream>
#include <algorithm>

using namespace OpenMM;
using namespace std;

ConstantPotentialForceImpl::ConstantPotentialForceImpl(const ConstantPotentialForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

ConstantPotentialForceImpl::~ConstantPotentialForceImpl() {
}

void ConstantPotentialForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcConstantPotentialForceKernel::Name(), context);

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("ConstantPotentialForce must have exactly as many particles as the System it belongs to.");

    // Check for errors in the specification of exceptions.
    vector<set<int> > exceptions(owner.getNumParticles());
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle[2];
        double chargeProd;
        owner.getExceptionParameters(i, particle[0], particle[1], chargeProd);
        int minp = min(particle[0], particle[1]);
        int maxp = max(particle[0], particle[1]);
        for (int j = 0; j < 2; j++) {
            if (particle[j] < 0 || particle[j] >= owner.getNumParticles()) {
                stringstream msg;
                msg << "ConstantPotentialForce: Illegal particle index for an exception: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (exceptions[minp].count(maxp) > 0) {
            stringstream msg;
            msg << "ConstantPotentialForce: Multiple exceptions are specified for particles ";
            msg << particle[0];
            msg << " and ";
            msg << particle[1];
            throw OpenMMException(msg.str());
        }
        exceptions[minp].insert(maxp);
    }

    // Check for problems with the periodic box vectors.
    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    double cutoff = owner.getCutoffDistance();
    if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
        throw OpenMMException("ConstantPotentialForce: "+Messages::cutoffTooLarge);

    // Check for errors in the specification of electrodes.
    vector<int> electrodeIndexMap(owner.getNumParticles(), -1);
    vector<int> electrodeIndices;
    for (int electrode = 0; electrode < owner.getNumElectrodes(); electrode++) {
        std::set<int> particles;
        double potential;
        double gaussianWidth;
        double thomasFermiScale;
        owner.getElectrodeParameters(electrode, particles, potential, gaussianWidth, thomasFermiScale);

        for (int particle : particles) {
            if (electrodeIndexMap[particle] != -1) {
                stringstream msg;
                msg << "ConstantPotentialForce: Particle " << particle << " belongs to electrodes " << electrodeIndexMap[particle] << " and " << electrode;
                throw OpenMMException(msg.str());
            }
            electrodeIndexMap[particle] = electrode;
            electrodeIndices.push_back(particle);
        }

        if (gaussianWidth < 0) {
            stringstream msg;
            msg << "ConstantPotentialForce: Electrode " << electrode << " has negative Gaussian width";
            throw OpenMMException(msg.str());
        }
        if (thomasFermiScale < 0) {
            stringstream msg;
            msg << "ConstantPotentialForce: Electrode " << electrode << " has negative Thomas-Fermi scale";
            throw OpenMMException(msg.str());
        }
    }

    // Make sure no exceptions involve electrode particles.
    for (int i = 0; i < owner.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd;
        owner.getExceptionParameters(i, particle1, particle2, chargeProd);
        int electrode1 = electrodeIndexMap[particle1];
        int electrode2 = electrodeIndexMap[particle2];
        if (electrode1 != -1) {
            stringstream msg;
            msg << "ConstantPotentialForce: Particle " << particle1 << " belongs to exception " << i << " and electrode " << electrode1;
            throw OpenMMException(msg.str());
        }
        if (electrode2 != -1) {
            stringstream msg;
            msg << "ConstantPotentialForce: Particle " << particle2 << " belongs to exception " << i << " and electrode " << electrode2;
            throw OpenMMException(msg.str());
        }
    }

    // Make sure that we do not try to apply a constraint when zero electrode
    // particles are present.
    if (!electrodeIndices.size() && owner.getUseChargeConstraint()) {
        throw OpenMMException("At least one electrode particle must exist to apply a total charge constraint");
    }

    // If we are precomputing a matrix, electrode particles must be frozen.
    if (owner.getConstantPotentialMethod() == ConstantPotentialForce::ConstantPotentialMethod::Matrix) {
        for (int i = 0; i < electrodeIndices.size(); i++) {
            if (system.getParticleMass(electrodeIndices[i]) != 0.0) {
                throw OpenMMException("Electrode particles must have zero mass for the matrix method");
            }
        }
    }

    kernel.getAs<CalcConstantPotentialForceKernel>().initialize(context.getSystem(), owner);
}

double ConstantPotentialForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups & (1 << forceGroup)) != 0) {
        return kernel.getAs<CalcConstantPotentialForceKernel>().execute(context, includeForces, includeEnergy);
    }
    return 0.0;
}

std::vector<std::string> ConstantPotentialForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcConstantPotentialForceKernel::Name());
    return names;
}

void ConstantPotentialForceImpl::updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle, int firstException, int lastException, int firstElectrode, int lastElectrode) {
    kernel.getAs<CalcConstantPotentialForceKernel>().copyParametersToContext(context, owner, firstParticle, lastParticle, firstException, lastException, firstElectrode, lastElectrode);
    context.systemChanged();
}

void ConstantPotentialForceImpl::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    kernel.getAs<CalcConstantPotentialForceKernel>().getPMEParameters(alpha, nx, ny, nz);
}

void ConstantPotentialForceImpl::getCharges(ContextImpl& context, std::vector<double>& charges) {
    kernel.getAs<CalcConstantPotentialForceKernel>().getCharges(context, charges);
}

void ConstantPotentialForceImpl::calcPMEParameters(const System& system, const ConstantPotentialForce& force, double& alpha, int& xsize, int& ysize, int& zsize) {
    force.getPMEParameters(alpha, xsize, ysize, zsize);
    if (alpha == 0.0) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double tol = force.getEwaldErrorTolerance();
        alpha = (1.0/force.getCutoffDistance())*std::sqrt(-log(2.0*tol));
        xsize = max((int) ceil(2*alpha*boxVectors[0][0]/(3*pow(tol, 0.2))), 6);
        ysize = max((int) ceil(2*alpha*boxVectors[1][1]/(3*pow(tol, 0.2))), 6);
        zsize = max((int) ceil(2*alpha*boxVectors[2][2]/(3*pow(tol, 0.2))), 6);
    }
}
