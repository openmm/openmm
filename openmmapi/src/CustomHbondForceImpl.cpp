/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

CustomHbondForceImpl::CustomHbondForceImpl(CustomHbondForce& owner) : owner(owner) {
}

CustomHbondForceImpl::~CustomHbondForceImpl() {
}

void CustomHbondForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomHbondForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.

    System& system = context.getSystem();
    vector<set<int> > exclusions(owner.getNumDonors());
    vector<double> parameters;
    int numDonorParameters = owner.getNumPerDonorParameters();
    for (int i = 0; i < owner.getNumDonors(); i++) {
        int primary, secondary;
        owner.getDonorParameters(i, primary, secondary, parameters);
        if (primary < 0 || primary >= system.getNumParticles()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal particle index for a donor: ";
            msg << primary;
            throw OpenMMException(msg.str());
        }
        if (secondary < 0 || secondary >= system.getNumParticles()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal particle index for a donor: ";
            msg << secondary;
            throw OpenMMException(msg.str());
        }
        if (parameters.size() != numDonorParameters) {
            stringstream msg;
            msg << "CustomHbondForce: Wrong number of parameters for donor ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    int numAcceptorParameters = owner.getNumPerAcceptorParameters();
    for (int i = 0; i < owner.getNumAcceptors(); i++) {
        int primary, secondary;
        owner.getAcceptorParameters(i, primary, secondary, parameters);
        if (primary < 0 || primary >= system.getNumParticles()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal particle index for an acceptor: ";
            msg << primary;
            throw OpenMMException(msg.str());
        }
        if (secondary < 0 || secondary >= system.getNumParticles()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal particle index for an acceptor: ";
            msg << secondary;
            throw OpenMMException(msg.str());
        }
        if (parameters.size() != numAcceptorParameters) {
            stringstream msg;
            msg << "CustomHbondForce: Wrong number of parameters for acceptor ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    for (int i = 0; i < owner.getNumExclusions(); i++) {
        int donor, acceptor;
        owner.getExclusionParticles(i, donor, acceptor);
        if (donor < 0 || donor >= owner.getNumDonors()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal donor index for an exclusion: ";
            msg << donor;
            throw OpenMMException(msg.str());
        }
        if (acceptor < 0 || acceptor >= owner.getNumAcceptors()) {
            stringstream msg;
            msg << "CustomHbondForce: Illegal acceptor index for an exclusion: ";
            msg << acceptor;
            throw OpenMMException(msg.str());
        }
        if (exclusions[donor].count(acceptor) > 0) {
            stringstream msg;
            msg << "CustomHbondForce: Multiple exclusions are specified for donor ";
            msg << donor;
            msg << " and acceptor ";
            msg << acceptor;
            throw OpenMMException(msg.str());
        }
        exclusions[donor].insert(acceptor);
    }
    if (owner.getNonbondedMethod() == CustomHbondForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("CustomHbondForce: The cutoff distance cannot be greater than half the periodic box size.");
    }
    dynamic_cast<CalcCustomHbondForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner);
}

void CustomHbondForceImpl::calcForces(ContextImpl& context) {
    dynamic_cast<CalcCustomHbondForceKernel&>(kernel.getImpl()).executeForces(context);
}

double CustomHbondForceImpl::calcEnergy(ContextImpl& context) {
    return dynamic_cast<CalcCustomHbondForceKernel&>(kernel.getImpl()).executeEnergy(context);
}

vector<string> CustomHbondForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomHbondForceKernel::Name());
    return names;
}

map<string, double> CustomHbondForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}
