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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomBondForceImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

CustomBondForceImpl::CustomBondForceImpl(const CustomBondForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

CustomBondForceImpl::~CustomBondForceImpl() {
}

void CustomBondForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomBondForceKernel::Name(), context);

    // Check for errors in the specification of bonds.

    const System& system = context.getSystem();
    vector<double> parameters;
    int numParameters = owner.getNumPerBondParameters();
    for (int i = 0; i < owner.getNumBonds(); i++) {
        int particle[2];
        owner.getBondParameters(i, particle[0], particle[1], parameters);
        for (int j = 0; j < 2; j++) {
            if (particle[j] < 0 || particle[j] >= system.getNumParticles()) {
                stringstream msg;
                msg << "CustomBondForce: Illegal particle index for a bond: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (parameters.size() != numParameters) {
            stringstream msg;
            msg << "CustomBondForce: Wrong number of parameters for bond ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    kernel.getAs<CalcCustomBondForceKernel>().initialize(context.getSystem(), owner);
}

double CustomBondForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcCustomBondForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomBondForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomBondForceKernel::Name());
    return names;
}

map<string, double> CustomBondForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

vector<pair<int, int> > CustomBondForceImpl::getBondedParticles() const {
    int numBonds = owner.getNumBonds();
    vector<pair<int, int> > bonds(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        owner.getBondParameters(i, bonds[i].first, bonds[i].second, parameters);
    }
    return bonds;
}

void CustomBondForceImpl::updateParametersInContext(ContextImpl& context, int firstBond, int lastBond) {
    kernel.getAs<CalcCustomBondForceKernel>().copyParametersToContext(context, owner, firstBond, lastBond);
    context.systemChanged();
}
