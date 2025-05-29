/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2024 Stanford University and the Authors.      *
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
#include "openmm/internal/CustomTorsionForceImpl.h"
#include "openmm/kernels.h"
#include <sstream>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

CustomTorsionForceImpl::CustomTorsionForceImpl(const CustomTorsionForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
}

CustomTorsionForceImpl::~CustomTorsionForceImpl() {
}

void CustomTorsionForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomTorsionForceKernel::Name(), context);

    // Check for errors in the specification of bonds.

    const System& system = context.getSystem();
    vector<double> parameters;
    int numParameters = owner.getNumPerTorsionParameters();
    for (int i = 0; i < owner.getNumTorsions(); i++) {
        int particle[4];
        owner.getTorsionParameters(i, particle[0], particle[1], particle[2], particle[3], parameters);
        for (int j = 0; j < 4; j++) {
            if (particle[j] < 0 || particle[j] >= system.getNumParticles()) {
                stringstream msg;
                msg << "CustomTorsionForce: Illegal particle index for a torsion: ";
                msg << particle[j];
                throw OpenMMException(msg.str());
            }
        }
        if (parameters.size() != numParameters) {
            stringstream msg;
            msg << "CustomTorsionForce: Wrong number of parameters for torsion ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    kernel.getAs<CalcCustomTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double CustomTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
        return kernel.getAs<CalcCustomTorsionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomTorsionForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomTorsionForceKernel::Name());
    return names;
}

map<string, double> CustomTorsionForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

void CustomTorsionForceImpl::updateParametersInContext(ContextImpl& context, int firstTorsion, int lastTorsion) {
    kernel.getAs<CalcCustomTorsionForceKernel>().copyParametersToContext(context, owner, firstTorsion, lastTorsion);
    context.systemChanged();
}
