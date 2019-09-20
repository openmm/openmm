/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2019 Stanford University and the Authors.      *
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

#include "openmm/internal/MonteCarloBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/kernels.h"
#include "openmm/OpenMMException.h"
#include "SimTKOpenMMUtilities.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace OpenMM_SFMT;
using std::vector;

MonteCarloBarostatImpl::MonteCarloBarostatImpl(const MonteCarloBarostat& owner) : owner(owner), step(0) {
}

void MonteCarloBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner);
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    volumeScale = 0.01*volume;
    numAttempted = 0;
    numAccepted = 0;
    SimTKOpenMMUtilities::setRandomNumberSeed(owner.getRandomNumberSeed());
}

void MonteCarloBarostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    step = 0;

    // Compute the current potential energy.

    double initialEnergy = context.getOwner().getState(State::Energy).getPotentialEnergy();

    // Modify the periodic box size.

    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double deltaVolume = volumeScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    double newVolume = volume+deltaVolume;
    double lengthScale = std::pow(newVolume/volume, 1.0/3.0);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale, lengthScale, lengthScale);
    context.getOwner().setPeriodicBoxVectors(box[0]*lengthScale, box[1]*lengthScale, box[2]*lengthScale);

    // Compute the energy of the modified system.
    
    double finalEnergy = context.getOwner().getState(State::Energy).getPotentialEnergy();
    double pressure = context.getParameter(MonteCarloBarostat::Pressure())*(AVOGADRO*1e-25);
    double kT = BOLTZ*context.getParameter(MonteCarloBarostat::Temperature());
    double w = finalEnergy-initialEnergy + pressure*deltaVolume - context.getMolecules().size()*kT*std::log(newVolume/volume);
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > std::exp(-w/kT)) {
        // Reject the step.

        kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        volume = newVolume;
    }
    else {
        numAccepted++;
        forcesInvalid = true;
    }
    numAttempted++;
    if (numAttempted >= 10) {
        if (numAccepted < 0.25*numAttempted) {
            volumeScale /= 1.1;
            numAttempted = 0;
            numAccepted = 0;
        }
        else if (numAccepted > 0.75*numAttempted) {
            volumeScale = std::min(volumeScale*1.1, volume*0.3);
            numAttempted = 0;
            numAccepted = 0;
        }
    }
}

std::map<std::string, double> MonteCarloBarostatImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    parameters[MonteCarloBarostat::Pressure()] = getOwner().getDefaultPressure();
    parameters[MonteCarloBarostat::Temperature()] = getOwner().getDefaultTemperature();
    return parameters;
}

std::vector<std::string> MonteCarloBarostatImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}

