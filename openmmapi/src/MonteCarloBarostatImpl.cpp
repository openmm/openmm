/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
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
using namespace std;

MonteCarloBarostatImpl::MonteCarloBarostatImpl(const MonteCarloBarostat& owner) : owner(owner), step(0) {
}

void MonteCarloBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner, 1);
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

    int groups = context.getIntegrator().getIntegrationForceGroups();
    double initialEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

    // Modify the periodic box size.

    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double deltaVolume = volumeScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    double newVolume = volume+deltaVolume;
    double lengthScale = pow(newVolume/volume, 1.0/3.0);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().saveCoordinates(context);
    context.getOwner().setPeriodicBoxVectors(box[0]*lengthScale, box[1]*lengthScale, box[2]*lengthScale);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale, lengthScale, lengthScale);

    // Compute the energy of the modified system.
    
    double finalEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double pressure = context.getParameter(MonteCarloBarostat::Pressure())*(AVOGADRO*1e-25);
    double kT = BOLTZ*context.getParameter(MonteCarloBarostat::Temperature());
    double w = finalEnergy-initialEnergy + pressure*deltaVolume - context.getMolecules().size()*kT*log(newVolume/volume);
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > exp(-w/kT)) {
        // Reject the step.

        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
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
            volumeScale = min(volumeScale*1.1, volume*0.3);
            numAttempted = 0;
            numAccepted = 0;
        }
    }
}

map<string, double> MonteCarloBarostatImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters[MonteCarloBarostat::Pressure()] = getOwner().getDefaultPressure();
    parameters[MonteCarloBarostat::Temperature()] = getOwner().getDefaultTemperature();
    return parameters;
}

vector<string> MonteCarloBarostatImpl::getKernelNames() {
    vector<string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}

double MonteCarloBarostatImpl::computeCurrentPressure(ContextImpl& context) {
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double delta = 1e-3;
    int groups = context.getIntegrator().getIntegrationForceGroups();
    kernel.getAs<ApplyMonteCarloBarostatKernel>().saveCoordinates(context);

    // Compute the first energy.

    double scale1 = 1.0+delta;
    context.getOwner().setPeriodicBoxVectors(box[0]*scale1, box[1]*scale1, box[2]*scale1);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, scale1, scale1, scale1);
    double energy1 = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

    // Compute the second energy.

    double scale2 = 1.0-delta;
    context.getOwner().setPeriodicBoxVectors(box[0]*scale2, box[1]*scale2, box[2]*scale2);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, scale2/scale1, scale2/scale1, scale2/scale1);
    double energy2 = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

    // Restore the context to its original state.

    context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);

    // Compute the pressure.

    vector<double> ke;
    kernel.getAs<ApplyMonteCarloBarostatKernel>().computeKineticEnergy(context, ke);
    double deltaVolume = volume*(scale1*scale1*scale1 - scale2*scale2*scale2);
    double pressure = (2.0/3.0)*ke[0]/volume - (energy1-energy2)/deltaVolume;
    return pressure/(AVOGADRO*1e-25);
}
