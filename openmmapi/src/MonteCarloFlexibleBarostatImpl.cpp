/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Sander Vandenhaute                                 *
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

#include "openmm/internal/MonteCarloFlexibleBarostatImpl.h"
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

MonteCarloFlexibleBarostatImpl::MonteCarloFlexibleBarostatImpl(const MonteCarloFlexibleBarostat& owner) : owner(owner), step(0) {
}

void MonteCarloFlexibleBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner, owner.getScaleMoleculesAsRigid());
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    lengthScale = 0.01*pow(volume, 1.0/3.0);
    numAttempted = 0;
    numAccepted = 0;
    SimTKOpenMMUtilities::setRandomNumberSeed(owner.getRandomNumberSeed());
}

void MonteCarloFlexibleBarostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    step = 0;

    // Compute the current potential energy.

    int groups = context.getIntegrator().getIntegrationForceGroups();
    double initialEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double pressure = context.getParameter(MonteCarloFlexibleBarostat::Pressure())*(AVOGADRO*1e-25);

    // Generate trial box vectors

    Vec3 box[3], trial[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    trial[0][0] = box[0][0] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    trial[1][0] = box[1][0] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    trial[1][1] = box[1][1] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    trial[2][0] = box[2][0] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    trial[2][1] = box[2][1] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    trial[2][2] = box[2][2] + lengthScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);

    // Recompute reduced form by flipping/linear combinations.

    for (auto i = 0; i < 3; i++)
        if (trial[i][i] < 0)
            trial[i] = -trial[i];
    trial[2] -= trial[1]*round(trial[2][1]/trial[1][1]);
    trial[2] -= trial[0]*round(trial[2][0]/trial[0][0]);
    trial[1] -= trial[0]*round(trial[1][0]/trial[0][0]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double newVolume = trial[0][0]*trial[1][1]*trial[2][2];

    // Scale particle coordinates and update box vectors in context.

    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, trial[0][0]/box[0][0], trial[1][1]/box[1][1], trial[2][2]/box[2][2]);
    context.getOwner().setPeriodicBoxVectors(trial[0], trial[1], trial[2]);

    // Compute the energy of the modified system.

    double numberOfScaledParticles;
    if (owner.getScaleMoleculesAsRigid())
        numberOfScaledParticles = context.getMolecules().size();
    else
        numberOfScaledParticles = context.getSystem().getNumParticles();
    double finalEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double kT = BOLTZ*context.getParameter(MonteCarloFlexibleBarostat::Temperature());
    double w0 = finalEnergy-initialEnergy;
    double w1 = pressure*(newVolume-volume);
    double w2 = -(numberOfScaledParticles-2)*kT*log(newVolume/volume);
    double w3 = -kT*log((trial[0][0]*trial[0][0]*trial[1][1])/(box[0][0]*box[0][0]*box[1][1]));
    double w  = w0+w1+w2+w3;
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > exp(-w/kT)) {
        // Reject the step.

        kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
    }
    else {
        numAccepted++;
        forcesInvalid = true;
    }
    numAttempted++;
    if (numAttempted >= 10) {
        if (numAccepted < 0.25*numAttempted) {
            lengthScale /= 1.1;
            numAttempted = 0;
            numAccepted = 0;
        }
        else if (numAccepted > 0.75*numAttempted) {
            lengthScale = min(lengthScale*1.1, pow(newVolume, 1.0/3.0)*0.3);
            numAttempted = 0;
            numAccepted = 0;
        }
    }

}

map<string, double> MonteCarloFlexibleBarostatImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters[MonteCarloFlexibleBarostat::Pressure()] = getOwner().getDefaultPressure();
    parameters[MonteCarloFlexibleBarostat::Temperature()] = getOwner().getDefaultTemperature();
    return parameters;
}

vector<string> MonteCarloFlexibleBarostatImpl::getKernelNames() {
    vector<string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}

