/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Lee-Ping Wang                                      *
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

#include "openmm/internal/MonteCarloAnisotropicBarostatImpl.h"
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

MonteCarloAnisotropicBarostatImpl::MonteCarloAnisotropicBarostatImpl(const MonteCarloAnisotropicBarostat& owner) : owner(owner), step(0) {
}

void MonteCarloAnisotropicBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner, 3);
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    for (int i=0; i<3; i++) {
        volumeScale[i] = 0.01*volume;
        numAttempted[i] = 0;
        numAccepted[i] = 0;
    }
    SimTKOpenMMUtilities::setRandomNumberSeed(owner.getRandomNumberSeed());
}

void MonteCarloAnisotropicBarostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    if (!owner.getScaleX() && !owner.getScaleY() && !owner.getScaleZ())
        return;
    step = 0;
    
    // Compute the current potential energy.
    
    int groups = context.getIntegrator().getIntegrationForceGroups();
    double initialEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double pressure;
    
    // Choose which axis to modify at random.
    int axis;
    while (true) {
        double rnd = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()*3.0;
        if (rnd < 1.0) {
            if (owner.getScaleX()) {
                axis = 0;
                pressure = context.getParameter(MonteCarloAnisotropicBarostat::PressureX())*(AVOGADRO*1e-25);
                break;
            }
        } else if (rnd < 2.0) {
            if (owner.getScaleY()) {
                axis = 1;
                pressure = context.getParameter(MonteCarloAnisotropicBarostat::PressureY())*(AVOGADRO*1e-25);
                break;
            }
        } else if (owner.getScaleZ()) {
            axis = 2;
            pressure = context.getParameter(MonteCarloAnisotropicBarostat::PressureZ())*(AVOGADRO*1e-25);
            break;
        }
    }
    
    // Modify the periodic box size.
    
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double deltaVolume = volumeScale[axis]*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    double newVolume = volume+deltaVolume;
    Vec3 lengthScale(1.0, 1.0, 1.0);
    lengthScale[axis] = newVolume/volume;
    kernel.getAs<ApplyMonteCarloBarostatKernel>().saveCoordinates(context);
    context.getOwner().setPeriodicBoxVectors(Vec3(box[0][0]*lengthScale[0], box[0][1]*lengthScale[1], box[0][2]*lengthScale[2]),
                                             Vec3(box[1][0]*lengthScale[0], box[1][1]*lengthScale[1], box[1][2]*lengthScale[2]),
                                             Vec3(box[2][0]*lengthScale[0], box[2][1]*lengthScale[1], box[2][2]*lengthScale[2]));
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale[0], lengthScale[1], lengthScale[2]);

    // Compute the energy of the modified system.
    
    double finalEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double kT = BOLTZ*context.getParameter(MonteCarloAnisotropicBarostat::Temperature());
    double w = finalEnergy-initialEnergy + pressure*deltaVolume - context.getMolecules().size()*kT*log(newVolume/volume);
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > exp(-w/kT)) {
        // Reject the step.
        
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
    }
    else {
        numAccepted[axis]++;
        forcesInvalid = true;
    }
    numAttempted[axis]++;
    if (numAttempted[axis] >= 10) {
        if (numAccepted[axis] < 0.25*numAttempted[axis]) {
            volumeScale[axis] /= 1.1;
            numAttempted[axis] = 0;
            numAccepted[axis] = 0;
        }
        else if (numAccepted[axis] > 0.75*numAttempted[axis]) {
            volumeScale[axis] = min(volumeScale[axis]*1.1, volume*0.3);
            numAttempted[axis] = 0;
            numAccepted[axis] = 0;
        }
    }
}

map<string, double> MonteCarloAnisotropicBarostatImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters[MonteCarloAnisotropicBarostat::PressureX()] = getOwner().getDefaultPressure()[0];
    parameters[MonteCarloAnisotropicBarostat::PressureY()] = getOwner().getDefaultPressure()[1];
    parameters[MonteCarloAnisotropicBarostat::PressureZ()] = getOwner().getDefaultPressure()[2];
    parameters[MonteCarloAnisotropicBarostat::Temperature()] = getOwner().getDefaultTemperature();
    return parameters;
}

vector<string> MonteCarloAnisotropicBarostatImpl::getKernelNames() {
    vector<string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}

Vec3 MonteCarloAnisotropicBarostatImpl::computeCurrentPressure(ContextImpl& context) {
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double delta = 1e-3;
    int groups = context.getIntegrator().getIntegrationForceGroups();
    kernel.getAs<ApplyMonteCarloBarostatKernel>().saveCoordinates(context);
    vector<double> ke;
    kernel.getAs<ApplyMonteCarloBarostatKernel>().computeKineticEnergy(context, ke);
    double deltaVolume = volume*2*delta;
    Vec3 pressure;

    // Compute the pressure along each axis.

    for (int axis = 0; axis < 3; axis++) {
        // Compute the first energy.

        Vec3 scale1(1, 1, 1);
        scale1[axis] = 1.0+delta;
        context.getOwner().setPeriodicBoxVectors(box[0]*scale1[0], box[1]*scale1[1], box[2]*scale1[2]);
        kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, scale1[0], scale1[1], scale1[2]);
        double energy1 = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

        // Compute the second energy.

        Vec3 scale2(1, 1, 1);
        scale2[axis] = 1.0-delta;
        context.getOwner().setPeriodicBoxVectors(box[0]*scale2[0], box[1]*scale2[1], box[2]*scale2[2]);
        kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, scale2[0]/scale1[0], scale2[1]/scale1[1], scale2[2]/scale1[2]);
        double energy2 = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

        // Reset the box shape.

        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, 1/scale2[0], 1/scale2[1], 1/scale2[2]);

        // Compute the pressure.

        pressure[axis] = (2.0*ke[axis]/volume - (energy1-energy2)/deltaVolume)/(AVOGADRO*1e-25);
    }

    // Restore the context to its original state.

    kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
    return pressure;
}
