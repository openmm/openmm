/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2026 Stanford University and the Authors.      *
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

#include "openmm/internal/RPMDMonteCarloAnisotropicBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Context.h"
#include "openmm/kernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/RPMDIntegrator.h"
#include "SimTKOpenMMUtilities.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace std;

RPMDMonteCarloAnisotropicBarostatImpl::RPMDMonteCarloAnisotropicBarostatImpl(const RPMDMonteCarloAnisotropicBarostat& owner) : owner(owner), step(0) {
}

void RPMDMonteCarloAnisotropicBarostatImpl::initialize(ContextImpl& context) {
    RPMDIntegrator* integrator = dynamic_cast<RPMDIntegrator*>(&context.getIntegrator());
    if (integrator == NULL)
        throw OpenMMException("RPMDMonteCarloAnisotropicBarostat must be used with an RPMDIntegrator");;
    if (!integrator->getApplyThermostat())
        throw OpenMMException("RPMDMonteCarloAnisotropicBarostat requires the integrator's thermostat to be enabled");;
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner, 1, owner.getScaleMoleculesAsRigid());
    savedPositions.resize(integrator->getNumCopies());
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    for (int i = 0; i < 3; i++) {
        volumeScale[i] = 0.01*volume;
        numAttempted[i] = 0;
        numAccepted[i] = 0;
    }
    SimTKOpenMMUtilities::setRandomNumberSeed(owner.getRandomNumberSeed());
}

void RPMDMonteCarloAnisotropicBarostatImpl::updateRPMDState(ContextImpl& context) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    if (!owner.getScaleX() && !owner.getScaleY() && !owner.getScaleZ())
        return;
    step = 0;

    // Record the initial positions and energy

    RPMDIntegrator& integrator = dynamic_cast<RPMDIntegrator&>(context.getIntegrator());
    int groups = integrator.getIntegrationForceGroups();
    double initialEnergy = 0;
    int numCopies = integrator.getNumCopies();
    for (int i = 0; i < numCopies; i++) {
        State state = integrator.getState(i, State::Positions | State::Energy, false, groups);
        savedPositions[i] = state.getPositions();
        initialEnergy += state.getPotentialEnergy();
    }

    // Compute the centroid.

    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3> centroid(numParticles, Vec3());
    for (int i = 0; i < numParticles; i++) {
        for (int j = 0; j < numCopies; j++)
            centroid[i] += savedPositions[j][i];
        centroid[i] *= 1.0/numCopies;
    }

    // Choose which axis to modify at random.

    double pressure;
    int axis;
    while (true) {
        double rnd = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()*3.0;
        if (rnd < 1.0) {
            if (owner.getScaleX()) {
                axis = 0;
                pressure = context.getParameter(RPMDMonteCarloAnisotropicBarostat::PressureX())*(AVOGADRO*1e-25);
                break;
            }
        }
        else if (rnd < 2.0) {
            if (owner.getScaleY()) {
                axis = 1;
                pressure = context.getParameter(RPMDMonteCarloAnisotropicBarostat::PressureY())*(AVOGADRO*1e-25);
                break;
            }
        }
        else if (owner.getScaleZ()) {
            axis = 2;
            pressure = context.getParameter(RPMDMonteCarloAnisotropicBarostat::PressureZ())*(AVOGADRO*1e-25);
            break;
        }
    }

    // Modify the periodic box size and scale the coordinates of the centroid.

    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double deltaVolume = volumeScale[axis]*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    double newVolume = volume+deltaVolume;
    Vec3 lengthScale(1.0, 1.0, 1.0);
    lengthScale[axis] = newVolume/volume;
    context.setPositions(centroid);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().saveCoordinates(context);
    context.getOwner().setPeriodicBoxVectors(Vec3(box[0][0]*lengthScale[0], box[0][1]*lengthScale[1], box[0][2]*lengthScale[2]),
                                             Vec3(box[1][0]*lengthScale[0], box[1][1]*lengthScale[1], box[1][2]*lengthScale[2]),
                                             Vec3(box[2][0]*lengthScale[0], box[2][1]*lengthScale[1], box[2][2]*lengthScale[2]));
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale[0], lengthScale[1], lengthScale[2]);
    State scaledState = context.getOwner().getState(State::Positions);

    // Now apply the same offset to all the copies.

    vector<Vec3> delta(numParticles);
    for (int i = 0; i < numParticles; i++)
        delta[i] = scaledState.getPositions()[i]-centroid[i];
    double finalEnergy = 0;
    vector<Vec3> positions(numParticles);
    for (int copy = 0; copy < numCopies; copy++) {
        for (int i = 0; i < numParticles; i++)
            positions[i] = savedPositions[copy][i]+delta[i];
        integrator.setPositions(copy, positions);
        finalEnergy += integrator.getState(copy, State::Energy, false, groups).getPotentialEnergy();
    }

    // Compute the energy of the modified system.

    double numberOfScaledParticles;
    if (owner.getScaleMoleculesAsRigid())
        numberOfScaledParticles = context.getMolecules().size();
    else
        numberOfScaledParticles = context.getSystem().getNumParticles();
    double kT = BOLTZ*integrator.getTemperature();
    double w = (finalEnergy-initialEnergy)/numCopies + pressure*deltaVolume - numberOfScaledParticles*kT*log(newVolume/volume);
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > exp(-w/kT)) {
        // Reject the step.

        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        for (int copy = 0; copy < numCopies; copy++)
            integrator.setPositions(copy, savedPositions[copy]);
    }
    else
        numAccepted[axis]++;
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

map<string, double> RPMDMonteCarloAnisotropicBarostatImpl::getDefaultParameters() {
    return {{RPMDMonteCarloAnisotropicBarostat::PressureX(), getOwner().getDefaultPressure()[0]},
            {RPMDMonteCarloAnisotropicBarostat::PressureY(), getOwner().getDefaultPressure()[1]},
            {RPMDMonteCarloAnisotropicBarostat::PressureZ(), getOwner().getDefaultPressure()[2]}};
}

vector<string> RPMDMonteCarloAnisotropicBarostatImpl::getKernelNames() {
    return {ApplyMonteCarloBarostatKernel::Name()};
}
