/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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

#include "openmm/internal/RPMDMonteCarloBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/kernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/RPMDIntegrator.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace OpenMM;
using namespace OpenMM_SFMT;
using std::vector;

const float BOLTZMANN = 1.380658e-23f; // (J/K)
const float AVOGADRO = 6.0221367e23f;
const float RGAS = BOLTZMANN*AVOGADRO; // (J/(mol K))
const float BOLTZ = RGAS/1000;         // (kJ/(mol K))

RPMDMonteCarloBarostatImpl::RPMDMonteCarloBarostatImpl(const RPMDMonteCarloBarostat& owner) : owner(owner), step(0) {
}

void RPMDMonteCarloBarostatImpl::initialize(ContextImpl& context) {
    RPMDIntegrator* integrator = dynamic_cast<RPMDIntegrator*>(&context.getIntegrator());
    if (integrator == NULL)
        throw OpenMMException("RPMDMonteCarloBarostat must be used with an RPMDIntegrator");;
    if (!integrator->getApplyThermostat())
        throw OpenMMException("RPMDMonteCarloBarostat requires the integrator's thermostat to be enabled");;
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner);
    savedPositions.resize(integrator->getNumCopies());
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    volumeScale = 0.01*volume;
    numAttempted = 0;
    numAccepted = 0;
    int randSeed = owner.getRandomNumberSeed();
    // A random seed of 0 means use a unique one
    if (randSeed == 0) randSeed = osrngseed();
    init_gen_rand(randSeed, random);
}

void RPMDMonteCarloBarostatImpl::updateRPMDState(ContextImpl& context) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    step = 0;

    // Compute the current potential energy.

    RPMDIntegrator& integrator = dynamic_cast<RPMDIntegrator&>(context.getIntegrator());

    // Record the initial positions and energy

    double initialEnergy = 0;
    int numCopies = integrator.getNumCopies();
    for (int i = 0; i < numCopies; i++) {
        State state = integrator.getState(i, State::Positions | State::Energy);
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

    // Modify the periodic box size and scale the coordinates of the centroid.

    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    double deltaVolume = volumeScale*2*(genrand_real2(random)-0.5);
    double newVolume = volume+deltaVolume;
    double lengthScale = std::pow(newVolume/volume, 1.0/3.0);
    context.setPositions(centroid);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale, lengthScale, lengthScale);
    context.getOwner().setPeriodicBoxVectors(box[0]*lengthScale, box[1]*lengthScale, box[2]*lengthScale);
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
        finalEnergy += integrator.getState(copy, State::Energy).getPotentialEnergy();
    }

    // Compute the energy of the modified system.

    double pressure = context.getParameter(RPMDMonteCarloBarostat::Pressure())*(AVOGADRO*1e-25);
    double kT = BOLTZ*integrator.getTemperature();
    double w = (finalEnergy-initialEnergy)/numCopies + pressure*deltaVolume - context.getMolecules().size()*kT*std::log(newVolume/volume);
    if (w > 0 && genrand_real2(random) > std::exp(-w/kT)) {
        // Reject the step.

        for (int copy = 0; copy < numCopies; copy++)
            integrator.setPositions(copy, savedPositions[copy]);
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        volume = newVolume;
    }
    else
        numAccepted++;
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

std::map<std::string, double> RPMDMonteCarloBarostatImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    parameters[RPMDMonteCarloBarostat::Pressure()] = getOwner().getDefaultPressure();
    return parameters;
}

std::vector<std::string> RPMDMonteCarloBarostatImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}
