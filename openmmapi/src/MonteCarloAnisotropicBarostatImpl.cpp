/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2013 Stanford University and the Authors.      *
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
#include "openmm/Context.h"
#include "openmm/kernels.h"
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

MonteCarloAnisotropicBarostatImpl::MonteCarloAnisotropicBarostatImpl(const MonteCarloAnisotropicBarostat& owner) : owner(owner), step(0) {
}

void MonteCarloAnisotropicBarostatImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(ApplyMonteCarloBarostatKernel::Name(), context);
    kernel.getAs<ApplyMonteCarloBarostatKernel>().initialize(context.getSystem(), owner);
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    for (int i=0; i<3; i++) {
        volumeScale[i] = 0.01*volume;
        numAttempted[i] = 0;
        numAccepted[i] = 0;
    }
    init_gen_rand(owner.getRandomNumberSeed(), random);
}

void MonteCarloAnisotropicBarostatImpl::updateContextState(ContextImpl& context) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    if (!owner.getScaleX() && !owner.getScaleY() && !owner.getScaleZ())
        return;
    step = 0;
    
    // Compute the current potential energy.
    
    double initialEnergy = context.getOwner().getState(State::Energy).getPotentialEnergy();
    double pressure;
    
    // Choose which axis to modify at random.
    int axis;
    while (true) {
        double rnd = genrand_real2(random)*3.0;
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
    double deltaVolume = volumeScale[axis]*2*(genrand_real2(random)-0.5);
    double newVolume = volume+deltaVolume;
    Vec3 lengthScale(1.0, 1.0, 1.0);
    lengthScale[axis] = newVolume/volume;
    kernel.getAs<ApplyMonteCarloBarostatKernel>().scaleCoordinates(context, lengthScale[0], lengthScale[1], lengthScale[2]);
    context.getOwner().setPeriodicBoxVectors(box[0]*lengthScale[0], box[1]*lengthScale[1], box[2]*lengthScale[2]);
    
    // Compute the energy of the modified system.
    
    double finalEnergy = context.getOwner().getState(State::Energy).getPotentialEnergy();
    double kT = BOLTZ*owner.getTemperature();
    double w = finalEnergy-initialEnergy + pressure*deltaVolume - context.getMolecules().size()*kT*std::log(newVolume/volume);
    if (w > 0 && genrand_real2(random) > std::exp(-w/kT)) {
        // Reject the step.
        
        kernel.getAs<ApplyMonteCarloBarostatKernel>().restoreCoordinates(context);
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
        volume = newVolume;
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
            volumeScale[axis] = std::min(volumeScale[axis]*1.1, volume*0.3);
            numAttempted[axis] = 0;
            numAccepted[axis] = 0;
        }
    }
}

std::map<std::string, double> MonteCarloAnisotropicBarostatImpl::getDefaultParameters() {
    std::map<std::string, double> parameters;
    parameters[MonteCarloAnisotropicBarostat::PressureX()] = getOwner().getDefaultPressure()[0];
    parameters[MonteCarloAnisotropicBarostat::PressureY()] = getOwner().getDefaultPressure()[1];
    parameters[MonteCarloAnisotropicBarostat::PressureZ()] = getOwner().getDefaultPressure()[2];
    return parameters;
}

std::vector<std::string> MonteCarloAnisotropicBarostatImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(ApplyMonteCarloBarostatKernel::Name());
    return names;
}

