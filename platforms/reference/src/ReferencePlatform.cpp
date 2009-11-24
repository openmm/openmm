/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "ReferencePlatform.h"
#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMRealType.h"

using namespace OpenMM;

ReferencePlatform::ReferencePlatform() {
    ReferenceKernelFactory* factory = new ReferenceKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcGBVIForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(CalcKineticEnergyKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
}

bool ReferencePlatform::supportsDoublePrecision() const {
    return (sizeof(RealOpenMM) >= sizeof(double));
}

void ReferencePlatform::contextCreated(ContextImpl& context) const {
    context.setPlatformData(new PlatformData(context.getSystem().getNumParticles()));
}

void ReferencePlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

ReferencePlatform::PlatformData::PlatformData(int numParticles) : time(0.0), stepCount(0), numParticles(numParticles) {
    RealOpenMM** positions = new RealOpenMM*[numParticles];
    RealOpenMM** velocities = new RealOpenMM*[numParticles];
    RealOpenMM** forces = new RealOpenMM*[numParticles];
    for (int i = 0; i < numParticles; ++i) {
        positions[i] = new RealOpenMM[3];
        velocities[i] = new RealOpenMM[3];
        forces[i] = new RealOpenMM[3];
    }
    this->positions = positions;
    this->velocities = velocities;
    this->forces = forces;
}

ReferencePlatform::PlatformData::~PlatformData() {
    RealOpenMM** positions = (RealOpenMM**) this->positions;
    RealOpenMM** velocities = (RealOpenMM**) this->velocities;
    RealOpenMM** forces = (RealOpenMM**) this->forces;
    for (int i = 0; i < numParticles; ++i) {
        delete[] positions[i];
        delete[] velocities[i];
        delete[] forces[i];
    }
    delete[] positions;
    delete[] velocities;
    delete[] forces;
}
