/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "ReferenceConstraints.h"
#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKOpenMMRealType.h"
#include "openmm/Vec3.h"
#include <map>
#include <vector>

using namespace OpenMM;
using namespace std;

ReferencePlatform::ReferencePlatform() {
    ReferenceKernelFactory* factory = new ReferenceKernelFactory();
    registerKernelFactory(CalcForcesAndEnergyKernel::Name(), factory);
    registerKernelFactory(UpdateStateDataKernel::Name(), factory);
    registerKernelFactory(ApplyConstraintsKernel::Name(), factory);
    registerKernelFactory(VirtualSitesKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomBondForceKernel::Name(), factory);
    registerKernelFactory(CalcHarmonicAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomAngleForceKernel::Name(), factory);
    registerKernelFactory(CalcPeriodicTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcRBTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCMAPTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomTorsionForceKernel::Name(), factory);
    registerKernelFactory(CalcNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomNonbondedForceKernel::Name(), factory);
    registerKernelFactory(CalcGBSAOBCForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomGBForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomExternalForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomHbondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCentroidBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCompoundBondForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomCVForceKernel::Name(), factory);
    registerKernelFactory(CalcRMSDForceKernel::Name(), factory);
    registerKernelFactory(CalcCustomManyParticleForceKernel::Name(), factory);
    registerKernelFactory(CalcGayBerneForceKernel::Name(), factory);
    registerKernelFactory(IntegrateVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateBrownianStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableLangevinStepKernel::Name(), factory);
    registerKernelFactory(IntegrateVariableVerletStepKernel::Name(), factory);
    registerKernelFactory(IntegrateCustomStepKernel::Name(), factory);
    registerKernelFactory(ApplyAndersenThermostatKernel::Name(), factory);
    registerKernelFactory(ApplyMonteCarloBarostatKernel::Name(), factory);
    registerKernelFactory(RemoveCMMotionKernel::Name(), factory);
}

double ReferencePlatform::getSpeed() const {
    return 1;
}

bool ReferencePlatform::supportsDoublePrecision() const {
    return true;
}

void ReferencePlatform::contextCreated(ContextImpl& context, const map<string, string>& properties) const {
    context.setPlatformData(new PlatformData(context.getSystem()));
}

void ReferencePlatform::contextDestroyed(ContextImpl& context) const {
    PlatformData* data = reinterpret_cast<PlatformData*>(context.getPlatformData());
    delete data;
}

ReferencePlatform::PlatformData::PlatformData(const System& system) : time(0.0), stepCount(0), numParticles(system.getNumParticles()) {
    positions = new vector<Vec3>(numParticles);
    velocities = new vector<Vec3>(numParticles);
    forces = new vector<Vec3>(numParticles);
    periodicBoxSize = new Vec3();
    periodicBoxVectors = new Vec3[3];
    constraints = new ReferenceConstraints(system);
    energyParameterDerivatives = new map<string, double>();
}

ReferencePlatform::PlatformData::~PlatformData() {
    delete (vector<Vec3>*) positions;
    delete (vector<Vec3>*) velocities;
    delete (vector<Vec3>*) forces;
    delete (Vec3*) periodicBoxSize;
    delete[] (Vec3*) periodicBoxVectors;
    delete (ReferenceConstraints*) constraints;
    delete (map<string, double>*) energyParameterDerivatives;
}
