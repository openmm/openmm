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

#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

KernelImpl* ReferenceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcForcesAndEnergyKernel::Name())
        return new ReferenceCalcForcesAndEnergyKernel(name, platform);
    if (name == UpdateStateDataKernel::Name())
        return new ReferenceUpdateStateDataKernel(name, platform, data);
    if (name == ApplyConstraintsKernel::Name())
        return new ReferenceApplyConstraintsKernel(name, platform, data);
    if (name == VirtualSitesKernel::Name())
        return new ReferenceVirtualSitesKernel(name, platform);
    if (name == CalcNonbondedForceKernel::Name())
        return new ReferenceCalcNonbondedForceKernel(name, platform);
    if (name == CalcCustomNonbondedForceKernel::Name())
        return new ReferenceCalcCustomNonbondedForceKernel(name, platform);
    if (name == CalcHarmonicBondForceKernel::Name())
        return new ReferenceCalcHarmonicBondForceKernel(name, platform);
    if (name == CalcCustomBondForceKernel::Name())
        return new ReferenceCalcCustomBondForceKernel(name, platform);
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new ReferenceCalcHarmonicAngleForceKernel(name, platform);
    if (name == CalcCustomAngleForceKernel::Name())
        return new ReferenceCalcCustomAngleForceKernel(name, platform);
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new ReferenceCalcPeriodicTorsionForceKernel(name, platform);
    if (name == CalcRBTorsionForceKernel::Name())
        return new ReferenceCalcRBTorsionForceKernel(name, platform);
    if (name == CalcCMAPTorsionForceKernel::Name())
        return new ReferenceCalcCMAPTorsionForceKernel(name, platform);
    if (name == CalcCustomTorsionForceKernel::Name())
        return new ReferenceCalcCustomTorsionForceKernel(name, platform);
    if (name == CalcGBSAOBCForceKernel::Name())
        return new ReferenceCalcGBSAOBCForceKernel(name, platform);
    if (name == CalcCustomGBForceKernel::Name())
        return new ReferenceCalcCustomGBForceKernel(name, platform);
    if (name == CalcCustomExternalForceKernel::Name())
        return new ReferenceCalcCustomExternalForceKernel(name, platform);
    if (name == CalcCustomHbondForceKernel::Name())
        return new ReferenceCalcCustomHbondForceKernel(name, platform);
    if (name == CalcCustomCentroidBondForceKernel::Name())
        return new ReferenceCalcCustomCentroidBondForceKernel(name, platform);
    if (name == CalcCustomCompoundBondForceKernel::Name())
        return new ReferenceCalcCustomCompoundBondForceKernel(name, platform);
    if (name == CalcCustomCVForceKernel::Name())
        return new ReferenceCalcCustomCVForceKernel(name, platform);
    if (name == CalcRMSDForceKernel::Name())
        return new ReferenceCalcRMSDForceKernel(name, platform);
    if (name == CalcCustomManyParticleForceKernel::Name())
        return new ReferenceCalcCustomManyParticleForceKernel(name, platform);
    if (name == CalcGayBerneForceKernel::Name())
        return new ReferenceCalcGayBerneForceKernel(name, platform);
    if (name == IntegrateVerletStepKernel::Name())
        return new ReferenceIntegrateVerletStepKernel(name, platform, data);
    if (name == IntegrateLangevinStepKernel::Name())
        return new ReferenceIntegrateLangevinStepKernel(name, platform, data);
    if (name == IntegrateBrownianStepKernel::Name())
        return new ReferenceIntegrateBrownianStepKernel(name, platform, data);
    if (name == IntegrateVariableLangevinStepKernel::Name())
        return new ReferenceIntegrateVariableLangevinStepKernel(name, platform, data);
    if (name == IntegrateVariableVerletStepKernel::Name())
        return new ReferenceIntegrateVariableVerletStepKernel(name, platform, data);
    if (name == IntegrateCustomStepKernel::Name())
        return new ReferenceIntegrateCustomStepKernel(name, platform, data);
    if (name == ApplyAndersenThermostatKernel::Name())
        return new ReferenceApplyAndersenThermostatKernel(name, platform);
    if (name == ApplyMonteCarloBarostatKernel::Name())
        return new ReferenceApplyMonteCarloBarostatKernel(name, platform);
    if (name == RemoveCMMotionKernel::Name())
        return new ReferenceRemoveCMMotionKernel(name, platform, data);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str());
}
