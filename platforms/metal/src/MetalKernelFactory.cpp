/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/src/OpenCLKernelFactory.cpp                      *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalKernelFactory.h"
#include "MetalKernels.h"
#include "MetalPlatform.h"
#include "openmm/common/CommonCalcCustomGBForceKernel.h"
#include "openmm/common/CommonCalcCustomHbondForceKernel.h"
#include "openmm/common/CommonCalcCustomManyParticleForceKernel.h"
#include "openmm/common/CommonCalcCustomNonbondedForceKernel.h"
#include "openmm/common/CommonIntegrateCustomStepKernel.h"
#include "openmm/common/CommonIntegrateNoseHooverStepKernel.h"
#include "openmm/common/CommonIntegrateQTBStepKernel.h"
#include "openmm/common/CommonMinimizeKernel.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

KernelImpl* MetalKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    MetalPlatform::PlatformData& data = *static_cast<MetalPlatform::PlatformData*>(context.getPlatformData());
    MetalContext& metal = *data.computeContext;
    if (name == CalcForcesAndEnergyKernel::Name())
        return new MetalCalcForcesAndEnergyKernel(name, platform, metal);
    if (name == UpdateStateDataKernel::Name())
        return new MetalUpdateStateDataKernel(name, platform, metal);
    if (name == ApplyConstraintsKernel::Name())
        return new CommonApplyConstraintsKernel(name, platform, metal);
    if (name == VirtualSitesKernel::Name())
        return new CommonVirtualSitesKernel(name, platform, metal);
    if (name == MinimizeKernel::Name())
        return new MetalMinimizeKernel(name, platform, metal);
    if (name == CalcHarmonicBondForceKernel::Name())
        return new CommonCalcHarmonicBondForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomBondForceKernel::Name())
        return new CommonCalcCustomBondForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new CommonCalcHarmonicAngleForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomAngleForceKernel::Name())
        return new CommonCalcCustomAngleForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new CommonCalcPeriodicTorsionForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcRBTorsionForceKernel::Name())
        return new CommonCalcRBTorsionForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCMAPTorsionForceKernel::Name())
        return new CommonCalcCMAPTorsionForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomTorsionForceKernel::Name())
        return new CommonCalcCustomTorsionForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcNonbondedForceKernel::Name())
        return new MetalCalcNonbondedForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcConstantPotentialForceKernel::Name())
        return new MetalCalcConstantPotentialForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomNonbondedForceKernel::Name())
        return new CommonCalcCustomNonbondedForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcGBSAOBCForceKernel::Name())
        return new CommonCalcGBSAOBCForceKernel(name, platform, metal);
    if (name == CalcCustomGBForceKernel::Name())
        return new CommonCalcCustomGBForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomExternalForceKernel::Name())
        return new CommonCalcCustomExternalForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomHbondForceKernel::Name())
        return new CommonCalcCustomHbondForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomCentroidBondForceKernel::Name())
        return new CommonCalcCustomCentroidBondForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomCompoundBondForceKernel::Name())
        return new CommonCalcCustomCompoundBondForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcCustomCVForceKernel::Name())
        return new MetalCalcCustomCVForceKernel(name, platform, metal);
    if (name == CalcATMForceKernel::Name())
        return new MetalCalcATMForceKernel(name, platform, metal);
    if (name == CalcCustomCPPForceKernel::Name())
        return new CommonCalcCustomCPPForceKernel(name, platform, context, metal);
    if (name == CalcOrientationRestraintForceKernel::Name())
        return new CommonCalcOrientationRestraintForceKernel(name, platform, metal);
    if (name == CalcPythonForceKernel::Name())
        return new CommonCalcPythonForceKernel(name, platform, context, metal);
    if (name == CalcRGForceKernel::Name())
        return new CommonCalcRGForceKernel(name, platform, metal);
    if (name == CalcRMSDForceKernel::Name())
        return new CommonCalcRMSDForceKernel(name, platform, metal);
    if (name == CalcCustomManyParticleForceKernel::Name())
        return new CommonCalcCustomManyParticleForceKernel(name, platform, metal, context.getSystem());
    if (name == CalcGayBerneForceKernel::Name())
        return new CommonCalcGayBerneForceKernel(name, platform, metal);
    if (name == CalcLCPOForceKernel::Name())
        return new CommonCalcLCPOForceKernel(name, platform, metal);
    if (name == IntegrateVerletStepKernel::Name())
        return new CommonIntegrateVerletStepKernel(name, platform, metal);
    if (name == IntegrateLangevinMiddleStepKernel::Name())
        return new CommonIntegrateLangevinMiddleStepKernel(name, platform, metal);
    if (name == IntegrateBrownianStepKernel::Name())
        return new CommonIntegrateBrownianStepKernel(name, platform, metal);
    if (name == IntegrateVariableVerletStepKernel::Name())
        return new CommonIntegrateVariableVerletStepKernel(name, platform, metal);
    if (name == IntegrateVariableLangevinStepKernel::Name())
        return new CommonIntegrateVariableLangevinStepKernel(name, platform, metal);
    if (name == IntegrateCustomStepKernel::Name())
        return new CommonIntegrateCustomStepKernel(name, platform, metal);
    if (name == IntegrateDPDStepKernel::Name())
        return new CommonIntegrateDPDStepKernel(name, platform, metal);
    if (name == IntegrateQTBStepKernel::Name())
        return new CommonIntegrateQTBStepKernel(name, platform, metal);
    if (name == ApplyAndersenThermostatKernel::Name())
        return new CommonApplyAndersenThermostatKernel(name, platform, metal);
    if (name == IntegrateNoseHooverStepKernel::Name())
        return new CommonIntegrateNoseHooverStepKernel(name, platform, metal);
    if (name == ApplyMonteCarloBarostatKernel::Name())
        return new CommonApplyMonteCarloBarostatKernel(name, platform, metal);
    if (name == RemoveCMMotionKernel::Name())
        return new CommonRemoveCMMotionKernel(name, platform, metal);
    throw OpenMMException("Tried to create kernel with illegal kernel name '"+name+"'");
}
