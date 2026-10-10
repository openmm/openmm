/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "MetalKernelFactory.h"
#include "MetalKernels.h"
#include "MetalPlatform.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/common/CommonParallelKernels.h"
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
    MetalContext& mc = *data.contexts[0];
    if (name == CalcForcesAndEnergyKernel::Name())
        return new MetalCalcForcesAndEnergyKernel(name, platform, mc);
    if (name == UpdateStateDataKernel::Name())
        return new MetalUpdateStateDataKernel(name, platform, mc);
    if (name == ApplyConstraintsKernel::Name())
        return new CommonApplyConstraintsKernel(name, platform, mc);
    if (name == VirtualSitesKernel::Name())
        return new CommonVirtualSitesKernel(name, platform, mc);
    if (name == MinimizeKernel::Name())
        return new CommonMinimizeKernel(name, platform, mc);
    if (name == CalcHarmonicBondForceKernel::Name())
        return new CommonCalcHarmonicBondForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomBondForceKernel::Name())
        return new CommonCalcCustomBondForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new CommonCalcHarmonicAngleForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomAngleForceKernel::Name())
        return new CommonCalcCustomAngleForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new CommonCalcPeriodicTorsionForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcRBTorsionForceKernel::Name())
        return new CommonCalcRBTorsionForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCMAPTorsionForceKernel::Name())
        return new CommonCalcCMAPTorsionForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomTorsionForceKernel::Name())
        return new CommonCalcCustomTorsionForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcNonbondedForceKernel::Name())
        return new MetalCalcNonbondedForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcConstantPotentialForceKernel::Name())
        return new MetalCalcConstantPotentialForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomNonbondedForceKernel::Name())
        return new CommonCalcCustomNonbondedForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcGBSAOBCForceKernel::Name())
        return new CommonCalcGBSAOBCForceKernel(name, platform, mc);
    if (name == CalcCustomGBForceKernel::Name())
        return new CommonCalcCustomGBForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomExternalForceKernel::Name())
        return new CommonCalcCustomExternalForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomHbondForceKernel::Name())
        return new CommonCalcCustomHbondForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomCentroidBondForceKernel::Name())
        return new CommonCalcCustomCentroidBondForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomCompoundBondForceKernel::Name())
        return new CommonCalcCustomCompoundBondForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcCustomCVForceKernel::Name())
        return new MetalCalcCustomCVForceKernel(name, platform, mc);
    if (name == CalcATMForceKernel::Name())
        return new MetalCalcATMForceKernel(name, platform, mc);
    if (name == CalcCustomCPPForceKernel::Name())
        return new CommonCalcCustomCPPForceKernel(name, platform, context, mc);
    if (name == CalcOrientationRestraintForceKernel::Name())
        return new CommonCalcOrientationRestraintForceKernel(name, platform, mc);
    if (name == CalcPythonForceKernel::Name())
        return new CommonCalcPythonForceKernel(name, platform, context, mc);
    if (name == CalcRGForceKernel::Name())
        return new CommonCalcRGForceKernel(name, platform, mc);
    if (name == CalcRMSDForceKernel::Name())
        return new CommonCalcRMSDForceKernel(name, platform, mc);
    if (name == CalcCustomManyParticleForceKernel::Name())
        return new CommonCalcCustomManyParticleForceKernel(name, platform, mc, context.getSystem());
    if (name == CalcGayBerneForceKernel::Name())
        return new CommonCalcGayBerneForceKernel(name, platform, mc);
    if (name == CalcLCPOForceKernel::Name())
        return new CommonCalcLCPOForceKernel(name, platform, mc);
    if (name == IntegrateVerletStepKernel::Name())
        return new CommonIntegrateVerletStepKernel(name, platform, mc);
    if (name == IntegrateLangevinMiddleStepKernel::Name())
        return new CommonIntegrateLangevinMiddleStepKernel(name, platform, mc);
    if (name == IntegrateBrownianStepKernel::Name())
        return new CommonIntegrateBrownianStepKernel(name, platform, mc);
    if (name == IntegrateVariableVerletStepKernel::Name())
        return new CommonIntegrateVariableVerletStepKernel(name, platform, mc);
    if (name == IntegrateVariableLangevinStepKernel::Name())
        return new CommonIntegrateVariableLangevinStepKernel(name, platform, mc);
    if (name == IntegrateCustomStepKernel::Name())
        return new CommonIntegrateCustomStepKernel(name, platform, mc);
    if (name == IntegrateDPDStepKernel::Name())
        return new CommonIntegrateDPDStepKernel(name, platform, mc);
    if (name == IntegrateQTBStepKernel::Name())
        return new CommonIntegrateQTBStepKernel(name, platform, mc);
    if (name == ApplyAndersenThermostatKernel::Name())
        return new CommonApplyAndersenThermostatKernel(name, platform, mc);
    if (name == IntegrateNoseHooverStepKernel::Name())
        return new CommonIntegrateNoseHooverStepKernel(name, platform, mc);
    if (name == ApplyMonteCarloBarostatKernel::Name())
        return new CommonApplyMonteCarloBarostatKernel(name, platform, mc);
    if (name == RemoveCMMotionKernel::Name())
        return new CommonRemoveCMMotionKernel(name, platform, mc);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
