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

#include "OpenCLKernelFactory.h"
#include "OpenCLKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

KernelImpl* OpenCLKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    OpenCLPlatform::PlatformData& data = *static_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcForcesAndEnergyKernel::Name())
        return new OpenCLCalcForcesAndEnergyKernel(name, platform, data);
    if (name == UpdateStateDataKernel::Name())
        return new OpenCLUpdateStateDataKernel(name, platform, data);
    if (name == CalcHarmonicBondForceKernel::Name())
        return new OpenCLCalcHarmonicBondForceKernel(name, platform, data, context.getSystem());
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new OpenCLCalcHarmonicAngleForceKernel(name, platform, data, context.getSystem());
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new OpenCLCalcPeriodicTorsionForceKernel(name, platform, data, context.getSystem());
    if (name == CalcRBTorsionForceKernel::Name())
        return new OpenCLCalcRBTorsionForceKernel(name, platform, data, context.getSystem());
//    if (name == CalcNonbondedForceKernel::Name())
//        return new OpenCLCalcNonbondedForceKernel(name, platform, data, context.getSystem());
//    if (name == CalcCustomNonbondedForceKernel::Name())
//        return new OpenCLCalcCustomNonbondedForceKernel(name, platform, data, context.getSystem());
//    if (name == CalcGBSAOBCForceKernel::Name())
//        return new OpenCLCalcGBSAOBCForceKernel(name, platform, data);
    if (name == IntegrateVerletStepKernel::Name())
        return new OpenCLIntegrateVerletStepKernel(name, platform, data);
//    if (name == IntegrateLangevinStepKernel::Name())
//        return new OpenCLIntegrateLangevinStepKernel(name, platform, data);
//    if (name == IntegrateBrownianStepKernel::Name())
//        return new OpenCLIntegrateBrownianStepKernel(name, platform, data);
//    if (name == IntegrateVariableVerletStepKernel::Name())
//        return new OpenCLIntegrateVariableVerletStepKernel(name, platform, data);
//    if (name == IntegrateVariableLangevinStepKernel::Name())
//        return new OpenCLIntegrateVariableLangevinStepKernel(name, platform, data);
//    if (name == ApplyAndersenThermostatKernel::Name())
//        return new OpenCLApplyAndersenThermostatKernel(name, platform, data);
    if (name == CalcKineticEnergyKernel::Name())
        return new OpenCLCalcKineticEnergyKernel(name, platform, data);
//    if (name == RemoveCMMotionKernel::Name())
//        return new OpenCLRemoveCMMotionKernel(name, platform, data);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
