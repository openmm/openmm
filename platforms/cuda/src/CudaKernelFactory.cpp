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

#include "CudaKernelFactory.h"
#include "CudaKernels.h"
#include "CudaParallelKernels.h"
#include "CudaPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

KernelImpl* CudaKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    if (data.contexts.size() > 1) {
        // We are running in parallel on multiple devices, so we may want to create a parallel kernel.
        
        if (name == CalcForcesAndEnergyKernel::Name())
            return new CudaParallelCalcForcesAndEnergyKernel(name, platform, data);
        if (name == CalcHarmonicBondForceKernel::Name())
            return new CudaParallelCalcHarmonicBondForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomBondForceKernel::Name())
            return new CudaParallelCalcCustomBondForceKernel(name, platform, data, context.getSystem());
        if (name == CalcHarmonicAngleForceKernel::Name())
            return new CudaParallelCalcHarmonicAngleForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomAngleForceKernel::Name())
            return new CudaParallelCalcCustomAngleForceKernel(name, platform, data, context.getSystem());
        if (name == CalcPeriodicTorsionForceKernel::Name())
            return new CudaParallelCalcPeriodicTorsionForceKernel(name, platform, data, context.getSystem());
        if (name == CalcRBTorsionForceKernel::Name())
            return new CudaParallelCalcRBTorsionForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCMAPTorsionForceKernel::Name())
            return new CudaParallelCalcCMAPTorsionForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomTorsionForceKernel::Name())
            return new CudaParallelCalcCustomTorsionForceKernel(name, platform, data, context.getSystem());
        if (name == CalcNonbondedForceKernel::Name())
            return new CudaParallelCalcNonbondedForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomNonbondedForceKernel::Name())
            return new CudaParallelCalcCustomNonbondedForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomExternalForceKernel::Name())
            return new CudaParallelCalcCustomExternalForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomHbondForceKernel::Name())
            return new CudaParallelCalcCustomHbondForceKernel(name, platform, data, context.getSystem());
        if (name == CalcCustomCompoundBondForceKernel::Name())
            return new CudaParallelCalcCustomCompoundBondForceKernel(name, platform, data, context.getSystem());
    }
    CudaContext& cu = *data.contexts[0];
    if (name == CalcForcesAndEnergyKernel::Name())
        return new CudaCalcForcesAndEnergyKernel(name, platform, cu);
    if (name == UpdateStateDataKernel::Name())
        return new CudaUpdateStateDataKernel(name, platform, cu);
    if (name == ApplyConstraintsKernel::Name())
        return new CudaApplyConstraintsKernel(name, platform, cu);
    if (name == VirtualSitesKernel::Name())
        return new CudaVirtualSitesKernel(name, platform, cu);
    if (name == CalcHarmonicBondForceKernel::Name())
        return new CudaCalcHarmonicBondForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomBondForceKernel::Name())
        return new CudaCalcCustomBondForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcHarmonicAngleForceKernel::Name())
        return new CudaCalcHarmonicAngleForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomAngleForceKernel::Name())
        return new CudaCalcCustomAngleForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcPeriodicTorsionForceKernel::Name())
        return new CudaCalcPeriodicTorsionForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcRBTorsionForceKernel::Name())
        return new CudaCalcRBTorsionForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCMAPTorsionForceKernel::Name())
        return new CudaCalcCMAPTorsionForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomTorsionForceKernel::Name())
        return new CudaCalcCustomTorsionForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcNonbondedForceKernel::Name())
        return new CudaCalcNonbondedForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomNonbondedForceKernel::Name())
        return new CudaCalcCustomNonbondedForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcGBSAOBCForceKernel::Name())
        return new CudaCalcGBSAOBCForceKernel(name, platform, cu);
    if (name == CalcCustomGBForceKernel::Name())
        return new CudaCalcCustomGBForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomExternalForceKernel::Name())
        return new CudaCalcCustomExternalForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomHbondForceKernel::Name())
        return new CudaCalcCustomHbondForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomCentroidBondForceKernel::Name())
        return new CudaCalcCustomCentroidBondForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomCompoundBondForceKernel::Name())
        return new CudaCalcCustomCompoundBondForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcCustomCVForceKernel::Name())
        return new CudaCalcCustomCVForceKernel(name, platform, cu);
    if (name == CalcRMSDForceKernel::Name())
        return new CudaCalcRMSDForceKernel(name, platform, cu);
    if (name == CalcCustomManyParticleForceKernel::Name())
        return new CudaCalcCustomManyParticleForceKernel(name, platform, cu, context.getSystem());
    if (name == CalcGayBerneForceKernel::Name())
        return new CudaCalcGayBerneForceKernel(name, platform, cu);
    if (name == IntegrateVerletStepKernel::Name())
        return new CudaIntegrateVerletStepKernel(name, platform, cu);
    if (name == IntegrateLangevinStepKernel::Name())
        return new CudaIntegrateLangevinStepKernel(name, platform, cu);
    if (name == IntegrateBrownianStepKernel::Name())
        return new CudaIntegrateBrownianStepKernel(name, platform, cu);
    if (name == IntegrateVariableVerletStepKernel::Name())
        return new CudaIntegrateVariableVerletStepKernel(name, platform, cu);
    if (name == IntegrateVariableLangevinStepKernel::Name())
        return new CudaIntegrateVariableLangevinStepKernel(name, platform, cu);
    if (name == IntegrateCustomStepKernel::Name())
        return new CudaIntegrateCustomStepKernel(name, platform, cu);
    if (name == ApplyAndersenThermostatKernel::Name())
        return new CudaApplyAndersenThermostatKernel(name, platform, cu);
    if (name == ApplyMonteCarloBarostatKernel::Name())
        return new CudaApplyMonteCarloBarostatKernel(name, platform, cu);
    if (name == RemoveCMMotionKernel::Name())
        return new CudaRemoveCMMotionKernel(name, platform, cu);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
