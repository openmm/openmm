/* -------------------------------------------------------------------------- *
 *                              OpenMMDrude                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2019 Stanford University and the Authors.      *
 * Portions copyright (c) 2021 Advanced Micro Devices, Inc.                   *
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

#include <exception>

#include "HipDrudeKernelFactory.h"
#include "CommonDrudeKernels.h"
#include "HipContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("HIP");
        HipDrudeKernelFactory* factory = new HipDrudeKernelFactory();
        platform.registerKernelFactory(CalcDrudeForceKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeLangevinStepKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeSCFStepKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerDrudeHipKernelFactories() {
    try {
        Platform::getPlatformByName("HIP");
    }
    catch (...) {
        Platform::registerPlatform(new HipPlatform());
    }
    registerKernelFactories();
}

KernelImpl* HipDrudeKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    HipContext& cu = *static_cast<HipPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == CalcDrudeForceKernel::Name())
        return new CommonCalcDrudeForceKernel(name, platform, cu);
    if (name == IntegrateDrudeLangevinStepKernel::Name())
        return new CommonIntegrateDrudeLangevinStepKernel(name, platform, cu);
    if (name == IntegrateDrudeSCFStepKernel::Name())
        return new CommonIntegrateDrudeSCFStepKernel(name, platform, cu);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
