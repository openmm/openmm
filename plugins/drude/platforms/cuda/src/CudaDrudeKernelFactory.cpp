/* -------------------------------------------------------------------------- *
 *                              OpenMMDrude                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2012 Stanford University and the Authors.      *
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

#include "CudaDrudeKernelFactory.h"
#include "CudaDrudeKernels.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        CudaDrudeKernelFactory* factory = new CudaDrudeKernelFactory();
        platform.registerKernelFactory(CalcDrudeForceKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeLangevinStepKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeSCFStepKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerDrudeCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* CudaDrudeKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaContext& cu = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == CalcDrudeForceKernel::Name())
        return new CudaCalcDrudeForceKernel(name, platform, cu);
    if (name == IntegrateDrudeLangevinStepKernel::Name())
        return new CudaIntegrateDrudeLangevinStepKernel(name, platform, cu);
    if (name == IntegrateDrudeSCFStepKernel::Name())
        return new CudaIntegrateDrudeSCFStepKernel(name, platform, cu);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
