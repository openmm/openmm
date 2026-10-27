/* -------------------------------------------------------------------------- *
 *                              OpenMMDrude                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2011-2026 Stanford University and the Authors.      *
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

#include "MetalDrudeKernelFactory.h"
#include "CommonDrudeKernels.h"
#include "MetalContext.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" void registerPlatforms() {
}

extern "C" void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("Metal");
        MetalDrudeKernelFactory* factory = new MetalDrudeKernelFactory();
        platform.registerKernelFactory(CalcDrudeForceKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeLangevinStepKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateDrudeSCFStepKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" void registerDrudeMetalKernelFactories() {
    try {
        Platform::getPlatformByName("Metal");
    }
    catch (...) {
        Platform::registerPlatform(new MetalPlatform());
    }
    registerKernelFactories();
}

KernelImpl* MetalDrudeKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    MetalContext& mc = *static_cast<MetalPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == CalcDrudeForceKernel::Name())
        return new CommonCalcDrudeForceKernel(name, platform, mc);
    if (name == IntegrateDrudeLangevinStepKernel::Name())
        return new CommonIntegrateDrudeLangevinStepKernel(name, platform, mc);
    if (name == IntegrateDrudeSCFStepKernel::Name())
        return new CommonIntegrateDrudeSCFStepKernel(name, platform, mc);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
