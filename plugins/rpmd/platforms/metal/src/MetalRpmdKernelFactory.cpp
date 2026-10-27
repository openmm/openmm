/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
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

#include "MetalRpmdKernelFactory.h"
#include "CommonRpmdKernels.h"
#include "MetalContext.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" void registerPlatforms() {
}

extern "C" void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("Metal");
        MetalRpmdKernelFactory* factory = new MetalRpmdKernelFactory();
        platform.registerKernelFactory(IntegrateRPMDStepKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" void registerRPMDMetalKernelFactories() {
    try {
        Platform::getPlatformByName("Metal");
    }
    catch (...) {
        Platform::registerPlatform(new MetalPlatform());
    }
    registerKernelFactories();
}

KernelImpl* MetalRpmdKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    MetalContext& mc = *static_cast<MetalPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == IntegrateRPMDStepKernel::Name())
        return new CommonIntegrateRPMDStepKernel(name, platform, mc);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
