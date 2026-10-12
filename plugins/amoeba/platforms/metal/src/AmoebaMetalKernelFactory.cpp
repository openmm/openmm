/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "AmoebaMetalKernelFactory.h"
#include "AmoebaMetalKernels.h"
#include "AmoebaCommonKernels.h"
#include "MetalPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerPlatforms() {
#else
extern "C" void registerPlatforms() {
#endif
}

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerKernelFactories() {
#else
extern "C" void registerKernelFactories() {
#endif
    try {
        Platform& platform = Platform::getPlatformByName("Metal");
        AmoebaMetalKernelFactory* factory = new AmoebaMetalKernelFactory();
        platform.registerKernelFactory(CalcAmoebaTorsionTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaMultipoleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaGeneralizedKirkwoodForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaVdwForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaWcaDispersionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcHippoNonbondedForceKernel::Name(), factory);
    }
    catch (...) {
        // Ignore.  The Metal platform isn't available.
    }
}

extern "C" void registerAmoebaMetalKernelFactories() {
    try {
        Platform::getPlatformByName("Metal");
    }
    catch (...) {
        Platform::registerPlatform(new MetalPlatform());
    }
    registerKernelFactories();
}

KernelImpl* AmoebaMetalKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    MetalPlatform::PlatformData& data = *static_cast<MetalPlatform::PlatformData*>(context.getPlatformData());
    MetalContext& mc = *data.contexts[0];

    if (name == CalcAmoebaTorsionTorsionForceKernel::Name())
        return new CommonCalcAmoebaTorsionTorsionForceKernel(name, platform, mc, context.getSystem());

    if (name == CalcAmoebaMultipoleForceKernel::Name())
        return new MetalCalcAmoebaMultipoleForceKernel(name, platform, mc, context.getSystem());

    if (name == CalcAmoebaGeneralizedKirkwoodForceKernel::Name())
        return new CommonCalcAmoebaGeneralizedKirkwoodForceKernel(name, platform, mc, context.getSystem());

    if (name == CalcAmoebaVdwForceKernel::Name())
        return new CommonCalcAmoebaVdwForceKernel(name, platform, mc, context.getSystem());

    if (name == CalcAmoebaWcaDispersionForceKernel::Name())
        return new CommonCalcAmoebaWcaDispersionForceKernel(name, platform, mc, context.getSystem());

    if (name == CalcHippoNonbondedForceKernel::Name())
        return new MetalCalcHippoNonbondedForceKernel(name, platform, mc, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}